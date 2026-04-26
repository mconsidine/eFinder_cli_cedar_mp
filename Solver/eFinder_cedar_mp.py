#!/usr/bin/python3

# eFinder_cedar_mp — electronic finder scope, plate-solving over LX200/WiFi
# Derived from original work Copyright (C) 2025 Keith Venables (GPL v3)
# Simplified: direct picamera2, no Nexus, no GPIO, no LED, no WiFi switching
#
# Cedar edition — replaces tetra3rs (Rust, in-process) with:
#   cedar-detect : Rust binary invoked via cedar_detect_client (bundled in cedar-solve)
#   cedar-solve  : Python fork of tetra3, imported as `tetra3`
#
# From migrating.rst in the cedar-solve repo, the canonical usage is:
#
#   from tetra3 import Tetra3, cedar_detect_client
#   t3 = Tetra3('t3_fov14_mag8')
#   cedar_detect = cedar_detect_client.CedarDetectClient()
#   centroids = cedar_detect.extract_centroids(image, sigma=8, max_size=10, use_binned=True)
#   solve_dict = t3.solve_from_centroids(centroids, fov_estimate=13.5)
#
# CedarDetectClient() starts and manages the cedar-detect-server binary itself;
# no subprocess or gRPC plumbing is needed in user code.
# The binary must be at <cedar-solve-install>/tetra3/bin/cedar-detect-server.
#
# solve_from_centroids returns a dict with keys:
#   'RA', 'Dec', 'Roll', 'FOV', 'matched_stars', ...
# or {} on failure.
#
# Coordinate convention: cedar-detect and cedar-solve both use (y, x)
# top-left-origin pixel coordinates — no centred-space conversion.
#
# Process layout (unchanged from tetra3rs_mp):
#   Process 0 (main):    spawns workers, monitors health.
#   Process 1 (camera):  picamera2 capture loop -> triple-buffered shared memory.
#   Process 2 (solver):  cedar-detect + cedar-solve -> shared Values + JSON.
#   Process 3 (lx200):   LX200/WiFi server — reads Values directly.

import os
import sys
import math
import socket
import time
import csv
import json
import ctypes
from datetime import datetime
from pathlib import Path
from multiprocessing import (
    Process, Queue, Event, Value,
    shared_memory, set_start_method
)

import numpy as np

# ---------------------------------------------------------------------------
# Shared constants
# ---------------------------------------------------------------------------
home_path   = str(Path.home())
version     = "6.6-cedar-mp-c2"
config_path = os.path.join(home_path, "Solver/eFinder.config")
solver_path = os.path.join(home_path, "Solver")

FRAME_H  = 760
FRAME_W  = 960
FRAME_SZ = FRAME_H * FRAME_W

CAM_ARCSEC_PX = 50.8
CAM_FOV_DEG   = 13.5

STATE_FILE = '/dev/shm/efinder_state.json'
LIVE_IMAGE = '/dev/shm/efinder_live.jpg'

# cedar-solve/detect coordinate convention: (y, x) top-left origin.
CENTRE_X = FRAME_W / 2.0
CENTRE_Y = FRAME_H / 2.0

N_FRAME_SLOTS = 3

# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------
def load_param():
    param = {}
    if os.path.exists(config_path):
        with open(config_path) as h:
            for line in h:
                parts = line.strip("\n").split(":")
                if len(parts) == 2:
                    param[parts[0]] = str(parts[1])
    return param

def save_param(param):
    with open(config_path, "w") as h:
        for key, value in param.items():
            h.write("%s:%s\n" % (key, value))

# ---------------------------------------------------------------------------
# Coordinate helpers
# ---------------------------------------------------------------------------
class Coordinates:
    def __init__(self):
        self._update_precession_constants()

    def _update_precession_constants(self):
        now  = datetime.now()
        decY = now.year + int(now.strftime('%j')) / 365.25
        self.t = decY - 2000
        self.T = self.t / 100
        self.m  = 3.07496 + 0.00186 * self.T
        self.n2 = 20.0431 - 0.0085  * self.T
        self.n1 = 1.33621 - 0.00057 * self.T

    def dateSet(self, timeOffset, timeStr, dateStr):
        days = 0
        sg   = float(timeOffset)
        hours, minutes, seconds = timeStr.split(':')
        hours = int(hours) + sg
        if hours >= 24:
            hours = str(int(hours - 24)); days = 1
        elif hours < 0:
            hours = str(int(hours + 24)); days = -1
        else:
            hours = str(int(hours))
        timeStr = hours + ':' + minutes + ':' + seconds
        month, day, year = dateStr.split('/')
        day = str(int(day) + days)
        dateStr = month + '/' + day + '/20' + year
        dt_str = dateStr + ' ' + timeStr
        print('Calculated UTC', dt_str)
        os.system('sudo date -u --set "%s"' % dt_str + '.000Z')
        self._update_precession_constants()

    def precess(self, r, d):
        dR = self.m + self.n1 * math.sin(math.radians(r)) * math.tan(math.radians(d))
        dD = self.n2 * math.cos(math.radians(r))
        return r + dR / 240 * self.t, d + dD / 3600 * self.t

    def hh2dms(self, dd):
        minutes, seconds = divmod(abs(dd) * 3600, 60)
        degrees, minutes = divmod(minutes, 60)
        return '%02d:%02d:%02d' % (degrees, minutes, seconds)

    def dd2aligndms(self, dd):
        sign = '+' if dd >= 0 else '-'
        minutes, seconds = divmod(abs(dd) * 3600, 60)
        degrees, minutes = divmod(minutes, 60)
        return '%s%02d*%02d:%02d' % (sign, degrees, minutes, seconds)

    def dd2dms(self, dd):
        sign = '+' if dd >= 0 else '-'
        minutes, seconds = divmod(abs(dd) * 3600, 60)
        degrees, minutes = divmod(minutes, 60)
        return '%s%02d:%02d:%02d' % (sign, degrees, minutes, seconds)

# ---------------------------------------------------------------------------
# Offset helpers — cedar-solve uses (y, x) top-left-origin pixel coordinates.
# d_x / d_y stored in eFinder.config as arcminutes on sky.
# solve_from_centroids accepts target_pixel=(y, x) to report RA/Dec for a
# pixel other than image centre.
# ---------------------------------------------------------------------------
def _offset_target_pixel(param):
    """Return (y, x) image pixel corresponding to configured sky offset."""
    dx_deg = float(param.get("d_x", "0")) / 60.0   # +X = right on sky
    dy_deg = float(param.get("d_y", "0")) / 60.0   # +Y = up on sky
    col = CENTRE_X + dx_deg * 3600.0 / CAM_ARCSEC_PX
    row = CENTRE_Y - dy_deg * 3600.0 / CAM_ARCSEC_PX
    return (row, col)

# ---------------------------------------------------------------------------
# CPU affinity helper
# ---------------------------------------------------------------------------
CPU_PINNING = {
    'main':   {0},
    'camera': {0},
    'lx200':  {1},
    'solver': {2, 3},
}

def _pin_cpu(label):
    cores = CPU_PINNING.get(label)
    if not cores:
        return
    try:
        avail = os.sched_getaffinity(0)
        want  = cores & avail
        if not want:
            print('[%s] cpu pin skipped: requested %s not in available %s' %
                  (label, sorted(cores), sorted(avail)))
            return
        os.sched_setaffinity(0, want)
        got = os.sched_getaffinity(0)
        if got == want:
            print('[%s] cpu pin: cores %s' % (label, sorted(got)))
        else:
            print('[%s] cpu pin partial: asked %s, got %s' %
                  (label, sorted(want), sorted(got)))
    except (OSError, AttributeError) as e:
        print('[%s] cpu pin not supported: %s' % (label, e))

# ===========================================================================
# PROCESS 1 - Camera  (unchanged from tetra3rs_mp)
# ===========================================================================
def camera_process(shm_names, frame_ready, cam_cmd_q, cam_result_q,
                   test_mode, latest_slot, frame_seq):
    _pin_cpu('camera')
    from picamera2 import Picamera2

    param = load_param()

    shms      = [shared_memory.SharedMemory(name=n) for n in shm_names]
    slot_bufs = [np.ndarray((FRAME_H, FRAME_W), dtype=np.uint8, buffer=s.buf)
                 for s in shms]

    picam2 = Picamera2()
    cfg = picam2.create_still_configuration(
        main={"size": (FRAME_W, FRAME_H), "format": "YUV420"},
        sensor={"output_size": (2028, 1520)},
        buffer_count=2,
    )
    picam2.configure(cfg)

    def _apply(exp_s, gain):
        picam2.stop()
        picam2.set_controls({
            "AeEnable":     False,
            "AwbEnable":    False,
            "ExposureTime": int(float(exp_s) * 1_000_000),
            "AnalogueGain": int(float(gain)),
        })
        picam2.start()

    _apply(param.get("Exposure", "0.2"), param.get("Gain", "20"))
    test_path = os.path.join(home_path, "Solver/test.npy")

    def _capture():
        if test_mode.value and os.path.exists(test_path):
            return np.load(test_path)
        arr = np.array(picam2.capture_array())
        return arr[0:FRAME_H, 0:FRAME_W]

    def _pick_write_slot(last_published):
        return (last_published + 1) % N_FRAME_SLOTS

    last_published = -1

    print('[camera] ready (triple-buffered)')

    while True:
        try:
            while True:
                cmd, a, b = cam_cmd_q.get_nowait()
                if cmd == 'set_exp':
                    _apply(a, b)
                elif cmd == 'capture_once':
                    cam_result_q.put(('frame', _capture().copy()))
                elif cmd == 'stop':
                    picam2.stop()
                    for s in shms: s.close()
                    return
        except Exception:
            pass

        write_idx = _pick_write_slot(last_published)
        arr = _capture()
        np.copyto(slot_bufs[write_idx], arr)

        with latest_slot.get_lock():
            latest_slot.value = write_idx
        with frame_seq.get_lock():
            frame_seq.value  += 1

        last_published = write_idx
        frame_ready.set()
        time.sleep(0.05)

# ===========================================================================
# PROCESS 2 - Solver (cedar-detect + cedar-solve)
# ===========================================================================
def solver_process(shm_names, frame_ready, cam_cmd_q, cam_result_q,
                   lx200_cmd_q, lx200_result_q,
                   shared_ra, shared_dec, offset_flag, test_mode,
                   latest_slot, frame_seq):
    """
    Uses cedar-solve (imported as `tetra3`) for plate solving and
    cedar_detect_client.CedarDetectClient() for centroid extraction.

    CedarDetectClient manages the cedar-detect-server binary itself —
    no subprocess or gRPC setup needed here.  The binary must be at:
        <cedar-solve-install>/tetra3/bin/cedar-detect-server

    extract_centroids() returns an Nx2 numpy array of (y, x) centroids,
    brightest first, in top-left-origin pixel coordinates — exactly what
    solve_from_centroids() expects.

    solve_from_centroids() returns a dict {'RA', 'Dec', 'Roll', 'FOV',
    'matched_stars', ...} or {} on failure.
    """
    _pin_cpu('solver')
    import serial as _pyserial
    from threading import Thread as _Thread, Lock as _Lock
    from queue import Queue as _ThreadQueue, Empty as _QueueEmpty, Full as _QueueFull
    from PIL import Image, ImageDraw, ImageFont, ImageEnhance, ImageOps

    param = load_param()
    coordinates = Coordinates()

    # Attach to triple-buffer slots
    shms      = [shared_memory.SharedMemory(name=n) for n in shm_names]
    slot_bufs = [np.ndarray((FRAME_H, FRAME_W), dtype=np.uint8, buffer=s.buf)
                 for s in shms]

    last_solved_seq = 0

    # ------------------------------------------------------------------
    # cedar-solve + cedar-detect setup
    # ------------------------------------------------------------------
    from tetra3 import Tetra3, cedar_detect_client

    # 't3_fov14_mag8' is a bare name; Tetra3 resolves it to
    # <tetra3_pkg>/data/t3_fov14_mag8.npz in the editable install.
    print('[solver] loading cedar-solve database...')
    t3 = Tetra3('t3_fov14_mag8')
    print('[solver] cedar-solve database loaded')

    # CedarDetectClient() locates and starts cedar-detect-server automatically.
    cedar_detect = cedar_detect_client.CedarDetectClient()
    print('[solver] cedar-detect client ready')

    try:
        fnt = ImageFont.truetype(os.path.join(home_path, "Solver/text.ttf"), 16)
    except Exception:
        fnt = ImageFont.load_default()

    # --- FOV calibration ---
    _fov_samples  = []
    _FOV_MIN, _FOV_MAX = 5, 20
    _fov_measured = float(param.get('fov_measured', '0'))

    def _get_fov():
        if _fov_measured > 0 and len(_fov_samples) >= _FOV_MIN:
            return _fov_measured, 0.3
        return CAM_FOV_DEG, 1.0

    def _update_fov(fov_deg):
        nonlocal _fov_measured
        if not fov_deg or fov_deg <= 0: return
        _fov_samples.append(fov_deg)
        if len(_fov_samples) > _FOV_MAX: _fov_samples.pop(0)
        if len(_fov_samples) < _FOV_MIN: return
        avg = sum(_fov_samples) / len(_fov_samples)
        if abs(avg - _fov_measured) > 0.05:
            _fov_measured = avg
            param['fov_measured'] = '%.4f' % avg
            save_param(param)
            print('[solver] FOV calibrated: %.3f deg' % avg)

    SOLVE_TIMEOUT_MS = 5000

    # solver state
    solve          = False
    solved_radec   = (0.0, 0.0)
    result_last    = {}
    centroids_last = None   # Nx2 (y, x) top-left
    stars          = '0'
    peak           = '0'
    eTime          = '00.00'
    keep           = False
    frame_n        = 0
    offset_str     = '%1.3f,%1.3f' % (0.0, 0.0)

    # --- background state-writer ---
    _state_snapshot = {}
    _state_lock     = _Lock()
    _state_dirty    = False

    def _snapshot_state():
        nonlocal _state_dirty
        try:
            with open('/sys/class/thermal/thermal_zone0/temp') as f:
                cpu_temp = int(f.read()) / 1000.0
        except Exception:
            cpu_temp = 0.0
        try:
            with open('/proc/self/status') as f:
                mem_kb = next(l for l in f if l.startswith('VmRSS:'))
            memory_mb = int(mem_kb.split()[1]) // 1024
        except Exception:
            memory_mb = 0
        s = {
            'ra':              solved_radec[0] / 15.0,
            'dec':             solved_radec[1],
            'solve_status':    'Solved' if solve else 'No solve',
            'solve_timestamp': int(time.time()),
            'stars':           stars,
            'peak':            peak,
            'exposure':        param.get('Exposure', '?'),
            'gain':            param.get('Gain', '?'),
            'solve_time':      eTime,
            'version':         version,
            'fov_measured':    round(_fov_measured, 3) if _fov_measured > 0 else None,
            'fov_samples':     len(_fov_samples),
            'offset_str':      offset_str,
            'cpu_temp':        round(cpu_temp, 1),
            'memory_usage':    memory_mb,
        }
        with _state_lock:
            _state_snapshot.clear()
            _state_snapshot.update(s)
            _state_dirty = True

    def _write_state():
        _snapshot_state()

    def _state_writer_thread():
        nonlocal _state_dirty
        while True:
            time.sleep(0.5)
            with _state_lock:
                if not _state_dirty:
                    continue
                snap = dict(_state_snapshot)
                _state_dirty = False
            try:
                tmp = STATE_FILE + '.tmp'
                with open(tmp, 'w') as f:
                    json.dump(snap, f)
                os.replace(tmp, STATE_FILE)
            except Exception as e:
                print('[solver] state flush failed:', e)

    _Thread(target=_state_writer_thread, daemon=True,
            name='eFinder-state').start()

    # --- background live-image writer ---
    _live_q = _ThreadQueue(maxsize=1)

    def _live_writer_thread():
        while True:
            try:
                arr, overlay = _live_q.get()
            except Exception:
                continue
            if arr is None:
                continue
            try:
                img  = Image.fromarray(arr)
                img2 = ImageEnhance.Contrast(img).enhance(5)
                img2 = img2.rotate(angle=180)
                if overlay is not None:
                    d = ImageDraw.Draw(img2)
                    d.text((5, 5), overlay, font=fnt, fill='white')
                tmp = LIVE_IMAGE + '.tmp'
                img2.save(tmp, format='JPEG')
                os.replace(tmp, LIVE_IMAGE)
            except Exception as e:
                print('[solver] live image write failed:', e)

    _Thread(target=_live_writer_thread, daemon=True,
            name='eFinder-live').start()

    def _write_live(arr):
        if solve and result_last:
            overlay = 'RA %s  Dec %s  Stars %s  %.2fs' % (
                coordinates.hh2dms(solved_radec[0] / 15),
                coordinates.dd2aligndms(solved_radec[1]),
                stars, float(eTime))
        else:
            overlay = None
        try:
            _live_q.put_nowait((arr, overlay))
        except _QueueFull:
            try:
                _live_q.get_nowait()
            except _QueueEmpty:
                pass
            try:
                _live_q.put_nowait((arr, overlay))
            except _QueueFull:
                pass

    def _save_debug(arr, txt):
        nonlocal frame_n, keep
        frame_n += 1
        img  = Image.fromarray(arr)
        img2 = ImageEnhance.Contrast(img).enhance(5)
        img2 = img2.rotate(angle=180)
        d    = ImageDraw.Draw(img2)
        d.text((70, 5), txt + "      Frame %d" % frame_n, font=fnt, fill='white')
        img2 = ImageOps.expand(img2, border=5, fill='red')
        img2.save(os.path.join(home_path, 'Solver/images/capture.jpg'))
        if frame_n > 1100:
            keep = False; frame_n = 0

    # ------------------------------------------------------------------
    # Centroid extraction via CedarDetectClient
    # ------------------------------------------------------------------
    def _extract_centroids(np_img):
        """Return (centroids_Nx2_yx, img_peak_int).

        cedar_detect.extract_centroids() takes a numpy uint8 array and
        returns an Nx2 float array [[y0,x0],[y1,x1],...] sorted brightest
        first, top-left origin.

        Parameters from migrating.rst example:
            sigma=8          — detection threshold in noise units
            max_size=10      — max star radius in pixels
            use_binned=True  — 2x2 binning before detection (faster, less noise)
        """
        try:
            centroids = cedar_detect.extract_centroids(
                np_img, sigma=8, max_size=10, use_binned=True)
            if centroids is None or len(centroids) == 0:
                return np.empty((0, 2), dtype=np.float32), 0
            # Peak: 5x5 window around brightest centroid
            row0 = int(round(centroids[0, 0]))
            col0 = int(round(centroids[0, 1]))
            r0, r1 = max(0, row0 - 2), min(FRAME_H, row0 + 3)
            c0, c1 = max(0, col0 - 2), min(FRAME_W, col0 + 3)
            img_peak = int(np_img[r0:r1, c0:c1].max()) if r1 > r0 and c1 > c0 else 0
            return centroids.astype(np.float32), img_peak
        except Exception as e:
            print('[solver] cedar-detect error:', e)
            return np.empty((0, 2), dtype=np.float32), 0

    # ------------------------------------------------------------------
    # Main solve function
    # ------------------------------------------------------------------
    def _do_solve(img):
        nonlocal solve, solved_radec, result_last, centroids_last
        nonlocal stars, peak, eTime

        t0 = time.time()
        np_img = img if img.dtype == np.uint8 else img.astype(np.uint8)

        centroids, img_peak = _extract_centroids(np_img)
        n_stars = len(centroids)

        print('[solver] centroids=%d  peak=%d' % (n_stars, img_peak))

        if n_stars < 15:
            solve = False
            if keep:
                _save_debug(img, "Bad image - %d stars  Exp=%ss Gain=%s" % (
                    n_stars, param.get('Exposure','?'), param.get('Gain','?')))
            return False

        stars = '%4d' % n_stars
        peak  = '%3d' % img_peak

        fov_est, fov_err = _get_fov()

        # Build target pixel from configured offset (image centre when d_x=d_y=0)
        tgt = _offset_target_pixel(param)
        use_target = (abs(tgt[0] - CENTRE_Y) > 0.5 or abs(tgt[1] - CENTRE_X) > 0.5)

        kwargs = dict(
            size          = (FRAME_H, FRAME_W),
            fov_estimate  = fov_est,
            fov_max_error = fov_err,
            solve_timeout = SOLVE_TIMEOUT_MS,
        )
        if use_target:
            kwargs['target_pixel'] = tgt

        result = t3.solve_from_centroids(centroids, **kwargs)

        eTime = ('%2.2f' % (time.time() - t0)).zfill(5)

        if not result or 'RA' not in result:
            solve = False
            if keep:
                _save_debug(img, "Not Solved - %s stars  Exp=%ss Gain=%s" % (
                    stars, param.get('Exposure','?'), param.get('Gain','?')))
            return False

        result_last    = result
        centroids_last = centroids

        ra_j2000  = result['RA']
        dec_j2000 = result['Dec']
        _update_fov(result.get('FOV'))

        if keep:
            _save_debug(img, "Peak=%d  Stars=%s  Exp=%ss Gain=%s" % (
                img_peak, stars, param.get('Exposure','?'), param.get('Gain','?')))

        # Precess J2000 -> JNow
        ra, dec = coordinates.precess(ra_j2000, dec_j2000)

        solved_radec = (ra, dec)
        solve = True

        shared_ra.value  = ra
        shared_dec.value = dec

        matched = result.get('matched_stars', '?')
        print('[solver] JNow', coordinates.hh2dms(ra / 15),
              coordinates.dd2aligndms(dec),
              '  matched=%s  solve=%.1fms' % (matched, float(eTime) * 1000))

        _snapshot_state()

        if _MOUNT_MODE != 'none':
            _Thread(target=_push_mount, args=(ra, dec), daemon=True).start()
        return True

    # --- mount push ---
    _MOUNT_MODE   = param.get('mount_mode', 'none').lower().strip()
    _MOUNT_HOST   = param.get('mount_host', '192.168.0.1').strip()
    _MOUNT_PORT   = int(param.get('mount_port', '9999'))
    _MOUNT_SERIAL = param.get('mount_serial', '/dev/ttyAMA0').strip()
    _MOUNT_BAUD   = int(param.get('mount_baud', '9600'))
    _mount_lock   = _Lock()
    _mount_sock   = None
    _mount_ser    = None

    def _fmt_ra(deg):
        h  = (deg / 15.0) % 24.0
        hh = int(h); mm = int((h - hh) * 60); ss = int(((h - hh) * 60 - mm) * 60)
        return '%02d:%02d:%02d' % (hh, mm, ss)

    def _fmt_dec(deg):
        s = '+' if deg >= 0 else '-'; d = abs(deg)
        dd = int(d); mm = int((d - dd) * 60); ss = int(((d - dd) * 60 - mm) * 60)
        return '%s%02d*%02d:%02d' % (s, dd, mm, ss)

    def _connect_mount():
        nonlocal _mount_sock, _mount_ser
        if _MOUNT_MODE == 'wifi':
            try:
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                s.settimeout(3.0); s.connect((_MOUNT_HOST, _MOUNT_PORT))
                s.settimeout(1.0)
                _mount_sock = s
                print('[solver] mount (WiFi) connected')
            except Exception as e:
                print('[solver] mount not reachable:', e)
        elif _MOUNT_MODE == 'serial':
            try:
                _mount_ser = _pyserial.Serial(
                    _MOUNT_SERIAL, _MOUNT_BAUD, timeout=1.0, write_timeout=1.0)
                print('[solver] mount (serial) connected')
            except Exception as e:
                print('[solver] mount serial not available:', e)

    def _push_mount(ra, dec):
        nonlocal _mount_sock, _mount_ser
        ra_s = _fmt_ra(ra); dec_s = _fmt_dec(dec)
        with _mount_lock:
            try:
                if _MOUNT_MODE == 'wifi' and _mount_sock:
                    for cmd in [':Sr%s#' % ra_s, ':Sd%s#' % dec_s, ':CM#']:
                        _mount_sock.sendall(cmd.encode('ascii'))
                        _mount_sock.recv(64)
                elif _MOUNT_MODE == 'serial' and _mount_ser:
                    for cmd in [':Sr%s#' % ra_s, ':Sd%s#' % dec_s, ':CM#']:
                        _mount_ser.reset_input_buffer()
                        _mount_ser.write(cmd.encode('ascii'))
                        time.sleep(0.1)
            except Exception as e:
                print('[solver] mount sync failed:', e)
                try:
                    if _mount_sock: _mount_sock.close()
                    if _mount_ser:  _mount_ser.close()
                except Exception: pass
                _mount_sock = _mount_ser = None
                _Thread(target=_connect_mount, daemon=True).start()

    if _MOUNT_MODE != 'none':
        _connect_mount()

    # --- camera helpers ---
    def _request_capture():
        cam_cmd_q.put(('capture_once', None, None))
        try:
            _, arr = cam_result_q.get(timeout=10.0)
            return arr
        except Exception:
            return slot_bufs[latest_slot.value].copy()

    def _set_camera(exp, gain):
        param['Exposure'] = str(exp)
        param['Gain']     = str(gain)
        save_param(param)
        cam_cmd_q.put(('set_exp', exp, gain))

    def _hip_lookup(hipId):
        name = sn = ''
        try:
            with open(os.path.join(home_path, 'Solver/starnames.csv')) as f:
                for row in csv.reader(f):
                    if str(row[1]) == str(hipId):
                        name = row[0].strip()
                        sn = (' (%s)' % row[2].strip()) if row[2].strip() else ''
                        break
        except Exception: pass
        return name, sn, str(hipId)

    # --- on-demand command handler ---
    def _handle(cmd, a, b):
        nonlocal solve, solved_radec, offset_str
        nonlocal keep, frame_n

        if cmd == 'adj_exp':
            new_exp = '%.1f' % max(0.1, float(param.get('Exposure','0.2'))
                                   + float(a) * 0.1)
            _set_camera(new_exp, param.get('Gain','20'))
            return new_exp

        elif cmd == 'adj_gain':
            g = max(0, min(50, float(param.get('Gain','20')) + float(a) * 5))
            _set_camera(param.get('Exposure','0.2'), '%.1f' % g)
            return '%.1f' % g

        elif cmd == 'select_exp':
            _set_camera(a, b); return '1'

        elif cmd == 'set_exp':
            _set_camera(float(a), param.get('Gain','20')); return '1'

        elif cmd == 'auto_exp':
            exp = float(param.get('Exposure','0.2'))
            _set_camera(exp, param.get('Gain','20'))
            img = _request_capture()
            for _ in range(20):
                cens, pk = _extract_centroids(img)
                nc = len(cens)
                print('[solver] auto_exp: %d stars %d peak' % (nc, pk))
                if nc < 20:
                    exp *= 2
                elif nc > 50 and pk > 250:
                    exp = int((exp / 2) * 10) / 10
                else:
                    break
                _set_camera(exp, param.get('Gain','20'))
                img = _request_capture()
            return str(exp)

        elif cmd == 'go_solve':
            img = _request_capture()
            return '1' if _do_solve(img) else '0'

        elif cmd == 'measure_offset':
            offset_flag.value = True
            ok = _do_solve(_request_capture())
            if not ok:
                offset_flag.value = False; return 'fail'

            pk = int(peak.strip()) if peak.strip().isdigit() else 0
            exp = float(param.get('Exposure','0.2'))
            while pk >= 220:   # >220 risks saturation on subsequent solve
                exp *= 0.75
                _set_camera(exp, param.get('Gain','20'))
                ok = _do_solve(_request_capture())
                pk = int(peak.strip()) if peak.strip().isdigit() else 0
            if not ok:
                offset_flag.value = False; return 'fail'

            if centroids_last is None or len(centroids_last) == 0:
                offset_flag.value = False; return 'fail'

            # centroids_last[0] is brightest: (y, x) top-left origin.
            cy_px, cx_px = centroids_last[0]
            dx_deg = (cx_px - CENTRE_X) * CAM_ARCSEC_PX / 3600.0
            dy_deg = (CENTRE_Y - cy_px) * CAM_ARCSEC_PX / 3600.0

            param['d_x'] = '{: .2f}'.format(dx_deg * 60)
            param['d_y'] = '{: .2f}'.format(dy_deg * 60)
            save_param(param)
            offset_str = '%1.3f,%1.3f' % (dx_deg, dy_deg)

            # Star name: cedar-solve may return matched_catalog_stars
            name, sn, hipId = '', '', '0'
            try:
                cat = result_last.get('matched_catalog_stars')
                if cat is not None and len(cat) > 0:
                    hipId = str(int(cat[0].get('hip_id', 0)))
                    name, sn, hipId = _hip_lookup(hipId)
            except Exception:
                pass

            offset_flag.value = False
            _write_state()
            return name + sn + ',HIP' + hipId + ',' + offset_str

        elif cmd == 'reset_offset':
            param['d_x'] = 0; param['d_y'] = 0
            offset_str = '%1.3f,%1.3f' % (0.0, 0.0)
            save_param(param); _write_state(); return '1'

        elif cmd == 'start_images':
            keep = (a == '1'); frame_n = 0 if keep else frame_n
            print('[solver] image saving:', 'on' if keep else 'off')
            return '1'

        elif cmd == 'date_set':
            coordinates.dateSet(*a); return '1'

        return 'ok'

    print('[solver] ready, entering main loop')

    while True:
        try:
            while True:
                cmd, a, b = lx200_cmd_q.get_nowait()
                result = _handle(cmd, a, b)
                lx200_result_q.put((cmd, result))
        except Exception:
            pass

        if frame_ready.wait(timeout=0.5) and not offset_flag.value:
            frame_ready.clear()

            seq_before  = frame_seq.value
            slot_before = latest_slot.value
            img = slot_bufs[slot_before].copy()
            seq_after   = frame_seq.value

            if seq_before == last_solved_seq:
                continue

            skipped = seq_before - last_solved_seq - 1
            if skipped > 0:
                print('[solver] skipped %d frame(s) (seq %d -> %d)' %
                      (skipped, last_solved_seq, seq_before))
            if seq_after != seq_before:
                print('[solver] frame seq changed during copy (%d -> %d) '
                      '— still safe via triple-buffer' % (seq_before, seq_after))

            last_solved_seq = seq_before

            _do_solve(img)
            _write_live(img)
            print('[solver] ****************')

# ===========================================================================
# PROCESS 3 - LX200 / WiFi server  (unchanged from tetra3rs_mp)
# ===========================================================================
def lx200_process(lx200_cmd_q, lx200_result_q,
                  shared_ra, shared_dec, offset_flag, test_mode):
    _pin_cpu('lx200')
    coordinates = Coordinates()

    def _read_state(key, default=''):
        try:
            with open(STATE_FILE) as f:
                return str(json.load(f).get(key, default))
        except Exception:
            return default

    def _cmd(cmd, a=None, b=None, timeout=15.0):
        lx200_cmd_q.put((cmd, a, b))
        try:
            _, result = lx200_result_q.get(timeout=timeout)
            return str(result)
        except Exception:
            return 'err'

    print('[lx200] starting on port 4060')
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    s.bind(('', 4060)); s.listen(50)
    raStr = decStr = ''
    timeOffset = '0'; timeStr = '23:00:00'

    while True:
        try:
            client, address = s.accept()
            try:
                client.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
            except Exception: pass
            print('[lx200] SkySafari connected from', address)
            while True:
                data = client.recv(1024)
                if not data: break
                pkt = data.decode('utf-8', 'ignore')
                time.sleep(0.001)

                ra  = shared_ra.value
                dec = shared_dec.value
                raPacket  = coordinates.hh2dms(ra / 15) + '#'
                decPacket = coordinates.dd2aligndms(dec) + '#'

                for x in pkt.split('#'):
                    if not x: continue
                    cmd = x[1:3]

                    if   x == ':GR':  client.send(raPacket.encode('ascii'))
                    elif x == ':GD':  client.send(decPacket.encode('ascii'))
                    elif cmd == 'St': client.send(b'1')
                    elif cmd == 'Sg': client.send(b'1')
                    elif cmd == 'SG':
                        if len(x) > 5:
                            client.send(b'1'); timeOffset = x[3:]
                        else:
                            res = _cmd('adj_gain', x[3:5])
                            client.send((':SG' + res + '#').encode('ascii'))
                    elif cmd == 'SL':
                        client.send(b'1'); timeStr = x[3:]
                    elif cmd == 'SC':
                        client.send(b'Updating Planetary Data#                              #')
                        _cmd('date_set', (timeOffset, timeStr, x[3:]))
                    elif cmd == 'RG': _cmd('select_exp', 0.1, 10)
                    elif cmd == 'RC': _cmd('select_exp', 0.1, 20)
                    elif cmd == 'RM': _cmd('select_exp', 0.2, 20)
                    elif cmd == 'RS': _cmd('select_exp', 0.5, 30)
                    elif cmd == 'Sr': raStr = x[3:]; client.send(b'1')
                    elif cmd == 'Sd': decStr = x[3:]; client.send(b'1')
                    elif cmd == 'MS': client.send(b'0')
                    elif cmd == 'Ms': _cmd('adj_exp', -1)
                    elif cmd == 'Mn': _cmd('adj_exp',  1)
                    elif cmd == 'Mw': _cmd('start_images', '0')
                    elif cmd == 'Me': _cmd('start_images', '1')
                    elif cmd == 'CM':
                        client.send(b'0')
                        _cmd('measure_offset', timeout=30.0)
                        try:
                            rp = raStr.split(':')
                            targetRa = int(rp[0]) + int(rp[1])/60 + int(rp[2])/3600
                            dp = decStr.split('*'); dd = dp[1].split(':')
                            targetDec = int(dp[0]) + math.copysign(
                                int(dd[0])/60 + int(dd[1])/3600, float(dp[0]))
                            print('[lx200] align target:', targetRa, targetDec)
                        except Exception: pass
                    elif x and x[-1] == 'Q':
                        _cmd('start_images', '0')
                    elif cmd == 'PS':
                        res = _cmd('go_solve')
                        client.send((':PS' + res + '#').encode('ascii'))
                    elif cmd == 'OF':
                        res = _cmd('measure_offset', timeout=30.0)
                        client.send((':OF' + res + '#').encode('ascii'))
                    elif cmd == 'GV':
                        client.send((':GV' + version + '#').encode('ascii'))
                    elif cmd == 'GO':
                        client.send((':GO' + _read_state('offset_str','0,0') + '#').encode('ascii'))
                    elif cmd == 'SO':
                        res = _cmd('reset_offset')
                        client.send((':SO' + res + '#').encode('ascii'))
                    elif cmd == 'GS':
                        client.send((':GS' + _read_state('stars','0') + '#').encode('ascii'))
                    elif cmd == 'GK':
                        client.send((':GK' + _read_state('peak','0') + '#').encode('ascii'))
                    elif cmd == 'Gt':
                        client.send((':Gt' + _read_state('solve_time','00.00') + '#').encode('ascii'))
                    elif cmd == 'SE':
                        res = _cmd('adj_exp', x[3:5])
                        client.send((':SE' + res + '#').encode('ascii'))
                    elif cmd == 'SX':
                        res = _cmd('set_exp', x.strip('#')[3:])
                        client.send((':SX' + res + '#').encode('ascii'))
                    elif cmd == 'GX':
                        res = _cmd('auto_exp', timeout=60.0)
                        client.send((':GX' + res + '#').encode('ascii'))
                    elif cmd == 'IM':
                        res = _cmd('start_images', x.strip('#')[3:4])
                        client.send((':IM' + res + '#').encode('ascii'))
                    elif cmd == 'TS':
                        test_mode.value = True
                        client.send(b':TS1#')
                    elif cmd == 'TO':
                        test_mode.value = False
                        client.send(b':TO1#')

            print('[lx200] SkySafari disconnected')
        except Exception as e:
            print('[lx200] server error:', e)
            try: s.close()
            except Exception: pass
            s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
            s.bind(('', 4060)); s.listen(50)

# ===========================================================================
# MAIN
# ===========================================================================
def main():
    if len(sys.argv) > 1:
        print('Killing running version')
        os.system('pkill -9 -f eFinder_cedar_mp.py')
        time.sleep(1)

    _pin_cpu('main')
    print('eFinder version', version)

    # Triple-buffered frame slots
    shms = [shared_memory.SharedMemory(create=True, size=FRAME_SZ)
            for _ in range(N_FRAME_SLOTS)]
    shm_names = [s.name for s in shms]
    print('Frame slots:', shm_names, '(%d bytes each)' % FRAME_SZ)

    latest_slot = Value(ctypes.c_int, 0)
    frame_seq   = Value(ctypes.c_uint64, 0)

    shared_ra   = Value(ctypes.c_double, 0.0)
    shared_dec  = Value(ctypes.c_double, 0.0)
    offset_flag = Value(ctypes.c_bool, False)
    test_mode   = Value(ctypes.c_bool, False)

    frame_ready    = Event()
    cam_cmd_q      = Queue()
    cam_result_q   = Queue()
    lx200_cmd_q    = Queue()
    lx200_result_q = Queue()

    proc_specs = {
        'camera': dict(
            target=camera_process,
            args=(shm_names, frame_ready, cam_cmd_q, cam_result_q,
                  test_mode, latest_slot, frame_seq)),
        'solver': dict(
            target=solver_process,
            args=(shm_names, frame_ready, cam_cmd_q, cam_result_q,
                  lx200_cmd_q, lx200_result_q,
                  shared_ra, shared_dec, offset_flag, test_mode,
                  latest_slot, frame_seq)),
        'lx200': dict(
            target=lx200_process,
            args=(lx200_cmd_q, lx200_result_q,
                  shared_ra, shared_dec, offset_flag, test_mode)),
    }

    procs = {}
    for name, spec in proc_specs.items():
        p = Process(target=spec['target'], args=spec['args'],
                    name='eFinder-' + name, daemon=True)
        p.start()
        procs[name] = p
        print('Started %s (pid %d)' % (name, p.pid))

    time.sleep(2.0)
    print('eFinder running — SkySafari -> port 4060')

    try:
        while True:
            time.sleep(30)
            for name, p in list(procs.items()):
                if not p.is_alive():
                    print('[main] %s died (exit %s) — restarting' % (name, p.exitcode))
                    spec = proc_specs[name]
                    new_p = Process(target=spec['target'], args=spec['args'],
                                    name='eFinder-' + name, daemon=True)
                    new_p.start()
                    procs[name] = new_p
                    print('[main] %s restarted (pid %d)' % (name, new_p.pid))
    except KeyboardInterrupt:
        print('eFinder stopped.')
    finally:
        for p in procs.values():
            p.terminate()
        for s in shms:
            try: s.close()
            except Exception: pass
            try: s.unlink()
            except Exception: pass

if __name__ == '__main__':
    set_start_method('fork')
    main()
