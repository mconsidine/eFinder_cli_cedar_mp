#!/usr/bin/python3

# eFinder_tetra3rs_mp — electronic finder scope, plate-solving over LX200/WiFi
# Derived from original work Copyright (C) 2025 Keith Venables (GPL v3)
# Simplified: direct picamera2, no Nexus, no GPIO, no LED, no WiFi switching
#
# Merged edition — combines the multiprocess architecture of the tiny_img
# branch with the tetra3rs Rust solver of the tetra3rs branch. Result:
#   - No cedar-detect gRPC server process
#   - No cedar-solve Python dependency
#   - In-process Rust centroid extraction + plate solve via tetra3rs
#   - Preserves tiny_img's mount-push, seeded-solve with blind fallback,
#     proc_specs-based restart pattern, and OTA update infrastructure.
#
#   Process 0 (main):    spawns worker processes; monitors health.
#   Process 1 (camera):  picamera2 capture loop -> shared memory frame slot.
#   Process 2 (solver):  tetra3rs extract + solve -> shared Values + JSON.
#   Process 3 (lx200):   LX200/WiFi server — reads Values directly, no IPC lag.
#
# Inter-process communication:
#   frame_shm      — SharedMemory  (760x960 uint8, camera -> solver)
#   frame_ready    — Event         (camera signals solver: new frame ready)
#   shared_ra      — Value(c_double) solver writes, lx200 reads — zero-copy
#   shared_dec     — Value(c_double) solver writes, lx200 reads — zero-copy
#   offset_flag    — Value(c_bool)   lx200/solver coordinate during offset meas.
#   test_mode      — Value(c_bool)   lx200 sets, camera reads
#   cmd_q          — Queue  lx200 -> solver: tuning/on-demand commands
#   result_q       — Queue  solver -> lx200: command results
#   cam_cmd_q      — Queue  solver -> camera: set_exp, capture_once
#   cam_result_q   — Queue  camera -> solver: captured frames
#
# /dev/shm/efinder_state.json  written by solver after each solve — slow
#   telemetry (stars, peak, eTime, etc.) read by lx200 only for diagnostic
#   commands (:GS, :GK, :Gt) which are not timing-critical.

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
version     = "6.6-tetra3rs-mp-tb8"
config_path = os.path.join(home_path, "Solver/eFinder.config")
solver_path = os.path.join(home_path, "Solver")

FRAME_H  = 760
FRAME_W  = 960
FRAME_SZ = FRAME_H * FRAME_W

CAM_ARCSEC_PX = 50.8
CAM_FOV_DEG   = 13.5

STATE_FILE = '/dev/shm/efinder_state.json'
LIVE_IMAGE = '/dev/shm/efinder_live.jpg'

# tetra3rs centroid convention: origin at image centre, +X right, +Y down.
# Offset stored in eFinder.config as arcminutes on sky; converted to centred
# pixel space for solver queries.
CENTRE_X = FRAME_W / 2.0
CENTRE_Y = FRAME_H / 2.0

# Triple-buffered camera -> solver frame handoff (items 1 + 3).
# Three SharedMemory slots; camera rotates writes across the two slots
# that the solver isn't currently reading. latest_slot and frame_seq are
# published atomically as the final step of each write.
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
# Offset / pixel conversion (tetra3rs centred convention)
# d_x / d_y stored in eFinder.config in arcminutes on sky.
# ---------------------------------------------------------------------------
def dxdy2centred(dx_deg, dy_deg):
    """Convert angular offset (degrees) to centred pixel coords."""
    cx =  dx_deg * 3600 / CAM_ARCSEC_PX      # +X right
    cy = -dy_deg * 3600 / CAM_ARCSEC_PX      # +Y down (dy up = negative cy)
    return cx, cy

def centred2dxdy(cx, cy):
    """Convert centred pixel coords to angular offset (degrees)."""
    dx_deg =  cx * CAM_ARCSEC_PX / 3600
    dy_deg = -cy * CAM_ARCSEC_PX / 3600
    return dx_deg, dy_deg

# ---------------------------------------------------------------------------
# CPU affinity helper (item 5)
# ---------------------------------------------------------------------------
# The Pi Zero 2W has 4x Cortex-A53. We pin each worker to specific cores
# to keep the solver's hot numeric loop from being interrupted by camera
# DMA completion handling or by SkySafari socket traffic, and to give the
# solver's two background threads (state writer, live JPEG renderer) room
# to run without stealing cycles from the main solve loop.
#
# Mapping (see CPU_PINNING dict). Tuneable if field measurements suggest
# a different assignment works better.
#
# sched_setaffinity can fail on restricted kernels, inside containers, or
# when systemd CPUAffinity= is already set; we log and carry on. Losing
# the tuning doesn't break correctness.

CPU_PINNING = {
    'main':   {0},      # supervisor — idle most of the time
    'camera': {0},      # light work, shares with USB/SDIO IRQs
    'lx200':  {1},      # isolated from solver so replies never queue
    'solver': {2, 3},   # heavy: main loop + state thread + live-jpeg thread
}

def _pin_cpu(label):
    """Best-effort CPU affinity pin. Logs result, never raises."""
    cores = CPU_PINNING.get(label)
    if not cores:
        return
    try:
        # Verify requested cores actually exist on this box. On non-Pi
        # hardware or odd kernels we just log and skip rather than pin
        # to a phantom core.
        avail = os.sched_getaffinity(0)
        want  = cores & avail
        if not want:
            print('[%s] cpu pin skipped: requested %s not in available %s' %
                  (label, sorted(cores), sorted(avail)))
            return
        os.sched_setaffinity(0, want)
        # Read back what actually stuck (cgroups can narrow the set).
        got = os.sched_getaffinity(0)
        if got == want:
            print('[%s] cpu pin: cores %s' % (label, sorted(got)))
        else:
            print('[%s] cpu pin partial: asked %s, got %s' %
                  (label, sorted(want), sorted(got)))
    except (OSError, AttributeError) as e:
        # AttributeError handles the unlikely case of a Python build
        # without sched_setaffinity (not Linux, not CPython on Linux, etc.)
        print('[%s] cpu pin not supported: %s' % (label, e))

# ===========================================================================
# PROCESS 1 - Camera
# ===========================================================================
def camera_process(shm_names, frame_ready, cam_cmd_q, cam_result_q,
                   test_mode, latest_slot, frame_seq):
    """
    Owns Picamera2. Rotates writes across three SharedMemory slots so the
    solver can read the most-recently-published slot without racing the
    writer. Publishes (latest_slot, frame_seq) atomically after each write.

    shm_names   : list of 3 SharedMemory names (created by main)
    frame_ready : Event — set after each completed publish
    latest_slot : Value(c_int)  — index (0/1/2) of last fully-written slot
    frame_seq   : Value(c_uint64) — monotonic frame counter

    Commands on cam_cmd_q:
        ('set_exp', exposure, gain)
        ('capture_once', None, None)  -> puts ('frame', array) on cam_result_q
        ('stop', None, None)
    """
    _pin_cpu('camera')
    from picamera2 import Picamera2

    param = load_param()

    # Attach to all three SHM slots; wrap each as a numpy view so np.copyto
    # writes straight into shared memory with no intermediate allocation.
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
        """
        Return a slot index that the solver is guaranteed not to be reading
        right now. With 3 slots and the solver always reading whichever slot
        was most recently published, we just pick any slot other than the
        published one. We prefer to rotate so we don't write to the same
        slot twice in a row — gives the solver a full extra cycle of
        read-safety headroom if it gets preempted mid-copy.
        """
        # Rotate: (published + 1) mod N. The solver may still be reading
        # `published` but not this next slot.
        return (last_published + 1) % N_FRAME_SLOTS

    last_published = -1   # sentinel: first pick goes to slot 0

    print('[camera] ready (triple-buffered)')

    while True:
        # drain on-demand commands
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

        # continuous capture into the next rotation slot
        write_idx = _pick_write_slot(last_published)
        arr = _capture()
        np.copyto(slot_bufs[write_idx], arr)

        # Publish: update latest_slot and bump frame_seq. The order matters:
        # incrementing the sequence *after* updating the slot index means
        # the solver's "read seq, copy, re-read seq" consistency check will
        # correctly detect a mid-copy write.
        with latest_slot.get_lock():
            latest_slot.value = write_idx
        with frame_seq.get_lock():
            frame_seq.value  += 1

        last_published = write_idx
        frame_ready.set()
        time.sleep(0.05)

# ===========================================================================
# PROCESS 2 - Solver (tetra3rs)
# ===========================================================================
def solver_process(shm_names, frame_ready, cam_cmd_q, cam_result_q,
                   lx200_cmd_q, lx200_result_q,
                   shared_ra, shared_dec, offset_flag, test_mode,
                   latest_slot, frame_seq):
    """
    Owns tetra3rs database.
    Continuous loop: wait for frame -> solve -> update shared_ra/shared_dec.
    Also handles on-demand commands from lx200_cmd_q.

    Reads frames from one of three SHM slots; `latest_slot` tells which
    slot is current, `frame_seq` is the monotonic frame number. We track
    `last_solved_seq` so we never redundantly solve the same frame.
    """
    _pin_cpu('solver')
    import tetra3rs
    import serial as _pyserial
    from threading import Thread as _Thread, Lock as _Lock
    from queue import Queue as _ThreadQueue, Empty as _QueueEmpty, Full as _QueueFull
    from PIL import Image, ImageDraw, ImageFont, ImageEnhance, ImageOps

    param = load_param()
    coordinates = Coordinates()

    # Attach to all three SHM slots — we pick which one to copy from at
    # each iteration based on latest_slot's value.
    shms      = [shared_memory.SharedMemory(name=n) for n in shm_names]
    slot_bufs = [np.ndarray((FRAME_H, FRAME_W), dtype=np.uint8, buffer=s.buf)
                 for s in shms]

    # Track last solved frame sequence so we don't re-solve the same frame
    # if the main loop wakes up more often than the camera publishes.
    last_solved_seq = 0

    # --- load database ---
    DB_PATH = os.path.join(home_path, 'Solver/efinder-tetra-database.bin')
    print('[solver] loading tetra3rs database...')
    db = tetra3rs.SolverDatabase.load_from_file(DB_PATH)
    print('[solver] tetra3rs ready  stars=%d  patterns=%d  max_fov=%.1f°' % (
        db.num_stars, db.num_patterns, db.max_fov_deg))

    # --- prewarm (item 9) ---
    # Force Rust-side scratch buffer allocation before the first real frame.
    # Synthesise a tiny list of centroids; we don't care if the solve itself
    # succeeds, only that the extraction/solve code paths have been exercised.
    try:
        _prewarm_img = np.zeros((FRAME_H, FRAME_W), dtype=np.uint8)
        # A handful of bright pixels so extract_centroids has something to find.
        for (y, x) in [(100, 120), (200, 300), (300, 500), (400, 200),
                       (500, 700), (600, 400), (250, 800), (350, 600)]:
            _prewarm_img[y-1:y+2, x-1:x+2] = 255
        _ext = tetra3rs.extract_centroids(_prewarm_img, sigma_threshold=10.0,
                                          max_centroids=20)
        try:
            db.solve_from_centroids(_ext.centroids,
                                    fov_estimate_deg=CAM_FOV_DEG,
                                    fov_max_error_deg=1.0,
                                    image_shape=(FRAME_H, FRAME_W),
                                    solve_timeout_ms=500)
        except Exception:
            pass   # a failed synthetic solve still allocates scratch
        del _prewarm_img, _ext
        print('[solver] prewarm complete')
    except Exception as e:
        print('[solver] prewarm skipped:', e)

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

    # --- offset in centred pixel space (tetra3rs convention) ---
    # d_x / d_y stored as arcminutes in config (on-sky).
    def _build_offset():
        dx_deg = float(param.get("d_x", "0")) / 60.0
        dy_deg = float(param.get("d_y", "0")) / 60.0
        return dxdy2centred(dx_deg, dy_deg)   # (cx, cy) centred pixels

    offset_cx, offset_cy = _build_offset()

    MAX_CENTROIDS          = 20
    SEEDED_TIMEOUT_MS      = 2000
    BLIND_TIMEOUT_MS       = 5000
    RESEED_TOLERANCE_DEG   = 4.0   # how close a new solve must be to trust the seed

    # --- Centroid peak-value attribute probe ---
    # Documented attribute on tetra3rs Centroid: `brightness` (integrated
    # intensity above background). Earlier experimental bindings exposed
    # it under other names; we probe a short list so the code works across
    # versions. Note: `brightness` / `mass` are INTEGRATED intensity, not
    # peak. `img_peak` in _do_solve uses the windowed np.max path instead,
    # which is the correct 0-255 peak. This attribute is only the fallback
    # when an image isn't passed to _centroid_peak.
    _PEAK_ATTR = None
    try:
        # Synthetic image with a few bright patches so extract_centroids
        # reliably returns at least one Centroid. Use 3x3 blocks rather
        # than single pixels — some centroid extractors require a local
        # neighborhood above threshold, not just one hot pixel.
        _probe_img = np.zeros((40, 40), dtype=np.uint8)
        for (_y, _x) in [(10, 10), (20, 30), (30, 15)]:
            _probe_img[_y-1:_y+2, _x-1:_x+2] = 255
        _probe_extraction = tetra3rs.extract_centroids(
            _probe_img, sigma_threshold=3.0, max_centroids=5,
        )
        if len(_probe_extraction.centroids) > 0:
            _c = _probe_extraction.centroids[0]
            for _name in ('brightness', 'mass', 'peak_val', 'peak_value', 'peak', 'intensity'):
                if hasattr(_c, _name) and getattr(_c, _name) is not None:
                    _PEAK_ATTR = _name
                    break
            print('[solver] centroid peak attribute: %r' % _PEAK_ATTR)
        else:
            print('[solver] centroid peak probe: no centroids from synthetic image')
    except Exception as _e:
        print('[solver] centroid peak probe failed:', _e)

    def _centroid_peak(centroid_list, np_img=None):
        """Peak intensity (0-255) of the brightest centroid.

        Primary path: read a 5x5 window around the brightest centroid's
        pixel location from np_img and take the max. This gives actual
        8-bit pixel intensity — what the log messages and auto_exp
        saturation check were designed around — without scanning the
        whole frame (25 pixels instead of 730K).

        Fallback: the cached centroid attribute (_PEAK_ATTR). Note that
        on tetra3rs 0.6.0 this resolves to `mass`, which is integrated
        intensity, not peak. Only useful as a relative brightness proxy
        when np_img isn't available.

        Returns 0 if the list is empty or neither path works.
        """
        if not centroid_list:
            return 0
        if np_img is not None:
            try:
                # centroid coords are centred (origin = image centre,
                # +X right, +Y down — tetra3rs convention). Convert to
                # top-left-origin indexing for the numpy array.
                cx = centroid_list[0].x
                cy = centroid_list[0].y
                row = int(round(CENTRE_Y + cy))
                col = int(round(CENTRE_X + cx))
                r0, r1 = max(0, row - 2), min(FRAME_H, row + 3)
                c0, c1 = max(0, col - 2), min(FRAME_W, col + 3)
                if r1 > r0 and c1 > c0:
                    return int(np_img[r0:r1, c0:c1].max())
            except Exception:
                pass   # fall through to attribute path
        if _PEAK_ATTR is not None:
            try:
                return int(getattr(centroid_list[0], _PEAK_ATTR))
            except Exception:
                pass
        return 0

    # solver state
    solve         = False
    solved_radec  = (0.0, 0.0)
    solution      = None
    centroids_last = []
    firstCentroid = None
    stars         = '0'
    peak          = '0'
    eTime         = '00.00'
    keep          = False
    frame_n       = 0
    offset_str    = '%1.3f,%1.3f' % (0.0, 0.0)

    # --- background state-writer (item 10) ---
    # Collect state into an in-memory dict on the hot path; flush it to
    # /dev/shm at 2 Hz from a daemon thread. Saves the JSON serialize +
    # write + fsync cost from the solve loop.
    _state_snapshot = {}
    _state_lock     = _Lock()
    _state_dirty   = False

    def _snapshot_state():
        """Hot-path: just stage the state dict, don't write."""
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
        """Back-compat wrapper kept for the measure_offset / reset_offset
        paths that expect a promptly-visible state file. Equivalent to
        the in-memory snapshot pattern; the background thread will flush
        within ~0.5 s, which is fine for those non-hot-path callers."""
        _snapshot_state()

    def _state_writer_thread():
        nonlocal _state_dirty
        last_flush = 0.0
        while True:
            time.sleep(0.5)   # 2 Hz flush cadence
            with _state_lock:
                if not _state_dirty:
                    continue
                snap = dict(_state_snapshot)
                _state_dirty = False
            try:
                # Atomic-ish: write to temp then rename. /dev/shm is tmpfs
                # so rename is cheap and avoids readers seeing half a file.
                tmp = STATE_FILE + '.tmp'
                with open(tmp, 'w') as f:
                    json.dump(snap, f)
                os.replace(tmp, STATE_FILE)
                last_flush = time.time()
            except Exception as e:
                print('[solver] state flush failed:', e)

    _Thread(target=_state_writer_thread, daemon=True,
            name='eFinder-state').start()

    # --- background live-image writer (item 11) ---
    # The live JPEG render is the single largest hidden cost in the solve
    # loop (Pillow contrast + rotate + JPEG encode on a 960x760 array is
    # easily 40-80ms on a Pi Zero 2W). Offload it to a thread with a
    # size-1 queue so the solver never blocks, and the renderer always
    # works on the newest available frame.
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
        """Hot-path: hand the frame to the render thread. Drops any backlog
        so the renderer always operates on the newest frame; if it falls
        behind we simply skip the older one rather than queuing it up."""
        if solve and solution is not None:
            overlay = 'RA %s  Dec %s  Stars %s  %.2fs' % (
                coordinates.hh2dms(solved_radec[0] / 15),
                coordinates.dd2aligndms(solved_radec[1]),
                stars, float(eTime))
        else:
            overlay = None
        # arr is already a copy (from the triple-buffer read in the main
        # loop), so passing it across the thread boundary is safe.
        try:
            _live_q.put_nowait((arr, overlay))
        except _QueueFull:
            # Drop the previous pending frame and push the new one.
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

    def _do_solve(img):
        nonlocal solve, solved_radec, solution, firstCentroid, centroids_last
        nonlocal stars, peak, eTime

        t0 = time.time()
        np_img = img if img.dtype == np.uint8 else img.astype(np.uint8)

        # tetra3rs centroid extraction — returns ExtractionResult with
        # centroids already in centred pixel space, brightest first.
        extraction = tetra3rs.extract_centroids(
            np_img,
            sigma_threshold=10.0,
            max_centroids=MAX_CENTROIDS,
        )
        centroid_list = extraction.centroids
        # Peak brightness is read as the max of a 5x5 window around the
        # brightest centroid — actual 8-bit pixel intensity, scoped to 25
        # pixels instead of a full-frame np.max() scan (~730K elements).
        img_peak = _centroid_peak(centroid_list, np_img)

        print('[solver] centroids=%d  peak=%d' % (len(centroid_list), img_peak))

        if len(centroid_list) < 15:
            solve = False
            if keep:
                _save_debug(img, "Bad image - %d stars  Exp=%ss Gain=%s" % (
                    len(centroid_list), param['Exposure'], param['Gain']))
            return False

        stars = '%4d' % len(centroid_list)
        peak  = '%3d' % img_peak

        fov_est, fov_err = _get_fov()

        # Tracking mode: when we have a recent successful solve, pass its
        # quaternion as `attitude_hint`. The solver then skips the expensive
        # 4-star pattern-hash search and does direct correspondence matching
        # against nearby catalog stars — typically 10-30 ms vs 100-500 ms
        # for lost-in-space. If the hint is stale (big slew, lost stars,
        # dome rotation lag), tetra3rs automatically falls back to
        # lost-in-space since strict_hint=False by default. No manual
        # retry needed.
        #
        # hint_uncertainty_deg sizes the cone search around the hint. For a
        # well-tracking mount between 200ms frames, actual motion is a few
        # arcsec; 1° is generous and keeps the cone small. Too tight and a
        # momentary gust of wind could miss; too wide and we lose the
        # speedup. 1° is a good starting point.
        kwargs = dict(
            fov_estimate_deg  = fov_est,
            fov_max_error_deg = fov_err,
            image_shape       = (FRAME_H, FRAME_W),
            solve_timeout_ms  = BLIND_TIMEOUT_MS,
        )
        if solve and solution is not None:
            try:
                kwargs['attitude_hint']        = solution.quaternion
                kwargs['hint_uncertainty_deg'] = 1.0
                kwargs['solve_timeout_ms']     = SEEDED_TIMEOUT_MS
            except Exception as _hint_e:
                # Older tetra3rs without tracking mode — ignore the hint and
                # fall back to lost-in-space. Logged once at startup would
                # be nicer than per-frame, but the cost is negligible.
                print('[solver] attitude_hint not supported:', _hint_e)
                kwargs.pop('attitude_hint', None)
                kwargs.pop('hint_uncertainty_deg', None)

        sol = db.solve_from_centroids(centroid_list, **kwargs)

        eTime = ('%2.2f' % (time.time() - t0)).zfill(5)

        if sol is None:
            solve = False
            if keep:
                _save_debug(img, "Not Solved - %s stars  Exp=%ss Gain=%s" % (
                    stars, param['Exposure'], param['Gain']))
            return False

        solution       = sol
        centroids_last = centroid_list
        firstCentroid  = centroid_list[0]   # brightest

        # sol.ra_deg / sol.dec_deg are the J2000 boresight.
        # pixel_to_world maps any pixel (centred coords) into J2000 RA/Dec —
        # this is the honest way to apply the offset, since tetra3rs does
        # the focal-plane projection (with parity_flip correction) internally.
        try:
            target = sol.pixel_to_world(offset_cx, offset_cy)
        except Exception:
            target = None
        if target is None:
            target_ra_deg, target_dec_deg = sol.ra_deg, sol.dec_deg
        else:
            target_ra_deg, target_dec_deg = target

        _update_fov(sol.fov_deg)

        # Precess J2000 -> JNow
        ra, dec = coordinates.precess(target_ra_deg, target_dec_deg)

        if keep:
            _save_debug(img, "Peak=%d  Stars=%s  Exp=%ss Gain=%s" % (
                img_peak, stars, param['Exposure'], param['Gain']))

        solved_radec = (ra, dec)
        solve = True

        # Zero-copy hot path — SkySafari :GR/:GD read these directly.
        shared_ra.value  = ra
        shared_dec.value = dec

        try:
            print('[solver] JNow', coordinates.hh2dms(ra / 15),
                  coordinates.dd2aligndms(dec),
                  '  solve=%.1fms' % sol.solve_time_ms)
        except Exception:
            print('[solver] JNow', coordinates.hh2dms(ra / 15),
                  coordinates.dd2aligndms(dec))
        _snapshot_state()   # hot path — in-memory only, bg thread flushes

        if _MOUNT_MODE != 'none':
            _Thread(target=_push_mount, args=(ra, dec), daemon=True).start()
        return True

    # --- mount push (ported from tiny_img) ---
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
            # Fallback: grab whichever slot camera published most recently.
            # Less fresh than an on-demand capture, but always available.
            return slot_bufs[latest_slot.value].copy()

    def _set_camera(exp, gain):
        param['Exposure'] = str(exp)
        param['Gain']     = str(gain)
        save_param(param)
        cam_cmd_q.put(('set_exp', exp, gain))

    # --- star name lookup ---
    def _lookup_star(catalog_id):
        hipId = str(abs(int(catalog_id)))   # negative IDs = Hipparcos gap-fill
        name = sn = ''
        try:
            with open(os.path.join(home_path, 'Solver/starnames.csv')) as f:
                for row in csv.reader(f):
                    if str(row[1]) == hipId:
                        name = row[0].strip()
                        sn = (' (%s)' % row[2].strip()) if row[2].strip() else ''
                        break
        except Exception: pass
        return name, sn, hipId

    # --- on-demand command handler ---
    def _handle(cmd, a, b):
        nonlocal solve, solved_radec, offset_cx, offset_cy, offset_str
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
                ext = tetra3rs.extract_centroids(img, sigma_threshold=10.0,
                                                 max_centroids=MAX_CENTROIDS)
                nc  = len(ext.centroids)
                pk  = _centroid_peak(ext.centroids, img)
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
            exp = float(param.get('Exposure','0.2'))
            while float(peak) > 255:
                exp *= 0.75
                _set_camera(exp, param.get('Gain','20'))
                ok = _do_solve(_request_capture())
            if not ok:
                offset_flag.value = False; return 'fail'

            # firstCentroid.x / .y are centred pixel coordinates.
            cx = firstCentroid.x
            cy = firstCentroid.y
            offset_cx = cx
            offset_cy = cy
            dx_deg, dy_deg = centred2dxdy(cx, cy)
            param['d_x'] = '{: .2f}'.format(float(60 * dx_deg))
            param['d_y'] = '{: .2f}'.format(float(60 * dy_deg))
            save_param(param)
            offset_str = '%1.3f,%1.3f' % (dx_deg, dy_deg)

            cat_id = solution.matched_catalog_ids[0] if \
                len(solution.matched_catalog_ids) > 0 else 0
            name, sn, hipId = _lookup_star(cat_id)

            offset_flag.value = False
            _write_state()
            return name + sn + ',HIP' + hipId + ',' + offset_str

        elif cmd == 'reset_offset':
            offset_cx = 0.0
            offset_cy = 0.0
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
        # drain on-demand commands first
        try:
            while True:
                cmd, a, b = lx200_cmd_q.get_nowait()
                result = _handle(cmd, a, b)
                lx200_result_q.put((cmd, result))
        except Exception:
            pass

        # wait for a new frame (up to 0.5 s)
        if frame_ready.wait(timeout=0.5) and not offset_flag.value:
            frame_ready.clear()

            # Triple-buffered read with sequence consistency check.
            # 1. Snapshot (slot, seq) — these two may be updated by camera
            #    between reads, but we only care that the SEQ we record
            #    belongs to the SLOT we actually copied from.
            # 2. Copy that slot's contents.
            # 3. Re-read seq. If it changed, camera published a new frame
            #    during our copy — we may have read a mix of old+new bytes.
            #    With 3 slots and camera rotating, the slot we copied from
            #    is guaranteed not to be the one camera just wrote, so the
            #    copy is actually consistent; the only thing we missed is
            #    freshness. Accept the copy but log the skip.
            seq_before  = frame_seq.value
            slot_before = latest_slot.value
            img = slot_bufs[slot_before].copy()
            seq_after   = frame_seq.value

            if seq_before == last_solved_seq:
                # Same frame we already solved; woke up spuriously or the
                # solver outran the camera. Skip without logging noise.
                continue

            skipped = seq_before - last_solved_seq - 1
            if skipped > 0:
                # Camera published more frames than we could solve. This is
                # normal when a solve takes longer than one exposure — log
                # at low volume so it shows up in profiling but doesn't
                # spam the journal during steady-state operation.
                print('[solver] skipped %d frame(s) (seq %d -> %d)' %
                      (skipped, last_solved_seq, seq_before))
            if seq_after != seq_before:
                print('[solver] frame seq changed during copy (%d -> %d) '
                      '— copy still safe via triple-buffer' %
                      (seq_before, seq_after))

            last_solved_seq = seq_before

            _do_solve(img)
            _write_live(img)
            print('[solver] ****************')

# ===========================================================================
# PROCESS 3 - LX200 / WiFi server
# ===========================================================================
def lx200_process(lx200_cmd_q, lx200_result_q,
                  shared_ra, shared_dec, offset_flag, test_mode):
    """
    Serves SkySafari on port 4060.
    Reads ra/dec directly from shared Values - no IPC latency on hot path.
    Reads other telemetry from /dev/shm/efinder_state.json (non-critical path).
    Sends tuning/on-demand commands to solver via lx200_cmd_q.
    """
    _pin_cpu('lx200')
    coordinates = Coordinates()

    def _read_state(key, default=''):
        """Read a single key from the JSON state file — non-critical path only."""
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
            # TCP-level tuning for SkySafari's many small :GR/:GD polls:
            # disable Nagle so each reply flushes immediately rather than
            # waiting 40 ms for a coalescing partner.
            try:
                client.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
            except Exception: pass
            print('[lx200] SkySafari connected from', address)
            while True:
                data = client.recv(1024)
                if not data: break
                pkt = data.decode('utf-8', 'ignore')
                # Item 7: previously time.sleep(0.02) — 20 ms of deliberate
                # latency on every SkySafari poll. Reduced to 1 ms as a
                # conservative first landing; can be removed entirely once
                # field-tested. TCP_NODELAY (set above) already handles
                # Nagle coalescing.
                time.sleep(0.001)

                # Hot path: read directly from shared memory — no IPC.
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
        os.system('pkill -9 -f eFinder_tetra3rs_mp.py')
        time.sleep(1)

    _pin_cpu('main')
    print('eFinder version', version)

    # Triple-buffered frame slots. Each is a separate SharedMemory region
    # of FRAME_SZ bytes; camera rotates writes across them, solver reads
    # whichever one latest_slot indicates.
    shms = [shared_memory.SharedMemory(create=True, size=FRAME_SZ)
            for _ in range(N_FRAME_SLOTS)]
    shm_names = [s.name for s in shms]
    print('Frame slots:', shm_names, '(%d bytes each)' % FRAME_SZ)

    # Atomic handoff Values (items 1 + 3).
    # latest_slot: which slot was most recently fully written.
    # frame_seq:   monotonic frame counter; solver uses it to detect
    #              skipped frames and to avoid re-solving the same frame.
    latest_slot = Value(ctypes.c_int, 0)
    frame_seq   = Value(ctypes.c_uint64, 0)

    # shared Values for hot-path ra/dec — direct memory read, zero IPC
    shared_ra   = Value(ctypes.c_double, 0.0)
    shared_dec  = Value(ctypes.c_double, 0.0)
    offset_flag = Value(ctypes.c_bool, False)
    test_mode   = Value(ctypes.c_bool, False)

    # queues and events
    frame_ready    = Event()
    cam_cmd_q      = Queue()
    cam_result_q   = Queue()
    lx200_cmd_q    = Queue()
    lx200_result_q = Queue()

    # Store target/args alongside each Process so the restart logic can
    # reconstruct them without relying on private _target/_args attributes.
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
