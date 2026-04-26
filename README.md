# eFinder — tetra3rs

this needs to be rewritten 

Plate-solving electronic finder scope for the Raspberry Pi Zero 2W.  
Connects to **SkySafari** over WiFi using the LX200 protocol on TCP port 4060.  
Uses an **IMX477 camera** (Raspberry Pi HQ Camera or Arducam IMX477) and the **Tetra3rs** star-pattern matching library for fast, reliable plate solving.

---

## Hardware Required

- Raspberry Pi Zero 2W
- Raspberry Pi HQ Camera or Arducam IMX477 board
- MicroSD card (16 GB or larger, Class 10 / A1 recommended)
- Micro-USB cable (data-capable, not charge-only) — for the USB data port
- Power supply via the PWR micro-USB port
- Optional: ADXL343 accelerometer on I2C for altitude readout

---

## Part 1 — Burn the SD Card

### 0.1 Download OS image

Under releases there should be an archive of a pre-built Pi
OS.  Download it and if necessary extract it from its zip 
file so that a .xz version is available

### 1.1 Download Raspberry Pi Imager

Download and install **Raspberry Pi Imager** from:  
https://www.raspberrypi.com/software/

### 1.2 Choose OS

- Click **Choose OS**
- Select **Use Custom**
- Select the local .xz file extracted in step 0.1 

### 1.3 Choose Storage

Insert your MicroSD card and select it.

### 1.4 DO NOT Configure OS Settings

The Pi needs internet access during installation to download packages. After installation you will access it and
switch it to STA (station) mode so that it has internet access.  Afterwards, it will switch to AP (hotspot) mode automatically. 
Your home WiFi credentials are only used once.

### 1.5 Write the Image

Click **Write** and wait for the process to complete. Eject the card safely when done.

---

## Part 2 — First Boot without a data tether

### 2.1 Insert the Card and Power On

Insert the MicroSD card into the Pi Zero 2W.  
Connect power to the **PWR** micro-USB port (the outer port, labelled PWR IN).  
Wait approximately 60–90 seconds for first boot to complete.

### 2.2 Find the Pi on Your Network

**macOS / Linux** — open Terminal and connect:

```bash
ssh efinder@efinder.local
```

If that does not resolve, find the IP address from your router's device list and connect directly:

```bash
ssh efinder@<ip-address>
```

**Windows** — SSH is built into Windows 10 and 11. Open **PowerShell** or **Command Prompt** and run the same command:

```powershell
ssh efinder@efinder.local
```

If `efinder.local` does not resolve, find the IP from your router and use it directly. Alternatively install **PuTTY** (https://www.putty.org) for a graphical SSH client — enter hostname `efinder.local` or the IP address, port `22`, and connect.

Log in with the password 12345678.

## OR

## Part 2 — First Boot with a data tether

### 2.1 Insert the Card and Power On

Insert the MicroSD card into the Pi Zero 2W.  
Connect power to the **DATA/PWr** micro-USB port (the inner port).  
Wait approximately 60–90 seconds for first boot to complete.

The use a serial terminal program or 'screen' (on linux) to access
the Pi at something like /dev/ttyUSB0 or /dev/ttyACM0 (linux)

Login with username efinder and the password 12345678

### 
---

## Part 3 — Install eFinder

All commands below are run on the Pi over SSH.

### 3.1 Put device into Station mode for internet access

```bash
sudo bash ~/station.sh
```

### 3.2 Run the Installer

```bash
sudo bash ~/install.sh
```

The installer will work through a number of steps automatically, among
which are:

| Step | What happens |
|------|-------------|
| a | Working directories created, support files deployed |
| b | Tetra3 plate-solving library installed and star databases downloaded |
| c | Samba file share configured |
| d | Apache/PHP web server configured for OTA updates |
| e | `/boot/firmware/config.txt` updated for IMX477 camera and USB serial gadget |
| f | SSH, I2C, USB serial console enabled |
| g | `eFinder.py` deployed |
| h | systemd services installed for auto-start on boot |

> **Duration:** The Tetra3 database download is the longest step — allow 10–20 minutes depending on your connection speed. Do not interrupt the installer once it has started.

### 3.3 Reboot

When the installer finishes it will print a summary and ask:

```
Reboot now? [y/N]
```

Type `y` and press Enter. The Pi will reboot.

---

## Part 4 — After Reboot

After rebooting the Pi will **no longer connect to your home WiFi**. It now runs as a WiFi access point.

### 4.1 Connect to the eFinder WiFi Network

On your phone, tablet, or laptop:

- Open WiFi settings
- Connect to the network named **`efinder`** followed by 4 hex characters, e.g. `efinder3a2f`
- Password: **`12345678`**

The SSID and password are also saved on the Pi at:  
`/home/efinder/Solver/default_hotspot.txt`

### 4.2 Verify the eFinder is Running

From a device connected to the eFinder WiFi, open a terminal and SSH in:

```bash
ssh efinder@192.168.50.1
```
or
```bash
ssh efinder@efinder.local
```

Check the service status:

```bash
sudo systemctl status efinder
```

You should see `active (running)`. To watch the live log:

```bash
journalctl -u efinder -f
```

A healthy startup has logging info like:

```
eFinder version 6.6
Coordinates ready
No accelerometer fitted
Loading Tetra3 database…
Tetra3 ready
Offset: (380.0, 480.0)
Starting solve loop…
Starting WiFi/LX200 server…
eFinder running — waiting for SkySafari connection on port 4060
```

---

## Part 5 — Connect SkySafari

### 5.1 Configure SkySafari

In SkySafari, go to **Settings → Telescope → Setup**:

| Setting | Value |
|---------|-------|
| Scope Type | Meade LX-200 GPS |
| Mount Type | Alt-Az (or match your mount) |
| IP Address | `192.168.50.1` or `efinder.local`|
| Port | `4060` |

### 5.2 Connect

Tap **Connect** in SkySafari. The eFinder log will show:

```
SkySafari connected from ('192.168.50.1', <port>)
```

SkySafari will send its date and time to the Pi on connection. The Pi sets its system clock from this — no internet or NTP is required in the field.

---

## Part 6 — Normal Use

Power on the Pi. Wait 30–60 seconds. Connect your phone to the eFinder WiFi. Connect SkySafari. The eFinder plate-solves continuously and reports the current position to SkySafari automatically.

No interaction with the Pi is needed during normal use.

---

## Part 7 — Focusing

TBD

### Focus Targets

**Best:** Stars at night. Point at a bright star (magnitude 1–3), use `--shutter 50000 --gain 10` as a starting point to avoid saturation, and adjust until the star disc is a tight point.

**Good:** The Moon. Bright, high-contrast, genuinely at infinity. Adjust `--shutter` down to `1000`–`5000` to avoid overexposure.

**Acceptable for bench setup:** An artificial star — a pinhole (0.1–0.5 mm) over a light source at 10 metres or more. Closer than 10 m will not be at true infinity focus for a 25 mm focal length lens.

**Avoid:** Daytime terrestrial targets at less than ~500 m. They are close enough that infinity focus will be slightly off, and the eFinder operates exclusively on star fields.

### Verifying Focus with the App

Once focused, restart the app and check the solve quality in the log:

```bash
sudo systemctl start efinder
journalctl -u efinder -f
```

A well-focused image on a clear night should show **30–80 centroids** and a **peak pixel value** in the range 150–230. Too few centroids means underexposed or out of focus. Peak value consistently above 240 means overexposed — reduce `Exposure` in `eFinder.config`.

---

## Maintenance Access

### SSH over WiFi (primary)

When connected to the eFinder WiFi network, open Terminal (macOS/Linux) or PowerShell/Command Prompt (Windows):

```bash
ssh efinder@192.168.50.1
```
or
```bash
ssh efinder@efinder.local
```

Windows users can also use **PuTTY** — enter IP `192.168.50.1`, port `22`.

### USB Serial Console (recovery / no WiFi)

Connect a data-capable USB cable to the **USB** port (the inner port, not PWR).

**macOS** — the device appears as `/dev/cu.usbmodem*`. Use **CoolTerm** or screen:

```bash
screen /dev/cu.usbmodem* 115200
```

**Linux** — the device appears as `/dev/ttyACM0`. Use screen or minicom:

```bash
screen /dev/ttyACM0 115200
```

**Windows** — the device appears as a **COM port** (e.g. `COM3`) in Device Manager under *Ports (COM & LPT)*. Use **CoolTerm** (https://freeware.the-meiers.org) or **PuTTY** — select Serial, enter the COM port number, speed 115200.

All platforms: connect at **115200 baud, 8N1**. This gives a login prompt directly — useful if the WiFi AP is not starting correctly.

### Samba File Share

**macOS** — in Finder press **Cmd+K** and enter:

```
smb://192.168.50.1/efindershare
```

**Linux** — open a file manager and connect to:

```
smb://192.168.50.1/efindershare
```

Or mount from the terminal:

```bash
sudo mount -t cifs //192.168.50.1/efindershare /mnt/efinder -o username=efinder
```

**Windows** — open **File Explorer**, click in the address bar and enter:

```
\\192.168.50.1\efindershare
```

Or map it as a network drive: right-click **This PC → Map network drive**, enter `\\192.168.50.1\efindershare`, and tick *Connect using different credentials*.

All platforms: Username `efinder` / Password `efinder`.

This gives read/write access to `/home/efinder` — useful for copying config files or captured images without SSH.

---

## OTA Updates

TBD

---

## Troubleshooting

**eFinder service fails to start**

```bash
journalctl -u efinder -b --no-pager
```

Common causes:
- Camera not detected — check `libcamera-hello --list-cameras`; verify the ribbon cable is seated; confirm `dtoverlay=imx477` is in `/boot/firmware/config.txt`
- Tetra3 database missing — check `~/venv-efinder/lib/python*/site-packages/tetra3/data/` contains `.npz` files

**SkySafari cannot connect**

- Confirm your device is on the eFinder WiFi, not your home network
- Confirm port is set to `4060` not `4061`
- Check `sudo systemctl status efinder` — the service must be running before SkySafari connects

**WiFi AP not appearing after reboot**

```bash
nmcli connection show
nmcli connection up efinder-ap
```

**Restore home WiFi temporarily** (to re-run installer or download updates):

```bash
sudo nmcli connection modify preconfigured autoconnect yes
sudo nmcli connection up preconfigured
```

To return to AP-only mode afterward:

```bash
sudo nmcli connection modify preconfigured autoconnect no
sudo reboot
```

**Check camera overlay**

```bash
vcgencmd get_config dtoverlay
libcamera-hello --list-cameras
```

---

## Diagnostic Commands over TCP

Any TCP client (e.g. `nc`, a Python script, or a custom app) connected to port 4060 can send these diagnostic commands while SkySafari is not connected:

| Command | Response | Description |
|---------|----------|-------------|
| `:PS#` | `:PS1#` or `:PS0#` | Trigger a plate solve |
| `:GV#` | `:GV6.6#` | Get software version |
| `:GS#` | `:GS  42#` | Get star count from last solve |
| `:GK#` | `:GK 187#` | Get peak pixel value from last solve |
| `:Gt#` | `:Gt02.34#` | Get elapsed solve time (seconds) |
| `:GO#` | `:GO0.012,0.003#` | Get current pointing offset |
| `:SO#` | `:SO1#` | Reset offset to image centre |
| `:OF#` | `:OFAlbireo,HIP98110,...#` | Measure offset, return alignment star |
| `:GX#` | `:GX0.2#` | Auto-expose, return chosen exposure |
| `:SX0.3#` | `:SX1#` | Set exposure to 0.3 seconds |
| `:GA#` | `:GA45#` | Get scope altitude (requires accelerometer) |
| `:TS#` | `:TS1#` | Enable test mode (uses test.npy instead of camera) |
| `:TO#` | `:TO1#` | Disable test mode |
| `:IM1#` | `:IM1#` | Start saving debug images to `Solver/images/` |
| `:IM0#` | `:IM1#` | Stop saving debug images |

Example using netcat (macOS/Linux):

```bash
echo -n ':GV#' | nc 192.168.50.1 4060
```

**Windows** — use PowerShell's built-in TCP client:

```powershell
$tcp = New-Object System.Net.Sockets.TcpClient('192.168.50.1', 4060)
$stream = $tcp.GetStream()
$bytes = [System.Text.Encoding]::ASCII.GetBytes(':GV#')
$stream.Write($bytes, 0, $bytes.Length)
Start-Sleep -Milliseconds 200
$buf = New-Object byte[] 1024
$len = $stream.Read($buf, 0, 1024)
[System.Text.Encoding]::ASCII.GetString($buf, 0, $len)
$tcp.Close()
```

Or install **netcat for Windows** via winget: `winget install netcat` and use the same syntax as macOS/Linux.

---

## File Layout on the Pi

```
/home/efinder/
├── venv-efinder/          Python virtual environment
├── tetra3/                Tetra3 source
├── eFinder_cli/           Repository clone (tinySS branch)
├── uploads/               OTA update zip drop location
└── Solver/
    ├── eFinder.py         Main application
    ├── eFinder.config     Saved exposure, gain, and offset settings
    ├── starnames.csv      HIP star name lookup for offset measurement
    ├── text.ttf           Font for debug image annotation
    ├── test.npy           Test image for test mode (northern hemisphere)
    ├── images/            Debug capture output (tmpfs — cleared on reboot)
    ├── default_hotspot.txt  AP SSID and password
    └── www/               Web UI files for OTA updater
```

---

## License

Derived from original eFinder work Copyright (C) 2025 Keith Venables, licensed under the GNU General Public License v3.  
See https://github.com/AstroKeith/eFinder_cli for the original project.
