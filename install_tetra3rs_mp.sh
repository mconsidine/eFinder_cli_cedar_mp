#!/bin/bash
# =============================================================================
# eFinder tetra3rs_mp — post-flash LAN connection script
#
# Everything except connecting to a local area network has been moved into
# the CI image build process. Flash the image produced by CI, power on,
# and the eFinder AP (efinder<last4> / 192.168.50.1) is ready without
# running this script.
#
# Run this script ONLY if you want to connect the Pi to an existing WiFi
# network (station mode) in addition to the AP, e.g. for OTA updates or
# to reach it via your home LAN.
#
# Usage:
#   Interactive : sudo bash install_tetra3rs_mp.sh
#   Automated   : WIFI_SSID=myssid WIFI_PASS=mypassword \
#                 bash install_tetra3rs_mp.sh --non-interactive
# =============================================================================
set -eo pipefail

EFINDER_USER=efinder
NON_INTERACTIVE=false
[[ "$1" == "--non-interactive" ]] && NON_INTERACTIVE=true

# ---------------------------------------------------------------------------
# Guard: must run as efinder user (with sudo)
# ---------------------------------------------------------------------------
if [ "$NON_INTERACTIVE" = false ]; then
    RUNNING_AS="${SUDO_USER:-$(logname 2>/dev/null)}"
    if [ "$RUNNING_AS" != "$EFINDER_USER" ]; then
        echo "ERROR: Run as '$EFINDER_USER': sudo bash install_tetra3rs_mp.sh"
        exit 1
    fi
fi

# ---------------------------------------------------------------------------
# Guard: hardware model
# ---------------------------------------------------------------------------
PI_MODEL=$(tr -d '\0' < /proc/device-tree/model 2>/dev/null || echo "unknown")
if [[ "$PI_MODEL" != *"Zero 2"* ]]; then
    echo "WARNING: Designed for Raspberry Pi Zero 2W. Detected: $PI_MODEL"
    if [ "$NON_INTERACTIVE" = false ]; then
        read -rp "Continue anyway? [y/N] " confirm
        [[ "$confirm" =~ ^[Yy]$ ]] || exit 1
    fi
fi

echo "============================================================================="
echo " eFinder LAN connection setup"
echo " Device : $PI_MODEL"
echo " Mode   : $( [ "$NON_INTERACTIVE" = true ] && echo non-interactive || echo interactive )"
echo "============================================================================="
echo ""
echo " NOTE: The eFinder AP (efinder<last4>, 192.168.50.1) is already active."
echo " This script adds a station-mode connection to a local WiFi network."
echo " Run ~/ap.sh at any time to switch back to AP-only mode."
echo "============================================================================="
echo ""

# ---------------------------------------------------------------------------
# Collect credentials
# ---------------------------------------------------------------------------
if [ "$NON_INTERACTIVE" = true ]; then
    : "${WIFI_SSID:?ERROR: WIFI_SSID must be exported}"
    : "${WIFI_PASS:?ERROR: WIFI_PASS must be exported}"
else
    # Show available networks to help the user pick the right SSID
    echo "Scanning for available networks..."
    nmcli device wifi rescan 2>/dev/null || true
    sleep 2
    nmcli -f SSID,SIGNAL,SECURITY device wifi list --rescan no 2>/dev/null \
        | head -20 || true
    echo ""
    read -rp  "Enter your WiFi SSID: " WIFI_SSID
    read -rsp "Enter your WiFi password (blank=open network): " WIFI_PASS; echo ""
fi

[ -z "$WIFI_SSID" ] && { echo "ERROR: SSID cannot be empty."; exit 1; }

# ---------------------------------------------------------------------------
# Determine security type from scan
# ---------------------------------------------------------------------------
SEC=$(nmcli -t -f SSID,SECURITY device wifi list --rescan no \
      | awk -F: -v s="$WIFI_SSID" '$1==s {print $2; exit}')

if [ -z "$SEC" ]; then
    echo ""
    echo "WARNING: '$WIFI_SSID' not found in current scan."
    echo "Visible networks:"
    nmcli -f SSID,SIGNAL,SECURITY device wifi list --rescan no | head -20
    if [ "$NON_INTERACTIVE" = false ]; then
        read -rp "Try to connect anyway? [y/N] " confirm
        [[ "$confirm" =~ ^[Yy]$ ]] || exit 1
        SEC="WPA"   # assume WPA — most common case
    else
        echo "ERROR: SSID not visible in scan. Aborting."
        exit 1
    fi
fi

echo ""
echo "Connecting to: $WIFI_SSID  (security: ${SEC:-(open)})"

# Remove any stale profile with the same name
nmcli connection delete "$WIFI_SSID" 2>/dev/null || true

# Build the connection profile
case "$SEC" in
    *SAE*|*WPA3*)  KEY_MGMT="sae"     ;;
    *WPA*|*PSK*)   KEY_MGMT="wpa-psk" ;;
    ""|--)         KEY_MGMT=""         ;;
    *)
        echo "WARNING: unrecognised security '$SEC', trying wpa-psk"
        KEY_MGMT="wpa-psk"
        ;;
esac

if [ -n "$KEY_MGMT" ]; then
    [ -z "$WIFI_PASS" ] && {
        echo "ERROR: network '$WIFI_SSID' requires a password ($SEC)."
        exit 1
    }
    nmcli connection add type wifi ifname wlan0 con-name "$WIFI_SSID" \
        ssid "$WIFI_SSID" \
        wifi-sec.key-mgmt "$KEY_MGMT" \
        wifi-sec.psk "$WIFI_PASS" \
        wifi-sec.psk-flags 0 \
        connection.autoconnect no \
        > /dev/null
else
    nmcli connection add type wifi ifname wlan0 con-name "$WIFI_SSID" \
        ssid "$WIFI_SSID" \
        connection.autoconnect no \
        > /dev/null
fi

if ! nmcli connection up "$WIFI_SSID"; then
    echo "ERROR: Failed to connect to '$WIFI_SSID'."
    nmcli connection delete "$WIFI_SSID" 2>/dev/null || true
    exit 1
fi

# Wait for an IP address (exclude the AP address range)
IP=""
for i in $(seq 1 20); do
    IP=$(ip -4 addr show wlan0 \
         | grep -oP "(?<=inet )[\d.]+" | grep -v "192\.168\.50\." | head -1)
    [ -n "$IP" ] && break
    sleep 1
done
[ -z "$IP" ] && IP=$(hostname -I | awk '{print $1}')

echo ""
echo "============================================================================="
echo " Connected to LAN."
echo ""
echo "   SSID     : $WIFI_SSID"
echo "   LAN IP   : $IP"
echo "   AP IP    : 192.168.50.1  (still active)"
echo ""
echo "   SSH (LAN): ssh efinder@$IP"
echo "   SSH (AP) : ssh efinder@192.168.50.1"
echo ""
echo "   Run ~/ap.sh to drop LAN and return to AP-only mode."
echo "============================================================================="
