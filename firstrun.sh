#!/bin/bash
# firstrun.sh — staged into the image at build time, runs once on first boot.
# Waits for internet, clones the repo, runs install.sh, then reboots.
set -e
LOG=/home/efinder/firstrun.log
exec > >(tee -a "$LOG") 2>&1

echo "=== eFinder first-boot setup: $(date) ==="

# Wait for internet (needed to clone repo and run apt)
echo "Waiting for internet connectivity..."
connected=false
for i in $(seq 1 30); do
    if ping -c1 -W2 github.com &>/dev/null; then
        connected=true
        echo "  Internet reachable after attempt $i"
        break
    fi
    echo "  attempt $i/30 -- retrying in 10s..."
    sleep 10
done

if [ "$connected" != "true" ]; then
    echo "ERROR: No internet after 5 minutes -- aborting first-boot setup."
    echo "       Ensure the Pi has WiFi or USB ethernet on first boot."
    exit 1
fi

# Fix ownership now that the efinder user exists at runtime
chown -R efinder:efinder /home/efinder

# Load baked-in credentials from config.env
source /home/efinder/config.env
export WIFI_PASS SAMBA_PASS

# Clone the repo (tiny_img branch)
cd /home/efinder
if [ -d eFinder_cli ]; then
    echo "Repo already present -- pulling latest..."
    git -C eFinder_cli fetch origin
    git -C eFinder_cli reset --hard origin/tiny_img
else
    git clone --branch tiny_img \
        https://github.com/mconsidine/eFinder_cli.git eFinder_cli
fi
chown -R efinder:efinder /home/efinder/eFinder_cli

# Run installer non-interactively
bash /home/efinder/eFinder_cli/install.sh --non-interactive

# Remove trigger file so firstrun.service skips on all subsequent boots.
# To re-run: touch /home/efinder/.run_firstrun && sudo reboot
rm -f /home/efinder/.run_firstrun

echo "=== eFinder first-boot setup complete: $(date) ==="
reboot
