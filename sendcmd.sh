#!/bin/bash
# sendcmd.sh -- send an LX200 command to the eFinder server
# Usage: ~/sendcmd.sh <command>
# Example: ~/sendcmd.sh GV    (get version)
#          ~/sendcmd.sh SO    (reset offset)
#          ~/sendcmd.sh GS    (get star count)
#          ~/sendcmd.sh GR    (get current RA)
#          ~/sendcmd.sh GD    (get current Dec)

[ -z "$1" ] && { echo "Usage: $0 <command>  e.g. $0 GV"; exit 1; }

echo -n ":${1}#" | nc -w1 localhost 4060
echo ""
