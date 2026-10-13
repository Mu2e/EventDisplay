#!/bin/bash

# Initialize variables
PORT_NUM=""
MACHINE_ID=""
USERNAME=""

# --- Function to open browser based on OS ---
open_browser() {
    local url=$1
    case "$(uname -s)" in
        Darwin)
            # macOS
            echo "Detected macOS. Opening Google Chrome..."
            open -a "Google Chrome" "$url"
            ;;
        Linux)
            # Linux variants
            echo "Detected Linux. Attempting to open browser..."
            # Try specific commands for Chrome/Chromium, fall back to xdg-open
            if command -v google-chrome >/dev/null; then
                google-chrome "$url"
            elif command -v chromium-browser >/dev/null; then
                chromium-browser "$url"
            else
                xdg-open "$url"
            fi
            ;;
        CYGWIN*|MINGW*|MSYS*)
            # Windows (Git Bash/Cygwin)
            echo "Detected Windows. Opening browser..."
            cmd.exe /c start "$url"
            ;;
        *)
            echo "Unsupported OS. Please open the URL manually: $url"
            ;;
    esac
}

# --- Parse command-line arguments ---
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --port)
            PORT_NUM="$2"
            shift; shift
            ;;
        --machine)
            MACHINE_ID="$2"
            shift; shift
            ;;
        --user)
            USERNAME="$2"
            shift; shift
            ;;
        *)
            echo "Unknown parameter passed: $1"
            echo "Usage: $0 --port <number> --machine <id> --user <name>"
            exit 1
            ;;
    esac
done

# --- Validation ---
if [ -z "$PORT_NUM" ] || [ -z "$MACHINE_ID" ] || [ -z "$USERNAME" ]; then
    echo "Error: Missing required arguments."
    echo "Usage: $0 --port <number> --machine <id> --user <name>"
    echo "Example: $0 --port 1234 --machine 05 --user sophie"
    exit 1
fi

# --- Execution ---
HOSTNAME="${MACHINE_ID}.fnal.gov"
URL="http://localhost:${PORT_NUM}/win1/"

open_browser "$URL"

echo "Establishing SSH tunnel for ${USERNAME}@${HOSTNAME} on port ${PORT_NUM}"
ssh -KXY -L "${PORT_NUM}:localhost:${PORT_NUM}" "${USERNAME}@${HOSTNAME}"
