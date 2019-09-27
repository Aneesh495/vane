#!/usr/bin/env bash
set -euo pipefail

# Scans all tracked and working files for em dashes (U+2014) and en dashes (U+2013)
FOUND=0
while IFS= read -r file; do
    if [ -f "$file" ]; then
        if grep -Pn "[\x{2014}\x{2013}]" "$file" 2>/dev/null; then
            echo "EM-DASH VIOLATION: $file"
            FOUND=1
        fi
    fi
done < <(find . -type f ! -path "./.git/*" ! -path "./target/*" ! -name "*.sh")

if [ "$FOUND" -ne 0 ]; then
    echo "ERROR: Disallowed dashes detected. Replace with hyphens (-)."
    exit 1
else
    echo "SUCCESS: Zero em dashes or en dashes found in codebase."
    exit 0
fi
