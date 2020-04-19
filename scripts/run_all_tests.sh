#!/usr/bin/env bash
# vane - Unified Multi-Language Test Runner
# Executes all zero-em-dash checks, Rust tests, BQN suites, and Python workflows.
# Note: Zero em dashes are used in this file per coding standards.

set -e

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

echo "========================================================================"
echo "                   vane: Full System Verification Suite                 "
echo "========================================================================"

echo ""
echo "[Step 1/5] Verifying zero em dashes and en dashes across repository..."
./scripts/verify_no_em_dash.sh

echo ""
echo "[Step 2/5] Running Rust core unit and integration test suite..."
cargo test --all

echo ""
echo "[Step 3/5] Executing BQN array-native numerical and physics test suites..."
bqn tests/bqn/run_all_bqn_tests.bqn

echo ""
echo "[Step 4/5] Running Python FFI wind tunnel simulation and visualization..."
python3 examples/python_wind_tunnel.py

echo ""
echo "[Step 5/5] Running vane-cli terminal simulation workstation..."
cargo run -p vane-cli -- run --chord 30 --aoa 4 --steps 100 --csv output/cli_smoke.csv --vtk output/cli_smoke.vti
cargo run -p vane-cli -- verify

echo ""
echo "========================================================================"
echo "  SUCCESS: All test stages passed with 100% verification fidelity!      "
echo "========================================================================"
