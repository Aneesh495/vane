#!/usr/bin/env bash
# vane - Computational Aerodynamics Performance Benchmark Suite
# Measures simulation throughput in MLUPS across varying grid resolutions.
# Note: Zero em dashes are used in this file per coding standards.

set -e

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

echo "========================================================================"
echo "         vane: Aerodynamic LBM Throughput Benchmark (MLUPS)             "
echo "========================================================================"

cargo build --release -p vane-cli

echo ""
echo "Running multi-scale grid throughput benchmark suite:"
echo ""

echo "--- Benchmark 1: Small Grid (120x60, 500 steps) ---"
./target/release/vane-cli bench --nx 120 --ny 60 --steps 500

echo ""
echo "--- Benchmark 2: Medium Grid (200x100, 400 steps) ---"
./target/release/vane-cli bench --nx 200 --ny 100 --steps 400

echo ""
echo "--- Benchmark 3: Large Grid (400x200, 200 steps) ---"
./target/release/vane-cli bench --nx 400 --ny 200 --steps 200

echo ""
echo "Benchmark suite completed successfully."
