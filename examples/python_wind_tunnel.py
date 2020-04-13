#!/usr/bin/env python3
"""
vane - Python Wind Tunnel Example

Demonstrates using Python bindings to configure, simulate, and visualize
transonic flow over a NACA 2412 cambered airfoil at 5 degrees angle of attack.
"""

import os
import sys
from pathlib import Path

# Add bindings to sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bindings" / "python"))

from vane import AsciiFlowRenderer, FlowPlotter, VaneSimulation


def main():
    print("=" * 72)
    print("  vane: Computational Aerodynamics Wind Tunnel (Python Client)")
    print("=" * 72)

    nx, ny = 160, 80
    chord = 40.0
    aoa = 5.0
    u_inlet = 0.05
    tau = 0.54

    print(f"Lattice domain: {nx} x {ny} nodes")
    print(f"Airfoil: NACA 2412, chord={chord} lu, AoA={aoa} deg")
    print(f"Relaxation time: tau={tau}, Inlet velocity: {u_inlet} lu/ts")

    with VaneSimulation(nx=nx, ny=ny, tau=tau, u_inlet=u_inlet, collision_type=2) as sim:
        print(f"Native engine version: {sim.version}")

        # Rasterize NACA 2412
        sim.set_naca_airfoil(
            code="2412",
            chord=chord,
            aoa_deg=aoa,
            x_lead=35.0,
            y_lead=40.0,
        )

        print("\nStepping wind tunnel simulation forward...")
        total_steps = 250
        step_chunk = 50

        for s in range(0, total_steps, step_chunk):
            sim.step(step_chunk)
            cl, cd, cm = sim.get_forces()
            print(
                f"Step {sim.current_step:4d} | "
                f"Lift CL: {cl:+7.4f} | "
                f"Drag CD: {cd:+7.4f} | "
                f"Moment CM: {cm:+7.4f}"
            )

        print("\nExtracting macroscopic velocity field...")
        speed = sim.get_field("speed")
        mask = sim.get_solid_mask()

        print("\nTerminal ASCII Flow Speed Visualization:")
        print("-" * 72)
        ascii_vis = AsciiFlowRenderer.render_field(
            speed,
            mask,
            out_width=72,
            out_height=20,
        )
        print(ascii_vis)
        print("-" * 72)

        # Output artifact directories
        out_dir = Path("output")
        out_dir.mkdir(exist_ok=True)

        vtk_file = out_dir / "naca2412_wind_tunnel.vti"
        csv_file = out_dir / "naca2412_telemetry.csv"
        png_file = out_dir / "naca2412_flow.png"

        sim.save_vtk(str(vtk_file))
        sim.save_csv(str(csv_file))
        print(f"Exported VTK field to: {vtk_file}")
        print(f"Exported Telemetry CSV to: {csv_file}")

        saved_png = FlowPlotter.save_image_pillow(speed, mask, str(png_file))
        if saved_png:
            print(f"Exported PNG visualization to: {png_file}")

    print("\nSimulation complete. Wind tunnel experiment passed successfully.")


if __name__ == "__main__":
    main()
