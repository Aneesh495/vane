"""
vane: High-Performance Computational Fluid Dynamics and Aerodynamics Workstation

Python bindings interfacing with the vane-bridge C-ABI native shared library
and BQN array engine.
"""

from .bridge import load_vane_bridge
from .simulation import VaneSimulation
from .visualize import AsciiFlowRenderer, FlowPlotter

__version__ = "0.1.0"
__all__ = [
    "load_vane_bridge",
    "VaneSimulation",
    "AsciiFlowRenderer",
    "FlowPlotter",
]
