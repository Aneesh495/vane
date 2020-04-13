"""
Flow visualization utilities for vane aerodynamic simulations.

Provides ASCII terminal rendering and high-resolution raster export.
"""

import math
from typing import List, Optional, Tuple, Union


class AsciiFlowRenderer:
    """
    Renders 2D macroscopic flow fields directly to ASCII characters in the terminal.
    """

    CHAR_RAMP = " .:-=+*#%@"

    @classmethod
    def render_field(
        cls,
        field: any,
        solid_mask: any,
        out_width: int = 80,
        out_height: int = 24,
        vmin: Optional[float] = None,
        vmax: Optional[float] = None,
    ) -> str:
        """
        Downsample and format a 2D field into an ASCII string grid.

        Parameters
        ----------
        field : 2D list or numpy array [ny, nx]
            Field scalar values to visualize.
        solid_mask : 2D list or numpy array [ny, nx]
            Obstacle mask (1 = solid, 0 = fluid).
        out_width : int
            Number of character columns in output.
        out_height : int
            Number of character rows in output.
        """
        try:
            ny = len(field)
            nx = len(field[0])
        except (TypeError, IndexError):
            return ""

        if nx <= 0 or ny <= 0:
            return ""

        # Flatten or scan min/max if not supplied
        if vmin is None or vmax is None:
            auto_min = 1e30
            auto_max = -1e30
            for y in range(ny):
                for x in range(nx):
                    if solid_mask[y][x] == 0:
                        v = field[y][x]
                        if not math.isnan(v):
                            if v < auto_min:
                                auto_min = v
                            if v > auto_max:
                                auto_max = v
            if vmin is None:
                vmin = auto_min if auto_min < auto_max else 0.0
            if vmax is None:
                vmax = auto_max if auto_max > auto_min else 1.0

        vrange = vmax - vmin if vmax > vmin else 1e-6
        num_chars = len(cls.CHAR_RAMP)

        lines: List[str] = []
        # In aerodynamics, y=0 is bottom, y=ny-1 is top. Terminal prints top to bottom.
        for row in range(out_height - 1, -1, -1):
            y_src = int(row * (ny - 1) / max(1, out_height - 1))
            line_chars = []
            for col in range(out_width):
                x_src = int(col * (nx - 1) / max(1, out_width - 1))

                if solid_mask[y_src][x_src] != 0:
                    line_chars.append("#")
                else:
                    val = field[y_src][x_src]
                    norm = (val - vmin) / vrange
                    clamped = max(0.0, min(0.999, norm))
                    char_idx = int(clamped * num_chars)
                    line_chars.append(cls.CHAR_RAMP[char_idx])
            lines.append("".join(line_chars))

        return "\n".join(lines)


class FlowPlotter:
    """
    Export raster plots of flow fields, vortex cores, and aerodynamic forces.
    """

    @staticmethod
    def save_image_pillow(
        field: any,
        solid_mask: any,
        filepath: str,
        colormap: str = "viridis",
    ) -> bool:
        """
        Render 2D field to PNG using Pillow.
        """
        try:
            from PIL import Image
        except ImportError:
            return False

        ny = len(field)
        nx = len(field[0])

        # Compute field min and max on fluid nodes
        vmin = 1e30
        vmax = -1e30
        for y in range(ny):
            for x in range(nx):
                if solid_mask[y][x] == 0:
                    v = field[y][x]
                    if not math.isnan(v):
                        if v < vmin:
                            vmin = v
                        if v > vmax:
                            vmax = v
        vrange = vmax - vmin if vmax > vmin else 1e-6

        img = Image.new("RGB", (nx, ny))
        pixels = img.load()

        for y in range(ny):
            img_y = ny - 1 - y  # Flip vertically so y=0 is bottom
            for x in range(nx):
                if solid_mask[y][x] != 0:
                    # Solid obstacle rendered as dark gray
                    pixels[x, img_y] = (40, 40, 40)
                else:
                    norm = max(0.0, min(1.0, (field[y][x] - vmin) / vrange))
                    if colormap == "inferno":
                        r = int(255 * math.sqrt(norm))
                        g = int(255 * (norm**2))
                        b = int(255 * (1.0 - norm) * 0.5)
                    else:  # Cool to warm
                        r = int(255 * norm)
                        g = int(180 * (1.0 - abs(norm - 0.5) * 2))
                        b = int(255 * (1.0 - norm))
                    pixels[x, img_y] = (r, g, b)

        img.save(filepath)
        return True
