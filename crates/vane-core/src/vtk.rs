//! ParaView VTK Structured XML (.vti) exporter.
//!
//! Formats simulation grid data arrays into standardized XML ImageData
//! files readable by ParaView, VisIt, and Blender.
//!
//! Note: Zero em dashes are used in this file per coding standards.

use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

/// Export structured simulation fields to XML VTK (.vti) file.
pub fn export_vti<P: AsRef<Path>>(
    path: P,
    nx: usize,
    ny: usize,
    rho: &[f64],
    ux: &[f64],
    uy: &[f64],
    vorticity: &[f64],
    mask: &[bool],
) -> io::Result<()> {
    let mut file = File::create(path)?;

    let extent_str = format!("0 {} 0 {} 0 0", nx - 1, ny - 1);

    writeln!(file, "<?xml version=\"1.0\"?>")?;
    writeln!(file, "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">")?;
    writeln!(file, "  <ImageData WholeExtent=\"{}\" Origin=\"0 0 0\" Spacing=\"1 1 1\">", extent_str)?;
    writeln!(file, "    <Piece Extent=\"{}\">", extent_str)?;
    writeln!(file, "      <PointData Scalars=\"Density\" Vectors=\"Velocity\">")?;

    // Density Array
    writeln!(file, "        <DataArray type=\"Float32\" Name=\"Density\" NumberOfComponents=\"1\" format=\"ascii\">")?;
    for val in rho {
        write!(file, "{:.6} ", val)?;
    }
    writeln!(file, "\n        </DataArray>")?;

    // Pressure Array (p = cs2 * rho = rho / 3)
    writeln!(file, "        <DataArray type=\"Float32\" Name=\"Pressure\" NumberOfComponents=\"1\" format=\"ascii\">")?;
    for val in rho {
        write!(file, "{:.6} ", val / 3.0)?;
    }
    writeln!(file, "\n        </DataArray>")?;

    // Vorticity Array
    writeln!(file, "        <DataArray type=\"Float32\" Name=\"Vorticity\" NumberOfComponents=\"1\" format=\"ascii\">")?;
    for val in vorticity {
        write!(file, "{:.6} ", val)?;
    }
    writeln!(file, "\n        </DataArray>")?;

    // Solid Mask Array
    writeln!(file, "        <DataArray type=\"Int32\" Name=\"SolidMask\" NumberOfComponents=\"1\" format=\"ascii\">")?;
    for &is_solid in mask {
        write!(file, "{} ", if is_solid { 1 } else { 0 })?;
    }
    writeln!(file, "\n        </DataArray>")?;

    // Velocity 3D Vector Array
    writeln!(file, "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">")?;
    for i in 0..(nx * ny) {
        writeln!(file, "{:.6} {:.6} 0.0", ux[i], uy[i])?;
    }
    writeln!(file, "        </DataArray>")?;

    writeln!(file, "      </PointData>")?;
    writeln!(file, "      <CellData>")?;
    writeln!(file, "      </CellData>")?;
    writeln!(file, "    </Piece>")?;
    writeln!(file, "  </ImageData>")?;
    writeln!(file, "</VTKFile>")?;

    Ok(())
}
