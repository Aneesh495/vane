# VTK Structured XML (.vti) and Legacy (.vtk) Export Engine
#
# Generates scientific visualization files for ParaView, VisIt, and Blender:
#   - XML Structured ImageData (.vti) format with LittleEndian ASCII payload
#   - PointData fields: Density, Pressure, Vorticity, Velocity vectors [ux, uy, uz], SolidMask
#   - Legacy ASCII Structured Points (.vtk) format
#
# Note: Zero em dashes are used in this file per coding standards.

math ← •Import "../lib/math_utils.bqn"
stencils ← •Import "../lib/stencils.bqn"

# Newline string
nl ← ⟨@+10⟩

# Format flat float vector into space-delimited string
FormatVector ← {
  vec ← ⥊ 𝕩
  ∾ { (•Fmt 𝕩) ∾ " " }¨ vec
}

# Format 3D velocity vectors [ux, uy, 0.0] as triples
FormatVelocityTriples ← {
  ⟨ux, uy⟩ ← 𝕩
  flatUx ← ⥊ ux
  flatUy ← ⥊ uy
  n ← ≠flatUx

  triples ← {
    i ← 𝕩
    u ← i ⊑ flatUx
    v ← i ⊑ flatUy
    (•Fmt u) ∾ " " ∾ (•Fmt v) ∾ " 0.0" ∾ nl
  }¨ (↕n)

  ∾ triples
}

# Generate complete XML VTK (.vti) string
# Arguments:
#   rho: [Ny, Nx] fluid density
#   ux: [Ny, Nx] horizontal velocity
#   uy: [Ny, Nx] vertical velocity
#   solidMask: [Ny, Nx] solid obstacle mask (1 = solid, 0 = fluid)
#
# Returns:
#   vtiContent: ASCII XML string ready to save to .vti file
FormatVTI ⇐ {
  ⟨rho, ux, uy, solidMask⟩ ← 𝕩
  ⟨ny, nx⟩ ← ≢rho

  # Compute pressure and vorticity
  cs2 ← 1.0 ÷ 3.0
  pressure ← cs2 × rho
  vortZ ← 1.0‿1.0 stencils.Curl2D ⟨ux, uy⟩

  nxM1 ← •Fmt (nx - 1)
  nyM1 ← •Fmt (ny - 1)
  extentStr ← "0 " ∾ nxM1 ∾ " 0 " ∾ nyM1 ∾ " 0 0"

  header ← ∾ ⟨
    "<?xml version=""1.0""?>" ∾ nl,
    "<VTKFile type=""ImageData"" version=""0.1"" byte_order=""LittleEndian"">" ∾ nl,
    "  <ImageData WholeExtent=""" ∾ extentStr ∾ """ Origin=""0 0 0"" Spacing=""1 1 1"">" ∾ nl,
    "    <Piece Extent=""" ∾ extentStr ∾ """>" ∾ nl,
    "      <PointData Scalars=""Density"" Vectors=""Velocity"">" ∾ nl
  ⟩

  rhoBlock ← ∾ ⟨
    "        <DataArray type=""Float32"" Name=""Density"" NumberOfComponents=""1"" format=""ascii"">" ∾ nl,
    "          " ∾ (FormatVector rho) ∾ nl,
    "        </DataArray>" ∾ nl
  ⟩

  pBlock ← ∾ ⟨
    "        <DataArray type=""Float32"" Name=""Pressure"" NumberOfComponents=""1"" format=""ascii"">" ∾ nl,
    "          " ∾ (FormatVector pressure) ∾ nl,
    "        </DataArray>" ∾ nl
  ⟩

  vortBlock ← ∾ ⟨
    "        <DataArray type=""Float32"" Name=""Vorticity"" NumberOfComponents=""1"" format=""ascii"">" ∾ nl,
    "          " ∾ (FormatVector vortZ) ∾ nl,
    "        </DataArray>" ∾ nl
  ⟩

  maskBlock ← ∾ ⟨
    "        <DataArray type=""Int32"" Name=""SolidMask"" NumberOfComponents=""1"" format=""ascii"">" ∾ nl,
    "          " ∾ (FormatVector solidMask) ∾ nl,
    "        </DataArray>" ∾ nl
  ⟩

  velBlock ← ∾ ⟨
    "        <DataArray type=""Float32"" Name=""Velocity"" NumberOfComponents=""3"" format=""ascii"">" ∾ nl,
    (FormatVelocityTriples ⟨ux, uy⟩),
    "        </DataArray>" ∾ nl
  ⟩

  footer ← ∾ ⟨
    "      </PointData>" ∾ nl,
    "      <CellData>" ∾ nl,
    "      </CellData>" ∾ nl,
    "    </Piece>" ∾ nl,
    "  </ImageData>" ∾ nl,
    "</VTKFile>" ∾ nl
  ⟩

  header ∾ rhoBlock ∾ pBlock ∾ vortBlock ∾ maskBlock ∾ velBlock ∾ footer
}

# Export VTI string to file path on disk
# Arguments:
#   filePath: destination file path (e.g. "output.vti")
#   payload: ⟨rho, ux, uy, solidMask⟩
ExportVTI ⇐ {
  ⟨filePath, rho, ux, uy, solidMask⟩ ← 𝕩
  xmlStr ← FormatVTI ⟨rho, ux, uy, solidMask⟩
  filePath •FChars xmlStr
  filePath
}

# Generate Legacy ASCII (.vtk) format
FormatLegacyVTK ⇐ {
  ⟨rho, ux, uy, solidMask⟩ ← 𝕩
  ⟨ny, nx⟩ ← ≢rho
  nTotal ← •Fmt (ny × nx)

  cs2 ← 1.0 ÷ 3.0
  pressure ← cs2 × rho
  vortZ ← 1.0‿1.0 stencils.Curl2D ⟨ux, uy⟩

  header ← ∾ ⟨
    "# vtk DataFile Version 3.0" ∾ nl,
    "vane lattice boltzmann simulation" ∾ nl,
    "ASCII" ∾ nl,
    "DATASET STRUCTURED_POINTS" ∾ nl,
    "DIMENSIONS " ∾ (•Fmt nx) ∾ " " ∾ (•Fmt ny) ∾ " 1" ∾ nl,
    "ORIGIN 0 0 0" ∾ nl,
    "SPACING 1 1 1" ∾ nl,
    "POINT_DATA " ∾ nTotal ∾ nl
  ⟩

  rhoBlock ← ∾ ⟨ "SCALARS Density float 1" ∾ nl, "LOOKUP_TABLE default" ∾ nl, (FormatVector rho) ∾ nl ⟩
  pBlock ← ∾ ⟨ "SCALARS Pressure float 1" ∾ nl, "LOOKUP_TABLE default" ∾ nl, (FormatVector pressure) ∾ nl ⟩
  vortBlock ← ∾ ⟨ "SCALARS Vorticity float 1" ∾ nl, "LOOKUP_TABLE default" ∾ nl, (FormatVector vortZ) ∾ nl ⟩
  velBlock ← ∾ ⟨ "VECTORS Velocity float" ∾ nl, (FormatVelocityTriples ⟨ux, uy⟩) ⟩

  header ∾ rhoBlock ∾ pBlock ∾ vortBlock ∾ velBlock
}
