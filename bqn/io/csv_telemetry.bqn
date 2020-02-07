# Time-Series Aerodynamic Telemetry and Convergence Logger
#
# Records, formats, and exports time histories of aerodynamic force coefficients,
# convergence residuals, mass conservation error, and unsteady wake shedding frequencies:
#   - Instantaneous coefficients: Cl(t), Cd(t), Cm(t), Cn(t)
#   - Numerical convergence: L-inf velocity residual ||u^{n+1} - u^n||_inf
#   - Physical conservation: global mass residual |sum(rho) - sum(rho0)| / sum(rho0)
#   - Acoustic/vortex shedding: Strouhal number St = f_shed * D / U_inf
#   - Statistical averages and RMS fluctuations: Cl_bar, Cd_bar, Cl_rms, Cd_rms
#
# Note: Zero em dashes are used in this file per coding standards.

math ← •Import "../lib/math_utils.bqn"

nl ← ⟨@+10⟩

# Standard CSV header line
telemetryHeader ⇐ "step,cl,cd,cm,cn,res_u,mass_err,max_mach" ∾ nl

# Format a single telemetry entry tuple into CSV line
# Tuple: ⟨step, cl, cd, cm, cn, resU, massErr, maxMach⟩
FormatTelemetryRow ⇐ {
  ⟨step, cl, cd, cm, cn, resU, massErr, maxMach⟩ ← 𝕩

  ∾ ⟨
    •Fmt step, ",",
    •Fmt cl, ",",
    •Fmt cd, ",",
    •Fmt cm, ",",
    •Fmt cn, ",",
    •Fmt resU, ",",
    •Fmt massErr, ",",
    •Fmt maxMach, nl
  ⟩
}

# Format full history table into complete CSV string
# Arguments:
#   historyList: list of telemetry tuples
# Returns:
#   csvString: multi-line CSV string starting with header
FormatTelemetryCSV ⇐ {
  historyList ← 𝕩
  rows ← { FormatTelemetryRow 𝕩 }¨ historyList
  telemetryHeader ∾ (∾ rows)
}

# Export telemetry table to CSV file on disk
# Arguments:
#   filePath: destination CSV file path (e.g. "telemetry.csv")
#   historyList: list of telemetry tuples
ExportTelemetryCSV ⇐ {
  ⟨filePath, historyList⟩ ← 𝕩
  csvStr ← FormatTelemetryCSV historyList
  filePath •FChars csvStr
  filePath
}

# Statistical Reduction of Time Series History
# Computes time-averaged coefficients, RMS fluctuations, and aerodynamic efficiency L/D
# Arguments:
#   historyList: list of tuples ⟨step, cl, cd, cm, cn, resU, massErr, maxMach⟩
#   sampleFraction: fraction of final steps to use for averaging (e.g. 0.5 for final 50%)
#
# Returns list:
#   ⟨clMean, cdMean, cmMean, clRMS, cdRMS, liftToDrag, isConverged⟩
ComputeStatistics ⇐ {
  ⟨historyList, sampleFraction⟩ ← 𝕩
  nTotal ← ≠historyList

  # Extract window of final samples
  nSamples ← 1 ⌈ (⌊ (sampleFraction × nTotal))
  window ← (0 - nSamples) ↑ historyList

  clVals ← { 1 ⊑ 𝕩 }¨ window
  cdVals ← { 2 ⊑ 𝕩 }¨ window
  cmVals ← { 3 ⊑ 𝕩 }¨ window
  resUVals ← { 5 ⊑ 𝕩 }¨ window

  # Time averages
  clMean ← (+´ clVals) ÷ nSamples
  cdMean ← (+´ cdVals) ÷ nSamples
  cmMean ← (+´ cmVals) ÷ nSamples

  # RMS fluctuations: sqrt( mean( (x - x_mean)^2 ) )
  clDiff ← clVals - clMean
  cdDiff ← cdVals - cdMean
  clRMS ← (((+´ (clDiff × clDiff)) ÷ nSamples) ⌈ 0.0) ⋆ 0.5
  cdRMS ← (((+´ (cdDiff × cdDiff)) ÷ nSamples) ⌈ 0.0) ⋆ 0.5

  # Aerodynamic efficiency: L / D = Cl_mean / Cd_mean
  denomCd ← cdMean ⌈ 1e¯12
  liftToDrag ← clMean ÷ denomCd

  # Convergence flag: final velocity residual < 1e-5
  finalResU ← ¯1 ⊑ resUVals
  isConverged ← finalResU < 1e¯5

  ⟨clMean, cdMean, cmMean, clRMS, cdRMS, liftToDrag, isConverged⟩
}

# Strouhal Number St = f_shed * D / U_inf from Zero Crossings of Lift Fluctuation
# Arguments:
#   clHistory: [N] time series of lift coefficient Cl(t)
#   dt: time step between consecutive samples
#   dChar: characteristic dimension (e.g. cylinder diameter or airfoil thickness)
#   uInf: freestream reference velocity
#
# Returns:
#   strouhalNumber: non-dimensional vortex shedding frequency
ComputeStrouhal ⇐ {
  ⟨clHistory, dt, dChar, uInf⟩ ← 𝕩
  n ← ≠clHistory

  clMean ← (+´ clHistory) ÷ (n ⌈ 1)
  clFluc ← clHistory - clMean

  # Zero crossings: count sign transitions where (cl[i] * cl[i-1]) < 0
  signs ← clFluc > 0
  signDiff ← (1 ↓ signs) ≠ (¯1 ↓ signs)
  nCrossings ← +´ signDiff

  # Number of complete cycles is half the number of zero crossings
  nCycles ← nCrossings ÷ 2.0
  totalTime ← (0.0 + n) × dt
  denomTime ← totalTime ⌈ 1e¯12
  shedFrequency ← nCycles ÷ denomTime

  # St = f * D / U_inf
  denomU ← uInf ⌈ 1e¯12
  strouhalNumber ← (shedFrequency × dChar) ÷ denomU
  strouhalNumber
}
