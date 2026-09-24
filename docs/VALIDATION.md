# Packaging validation

Checked on 24 September 2026 using MATLAB R2024a on Windows.

- MATLAB Code Analyzer: no syntax errors in the eight MATLAB files; existing style/performance notices remain.
- Training data: 188 x 4 inputs; three 188 x 3179 response arrays; 3 x 3179 coordinates. All checked values are finite.
- Archived result coordinates match the supplied training coordinates; the saved difference equals SSIS Pf minus MCS Pf.
- `run_workflow('figures')`: completed; generated the comparison figure, two spatial maps, twelve conditional GMM isosurface plots, and the center-coordinate CSV. Representative figures were visually inspected.
- Reference comparison: R squared 0.9925, MAE approximately 4.99e-4, RMSE approximately 8.06e-4. MAPE is approximately 8.61% over the 3178 nodes with strictly positive MCS probability; one zero-probability node is excluded only from relative-error metrics.
- `run_workflow('benchmark')`: completed with seed 42 and generated the three-stage illustration. The corrected ordinary IS diagnostic is approximately 3.045e-3; this smoke run is not an independent convergence study.
- Source preservation: SHA-256 hashes of all 18 selected source files match the values recorded before packaging. Copied input data and the COMSOL model retain their original hashes.
- No Chinese characters or machine-specific drive paths remain in the packaged MATLAB source files.

The full surrogate-training/SSIS workflow and COMSOL forward simulations were not rerun. Successful plotting validates the packaged saved-data workflow; it does not establish exact regeneration of the archived scientific results. See `REPRODUCTION.md` for the manuscript/code differences.
