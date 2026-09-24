# Reproduction notes

## Selected sources

The package follows `paper writing/Manuscript.docx` and the associated `paper writing/draw pictures` scripts in the supplied research directory. Earlier Bayesian-inversion trials, PCA experiments, backups, unrelated UQLab examples, manuscript drafts, and reference papers are excluded. Source paths relative to that directory and SHA-256 checksums are recorded in `source_manifest.csv`.

The main calculation comes from `draw pictures/run_analysis.m`. The comparison figure uses the accuracy-analysis function from `create_visualizations_GMM2.m`; unused experimental functions were omitted. Conditional GMM plots come from `create_visualizations_GMM1_modify.m`, and the 2D illustration from `toy_example_complete2.m`. The COMSOL interface and model come from `test/modify20_subset_sobol`.

## Data and defaults

The MAT files preserve their original arrays and variable names. Rows correspond to samples and columns to nodes unless stated otherwise.

| File | Variable and shape | Meaning |
| --- | --- | --- |
| `X_data_sobol_188all1.mat` | `X`, 188 x 4 | E (Pa), nu, cohesion (Pa), friction angle (rad) |
| `Yall_data_sobol_188all1.mat` | `Yall`, 188 x 3179 | Vertical displacement (m) |
| `S1_sobol_188all1.mat` | `dataS1_all`, 188 x 3179 | Principal stress 1 (Pa) |
| `S3_sobol_188all1.mat` | `dataS3_all`, 188 x 3179 | Principal stress 3 (Pa) |
| `coordinate.mat` | `coordinate`, 3 x 3179 | Node coordinates (m); copied from `coordinate1.mat` |

Independent Gaussian inputs have means `[8.9e9, 0.25, 3e6, (3.14/180)*28.4]` and standard deviations `[2e7, 0.02, 5e5, 3.14/180]`. The original 3.14 approximation is retained. The analysis seed is 100; the 2D illustration uses 42. Analysis settings remain near the beginning of `analysis/run_analysis.m` and in its MCS/SSIS sections.

The supplied reference results contain 3 GMM components, 10,000 IS samples, and ESS approximately 2069.72. Their MCS evaluation-count field is a per-node vector with maximum 100,000. Newly generated results store that maximum as a scalar. The original result files are unchanged.

## Implementation details to reconcile with the manuscript

These observations describe the selected source version; the packaging does not change the scientific estimator to match the prose.

- **Training design:** the active script uses the first 100 of 188 archived rows. Related generation scripts use Sobol sampling, and the archived filenames contain `sobol`; the manuscript describes Latin hypercube training. The exact historical generation of all 188 rows has not been independently reconstructed.
- **Surrogates:** Y and S1 use LARS PCE with degrees 6-8, qNorm 0.6, and maximum interaction 2. S3 uses a constant model for relative range below 2%, linear regression below 5%, and PCE otherwise (degrees 2-7, qNorm 0.5). The manuscript describes PCE with qNorm 0.5 more generally. Old comments claiming that active S3 models are Kriging were corrected.
- **Error metrics:** the S3 field labeled LOO contains relative standard deviation for the constant model, mean leave-one-out relative absolute error for the linear model, or UQLab LOO for PCE. These are different metrics. Five-fold CV is skipped in the selected script.
- **Failure convention:** local failure is `g >= -1e-6`, with `g = -S3 - (-S1*t^2 + 2*c*t)` and `t = tan(phi/2 + pi/4)`. Subset simulation uses the rowwise 0.9 quantile of nodal g values; reported Pf arrays are nodewise estimates.
- **Importance sampling:** the selected main script uses self-normalized importance weights. This is not the ordinary unbiased IS estimator described in parts of the manuscript. GMM/wide-GMM/standard-normal weights are 0.6/0.2/0.2. `tau=2` multiplies covariance by `tau^2=4`, so it doubles standard deviation.
- **Cost comparison:** the 10,000 IS draws exclude the 3,000 pilot samples, subset-simulation evaluations, and surrogate training. They are not the total number of evaluations for the complete method.
- **Reference provenance:** archived results are bundled as supplied, but an exact match between them and a fresh execution of the selected source has not been established. The COMSOL model is the matching named candidate used by the available dataset-extension script; all historical training solves were not replayed.

## Packaging changes and checks

Chinese comments, messages, and plot labels were translated. Scripts were wrapped as functions and dispatched by `run_workflow`; the caller's current directory and MATLAB path are restored. Machine-specific output paths were replaced, repeated filename suffixes were unified, and a spatial-map function was added. GMM plots support more than three components and reserve space for axis labels.

Two diagnostic fixes were made: relative-error plots exclude zero MCS denominators (the original reported infinite MAPE), and the toy example restores the removed log-weight scale when printing the ordinary IS probability. The toy visualization still uses the original scaled weights, so its sampling/plotting procedure is unchanged. No archived probabilities or training arrays were modified.

Checks include MATLAB R2024a Code Analyzer, finite-value and dimension checks, coordinate/result consistency, and execution of the reference-result plotting and 2D example. Code Analyzer reports style/performance notices, but no syntax errors. Selected source-file hashes are checked against the originals. The full 3,179-node analysis and the COMSOL forward solves have not been executed during packaging.

## Optional forward evaluation

Start a COMSOL LiveLink session and run `run_workflow('forward')`. This evaluates the 188 archived inputs, rather than inventing a replacement training design. It calls solver sequences `sol5` and `sol6`, exports `w`, `solid.sp1`, and `solid.sp3` from `dset6`, and uses export columns 10-12 as in the source interface. An explicit check detects unexpected export dimensions. Verify the intended output time and node order when changing the model.

Forward results are written as `X_generated.mat`, `Yall_generated.mat`, `S1_generated.mat`, `S3_generated.mat`, and `coordinate.mat` in their own output directory. They do not replace the archived inputs used by the main analysis automatically.
