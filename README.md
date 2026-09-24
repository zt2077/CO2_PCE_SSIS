# PCE-assisted SSIS for CO2 injection risk

The MATLAB workflow combines sparse polynomial chaos expansion (PCE), subset simulation, Gaussian mixture models (GMM), and importance sampling to estimate injection-induced mechanical failure probabilities in the Johansen Formation.

## Contents

| Directory | Contents |
| --- | --- |
| `analysis/` | Surrogate training, MCS, and SSIS analysis |
| `visualization/` | Probability comparisons, spatial maps, and conditional GMM plots |
| `examples/` | Two-dimensional SS-GMM-IS illustration |
| `forward_model/` | COMSOL LiveLink evaluation interface |
| `models/` | COMSOL model used by the forward interface |
| `data/training/` | 188 input samples and responses at 3,179 nodes |
| `data/reference_results/` | Archived results for plotting without rerunning the analysis |
| `docs/` | Data definitions, provenance, and reproduction notes |

## Requirements and use

Use MATLAB with Statistics and Machine Learning Toolbox. Analysis additionally requires UQLab (the source installation is version 2.0). COMSOL 6.2 with LiveLink for MATLAB and the relevant physics modules is needed only to rerun the forward model. MATLAB R2024a was used for package checks.

Open this folder in MATLAB, then run:

```matlab
run_workflow('figures');       % Plot the bundled reference results
run_workflow('benchmark');     % Run the illustrative 2D example

% After installing UQLab and adding its core folder to the MATLAB path:
result_dir = run_workflow('analysis');
run_workflow('figures', result_dir);

% Optional: from a MATLAB session connected to COMSOL LiveLink:
run_workflow('forward');
```

Each run writes to a new directory under `outputs/`. Input data are not overwritten. The main analysis uses the first 100 archived training samples, an injection rate of 20 kg/s, up to 100,000 MCS samples, and 10,000 final IS samples.

The plotting and 2D workflows were executed successfully; the full surrogate analysis and COMSOL simulations were not rerun during packaging. Read [reproduction notes](docs/REPRODUCTION.md) for implementation details and differences from the manuscript before claiming exact reproduction.

## License and citation

The project code is released under the [MIT License](LICENSE). Third-party software and the COMSOL-derived model retain their respective terms; see [third-party notices](docs/THIRD_PARTY.md). Please cite the accompanying manuscript when using this work. Publication details and a repository DOI can be added after release.
