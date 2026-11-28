# PLB-NGRC

MATLAB code for parallel locality-blended next-generation reservoir computing (PLB-NGRC) and the baselines P-NGRC, LB-NGRC and parallel reservoir computing (P-RC), accompanying the manuscript

J. Hwang, T. Kim, Y. Do and Y.-C. Lai, *Dual locality for data-driven prediction of spatiotemporal chaos with next-generation reservoir computing*.

The result files of the manuscript (selected hyperparameters, Bayesian optimization histories, predictions and measures) are archived on Zenodo: [https://doi.org/10.5281/zenodo.24548630](https://doi.org/10.5281/zenodo.24548630).

## Quick start: figures from the archived results

1. Download the result archive from [https://doi.org/10.5281/zenodo.24548630](https://doi.org/10.5281/zenodo.24548630) and unzip it in the top folder of this repository. It fills `experiments/<data set>/data`, `experiments/<data set>/results` and `timing/training_time.csv`.
2. Run any figure script, e.g.

```matlab
run('figures/code/Fig_Ikeda.m')
```

The figure is written to `figures/img/figure4.png`.

To reproduce one figure from scratch instead, delete the corresponding `data` and `results` folders and run, e.g., `Run_model_system('Ikeda')` followed by `Fig_Ikeda.m`. `Run_model_system('Ikeda')` also runs the observational-noise and scaling experiments of this system.

## Requirements

- MATLAB R2024a or later
- Statistics and Machine Learning Toolbox (Bayesian optimization, `bayesopt`)
- Parallel Computing Toolbox (parallel Bayesian optimization)
- Mapping Toolbox (coastlines in Figs. 10 and 11 only)
- Internet access for the Z500 data (NCEP/NCAR Reanalysis 1, NOAA PSL), which are downloaded once and cached

## Run

```matlab
Run_all
```

`Run_all` runs the blocks below in order and writes a log to `Run_all_log.txt`. A block that fails does not stop the others. Every block can also be run on its own after adding the folders to the path, e.g. `Run_model_system('Ikeda')`. Interrupted Bayesian optimizations resume from their checkpoints, and finished results are not recomputed.

| Block | Function | Content |
|---|---|---|
| Lorenz-96 system | `Run_model_system('Lorenz96')` | main comparison, ablation, scaling with N |
| Coupled Ikeda map | `Run_model_system('Ikeda')` | as above, and observational noise |
| Thomas system | `Run_model_system('Thomas')` | as above, and observational noise |
| Kuramoto-Sivashinsky equation | `KS_run` | L = 22, 44, 100, 200 and power spectra |
| Z500 along 50N | `Z500_1D_run` | P-NGRC and PLB-NGRC |
| Z500 over 30-70N | `Z500_2D_run` | P-NGRC, PLB-NGRC and the ablation of the latitude and seasonal conditioning |
| Training time | `Measure_training_time` | `timing/training_time.csv` |
| Figures | `figures/code/Fig_*.m` | `figures/img/*.png` |

## Folders

```
core/          NGRC and P-RC models, Bayesian optimization, evaluation (VPT, W1, largest Lyapunov exponent)
experiments/   one folder per data set; each creates data/ and results/ when run
timing/        training time of the selected models
figures/       code/ (one script per figure) and img/ (output)
```

## Result files

Each experiment folder holds `data/` (input time series) and `results/` (results).

| File | Content |
|---|---|
| `data/<system>_N<sites>_data.mat` | normalized time series `S` with the training, validation and test indices (model systems; KS: `KS_L<domain length>_data.mat`; Z500: `Z500_1D_data.mat`, `Z500_2D_data.mat`) |
| `results/<system>_N<sites>_true_lyapunov.mat` | largest Lyapunov exponent of the true system from 10 initial conditions (KS: `KS_L<domain length>_true_lyapunov.mat`) |
| `results/<system>_N<sites>_<method>[_<variant>].mat` | result of one method, with method `PNGRC`, `PLBNGRC`, `LBNGRC` or `PRC` and variant `inball`, `weighted` (ablation) or `noise01`, `noise05`, `noise10` (observational noise in percent) |
| `results/KS_L<domain length>_<method>.mat`, `KS_spectrum.mat` | Kuramoto-Sivashinsky results and power spectra |
| `results/Z500_1D_<method>.mat` | Z500 along 50N |
| `results/Z500_2D_<variant>_best.mat`, `_test.mat` | Z500 over 30-70N: selected hyperparameters `p` and test measures; variants `PLBNGRC`, `PNGRC`, and for the ablation `season_only`, `latitude_only`, `neither`, `in_features` (hyperparameters selected again) or with the prefix `fixed_` (other hyperparameters fixed) |
| `*_bo.mat`, `*_log.txt` | Bayesian optimization checkpoint and console log |

A result file of the model systems holds a structure `res` with the selected hyperparameters (`best`), the validation objective (`val_obj`), the evaluations of the Bayesian optimization (`bo_x`, `bo_f`), the training time (`t_train`), the first long-term prediction (`P_long`), and the measures of the 10 predictions: VPT in Lyapunov times (`vpt`), W1 (`D_all`), largest Lyapunov exponent (`lle_all`), divergence (`div_all`) and the reference value of W1 (`D_floor`).

## Figures and tables

Each script writes the figure of the manuscript with the same number to `figures/img/`. In the manuscript, the panels of Figs. 2 and 4-7 are arranged in a different layout and the labels are typeset again.

| Figure | Script | Output |
|---|---|---|
| 2 Lorenz-96 prediction | `Fig_Lorenz96_prediction.m` | `figure2.png` |
| 3 Lorenz-96 measures | `Fig_Lorenz96_measures.m` | `figure3.png` |
| 4 Coupled Ikeda map | `Fig_Ikeda.m` | `figure4.png` |
| 5 Thomas system | `Fig_Thomas.m` | `figure5.png` |
| 6 Coupled Ikeda map with noise | `Fig_Ikeda_noise.m` | `figure6.png` |
| 7 Thomas system with noise | `Fig_Thomas_noise.m` | `figure7.png` |
| 8 Scaling with N | `Fig_scaling.m` | `figure8.png` |
| 9 Kuramoto-Sivashinsky | `Fig_KS.m` | `figure9.png` |
| 10 Z500 along 50N | `Fig_Z500_1D.m` | `figure10.png` |
| 11 Z500 over 30-70N | `Fig_Z500_2D.m` | `figure11.png` |

Figure 1 is a schematic. The settings, selected hyperparameters and training times (Table I) and the values in the other tables are stored in the result files (`res.best`) and in `timing/training_time.csv`.

## Data

The model-system data are generated with fixed random seeds. The Z500 data are the daily 500 hPa geopotential height of the NCEP/NCAR Reanalysis 1, provided by the NOAA PSL, Boulder, Colorado, USA, from https://psl.noaa.gov.

## License and citation

The code is released under the MIT License (`LICENSE`). If you use it, please cite this code (https://doi.org/10.5281/zenodo.24548886); the citation metadata of the code are in `CITATION.cff`.
