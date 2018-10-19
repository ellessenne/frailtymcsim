# Impact of model misspecification in shared frailty survival models

This repository contains the code required to simulate data, fit each model, and produce summary tables and figures for the manuscript titled _Impact of model misspecification in shared frailty survival models_ by Gasparini et al. In addition to that, this repository contains the Supporting Web Material and a dataset with the full results of the simulation study.

A preprint of the manuscript is available on [arXiv](https://arxiv.org/abs/1810.08140).

The following files are included in this repository:

* `make_data_dgms.R`, an R script file used to simulate the data under each data-generating mechanism;
* `make_models.R`, an R script file used to fit the models;
* `make_summaries.R`, an R script file used to aggregate the results and obtain summary statistics for each model and data-generating mechanism;
* `Results.RDS`, a serialised R object containing the raw results of the simulation study;
* `summary.RDS`, a serialised R object containing the summary statistics for each model and data-generating mechanism. `summary.RDS` was produced by `make_summaries.R`;
* `suppl.pdf`, a PDF document with the Supporting Web Material.
