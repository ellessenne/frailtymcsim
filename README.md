
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Impact of Model Misspecification in Shared Frailty Survival Models

This repository contains the code required to simulate data, fit each
model, and produce summary tables and figures for the manuscript titled
*Impact of model misspecification in shared frailty survival models* by
Gasparini et al. In addition to that, this repository contains the
Supporting Web Material and a dataset with the full results of the
simulation study.

The final, published manuscript is available from the Statistics in
Medicine website: <https://doi.org/10.1002/sim.8309>

A pre-prints history is available on
[arXiv](https://arxiv.org/abs/1810.08140); please get in touch via
[Twitter](https://twitter.com/ellessenne) or
[e-mail](mailto:ag475@leicester.ac.uk?subject=Preprint:%20Impact%20of%20Model%20Misspecification%20in%20Shared%20Frailty%20Survival%20Models)
if you have any comment.

The following files are included in this repository:

  - `make_data_dgms.R`, an R script file used to simulate the data under
    each data-generating mechanism;
  - `make_models.R`, an R script file used to fit the models;
  - `make_summaries.R`, an R script file used to aggregate the results
    and obtain summary statistics for each model and data-generating
    mechanism;
  - `results-nn.RDS`, `nn` serialised R objects containing the raw
    results of the simulation study. See the next chapter for
    information on how to combine them;
  - `summary.RDS`, a serialised R object containing the summary
    statistics for each model and data-generating mechanism.
    `summary.RDS` was produced by `make_summaries.R`;
  - `suppl.pdf`, a PDF document with the Supporting Web Material.

### Abstract

<p align="center">

<img src = "./abstract.png">

</p>

### Splitting/Combining Results

In order to upload the full results of the simulation study, we had to
split them in chunks as the single file `Results.RDS` was too large to
upload on GitHub. In this section, we illustrate how the file was split
in chunks and how to re-combine them if needed.

We first imported the single file with results, too large to store on
GitHub:

``` r
# done <- lapply(list.files("Results", pattern = "^out_|^multiout_", full.names = TRUE), readRDS)
# saveRDS(done, "Results.RDS")
done <- readRDS("Results.RDS")
```

As `done` is a list of `data.frame` objects, we combined it into a
single `data.frame`:

``` r
library(dplyr)
done <- dplyr::bind_rows(done)
```

We then split the single `data.frame` into chunks of roughly equal size:

``` r
results.01 <- done[1:ceiling(nrow(done) / 10), ]
results.02 <- done[(ceiling(nrow(done) / 10) + 1):ceiling(nrow(done) / 10 * 2), ]
results.03 <- done[(ceiling(nrow(done) / 10 * 2) + 1):ceiling(nrow(done) / 10 * 3), ]
results.04 <- done[(ceiling(nrow(done) / 10 * 3) + 1):ceiling(nrow(done) / 10 * 4), ]
results.05 <- done[(ceiling(nrow(done) / 10 * 4) + 1):ceiling(nrow(done) / 10 * 5), ]
results.06 <- done[(ceiling(nrow(done) / 10 * 5) + 1):ceiling(nrow(done) / 10 * 6), ]
results.07 <- done[(ceiling(nrow(done) / 10 * 6) + 1):ceiling(nrow(done) / 10 * 7), ]
results.08 <- done[(ceiling(nrow(done) / 10 * 7) + 1):ceiling(nrow(done) / 10 * 8), ]
results.09 <- done[(ceiling(nrow(done) / 10 * 8) + 1):ceiling(nrow(done) / 10 * 9), ]
results.10 <- done[(ceiling(nrow(done) / 10 * 9) + 1):nrow(done), ]
```

We test that by re-combining the chunks we obtain the same `data.frame`:

``` r
done.c <- dplyr::bind_rows(
  results.01,
  results.02,
  results.03,
  results.04,
  results.05,
  results.06,
  results.07,
  results.08,
  results.09,
  results.10
)
all.equal(current = done.c, target = done)
#> [1] TRUE
```

If the previous test passed, we export each chunk as a separate `.RDS`
file:

``` r
if (all.equal(current = done.c, target = done)) {
  saveRDS(results.01, file = "results-01.RDS")
  saveRDS(results.02, file = "results-02.RDS")
  saveRDS(results.03, file = "results-03.RDS")
  saveRDS(results.04, file = "results-04.RDS")
  saveRDS(results.05, file = "results-05.RDS")
  saveRDS(results.06, file = "results-06.RDS")
  saveRDS(results.07, file = "results-07.RDS")
  saveRDS(results.08, file = "results-08.RDS")
  saveRDS(results.09, file = "results-09.RDS")
  saveRDS(results.10, file = "results-10.RDS")
}
```

To recombine the chunks uploaded on GitHub, first download all the files
in a single directory. Then, run the following `R` code making sure that
the working directory is the directory where the files are saved:

``` r
done.r <- lapply(1:10, function(n) readRDS(file = paste0("results-", sprintf("%02d", n), ".RDS")))
done.r <- do.call(rbind.data.frame, done.r)
```

Test that we obtain the same results:

``` r
all.equal(current = done.r, target = done)
#> [1] TRUE
```
