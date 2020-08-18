# Genomic Architecture

This is code for [Ashraf and Lawson 2020](https://www.biorxiv.org/content/10.1101/2020.08.17.254110v1) (biorxiv link) "Genetic drift from the out-of-Africa bottleneck leads to biased estimation of genetic architecture and selection".

It implements a model of genetic drift that highlights some issues in how Genomic Architecture is understood in the context of selection.

## Installation of RStan

You need to install [Rstan](https://mc-stan.org/users/interfaces/rstan). Short version:

```R
install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)
```

Note that the C++ Toolchain needs to be installed, for which you may need to follow the advice in the [RStan Getting Started Guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

## Model compilation

Models need only be compiled once to be run with any input parameters. To compile, at the shell script, type:

```sh
Rscript sim_bayesinf.R compile
```

If you are on Windows, you may simply wish to load the script and run the corresponding lines of code from [sim_bayesinf.R](sim_bayesinf.R).

## Running the inference

The model comes with two options for how data are generated. Both of these use a function called `make_test_data`, which in the simple model simulates from the inference model.

To run inference using the simple model, simply run:

```sh
Rscript sim_bayesinf.R <S> <Fst> <seed>
```
e.g.
```sh
Rscript sim_bayesinf.R -1 0.1 1
```
The results will be saved to an RData object called `paste0(s_',S,'_Fst',Fst,'_seed',seed,'.RData')`

To run the inference using the 1000G frequencies, run:

```sh
Rscript sim_bayesinf_1kg.R <S> <pop> <Fst> <seed>
```
where `<pop>` is the name of the 1000G population of interest (afr,eur,sas,eas, or amr). The results will be saved to an RData object called `paste0('s_',S,'_target',target,'_Fst',Fst,'_seed',seed,'_1kg.RData')`.

## Notes on run time

We used 10000 MCMC iterations and 10000 SNPs in our paper, and the scripts are set to these values. These take over a day to run and so are not suitable for experimentation.

The model likelihood evaluation is linear in the number of SNPs, but the convergence time grows with this number, and Stan uses a gradient-based optimiser that is quadratic in the number of parameters. Therefore, you can obtain valid experimentation using say 10k iterations and 1k SNPs and this takes under a minute. However, many SNPs are needed to obtain a stable distribution and so you may need more SNPs to accurately replicate our infererences.

## Thoughts on practical usage

This is a development model and is not designed for serious deployment without care; it is trivial to introduce uncertainties in the true values of $\beta_i$ as well as weights for each SNP. However, moving into a dedicated MCMC framework may be necessary to obtain fast convergence for larger datasets.

## Contact

Daniel Lawson [dan.lawson@bristol.ac.uk](mailto:dan.lawson@bristol.ac.uk)
