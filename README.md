# genomic Architecture

Code for Ashraf and Lawson 2020 (to appear).

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
