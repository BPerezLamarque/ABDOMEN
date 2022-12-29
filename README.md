# ABDOMEN (A Brownian moDel Of Microbiota EvolutioN):





This tutorial explains how to use ABDOMEN (A Brownian moDel Of Microbiota EvolutioN) for a given host-microbiota system (i.e. a fixed, bifurcating host phylogeny and microbiota relative abundances for each extant host species; Fig. 1a). This phylogenetic comparative approach apply the multivariate Brownian motion process to compositional data. It also includes the widely-used Pagel's $\lambda$ tree transformation that quantifies phylosymbiosis by evaluating how much host phylogeny contributes to explaining interspecific variation in present-day microbiota composition (Fig. 1b). ABDOMEN assumes that, from ancestral values at the root $X_0$, the log-absolute abundances of the different microbial taxa change on the host phylogeny following a multivariate Brownian motion model with variance-covariance matrix $R$ (Fig. 1a). 

<p align="center">
    <img src="https://github.com/BPerezLamarque/ABDOMEN/blob/main/example/ABDOMEN.jpg" width="300">
</p>

<p align="center">
    <b>ABDOMEN: A comparative phylogenetic model for the dynamics of microbiota composition during host diversification.</b>
</p>


Noting $X_{ij}$ the unmeasured absolute abundance of microbial taxon $j$ in host $i$, $Y_i=\sum_j X_{ij}$ the unmeasured total microbial abundance in the microbiota of host $i$, $\tilde{Y}_i = Y_i Y_0$ its value relative to the unknown total microbial abundance at the root $Y_0$, and $Z_{ij}=X_{ij}/Y_i$  the measured relative abundance of microbial taxon $j$ in host $i$, we sample from the joint posterior distribution $P(\log⁡Z_0, R, \lambda, \log⁡\tilde{Y}̃_1,...,\log⁡\tilde{Y}_n | Z_{11},…,Z_{ij},…,Z_{np},C), where $Z_0$ is the vector of relative abundances at the root, $n$ is the number of host species, $p$ is the number of microbial taxa, and $C$ is the phylogenetic variance-covariance matrix. ABDOMEN perform this sampling using a No U-turn Hamiltonian Markov Chain Monte Carlo algorithm implemented via the Stan probabilistic programming language. 

ABDOMEN thus outputs  an estimate of phylosymbiosis (measured as Pagel’s $\lambda$), microbiota integration (reflected in the $R$ matrix), and ancestral microbiota composition (given by $Z_0$). 



**Citation:** 


Benoît Perez-Lamarque, Guilhem Sommeria-Klein, Loréna Duret, and Hélène Morlon
,
*Phylogenetic comparative approach reveals evolutionary conservatism, ancestral composition, and integration of vertebrate gut microbiota*, bioRxiv, [doi: https://doi.org/10.1101/XXXXX](https://www.biorxiv.org/).


**Contact:** Benoît Perez-Lamarque (benoit.perez.lamarque@gmail.com) and Guilhem Sommeria-Klein (guilhem.sk@gmail.com)




# Contents:
**[Installation](#installation)**\
**[Running ABDOMEN](#running-abdomen)**\
**[Plotting outputs](#plotting-outputs)**\


# Installation:


The scripts to run ABDOMEN ("ABDOMEN.R" and "ABDOMEN_v1.0.stan") must be downloaded from the folder ['script'](https://github.com/BPerezLamarque/ABDOMEN/tree/main/script/). 
In addition, the following packages must be installed:

```r
install.packages("ggplot2", "mvMORPH", "RPANDA", "rstan")

```



# Running ABDOMEN:



First, you can set your working directory and load the packages and functions:

```r

setwd("YOUR_WORKING_DIRECTORY")

library(ggplot2)
library(mvMORPH)
library(RPANDA)
library(rstan)
rstan_options(auto_write = TRUE)


source("ABDOMEN.R")

```

Then, you can load the example dataset of the gut microbiota of Certardiodactyla (from Song et al., 2020) that can be downloaded from the folder ['example'](https://github.com/BPerezLamarque/ABDOMEN/tree/main/example/).: 

```r

# Open the phylogenetic tree of the Cetartiodactyla
tree <- read.tree("tree_Cetartiodactyla.tre")

# Open the table with the mean abundances of the main bacterial orders in the gut microbiota of each Cetartiodactyla species
table <- t(read.table("table_bacterial_orders_Cetartiodactyla.txt", header=TRUE, sep="\t"))
# each row corresponds to a Cetartiodactyla species, and each column to a bacterial order

```


The following parameters must be specified to run ABDOMEN:

```r

name="run_Cetartiodactyla_bacterial_orders" # the name of the run

detection_threshold <- 1e-05 # the detection threshold: below this threshold, we assume that we cannot detect a given microbial taxa. Then, all relative abundances below this threshold are set to this threshold. 

seed <- 3 # seed for reproductibility

mean_prior_logY <- 0 # mean value for the Gaussian prior of logY (the latent variable that correspond to the total microbial abundances, relative to the ancestral ones)
sd_prior_logY <- 2  # standart deviation for the Gaussian prior of logY (the latent variable that correspond to the total microbial abundances, relative to the ancestral ones)

nb_cores <- 4 # number of cores to run the analyses
chains <-  4 # number of chains for the inference
warmup <-  500 # number of warmup iterations in STAN
iter <-  1000 # total number of iterations in STAN

```

Now you can run ABDOMEN: 

```r

fit_summary <- ABDOMEN(tree, table, name, nb_cores = 2, chains = 4, warmup = 1000, iter = 2000)

```

# Plotting outputs:

You can now automatically plot the estimated variances and covariances ($R$), the ancestral microbiota composition ($Z_0$), as well as the ancestral microbiota composition in each node of the phylogeny:

```r

ABDOMEN_process_output(tree, table, name, fit_summary)

```
