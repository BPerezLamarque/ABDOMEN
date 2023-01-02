# ABDOMEN (A Brownian moDel Of Microbiota EvolutioN):





This tutorial explains how to use ABDOMEN (A Brownian moDel Of Microbiota EvolutioN) for a given host-microbiota system, *i.e.* a fixed, bifurcating host phylogeny and microbiota relative abundances for each extant host species (Fig. 1a). 

This phylogenetic comparative approach apply the multivariate Brownian motion process to compositional data. It also includes the widely-used Pagel's $\lambda$ tree transformation that quantifies phylosymbiosis by evaluating how much host phylogeny contributes to explaining interspecific variation in present-day microbiota composition (Fig. 1b). ABDOMEN assumes that, from ancestral values at the root $X_0$, the log-absolute abundances of the different microbial taxa change on the host phylogeny following a multivariate Brownian motion model with variance-covariance matrix $R$ (Fig. 1a). 

<p align="center">
    <img src="https://github.com/BPerezLamarque/ABDOMEN/blob/main/example/ABDOMEN.png" width="500">
</p>

<p align="center">
    <b>Figure 1: ABDOMEN: A comparative phylogenetic model for the dynamics of microbiota composition during host diversification.</b>
</p>

Noting $X_{ij}$ the unmeasured absolute abundance of microbial taxon $j$ in host $i$, $Y_i=\sum_j X_{ij}$ the unmeasured total microbial abundance in the microbiota of host $i$, ${\tilde{Y_i}} = Y_i / Y_0$ its value relative to the unknown total microbial abundance at the root $Y_0$, and $Z_{ij}=X_{ij}/Y_i$ the measured relative abundance of microbial taxon $j$ in host $i$, we sample from the joint posterior distribution $P(\log Z_0, R, \lambda, \log {\tilde{Y_1}},...,\log {\tilde{Y_n}} | Z_{11},…,Z_{ij},…,Z_{np},C)$, where $Z_0$ is the vector of relative abundances at the root, $n$ is the number of host species, $p$ is the number of microbial taxa, and $C$ is the phylogenetic variance-covariance matrix. ABDOMEN performs this sampling using a No U-turn Hamiltonian Markov Chain Monte Carlo algorithm implemented via the Stan probabilistic programming language.

ABDOMEN thus outputs  an estimate of phylosymbiosis (measured as Pagel’s $\lambda$), microbiota integration (reflected in the $R$ matrix), and ancestral microbiota composition (given by $Z_0$). 



**Citation:** Benoît Perez-Lamarque, Guilhem Sommeria-Klein, Loréna Duret, and Hélène Morlon
,
*Phylogenetic comparative approach reveals evolutionary conservatism, ancestral composition, and integration of vertebrate gut microbiota*, bioRxiv, [doi: https://doi.org/10.1101/XXXXX](https://www.biorxiv.org/).


**Contact:** Benoît Perez-Lamarque (benoit.perez.lamarque@gmail.com) and Guilhem Sommeria-Klein (guilhem.sk@gmail.com)




# Contents:
**[Installation](#installation)**\
**[Running ABDOMEN](#running-abdomen)**\
**[Plotting outputs](#plotting-outputs)**


# Installation:


The scripts to run ABDOMEN ("ABDOMEN.R" and "ABDOMEN_v1.0.stan") must be downloaded from the folder ['script'](https://github.com/BPerezLamarque/ABDOMEN/tree/main/script/). 
In addition, the following packages must be installed:

```r
install.packages("ggplot2", "mvMORPH", "RPANDA", "rstan", "RColorBrewer")

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
library(RColorBrewer)


source("ABDOMEN.R")

```

Then, you can load the example dataset of the gut microbiota of Certardiodactyla (from Song et al., 2020) that can be downloaded from the folder ['example'](https://github.com/BPerezLamarque/ABDOMEN/tree/main/example/): 

```r

# Open the phylogenetic tree of the Cetartiodactyla
tree <- read.tree("tree_Cetartiodactyla.tre")

# Open the table with the mean abundances of the main bacterial orders in the gut microbiota of each Cetartiodactyla species
table <- t(read.table("table_bacterial_orders_Cetartiodactyla.txt", header=TRUE, sep="\t"))
# each row corresponds to a Cetartiodactyla species, and each column to a bacterial order (there are a total of 14 bacterial orders)

```

NBB: The Certardiodactyla dataset is only composed on 33 host species, which is not enought for having very robust and accurate parameters estimates (ideally, we recommend >50 host species and >5 microbial taxa), but ABDOMEN can be run quickly of this small dataset for demonstration. 


The following parameters must be specified to run ABDOMEN:

```r

name <- "run_Cetartiodactyla_bacterial_orders" # the name of the run

code_path <- getwd() # indicates where the stan codes are stored (here, there are directly stored in the working directory)

detection_threshold <- 1e-05 # the detection threshold: below this threshold, we assume that we cannot detect a given microbial taxa. Then, all relative abundances below this threshold are set to this threshold. 

seed <- 3 # seed for reproductibility

mean_prior_logY <- 0 # mean value for the Gaussian prior of logY (the latent variable that correspond to the total microbial abundances, relative to the ancestral ones)
sd_prior_logY <- 2  # standart deviation for the Gaussian prior of logY (the latent variable that correspond to the total microbial abundances, relative to the ancestral ones)

nb_cores <- 4 # number of cores to run the analyses
chains <-  4 # number of chains for the inference
warmup <-  1000 # number of warmup iterations in STAN
iter <-  2000 # total number of iterations in STAN

```

Now you can run ABDOMEN: 

```r

fit_summary <- ABDOMEN(tree, table, name, 
                       code_path = code_path,
                       detection_threshold = detection_threshold, seed = seed, 
                       mean_prior_logY = mean_prior_logY, sd_prior_logY = sd_prior_logY,
                       nb_cores = nb_cores, chains = chains, warmup = warmup, iter = iter)

```

# Plotting outputs:

You can now automatically plot the estimated variances and covariances ($R$), the ancestral microbiota composition ($Z_0$), as well as the ancestral microbiota composition in each node of the phylogeny:

```r

ABDOMEN_process_output(tree, table, name, fit_summary)

```


<p align="center">
    <img src="https://github.com/BPerezLamarque/ABDOMEN/blob/main/example/plot_ABDOMEN/convergence_chains_run_Cetartiodactyla_bacterial_orders_lambda.png" width="350">
</p>
<p align="center">
    <b>Figure 2: Convergence of the 4 chains of ABDOMEN for the estimation of Pagel's $\lambda$, $i.e.$ the measure of phylosymbiosis. Here, $\lambda$ is estimated at 0.23 (95% CI: [0.1; 0.36]).</b>
</p>


<p align="center">
    <img src="https://github.com/BPerezLamarque/ABDOMEN/blob/main/example/plot_ABDOMEN/results_Z0_run_Cetartiodactyla_bacterial_orders.png" width="350">
</p>
<p align="center">
    <b>Figure 3: Estimated ancestral microbiota composition of Cetartiodactyla (mean of the posterior distribution at the root).</b>
</p>


<p align="center">
    <img src="https://github.com/BPerezLamarque/ABDOMEN/blob/main/example/plot_ABDOMEN/results_Z0_run_Cetartiodactyla_bacterial_orders.png" width="350">
</p>
<p align="center">
    <b>Figure 4: Estimated ancestral microbiota composition of Cetartiodactyla (mean of the posterior distribution at the root and generalized least squares estimates at internal nodes).</b>
</p>


<p align="center">
    <img src="https://github.com/BPerezLamarque/ABDOMEN/blob/main/example/plot_ABDOMEN/results_variances_run_Cetartiodactyla_bacterial_orders.png" width="350">
</p>
<p align="center">
    <b>Figure 5: Estimated variances for each bacterial order in the gut microbiota of Cetartiodactyla (mean of the posterior distribution).</b>
</p>

<p align="center">
    <img src="https://github.com/BPerezLamarque/ABDOMEN/blob/main/example/plot_ABDOMEN/results_variances_run_Cetartiodactyla_bacterial_orders.png" width="350">
</p>
<p align="center">
    <b>Figure 6: Estimated covariances between bacterial order in the gut microbiota of Cetartiodactyla (mean of the posterior distribution) All covariances are represented in (a), while only significant ones are in (b; a significant covariance means that 0 is not in its 95% CI - many covariances are not significant here because the dataset is too small).</b>
</p>




# Assessing the significance of phylosymbiosis:

To assess the significance of phylosymbiosis, one can perform permutatations. Here, we will randomly permutate all the host species:

```r

name_random <- "run_Cetartiodactyla_bacterial_orders_permutation_1"

seed <- 100
set.seed(seed)

table_random <- table[sample(tree$tip.label),] # randomly permutates all the Cetartiodactyla species:

fit_summary <- ABDOMEN(tree, table_random, name_random, 
                       code_path = code_path,
                       detection_threshold = detection_threshold, seed = seed, 
                       mean_prior_logY = mean_prior_logY, sd_prior_logY = sd_prior_logY,
                       nb_cores = nb_cores, chains = chains, warmup = warmup, iter = iter)

ABDOMEN_process_output(tree, table_random, name, fit_summary)

# this step must we replicated a large number of times (e.g. 100 with different seeds) to compare the original lambda values with the ones obtained when permuttatting the dataset. 

```



NB: More constraint permutations can be done, *e.g.* only permuting species having a similar diet to test the effect of diet conservatism on phylosymbiosis. 


