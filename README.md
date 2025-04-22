# ABDOMEN (A Brownian moDel Of Microbiota EvolutioN):





This tutorial explains how to use ABDOMEN (A Brownian moDel Of Microbiota EvolutioN) using R for a given host-microbiota system, *i.e.* a fixed, bifurcating, calibrated host phylogeny and microbiota relative abundances for each extant host species typically measured using metabarcoding techniques (Figure 1a). 

This phylogenetic comparative approach applies the multivariate Brownian motion process to compositional data. It also includes the widely-used Pagel's $\lambda$ tree transformation that quantifies phylosymbiosis by evaluating how much host phylogeny contributes to explaining interspecific variation in present-day microbiota composition (Figure 1b). ABDOMEN assumes that, from ancestral values at the root $X_0$, the log-absolute abundances of the different microbial taxa change on the host phylogeny following a multivariate Brownian motion model with a variance-covariance matrix $R$ (Figure 1a). 

<p align="center">
    <img src="https://github.com/BPerezLamarque/ABDOMEN/blob/main/example/ABDOMEN.png" width="600">
</p>

<p align="center">
    <b>Figure 1: ABDOMEN: A comparative phylogenetic model for the dynamics of microbiota composition during host diversification.</b>
</p>

Noting $X_{ij}$ the unmeasured absolute abundance of microbial taxon $j$ in host $i$, $Y_i=\sum_j X_{ij}$ the unmeasured total microbial abundance in the microbiota of host $i$, ${\tilde{Y_i}} = Y_i / Y_0$ its value relative to the unknown total microbial abundance at the root $Y_0$, and $Z_{ij}=X_{ij}/Y_i$ the measured relative abundance of microbial taxon $j$ in host $i$, we sample from the joint posterior distribution $P(\log Z_0, R, \lambda, \log {\tilde{Y_1}},...,\log {\tilde{Y_n}} | Z_{11},…,Z_{ij},…,Z_{np},C)$, where $Z_0$ is the vector of relative abundances at the root, $n$ is the number of host species, $p$ is the number of microbial taxa, and $C$ is the phylogenetic variance-covariance matrix. ABDOMEN performs this sampling using a No U-turn Hamiltonian Markov Chain Monte Carlo algorithm implemented via the Stan probabilistic programming language.

ABDOMEN thus outputs  an estimate of phylosymbiosis (measured as Pagel’s $\lambda$), microbiota integration (covariances between microbial taxa reflected in the $R$ matrix), and ancestral microbiota composition (given by $Z_0$). 



**Citation:** Benoît Perez-Lamarque, Guilhem Sommeria-Klein, Loréna Duret, and Hélène Morlon, *Phylogenetic comparative approach reveals evolutionary conservatism, ancestral composition, and integration of vertebrate gut microbiota*, Molecular Biology and Evolution, Volume 40, Issue 7, July 2023, msad144, [doi: https://doi.org/10.1093/molbev/msad144](https://doi.org/10.1093/molbev/msad144).


**Contact:** Benoît Perez-Lamarque (benoit.perez.lamarque@gmail.com) and Guilhem Sommeria-Klein (guilhem.sk@gmail.com)




# Contents:
**[Installation](#installation)**\
**[Running ABDOMEN](#running-abdomen)**\
**[Plotting outputs](#plotting-outputs)**\
**[Assessing the significance of phylosymbiosis](#assessing-the-significance-of-phylosymbiosis)**




# Installation:


The R and Stan scripts to run ABDOMEN ("ABDOMEN.R" and "ABDOMEN_v1.0.stan") must be downloaded from the folder ['script'](https://github.com/BPerezLamarque/ABDOMEN/tree/main/script/) and stored in your R working directory. 
In addition, the following packages must be installed:

```r
install.packages("ggplot2", "mvMORPH", "RPANDA", "rstan", "RColorBrewer")

```



# Running ABDOMEN:



First, you can set your working directory and **load the packages and functions**:

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

Then, you can **load the example dataset** of the gut microbiota of Certardiodactyla (from Song et al., 2020) that can be downloaded from the folder ['example'](https://github.com/BPerezLamarque/ABDOMEN/tree/main/example/): 

```r

# Open the phylogenetic tree of the Cetartiodactyla
tree <- read.tree("tree_Cetartiodactyla.tre")
# Please note that the phylogenetic tree must be rooted, binary, and ultrametric to run ABDOMEN. 

# Open the table with the mean abundances of the main bacterial orders in the gut microbiota of each Cetartiodactyla species
table <- t(read.table("table_bacterial_orders_Cetartiodactyla.txt", header=TRUE, sep="\t"))
# each row corresponds to a Cetartiodactyla species and each column to a bacterial order (there are a total of 14 bacterial orders)

```

NB: The Certardiodactyla dataset is only composed of 33 host species, which is not enough for having very robust and accurate parameters estimates (ideally, we recommend >50 host species and >5 microbial taxa), but ABDOMEN can be run quickly on this small dataset for demonstration. 


The following **parameters must be specified** to run ABDOMEN:

```r

name <- "run_Cetartiodactyla_bacterial_orders" # the name of the run

code_path <- getwd() # indicates where the stan codes are stored (here, there are directly stored in the working directory) and where the ABDOMEN plots will be generated.

detection_threshold <- 1e-05 # the detection threshold: below this threshold, we assume that we cannot detect the microbial taxa (either because it is not present or because we cannot detect very rare taxa with metabarcoding techniques). Then, all relative abundances below this threshold are set to this threshold. 

seed <- 3 # seed for reproductibility

mean_prior_logY <- 0 # mean value for the Gaussian prior of logY (the latent variables that correspond to the total microbial abundances, relative to the ancestral ones)
sd_prior_logY <- 2  # standard deviation for the Gaussian prior of logY (the latent variables that correspond to the total microbial abundances, relative to the ancestral ones)

nb_cores <- 4 # number of cores to run the analyses
chains <-  4 # number of chains for the inference
warmup <-  1000 # number of warmup iterations in STAN
iter <-  2000 # total number of iterations in STAN

```

Now you can **run ABDOMEN**: 

```r

fit_summary <- ABDOMEN(tree, table, name, 
                       code_path = code_path,
                       detection_threshold = detection_threshold, seed = seed, 
                       mean_prior_logY = mean_prior_logY, sd_prior_logY = sd_prior_logY,
                       nb_cores = nb_cores, chains = chains, warmup = warmup, iter = iter)

```

# Plotting outputs:

You can now automatically check the convergence of the chains and plot the estimated variances and covariances between microbial taxa ($R$), the ancestral microbiota composition ($Z_0$), as well as the ancestral microbiota composition in each node of the phylogeny:

All plots are saved as PDF files in the 'plot_ABDOMEN/' folder, which is automatically created at the location specified by the variable ``code_path`` (default: the working directory). To customize the appearance of the plots, you can directly modify the script of the function [`ABDOMEN_process_output`](https://github.com/BPerezLamarque/ABDOMEN/blob/main/script/ABDOMEN.R) function. The colors can be specified using the argument `list_colors`.


```r

# Specify the color for each microbial taxon:
list_colors <- scales::hue_pal()(ncol(table))
names(list_colors) <- colnames(table)


# Generate the PDF plots:    
ABDOMEN_process_output(tree, table, name, fit_summary, code_path = code_path, 
		detection_threshold = detection_threshold, list_colors=list_colors)

```


<p align="center">
    <img src="https://github.com/BPerezLamarque/ABDOMEN/blob/main/example/plot_ABDOMEN/convergence_chains_run_Cetartiodactyla_bacterial_orders_lambda.png" width="350">
</p>
<p align="center">
    <b>Figure 2: Convergence of the 4 chains of ABDOMEN for the estimation of Pagel's lambda, i.e. the measure of phylosymbiosis. Here, lambda is estimated at 0.23 (95% CI: [0.1; 0.36]).</b>
</p>


Besides the plots, you can also directly extract the estimated parameter values: 

```r

# Pagel's lambda (the measure of phylosymbiosis):
original_lambda <- ABDOMEN_extract_lambda(tree, table, fit_summary, detection_threshold = detection_threshold) # gives the mean estimated Pagel's lambda and its 95% CI
original_lambda

# The ancestral microbiota composition (Z0):
ABDOMEN_extract_Z0(tree, table, fit_summary, detection_threshold = detection_threshold) #  gives the mean ancestral microbiota composition and the 95% CI of each microbial taxa

# The ancestral microbiota composition estimated for each node of the host phylogenetic tree  (Z0_nodes):
Z0_nodes <- ABDOMEN_extract_Z0_nodes(tree, table, fit_summary, detection_threshold = detection_threshold) # Z0_nodes contains the relative abundances of each bacterial taxa, the node number in the host tree, and the list of all extant taxa descending from this node (column MRCA).


```


<p align="center">
    <img src="https://github.com/BPerezLamarque/ABDOMEN/blob/main/example/plot_ABDOMEN/results_Z0_run_Cetartiodactyla_bacterial_orders.png" width="350">
</p>
<p align="center">
    <b>Figure 3: Estimated ancestral microbiota composition of Cetartiodactyla (Z0; mean of the posterior distribution at the root).</b>
</p>


<p align="center">
    <img src="https://github.com/BPerezLamarque/ABDOMEN/blob/main/example/plot_ABDOMEN/results_ancestral_states_run_Cetartiodactyla_bacterial_orders.png" width="350">
</p>
<p align="center">
    <b>Figure 4: Estimated ancestral microbiota composition of Cetartiodactyla (mean of the posterior distribution at the root (Z0) and generalized least squares estimates at internal nodes).</b>
</p>

```r
# The microbial integration (R; variance-covariance matrix between microbial taxa):
R_matrices <- ABDOMEN_extract_R(tree, table, fit_summary, detection_threshold = detection_threshold) 
R_matrices$R #  gives the mean variance-covariance matrix R between microbial taxa
R_matrices$R_lower_bound #  gives the lower bound of the 95% CI of the variance-covariance matrix R 
R_matrices$R_upper_bound #  gives the upper bound of the 95% CI of the variance-covariance matrix R 
R_matrices$R_signif #  only outputs the significant covariances (i.e. covariances for which 0 is not in the 95% CI)

```

<p align="center">
    <img src="https://github.com/BPerezLamarque/ABDOMEN/blob/main/example/plot_ABDOMEN/results_variances_run_Cetartiodactyla_bacterial_orders.png" width="350">
</p>
<p align="center">
    <b>Figure 5: Estimated variances for each bacterial order in the gut microbiota of Cetartiodactyla (diagonal elements of R; mean of the posterior distribution).</b>
</p>

<p align="center">
    <img src="https://github.com/BPerezLamarque/ABDOMEN/blob/main/example/plot_ABDOMEN/results_covariances_run_Cetartiodactyla_bacterial_orders.png" width="350">
</p>
<p align="center">
    <img src="https://github.com/BPerezLamarque/ABDOMEN/blob/main/example/plot_ABDOMEN/results_covariances_run_Cetartiodactyla_bacterial_orders_signif.png" width="350">
</p>
<p align="center">
    <b>Figure 6: Estimated covariances between bacterial order in the gut microbiota of Cetartiodactyla (off-diagonal elements of R; mean of the posterior distribution). All covariances are represented in (a), while only significant ones are in (b). A significant covariance means that 0 is not in its 95% CI - many covariances are not significant here because the dataset is too small.</b>
</p>




# Assessing the significance of phylosymbiosis:

To assess the significance of phylosymbiosis, one can perform permutations. Here, we permutate all the host species at random:

```r

name_random <- "run_Cetartiodactyla_bacterial_orders_permutation_1"

seed <- 100
set.seed(seed)

table_random <- table[sample(tree$tip.label),] # randomly permutates all the Cetartiodactyla species:
rownames(table_random) <- rownames(table)

fit_summary_permut <- ABDOMEN(tree, table_random, name = name_random, 
                       code_path = code_path,
                       detection_threshold = detection_threshold, seed = seed, 
                       mean_prior_logY = mean_prior_logY, sd_prior_logY = sd_prior_logY,
                       nb_cores = nb_cores, chains = chains, warmup = warmup, iter = iter)

ABDOMEN_extract_lambda(tree, table_random, fit_summary_permut, detection_threshold = detection_threshold)

```

This step must be replicated a large number of times (e.g. 100 with different seeds) to compare the original lambda values with the ones obtained when permutating the dataset and obtain a p-value. For instance :

```r

nb_permutations <- 100
list_lambda_permutations <- c()
for (seed in 1:nb_permutations){
    set.seed(seed)
    name_random <- paste0("run_Cetartiodactyla_bacterial_orders_permutation_",seed)
    table_random <- table[sample(tree$tip.label),] # randomly permutates all the Cetartiodactyla species:
    rownames(table_random) <- rownames(table)
    fit_summary_permut <- ABDOMEN(tree, table_random, name = name_random, 
                       code_path = code_path,
                       detection_threshold = detection_threshold, seed = seed, 
                       mean_prior_logY = mean_prior_logY, sd_prior_logY = sd_prior_logY,
                       nb_cores = nb_cores, chains = chains, warmup = warmup, iter = iter)
    list_lambda_permutations <- rbind(list_lambda_permutations, ABDOMEN_extract_lambda(tree, table_random, fit_summary_permut, detection_threshold = detection_threshold))
}

# p-value for the significance of phylosymbiosis:
length(which(list_lambda_permutations[,1]>=original_lambda[1]))/nb_permutations 

```



In addition, more constrained permutations can be done, *e.g.* only permuting species having a similar diet or the same geographic location to test the effect of diet conservatism or biogeographic conservatism on phylosymbiosis:


```r

metadata <- read.table("metadata_Cetartiodactyla.csv", header=TRUE, sep=";")
colnames(metadata) # in this example, we will randomize the species based on their geographic location (indicated in the column "location").
random = "location" # name of the variable used for the constrained permutations. It can be done using any kind of qualitative variable. 

nb_permutations <- 100
list_lambda_permutations <- c()
for (seed in 1:nb_permutations){
    set.seed(seed)
    name_random <- paste0("run_Cetartiodactyla_bacterial_orders_permutation_diet_",seed)
    
    # randomly permutate the Cetartiodactyla species having the same geographic location 
    categories <- sort(unique(metadata[,random]))
    metadata$species_random <- metadata$species
    for (cat in categories){
      list_index <- metadata$species_random[which(metadata[,random]==cat)]
      metadata$species_random[which(metadata[,random]==cat)] <- list_index[sample(length(list_index))]
    }
    table_random <- table[metadata$species_random,]
    rownames(table_random) <- metadata$species
    
    fit_summary_permut <- ABDOMEN(tree, table_random, name = name_random, 
                       code_path = code_path,
                       detection_threshold = detection_threshold, seed = seed, 
                       mean_prior_logY = mean_prior_logY, sd_prior_logY = sd_prior_logY,
                       nb_cores = nb_cores, chains = chains, warmup = warmup, iter = iter)
    list_lambda_permutations <- rbind(list_lambda_permutations, ABDOMEN_extract_lambda(tree, table_random, fit_summary_permut, detection_threshold = detection_threshold))
}

# p-value for the significance of phylosymbiosis:
length(which(list_lambda_permutations[,1]>=original_lambda[1]))/nb_permutations 


```