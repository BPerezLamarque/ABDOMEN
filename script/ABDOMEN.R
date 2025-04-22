ABDOMEN <- function(tree, table, name, code_path = getwd(), detection_threshold=1e-05, seed=3, mean_prior_logY=0, sd_prior_logY=2,
                    nb_cores = 1, chains = 4, warmup = 500, iter = 1000, plot_chains=TRUE){
  
  # scale the tree
  tree$edge.length <- tree$edge.length/max(node.depth.edgelength(tree)) 
  
  if (!is.rooted(tree)){print("WARNING: Your tree is not rooted. Please root your tree before running ABDOMEN.")}
  if (!is.ultrametric(tree)){print("WARNING: Your tree is not ultrametric. Please calibrate your tree before running ABDOMEN. If and only if, the tree is non ultrametric because of numerical precisions, you can use the function 'force.ultrametric()'.")}
  if (!all(tree$tip.label %in% rownames(table))){print("WARNING: Some species in the tree are not in the OTU table. Please remove these species from the phylogenetic tree before running AB.")}
  
  # scale the abundance per row (each row sum must be equal to 1)
  for (i in 1:nrow(table)) {table[i,] <- table[i,]/sum(table[i,])}
  
  # reorder table as tree$tip.label
  table <- table[tree$tip.label,]
  
  while (length(table[which(table<detection_threshold)])>0){
    table[table<detection_threshold] <- detection_threshold
    for(i in 1:nrow(table)) {table[i,] <- table[i,]/sum(table[i,])}
  }
  
  
  # Other parameters
  n=nrow(table)
  p=ncol(table)
  C=vcv(tree)
  
  
  # Input data
  microbiotree_input_data = list(n = n, p = p, 
                                 logZ = log(table), 
                                 C = C,
                                 mean_prior_logY = mean_prior_logY,
                                 sd_prior_logY = sd_prior_logY)
  
  # Initiation 
  
  init_function = function(n, p, Posdef, tree, mvSIM){
    
    Z0 = runif(p, 0, 1)
    Z0 = Z0/sum(Z0)
    
    R = Posdef(p, rexp(p, 1/2))    
    
    logX=mvSIM(tree = tree, nsim=1, model="BM1", param = list(ntraits=p, sigma=R,theta=log(Z0)))
    logY=log(rowSums(exp(logX)))
    
    lambda = runif(1, 0, 1)
    
    return = list(Z0 = Z0, logY = logY, R = R, lambda=lambda)
  }
  
  # Prepare the working directory
  dir.create(file.path(code_path, "plot_ABDOMEN/"), showWarnings = FALSE)
  
  fit = stan(file=paste0(code_path,"/ABDOMEN_v1.0.stan"),
             data = microbiotree_input_data, cores=nb_cores,
             include = TRUE, pars=c("Z0","logY","logY_trans", "R", "lambda"),
             chains = chains, iter = iter, warmup = warmup, thin = 1, 
             init = function() init_function(n, p, Posdef, tree, mvSIM), seed = seed,
             control = list(max_treedepth = 15, adapt_delta = 0.99) )
  ### return an array of three dimensions: iterations, chains, parameters 
  
  
  # Process the output of the run
  
  if (plot_chains){
    
    a <- extract(fit, permuted = FALSE) 
    a_warmup <- extract(fit, permuted = FALSE, inc_warmup = TRUE)
    
    list_parameters <- c("lambda", paste0("Z0[",1:p,"]"), sort(apply(expand.grid(paste0("R[",1:p,","), paste0(1:p,"]")), 1, paste, collapse="")), paste0("logY[",1:n,"]"))
    
    # Check the convergence of the chain
    list_Rhat <- c()
    list_ESS_bulk <- c()
    list_ESS_tail <- c()
    
    pdf(paste0("plot_ABDOMEN/convergence_chains_",name,".pdf"),width=4.5, height=4)
    for (param in  list_parameters){
      min <- min(as.vector(c(a[,,param])))
      max <- max(as.vector(c(a[,,param])))
      
      Rhat <- Rhat(a_warmup[,,param])
      ESS_bulk <- ess_bulk(a_warmup[,,param]) 
      ESS_tail <- ess_tail(a_warmup[,,param])
      list_Rhat <- c(list_Rhat, Rhat)
      list_ESS_bulk <- c(list_ESS_bulk, ESS_bulk)
      list_ESS_tail <- c(list_ESS_tail, ESS_tail)
      
      thin=5 # only for plotting 
      chain=1
      (plot(a[seq(1,(iter-warmup)/2,thin),chain,param], type="l", ylim=c(min, max), ylab = param, xlab=paste0("thin=",thin*2), col=chain, main=paste0("Rhat=",round(Rhat,3), ", ESS=", round(ESS_bulk))))
      for (chain in 2:chains){
        par(new=TRUE)
        (plot(a[seq(1,(iter-warmup)/2,thin),chain,param], type="l", ylim=c(min, max), col=chain, ylab = "", xlab="", axes = F))
      }
    }
    dev.off()
  }
  
  fit_summary <- summary(fit)
  
  return(fit_summary)
  
}


ABDOMEN_process_output <- function(tree, table, name, fit_summary, code_path = getwd(),
                                   detection_threshold=1e-05, list_colors=NULL){
  
  # scale the tree
  tree$edge.length <- tree$edge.length/max(node.depth.edgelength(tree)) 
  
  # scale the abundance per row (each row sum must be equal to 1)
  for (i in 1:nrow(table)) {table[i,] <- table[i,]/sum(table[i,])}
  
  # reorder table as tree$tip.label
  table <- table[tree$tip.label,]
  
  while (length(table[which(table<detection_threshold)])>0){
    table[table<detection_threshold] <- detection_threshold
    for(i in 1:nrow(table)) {table[i,] <- table[i,]/sum(table[i,])}
  }
  
  
  # Determine the color for each microbial taxon
  if (is.null(list_colors)){
    list_colors <- scales::hue_pal()(ncol(table))
    names(list_colors) <- colnames(table)
  }
  
  
  
  p <- ncol(table)
  n <- nrow(table)
  
  Z0 <- fit_summary$summary[1:p,"mean"]
  
  Z0_2.5 <- fit_summary$summary[1:p,"2.5%"]
  Z0_97.5 <- fit_summary$summary[1:p,"97.5%"]
  
  names(Z0) <- colnames(table)
  names(Z0_2.5) <- colnames(table)
  names(Z0_97.5) <- colnames(table)
  
  
  R <- fit_summary$summary[(2*n+p+1):(2*n+p+p*p),1]
  R_mat <- matrix(R, nrow=p, byrow = TRUE)
  
  R_mat_2.5 <- matrix(fit_summary$summary[(2*n+p+1):(2*n+p+p*p),"2.5%"], nrow=p, byrow = TRUE)
  R_mat_97.5 <- matrix(fit_summary$summary[(2*n+p+1):(2*n+p+p*p),"97.5%"], nrow=p, byrow = TRUE)
  
  rownames(R_mat) <- colnames(R_mat) <- colnames(table)
  rownames(R_mat_2.5) <- colnames(R_mat_2.5) <- colnames(table)
  rownames(R_mat_97.5) <- colnames(R_mat_97.5) <- colnames(table)
  
  lambda <- round(fit_summary$summary[nrow(fit_summary$summary)-1,1],2)
  lambda_2.5 <- round(fit_summary$summary[nrow(fit_summary$summary)-1,4],2)
  lambda_97.5 <- round(fit_summary$summary[nrow(fit_summary$summary)-1,8],2)
  print(paste0("Pagel's lambda: ", lambda, ", 95% CI: [",lambda_2.5, "; ", lambda_97.5, "]"))
  
  dir.create(file.path(code_path, "plot_ABDOMEN/"), showWarnings = FALSE)
  
  pdf(paste0(code_path, "/plot_ABDOMEN/results_Z0_",name, ".pdf"), width=5, height=7)
  names(Z0) <- gsub("_", " - ", names(Z0))
  list_colors_plot <- list_colors
  names(list_colors_plot) <- gsub("_", " - ", names(list_colors_plot))
  df=data.frame(cbind(names(Z0),Z0))
  df$V1 <- as.factor(df$V1)
  df$V1 <-  factor(df$V1,levels = df$V1)
  df$Z0 <- as.numeric(df$Z0)
  print(ggplot(df, aes(x="", y=Z0, fill=V1)) +
          geom_bar(stat="identity", width=1) +  guides(fill=guide_legend(title=" "))+
          scale_fill_manual(values=list_colors_plot)+
          coord_polar("y", start=0, direction=-1) + theme_void())
  dev.off()
  
  # Plot Covariances
  
  pdf(paste0(code_path, "/plot_ABDOMEN/results_covariances_", name, ".pdf"), width=8, height=6)
  
  R_cov <- R_mat
  diag(R_cov) <- 0
  rownames(R_cov) <- colnames(R_cov) <- gsub("_", " - ", colnames(R_cov))
  
  ncol <- 1000
  my.cols <- suppressWarnings(colorRampPalette(brewer.pal(11, "RdBu")))(ncol)
  my.cols[(1000/2-1):(1000/2+1)] <- "#ffffff"
  print(levelplot(R_cov, at=seq(-max(abs(R_cov)), max(abs(R_cov)), len=ncol), col.regions = my.cols, xlab="", ylab="", scales=list(x=list(rot=45), alternating=2, tck = c(0,1))))
  
  # Plot only significant covariances
  R_cov[intersect(which(R_mat_2.5<0), which(R_mat_97.5>0))] <- 0
  
  if (!all(R_cov==0)){
    print(levelplot(R_cov, at=seq(-max(abs(R_cov)), max(abs(R_cov)), len=ncol), col.regions = my.cols, xlab="", ylab="", scales=list(x=list(rot=45), alternating=2, tck = c(0,1))))
  }
  dev.off()
  
  
  ## Plot variances Mat R 
  
  pdf(paste0(code_path, "/plot_ABDOMEN/results_variances_", name, ".pdf"), width=8, height=6)
  
  R_var <- R_mat*0
  diag(R_var) <- diag(R_mat)
  rownames(R_var) <- colnames(R_var) <- gsub("_", " - ", colnames(R_var))
  
  ncol <- 1000
  my.cols <- suppressWarnings(colorRampPalette(brewer.pal(11, "YlGn")))(ncol)
  my.cols[1] <- "transparent"
  
  if (!all(R_var==0)){
    colnames(R_var) <- paste0(colnames(R_var), " (",round(colSums(table)/sum(table)*100,1), "%)")
    print(levelplot(R_var, at=seq(0, max(c(max(abs(R_var)))), len=ncol), col.regions = my.cols, xlab="", ylab="", scales=list(x=list(rot=45), alternating=2, tck = c(0,1))))
    
  }
  dev.off()
  
  
  
  # Ancestral states at all nodes 
  
  # Build a matrix with tip and internal covariances
  vcvPhyloInternal <- function(tree){
    nbtip <- Ntip(tree)
    dis <- dist.nodes(tree)
    MRCA <- mrca(tree, full = TRUE)
    M <- dis[as.character(nbtip + 1), MRCA]
    dim(M) <- rep(sqrt(length(M)), 2)
    return(M)
  }
  
  
  logX=log(table*exp(fit_summary$summary[(p+1):(p+n)]))
  
  logX <- logX[tree$tip.label,] # should not change anything
  
  
  # transform the tree
  tree_lambda <- tree
  tree_lambda$edge.length <- tree_lambda$edge.length*lambda
  tree_lambda$edge.length[which(tree_lambda$edge[,2] %in% 1:n)] <- tree_lambda$edge.length[which(tree_lambda$edge[,2] %in% 1:n)] + (1-lambda)
  tree_lambda$edge.length <- tree_lambda$edge.length/max(node.depth.edgelength(tree_lambda)) # useless
  
  # covariance for the nodes
  V <- vcvPhyloInternal(tree_lambda)
  indice <- (1:n)
  AY <- V[-indice,indice]
  vY <- V[indice,indice]
  
  # Ancestral state at the root
  one <- t(rbind(rep(1,n)))
  logZ0 <- log(Z0)
  
  
  # states at the nodes
  if (Nnode(tree)==n-1){ # the tree must be rooted and binary
    
    state_nodes <- (AY%*%pseudoinverse(vY)%*%(logX-one%*%logZ0))+(one[1:(n-1),,drop=F]%*%logZ0)
    colnames(state_nodes) = colnames(table)
    rownames(state_nodes) = paste("node_",n+1:Nnode(tree), sep="")
    
    # plot ancestral states 
    
    p12 <- ggtree::ggtree(tree, ladderize = FALSE) + ggtree::theme_tree(bgcolor = "transparent") + ggtree::geom_treescale(x=1.15, y=1, width=0, color="transparent")
    
    Z0_nodes <- data.frame(exp(state_nodes))
    Z0_nodes$node <- n+1:Nnode(tree)
    pies <- ggtree::nodepie(Z0_nodes, cols=1:(ncol(Z0_nodes)-1), color=list_colors, alpha=1)
    p12_nodes <- ggtree::inset(p12, pies, width=0.1, height=0.1 )
    
    
    pdf(paste0(code_path, "/plot_ABDOMEN/results_ancestral_states_", name, ".pdf"), width=6, height=6)
    plot(p12_nodes+ggtree::geom_tiplab(size=0) + theme(plot.margin = margin(0.1,0.1,0.1,0.1, "cm")))
    plot(p12_nodes+ggtree::geom_tiplab(size=2) + theme(plot.margin = margin(0.1,0.1,0.1,0.1, "cm")))
    dev.off()
  }
  
}




ABDOMEN_extract_Z0 <- function(tree, table, fit_summary, detection_threshold=1e-05){
  
  # scale the tree
  tree$edge.length <- tree$edge.length/max(node.depth.edgelength(tree)) 
  
  # scale the abundance per row (each row sum must be equal to 1)
  for (i in 1:nrow(table)) {table[i,] <- table[i,]/sum(table[i,])}
  
  # reorder table as tree$tip.label
  table <- table[tree$tip.label,]
  
  while (length(table[which(table<detection_threshold)])>0){
    table[table<detection_threshold] <- detection_threshold
    for(i in 1:nrow(table)) {table[i,] <- table[i,]/sum(table[i,])}
  }
  
  p <- ncol(table)
  n <- nrow(table)
  
  Z0 <- fit_summary$summary[1:p,"mean"]
  
  Z0_2.5 <- fit_summary$summary[1:p,"2.5%"]
  Z0_97.5 <- fit_summary$summary[1:p,"97.5%"]
  
  names(Z0) <- colnames(table)
  names(Z0_2.5) <- colnames(table)
  names(Z0_97.5) <- colnames(table)
  
  all_Z0 <- data.frame(Z0=(Z0), lower_bound=(Z0_2.5), upper_bound=(Z0_97.5))
  
  return(all_Z0)
  
}


ABDOMEN_extract_Z0_nodes <- function(tree, table, fit_summary, detection_threshold=1e-05){
  
  # scale the tree
  tree$edge.length <- tree$edge.length/max(node.depth.edgelength(tree)) 
  
  # scale the abundance per row (each row sum must be equal to 1)
  for (i in 1:nrow(table)) {table[i,] <- table[i,]/sum(table[i,])}
  
  # reorder table as tree$tip.label
  table <- table[tree$tip.label,]
  
  while (length(table[which(table<detection_threshold)])>0){
    table[table<detection_threshold] <- detection_threshold
    for(i in 1:nrow(table)) {table[i,] <- table[i,]/sum(table[i,])}
  }
  
  p <- ncol(table)
  n <- nrow(table)
  
  Z0 <- fit_summary$summary[1:p,"mean"]
  names(Z0) <- colnames(table)
  
  lambda <- round(fit_summary$summary[nrow(fit_summary$summary)-1,1],2)

  # Ancestral states at all nodes 
  
  # Build a matrix with tip and internal covariances
  vcvPhyloInternal <- function(tree){
    nbtip <- Ntip(tree)
    dis <- dist.nodes(tree)
    MRCA <- mrca(tree, full = TRUE)
    M <- dis[as.character(nbtip + 1), MRCA]
    dim(M) <- rep(sqrt(length(M)), 2)
    return(M)
  }
  
  logX=log(table*exp(fit_summary$summary[(p+1):(p+n)]))
  
  logX <- logX[tree$tip.label,] # should not change anything
  
  # transform the tree
  tree_lambda <- tree
  tree_lambda$edge.length <- tree_lambda$edge.length*lambda
  tree_lambda$edge.length[which(tree_lambda$edge[,2] %in% 1:n)] <- tree_lambda$edge.length[which(tree_lambda$edge[,2] %in% 1:n)] + (1-lambda)
  tree_lambda$edge.length <- tree_lambda$edge.length/max(node.depth.edgelength(tree_lambda)) # useless
  
  # covariance for the nodes
  V <- vcvPhyloInternal(tree_lambda)
  indice <- (1:n)
  AY <- V[-indice,indice]
  vY <- V[indice,indice]
  
  # Ancestral state at the root
  one <- t(rbind(rep(1,n)))
  logZ0 <- log(Z0)
  
  # states at the nodes
  if (Nnode(tree)==n-1){ # the tree must be rooted and binary
    
    state_nodes <- (AY%*%pseudoinverse(vY)%*%(logX-one%*%logZ0))+(one[1:(n-1),,drop=F]%*%logZ0)
    colnames(state_nodes) = colnames(table)
    rownames(state_nodes) = paste("node_",n+1:Nnode(tree), sep="")
    
    Z0_nodes <- data.frame(exp(state_nodes))
    Z0_nodes$node <- n+1:Nnode(tree)
    
    Z0_nodes$MRCA <- NA
    
    for (i in 1:nrow(Z0_nodes)){
      Z0_nodes[i, 1:p] <- Z0_nodes[i, 1:p]/sum(Z0_nodes[i, 1:p])
      Z0_nodes$MRCA[i] <- paste(sort(extract.clade(tree, node=Z0_nodes$node[i])$tip.label),collapse = "-")
    }
  
    return(Z0_nodes)

  } else {print("Warning: the tree must be rooted and binary.")}
}



ABDOMEN_extract_lambda <- function(tree, table, fit_summary, detection_threshold=1e-05){
  
  # scale the tree
  tree$edge.length <- tree$edge.length/max(node.depth.edgelength(tree)) 
  
  # scale the abundance per row (each row sum must be equal to 1)
  for (i in 1:nrow(table)) {table[i,] <- table[i,]/sum(table[i,])}
  
  # reorder table as tree$tip.label
  table <- table[tree$tip.label,]
  
  while (length(table[which(table<detection_threshold)])>0){
    table[table<detection_threshold] <- detection_threshold
    for(i in 1:nrow(table)) {table[i,] <- table[i,]/sum(table[i,])}
  }
  
  p <- ncol(table)
  n <- nrow(table)
  
  lambda <- round(fit_summary$summary[nrow(fit_summary$summary)-1,1],2)
  lambda_2.5 <- round(fit_summary$summary[nrow(fit_summary$summary)-1,4],2)
  lambda_97.5 <- round(fit_summary$summary[nrow(fit_summary$summary)-1,8],2)
  
  all_lambda <- c(lambda,lambda_2.5, lambda_97.5)
  names(all_lambda) <- c("Pagels_lambda", "lower_bound", "upper_bound")
  
  return(all_lambda)
  
}

ABDOMEN_extract_R <- function(tree, table, fit_summary, detection_threshold=1e-05){
  
  # scale the tree
  tree$edge.length <- tree$edge.length/max(node.depth.edgelength(tree)) 
  
  # scale the abundance per row (each row sum must be equal to 1)
  for (i in 1:nrow(table)) {table[i,] <- table[i,]/sum(table[i,])}
  
  # reorder table as tree$tip.label
  table <- table[tree$tip.label,]
  
  while (length(table[which(table<detection_threshold)])>0){
    table[table<detection_threshold] <- detection_threshold
    for(i in 1:nrow(table)) {table[i,] <- table[i,]/sum(table[i,])}
  }
  
  p <- ncol(table)
  n <- nrow(table)
  
  R <- fit_summary$summary[(2*n+p+1):(2*n+p+p*p),1]
  R_mat <- matrix(R, nrow=p, byrow = T)
  
  R_mat_2.5 <- matrix(fit_summary$summary[(2*n+p+1):(2*n+p+p*p),"2.5%"], nrow=p, byrow = T)
  R_mat_97.5 <- matrix(fit_summary$summary[(2*n+p+1):(2*n+p+p*p),"97.5%"], nrow=p, byrow = T)
  
  rownames(R_mat) <- colnames(R_mat) <- colnames(table)
  rownames(R_mat_2.5) <- colnames(R_mat_2.5) <- colnames(table)
  rownames(R_mat_97.5) <- colnames(R_mat_97.5) <- colnames(table)
  
  R_signif <- R_mat
  R_signif[intersect(which(R_mat_2.5<0), which(R_mat_97.5>0))] <- 0
  
  
  all_R <- list(R=R_mat, R_lower_bound=R_mat_2.5, R_upper_bound=R_mat_97.5, R_signif=R_signif)
  
  return(all_R)
  
}

