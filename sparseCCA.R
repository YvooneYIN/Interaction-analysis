#!/usr/bin/env Rscript
#########
#setwd("/Users/yinwen/myfiles/Projects/2NAFLD/0.8sparseCCA/")
rm(list= ls())
options(warn = -1)
library(optparse)
option_list <- list(
  make_option(c("-s","--script_path"),type = 'character',default = NULL,
              help = 'Path to the script', metavar = "/SCRIPTS/PATH"),
  make_option(c("--input_predictor_file"),type = 'character',default = NULL,
              help = 'Name of  the input predictor(x) file (csv file,column:taxa name, row:sample)', metavar = "/INPUT/X/FILE/NAME"),
  make_option(c("--input_independent_file"),type = 'character',default = NULL,
              help = 'Name of the input independent(y) file (csv file,column:gene name, row:sample)', metavar = "/INPUT/Y/FILE/NAME"),
  make_option(c("-n","--thread"),type = 'integer',default = 12,
              help = 'Number of threads', metavar = "12"),
  make_option(c("-o","--output_path"),type = 'character',default = NULL,
              help = 'Path to the output', metavar = "/OUTPUT/PATH"),
  make_option(c("--name"),type = 'character',default = NULL,
              help = 'Name of output', metavar = "Control_S")
)

opt_parser <- OptionParser(usage = "./sparseCCA.R --script_path /SCRIPTS/PATH --input_predictor_file /INPUT/X/FILE/NAME --input_independent_file /INPUT/Y/FILE/NAME --thread 12 --output_path /OUTPUT/PATH --name Control_S",option_list = option_list)
opt <- parse_args(opt_parser)


script_path <- opt$script_path
#script_path <- '../0.8sparseCCA/'
input_independent_file <- opt$input_independent_file
#input_independent_file <- '../0.7lasso/01Input/gene_control.csv'
input_predictor_file <- opt$input_predictor_file
#input_predictor_file <- '../0.7lasso/01Input/taxa_S_control.csv'
output_path <- opt$output_path
#output_path <- '../0.8sparseCCA/Output'
name <- opt$name
#name <- 'Control_S'
thread <- opt$thread
#thread <- 2
suppressPackageStartupMessages({
  library(data.table)
  library(PMA)
  library(parallel)
})

################# Functions ################

run_sparseCCA <- function(X, Z, CCA.K, penaltyX, penaltyZ, vInit=NULL, outputFile=NULL){
  CCA.out <-  CCA(X,Z,typex="standard",typez="standard",K=CCA.K,
                  penaltyx=penaltyX,penaltyz=penaltyZ,
                  v=vInit)
  if(!is.null(outputFile)){
    sink(outputFile)
    print(CCA.out)
    sink()
  }
  
  ## add rownames to output factors
  rownames(CCA.out$u) <- colnames(X)
  rownames(CCA.out$v) <- colnames(Z)
  ## compute contribution of selected features to each of the samples.
  CCA_var_genes <- X %*% CCA.out$u ## canonical variance for genes 
  CCA_var_microbes <- Z %*% CCA.out$v ## canonical variance for microbes
  
  return(list(CCA.out, CCA_var_genes, CCA_var_microbes))
  
}

get_avg_features <- function(cca_cov, CCA.K){
  num_features <- 0
  for(k in 1:CCA.K){
    num_features <- num_features + length(which(cca_cov[,k]!=0))
  }
  avg_features <- num_features/CCA.K
}

save_CCA_components <- function(CCA.out, CCA.K, dirname){
  ## Print canonical covariates in files 
  for(i in CCA.K){
    print(paste0("Writing significant component = ", i))
    selected_X <- which(CCA.out$u[,i]!=0) 
    selected_X <- rownames(CCA.out$u)[selected_X]
    coeff_X <- unname(CCA.out$u[selected_X,i])
    selected_Z <- which(CCA.out$v[,i]!=0)
    selected_Z <- rownames(CCA.out$v)[selected_Z]
    coeff_Z <- unname(CCA.out$v[selected_Z,i])
    ## Make all vectors of same length to avoid repetition of elements from shorter vectors.
    n <- max(length(selected_X), length(selected_Z))
    length(selected_X) <- n                      
    length(selected_Z) <- n
    length(coeff_X) <- n
    length(coeff_Z) <- n
    selected_XZ <- as.data.frame(cbind(gene = selected_X, gene_coeff = coeff_X,
                                       taxa = selected_Z, taxa_coeff = coeff_Z))
    write.table(selected_XZ, file=paste0(dirname,"gene_taxa_component_",i,".txt"), sep = "\t", col.names = NA,quote = F)
  }
  
}

tune_params_grid_search <- function( X, Y){
  scoreXcv <- c()
  scoreYcv <- c()
  penaltyX <- seq(0.05,1,length=100)
  penaltyY <- seq(0.05,1,length=100)
  corr_demo <- matrix(nrow = length(penaltyX), ncol =  length(penaltyY))
  num_samples <- nrow(Y)
  start_time <- Sys.time()
  for( i in 1:length(penaltyX)){
    for(j in 1:length(penaltyY)){
      # print(paste0("Index: i = ",i,", j =", j)); flush.console()
      for(k in 1:num_samples){
        #i = 1
        #j = 1
        #k = 1
        print(paste0("Index: i = ",i,", j =", j," k = ",k)); flush.console()
        #compute weights with sample k held out:
        #Default niter = 15 edited to 5 to speed this up.
        res <- CCA(X[-k,],Y[-k,], penaltyx = penaltyX[i], penaltyz = penaltyY[j], K=1, niter = 5, trace = F, standardize = T)
        ## Compute scores for k'th sample for first pair of canonical variables
        ## Take weight of features (res$u and res$v) computed using all except 
        ## the kth sample and multiple it by values for the kth sample in the 
        ## feature matrix X and Y.
        scoreXcv[k] <- X[k,]%*%res$u ## single value
        scoreYcv[k] <- Y[k,]%*%res$v ## single value
      }
      ## correlation between scores for X and Y for all held out samples.
      corr_demo[i,j] = cor(scoreXcv,scoreYcv) 
    }
  }
  end_time <- Sys.time()
  time_elapsed <- end_time - start_time
  print(paste0("Time elapsed for param tuning = ", time_elapsed)); flush.console()
  
  row.names(corr_demo) <- as.character(penaltyX)
  colnames(corr_demo) <- as.character(penaltyY)
  
  corr_demo_df <- as.data.frame(corr_demo)
  rownames(corr_demo_df)
  colnames(corr_demo_df)
  
  ##identify best penalty parameters
  # find index with max absolute corr
  bestpenalty <- which(abs(corr_demo) == max(abs(corr_demo)), arr.ind = TRUE)
  bestpenalty
  bestpenaltyX <- penaltyX[bestpenalty[1]]
  bestpenaltyY <- penaltyY[bestpenalty[2]]
  
  return (c(bestpenaltyX,bestpenaltyY))
}

test_significance_LOOCV <- function(X, Y, bestpenaltyX, bestpenaltyY, num_components){
  cca.k = num_components
  scoresXcv <- matrix(nrow = nrow(X), ncol = cca.k)
  scoresYcv <-  matrix(nrow = nrow(Y), ncol = cca.k)
  corr_pval <- c()
  corr_r <- c()
  for(i in 1:nrow(Y)){ #n = no. of samples
    #compute weights with sample i held out:
    res <- CCA(X[-i,],Y[-i,], penaltyx=bestpenaltyX, penaltyz=bestpenaltyY, K=cca.k, trace = F) ## default niter = 15 which is spit out when trace = T (default)
    ###compute scores for i'th sample for each component (pair of canonical variables)
    for(j in 1:cca.k){
      #print(paste0("i = ", i," K = ", j)); flush.console()
      scoresXcv[i,j] <- X[i,]%*%res$u[,j]
      scoresYcv[i,j] <- Y[i,]%*%res$v[,j]
    }
  }
  ## Test for each components
  for(j in 1:cca.k){
    corr <- cor.test(scoresXcv[,j],scoresYcv[,j]) ## Pearson correlation.
    corr_pval[j] <- corr$p.value
  }
  corr_pval
}

tune_params_grid_search_parallel <- function(X_matrix,Y_matrix,thread,output_path,penaltyX,penaltyY,num_samples,name) {
  corr_demo <- matrix(nrow = length(penaltyX), ncol = length(penaltyY))
  no_cores <- thread
  cl <- makeCluster(no_cores)
  
  # 导出所需变量到每个核心
  clusterExport(cl, varlist = c("X_matrix", "Y_matrix", "penaltyX", "penaltyY", "num_samples", "CCA",'output_path','sprintf','name'))
  
  start_time <- Sys.time()
  
  # 使用并行化处理参数搜索
  results <- parLapply(cl, seq_len(length(penaltyX) * length(penaltyY)), function(idx) {
    i <- ((idx - 1) %% length(penaltyX)) + 1
    j <- ((idx - 1) %/% length(penaltyX)) + 1
    scoreXcv <- numeric(num_samples)
    scoreYcv <- numeric(num_samples)
    
    for (k in 1:num_samples) {
      cat(sprintf("Index: i = %d, j = %d, k = %d\n", i, j, k), file = paste0(output_path,'/' ,name,'/log.txt'), append = TRUE)
      res <- CCA(X_matrix[-k,], Y_matrix[-k,], penaltyx = penaltyX[i], penaltyz = penaltyY[j], K = 1, niter = 5, trace = F, standardize = T)
      scoreXcv[k] <- X_matrix[k,] %*% res$u
      scoreYcv[k] <- Y_matrix[k,] %*% res$v
    }
    
    cor_val <- cor(scoreXcv, scoreYcv)
    return(list(index = c(i, j), cor_val = cor_val))
  })
  
  stopCluster(cl)
  
  # 处理并行计算的结果
  for (res in results) {
    corr_demo[res$index[1], res$index[2]] <- res$cor_val
  }
  
  end_time <- Sys.time()
  print(paste0("Time elapsed for param tuning = ", end_time - start_time))
  
  # 确定最佳惩罚参数
  bestpenalty <- which(abs(corr_demo) == max(abs(corr_demo)), arr.ind = TRUE)
  bestpenaltyX <- penaltyX[bestpenalty[1]]
  bestpenaltyY <- penaltyY[bestpenalty[2]]
  
  return(c(bestpenaltyX, bestpenaltyY))
}

########## Tune and run sparse CCA #########

# #### load data

independent_df <- as.matrix(read.csv(paste0(script_path,'/',input_independent_file),header = T,row.names = 1,quote = '',check.names = F))
predictor_df <- as.matrix(read.csv(paste0(script_path,'/',input_predictor_file),header = T,row.names = 1,check.names = F,quote = ''))

# ## Ensure same sampleIDs in both genes and microbes data before sparse CCA
tryCatch({
  stopifnot(all(rownames(independent_df) == rownames(predictor_df)))
},
error =  function(e){
  print('Ensure same sampleIDs in both independent_df and predictor_df before sparse CCA!')
}
)
if (!dir.exists(paste0(output_path,'/',name))){
  dir.create(paste0(output_path,'/',name),recursive = T)
}

# 
# ## set penalty parameters
# bestpenaltyX <- 0.05
# bestpenaltyY <- 0.3222
# 
# ## SKIP if using pre-computed values above
# ## select tuning parameters
X_matrix = independent_df
Y_matrix = predictor_df
penaltyX <- seq(0.05, 0.4, length = 5)
penaltyY <- seq(0.05, 0.4, length = 5)

num_samples <- nrow(Y_matrix)

bestPenalty <- tune_params_grid_search_parallel(X_matrix,Y_matrix,thread,output_path,penaltyX,penaltyY,num_samples,name)
bestpenaltyX <- bestPenalty[1]
bestpenaltyY <- bestPenalty[2]
# 
# #### Run sparse CCA
# 
# ## Set the number of desired components
cca.k = 10


# 
# ## Run sparse CCA using selected tuning param using permutation search
cca <- run_sparseCCA( predictor_df,independent_df, cca.k, bestpenaltyX, bestpenaltyY,
                     outputFile=paste0(output_path,"/",name,'/',bestpenaltyX,"_",bestpenaltyY,".txt"))
# 
# ## average number of genes and microbes in resulting components
avg_independent <- get_avg_features(cca[[1]]$u, cca.k)
avg_independent
# 
avg_predictor <- get_avg_features(cca[[1]]$v, cca.k)
avg_predictor
# 
# #### Test significance of components using LOOCV
CCA_pval <- test_significance_LOOCV(as.matrix(independent_df), as.matrix(predictor_df), bestpenaltyX, bestpenaltyY, cca.k)
# 
length(which(CCA_pval < 0.1))
which(CCA_pval < 0.1)
# 
CCA_padj <- p.adjust(CCA_pval, method = "BH")
CCA_padj
# 
length(which(CCA_padj < 0.1))
which(CCA_padj < 0.1)
# 
# #### Output significant components
sig_cutoff <- 0.05
sig <- which(CCA_padj < sig_cutoff)
dirname <- paste0(output_path,'/',name,'/','components','/')
# ## This will return FALSE if the directory already exists or is uncreatable,
# ## and TRUE if it didn't exist but was succesfully created.
ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)
save_CCA_components(cca[[1]],sig,dirname)




