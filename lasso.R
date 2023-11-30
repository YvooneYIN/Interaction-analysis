#!/usr/bin/env Rscript
rm(list = ls())
options(warn = -1)
library(optparse)
############
option_list <- list(
  make_option(c("-s","--script_path"),type = 'character',default = NULL,
              help = 'Path to the script', metavar = "/SCRIPTS/PATH"),
  make_option(c("--input_predictor_file"),type = 'character',default = NULL,
              help = 'Name of  the input predictor(x) file (csv file,column:taxa name, row:sample)', metavar = "/INPUT/X/FILE/NAME"),
  make_option(c("--input_independent_file"),type = 'character',default = NULL,
              help = 'Name of the input independent(y) file (csv file,column:gene name, row:sample)', metavar = "/INPUT/Y/FILE/NAME"),
  make_option(c("-n","--thread"),type = 'integer',default = 6,
              help = 'Number of threads', metavar = "6"),
  make_option(c("-o","--output_path"),type = 'character',default = NULL,
              help = 'Path to the output', metavar = "/OUTPUT/PATH"),
  make_option(c("-x","--predictor"),type = 'character',default = NULL,
              help = 'Name of predictor variable', metavar = "taxa"),
  make_option(c("-y","--independent"),type = 'character',default = NULL,
              help = 'Name of independent variable', metavar = "gene")
)
opt_parser <- OptionParser(usage = "lasso.R -s . --input_predictor_file /INPUT/X/FILE/NAME --input_independent_file /INPUT/Y/FILE/NAME -n 6 -o  /OUTPUT/PATH -x taxa -y gene",
                           option_list = option_list)
opt <- parse_args(opt_parser)

script_path <- opt$script_path
#script_path <- '../0.7lasso/'
input_independent_file <- opt$input_independent_file
#input_independent_file <- '../0.7lasso/01Input/gene_control.csv'
input_predictor_file <- opt$input_predictor_file
#input_predictor_file <- '../0.7lasso/01Input/taxa_S_control.csv'
thread <- opt$thread
#thread <- 3
output_path <- opt$output_path
#output_path <- '../0.7lasso/02Output'
predictor <- opt$predictor
#predictor <- 'taxa'
independent <- opt$independent
#independent <- 'gene'
check_file_existence <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
}
check_file_existence(input_independent_file)
check_file_existence(input_predictor_file)


############
suppressPackageStartupMessages({
  library(parallel)
  library(doParallel)
  library(data.table)
  library(glmnet)
  library(hdi)
  library(stabs)
  library(utils)
})
############
cl <- makeCluster(thread)
clusterEvalQ(cl,library(glmnet))
clusterEvalQ(cl,library(data.table))
clusterEvalQ(cl,library(hdi))
clusterEvalQ(cl,library(stabs))
safe_file_read <- function(file_path) {
  tryCatch(
    {
      read.delim(file_path, sep=',', head=T, row.names = 1, check.names = F, stringsAsFactors = F)
    },
    error = function(e) {
      cat("Error in reading file:", e$message, "\n")
      NULL  # 返回 NULL 或其他适当的值
    }
  )
}
fit.cv.lasso <- function(x, y_i, kfold){
 # fit.model <- fit.cv.lasso(x, y_i,  kfold = length(y_i))
  lambdas = NULL
  r.sqr.final <- numeric()
  r.sqr.final.adj <- numeric()
  r.sqr.CV.test <- numeric()
  lasso.cv.list <- list()
  ## glmnet CV
  cv.fit <- cv.glmnet(x, y_i, alpha=1, nfolds=kfold, type.measure = "mse", keep =TRUE, grouped=FALSE, standardize = T)  
  lambdas = data.frame(cv.fit$lambda,cv.fit$cvm)
  ## get best lambda -- lambda that gives min cvm
  bestlambda <- cv.fit$lambda.min
  bestlambda_index <- which(cv.fit$lambda == bestlambda)
  ## Get R^2 of final model
  final_model <- cv.fit$glmnet.fit
  r_sqr_final_model <- cv.fit$glmnet.fit$dev.ratio[bestlambda_index]
  ## Get adjusted R^2
  r_sqr_final_adj <- adj_r_squared(r_sqr_final_model, n = nrow(x), 
                                   p = sum(as.vector(coef(cv.fit$glmnet.fit, 
                                                          s = cv.fit$lambda.min)) > 0))
  return(list(bestlambda = bestlambda, r.sqr = r_sqr_final_model, 
              r.sqr.adj = r_sqr_final_adj
  ))
}
estimate.sigma.loocv <- function(x, y_i, bestlambda, tol) {
  ## Fit a lasso object
  lasso.fit = glmnet(x,y_i,alpha = 1) ## this is same as cv.fit$glmnet.fit from loocv code below.
  beta <- as.vector(coef(lasso.fit, s = bestlambda)) ## This gives coefficients of fitted model, not predicted coeff.
  # try(if(length(which(abs(beta) > tol)) > n) stop(" selected predictors more than number of samples! Abort function"))
  y = as.vector(y_i)
  yhat = as.vector(predict(lasso.fit, newx = x, s = bestlambda))
  ## predicted coefficients, same as coefficient of fitted model lasso. Either one is fine.
  # beta = predict(lasso.fit,s=bestlambda, type="coef")
  df = sum(abs(beta) > tol) ## Number of non-zero coeff. Floating-point precision/tolerance used instead of checking !=0
  n = length(y_i)
  ss_res = sum((y - yhat)^2)
  if((n-df-1) >= 1) {
    sigma = sqrt(ss_res / (n-df-1))
    sigma.flag = 0
  } else{
    sigma = 1 ## conservative option
    # sigma = ss_res ## lenient option
    sigma.flag = 2
  }
  
  
  return(list(sigmahat = sigma, sigmaflag = sigma.flag, betahat = beta)) ## we return beta to be used later in hdi function.
  
}
r_squared <- function(y, yhat) {
  ybar <- mean(y)
  ## Total SS
  ss_tot <- sum((y - ybar)^2)
  ## Residual SS
  ss_res <- sum((y - yhat)^2)
  ## R^2 = 1 - ss_res/ ss_tot
  1 - (ss_res / ss_tot)
}
adj_r_squared <- function(r_squared, n, p) {
  1 - (1 - r_squared) * (n - 1) / (n - p - 1)
}
clusterExport(cl, "fit.cv.lasso")
clusterExport(cl, "estimate.sigma.loocv")
clusterExport(cl, "r_squared")
clusterExport(cl, "adj_r_squared")
###################
independent_file <- as.matrix(safe_file_read(input_independent_file))
#independent_file <- t(independent_file)
if (is.null(independent_file)) {
  stop("Failed to read input gene file.")
}
predictor_file <- as.matrix(safe_file_read(input_predictor_file))
if (is.null(predictor_file)) {
  stop("Failed to read input micro file.")
}
clusterExport(cl, "independent_file")
clusterExport(cl, "predictor_file")
tryCatch({
  stopifnot(rownames(independent_file) == rownames(predictor_file))
  },
  error = function(e){
    print('Please check the form of input files, make sure the rows correspond to samples!')
  }
)

##
registerDoParallel(cl)
start <- Sys.time()
print(paste0("Number of threads: ", getDoParWorkers()))
print('Start analysing ...')

lasso_function <- function(i,independent_file,predictor_file,seed,output_path,script_path,ncol,independent,predictor){
 # i = 1
  #seed = 1129
  if (!is.null(seed)) {
    set.seed(seed + i)
  }
  ##for each independent variable, find the best x to predict y
  x <- predictor_file
  y <- independent_file
  y_i <- y[,i]
  y_name <- colnames(y)[i]
  #length(y_i) is the same as the sample number
  fit.model <- fit.cv.lasso(x, y_i,  kfold = length(y_i))
  bestlambda <- fit.model$bestlambda
  r.sqr <- fit.model$r.sqr 
  r.sqr.adj <- fit.model$r.sqr.adj
  sigma.myfun <- estimate.sigma.loocv(x, y_i, bestlambda, tol=1e-4)
  sigma <- sigma.myfun$sigmahat
  beta <- as.vector(sigma.myfun$betahat)[-1]
  sigma.flag <- sigma.myfun$sigmaflag
  lasso.proj.fit <- lasso.proj(x, y_i, multiplecorr.method = "BH", betainit = beta, sigma = sigma, suppress.grouptesting = T)
  lasso.ci <- as.data.frame(confint(lasso.proj.fit, level = 0.95))
  lasso.df <- data.frame(independent = rep(y_name, length(lasso.proj.fit$pval)),
                         predictor = names(lasso.proj.fit$pval.corr),
                         r.sqr = r.sqr,r.sqr.adj = r.sqr.adj,
                         pval = lasso.proj.fit$pval,padj = lasso.proj.fit$pval.corr,
                         ci.lower = lasso.ci$lower, ci.upper = lasso.ci$upper,
                         sigma = sigma, sigma.flag = sigma.flag,
                         row.names=NULL)
  lasso.df <- lasso.df[order(lasso.df$pval),]
  stab.glmnet <- stabsel(x = x, y = y_i,
                         fitfun = glmnet.lasso, cutoff = 0.6,
                         PFER = 1)
  predictor.selected <- names(stab.glmnet$selected)
  lasso.df$selected <- ifelse(lasso.df$predictor %in% predictor.selected,'yes','no')
  colnames(lasso.df)[1:2] <- c(independent,predictor)
  write.table(lasso.df,paste0(output_path,"/Rawfile/",colnames(y)[i],".txt"),quote = F)
  cat(sprintf("完成 %.2f%%\n", round(100* i/ncol,4)), file = paste0(output_path,"progress_log.txt"), append = TRUE)
  return(lasso.df)
}

if (!dir.exists(output_path)) {
  dir.create(output_path,recursive = T)
}
if (!dir.exists(paste0(output_path,'/Rawfile'))) {
  dir.create(paste0(output_path,'/Rawfile'),recursive = T)
}

ncol <- ncol(independent_file)
#nrow <- 3
system.time({
  res <- parLapply(cl,1:ncol,lasso_function,independent_file,predictor_file,seed = 1129,output_path,script_path,ncol,independent,predictor)
})
print('Finish analysed!')
stopCluster(cl)

merged_overlap <- rbindlist(res)
end <- Sys.time()
print(paste0('Spent time: ',end - start))

print('Start correcting for p value ...')
merged_overlap$FDR <- p.adjust(merged_overlap$pval,"BH")
print('Start choosing for stable results ...')
stable_result <- merged_overlap[merged_overlap$selected == 'yes',]
stable_result_sig <- stable_result[stable_result$FDR < 0.05,]
print('Writing output file ...')
write.csv(stable_result_sig,paste0(output_path,'lasso_results_FDR_0.05_stable.csv'),quote = F)
print('Lasso analysis finished!')
