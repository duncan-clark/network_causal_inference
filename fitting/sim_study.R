# This script is the basis for the reproducible example of the sim study for the paper:
# For fuller sim study see sim_study_5.R


# EFFECTS of INTEREST
# - Main effect of outcome
# - k-peer treatment effect (for treated,and untreated)
# - k-peer outcome effect

# DGPs Considered
# - ERNM with transitivity + geodegree + homophily on both
# - ERNM with homophily on both

# Models Considered for each DGP
# - ERNM with transitivity + geodegree + homophily on both
# - ERGM with separated logistic regression
# - MRF with neighbour outcome and treatment counts
# - Logistic regression with neighbour treatment and outcome counts 

library(statnet)
library(ernm)
library(mvtnorm)
library(coda)
library(parallel)
library(doParallel)
library(rstanarm)
source("functions/utils.R")
source("functions/fire.R")

# =============================================================================================
# - Setup + True Effects
# =============================================================================================

# use 100 node network, 50% treatment density is the only thing that is used, since simulations will change the outcomes.
set.seed(1)
net <- network(matrix(rbinom(100**2,1,0.01),nrow =100),directed  = F)
set.vertex.attribute(net,"outcome",as.character(rbinom(100,1,0.5)))
set.vertex.attribute(net,"treatment",as.character(rbinom(100,1,0.5)))
net <- add_treated_neighs(net,"treatment")

#Use uniform prior that drops off at -20 and 20
unif_prior <- function(x){
  if(any(abs(x) > 20) ){
    return(0)
  }
  else{return(1)}
}

# Get the simulated networks:
DGP_formula <- as.formula(net  ~  edges() +    
                            gwesp(0.5,"outcome") +
                            gwdegree(0.5) +
                            homophily("outcome") +
                            nodeCount("outcome") +
                            logistic("outcome","treatment") +
                            logisticNeighbors("outcome","treatment") +
                            logisticNeighborsTopLevel("outcome","outcome")
                          | outcome)
DGP <- ernm(DGP_formula,fitModel = F)
DGP$setThetas(c(-4.5,rep(1,5),0.1,0.1))

DGP_sims <- DGP$sampler$generateSample(10000,1000,1000)
DGP_sims_R <- lapply(DGP_sims[1:100],binaryNet_to_network)

# Infer the "true" main effect
# Note that after any given network is observed the target is the posterior distribution of the true effect not the actual true effect.
# Thus the "true" effect is not of immediate interest - though we hope we are not "far" from it
# Note that k-peer treatment effect for treated and untreated is expected to be quite different due to homophily.

DGP_effects_treated <- extract_effects_from_sims(
  sims = DGP_sims,
  outcome_var = "outcome",
  treatment_var = "treatment",
  treatment_indicator = 0,
  node_treatment_statuses = c(1) # the statuses of the nodes to be considered
)

DGP_effects_not_treated <- extract_effects_from_sims(
  sims = DGP_sims,
  outcome_var = "outcome",
  treatment_var = "treatment",
  treatment_indicator = 0,
  node_treatment_statuses = c(0) # the statuses of the nodes to be considered
)

DGP_effects_both <- extract_effects_from_sims(
  sims = DGP_sims,
  outcome_var = "outcome",
  treatment_var = "treatment",
  treatment_indicator = 0,
  node_treatment_statuses = c(0,1) # the statuses of the nodes to be considered
)

do.call(rbind,list(unlist(DGP_effects_both),unlist(DGP_effects_treated),unlist(DGP_effects_not_treated)))

# =============================================================================================
# - ERNM
# =============================================================================================
# Use theta0 based on ERNM package - just gets the edge density in the ball-park
# Pilot takes ~20 mins on our server
theta0 <-  ernm(DGP_formula,fitModel = F)$thetas()
theta0 <- c(-4.5,rep(1,5),0.1,0.1)
ERNM_pilot <- fire(DGP_formula,
                   prior_dens = unif_prior,
                   theta0 = theta0,
                   ernm_burnin = 10000,
                   ernm_steps = 100,
                   burn_in = 100,
                   samples = 10,
                   n_sims = 500,
                   gamma = 0.2,
                   run_type = "pilot",
                   cores = 9,
                   pilot_chains = 8)

ERNM_pilot$time
ERNM_pilot$median_acceptance
# looks burnt in after around ~ 5000 steps

# make function that can be parallelised easily
ERNM_parallel_func <- function(nw,burn_in,samples,gamma = 0.2,n_sims = 200,timeout = 360){
  form <- DGP_formula
  form[[2]] <- nw
  theta0 <-  ernm(form,fitModel = F)$thetas()
  theta0 <- c(-4.5,rep(1,5),0.1,0.1)
  model <- tryCatch(
    withTimeout(
      fire(form,
           prior_dens = unif_prior,
           theta0 = theta0,
           ernm_burnin = 10000,
           ernm_steps = 100,
           burn_in = burn_in,
           samples = samples,
           n_sims = n_sims,
           gamma = gamma,
           run_type = "main",
           verbose = F,
           cores = NULL,
           pilot_chains = NULL,
           small_results = T),
      timeout = timeout),
    error = function(e){return(NA)}
  )
  return(model)
}

# get the posteriors for each of the networks
# problem is if we fall off one simulation can take a long time ! 
t <- proc.time()
cl <- parallel::makeForkCluster(20)
clusterSetRNGStream(cl, iseed = 1234)
ERNM_models <- parallel::parLapply(cl = cl,
                                   DGP_sims_R,
                                   fun = ERNM_parallel_func,
                                   burn_in = 10000,
                                   samples = 1000,
                                   gamma = 0.05,
                                   n_sims = 500,
                                   timeout = 7200
)
stopCluster(cl)
rm(cl)
t <- proc.time()-t
print(paste("ERNM modelling took ", t[3], " seconds"))

# =============================================================================================
# - MRF
# =============================================================================================
MRF_formula <- as.formula(net  ~ 
                            homophily("outcome") +
                            nodeCount("outcome") +
                            logistic("outcome","treatment") +
                            logisticNeighbors("outcome","treatment") +
                            logisticNeighborsTopLevel("outcome","outcome")
                          | outcome)
theta0 <-  ernm(MRF_formula,fitModel = F)$thetas()

# MRF_pilot <- fire(MRF_formula,
#                    prior_dens = unif_prior,
#                    theta0 = theta0,
#                    ernm_burnin = 10000,
#                    ernm_steps = 100,
#                    burn_in = 10000,
#                    samples = 10,
#                    n_sims = 100,
#                    gamma = 0.3,
#                    run_type = "pilot",
#                    cores = 9,
#                    pilot_chains = 8,
#                    fixed_edges = T,
#                    verbose = T)
# 
# MRF_pilot$median_acceptance
# MRF_pilot$time
# looks burnt in after around ~ 5000 steps

MRF_parallel_func <- function(nw,burn_in,samples,gamma = 0.2,n_sims = 200,timeout = 360){
  form <- MRF_formula
  form[[2]] <- nw
  theta0 <-  c(ernm(form,fitModel = F)$thetas())
  
  model <- tryCatch(
    withTimeout(
      fire(form,
           prior_dens = unif_prior,
           theta0 = theta0,
           ernm_burnin = 10000,
           ernm_steps = 100,
           burn_in = burn_in,
           samples = samples,
           n_sims = n_sims,
           gamma = gamma,
           run_type = "main",
           cores = NULL,
           pilot_chains = NULL,
           fixed_edges = T,
           verbose = F,
           small_results = T),
      timeout = timeout),
    error = function(e){return(NA)}
  )
  return(model)
}

t <- proc.time()
cl <- parallel::makeForkCluster(20)
clusterSetRNGStream(cl, iseed = 1234)
MRF_models <- parallel::parLapply(cl = cl,
                                  DGP_sims_R,
                                  fun = MRF_parallel_func,
                                  burn_in = 10000,
                                  samples = 1000,
                                  gamma = 0.3,
                                  n_sims = 100,
                                  timeout = 3600)
stopCluster(cl)
rm(cl)
t <- proc.time()-t
print(paste("MRF modelling took ", t[3], " seconds"))

# =============================================================================================
# - ERGM + Logistic
# =============================================================================================
ergm_formula <- as.formula(net ~ edges() +    
                             gwdegree(0.5)+
                             gwesp(0.5) +
                             homophily("outcome") +
                             logisticNeighbors("outcome","treatment") +
                             logisticNeighborsTopLevel("outcome","outcome"))

theta0 <-  ernm(ergm_formula,fitModel = F)$thetas()
# ERGM_pilot <- fire(ergm_formula,
#                    prior_dens = unif_prior,
#                    theta0 = theta0,
#                    ernm_burnin = 10000,
#                    ernm_steps = 100,
#                    burn_in = 10000,
#                    samples = 10,
#                    n_sims = 200,
#                    gamma = 0.2,
#                    run_type = "pilot",
#                    cores = 9,
#                    pilot_chains = 8,
#                    fixed_nodes = T)
# 
# ERGM_pilot$time
# ERGM_pilot$median_acceptance
# burnt in after ~ 5000

ERGM_logistic_parallel_func <- function(nw,burn_in,samples,gamma = 0.2,n_sims = 200,timeout = 360){
  form <- ergm_formula
  form[[2]] <- nw
  theta0 <-  ernm(form,fitModel = F)$thetas()
  
  ergm_model <- tryCatch(
    withTimeout(
      fire(form,
           prior_dens = unif_prior,
           theta0 = theta0,
           ernm_burnin = 10000,
           ernm_steps = 100,
           burn_in = burn_in,
           samples = samples,
           n_sims = n_sims,
           gamma = gamma,
           run_type = "main",
           verbose = F,
           cores = NULL,
           pilot_chains = NULL,
           small_results = T,
           fixed_nodes = T),
      timeout = timeout),
    error = function(e){return(NA)}
  )
  
  if(is.na(ergm_model)){return(NA)}
  
  outcome <- binary_convert(as.numeric(nw %v% "outcome"))
  treatment <-  binary_convert(as.numeric(nw %v% "treatment"))
  treatment_neighbors <- as.numeric(nw %v% "treatment_neighbors")
  nw <- add_treated_neighs(nw,"outcome")
  outcome_neighbors <- as.numeric(nw %v% "outcome_neighbors")
  
  model <- rstanarm::stan_glm(outcome ~ treatment + treatment_neighbors + outcome_neighbors,
                              data = data.frame(outcome = outcome,
                                                treatment = treatment,
                                                treatment_neighbors= treatment_neighbors,
                                                outcome_neighbors = outcome_neighbors),
                              family=binomial(link="logit"))
  
  return(list(ergm = ergm_model,glm = model))
}

tmp <- ERGM_logistic_parallel_func(DGP_sims_R[[1]],10,10,0.001,10)

t <- proc.time()
cl <- parallel::makeForkCluster(20)
clusterSetRNGStream(cl, iseed = 1234)
ERGM_logistic_models <- parallel::parLapply(cl = cl,
                                            DGP_sims_R,
                                            fun = ERGM_logistic_parallel_func,
                                            burn_in = 10000,
                                            samples = 1000,
                                            gamma = 0.2,
                                            n_sims = 200,
                                            timeout = 3600)
stopCluster(cl)
rm(cl)
t <- proc.time()-t
print(paste("ERGM Logistic modelling took ", t[3], " seconds"))

# =============================================================================================
# - Logistic
# =============================================================================================
regression_formula = "treatment + treatment_neighbors + outcome_neighbors" 
logistic_parallel_func <- function(nw){
  
  # need to do the glm fit also:
  # swap the treatment so 0 is the "treatment
  # set.vertex.attribute(nw,"treatment", as.character(1-as.numeric(get.vertex.attribute(nw,"treatment"))))
  
  outcome <- binary_convert(as.numeric(nw %v% "outcome"))
  treatment <-  binary_convert(as.numeric(nw %v% "treatment"))
  treatment_neighbors <- as.numeric(nw %v% "treatment_neighbors")
  # make sure we have outcome neighbors since these chagne with the sim
  nw <- add_treated_neighs(nw,"outcome")
  outcome_neighbors <- as.numeric(nw %v% "outcome_neighbors")
  
  glm_data = data = data.frame(outcome = outcome,
                               treatment = treatment,
                               treatment_neighbors= treatment_neighbors,
                               outcome_neighbors = outcome_neighbors)
  
  model <- rstanarm::stan_glm(as.formula(paste("outcome ~",regression_formula)),
                              data = glm_data,
                              family=binomial(link="logit"))
  return(model)
}

tmp <- logistic_parallel_func(DGP_sims_R[[1]])

t <- proc.time()
cl <- parallel::makeForkCluster(20)
clusterSetRNGStream(cl, iseed = 1234)
logistic_models <- parallel::parLapply(cl = cl,
                                       DGP_sims_R,
                                       fun = logistic_parallel_func)
stopCluster(cl)
rm(cl)
t <- proc.time() -t
print(paste("Logistic modelling took ", t[3], " seconds"))


# =============================================================================================
# - Comparison
# =============================================================================================
# function to compare the distributions that can be parallised
# very memory inefficient since passes lots of things to clusters that are unnecessary 
model_names <- c("ERNM_models","MRF_models","ERGM_models","logistic_models")
compare_posteriors <- function(j, # which one of the ERNM models to calculate for
                               ERNM_models,
                               MRF_models,
                               ERGM_models,
                               logistic_models,
                               ERNM_formulas,
                               MRF_formulas,
                               ERGM_formulas,
                               ERGM_nsims,
                               ERNM_nsims,
                               ERNM_draws,
                               ERGM_draws,
                               glm_theta_samples,
                               treatment_var,
                               outcome_var,
                               node_treatment_statuses = c(0,1),
                               treatment_indicator = 0,
                               model_names,
                               small = F,
                               keep = NULL,
                               reldist_method = "gam",
                               return_sims = F){
  
  
  #t<-proc.time()
  # Do the ernm_posteriors
  if(length(ERNM_models) != 0){
    ERNM_posteriors <- tryCatch(extract_posteriors_fast(model = ERNM_models[[j]],
                                                        formula = ERNM_formulas[[j]],
                                                        theta_draws = ERNM_draws,
                                                        fixed_edges = F,
                                                        nsims = ERNM_nsims,
                                                        treatment_var = treatment_var,
                                                        treatment_indicator = treatment_indicator,
                                                        outcome_var = outcome_var,
                                                        node_treatment_statuses = node_treatment_statuses,
                                                        cores = NULL,
                                                        return_sims = return_sims),
                                error = function(e){
                                  print("error with ERNM")
                                  print(j)
                                  return(NA)})
  }
  else{
    ERNM_posteriors = list()
  }
  
  if(length(MRF_models) != 0){
    MRF_posteriors <- tryCatch(extract_posteriors_fast(model = MRF_models[[j]],
                                                       formula = MRF_formulas[[j]],
                                                       theta_draws = ERNM_draws,
                                                       fixed_edges = T, # ONLY change from full ERNM above
                                                       nsims = ERNM_nsims,
                                                       treatment_var = treatment_var,
                                                       treatment_indicator = treatment_indicator,
                                                       outcome_var = outcome_var,
                                                       node_treatment_statuses = node_treatment_statuses,
                                                       cores = NULL,
                                                       return_sims = return_sims),
                               error = function(e){
                                 print("error with MRF")
                                 print(j)
                                 return(NA)})
  }else{
    MRF_posteriors = list()
  }
  
  if(length(ERGM_models) != 0){
    if(length(node_treatment_statuses) !=1){
      tmp <- NULL
    }else{
      tmp <- node_treatment_statuses
    }
    
    ERGM_posteriors <- tryCatch(extract_posteriors_ERGM_glm(model = ERGM_models[[j]],
                                                            formula = ERGM_formulas[[j]],
                                                            nsims = ERGM_nsims,
                                                            ergm_theta_samples = ERGM_draws,
                                                            glm_theta_samples = glm_theta_samples,
                                                            main_index=1,
                                                            k_treat_index=2,
                                                            k_out_index=3,
                                                            treatment_var = treatment_var,
                                                            outcome_var = outcome_var,
                                                            cores = NULL,
                                                            treatment_indicator = treatment_indicator,
                                                            node_treatment_status = tmp),
                                error = function(e){
                                  print("error ERGM")
                                  print(j)
                                  return(NA)})
  }else{
    ERGM_posteriors = list()
  }
  
  
  if(length(ERNM_models) != 0){
    if(length(node_treatment_statuses) !=1){
      tmp <- NULL
    }else{
      tmp <- node_treatment_statuses
    }
    logistic_posteriors <- tryCatch(extract_posteriors_logistic_directly(logistic_models[[j]],
                                                                         theta_samples = glm_theta_samples,
                                                                         main_index = 1,
                                                                         k_treat_index = 2,
                                                                         k_out_index = 3,
                                                                         node_treatment_status = tmp),
                                    error = function(e){
                                      print("error with logistic")
                                      print(j)
                                      return(NA)})
  }else{
    logistic_posteriors = list()
  }
  
  posterior_list <- list(ERNM_posteriors,MRF_posteriors,ERGM_posteriors,logistic_posteriors)
  posterior_list <- posterior_list[unlist(sapply(posterior_list,function(x){length(x) !=0}))]
  means <- data.frame(do.call(rbind,lapply(posterior_list,function(x){sapply(x,mean,na.rm = T)})))
  if(dim(means)[1] == length(model_names)){rownames(means) <- model_names}
  
  KL_dists <- tryCatch(KL_distances_rel_rank(posterior_list,
                                             method = reldist_method),
                       error = function(e){
                         print("error with KL distances")
                         return(NA)
                       })
  if(!return_sims){
    return(list(means = means,
                KL_dists = KL_dists,
                CIs = lapply(posterior_list,function(x){posterior_CIs(x,0.05)})))
  }else{
    return(list(means = means,
                KL_dists = KL_dists,
                CIs = lapply(posterior_list,function(x){posterior_CIs(x$result,0.05)}),
                posteriors = posterior_list))
  }
}

# Function to actually extract the things we want:

# Freq results comparisons that don't really make sense:
freq_results <- function(compare_output,true_effects){
  # show the true effect and the mean mean-a-posteriori for information:
  mean_mean <- lapply(compare_output,function(x){
    return(x$mean)
  })
  # remove any estimates that are NA
  mean_mean <- mean_mean[sapply(mean_mean,function(x){dim(x)[2] != 1})]
  
  # for each estimand of each position count the number of missing values
  num_nas <- lapply(mean_mean,function(x){
    is.na(Matrix::Matrix(as.matrix(x)))})
  num_nas <- Reduce("+",num_nas)
  mean_mean <- lapply(mean_mean,function(x){
    x[is.na(x)] <- 0
    return(x)
  })
  
  mean_mean <- Matrix::Matrix(as.matrix(Reduce("+",mean_mean)))/(length(compare_output)-num_nas)
  mean_mean <- t(round(as.matrix(mean_mean),3))
  mean_mean <- as.data.frame(mean_mean)
  mean_mean$true <- round(true_effects,3)
  
  # percent of bayesian credible intervals that contain the true posterior mean:
  # note need to change index to 3 in the CIs cbind line
  coverage <- data.frame()
  for(i in 1:length(compare_output[[1]]$CIs)){
    tmp <- lapply(compare_output,
                  FUN = function(z){
                    # get the ith result
                    x = z$CIs
                    CIs <- x[[i]]
                    # add the truth
                    CIs <- cbind(CIs,true_effects)
                    coverage = apply(CIs,1,function(x){median(x,na.rm = F) == x[3]})*1
                    coverage[is.na(coverage)] <- 0 # get no coverage if returns NA
                    return(coverage)
                  })
    tmp <- do.call(cbind,tmp)
    tmp <- apply(tmp,1,mean,na.rm =T)
    tmp <- data.frame(t(tmp))
    
    #Make sure names are okay
    if(i!=1){names(tmp) <- names(coverage)}
    coverage <- rbind(coverage,tmp)
  }
  coverage <- t(coverage)
  colnames(coverage) <- names(mean_mean[1:(dim(mean_mean)[2] -1)])
  
  return(list(mean_mean = mean_mean,
              coverage = coverage))
}

# Bayes results that make more sense though are harder to understand
bayes_results  <- function(compare_output,ref_dist = 1){
  
  # remove any NAs
  keep <- sapply(compare_output,function(y){sum(is.na(y$means)) != length(y$means[[1]])})
  
  KL_dists <- data.frame()
  for(i in 1:length(compare_output[[1]]$CIs)){
    tmp <- lapply(compare_output[keep],function(x){
      dists <- sapply(x$KL_dists,function(y){y[ref_dist,i]})
      return(dists)})
    tmp <- apply(do.call(rbind,tmp),2,mean,na.rm = T)
    KL_dists <- rbind(KL_dists,tmp)
  }
  
  names(KL_dists) <- names(compare_output[[1]]$means)
  rownames(KL_dists) <- rownames(compare_output[[1]]$means)
  return(as.data.frame(t(KL_dists)))
}

# Check how long it takes !
t <- proc.time()
tmp <- compare_posteriors(j = 98,
                          ERNM_models = ERNM_models,
                          MRF_models = MRF_models,
                          ERGM_models = ERGM_logistic_models,
                          logistic_models = logistic_models,
                          ERNM_formulas = lapply(1:length(ERNM_models),function(x){DGP_formula}),
                          MRF_formulas = lapply(1:length(MRF_models),function(x){MRF_formula}),
                          ERGM_formulas = lapply(1:length(ERGM_logistic_models),function(x){ergm_formula}),
                          ERGM_nsims = 10,
                          ERNM_nsims = 10,
                          ERNM_draws = 10,
                          ERGM_draws = 10,
                          glm_theta_samples = 1000,
                          treatment_var = "treatment",
                          outcome_var = "outcome",
                          node_treatment_statuses = c(0,1),
                          treatment_indicator = 0,
                          model_names = c("ERNM","MRF","ERGM","logistic"),
                          small = F,
                          keep = NULL,
                          #reldist_method = "gam"
                          reldist_method = "quick"
)
print(proc.time() - t)

# note ERGM nsims must be kept low for memroy considerations

# remove any networks that we were too slow to fit on
keep <- 1:length(DGP_sims_R)
keep <- keep[!is.na(ERNM_models) & !is.na(MRF_models) & !is.na(ERGM_logistic_models)]

t <- proc.time()
cl = parallel::makeForkCluster(10)
compare_results <- parLapply(keep,
                             fun = compare_posteriors,
                             cl = cl,
                             ERNM_models = ERNM_models[keep],
                             MRF_models = MRF_models[keep],
                             ERGM_models = ERGM_logistic_models[keep],
                             logistic_models = logistic_models[keep],
                             ERNM_formulas = lapply(1:length(keep),function(x){DGP_formula}),
                             MRF_formulas = lapply(1:length(keep),function(x){MRF_formula}),
                             ERGM_formulas = lapply(1:length(keep),function(x){ergm_formula}),
                             ERGM_nsims = 50,
                             ERNM_nsims = 100,
                             ERNM_draws = 100,
                             ERGM_draws = 50,
                             glm_theta_samples = 1000,
                             treatment_var = "treatment",
                             outcome_var = "outcome",
                             node_treatment_statuses = c(0,1),
                             treatment_indicator = 0,
                             model_names = c("ERNM","MRF","ERGM","logistic"),
                             small = F,
                             keep = NULL,
                             reldist_method = "gam")
stopCluster(cl)
rm(cl)
print("Comparing the posteriors took this long:")
print(proc.time()-t)
# comparing the posteriors with 10 cores takes ~ 20 hours - need o use only 20 cores due to memory constraints


f_results = freq_results(compare_results,unlist(do.call(c,DGP_effects_both)))
b_results = bayes_results(compare_results,ref_dist = 1)

# =====================
#  Table for paper
# ====================

mean_mean = f_results$mean_mean
mean_mean[mean_mean==0] <- NA
mean_mean <- as.data.frame(as.matrix(mean_mean))
tmp <- round(mean_mean,2)
tmp <- tmp[,c(5,1,2,3,4)]
coverage <- f_results$coverage

#add the coverages:
for(j in 2:dim(tmp)[2]){
  for(i in 1:dim(tmp)[1]){
    if(is.na(tmp[i,j])){
      tmp[i,j] <- as.character(tmp[i,j])
      next
    }
    tmp[i,j] <- paste(as.character(tmp[i,j])," (", as.character(round(as.numeric(coverage[i,(j-1)])*100,0)),"\\","%)",sep ="")
  }
}

tmp <- tmp[c(1:6,12,13,14,15,16),]

rownames(tmp) <- c("main",
                   "1-peer-out","2-peer-out","3-peer-out","4-peer-out","5-peer-out",
                   "1-peer-treat","2-peer-treat","3-peer-treat","4-peer-treat","5-peer-treat"
)
tmp <- as.data.frame(tmp)
names(tmp) <- c("True","ERNM","MRF","ERGM+Logistic","Logistic")

#======================
# Save the results
#======================

# Save the results - i.e. without all the models
setwd("results/simulation_study")
tmp <- c("f_results","b_results","compare_results")
save(file = "results_5_small.RData",list = tmp)
#Save the full_results i.e. with all the models
tmp <- c("ERNM_models","MRF_models","ERGM_logistic_models","logistic_models")
save(file = "results_5.RData",list = tmp)

for(i in ls()){
  size = object.size(eval(parse(text = i)))*(10**(-6))
  if(size >10){
    print(i)
    print(size)
  }
}





