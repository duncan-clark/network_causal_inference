# Modified 
# BERNM model = fire model

fire <- function(formula,
                 prior_dens,
                 run_type = "main", #either "main" or "pilot", if pilot will run chains together to give estimate of burnin time
                 burn_in = 100,
                 samples = 1000,
                 ernm_burnin = 100000,
                 ernm_steps = 10000,
                 n_sims = 10, #number of sims used to estimate Fisher mat
                 gamma_walk = 0.1, #multiplier for random walk portion of proposal step
                 theta0 = NULL,
                 pilot_chains = 8,
                 verbose  = TRUE,
                 cores = NULL,
                 fixed_edges = F,
                 fixed_nodes = F,
                 small_results = T){
  
  if(run_type == "pilot"){
    verbose = F
    small_results = F
  }
  
  y <- ergm.getnetwork(formula)
  model <- ernm(formula,fitModel = F)
  
  t_start <- proc.time()
  
  if (is.null(theta0)) {
    theta0 <- ernm(formula,fitModel = F)$thetas()
  }
  
  model_dim <- length(theta0)

  do_chain <- function(theta){
    
    theta0 <- theta
    theta1 <- rep(NA, length(theta0))
    theta <- matrix(rep(NA,length(theta0)*(burn_in+samples)),
                    byrow = T,
                    ncol = length(theta0))
    
    iters <- burn_in + samples
    # keeps track of current point in each chain, each row is a theta
    # keeps track of whole chain
    
    t_start <- Sys.time()
    t_sims <- proc.time()-proc.time()
    
    
    # EXPERIMENTAL FOR MRF
    if(fixed_edges){
      tmp_model <- ernm(formula,,fitModel = F,nodeSamplingPercentage = 1)
    }else if(fixed_nodes){
      tmp_model <- ernm(formula,,fitModel = F,nodeSamplingPercentage = 0)
    }
    else{
      tmp_model <- ernm(formula,,fitModel = F)
    }
    tmp_model$setThetas(theta0)
    stats <- ernm::calculateStatistics(formula)
    # sims does 10000 burnin, only 10 between samples
    # only use one sample for actual auxillary sampling
    # others are used for estimating fisher info matrix
    sims <- tmp_model$sampler$generateSampleStatistics(ernm_burnin,ernm_steps,n_sims)
    # keep track of acceptances:
    probs <- list()
    length(probs) <- iters
    for(k in 1:iters){
      # Use the previous simulations to approximate the fisher information matrix
      # to give a better sample step as per RMHMCMC
      # demean the sims first:
      sims_demean <- scale(sims,center = T,scale =F)
      info_estim <- (1/n_sims)*(t(sims_demean)%*%sims_demean)
      
      if (abs(det(info_estim)) < 10**(-6)) {
        for (i in 10**seq(-6, 20, 1)) {
          if (abs(det(info_estim)) > 10**(-6)) {
            break
          }
          else {
            tmp <- diag(info_estim)
            tmp <- min(tmp[tmp!=0])
            info_estim <- info_estim + diag(tmp * i, dim(info_estim)[1])
          }
        }
      }
      info_inv <- tryCatch(solve(info_estim),error =function(e){return(diag(1,length(theta0)))})
      # numerical reasons may have made some entries < 0
      info_inv[info_inv<0] <- 0
      #make sure remains symetric
      info_inv[upper.tri(info_inv)] <- t(info_inv)[upper.tri(info_inv)]
      
      theta1 <- theta0 + gamma_walk * rmvnorm(1, sigma = info_inv)[1,]
      
      theta1 <- as.vector(theta1)
      
      #change to sample from ERNM instead of and ERGM
      # delta <- ergm_MCMC_slave(Clist = Clist, proposal = proposal,
      #                          eta = theta1, control = control, verbose = FALSE)$s
      
      tmp_model$setThetas(theta1)
      # generate sample for info estim and the auxillary sampling
      t <- proc.time()
      sims <- tmp_model$sampler$generateSampleStatistics(ernm_burnin,ernm_steps,n_sims)
      t_sims <- t_sims + (proc.time() - t)
      
      sim <- sims[n_sims,] - stats
      
      if(verbose){
        print("old theta is:")
        print(theta0)
        print("proposal is :")
        print(theta1)
        print("Difference between sim and stats is")
        print(sim)
      }

      
      p0 <- prior_dens(theta0)
      p1 <- prior_dens(theta1)
      
      beta <- (theta0 - theta1) %*% sim
      beta - as.numeric(beta[1,1])
      if(p1 !=0){
        log_prior_ratio <- log(p1)- log(p0)
        probs[[k]] <- min(exp(beta + log_prior_ratio),1)
      }else{
        probs[[k]] <- 0
      }

      
      if(verbose){
        print("probability of step is :")
        print(exp(beta))
        print('Prior on old step is:')
        print(p0)
        print('Prior on new step is:')
        print(p1)
        print("iter is:")
        print(k)
      }

      if(runif(1) <= probs[[k]]){
        theta[k,] <- theta1}
      else{
        theta[k,] <- theta0
      }
      theta0 <- theta[k,]
    }
    t_end <- Sys.time()
    if(small_results){
      out = list(time = t_end - t_start,
                 time_sims = t_sims,
                 specs = tmp_model$coef.names,
                 theta_sample = theta[burn_in:(burn_in + samples),])
    }else{
      out = list(time = t_end - t_start,
                 time_sims = t_sims,
                 specs = tmp_model$coef.names,
                 theta_sample = theta[burn_in:(burn_in + samples),],
                 theta_burnin = theta[1:burn_in,],
                 probs = probs)
    }

    class(out) <- "fire"
    
    return(out)
  }
  # debug(do_chain)
  
  if(run_type == "main"){
    out <- do_chain(theta0)
  }
  if(run_type == "pilot"){
    verbose = F
    small_results = F
    # estimate the fisher info
    if(fixed_edges){
      tmp_model <- ernm(formula,fitModel = F,nodeSamplingPercentage = 1)
    }else if(fixed_nodes){
      tmp_model <- ernm(formula,fitModel = F,nodeSamplingPercentage = 1)
    }
    else{
      tmp_model <- ernm(formula,fitModel = F)
    }
    tmp_model$setThetas(theta0)
    stats <- ernm::calculateStatistics(formula)
    # sims does 10000 burnin, only 10 between samples
    # only use one sample for actual auxillary sampling
    # others are used for estimating fisher info matrix
    sims <- tmp_model$sampler$generateSampleStatistics(ernm_burnin,ernm_steps,n_sims)
    sims_demean <- scale(sims,center = T,scale =F)
    info_estim <- (1/n_sims)*(t(sims_demean)%*%sims_demean)
    
    if (abs(det(info_estim)) < 10**(-6)) {
      for (i in 10**seq(-6, 20, 1)) {
        if (abs(det(info_estim)) > 10**(-6)) {
          break
        }
        else {
          tmp <- diag(info_estim)
          tmp <- min(tmp[tmp!=0])
          info_estim <- info_estim + diag(tmp * i, dim(info_estim)[1])
        }
      }
    }
    info_inv <- tryCatch(solve(info_estim),error =function(e){return(diag(1,length(theta0)))})
    # numerical reasons may have made some entries < 0
    info_inv[info_inv<0] <- 0
    #make sure remains symetric
    info_inv[upper.tri(info_inv)] <- t(info_inv)[upper.tri(info_inv)]
    
    thetas <- rmvnorm(n=pilot_chains,mean = theta0, sigma = info_inv)
    
    if(!is.null(cores)){
      cl <- parallel::makeCluster(cores,type = "FORK")
      parallel::clusterSetRNGStream(cl, iseed = 1)
      parallel::clusterEvalQ(cl,{library(network)})
      parallel::clusterEvalQ(cl,{library(ernm)})
      parallel::clusterEvalQ(cl,{library(mvtnorm)})
      parallel::clusterExport(cl,c("do_chain",as.character(formula[[2]]),"burn_in","samples","n_sims","ernm"),envir = environment())
      chains <- parallel::parApply(cl = cl, thetas, 1, FUN = do_chain)
      parallel::stopCluster(cl)
    }else{
      chains <- apply(thetas, 1, FUN = do_chain)
    }

    # do the Gelman stuff to judge the convergence:
    MC_chains <- lapply(chains,function(c){
     as.mcmc(c$theta_burnin)
    })
      
    # Calc median acceptance prob
    median_acceptance <- lapply(chains, function(c) {
      median(unlist(c$probs))
    })
    
    # plot the Gelman diagnostics
    len <- length(MC_chains[[1]][[1]]) / pilot_chains
    
    diagnostics_list <- list()
    length(diagnostics_list) <- model_dim
    plot_list <- list()
    length(plot_list) <- model_dim
    
    for (i in 1:model_dim){
      MCMC_list <- lapply(MC_chains,function(c){
        coda::as.mcmc(apply(c,1,FUN = function(x){x[i]}))
      })
      
      MCMC_list <- mcmc.list(MCMC_list)
      plot(MCMC_list,main = paste("MCMC Plot for",names(theta0)[i]))
      diagnostics <- gelman.diag(MCMC_list)
      diagnostics_list[[i]] <- diagnostics
      plot <- gelman.plot(MCMC_list,main = paste("Shrink Factor Plot for",names(theta0)[i]))
      plot_list[[i]] <- plot
    }
    
    out <- list(MC_chains = MC_chains, 
                chains = chains,
                median_acceptance = median_acceptance,
                diagnostics = diagnostics_list,
                thetas_init = thetas,
                time = proc.time() - t_start,
                fixed_edges = fixed_edges,
                fixed_nodes = fixed_nodes)
  }
  return(out)
}

 