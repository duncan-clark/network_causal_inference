# Contains utilities 

convert_factor <- function(x){as.character(levels(x)[x])}

# function to convert vector of two distinct numbers to binary
binary_convert <- function(x){
  if(length(unique(x) <=2)){
    max <- max(x)
    ones <- (x==max)
    zeros <- (x!=max)
    x[ones] <- 1
    x[zeros] <- 0
    return(x)
  }else{
    warning("Binary convert did not run since there are more than 2 categories")
    return(x)
  }
}

binary_convert <- function(x){
  (x==max(x))*1
}

#Function to convert ERNM binaryNet to network object since ERNM does not have that functionality
binaryNet_to_network <- function(x){
  if(class(x) == "Rcpp_UndirectedNet"){
    el <- x$edges()
    attr(el, "n") <- n <- x$size()
    
    if (nrow(el) > 0)
      nw <- network(el,matrix.type = "edgelist", directed = FALSE)
    else
      nw <- network.initialize(n, directed = FALSE)
    
    for (i in which(x$nMissing(1:n) > 0)) {
      nas <- which(is.na(x[i, 1:n]))
      nw[i, nas] <- NA
    }
    
    vn <- x$variableNames(TRUE)
    if (length(vn) > 0) {
      for (i in 1:length(vn)) {
        vals <- x[[vn[i]]]
        if (vn[i] == "vertex.names") {
          network.vertex.names(nw) <- as.character(vals)
        } else if (vn[i] == "na") {
          
        } else{
          nw %v% vn[i] <-
            if (is.factor(vals))
              as.character(vals)
          else
            as.vector(vals)
        }
      }
    }
  }
  else if(class(x) == "Rcpp_DirectedNet"){
    el <- x$edges()
    attr(el, "n") <- n <- x$size()
    
    if (nrow(el) > 0)
      nw <- network(el, directed = TRUE,matrix.type = "edgelist")
    else
      nw <- network.initialize(n, directed = TRUE)
    
    for (i in which(x$nMissing(1:n) > 0)) {
      nas <- which(is.na(x[i, 1:n]))
      nw[i, nas] <- NA
    }
    
    vn <- x$variableNames(TRUE)
    if (length(vn) > 0) {
      for (i in 1:length(vn)) {
        vals <- x[[vn[i]]]
        if (vn[i] == "vertex.names") {
          network.vertex.names(nw) <- as.character(vals)
        } else if (vn[i] == "na") {
          
        } else{
          nw %v% vn[i] <-
            if (is.factor(vals))
              as.character(vals)
          else
            as.vector(vals)
        }
      }
    }
  }
  
  return(nw)
}

#Function to add the number of treated neighbors to network or binary net object
add_treated_neighs <- function(net,
                               treatment_var){
  tmp <- as.numeric(get.vertex.attribute(net,treatment_var))
  set.vertex.attribute(net,paste(treatment_var,"_neighbors",sep=""),as.character(sapply(1:(net%n%'n'),function(i){sum(tmp[get.neighborhood(net,i)])})))
  return(net)
}

add_treated_neighs_ind <- function(net,
                               treatment_var,
                               treatment_indicator =1){
  tmp <- (as.numeric(get.vertex.attribute(net,treatment_var)) ==treatment_indicator)*1 
  set.vertex.attribute(net,paste(treatment_var,"_neighbors",sep=""),as.character(sapply(1:(net%n%'n'),function(i){sum(tmp[get.neighborhood(net,i)])})))
  return(net)
}


add_treated_neighs_binaryNet <- function(net,treatment_var){
  node_treats <- binary_convert(as.numeric(net$getVariable(treatment_var)))
  neighbor_treats <- sapply(1:net$size(),function(i){
    tmp <- net$neighbors(i)[[1]]
    return(sum(node_treats[tmp]))
  })
  net[[paste(treatment_var,"_neighbors",sep="")]] <- neighbor_treats
  return(net)
}


# function to use simulations from a network model
# to infer distribution of missing potential outcomes
# so we can infer peers effects, overall effects etc
# this is P.O inference in the Bayesian sense:

PO_infer <- function(sims, #sims from model (should be bayesian posterior sims)
                     outcome_var, #needs to be coded 0/1
                     treatment_var, #needs to be coded 0/1 
                     treatment_indicator =1, #whether to use 1 or 0 as treatment indicator
                     num_peers = NULL # if null will calc P.O for ATE, else will calc for P.O for k peer effect of treatment as per Toulis 2013 
){
  # for each simulation for each vertex record treatment and outcome
  n <- get.network.attribute(sims[[1]],"n")
  
  # if we are doing peer treatment we are interested
  # only in nodes that have exactly k treated neighbors and are not treated themselves
  # control is nodes that have no treated neighbors and no self treatment
  if(is.null(num_peers)){
    treatments <- lapply(sims,function(net){
      (as.numeric(get.vertex.attribute(net,treatment_var)) == treatment_indicator)*1
    })
  }else{
    treatments <- lapply(sims,function(net){
      # get the treatment of the nodes
      node_treats <- (as.numeric(get.vertex.attribute(net,treatment_var)) == treatment_indicator)*1
      # get the numberof treated neighbors
      neighbor_treats <- sapply(1:n,function(i){
        tmp <- get.neighborhood(net,i)
        return(sum(node_treats[tmp]))
      })
      
      # return 1 if correct treatment
      # return 0 if correct control treatment
      # else return NA
      # if treatment variable is the same outcome variable, don't require 0 treatment on node
      if(outcome_var != treatment_var){
        treats <- (node_treats==0)*(neighbor_treats==num_peers)
        controls <- (node_treats==0)*(neighbor_treats ==0)
      }else{
        treats <- (neighbor_treats==num_peers)
        controls <- (neighbor_treats ==0)
      }
      
      #controls <- (node_treats==0)*(neighbor_treats ==0)
      others <- (treats ==0)*(controls==0)
      
      result = rep(0,n)
      
      result[as.logical(treats)] <- 1
      result[as.logical(controls)] <- 0 
      result[as.logical(others)] <- NA
      
      return(result)
      
    })
  }
  
  outcomes <- lapply(sims,function(net){
    as.numeric(get.vertex.attribute(net,outcome_var))
  })
  
  ATEs <- mapply(treatments,outcomes,FUN = function(x,y){
    treated <- which(x == 1)
    controls <- which(x == 0)
    if(length(treated) !=0 & length(controls) !=0){
      mean_outcomes_treated <- mean(y[treated])
      mean_outcomes_controls <- mean(y[controls])
      return(mean_outcomes_treated - mean_outcomes_controls)
    }else{return(NA)}
  })
  
  # Extract our control and treated mean outcomes
  # NOTE NO KIND OF MATCHING ETC IS DONE HERE
  # combine the treaments and outcomes, since network structure not important:
  treatments <- do.call(c,treatments)
  outcomes <- do.call(c,outcomes)
  treated <- which(treatments == 1)
  controls <- which(treatments ==0 )
  mean_outcomes_treated <- mean(outcomes[treated])
  mean_outcomes_controls <- mean(outcomes[controls])
  
  return(list(ATE = mean_outcomes_treated - mean_outcomes_controls,
              ATEs = ATEs))
}

PO_infer_binaryNet <- function(sims, #should be in binary net format
                                outcome_var,
                                treatment_var,
                                treatment_indicator = 1, #whether to use 1 or 0 as treatment indicator
                                num_peers = NULL
){
  
  # for each simulation for each vertex record treatment and outcome
  n <- sims[[1]]$size()
  
  # if we are doing peer treatment we are interested
  # only in nodes that have exactly k treated neighbors and are not treated themselves
  # control is nodes that have no treated neighbors and no self treatment
  if(is.null(num_peers)){
    treatments <- lapply(sims,function(net){
      (as.numeric(convert_factor(net$getVariable(treatment_var))) == treatment_indicator)*1
    })
  }else{
    treatments <- lapply(sims,function(net){
      # get the treatment of the nodes
      node_treats <- (as.numeric(convert_factor(net$getVariable(treatment_var))) == treatment_indicator)*1
      # get the number of treated neighbors
      neighbor_treats <- sapply(1:n,function(i){
        tmp <- net$neighbors(i)[[1]]
        return(sum(node_treats[tmp]))
      })
      
      # return 1 if correct treatment
      # return 0 if correct control treatment
      # else return NA
      # if treatment variable is the same outcome variable, don't require 0 treatment on node
      if(outcome_var != treatment_var){
        treats <- (node_treats==0)*(neighbor_treats==num_peers)
        controls <- (node_treats==0)*(neighbor_treats ==0)
      }else{
        treats <- (neighbor_treats==num_peers)
        controls <- (neighbor_treats ==0)
      }
      
      #controls <- (node_treats==0)*(neighbor_treats ==0)
      others <- (treats ==0)*(controls==0)
      
      result = rep(0,n)
      
      result[as.logical(treats)] <- 1
      result[as.logical(controls)] <- 0 
      result[as.logical(others)] <- NA
      
      return(result)
      
    })
  }
  
  outcomes <- lapply(sims,function(net){
    binary_convert(as.numeric(convert_factor(net$getVariable(outcome_var))))
  })
  
  ATEs <- mapply(treatments,outcomes,FUN = function(x,y){
    treated <- which(x == 1)
    controls <- which(x == 0)
    if(length(treated) !=0 & length(controls) !=0){
      mean_outcomes_treated <- mean(y[treated])
      mean_outcomes_controls <- mean(y[controls])
      return(mean_outcomes_treated - mean_outcomes_controls)
    }else{return(NA)}
  })
  
  # Extract our control and treated mean outcomes
  # NOTE NO KIND OF MATCHING ETC IS DONE HERE
  # combine the treaments and outcomes, since network structure not important:
  treatments <- do.call(c,treatments)
  outcomes <- do.call(c,outcomes)
  treated <- which(treatments == 1)
  controls <- which(treatments == 0)
  mean_outcomes_treated <- mean(outcomes[treated])
  mean_outcomes_controls <- mean(outcomes[controls])
  
  
  
  return(list(treated_outcomes = outcomes[treated],
              control_outcomes = outcomes[controls],
              ATE = mean(outcomes[treated]) - mean(outcomes[controls]),
              net_ATEs = ATEs))
}

# Special internal function for the extract posterior function
# Idea is to avoid accsssing things for the same sim over and over again to save time
PO_infer_fast <- function(sims, #should be in binary net format
                          outcome_var,
                          treatment_var,
                          treatment_indicator = 1, #whether to use 1 or 0 as treatment indicator
                          num_peers = NULL,
                          neighbor_treats, #list of vectors of sums of neighborhood treatments - so can be reused
                          node_treatment_statuses = c(0,1) # the statuses of the nodes to be considered
                          ){  #for each simulation for each vertex record treatment and outcome
  n <- sims[[1]]$size()
  
  # if we are doing peer treatment we are interested
  # only in nodes that have exactly k treated neighbors and are not treated themselves
  # control is nodes that have no treated neighbors and no self treatment
  if(is.null(num_peers)){
    treatments <- lapply(sims,function(net){
      node_treats = (as.numeric(convert_factor(net$getVariable(treatment_var))) == treatment_indicator)*1
    })
  }else{
    treatments <- mapply(sims,neighbor_treats,FUN = function(net,neighbor_treats){
      # get the treatment of the nodes
      node_treats <- (as.numeric(convert_factor(net$getVariable(treatment_var))) == treatment_indicator)*1
      # get the number of treated neighbors
      # no longer needed as is supplied
      # neighbor_treats <- sapply(1:n,function(i){
      #   tmp <- net$neighbors(i)[[1]]
      #   return(sum(node_treats[tmp]))
      # })
      
      # return 1 if correct treatment
      # return 0 if correct control treatment
      # else return NA
      # if treatment variable is the same outcome variable, don't require 0 treatment on node
      if(outcome_var != treatment_var){
        treats <- (node_treats %in% node_treatment_statuses)*(neighbor_treats == num_peers)
        controls <- (node_treats %in% node_treatment_statuses)*(neighbor_treats == 0)
      }else{
        treats <- (neighbor_treats == num_peers)
        controls <- (neighbor_treats == 0)
      }
      
      #controls <- (node_treats==0)*(neighbor_treats ==0)
      others <- (treats ==0)*(controls==0)
      
      result = rep(0,n)
      
      result[as.logical(treats)] <- 1
      result[as.logical(controls)] <- 0 
      result[as.logical(others)] <- NA
      
      return(result)
      
    },SIMPLIFY = F)
  }
  
  outcomes <- lapply(sims,function(net){
    binary_convert(as.numeric(convert_factor(net$getVariable(outcome_var))))
  })
  
  ATEs <- mapply(treatments,outcomes,FUN = function(x,y){
    treated <- which(x == 1)
    controls <- which(x == 0)
    if(length(treated) !=0 & length(controls) !=0){
      mean_outcomes_treated <- mean(y[treated])
      mean_outcomes_controls <- mean(y[controls])
      return(mean_outcomes_treated - mean_outcomes_controls)
    }else{return(NA)}
  })
  
  # Extract our control and treated mean outcomes
  # NOTE NO KIND OF MATCHING ETC IS DONE HERE
  # combine the treaments and outcomes, since network structure not important:
  treatments <- do.call(c,treatments)
  outcomes <- do.call(c,outcomes)
  treated <- which(treatments == 1)
  controls <- which(treatments == 0)
  mean_outcomes_treated <- mean(outcomes[treated])
  mean_outcomes_controls <- mean(outcomes[controls])
  
  
  
  return(list(treated_outcomes = outcomes[treated],
              control_outcomes = outcomes[controls],
              ATE = mean(outcomes[treated]) - mean(outcomes[controls]),
              net_ATEs = ATEs))
}

# Idea is to estimate the effect for each node and then take the mean of this,
# not the mean of all nodes that achieve the treatment
# this effectively weights all nodes evenly, whereas the non test version effectively weights, by how often a node received treatment
# this is the "correct" version
PO_infer_fast_ind <- function(sims, #should be in binary net format
                          outcome_var,
                          treatment_var,
                          treatment_indicator = 1, # whether to use 1 or 0 as treatment indicator
                          num_peers = NULL,
                          neighbor_treats, # list of vectors of sums of neighborhood treatments - so can be reused
                          node_treatment_statuses = c(0,1) # the statuses of the nodes to be considered
){  #for each simulation for each vertex record treatment and outcome
  n <- sims[[1]]$size()
  
  # if we are doing peer treatment we are interested
  # only in nodes that have exactly k treated neighbors and are not treated themselves
  # control is nodes that have no treated neighbors and no self treatment
  if(is.null(num_peers)){
    treatments <- lapply(sims,function(net){
      node_treats = (as.numeric(convert_factor(net$getVariable(treatment_var))) == treatment_indicator)*1
    })
  }else{
    treatments <- mapply(sims,neighbor_treats,FUN = function(net,neighbor_treats){
      # get the treatment of the nodes
      node_treats <- (as.numeric(convert_factor(net$getVariable(treatment_var))) == treatment_indicator)*1
      # get the number of treated neighbors
      # no longer needed as is supplied
      # neighbor_treats <- sapply(1:n,function(i){
      #   tmp <- net$neighbors(i)[[1]]
      #   return(sum(node_treats[tmp]))
      # })
      
      # return 1 if correct treatment
      # return 0 if correct control treatment
      # else return NA
      # if treatment variable is the same outcome variable, don't require 0 treatment on node
      if(outcome_var != treatment_var){
        treats <- (node_treats %in% node_treatment_statuses)*(neighbor_treats == num_peers)
        controls <- (node_treats %in% node_treatment_statuses)*(neighbor_treats == 0)
      }else{
        treats <- (neighbor_treats == num_peers)
        controls <- (neighbor_treats == 0)
      }
      
      #controls <- (node_treats==0)*(neighbor_treats ==0)
      others <- (treats ==0)*(controls==0)
      
      result = rep(0,n)
      
      result[as.logical(treats)] <- 1
      result[as.logical(controls)] <- 0 
      result[as.logical(others)] <- NA
      
      return(result)
      
    },SIMPLIFY = F)
  }
  
  outcomes <- lapply(sims,function(net){
    binary_convert(as.numeric(convert_factor(net$getVariable(outcome_var))))
  })
  
  # node wise treatment means
  treat_means <- mapply(treatments,outcomes,FUN = function(x,y){
    tmp = y
    tmp[x !=1 | is.na(x)] = NA 
    return(tmp)
  },SIMPLIFY = F)
  treat_means = apply(do.call(rbind,treat_means),2,mean,na.rm = T)
  
  #nodewise control means
  control_means <- mapply(treatments,outcomes,FUN = function(x,y){
    tmp = y
    tmp[x !=0 | is.na(x)] = NA 
    return(tmp)
  },SIMPLIFY = F)
  control_means = apply(do.call(rbind,control_means),2,FUN = mean,na.rm = T)

  ATE = mean(treat_means - control_means,na.rm = T)
  
  return(list(control_means = control_means,
              treat_means = treat_means,
              ATE = ATE))
}


# Function to extract all effects from simulated networks
extract_effects_from_sims <- function(sims, #should be in binary net format
                                      outcome_var,
                                      treatment_var,
                                      treatment_indicator = 1, #whether to use 1 or 0 as treatment indicator
                                      node_treatment_statuses = c(0,1) # the statuses of the nodes to be considered
){
  n <- sims[[1]]$size()
  neighbor_treats <- lapply(sims,function(net){
    node_treats <- (as.numeric(convert_factor(net$getVariable(treatment_var))) == treatment_indicator)*1
    return(sapply(1:n,function(i){
      tmp <- net$neighbors(i)[[1]]
      return(sum(node_treats[tmp]))}))
  })
  neighbor_outs <- lapply(sims,function(net){
    node_outs <- binary_convert(as.numeric(net$getVariable(outcome_var)))
    return(sapply(1:n,function(i){
      tmp <- net$neighbors(i)[[1]]
      return(sum(node_outs[tmp]))}))
  })
  
  main_effect <- PO_infer_binaryNet(sims = sims,
                                    outcome_var = outcome_var,
                                    treatment_var = treatment_var,
                                    treatment_indicator = treatment_indicator,
                                    num_peers = NULL)
  
  k_peer_out <- list()
  length(k_peer_out) <- 10
  for(k in 1:10){
    k_peer_out[[k]] <- PO_infer_fast(sims = sims,
                                     outcome_var = outcome_var,
                                     treatment_var = outcome_var,
                                     neighbor_treats = neighbor_outs,
                                     treatment_indicator = 1,
                                     num_peers = k)
  }
  
  k_peer_treat <- list()
  length(k_peer_treat) <- 10
  for(k in 1:10){
    #if(k==2){browser()}
    k_peer_treat[[k]] <- PO_infer_fast(sims = sims,
                                       outcome_var = outcome_var,
                                       treatment_var = treatment_var,
                                       treatment_indicator = treatment_indicator,
                                       neighbor_treats = neighbor_treats,
                                       node_treatment_statuses = node_treatment_statuses,
                                       num_peers = k)
  }
  #return the estimated effects:
  return(list(main = main_effect$ATE,
              k_peer_out = lapply(k_peer_out,function(x){x$ATE}),
              k_peer_treat = lapply(k_peer_treat,function(x){x$ATE}))
  )
}

# function to calculate the effect on the response scale for logistic regression 
# beta should be a row vector
ATE_logistic_func <- function(beta,x_0,x_1){
  mean((1/(1+exp(-x_1%*%t(beta)))) - (1/(1+exp(-x_0%*%t(beta)))))
}

# Some utilities to extract from the fire models the things that we want
# posterior functions take in the output of the extract function

extract_posteriors <- function(model,theta_draws=100,fixed_edges = F,nsims,treatment_indicator = 1){
  #simulate from the posterior to get network that can be used to get the posterior distribution of the
  if(fixed_edges){
    tmp_model <- ernm(model$formula,
                      maxIter = 2,
                      mcmcSampleSize = 100,
                      mcmcBurnIn = 100,
                      nodeSamplingPercentage = 1)
  }else{
    tmp_model <- ernm(model$formula,
                      maxIter = 2,
                      mcmcSampleSize = 100,
                      mcmcBurnIn = 100)
  }
  
  # select the required number of draws:
  samp <- model$theta_sample[sample(dim(model$theta_sample)[1],theta_draws,replace = T),]
  
  estims <- lapply(1:dim(samp)[1],function(i){
    theta <- samp[i,]
    tmp_model$m$setThetas(theta)
    sims <- tmp_model$m$sampler$generateSample(10000,1000,nsims)
    
    main_effect <- PO_infer_binaryNet(sims = sims,
                                      outcome_var = "outcome",
                                      treatment_var = "treatment",
                                      treatment_indicator = treatment_indicator,
                                      num_peers = NULL)
    
    k_peer_out <- list()
    length(k_peer_out) <- 10
    for(k in 1:10){
      k_peer_out[[k]] <- PO_infer_binaryNet(sims = sims,
                                            outcome_var = "outcome",
                                            treatment_var = "outcome",
                                            num_peers = k)
    }
    
    k_peer_treat <- list()
    length(k_peer_treat) <- 10
    for(k in 1:10){
      k_peer_treat[[k]] <- PO_infer_binaryNet(sims = sims,
                                              outcome_var = "outcome",
                                              treatment_var = "treatment",
                                              treatment_indicator = treatment_indicator,
                                              num_peers = k)
    }
    
    #return the estimated effects:
    return(list(main = main_effect$ATE,
                k_peer_out = lapply(k_peer_out,function(x){x$ATE}),
                k_peer_treat = lapply(k_peer_treat,function(x){x$ATE}))
    )
  }
  )
  result <- list()
  result[[1]] <- sapply(estims,function(x){x$main})
  for(k in 2:11){
    result[[k]] <- sapply(estims,function(x){x$k_peer_out[[(k-1)]]})
  }
  for(k in 12:21){
    result[[k]] <- sapply(estims,function(x){x$k_peer_treat[[(k-11)]]})
  }
  names(result) <- c("main",
                     sapply(seq(1,10),function(i){paste("k_peer_outcome",i,sep="")}),
                     sapply(seq(1,10),function(i){paste("k_peer_treatment",i,sep="")}))
  return(result)
}
posterior_means <- function(extract_posterior_output){
  tmp <- extract_posterior_output
  return(sapply(tmp,mean,na.rm = T))
}
posterior_missingness <- function(extract_posterior_output){
  tmp <- extract_posterior_output
  return(sapply(tmp,function(x){sum(is.na(x))}))
}
posterior_CIs <- function(extract_posterior_output,p){
  tmp <- extract_posterior_output
  result <- lapply(tmp,function(x){
    return(c(quantile(x,p,na.rm = T), quantile(x,(1-p),na.rm = T)))
  })
  result = do.call(rbind,result)
  rownames(result) <- names(tmp)
  return(result)
}

extract_posteriors_fast <- function(model,
                                    formula,
                                    theta_draws=100,
                                    fixed_edges = F,
                                    nsims,
                                    treatment_var = "treatment",
                                    outcome_var ="outcome",
                                    treatment_indicator = 1,
                                    node_treatment_statuses = c(0,1),
                                    cores = NULL,
                                    return_sims = F){
  #simulate from the posterior to get network that can be used to get the posterior distribution of the
  if(fixed_edges){
    tmp_model <- ernm(formula,
                      maxIter = 2,
                      mcmcSampleSize = 100,
                      mcmcBurnIn = 100,
                      nodeSamplingPercentage = 1)
  }else{
    tmp_model <- ernm(formula,
                      maxIter = 2,
                      mcmcSampleSize = 100,
                      mcmcBurnIn = 100)
  }
  
  # select the required number of draws:
  samp <- model$theta_sample[sample(dim(model$theta_sample)[1],theta_draws,replace = T),]
  
  do_estimate = function(i,refit = T,big = return_sims){
    
    if(refit){
      if(fixed_edges){
        tmp_model <- ernm(formula,
                          maxIter = 2,
                          mcmcSampleSize = 100,
                          mcmcBurnIn = 100,
                          nodeSamplingPercentage = 1)
      }else{
        tmp_model <- ernm(formula,
                          maxIter = 2,
                          mcmcSampleSize = 100,
                          mcmcBurnIn = 100)
      }
    }
    theta <- samp[i,]
    tmp_model$m$setThetas(theta)
    sims <- tmp_model$m$sampler$generateSample(10000,1000,nsims)
    
    #get the number of treated neighbors
    n <- sims[[1]]$size()
    neighbor_treats <- lapply(sims,function(net){
      node_treats <- (as.numeric(convert_factor(net$getVariable(treatment_var))) == treatment_indicator)*1
      return(sapply(1:n,function(i){
        tmp <- net$neighbors(i)[[1]]
        return(sum(node_treats[tmp]))}))
    })
    neighbor_outs <- lapply(sims,function(net){
      node_outs <- binary_convert(as.numeric(net$getVariable(outcome_var)))
      return(sapply(1:n,function(i){
        tmp <- net$neighbors(i)[[1]]
        return(sum(node_outs[tmp]))}))
    })
    
    main_effect <- PO_infer_binaryNet(sims = sims,
                                      outcome_var = outcome_var,
                                      treatment_var = treatment_var,
                                      treatment_indicator = treatment_indicator,
                                      num_peers = NULL)
    
    k_peer_out <- list()
    length(k_peer_out) <- 10
    for(k in 1:10){
      k_peer_out[[k]] <- PO_infer_fast(sims = sims,
                                       outcome_var = outcome_var,
                                       treatment_var = outcome_var,
                                       neighbor_treats = neighbor_outs,
                                       treatment_indicator = 1,
                                       num_peers = k)
    }
    
    k_peer_treat <- list()
    length(k_peer_treat) <- 10
    for(k in 1:10){
      k_peer_treat[[k]] <- PO_infer_fast(sims = sims,
                                         outcome_var = outcome_var,
                                         treatment_var = treatment_var,
                                         treatment_indicator = treatment_indicator,
                                         neighbor_treats = neighbor_treats,
                                         node_treatment_statuses = node_treatment_statuses,
                                         num_peers = k)
    }
    if(big){
      return(list(main = main_effect$ATE,
                  k_peer_out = lapply(k_peer_out,function(x){x$ATE}),
                  k_peer_treat = lapply(k_peer_treat,function(x){x$ATE}),
                  sims = sims))}
    else{
      #return the estimated effects:
      return(list(main = main_effect$ATE,
                  k_peer_out = lapply(k_peer_out,function(x){x$ATE}),
                  k_peer_treat = lapply(k_peer_treat,function(x){x$ATE}))
      )
    }
  }
  
  if(!is.null(cores)){
    cl <- parallel::makeCluster(cores,type = "FORK")
    parallel::clusterSetRNGStream(cl, iseed = 1)
    parallel::clusterEvalQ(cl,{library(network)})
    parallel::clusterEvalQ(cl,{library(ernm)})
    parallel::clusterExport(cl,c("samp","PO_infer_binaryNet","PO_infer_fast","convert_factor","fixed_edges"),envir = environment())
    estims <- parallel::parLapply(cl = cl, 1:dim(samp)[1], fun = do_estimate)
    parallel::stopCluster(cl)
  }else{
    estims <- lapply(1:dim(samp)[1],do_estimate,refit = F)
  }

  
  result <- list()
  result[[1]] <- sapply(estims,function(x){x$main})
  for(k in 2:11){
    result[[k]] <- sapply(estims,function(x){x$k_peer_out[[(k-1)]]})
  }
  for(k in 12:21){
    result[[k]] <- sapply(estims,function(x){x$k_peer_treat[[(k-11)]]})
  }
  names(result) <- c("main",
                     sapply(seq(1,10),function(i){paste("k_peer_outcome",i,sep="")}),
                     sapply(seq(1,10),function(i){paste("k_peer_treatment",i,sep="")}))
  
  if(return_sims){
    return(list(result = result,
                sims = do.call(c,lapply(estims,function(x){x$sims}))))
  }
  return(result)
}

# ERGM  + logistic extract the posterior
# - simulate network, refit linear model to new data
# - from linear model directly infer treatment effects by simulation
# - these simulations don't really make sense - should be inferring effects directly from the models
# - but since the models can differ this is tricky...
extract_posteriors_ERGM_fast <- function(model,formula,theta_draws = 100,nsims=100,treatment_var = "treatment",outcome_var ="outcome",treatment_indicator = 1){
  #simulate from the posterior to get network that can be used to get the posterior distribution of the
  tmp_model <- ernm(formula,
                    maxIter = 3,
                    mcmcSampleSize = 100,
                    mcmcBurnIn = 100)
  
  # thin the sample
  samp <- model$ergm$theta_sample[sample(dim(model$ergm$theta_sample)[1],theta_draws,replace = T),]
  
  estims <- lapply(1:dim(samp)[1],function(i){
    theta <- samp[i,]
    tmp_model$m$setThetas(theta)
    sims <- tmp_model$m$sampler$generateSample(10000,1000,nsims)
    sims <- lapply(sims,function(net){
      
      # now for each sim do the glm addition bit:
      net <- add_treated_neighs_binaryNet(net,treatment_var)
      net <- add_treated_neighs_binaryNet(net,outcome_var)
      
      assign(outcome_var,binary_convert(as.numeric(net$getVariable(outcome_var))))
      assign(treatment_var,binary_convert(as.numeric(net$getVariable(treatment_var))))
      assign(paste(treatment_var,"_neighbors",sep=""),as.numeric(net$getVariable(paste(treatment_var,"_neighbors",sep=""))))
      assign(paste(outcome_var,"_neighbors",sep=""),as.numeric(net$getVariable(paste(outcome_var,"_neighbors",sep=""))))
      
      newdata <- list(eval(parse(text = outcome_var)),
                      eval(parse(text = treatment_var)),
                      eval(parse(text = paste(outcome_var,"_neighbors",sep=""))),
                      eval(parse(text = paste(treatment_var,"_neighbors",sep=""))))
      names(newdata) <- c(outcome_var,
                          treatment_var,
                          paste(outcome_var,"_neighbors",sep=""),
                          paste(treatment_var,"_neighbors",sep=""))
      newdata <- as.data.frame(newdata)
      
      net[[outcome_var]] <- as.factor(1*(predict(model$glm,newdata=newdata,type = "response")>0.5))
      return(net)}
    )
    
    # get the number of treated neighbors
    n <- sims[[1]]$size()
    neighbor_treats <- lapply(sims,function(net){
      node_treats <- (as.numeric(convert_factor(net$getVariable(treatment_var))) == treatment_indicator)*1
      return(sapply(1:n,function(i){
        tmp <- net$neighbors(i)[[1]]
        return(sum(node_treats[tmp]))}))
    })
    neighbor_outs <- lapply(sims,function(net){
      node_outs <- binary_convert(as.numeric(net$getVariable(outcome_var)))
      return(sapply(1:n,function(i){
        tmp <- net$neighbors(i)[[1]]
        return(sum(node_outs[tmp]))}))
    })
    main_effect <- PO_infer_binaryNet(sims = sims,
                                      outcome_var = outcome_var,
                                      treatment_var = treatment_var,
                                      treatment_indicator = treatment_indicator,
                                      num_peers = NULL)
    
    k_peer_out <- list()
    length(k_peer_out) <- 10
    for(k in 1:10){
      k_peer_out[[k]] <- PO_infer_fast(sims = sims,
                                       outcome_var = outcome_var,
                                       treatment_var = outcome_var,
                                       neighbor_treats = neighbor_outs,
                                       num_peers = k)
    }
    
    k_peer_treat <- list()
    length(k_peer_treat) <- 10
    for(k in 1:10){
      k_peer_treat[[k]] <- PO_infer_fast(sims = sims,
                                         outcome_var = outcome_var,
                                         treatment_var = treatment_var,
                                         treatment_indicator = treatment_indicator,
                                         neighbor_treats = neighbor_treats,
                                         num_peers = k)
    }
    
    #return the estimated effects:
    return(list(main = main_effect$ATE,
                k_peer_out = lapply(k_peer_out,function(x){x$ATE}),
                k_peer_treat = lapply(k_peer_treat,function(x){x$ATE}))
    )
  }
  )
  result <- list()
  result[[1]] <- sapply(estims,function(x){x$main})
  for(k in 2:11){
    result[[k]] <- sapply(estims,function(x){x$k_peer_out[[(k-1)]]})
  }
  for(k in 12:21){
    result[[k]] <- sapply(estims,function(x){x$k_peer_treat[[(k-11)]]})
  }
  names(result) <- c("main",
                     sapply(seq(1,10),function(i){paste("k_peer_outcome",i,sep="")}),
                     sapply(seq(1,10),function(i){paste("k_peer_treatment",i,sep="")}))
  return(result)
}
extract_posteriors_ERGM_fit <- function(model,formula,theta_draws = 100,nsims=100,treatment_var = "treatment",outcome_var ="outcome",treatment_indicator = 1){
  #simulate from the posterior to get network that can be used to get the posterior distribution of the
  tmp_model <- ernm(formula,
                    maxIter = 3,
                    mcmcSampleSize = 100,
                    mcmcBurnIn = 100)
  
  # thin the sample
  samp <- model$ergm$theta_sample[sample(dim(model$ergm$theta_sample)[1],theta_draws,replace = T),]
  
  estims <- lapply(1:dim(samp)[1],function(i){
    theta <- samp[i,]
    tmp_model$m$setThetas(theta)
    sims <- tmp_model$m$sampler$generateSample(10000,1000,nsims)
    sims <- lapply(sims,function(net){
      
      # now for each sim do the glm addition bit:
      net <- add_treated_neighs_binaryNet(net,treatment_var)
      net <- add_treated_neighs_binaryNet(net,outcome_var)
      
      assign(outcome_var,binary_convert(as.numeric(net$getVariable(outcome_var))))
      assign(treatment_var,binary_convert(as.numeric(net$getVariable(treatment_var))))
      assign(paste(treatment_var,"_neighbors",sep=""),as.numeric(net$getVariable(paste(treatment_var,"_neighbors",sep=""))))
      assign(paste(outcome_var,"_neighbors",sep=""),as.numeric(net$getVariable(paste(outcome_var,"_neighbors",sep=""))))
      
      newdata <- list(eval(parse(text = outcome_var)),
                      eval(parse(text = treatment_var)),
                      eval(parse(text = paste(outcome_var,"_neighbors",sep=""))),
                      eval(parse(text = paste(treatment_var,"_neighbors",sep=""))))
      names(newdata) <- c(outcome_var,
                          treatment_var,
                          paste(outcome_var,"_neighbors",sep=""),
                          paste(treatment_var,"_neighbors",sep=""))
      newdata <- as.data.frame(newdata)
      
      new_model <- model$glm
      new_model <- update(new_model,data=newdata)
      
      net[[outcome_var]] <- as.factor(1*(predict(new_model,newdata=newdata,type = "response")>0.5))
      rm(new_model)
      return(net)}
    )
    
    # get the number of treated neighbors
    n <- sims[[1]]$size()
    neighbor_treats <- lapply(sims,function(net){
      node_treats <- (as.numeric(convert_factor(net$getVariable(treatment_var))) == treatment_indicator)*1
      return(sapply(1:n,function(i){
        tmp <- net$neighbors(i)[[1]]
        return(sum(node_treats[tmp]))}))
    })
    neighbor_outs <- lapply(sims,function(net){
      node_outs <- binary_convert(as.numeric(net$getVariable(outcome_var)))
      return(sapply(1:n,function(i){
        tmp <- net$neighbors(i)[[1]]
        return(sum(node_outs[tmp]))}))
    })
    main_effect <- PO_infer_binaryNet(sims = sims,
                                      outcome_var = outcome_var,
                                      treatment_var = treatment_var,
                                      treatment_indicator = treatment_indicator,
                                      num_peers = NULL)
    
    k_peer_out <- list()
    length(k_peer_out) <- 10
    for(k in 1:10){
      k_peer_out[[k]] <- PO_infer_fast(sims = sims,
                                       outcome_var = outcome_var,
                                       treatment_var = outcome_var,
                                       neighbor_treats = neighbor_outs,
                                       num_peers = k)
    }
    
    k_peer_treat <- list()
    length(k_peer_treat) <- 10
    for(k in 1:10){
      k_peer_treat[[k]] <- PO_infer_fast(sims = sims,
                                         outcome_var = outcome_var,
                                         treatment_var = treatment_var,
                                         treatment_indicator = treatment_indicator,
                                         neighbor_treats = neighbor_treats,
                                         num_peers = k)
    }
    
    #return the estimated effects:
    return(list(main = main_effect$ATE,
                k_peer_out = lapply(k_peer_out,function(x){x$ATE}),
                k_peer_treat = lapply(k_peer_treat,function(x){x$ATE}))
    )
  }
  )
  result <- list()
  result[[1]] <- sapply(estims,function(x){x$main})
  for(k in 2:11){
    result[[k]] <- sapply(estims,function(x){x$k_peer_out[[(k-1)]]})
  }
  for(k in 12:21){
    result[[k]] <- sapply(estims,function(x){x$k_peer_treat[[(k-11)]]})
  }
  names(result) <- c("main",
                     sapply(seq(1,10),function(i){paste("k_peer_outcome",i,sep="")}),
                     sapply(seq(1,10),function(i){paste("k_peer_treatment",i,sep="")}))
  return(result)
}
extract_posteriors_logistic <- function(model,net,theta_draws = 100,nsims=100,treatment_var = "treatment",outcome_var ="outcome",treatment_indicator = 1){
  preds <- posterior_predict(model$glm)
  
  estims <- lapply(1:theta_draws,function(i){
    sims <- list(net)
    sims[[1]][[outcome_var]] <- preds[i,]
    
    # get the number of treated neighbors
    n <- sims[[1]]$size()
    neighbor_treats <- lapply(sims,function(net){
      node_treats <- (as.numeric(convert_factor(net$getVariable(treatment_var))) == treatment_indicator)*1
      return(sapply(1:n,function(i){
        tmp <- net$neighbors(i)[[1]]
        return(sum(node_treats[tmp]))}))
    })
    neighbor_outs <- lapply(sims,function(net){
      node_outs <- binary_convert(as.numeric(net$getVariable(outcome_var)))
      return(sapply(1:n,function(i){
        tmp <- net$neighbors(i)[[1]]
        return(sum(node_outs[tmp]))}))
    })
    main_effect <- PO_infer_binaryNet(sims = sims,
                                      outcome_var = outcome_var,
                                      treatment_var = treatment_var,
                                      treatment_indicator = treatment_indicator,
                                      num_peers = NULL)
    
    k_peer_out <- list()
    length(k_peer_out) <- 10
    for(k in 1:10){
      k_peer_out[[k]] <- PO_infer_fast(sims = sims,
                                       outcome_var = outcome_var,
                                       treatment_var = outcome_var,
                                       neighbor_treats = neighbor_outs,
                                       num_peers = k)
    }
    
    k_peer_treat <- list()
    length(k_peer_treat) <- 10
    for(k in 1:10){
      k_peer_treat[[k]] <- PO_infer_fast(sims = sims,
                                         outcome_var = outcome_var,
                                         treatment_var = treatment_var,
                                         treatment_indicator = treatment_indicator,
                                         neighbor_treats = neighbor_treats,
                                         num_peers = k)
    }
    
    #return the estimated effects:
    return(list(main = main_effect$ATE,
                k_peer_out = lapply(k_peer_out,function(x){x$ATE}),
                k_peer_treat = lapply(k_peer_treat,function(x){x$ATE}))
    )
  }
  )
  result <- list()
  result[[1]] <- sapply(estims,function(x){x$main})
  for(k in 2:11){
    result[[k]] <- sapply(estims,function(x){x$k_peer_out[[(k-1)]]})
  }
  for(k in 12:21){
    result[[k]] <- sapply(estims,function(x){x$k_peer_treat[[(k-11)]]})
  }
  names(result) <- c("main",
                     sapply(seq(1,10),function(i){paste("k_peer_outcome",i,sep="")}),
                     sapply(seq(1,10),function(i){paste("k_peer_treatment",i,sep="")}))
  return(result)
}

extract_posteriors_logistic_directly <- function(model,main_index,k_treat_index,k_out_index,theta_samples = 100,node_treatment_status = NULL){
  
  samp <- as.matrix(model)
  X = model$model
  
  # get the main effect
  x_0 <- as.matrix(cbind(1,X[,-1])) 
  x_1 <- x_0
  x_0[,main_index + 1] <- 0
  x_1[,main_index + 1] <- 1
  
  main <- sapply(sample(1:dim(samp)[1],theta_samples,replace = T),function(i){
    # construct X_0 and X_1 from the model 
    # add intercept and leave of outcome column

    return(ATE_logistic_func(t(as.matrix(samp[i,])),
                             x_0 = x_0,
                             x_1 = x_1))
  })
  
  k_peer_treat <- list()
  length(k_peer_treat) <- 10
  for(k in 1:10){
    x_0 <- as.matrix(cbind(1,X[,-1])) 
    x_1 <- x_0
    if(!is.null(node_treatment_status)){
      x_0[,main_index + 1] <- node_treatment_status
      x_1[,main_index + 1] <- node_treatment_status
    }
    
    # NOTE - plus 1 is due to the adding of the intercept
    x_0[,(k_treat_index+1)] <- 0
    x_1[,(k_treat_index+1)] <- k
    
    k_peer_treat[[k]] <- sapply(sample(1:dim(samp)[1],1000),function(i){
      
      return(ATE_logistic_func(t(as.matrix(samp[i,])),
                               x_0 = x_0,
                               x_1 = x_1))
    })
  }
  
  k_peer_out <- list()
  length(k_peer_out) <- 10
  for(k in 1:10){
    x_0 <- as.matrix(cbind(1,X[,-1])) 
    x_1 <- x_0
    x_0[,(k_out_index+1)] <- 0
    x_1[,(k_out_index+1)] <- k
    
    k_peer_out[[k]] <- sapply(sample(1:dim(samp)[1],1000),function(i){
      
      return(ATE_logistic_func(t(as.matrix(samp[i,])),
                               x_0 = x_0,
                               x_1 = x_1))
    })
  }
  logistic_posteriors <- c(list(main),k_peer_out,k_peer_treat)
}
# Use logistic infer function after redoing fit based on ERGM simulations
extract_posteriors_ERGM_glm <- function(model,
                                        formula,
                                        nsims=100,
                                        ergm_theta_samples,
                                        glm_theta_samples =100,
                                        main_index,
                                        k_treat_index,
                                        k_out_index,
                                        neighs_to_add = NULL,
                                        treatment_indicator = 1,
                                        node_treatment_status = NULL,
                                        treatment_var = "treatment",
                                        outcome_var ="outcome",
                                        cores = NULL){
  #simulate from the posterior to get network that can be used to get the posterior distribution of the
  tmp_model <- ernm(formula,
                    maxIter = 3,
                    mcmcSampleSize = 100,
                    mcmcBurnIn = 100)
  
  # thin the sample
  samp <- model$ergm$theta_sample[sample(dim(model$ergm$theta_sample)[1],ergm_theta_samples,replace = T),]
  
  do_estim <- function(i,refit = T,return_sims = return_sims){
    # need the refit since can't pass CPP models to parallel clusters
    if(refit){
      tmp_model <- ernm(formula,
                        maxIter = 3,
                        mcmcSampleSize = 100,
                        mcmcBurnIn = 100)
    }
    
    theta <- samp[i,]
    tmp_model$m$setThetas(theta)
    sims <- tmp_model$m$sampler$generateSample(10000,1000,nsims)
    models <- lapply(sims,function(net){
      
      # change the network treatment so that we have the correct treatment indicator:
      if(treatment_indicator == 0){
        tmp <- net[[treatment_var]]
        tmp <- as.factor(1-as.numeric(convert_factor(tmp)))
        net[[treatment_var]] <- tmp
      }
      
      # now for each sim do the glm addition bit:
      assign(outcome_var,binary_convert(as.numeric(net$getVariable(outcome_var))))
      assign(treatment_var,binary_convert(as.numeric(net$getVariable(treatment_var))))
      
      net <- add_treated_neighs_binaryNet(net,treatment_var)
      net <- add_treated_neighs_binaryNet(net,outcome_var)
      
      assign(paste(treatment_var,"_neighbors",sep=""),as.numeric(net$getVariable(paste(treatment_var,"_neighbors",sep=""))))
      assign(paste(outcome_var,"_neighbors",sep=""),as.numeric(net$getVariable(paste(outcome_var,"_neighbors",sep=""))))
      #browser()
      if(!is.null(neighs_to_add)){
        for(x in neighs_to_add){
          net <- add_treated_neighs_binaryNet(net,x)
          assign(paste(x,"_neighbors",sep=""),as.numeric(net$getVariable(paste(x,"_neighbors",sep=""))))
        }
      }
      
      newdata <- list(eval(parse(text = outcome_var)),
                       eval(parse(text = treatment_var)))
      names = c(outcome_var,
                treatment_var)
      x <- c(outcome_var,treatment_var,neighs_to_add)
      for(i in 1:length(x)){
        newdata[[i+2]] = eval(parse(text = paste(x[[i]],"_neighbors",sep="")))
        names = c(names,paste(x[[i]],"_neighbors",sep=""))
      }
      names(newdata) <- names
      newdata <- as.data.frame(newdata)
      # make the right order:
      reg_terms <- labels(terms(model$glm$formula))
      reg_terms <- c(reg_terms,outcome_var)
      matches <- match(reg_terms,names(newdata))
      if(sum(is.na(matches)!=0)){
        print("Need to add some more treated neighbors- will error out!")
      }
      newdata <- newdata[,matches]
      
      #if the nw data has any constant colummns add a very small amount of random noise:
      sd <- apply(newdata,2,function(x){sd(as.numeric(x),na.rm = T)})
      if(any(sd == 0)){
        for(i in 1:length(sd)){
          if(sd[i]==0){
            print("adding random noise to variable")
            print(names(newdata)[i])
            newdata[,i] <- rnorm(dim(newdata)[1],0,min(sd)*0.1)
          }
        }
      }

      # new_model <- model$glm
      # new_model <- update(new_model,data=newdata)
      new_model  <- rstanarm::stan_glm(formula = model$glm$formula,family = model$glm$family,data = newdata)
      #print(new_model$coefficients)
      return(new_model)
    }
    )
    
    # return the estimates for each of these models:
    #browser()
    return(lapply(models,function(m){
      #browser()
      # print("Printing mean and sd of the fitted model coefficients")
      # print(apply(as.matrix(m),2,mean))
      # print(apply(as.matrix(m,),2,sd))
      extract_posteriors_logistic_directly(model=m,
                                           main_index = main_index,
                                           k_treat_index = k_treat_index,
                                           k_out_index = k_out_index,
                                           theta_samples = glm_theta_samples,
                                           node_treatment_status = node_treatment_status)
    }
    )
    )
  }
  if(!is.null(cores)){
    cl <- parallel::makeCluster(cores,type = "FORK")
    parallel::clusterSetRNGStream(cl, iseed = 1)
    parallel::clusterEvalQ(cl,{library(network)})
    parallel::clusterEvalQ(cl,{library(ernm)})
    parallel::clusterEvalQ(cl,{library(rstanarm)})
    parallel::clusterExport(cl,c("samp","binary_convert","nsims","convert_factor","add_treated_neighs_binaryNet"),envir = environment())
    estims <- parallel::parLapply(cl = cl, 1:dim(samp)[1], fun = do_estim)
    parallel::stopCluster(cl)
  }else{
    estims <- lapply(1:dim(samp)[1],do_estim)
  }

  estims <- do.call(c,estims)
    
  result <- list()
  result[[1]] <- do.call(c,lapply(estims,function(x){x[[1]]}))
  for(k in 2:11){
    result[[k]] <-  do.call(c,lapply(estims,function(x){x[[k]]}))
  }
  for(k in 12:21){
    result[[k]] <-  do.call(c,lapply(estims,function(x){x[[k]]}))
  }
  names(result) <- c("main",
                     sapply(seq(1,10),function(i){paste("k_peer_outcome",i,sep="")}),
                     sapply(seq(1,10),function(i){paste("k_peer_treatment",i,sep="")}))
  return(result)
}

# compares the KL divergence for each of the 
KL_distances <- function(model_list,lower_bound = -10,upper_bound = 10,bandwidth_adjust=1,points = 1000){
  
  # compare the K-L divergence for each estimand for each model in the list:
  #estimate the densities for each of the models
  dens <- list()
  for(i in 1:length(model_list)){
    dens[[i]] <- lapply(model_list[[i]],function(x){
      x = x[!is.na(x)]
      if(length(x) >2){return(density(x,n=points,from = lower_bound,to = upper_bound,adjust = bandwidth_adjust)$y)}else{
        return(0)}
    })
  }
  
  #now get the KL distance
  KL_dists <- list()
  length(KL_dists) <- length(model_list[[1]])
  dims <- length(model_list)
  
  
  for(i in 1:length(KL_dists)){
    
    KL_dists[[i]] <- outer(1:dims,1:dims,FUN = function(x,y){
      
      
      x = lapply(x,function(z){dens[[z]][[i]]})
      y = lapply(y,function(z){dens[[z]][[i]]})
      
      tmp <- mapply(x,y,FUN = function(p,t){
        if((length(p)<2) | (length(t)<2)){
          return(NaN)
        }
        keep <- which(as.logical((!is.na(p))*(!is.na(t))*(!(p==Inf))*(!(t==Inf))))
        t = t[keep]
        p = p[keep]
        p = p[t!=0]
        t = t[t!=0]
        # print("these are some densities to test")
        # print(p)
        # print(t)
        if((length(p)<2) | (length(t)<2)){
          return(NaN)
        }
        return(LaplacesDemon::KLD(p,t)$sum.KLD.px.py)
        
        #return(sum(p*log(p/t)))
      },SIMPLIFY =F)
      
      return(unlist(tmp))
    })
  }
  names(KL_dists) <- names(model_list[[1]])
  
  return(KL_dists)
}

# compares KL distances using the relative rank method
KL_distances_rel_rank <- function(model_list,method = "gam"){
  
  #now get the KL distance
  KL_dists <- list()
  length(KL_dists) <- length(model_list[[1]])
  dims <- length(model_list)
  
  
  for(i in 1:length(KL_dists)){
    
    KL_dists[[i]] <- outer(1:dims,1:dims,FUN = function(x,y){
      # estimate the relative distribution for models x, and y  for causal estimand i
      # accees the correct samples
      x = lapply(x,function(x){model_list[[x]][[i]]})
      y = lapply(y,function(y){model_list[[y]][[i]]})
      tmp <- mapply(x,y,FUN = function(x,y){
        if(sum(is.na(y)) == length(y) | sum(is.na(x)) == length(x)){
          return(NA)
        }
        
        tmp <- reldist::reldist(y=y,yo=x,method = method,graph = F)
        KL <- mean(tmp$y*log(tmp$y))
        return(KL)
      },SIMPLIFY =F)
      tmp <- unlist(tmp)
      return(tmp)
    })
  }
  return(KL_dists)
}

# Function to get the various GOF measures from the sims

get_gof <- function(sims){
  
  degree <- lapply(sims,function(net){
    tmp <- ergm::summary_formula(net ~ degree(0:20))
    return(tmp)
  })
  
  degree <- do.call(rbind,degree)
  
  
  smoker_degree <- lapply(sims,function(net){
    tmp <- ergm::summary_formula(net ~ degree(0:20,by = "c.smoke"))[22:42]
    return(tmp)
  })
  
  smoker_degree <- do.call(rbind,smoker_degree)
  
  triads <- sapply(sims,function(net){
    tmp <- ergm::summary_formula(net ~ triangles)
    return(tmp)
  })
  
  smoker_triads  <- sapply(sims,function(net){
    tmp <- ergm::summary_formula(net ~ triangles("c.smoke"))
    return(tmp)
  })
  
  esp <- lapply(sims,function(net){
    tmp <- ergm::summary_formula(net ~ esp(0:20))
    return(tmp)
  })
  
  esp <- do.call(rbind,esp)
  
  geodist <- lapply(sims,function(net){
    ergm.geodistdist(net)[1:15]
  })
  
  geodist <- do.call(rbind,geodist)
  
  smoker_edges <- sapply(sims,function(net){
    tmp <- ergm::summary_formula(net ~ nodematch("c.smoke"))
    return(tmp)
  })
  
  edges <-  sapply(sims,function(net){
    tmp <- ergm::summary_formula(net ~ edges)
    return(tmp)
  })
  
  return(list(degree = degree,
              esp = esp,
              geodist = geodist,
              smoker_degree = smoker_degree,
              smoker_edges = smoker_edges,
              triads = triads,
              edges = edges,
              smoker_triads = smoker_triads))
}

mean_compare_plot <- function(means,sd,title){
  means <- as.data.frame(means)
  #build the long mean + sd dataset:
  means$method <- rownames(means)
  mean_long <- reshape2::melt(means,id.vars = "method")
  names(mean_long)[3] <- "mean"
  
  sd <- as.data.frame(sd)
  sd$method <- rownames(sd)
  sd_long <- reshape2::melt(sd,id.vars = "method")
  names(sd_long)[3] <- "sd"
  
  tmp <- mean_long %>%dplyr::full_join(sd_long,by = c("variable","method"))
  #remove main effects
  tmp <- tmp[tmp$variable != "main",]
  tmp$variable <- levels(tmp$variable)[tmp$variable]
  
  tmp$variable <- sapply(tmp$variable,function(x){as.numeric(substring(x,nchar(x),nchar(x)))}) 
  names(tmp)[2] <- "peers"
  tmp$peers[tmp$peers==0] <- 10
  
  plot <- ggplot(data = tmp,aes(x =peers,y = mean,col = method))+
    geom_line(position=position_dodge(0.1))+
    geom_errorbar(position=position_dodge(0.1),aes(ymin=mean-sd, ymax=mean+sd),
                  width=.2)+
    scale_x_continuous(breaks = seq(0,10))+
    ylab("k-peer-effect")+
    ggtitle(title)
  plot
}

# function to extract smoker status and proportion of gender homogenous friends from sims:
extract_prop_treated_friends <- function(sims,
                                 treatment_var,
                                 outcome_var,
                                 treatment_indicator){
  outs = lapply(sims,function(net){
    net <- add_treated_neighs_ind(net,treatment_var,treatment_indicator)
    # get the outcomes
    treats = (as.numeric(get.vertex.attribute(net,treatment_var)) == treatment_indicator)*1
    outs = as.numeric(get.vertex.attribute(net,outcome_var))
    neighs = sapply(1:(net%n%'n'),function(i){length(get.neighborhood(net,i))})
    prop_treated_neighs = as.numeric((net %v% paste(treatment_var,"_neighbors",sep="")))/neighs
    return(data.frame(
      treats = treats,
      outs = outs,
      neighs = neighs,
      prop_treated_neighs = prop_treated_neighs
    ))
  })
  return(do.call(rbind,outs))
}

lolog_boxplot <- function(sims_list,
                          names_list,
                          observed_list,
                          title = "",
                          title_list = "",
                          xlab = "",
                          ylim = NULL,
                          superimpose = FALSE, #if false will be facetted, else will be superimposed
                          superimpose_sims_list =NULL, #if null will superimpose all on one grid, otherwise will facet and superimpose this on all facets
                          superimpose_names_list = NULL,
                          superimpose_observed_list = NULL,
                          superimpose_title_list = NULL,
                          super_labels = c("Yes","No")
){
  
  list_process_func <- function(sims,names,observed,title){
    sims <- as.data.frame(sims)
    names(sims) <- names
    #make sims only as big as the names vector:
    sims <- sims[,1:length(names)]
    
    sims <- reshape2::melt(sims)
    names(sims) <- c("statistic","observed")
    sims$statistic <- as.factor(sims$statistic)
    
    observed <- as.data.frame(matrix(observed[1:length(names)],ncol = 1))
    names(observed) <- "observed"
    observed$statistic <- levels(sims[,1])
    observed$type <- "observed"
    return(list(sims =sims, observed = observed))
  }
  
  data <- mapply(FUN =list_process_func,sims_list,names_list, observed_list,title_list,SIMPLIFY =  FALSE)
  data <- mapply(data,title_list,FUN = function(x,y){
    x$sims$facet =y
    x$observed$facet = y
    return(x)},SIMPLIFY = FALSE)
  sims <- do.call(rbind,lapply(data,function(x){x$sims}))
  obs <- do.call(rbind,lapply(data,function(x){x$observed}))
  
  if(superimpose){
    if(is.null(superimpose_sims_list) | is.null(superimpose_observed_list) | is.null(superimpose_names_list)){
      names(sims)[which(names(sims)=="facet")] <- 'model'
      names(obs)[which(names(obs)=="facet")] <- 'model'
      
      gg <- ggplot(data = sims,aes(x = statistic, y = observed, group = interaction(statistic,model),fill = model)) +
        geom_boxplot(position = 'identity',alpha =0.5) +
        ggduncan::scale_fill_duncan() +
        geom_point(data = obs,aes(y=observed,x=statistic,group = interaction(statistic,model),col = type),cex = 3)+
        scale_colour_manual(labels = c("Observed Values"), values =  "Magenta",guide = FALSE)+
        xlab(xlab)+
        ylab("")+
        ggduncan::theme_duncan(axis.text.x = ggplot2::element_text(angle = 45,
                                                                   hjust = 1),
                               panel.grid.major.x = ggplot2::element_blank(),
                               panel.grid.minor.x = ggplot2::element_blank())+
        ggtitle(title)
      
    }else{
      #get the data for the layer to be superimposed on all:
      data_super <- mapply(FUN =list_process_func,superimpose_sims_list,superimpose_names_list, superimpose_observed_list,title_list,SIMPLIFY =  FALSE)
      data_super <- mapply(data_super,superimpose_title_list,FUN = function(x,y){
        x$sims$facet =y
        x$observed$facet = y
        return(x)},SIMPLIFY = FALSE)
      sims_super <- lapply(data_super,function(x){x$sims})
      obs_super <- lapply(data_super,function(x){x$observed})
      
      sims$super <- super_labels[2]
      sims_super <- lapply(sims_super,function(x){x$super <- super_labels[1]
      return(x)})
      sims <- rbind(sims,do.call(rbind,sims_super))
      names(sims)[which(names(sims)=="facet")] <- 'model'
      names(obs)[which(names(obs)=="facet")] <- 'model'
      
      gg <- ggplot(data = sims,aes(x = statistic, y = observed, group = interaction(statistic,super)))+
        facet_grid(rows = vars(model),scales = "free_y")+
        geom_boxplot(position = 'identity',alpha =0.5,aes(fill = super))+
        ggduncan::scale_fill_duncan()+
        geom_point(data = obs,aes(y=observed,x=statistic,group = interaction(statistic,model),col = type),cex = 3)+
        scale_colour_manual(labels = c("Observed Values"), values =  "Magenta",guide = FALSE)+
        xlab(xlab)+
        ylab("")+
        ggplot2::theme_bw() +
        ggduncan::theme_duncan(axis.text.x = ggplot2::element_text(angle = 45,
                                                                   hjust = 1),
                               panel.grid.major.x = ggplot2::element_blank(),
                               panel.grid.minor.x = ggplot2::element_blank())+
        ggtitle(title)
    }
    
  }else{
    
    gg <- ggplot(data = sims,aes(x = statistic, y = observed, group = statistic))+
      geom_boxplot()+
      ggduncan::scale_fill_duncan()+
      facet_grid(rows = vars(facet),scales = "free_y")+
      geom_point(data = obs,aes(y=observed,x=statistic,group = statistic,col = type),cex =3 )+
      scale_colour_manual(labels = c("Observed Values"), values =  "Magenta",guide = FALSE)+
      xlab(xlab)+
      ylab("")+
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                         hjust = 1), panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank())+
      ggduncan::theme_duncan()
    ggtitle(title)
    
  }
  
  if(is.null(ylim)){}else{
    gg <- gg + coord_cartesian(ylim = c(0,ylim))
  }
  return(gg)
}

#without facetting
lolog_boxplot_single <- function(sims,
                                 names,
                                 observed,
                                 title = "",
                                 xlab = ""){
  
  sims <- as.data.frame(sims)
  names(sims) <- names
  sims <- reshape2::melt(sims)
  names(sims) <- c("statistic","observed")
  sims$statistic <- as.factor(sims$statistic)
  
  observed <- as.data.frame(matrix(observed,ncol = 1))
  names(observed) <- "observed"
  observed$statistic <- levels(sims[,1])
  observed$type <- "observed"
  
  gg <- ggplot(data = sims,aes(x = statistic, y = observed, group = statistic))+
    geom_boxplot()+
    geom_point(data = observed,aes(y=observed,x=statistic,group = statistic,col = type))+
    scale_colour_manual(labels = c("Observed Values"), values =  "red",guide = FALSE)+
    xlab(xlab)+
    ylab("")+
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                                                       hjust = 1), panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank())
  ggtitle(title)
  return(gg)
}






