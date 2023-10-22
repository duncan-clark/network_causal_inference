
library(statnet)
library(ernm)
library(mvtnorm)
library(coda)
library(parallel)
library(doParallel)
library(rstanarm)
library(lolog.catelog.helper)
library(ggplot2)
library(scales)
source("functions/utils.R")
source("functions/fire.R")

# Configs:
SAVE_DIR <- ""

# ======================================================================================================================
# Setup
# ======================================================================================================================
load("data/P0AH_s270_weak_.RData")
str(pops)
nw <- network(pops,directed=FALSE)

n <- N <- network.size(nw)
vars <- data.frame(attr(pops,"vertex.attributes"))
d <- sna::degree(nw)/2
el <- as.data.frame(as.matrix(nw,matrix.type="edgelist"))
el[,3] <- 1
el <- as.matrix(el)
g <- vars$grade > 10
#
vn <- nw %v% "vertex.names"
vnm <- match(vn, vars[,"vertex.names"])
for(i in 1:ncol(vars)){
  set.vertex.attribute(nw, colnames(vars)[i], vars[vnm,i])
}

# make a categorical grade, sex and smoke variable, make into 0s and 1s
set.vertex.attribute(nw,"c.sex",as.character(get.vertex.attribute(nw,'sex')-1))
set.vertex.attribute(nw,"c.sex_female",as.character(get.vertex.attribute(nw,'sex')-1))
set.vertex.attribute(nw,"c.sex_male",as.character(1 - (get.vertex.attribute(nw,'sex')-1)))
set.vertex.attribute(nw,"c.grade",as.character(get.vertex.attribute(nw,'grade')))
set.vertex.attribute(nw,"c.smoke",as.character(get.vertex.attribute(nw,'smoke')))
set.vertex.attribute(nw,"c.intercept",as.character(1))

# functions to deal with treated neighbors
add_treated_neighs <- function(net,treatment_var){
  tmp <- as.numeric(get.vertex.attribute(net,treatment_var))
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
nw <- add_treated_neighs(nw,"c.smoke")
nw <- add_treated_neighs(nw,"c.sex")
nw <- add_treated_neighs(nw,"c.sex_female")
nw <- add_treated_neighs(nw,"c.sex_male")


regression_formula <- "c.sex + c.sex_female_neighbors + c.sex_male_neighbors + c.smoke_neighbors"
regression_formula_2 <- "c.sex + c.sex_female_neighbors + c.sex_male_neighbors"
regression_formula_3 <- "c.sex + c.sex_female_neighbors + c.smoke_neighbors"
ERGM_formula <- "edges() + 
                gwesp(0.5,homogenous = T,variableName = 'c.grade') + 
                gwdegree(0.5) + 
                homophily('c.grade') +
                homophily('c.sex') +
                homophily('c.smoke')"
ERNM_formula <- "edges() + 
                gwesp(0.5,homogenous = T,variableName = 'c.grade') + 
                gwdegree(0.5) + 
                homophily('c.grade') +
                homophily('c.sex') +
                homophily('c.smoke') +
                nodeCount('c.smoke') +
                logisticNeighbors('c.smoke','c.sex') +
                logisticNeighborsTopLevel('c.smoke','c.smoke') +
                logistic('c.smoke','c.sex') | c.smoke"
MRF_formula <- "homophily('c.smoke') +
                nodeCount('c.smoke') +
                logisticNeighbors('c.smoke','c.sex') +
                logisticNeighborsTopLevel('c.smoke','c.smoke') +
                logistic('c.smoke','c.sex') | c.smoke"

# Note that differential homophily for the smoker terms leads to problems

# Use uniform prior that drops off at -10 and 10
unif_prior <- function(x){
  if(any(abs(x) > 10) ){
    return(0)
  }
  else{return(1)}
}


#====================================================================================================================
# Pilot Runs 
#====================================================================================================================

# Use this theta to get a reasonable starting point - avoid starting in degenerate region
ernm_mle <- ernm(as.formula(paste("nw ~ " ,ERNM_formula)),maxIter = 5)
ernm_theta0 <- ernm_mle$theta

ernm_theta0 <- ernm(as.formula(paste("nw ~ " ,ERNM_formula)),fitModel = F)$theta()
t<-proc.time()[3]
ernm_pilot <- fire(as.formula(paste("nw ~ " ,ERNM_formula)),
                   prior_dens = unif_prior,
                   theta0 = ernm_theta0,
                   ernm_burnin = 10000,
                   ernm_steps = 1000,
                   burn_in = 100,
                   samples = 10,
                   n_sims = 200,
                   gamma = 0.1,
                   run_type = "pilot",
                   verbose = T,
                   cores = 9,
                   pilot_chains = 8)
print(proc.time()[3] - t)
ernm_pilot$time
ernm_pilot$median_acceptance
for(i in 1:8){
  print(summary(unlist(ernm_pilot$chains[[i]]$probs)))
}
# ernm_pilot = NULL
# 100 steps takes ~ 12 minutes


#1000 may not be sufficient for this one!
# - with 200 samples gamma = 0.1 we still get some "falling off"
# - try with 500 samples gamma  = 0.05
mrf_theta0 <- ernm(as.formula(paste("nw ~ " ,MRF_formula)),maxIter = 20,nodeSamplingPercentage = 1)$theta
mrf_pilot <- fire(as.formula(paste("nw ~ " ,MRF_formula)),
                   prior_dens = unif_prior,
                   theta0 = mrf_theta0,
                   ernm_burnin = 10000,
                   ernm_steps = 1000,
                   burn_in = 100,
                   samples = 10,
                   n_sims = 20,
                   gamma = 0.2,
                   run_type = "pilot",
                   verbose = T,
                   cores = 9,
                   pilot_chains = 8,
                   fixed_edges  = T)

mrf_pilot$time
mrf_pilot$median_acceptance
for(i in 1:8){
  print(summary(unlist(mrf_pilot$chains[[i]]$probs)))
}

ergm_theta0 <- ernm(as.formula(paste("nw ~ " ,ERGM_formula)),maxIter = 20,nodeSamplingPercentage = 0)$theta
ergm_pilot <- fire(as.formula(paste("nw ~ " ,ERGM_formula)),
                  prior_dens = unif_prior,
                  theta0 = ergm_theta0,
                  ernm_burnin = 10000,
                  ernm_steps = 1000,
                  burn_in = 100,
                  samples = 10,
                  n_sims = 50,
                  gamma = 0.2,
                  run_type = "pilot",
                  verbose = T,
                  cores = 9,
                  pilot_chains = 8,
                  fixed_edges  = F,
                  fixed_nodes = T)

ergm_pilot$time
ergm_pilot$median_acceptance
for(i in 1:8){
  print(summary(unlist(ergm_pilot$chains[[i]]$probs)))
}

#====================================================================================================================
# Model Fits 
#====================================================================================================================

# ERNM model
#ernm_mle <- ernm(as.formula(paste("nw ~ " ,ERNM_formula)),maxIter = 50)
#ernm_theta0 <- ernm_mle$theta
t<-proc.time()[3]
ernm_model <- fire(as.formula(paste("nw ~ " ,ERNM_formula)),
                   prior_dens = unif_prior,
                   theta0 = ernm_theta0,
                   ernm_burnin = 10000,
                   ernm_steps = 1000,
                   burn_in = 1000,
                   samples = 1000,
                   n_sims = 200,
                   gamma = 0.1,
                   verbose = F,
                   run_type = "main",
                   cores = NULL,
                   pilot_chains = NULL)
print("ERNM model took:")
print(proc.time()[3] - t)
# ernm model takes ~ 4 hours to  run
apply(ernm_model$theta_sample,2,mean)
apply(ernm_model$theta_sample,2,sd)

# MRF model:

mrf_model <- fire(as.formula(paste("nw ~ " ,MRF_formula)),
                   prior_dens = unif_prior,
                   theta0 = mrf_theta0,
                   ernm_burnin = 10000,
                   ernm_steps = 1000,
                   burn_in = 2000,
                   samples = 1000,
                   n_sims = 20,
                   gamma = 0.2,
                   verbose = F,
                   run_type = "main",
                   cores = NULL,
                   pilot_chains = NULL,
                   fixed_edges = T)

apply(mrf_model$theta_sample,2,mean)
apply(mrf_model$theta_sample,2,sd)


# MRF pseudo likelihood model

# make design matrix for MRF model
mrf_tmp_model <- ernm(as.formula(paste("nw ~ " ,MRF_formula)),
                      maxIter = 2,
                      nodeSamplingPercentage = 1)
mrf_tmp_model <- mrf_tmp_model$m$sampler$getModel()
# set the network after the fitting:
mrf_tmp_model$setNetwork(as.BinaryNet(nw))
mrf_tmp_model$calculate()
# note that change to 1 corresponds to change from "1" to "0"
# note that change to 2 corresponds to change from "0" to "1" 
# for some reason!

mrf_pseudo_data <- lapply(1:n,function(i){
  # calculate the change stat conditional on the rest of the network
  stats <- mrf_tmp_model$statistics()
  old <- (nw %v% "c.smoke")[i]
  new <- (old=="1")*1 + (old=="0")*2
  mrf_tmp_model$discreteVertexUpdate(i,'c.smoke',new)
  new_stats <- mrf_tmp_model$statistics()
  # change stat should be toggle on so need to adjust if smoker or not
  change_stats =  (new_stats- stats) * ((-1)**(old == "1"))
  #reset
  mrf_tmp_model$calculate()
  return(change_stats)
})
mrf_pseudo_data <- as.data.frame(do.call(rbind,mrf_pseudo_data))


mrf_pseudo_data$nodecount.c.smoke.1 = mrf_pseudo_data$nodecount.c.smoke.1*-1
names(mrf_pseudo_data) = c("homophily('c.smoke')",
                           "nodeCount('c.smoke')",
                           "logisticNeighbors('c.smoke','c.sex')",
                           "logisticNeighborsTopLevel('c.smoke','c.smoke')",
                           "logistic('c.smoke','c.sex')")
mrf_pseudo_data$'c.smoke' = as.numeric((nw %v% "c.smoke"))

mrf_pseudo_model <- glm(c.smoke ~ .,
                        data = mrf_pseudo_data[,-2],
                        family=binomial(link="logit"))

mrf_pseudo_model$coefficients
apply(mrf_model$theta_sample,2,mean)[c(2,1,3,4,5)]

mrf_pseudo_model <- list(theta_sample = matrix(mrf_pseudo_model$coefficients[c(2,1,3,4,5)],nrow=1))
mrf_pseudo_model[[1]][1,] <- mrf_pseudo_model[[1]][1,]*1 # do this since the MRF model formula is the other way around on node count


c.sex <- as.numeric(nw %v% "c.sex")
c.sex_female_neighbors <- as.numeric(nw %v% "c.sex_female_neighbors")
c.sex_male_neighbors <- as.numeric(nw %v% "c.sex_male_neighbors")
c.smoke_neighbors <- as.numeric(nw %v% "c.smoke_neighbors")
c.smoke <- as.numeric(nw %v% "c.smoke")
glm_data = data.frame(c.sex = c.sex,
                      c.sex_female_neighbors = c.sex_female_neighbors,
                      c.sex_male_neighbors = c.sex_male_neighbors,
                      c.smoke_neighbors = c.smoke_neighbors,
                      c.smoke = c.smoke)

logistic_model <- rstanarm::stan_glm(as.formula(paste("c.smoke ~",regression_formula)),
                                     data = glm_data,
                                     family=binomial(link="logit"))
logistic_model$coefficients

logistic_model <- rstanarm::stan_glm(as.formula(paste("c.smoke ~",regression_formula_2)),
                                     data = glm_data,
                                     family=binomial(link="logit"))
logistic_model$coefficients
diag(logistic_model$covmat)**0.5
cov2cor(logistic_model$covmat)

logistic_model_3 <- rstanarm::stan_glm(as.formula(paste("c.smoke ~",regression_formula_3)),
                                     data = glm_data,
                                     family=binomial(link="logit"))
logistic_model_3$coefficients
diag(logistic_model_3$covmat)**0.5
cov2cor(logistic_model_3$covmat)

# Comment: negative coefficient is result of smoker neigbors being confounded with female neighbors!


ergm_theta0 <- ernm(as.formula(paste("nw ~ " ,ERGM_formula)),maxIter = 20,nodeSamplingPercentage = 0)$theta
ergm_model <- fire(as.formula(paste("nw ~ " ,ERGM_formula)),
                   prior_dens = unif_prior,
                   theta0 = ergm_theta0,
                   ernm_burnin = 10000,
                   ernm_steps = 1000,
                   burn_in = 2000,
                   samples = 1000,
                   n_sims = 50,
                   gamma = 0.2,
                   run_type = "main",
                   verbose = F,
                   cores = NULL,
                   pilot_chains = NULL,
                   fixed_nodes = T)
ergm_model <- list(ergm=ergm_model)

ergm_model <- list(ergm = ergm_model$ergm,glm = logistic_model)

# =========================================================================================================================
# Estimate Causal Estimands
# =========================================================================================================================
# make a function to get the gender facetted data for each of the 4 models:
get_effects <- function(
  ERNM_models,
  logistic_models,
  ERGM_models,
  ERNM_formulas,
  ERNM_fixed,
  ERGM_formulas,
  names,
  ERNM_nsims,
  ERGM_nsims,
  ERNM_draws,
  ERGM_draws,
  logistic_draws,
  cores,
  treatment_var,
  outcome_var,
  neighs_to_add
){
  # get the results for the ERNM models
  ERNM_results <- mapply(SIMPLIFY = FALSE,ERNM_models,ERNM_formulas,ERNM_fixed,FUN = function(model,form,fixed){
    print("male_female")
    male_female <- extract_posteriors_fast(model = model,
                                           formula = as.formula(paste("nw ~ " ,form)),
                                           theta_draws = ERNM_draws,
                                           fixed_edges = fixed,
                                           nsims = ERNM_nsims,
                                           treatment_var = treatment_var,
                                           treatment_indicator = 1,
                                           outcome_var = outcome_var,
                                           node_treatment_statuses = c(0),
                                           cores = cores)
    print("female_female")
    female_female <- extract_posteriors_fast(model = model,
                                             formula = as.formula(paste("nw ~ " ,form)),
                                             theta_draws = ERNM_draws,
                                             fixed_edges = fixed,
                                             nsims = ERNM_nsims,
                                             treatment_var = treatment_var,
                                             treatment_indicator = 1,
                                             outcome_var = outcome_var,
                                             node_treatment_statuses = c(1),
                                             cores = cores)
    print("male_male")
    male_male <- extract_posteriors_fast(model = model,
                                         formula = as.formula(paste("nw ~ " ,form)),
                                         theta_draws = ERNM_draws,
                                         fixed_edges = fixed,
                                         nsims = ERNM_nsims,
                                         treatment_var = treatment_var,
                                         treatment_indicator = 0,
                                         node_treatment_statuses = c(0),
                                         outcome_var = outcome_var,
                                         cores = cores)
    print("female_male")
    female_male <- extract_posteriors_fast(model = model,
                                           formula = as.formula(paste("nw ~ " ,form)),
                                           theta_draws = ERNM_draws,
                                           fixed_edges = fixed,
                                           nsims = ERNM_nsims,
                                           treatment_var = treatment_var,
                                           treatment_indicator = 0,
                                           outcome_var = outcome_var,
                                           node_treatment_statuses = c(1),
                                           cores = cores)
    return(list(male_female = male_female,
                female_female = female_female,
                male_male = male_male,
                female_male = female_male
                ))
    })
  print("ERNM done")
  
  # get the results for any ERGM models
  ERGM_results <- mapply(SIMPLIFY = FALSE,ERGM_models,ERGM_formulas,FUN = function(model,form){
    print("male_female")
    male_female <- extract_posteriors_ERGM_glm(model = model,
                                               formula = as.formula(paste("nw ~ " ,form)),
                                               nsims = ERGM_nsims,
                                               ergm_theta_samples = ERGM_draws,
                                               glm_theta_samples = logistic_draws,
                                               main_index = 1,
                                               k_treat_index=  2,
                                               k_out_index =  4,
                                               neighs_to_add = neighs_to_add,
                                               treatment_var = treatment_var,
                                               outcome_var = outcome_var,
                                               cores = cores,
                                               treatment_indicator = 1,
                                               node_treatment_status = 0)
    print("female_female")
    female_female <- extract_posteriors_ERGM_glm(model = model,
                                               formula = as.formula(paste("nw ~ " ,form)),
                                               nsims = ERGM_nsims,
                                               ergm_theta_samples = ERGM_draws,
                                               glm_theta_samples = logistic_draws,
                                               main_index = 1,
                                               k_treat_index=  2,
                                               k_out_index =  4,
                                               neighs_to_add = neighs_to_add,
                                               treatment_var = treatment_var,
                                               outcome_var = outcome_var,
                                               cores = cores,
                                               treatment_indicator = 1,
                                               node_treatment_status = 1)
    print("male_male")
    male_male <- extract_posteriors_ERGM_glm(model = model,
                                               formula = as.formula(paste("nw ~ " ,form)),
                                               nsims = ERGM_nsims,
                                               ergm_theta_samples = ERGM_draws,
                                               glm_theta_samples = logistic_draws,
                                               main_index = 1,
                                               k_treat_index=  3,
                                               k_out_index =  4,
                                               neighs_to_add = neighs_to_add,
                                               treatment_var = treatment_var,
                                               outcome_var = outcome_var,
                                               cores = cores,
                                               treatment_indicator = 1,
                                               node_treatment_status = 1)
    
    print("female_male")
    female_male <- extract_posteriors_ERGM_glm(model = model,
                                                 formula = as.formula(paste("nw ~ " ,form)),
                                                 nsims = ERGM_nsims,
                                                 ergm_theta_samples = ERGM_draws,
                                                 glm_theta_samples = logistic_draws,
                                                 main_index = 1,
                                                 k_treat_index=  3,
                                                 k_out_index =  4,
                                                 neighs_to_add = neighs_to_add,
                                                 treatment_var = treatment_var,
                                                 outcome_var = outcome_var,
                                                 cores = cores,
                                                 treatment_indicator = 1,
                                                 node_treatment_status = 0)
    return(list(male_female = male_female,
                female_female = female_female,
                male_male = male_male,
                female_male = female_male))
  })
  print("ERGM done")
  
  # get results for any logistic models:
  logistic_results <- lapply(logistic_models,FUN = function(model){
    
    male_female <- extract_posteriors_logistic_directly(model,
                                                        theta_samples = logistic_draws,
                                                        main_index = 1,
                                                        k_treat_index = 2,
                                                        k_out_index = 4,
                                                        node_treatment_status = 0)
    female_female <- extract_posteriors_logistic_directly(model,
                                                        theta_samples = logistic_draws,
                                                        main_index = 1,
                                                        k_treat_index = 2,
                                                        k_out_index = 4,
                                                        node_treatment_status = 1)
    
    male_male <- extract_posteriors_logistic_directly(model,
                                                        theta_samples = logistic_draws,
                                                        main_index = 1,
                                                        k_treat_index = 3,
                                                        k_out_index = 4,
                                                        node_treatment_status = 1)
    female_male <- extract_posteriors_logistic_directly(model,
                                                          theta_samples = logistic_draws,
                                                          main_index = 1,
                                                          k_treat_index = 3,
                                                          k_out_index = 4,
                                                          node_treatment_status = 0)
    return(list(male_female = male_female,
                female_female = female_female,
                male_male = male_male,
                female_male = female_male
    ))
  })
  print("logistic done")
  
  # now combine to make a mega data set!
  # k - peer outcomes
  #ERNM_names = names[1:length(ERNM_results)]
  # the effects that are currently calculated:
  effect_names = c("main",
                   "k_peer_outcome1","k_peer_outcome2","k_peer_outcome3","k_peer_outcome4","k_peer_outcome5",   
                   "k_peer_outcome6","k_peer_outcome7","k_peer_outcome8","k_peer_outcome9","k_peer_outcome10",
                   "k_peer_treatment1","k_peer_treatment2","k_peer_treatment3","k_peer_treatment4","k_peer_treatment5",
                   "k_peer_treatment6","k_peer_treatment7","k_peer_treatment8","k_peer_treatment9","k_peer_treatment10",
                   "neigh_treat","node_treat","method","var")  
  tmp <- mapply(SIMPLIFY = F,c(ERNM_results,ERGM_results,logistic_results),names,FUN = function(x,y){
    # combine all the different results
    tmp <- data.frame()
    for(i in 1:length(x)){
      print(i)
      # get the neighbor and node treatment
      name = names(x)[i]
      where = regexpr("_",name)[1]
      neigh_treat = substring(name,where+1,nchar(name)) 
      node_treat = substring(name,1,where-1) 
      
      # extract the means:
      means <- as.data.frame(lapply(x[[i]],mean,na.rm = T))
      means$neigh_treat = neigh_treat
      means$node_treat = node_treat
      means$method = y
      means$var = "mean"
      names(means) = effect_names
      tmp = rbind(tmp,means)
      
      # extract the sds. 
      sd <- as.data.frame(lapply(x[[i]],sd,na.rm = T))
      sd$neigh_treat = neigh_treat
      sd$node_treat = node_treat
      sd$method = y
      sd$var = "sd"
      names(sd) = effect_names
      tmp = rbind(tmp,sd)
    }
    return(tmp)})
  
  tmp <- reshape2::melt(tmp,id.vars = list("neigh_treat","method","var","node_treat"))
  
  return(list(raw_results = list(
                            ERNM = ERNM_results,
                            ERGM = ERGM_results,
                            logistic = logistic_results),
         processed_results = tmp))
}

t<-proc.time()
all_effects <- get_effects(logistic_models = list(logistic_model),
                           ERGM_models = list(ergm_model),
                           ERNM_models = list(ernm_model,mrf_model,mrf_pseudo_model),
                           ERNM_formulas = list(ERNM_formula,MRF_formula,MRF_formula),
                           ERNM_fixed = list(F,T,T),
                           ERGM_formulas = list(ERGM_formula),
                           names = list("ERNM","MRF","pseudo_MRF","ERGM","logistic"),
                           ERGM_nsims = 50,
                           ERNM_nsims = 10,
                           ERNM_draws = 100,
                           ERGM_draws = 20,
                           logistic_draws = 1000,
                           cores = 20,
                           treatment_var = "c.sex",
                           outcome_var = "c.smoke",
                           neighs_to_add = c("c.sex_male","c.sex_female")
)
print(proc.time() -t)


# make the mega facet plot:
# parse the effects
tmp <- all_effects$processed_results
tmp$peers <- sapply(convert_factor(tmp$variable),function(x){
  #check if last two digits are number 
  peer = substring(x,nchar(x)-1,nchar(x))
  if(!is.na(as.numeric(peer))){return(as.numeric(peer))}
  # otherwise use the last one
  peer = substring(x,nchar(x),nchar(x))
  if(!is.na(as.numeric(peer))){return(as.numeric(peer))}
  # otherwise return 0
  return(0)
})
tmp$type <- mapply(convert_factor(tmp$variable),tmp$peers,FUN = function(x,y){
  # remove the peers from the end
  x = substring(x,1,nchar(x)-nchar(y))
  locate = gregexpr("_",x)[[1]][2]
  if(is.na(locate)){
    #it was the main effect
    return("treatment")
  }
  type = substring(x,locate+1,nchar(x))
  return(type)
})
tmp <- reshape2::dcast(tmp,... ~ var)
tmp$node_treat <- sapply(tmp$node_treat,function(x){paste(x,"_nodes")})
tmp$neigh_treat <- sapply(tmp$neigh_treat,function(x){paste(x,"_neighs")})


plot <- ggplot(data = tmp[tmp$type=="treatment",],aes(x =peers,y = mean,col = method))+
  facet_grid(neigh_treat ~ node_treat) +
  geom_line(position=position_dodge(0.1))+
  geom_errorbar(position=position_dodge(0.1),aes(ymin=mean-sd, ymax=mean+sd),
                width=.2)+
  scale_x_continuous(breaks = seq(0,10))+
  ylab("k-peer-effect")+
  ggtitle("Mega plot of treatments")
plot

# Logistic Model:

# This posterior distribution gives us the effects as required
logistic_posteriors <- extract_posteriors_logistic_directly(logistic_model,
                                                            theta_samples = 1000,
                                                            main_index = 1,
                                                            k_treat_index = 2,
                                                            k_out_index = 4,
                                                            node_treatment_status = NULL)

# MRF Model

# MRF_posteriors <- extract_posteriors_fast(mrf_model,theta_draws=100,fixed_edges = T,nsims = 100,treatment_var = "c.sex",outcome_var ="c.smoker")

# ERGM + Logistic Model:

# low level of ergm theta sample, and network sims since these slow it down a lot ...
t <- proc.time()
ERGM_posteriors <- extract_posteriors_ERGM_glm(model = ergm_model,
                                               formula = as.formula(paste("nw ~ " ,ERGM_formula)),
                                               nsims = 10,
                                               ergm_theta_samples = 20,
                                               glm_theta_samples = 1000,
                                               main_index=1,
                                               k_treat_index=2,
                                               k_out_index=4,
                                               neighs_to_add =  c("c.sex_male","c.sex_female"),
                                               treatment_var = "c.sex",
                                               outcome_var ="c.smoke",
                                               cores = 20,
                                               node_treatment_status = NULL)
print(proc.time() -t)
# Thins since is too big !
ERGM_posteriors <- lapply(ERGM_posteriors,function(x){sample(x,100000)})

# COMMENT:
# - The fully separable model below is naive and does not accout  for network uncertainty
# - Used the fit on the observed network on network that could have been observed 
# - Stay with the above multi fitting rstanarm PO inferer function since is in line with Toulis2013 approach.

t <- proc.time()
ERGM_posteriors <- extract_posteriors_ERGM_fast(ergm_model,
                                                formula = as.formula(paste("nw ~ " ,ERGM_formula)),
                                                theta_draws = 100,
                                                nsims = 500,
                                                treatment_var = "c.sex",
                                                outcome_var ="c.smoke",
                                                cores = 20)
print(proc.time() -t)

t <- proc.time()
MRF_posteriors <- extract_posteriors_fast(mrf_model,
                                          formula = as.formula(paste("nw ~ " ,MRF_formula)),
                                          theta_draws = 100,
                                          fixed_edges = T,
                                          nsims = 50,
                                          treatment_var = "c.sex",
                                          outcome_var ="c.smoke",
                                          node_treatment_statuses = c(0,1),
                                          cores = 20)
print(proc.time() -t)


t <- proc.time()
MRF_pseudo_posteriors <- extract_posteriors_fast(mrf_pseudo_model,
                                          formula = as.formula(paste("nw ~ " ,MRF_formula)),
                                          theta_draws = 100,
                                          fixed_edges = T,
                                          nsims = 50,
                                          treatment_var = "c.sex",
                                          outcome_var ="c.smoke",
                                          node_treatment_statuses = c(0,1),
                                          cores = 20)
print(proc.time() -t)


# ERNM model
t <- proc.time()
ERNM_posteriors <- extract_posteriors_fast(ernm_model,
                                           formula = as.formula(paste("nw ~ " ,ERNM_formula)),
                                           theta_draws = 100,
                                           fixed_edges = F,
                                           nsims = 50,
                                           treatment_var = "c.sex",
                                           treatment_indicator = 1,
                                           outcome_var ="c.smoke",
                                           node_treatment_statuses = c(0,1),
                                           cores = 20)
print(proc.time() -t)

#first denotes the nodes, second denotes the neighbours
t <- proc.time()
ERNM_posteriors_both_female <- extract_posteriors_fast(ernm_model,
                                           formula = as.formula(paste("nw ~ " ,ERNM_formula)),
                                           theta_draws = 100,
                                           fixed_edges = F,
                                           nsims = 50,
                                           treatment_var = "c.sex",
                                           treatment_indicator = 1,
                                           outcome_var ="c.smoke",
                                           cores = 20)
print(proc.time() -t)

t <- proc.time()
ERNM_posteriors_male_female <- extract_posteriors_fast(ernm_model,
                                                formula = as.formula(paste("nw ~ " ,ERNM_formula)),
                                                theta_draws = 100,
                                                fixed_edges = F,
                                                nsims = 50,
                                                treatment_var = "c.sex",
                                                treatment_indicator = 1,
                                                outcome_var ="c.smoke",
                                                node_treatment_statuses = c(0),
                                                cores = 20)
print(proc.time() -t)

t <- proc.time()
ERNM_posteriors_female_female <- extract_posteriors_fast(ernm_model,
                                                formula = as.formula(paste("nw ~ " ,ERNM_formula)),
                                                theta_draws = 100,
                                                fixed_edges = F,
                                                nsims = 50,
                                                treatment_var = "c.sex",
                                                treatment_indicator = 1,
                                                outcome_var ="c.smoke",
                                                node_treatment_statuses = c(1),
                                                cores = 20)
print(proc.time() -t)

t <- proc.time()
ERNM_posteriors_both_male <- extract_posteriors_fast(ernm_model,
                                                       formula = as.formula(paste("nw ~ " ,ERNM_formula)),
                                                       theta_draws = 100,
                                                       fixed_edges = F,
                                                       nsims = 50,
                                                       treatment_var = "c.sex",
                                                       treatment_indicator = 0,
                                                       outcome_var ="c.smoke",
                                                       cores = 20)
print(proc.time() -t)

t <- proc.time()
ERNM_posteriors_male_male <- extract_posteriors_fast(ernm_model,
                                                       formula = as.formula(paste("nw ~ " ,ERNM_formula)),
                                                       theta_draws = 100,
                                                       fixed_edges = F,
                                                       nsims = 50,
                                                       treatment_var = "c.sex",
                                                       treatment_indicator = 0,
                                                       outcome_var ="c.smoke",
                                                       node_treatment_statuses = c(0),
                                                       cores = 20)
print(proc.time() -t)

t <- proc.time()
ERNM_posteriors_female_male <- extract_posteriors_fast(ernm_model,
                                                         formula = as.formula(paste("nw ~ " ,ERNM_formula)),
                                                         theta_draws = 100,
                                                         fixed_edges = F,
                                                         nsims = 50,
                                                         treatment_var = "c.sex",
                                                         treatment_indicator = 0,
                                                         outcome_var ="c.smoke",
                                                         node_treatment_statuses = c(1),
                                                         cores = 20)
print(proc.time() -t)




# =========================================================================================================================
# Examine Results
# =========================================================================================================================
apply(ernm_model$theta_sample,2,mean)
apply(ernm_model$theta_sample,2,sd)

# k -peer outcomes
means_out <- do.call(rbind,lapply(list(logistic_posteriors[2:11],
                                       ERGM_posteriors[2:11],
                                       ERNM_posteriors[2:11],
                                       MRF_posteriors[2:11],
                                       MRF_pseudo_posteriors[2:11]),
                                  function(x){sapply(x,mean,na.rm = T)}))
rownames(means_out) <- c("logistic","ERGM+Logistic","ERNM","MRF","pseudo MRF")
means_out

sd_out <- do.call(rbind,lapply(list(logistic_posteriors[2:11],
                                       ERGM_posteriors[2:11],
                                       ERNM_posteriors[2:11],
                                       MRF_posteriors[2:11],
                                       MRF_pseudo_posteriors[2:11]),
                                  function(x){sapply(x,sd,na.rm = T)}))
rownames(sd_out) <- c("logistic","ERGM+Logistic","ERNM","MRF","pseudo MRF")
sd_out

# k-peer "treatment" (gender is treatment)
means_treat <- do.call(rbind,lapply(list(logistic_posteriors[12:21],
                                       ERGM_posteriors[12:21],
                                       ERNM_posteriors[12:21],
                                       MRF_posteriors[12:21],
                                       MRF_pseudo_posteriors[12:21]),
                                  function(x){sapply(x,mean,na.rm = T)}))
rownames(means_treat) <- c("logistic","ERGM+Logistic","ERNM","MRF","pseudo MRF")
means_treat

sd_treat <- do.call(rbind,lapply(list(logistic_posteriors[12:11],
                                    ERGM_posteriors[12:21],
                                    ERNM_posteriors[12:21],
                                    MRF_posteriors[12:21],
                                    MRF_pseudo_posteriors[12:21]),
                               function(x){sapply(x,sd,na.rm = T)}))
rownames(sd_treat) <- c("logistic","ERGM+Logistic","ERNM","MRF","pseudo MRF")
sd_treat

# compare male/female/both for ERNM for "female neighs"
means_female_neighs <- do.call(rbind,lapply(list(ERNM_posteriors_both_female[12:21],ERNM_posteriors_male_female[12:21],ERNM_posteriors_female_female[12:21]),function(x){sapply(x,mean,na.rm = T)}))
rownames(means_female_neighs) <- c("Both","Female","Male")
means_female_neighs

sd_female_neighs <- do.call(rbind,lapply(list(ERNM_posteriors_both_female[12:21],ERNM_posteriors_male_female[12:21],ERNM_posteriors_female_female[12:21]),function(x){sapply(x,sd,na.rm = T)}))
rownames(sd_female_neighs) <- c("Both","Female","Male")
sd_female_neighs

# compare male/female/both for ERNM
means_male_neighs <- do.call(rbind,lapply(list(ERNM_posteriors_both_male[12:21],ERNM_posteriors_male_male[12:21],ERNM_posteriors_female_male[12:21]),function(x){sapply(x,mean,na.rm = T)}))
rownames(means_male_neighs) <- c("Both","Female","Male")
means_male_neighs

sd_male_neighs <- do.call(rbind,lapply(list(ERNM_posteriors_both_male[12:21],ERNM_posteriors_male_male[12:21],ERNM_posteriors_female_male[12:21]),function(x){sapply(x,sd,na.rm = T)}))
rownames(sd_male_neighs) <- c("Both","Female","Male")
sd_male_neighs

# plot the results
mean_compare_plot(means_out,sd_out,"Comparison Plot of MRF, ERGM, Logistic, ERNM k-peer-outcome effects")
mean_compare_plot(means_treat,sd_treat,"Comparison Plot of MRF, ERGM, Logistic, ERNM k-peer-treatment effects")

mean_compare_plot(means_male_neighs,sd_male_neighs,"Comparison Plot of ERNM peer treatments for males,females, both, for male neigbhours")
mean_compare_plot(means_female_neighs,sd_female_neighs,"Comparison Plot of ERNM peer treatments for males,females, both, for female neighbours")

# Next Goal:
# - make facetted plot that shows the k-peer neighbor gender effect for male/female nodes, and male/female neighbors.



# ==========================================================================================================================
# Goodness of Fit :
# ======================================================================================================================
outcome_var <- "c.smoke"
samp <- ernm_model$theta_sample[sample(1:1000,200),]
tmp_model <- tmp_model <- ernm(as.formula(paste("nw ~",ERNM_formula)),
                               maxIter = 2,
                               mcmcSampleSize = 100,
                               mcmcBurnIn = 100)

ernm_sims <- lapply(1:200,function(i){
  theta <- samp[i,]
  tmp_model$m$setThetas(theta)
  sims <- tmp_model$m$sampler$generateSample(10000,1000,1)
  return(binaryNet_to_network(sims[[1]]))
})
print("ERNM sims done")

samp <- ergm_model$ergm$theta_sample[sample(1:1000,200),]
tmp_model <- tmp_model <- ernm(as.formula(paste("nw ~",ERGM_formula)),
                               maxIter = 2,
                               mcmcSampleSize = 100,
                               mcmcBurnIn = 100)

preds <- posterior_predict(logistic_model)
ergm_sims <- lapply(1:200,function(i){
  theta <- samp[i,]
  tmp_model$m$setThetas(theta)
  sims <- tmp_model$m$sampler$generateSample(10000,1000,1)
  net <- binaryNet_to_network(sims[[1]])
  set.vertex.attribute(net,outcome_var,preds[i,])
  return(net)
})
print("ERGM sims done")

preds <- posterior_predict(logistic_model)
logistic_sims <- lapply(1:200,function(i){
  net <- nw
  set.vertex.attribute(net,outcome_var,preds[i,])
  return(net)
})

mrf_tmp_model <- ernm(as.formula(paste("nw ~ " ,MRF_formula)),
                      maxIter = 2,
                      nodeSamplingPercentage = 1)

samp <- mrf_model$theta_sample[sample(1:1000,200),]
mrf_sims <- lapply(1:200,function(i){
  theta <- samp[i,]
  mrf_tmp_model$m$setThetas(theta)
  sims <- mrf_tmp_model$m$sampler$generateSample(10000,1000,1)
  net <- binaryNet_to_network(sims[[1]])
  return(net)
})
print("MRF sims done")
  
samp <- mrf_pseudo_model$theta_sample[rep(1,200),]
mrf_pseudo_sims <- lapply(1:200,function(i){
  theta <- samp[i,]
  mrf_tmp_model$m$setThetas(theta)
  sims <- mrf_tmp_model$m$sampler$generateSample(10000,1000,1)
  net <- binaryNet_to_network(sims[[1]])
  return(net)
})
print("pseudo MRF sims done")



# Get the GOF - function in utils.R file
gof_ERNM <- get_gof(ernm_sims)
gof_ERGM <- get_gof(ergm_sims)
gof_logistic <- get_gof(logistic_sims)
gof_mrf <- get_gof(mrf_sims)
gof_mrf_pseudo <- get_gof(mrf_pseudo_sims)
gof_obs <- get_gof(list(nw))

# see the effect of proportion of treated friends on smoking:
prop_treated_female_neighs <- lapply(list(ernm_sims,
                                   ergm_sims,
                                   logistic_sims,
                                   mrf_sims,
                                   mrf_pseudo_sims),
                                   FUN = function(sims){
                                     print("here")
                                     extract_prop_treated_friends(sims,
                                                                  treatment_var = "c.sex",
                                                                  outcome_var = "c.smoke",
                                                                  treatment_indicator =  1)})
prop_treated_male_neighs <- lapply(list(ernm_sims,
                                          ergm_sims,
                                          logistic_sims,
                                          mrf_sims,
                                          mrf_pseudo_sims),
                                     FUN = function(sims){
                                       print("here")
                                       extract_prop_treated_friends(sims,
                                                                    treatment_var = "c.sex",
                                                                    outcome_var = "c.smoke",
                                                                    treatment_indicator =  0)})
names(prop_treated_male_neighs) = c("ERNM","ERGM","Logistic","MRF","pseudo_MRF")
names(prop_treated_female_neighs) = c("ERNM","ERGM","Logistic","MRF","pseudo_MRF")

# make smaller
prop_treated_female_neighs = prop_treated_female_neighs[1]
prop_treated_male_neighs = prop_treated_male_neighs[1]

# function to make facetted plot for the each gender proportion treated.
prop_treated_hist <- function(output_female,
                              output_male){
  # remove neighs
  
  female_female = output_female[output_female$treats ==1,]
  female_female$neigh_treat = "Female Neighbours"
  female_female$node_treat = "Female Nodes"
  
  male_female = output_female[output_female$treats ==0,]
  male_female$neigh_treat = "Female Neighbours"
  male_female$node_treat = "Male Nodes"
  
  female_male = output_male[output_male$treats ==0,]
  female_male$neigh_treat = "Male Neighbours"
  female_male$node_treat = "Female Nodes"
  
  male_male = output_male[output_male$treats ==1,]
  male_male$neigh_treat = "Male Neighbours"
  male_male$node_treat = "Male Nodes"
  
  tmp <- rbind(female_female,
               male_female,
               female_male,
               male_male)
  
  plot = ggplot(data = tmp,aes(x = prop_treated_neighs,y = outs))+
    facet_grid(neigh_treat ~ node_treat) +
    stat_summary_bin(fun = "mean",
                     geom="bar",
                     binwidth=0.1)+
    xlab("Proportion of gender neighbours")+
    ylab("Proportion of smokers")
  plot
}

prop_treated_hist(prop_treated_female_neighs[[1]],
                prop_treated_male_neighs[[1]])



tmp = prop_treated_neighs[[1]]
plot = ggplot(data = tmp,aes(x = prop_treated_neighs,y = outs))+
  stat_summary_bin(fun = "mean",
                   geom="bar",
                   binwidth=0.1)
plot

plot = ggplot(data = tmp[tmp$treats ==0,],aes(x = prop_treated_neighs,y = outs))+
  stat_summary_bin(fun = "mean",
                   geom="bar",
                   binwidth=0.1)
plot

# # get the observed values
# geodist_obs <- ergm.geodistdist(nw)[1:15]
# esp_obs <- ergm::summary_formula(nw ~ esp(0:20))
# degree_obs <- ergm::summary_formula(nw ~ degree(0:20))
# smoker_degree_obs <- ergm::summary_formula(nw ~ degree(0:20,by = "c.smoke"))[22:42]
# triads_obs <- ergm::summary_formula(nw ~ triangles)
# edges_obs <- ergm::summary_formula(nw ~ edges)
# smoker_triads_obs <- ergm::summary_formula(nw ~ triangles("c.smoke"))
# smoker_edges_obs <- ergm::summary_formula(nw ~ nodematch("c.smoke"))


# plot the goodness of fit
# plot GOF boxplots:

lolog.catelog.helper::lolog_boxplot(sims_list = list(gof_ERNM$geodist,gof_ERGM$geodist),
                                    names_list = list(1:15,1:15),
                                    observed_list = list(gof_obs$geodist),
                                    title_list = list("ERNM","ERGM + Logistic"),
                                    superimpose =T,
                                    title = "Geodist GOF")
lolog.catelog.helper::lolog_boxplot(sims_list = list(gof_ERNM$esp,gof_ERGM$esp),
                                    names_list = list(0:14,0:14),
                                    observed_list = list(gof_obs$esp),
                                    title_list = list("ERNM","ERGM + Logistic"),
                                    superimpose =T,
                                    title = "ESP GOF")
lolog.catelog.helper::lolog_boxplot(sims_list = list(gof_ERNM$degree,gof_ERGM$degree),
                                    names_list = list(0:14,0:14),
                                    observed_list = list(gof_obs$degree),
                                    title_list = list("ERNM","ERGM + Logistic"),
                                    superimpose =T,
                                    title = "Degree GOF")

# Smoker triads and dyads:

par(mfcol = c(1,2))
hist(gof_ERNM$smoker_triads,main = "ERNM Smoker Triad Dist",xlab = "Triads",xlim = c(100,600))
abline(v=  gof_obs$smoker_triads,col = "red")
hist(gof_ERNM$smoker_edges,main = "ERNM Smoker Dyad Dist",xlab = "Dyads")
abline(v=  gof_obs$smoker_edges,col = "red")

# Triads and edges:

par(mfcol = c(1,2))
hist(gof_ERNM$triads,main = "ERNM Triad Dist",xlab = "Triads",xlim = c(300,1200))
abline(v=  gof_obs$triads,col = "red")
hist(gof_ERNM$edges,main = "ERNM edge Dist",xlab = "Triads")
abline(v=  gof_obs$edges,col = "red")

# Proportion of Smoker Triads

par(mfcol = c(1,2))
hist(gof_ERNM$smoker_triads/gof_ERNM$triads,main = "ERNM Smoker Triad Dist",xlab = "Triads")
abline(v=  gof_obs$smoker_triads/gof_obs$triads,col = "red")
hist(gof_ERNM$smoker_edges/gof_ERNM$edges,main = "ERNM Smoker Dyad Dist",xlab = "Dyads")
abline(v=  gof_obs$smoker_edges/gof_obs$edges,col = "red")


# =========================================================================================================================
# Save Results
# =========================================================================================================================
# thin the ERGM_posteriors

# make some null values that we don't have yet:
setwd(SAVEDIR)
tmp <- c("ERNM_posteriors","ERGM_posteriors","MRF_posteriors","MRF_pseudo_posteriors","logistic_posteriors",
         "all_effects",
         "prop_treated_female_neighs","prop_treated_male_neighs",
         "ernm_pilot","ernm_model","mrf_pilot","mrf_model","ergm_model","logistic_model",
         "gof_ERNM","gof_ERGM","gof_logistic","gof_obs","gof_mrf","gof_mrf_pseudo")
all_effects_small <- all_effects
all_effects_small$raw_results <- NULL
ERGM_posteriors_small <- lapply(ERGM_posteriors,function(x){sample(x,length(x)*0.01,replace = F)})
tmp_small <- c("ERNM_posteriors","ERGM_posteriors_small","MRF_posteriors","MRF_pseudo_posteriors","logistic_posteriors",
         "all_effects_small",
         "prop_treated_female_neighs","prop_treated_male_neighs",
         "ernm_pilot","ernm_model","mrf_pilot","mrf_model","mrf_pseudo_model","ergm_model","logistic_model",
         "gof_ERNM","gof_ERGM","gof_logistic","gof_mrf","gof_mrf_pseudo","gof_obs")

save(file = "results_add_health_fitting_2.RData",list = tmp)
save(file = "results_add_health_fitting_2_small.RData",list = tmp_small)

for(i in ls()){
  size = object.size(eval(parse(text = i)))*(10**(-6))
  if(size >2){
    print(i)
    print(size)
  }
}



