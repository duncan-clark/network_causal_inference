# This script carries out the triad census and makes the GOF plots contained in the paper:

library(statnet)
library(ernm)
library(mvtnorm)
library(coda)
library(parallel)
library(doParallel)
library(rstanarm)
source("functions/utils.R")
source("functions/fire.R")

load("/results/add_health/results_add_health_fitting.RData")

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
outcome_var <- "c.smoke"

# =======================================================================================
# 1) Logistic Regression Simulations
# =======================================================================================
preds <- posterior_predict(logistic_model)
sims_1 <- lapply(1:1000,function(i){
  net <- nw
  set.vertex.attribute(net,outcome_var,preds[i,])
  return(net)
})

# =======================================================================================
# 2) ERNM Simulations
# =======================================================================================
ERNM_formula <- "edges() + 
                gwesp(0.5,homogenous = T,variableName = 'c.grade') + 
                gwdegree(0.5) + 
                homophily('c.grade') +
                homophily('c.sex') +
                homophily('c.smoke') +
                nodeCount('c.smoke') +
                logisticNeighborsTopLevel('c.smoke','c.smoke') +
                logistic('c.smoke','c.sex') | c.smoke"

tmp_model <- ernm(as.formula(paste("nw ~", ERNM_formula)),
                  maxIter = 2,
                  mcmcSampleSize = 100,
                  mcmcBurnIn = 100)

samp <- ernm_model$theta_sample[sample(dim(ernm_model$theta_sample)[1],1000,replace = T),]

sims_2 <- lapply(1:1000,function(i){
  theta <- samp[i,]
  tmp_model$m$setThetas(theta)
  sims <- tmp_model$m$sampler$generateSample(10000,1000,1)
  return(binaryNet_to_network(sims[[1]]))
})


# =======================================================================================
# 3) ERGM + GLM Simulations
# =======================================================================================
ERGM_formula <- "edges() + gwesp(0.5) + gwdegree(0.5)  + homophily('c.grade') + homophily('c.sex') +  homophily('c.smoke')"

samp <- ergm_model$ergm$theta_sample[sample(dim(ergm_model$ergm$theta_sample)[1],1000,replace = T),]
tmp_model <- ernm(as.formula(paste("nw ~", ERGM_formula)),
                  maxIter = 2,
                  mcmcSampleSize = 100,
                  mcmcBurnIn = 100)

# do the network sims
sims_3 <- lapply(1:1000,function(i){
  theta <- samp[i,]
  tmp_model$m$setThetas(theta)
  sims <- tmp_model$m$sampler$generateSample(10000,1000,1)
  return(binaryNet_to_network(sims[[1]]))
})

preds <- posterior_predict(ergm_model$glm)
# do the logistic sims
sims_3 <- lapply(1:1000,function(i){
  net <- sims_3[[i]]
  set.vertex.attribute(net,outcome_var,preds[i,])
  return(net)
})

# Save the simulation:

# =======================================================================================
# Compare distribuitons
# =======================================================================================

# get the gof:
gof_logistic <- get_gof(sims_1)
gof_ernm <- get_gof(sims_2)
gof_ergm <- get_gof(sims_3)
gof_obs <- get_gof(list(nw))

par(mfcol = c(1,2))
hist(gof_ernm$triads,main = "ERNM Triad Dist",xlab = "Triads",xlim = c(200,1200))
abline(v = gof_obs$triads,col = "red")
hist(gof_ernm$triads,main = "ERGM + GLM Triad Dist",xlab = "Triads",xlim = c(200,1200))
abline(v =  gof_obs$triads,col = "red")

# smokers
par(mfcol = c(1,3))
hist(sapply(sims_1,function(net){sum(get.vertex.attribute(net,outcome_var))}),main = "Logistic Smoker Dist",xlab = "Smokers")
abline(v = sum(as.numeric(get.vertex.attribute(nw,outcome_var))),col = "red")
hist(sapply(sims_2,function(net){sum(as.numeric(get.vertex.attribute(net,outcome_var)))}),main = "ERNM Smoker Dist",xlab = "Smokers")
abline(v = sum(as.numeric(get.vertex.attribute(nw,outcome_var))),col = "red")
hist(sapply(sims_3,function(net){sum(get.vertex.attribute(net,outcome_var))}),main = "ERGM + GLM Smoker Dist",xlab = "Smokers")
abline(v = sum(as.numeric(get.vertex.attribute(nw,outcome_var))),col = "red")



# smoker  triads:
par(mfcol = c(1,2))
hist(gof_ernm$smoker_triads,main = "ERNM Triad Dist",xlab = "Triads",xlim = c(100,550))
abline(v=  gof_obs$smoker_triads,col = "red")
hist(gof_ergm$smoker_triads,main = "ERGM + GLM Triad Dist",xlab = "Triads",xlim = c(100,550))
abline(v=  gof_obs$smoker_triads,col = "red")

par(mfcol = c(1,2))
hist(gof_ernm$smoker_triads/gof_ernm$triads,main = "ERNM Triad Dist",xlab = "Triads",xlim = c(0,1))
abline(v=  gof_obs$smoker_triads/gof_obs$triads,col = "red")
hist(gof_ergm$smoker_triads/gof_ernm$triads,main = "ERGM + GLM Triad Dist",xlab = "Triads",xlim = c(0,1))
abline(v=  gof_obs$smoker_triads/gof_obs$triads,col = "red")


par(mfcol = c(1,2))
hist(gof_ernm$smoker_edges,main = "ERNM Triad Dist",xlab = "edges",xlim = c(200,550))
abline(v=  gof_obs$smoker_edges,col = "red")
hist(gof_ergm$smoker_edges,main = "ERGM + GLM Triad Dist",xlab = "edges",xlim = c(200,550))
abline(v=  gof_obs$smoker_edges,col = "red")


# Comments : 
# - ERNM is doing slightly better than ERGM + logistic regression
# - Still does not do well on the triad fit


# plot GOF boxplots:
library(lolog.catelog.helper)
library(ggplot2)
library(scales)
lolog_boxplot(sims_list = list(gof_ernm$esp,gof_ergm$esp),
                                    names_list = list(0:20,0:20),
                                    observed_list = list(gof_obs$esp,gof_obs$esp),
                                    title_list = list("ERNM","ERGM+Logistic"),
                                    superimpose =T,
                                    title = "ESP GOF")

# plot GOF boxplots:
library(lolog.catelog.helper)
library(ggplot2)
library(scales)
lolog_boxplot(sims_list = list(gof_ernm$degree,gof_ergm$degree),
                                    names_list = list(1:20,1:20),
                                    observed_list = list(gof_obs$degree,gof_obs$degree),
                                    title_list = list("ERNM","ERGM+Logistic"),
                                    superimpose =T,
                                    title = "Degree GOF")

# plot GOF boxplots:
library(lolog.catelog.helper)
library(ggplot2)
library(scales)
lolog_boxplot(sims_list = list(gof_ernm$smoker_degree,gof_ergm$smoker_degree),
                                    names_list = list(1:20,1:20),
                                    observed_list = list(gof_obs$smoker_degree,gof_obs$smoker_degree),
                                    title_list = list("ERNM","ERGM+Logistic"),
                                    superimpose =T,
                                    title = "Smoker degree GOF")
# plot GOF boxplots:
library(lolog.catelog.helper)
library(ggplot2)
library(scales)
lolog_boxplot(sims_list = list(gof_ernm$geodist,gof_ergm$geodist),
                                    names_list = list(1:15,1:15),
                                    observed_list = list(gof_obs$geodist,gof_obs$geodist),
                                    title_list = list("ERNM","ERGM+Logistic"),
                                    superimpose =T,
                                    title = "Geodesic Distance GOF")

save(file = "GOF_sims.RData",list = c("nw",
                                      "gof_ergm",
                                      "gof_ernm",
                                      "gof_logistic",
                                      "gof_obs"))
