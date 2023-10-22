# This script does ERNM Pseudo likelihood estimation:

ernm_MPLE <- function(formula){
  
  net <- eval(parse(text = (formula[[2]])))
  n <- net %n% "n"
  tmp_model <- ernm(formula,
                    maxIter = 2,
                    mcmcSampleSize = 1000,
                    mcmcBurnIn = 100)
  tmp_model <- tmp_model$m$sampler$getModel()
  
  obs_stats <- tmp_model$statistics()
  stats <- list()
  for(i in 1:n){
    for(j in 1:n){
      # update the network and calculate the change stats:
      tmp <- tryCatch({
        tmp_model$dyadUpdate(i,j)
        output <- (-1)**(net[i,j]) * (tmp_model$statistics() - obs_stats)
        tmp_model$calculate()
        output},
        error = function(e){
          return(NULL)
        }
      )
      if(!is.null(tmp)){
        stats[[(length(stats)+1)]] <- tmp
        stats[[(length(stats))]][(length(output)+1)] = net [i,j]
      }
    }
    if((n*(i-1) + j)>10000){
      break
    }
  }
  
  # return logistic regression on change stats:
  stats <- as.data.frame(do.call(rbind,stats))
  tmp <- which(names(stats) == "??")
  names(stats)[tmp] <- 1:length(tmp)
  names(stats)[dim(stats)[2]] <- "edge_present"
  stats <- stats[-which(names(stats) == "edges")]
  tmp <- glm(edge_present ~., data = stats)
}