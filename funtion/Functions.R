################################################################################
# N: size of grid
# M: the number of realization
# x.sample: initial data.independent
# y.sample: initial data.dependent
# g.sample: initial data.constraint
# g_constraint: constraint threshold
# d: d-dimension of input variable
# lb: lower bound - vector
# ub: upper bound - vector
################################################################################  
source("Scales.R")
library(lhs)
library(SLHD)
library(DiceKriging)
library(MASS)
library(mvtnorm)
################################################################################

################################################################################
##                           Constrained Problem                              ##
################################################################################

################### Constrained GP Optimal In Predictions ######################
gp_opt_mean_constraint <- function(x.sample, y.sample, g.sample, g_constraint, lb, ub, d){
  eps <- sqrt(.Machine$double.eps)
  d <- ncol(x.sample)
  xF_u <- scale_to_uni(x.sample, lb, ub)
  xF_o <- x.sample
  
  gpm <- km(design=data.frame(x=xF_u), 
            response=data.frame(y=y.sample), nugget = eps)
  
  gpm_g <- km(design=data.frame(x=xF_u), 
              response=data.frame(y=g.sample), nugget = eps)
  
  x.grid <- randomLHS(10000, d)
  prd <- predict(gpm, newdata=data.frame(x=x.grid), checkNames=FALSE, type="SK", cov.compute = TRUE)
  prd_g <- predict(gpm_g, newdata=data.frame(x=x.grid), checkNames=FALSE, type="SK", cov.compute = TRUE)
  y.gphat <- min(prd$mean[prd_g$mean >= g_constraint])
  index_y.gphat <- which(prd$mean == y.gphat & prd_g$mean > g_constraint)
  
  if (d == 1) {
    xu.gphat <- x.grid[index_y.gphat]
  } else {
    xu.gphat <- x.grid[index_y.gphat, ]
  }
  
  gp_mean_hat_x <- scale_to_org(xu.gphat, lb, ub)
  
  return(gp_mean_hat_x=gp_mean_hat_x)
  
}
################################################################################

#################### Constrained GP Optimal Decision Sets ######################
gp_opts_constraint <- function(x.sample, y.sample, g.sample, g_constraint, lb, ub, N, M, d){
  eps <- sqrt(.Machine$double.eps)
  d <- ncol(x.sample)
  xF_u <- scale_to_uni(x.sample, lb, ub)
  xF_o <- x.sample
  
  gpm <- km(design=data.frame(x=xF_u), 
            response=data.frame(y=y.sample), nugget = eps)
  
  gpm_g <- km(design=data.frame(x=xF_u), 
              response=data.frame(y=g.sample), nugget = eps)
  
  gpm_callable <- function(new_data) {
    new_data_u <- scale_to_uni(new_data, lb, ub)
    predict(gpm, newdata = data.frame(x = new_data_u), 
            type = "SK", checkNames = FALSE)
  }
  
  gpm_g_callable <- function(new_data) {
    new_data_u <- scale_to_uni(new_data, lb, ub)
    predict(gpm_g, newdata = data.frame(x = new_data_u), 
            type = "SK", checkNames = FALSE)
  }
  
  gp.xhat <- matrix(NA, ncol = d, nrow = M)
  gp.yhat <- c()
  gp.ghat <- c()
  
  for (i in 1:M) {
    x.grid <- randomLHS(N, d)
    prd <- predict(gpm, newdata=data.frame(x=x.grid), checkNames=FALSE, type="SK", cov.compute = TRUE)
    prd_g <- predict(gpm_g, newdata=data.frame(x=x.grid), checkNames=FALSE, type="SK", cov.compute = TRUE)
    posterior.y <- rmvnorm(1, mean = prd$mean, sigma = prd$cov)
    posterior.g <- rmvnorm(1, mean = prd_g$mean, sigma = prd_g$cov)
    y.gphat <- min(posterior.y[posterior.g >= g_constraint])
    index_y.gphat <- which(posterior.y == y.gphat & posterior.g > g_constraint)
    g.gphat <- posterior.g[index_y.gphat]
    
    if (d == 1) {
      gp.xhat[i] <- scale_to_org(x.grid[index_y.gphat], lb, ub)
      gp.yhat <- c(gp.yhat, y.gphat)
      gp.ghat <- c(gp.ghat, g.gphat)
    } else {
      gp.xhat[i, ] <- scale_to_org(x.grid[index_y.gphat, ], lb, ub)
      gp.yhat <- c(gp.yhat, y.gphat)
      gp.ghat <- c(gp.ghat, g.gphat)
    }
  }
  return(list(gpm_y = gpm_callable, gpm_g = gpm_g_callable, gp.xhat=gp.xhat, gp.yhat=unlist(gp.yhat), gp.ghat=unlist(gp.ghat)))
}
################################################################################

################################################################################
##                           Unconstrained Problem                            ##
################################################################################

################## Unconstrained GP Optimal In Predictions #####################
gp_opt_mean_unconstraint <- function(x.sample, y.sample, lb, ub, d){
  eps <- sqrt(.Machine$double.eps)
  d <- ncol(x.sample)
  xF_u <- scale_to_uni(x.sample, lb, ub)
  xF_o <- x.sample
  
  gpm <- km(design=data.frame(x=xF_u), 
            response=data.frame(y=y.sample), nugget = eps)
  
  x.grid <- randomLHS(10000, d)
  prd <- predict(gpm, newdata=data.frame(x=x.grid), checkNames=FALSE, type="SK", cov.compute = TRUE)
  y.gphat <- min(prd$mean)
  index_y.gphat <- which(prd$mean == y.gphat)
  
  if (d == 1) {
    xu.gphat <- x.grid[index_y.gphat]
  } else {
    xu.gphat <- x.grid[index_y.gphat, ]
  }
  
  gp_mean_hat_x <- scale_to_org(xu.gphat, lb, ub)
  
  return(gp_mean_hat_x=gp_mean_hat_x)
  
}
################################################################################

################### Unconstrained GP Optimal Decision Sets #####################
gp_opts_unconstraint <- function(x.sample, y.sample, lb, ub, N, M, d){
  eps <- sqrt(.Machine$double.eps)
  d <- ncol(x.sample)
  xF_u <- scale_to_uni(x.sample, lb, ub)
  xF_o <- x.sample
  
  gpm <- km(design=data.frame(x=xF_u), 
            response=data.frame(y=y.sample), nugget = eps)
  
  gpm_callable <- function(new_data) {
    new_data_u <- scale_to_uni(new_data, lb, ub)
    predict(gpm, newdata = data.frame(x = new_data_u), 
            type = "SK", checkNames = FALSE)
  }
  
  gp.xhat <- matrix(NA, ncol = d, nrow = M)
  gp.yhat <- c()
  
  for (i in 1:M) {
    x.grid <- randomLHS(N, d)
    prd <- predict(gpm, newdata=data.frame(x=x.grid), checkNames=FALSE, type="SK", cov.compute = TRUE)
    posterior.y <- rmvnorm(1, mean = prd$mean, sigma = prd$cov)
    y.gphat <- min(posterior.y)
    xhat.index <- which.min(posterior.y)
    
    if (d == 1) {
      gp.xhat[i] <- scale_to_org(x.grid[xhat.index], lb, ub)
      gp.yhat <- c(gp.yhat, y.gphat)
    } else {
      gp.xhat[i, ] <- scale_to_org(x.grid[xhat.index, ], lb, ub)
      gp.yhat <- c(gp.yhat, y.gphat)
    }
    
  }
  return(list(gpm_y = gpm_callable, gp.xhat=gp.xhat, gp.yhat=unlist(gp.yhat)))
}
################################################################################