fit_n_models <- function(N,...){
  
  set.seed(92749402)
  
  all_betas <- tribble(
    ~lap, ~method, ~var_id, ~var, ~value, ~is_zero, ~correct_value, ~correct_zero,  ~diff_value, ~is_correct
  )
  
  for(i in 1:N){
  dat <- sim_data(set_seed = FALSE, ...)
  cat("Lap ", i, "/", N, "\n", sep = "")
  # Regular LAD-lasso model
  fit_ladlasso <- ladlasso_cv(Y = dat$Y, X = dat$X,
                              lambda.min = 0.0001, lambda.max = 0.5, B = 100, 
                              print_progress = FALSE)
  # Adaptive LAD-lasso model
  fit_adladlasso <- adladlasso_cv(Y = dat$Y, X = dat$X,
                                  lambda.min = 0.0001, lambda.max = 0.5, B = 100,
                                  beta0 = fit_ladlasso$beta,
                                  print_progress = FALSE)
  
  # Collect betas into single matrix
  beta_mat <- rbind(dat$beta_orig, fit_ladlasso$beta[-1,], fit_adladlasso$beta[-1,])
  
  # Number of methods 
  n_methods <- nrow(beta_mat)/dat$n_beta
  
  # Make the matrix into tibble
  beta_tibble <- as_tibble(beta_mat) |>
    mutate(method = rep(c("original", "ladlasso", "ad_ladlasso"), each = dat$n_beta),
           var_id = rep(1:dat$n_beta, n_methods)) |> 
    gather("var", "value", -method,  -var_id) |> 
    mutate(is_zero = (value > -1e-05) & (value < 1e-05)) |> 
    add_column(lap = i, .before = 1)
  
  # Save the correct values to compare them to models values
  correct_zeroes <- beta_tibble |>
    filter(method == "original") |> 
    pull(is_zero)
  correct_values <- beta_tibble |>
    filter(method == "original") |> 
    pull(value)
  
  # Add new rows to indicate the difference to the real values
  beta_tibble <- 
    beta_tibble |> 
    # Add temporary columns to create other variables
    add_column(correct_value = rep(correct_values, n_methods),
               correct_zero = rep(correct_zeroes, n_methods)) |> 
    # Add difference to original and a variable indicating whether the zero is correct
    mutate(diff_value = abs(correct_value - value),
           is_correct = correct_zeroes == is_zero) #|> 
    # Remove temporary columns
    #select(-correct_zero, -correct_value)
    
  all_betas <- rbind(all_betas, beta_tibble)
  
  }
  
  return(all_betas)
}

beta_bar_all <- function(beta_tibble, bar_width = 0.5, x_names = paste("Beta", 1:dat$n_beta, sep = "")){
  box_x <- c(dat$n_effect + 0.5, dat$n_beta + 0.5)
  box_y <- c(-1,1) * 0.2
  
  beta_tibble |> 
    filter(var == "Y1") |> 
    ggplot() +
    geom_hline(yintercept = 0)+
    geom_bar(mapping = aes(x = factor(var_id), y = value, fill = factor(method, levels = c("original", "ladlasso", "ad_ladlasso"))),
             width = bar_width,
             color = "black",
             position = position_dodge(preserve = "single"), stat = "identity") +
    #geom_text(aes(label = round(value,1)), position = position_dodge(), vjust = -0.5) +
    labs(x = "Estimates", y = "Value", fill = "Method") +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.0,
        hjust = 0.0
      ))+
    scale_x_discrete(labels = x_names) +
    geom_path(data = tibble(x = box_x[c(1,1,2,2,1)], #c(3.5,3.5,10.5,10.5,3.5), #center_x[c(1,1,2,2,1)]
                            y = box_y[c(1,2,2,1,1)]),#c(-0.2,0.2,0.2,-0.2,-0.2)) #center_y[c(1,2,2,1,1)]
              aes(x = x,
                  y = y),
              color = "darkgrey",
              size = 1) 
}

beta_bar_eff <- function(beta_tibble, bar_width = 0.5, x_names = paste("Beta", 1:dat$n_beta, sep = "")){
  beta_tibble_eff <- beta_tibble |> 
    filter(var == "Y1", var_id %in% (dat$n_effect + 1):dat$n_beta)
  
  ggplot(beta_tibble_eff) +
    geom_hline(yintercept = 0)+
    geom_bar(mapping = aes(x = factor(var_id), y = value, fill = factor(method, levels = c("original", "ladlasso", "ad_ladlasso"))),
             color = "black",position = position_dodge(preserve = "single"), stat = "identity",
             width = bar_width) +
    geom_text(aes(x = factor(var_id), y = value, 
                  color = factor(method, levels = c("original", "ladlasso", "ad_ladlasso")), 
                  label = ifelse(is_zero, "0", "")),
              position = position_dodge(width = 0.55), 
              vjust = -0.4,
              size = 5,
              show.legend = FALSE) +
    geom_rect(data = beta_tibble_eff |> filter(is_zero), 
              aes(xmin = var_id - 3 - 0.1 + ifelse(method == "original", -0.19, ifelse(method == "ad_ladlasso", 0.19, 0)),
                  xmax = var_id - 3 + 0.1 + ifelse(method == "original", -0.19, ifelse(method == "ad_ladlasso", 0.19, 0)),
                  fill = method
              ),
              ymin = 0.00014,
              ymax = -0.00014,
              color = "black",
              show.legend = FALSE) +
    scale_x_discrete(labels = x_names[-(1:dat$n_effect)]) +
    scale_y_continuous(breaks = c(-0.015, -0.005, 0.005, 0.015))+
    labs(x = "Estimates",
         y = "",
         fill = "Method") 
}

plot_beta_bar <- function(beta_tibble){
  bar_width <- 0.5
  
  x_names <- paste("Beta", 1:dat$n_beta, sep = "")
  
  ret_plot <- beta_bar_all(beta_tibble, bar_width = bar_width, x_names = x_names) +
    beta_bar_eff(beta_tibble, bar_width = bar_width, x_names = x_names) +
    plot_layout(guides = 'collect', design =
                  "
      AAA#####
      AAABBBBB
      AAABBBBB
      AAA#####
      "
    )
  
  return(ret_plot)  
}

sim_data <- function(n = 100, n_beta = 10, n_effect = 3, n_response = 2, p_zeroes = 0.1, seed = 418630, set_seed = TRUE, error_type = "uniform"){
  # Set seed
  if(set_seed) set.seed(seed)
  
  # Covariates
  beta_orig <- c(rnorm(n_effect * n_response), rep(0, (n_beta - n_effect) * n_response)) |> 
    # Vector to matrix
    matrix(ncol = n_response, nrow = n_beta, byrow = TRUE)
  
  # Explanatory variables
  X <- rnorm(n * n_beta, mean = 3) |> 
    # Add zeroes with the proportion p_zeroes
    {\(x) ifelse(test = (runif(n * n_beta) < p_zeroes) | x < 0, yes = 0, no = x)}() |> 
    # Vector to matrix
    matrix(nrow = n, ncol = n_beta)
  
  # Response variables
  Y <- X %*% beta_orig + 
  # Error
  switch (error_type,
    "uniform" = matrix(runif(n * n_response), ncol = n_response, nrow = n),
    "laplace" = raml(n, mu = c(3,6), Sigma = diag(2)*0.5)
  )
  
  colnames(beta_orig) <- paste("Y", 1:n_response, sep = "")
  
  return(list(Y = Y, X = X, beta_orig = beta_orig, n_effect = n_effect, n_beta = n_beta))
}


## Input variables:
##   Y:         
##   X:         
##   beta0:     The initial value of beta
##   lambda.min: The minimum value in the grid
##   lambda.max: The maximum value in the grid
##   B:         The number of values in the grid
##   K:         The number of folds for cross validation, default is 5.
##   penalty:   The indices of the non-penalized variables
## Output variables:
##   gamma:     The weight vector
##   lambda:    The grid of lambdas
##   PE:        A sequence of prediction error corresponding to each possible value of lambda
##   min_ind:   The index of the lambda which minimizes PE
##   PE_min:    The minimum value of prediction error
##   lambda_op: The optimal lambda value corresponding to PE_min
##   beta:      The adaptive lad-lasso coefficient estimate using gamma and lambda
##   conv_measure1: The value of the convergence measure 1
##   conv_measure2: The value of the convergence measure 2

adladlasso_cv <- function(Y=Y300,X=X300,beta0=beta0,lambda.min=lambda.min,lambda.max=lambda.max,
                          B=50,k=5,penalty=c(1), print_progress = TRUE)
{
  require('cvTools')
  require('MNM')
  
  lambda<-exp(seq(log(lambda.min),log(lambda.max),length=B))
  X <- cbind(1,X)
  p<-ncol(Y)
  q<-ncol(X)
  n<-nrow(Y) 
  npen<-length(penalty)
  C <- cvFolds(n, K = k,type = c("random"))
  qn <- rep(0,B)
  PE <- rep(0,B)
  
  norms<-sqrt(diag(beta0%*%t(beta0)))
  gamma<-1/(norms+1/n)
  
  for (j in 1:B)
  {
    PEi <- rep(0,k)
    qni <- rep(0,k)
    for (i in 1:k)
    {
      tra <- C$subsets[C$which!=i]
      te <- C$subsets[C$which==i]
      Yte <- Y[te,]
      Ytra <- Y[tra,]
      Xte <- X[te,]
      Xtra <- X[tra,]
      ntra <- dim(Ytra)[1] 
      y1<-matrix(0,q-npen,p)
      ytra<-rbind(Ytra,y1)
      x1<-ntra*lambda[j]*diag(gamma)
      x1<-x1[-penalty,]
      xtra<-rbind(Xtra,x1)
      mod1 <-mv.l1lm(ytra~-1+xtra,score="s",stand="o",maxiter = 10000, eps = 1e-4, eps.S = 1e-6)
      beta <- coef(mod1)
      E <- Yte-Xte%*%beta
      PEi[i] <- mean(sqrt(diag(E%*%t(E))))
      norms<-sqrt(diag(beta%*%t(beta)))
      qni[i]<-sum(norms>1e-6)-1
    }
    PE[j] <- mean(PEi)
    qn[j] <- mean(qni)
    if(print_progress) print(paste("j=",j,"PE=",PE[j]))
  }
  
  PE_min <- min(PE)
  lambda_op <- lambda[which.min(PE)]
  min_ind <- which.min(PE)
  
  y1<-matrix(0,q-npen,p)
  y<-rbind(Y,y1)
  x1<-n*lambda_op*diag(gamma)
  x1<-x1[-penalty,]
  x<-rbind(X,x1)
  mod1 <-mv.l1lm(y~-1+x,score="s",stand="o",maxiter = 10000, eps = 1e-4, eps.S = 1e-6)
  beta <- coef(mod1)
  norms<-sqrt(diag(beta%*%t(beta)))
  
  diff<-beta-beta0
  osoittaja<-sqrt(sum(diag(t(diff)%*%diff)))
  nimittaja<-sqrt(sum(diag(t(beta0)%*%beta0)))
  conv_measure1<-osoittaja/nimittaja
  
  diff<-c(diff)
  osoittaja <- sqrt(sum(diff*diff))
  nimittaja <- sqrt(sum(beta0*beta0))
  conv_measure2<-osoittaja/nimittaja
  
  plot(lambda,PE,type="l",xlab="lambda",ylab="PE")
  abline(v=lambda_op,lty=2)
  
  #plot(norms[-1],type="h",xlab="marker",ylab=bquote(paste("Marker effect ||",beta[j],"||")))
  #points(c(50,75,100,150),c(-10,-10,-10,-10),type="h",col="red")
  
  return(list(gamma=gamma,lambda=lambda,PE=PE, min_ind=min_ind, PE_min=PE_min, 
              lambda_op = lambda_op, beta=beta, conv_measure1=conv_measure1, 
              conv_measure2=conv_measure2))  
}

## Input variables:
##   Y:         
##   X:         
##   lambda.min: The minimum value in the grid
##   lambda.max: The maximum value in the grid
##   len:         The number of values in the grid
##   penalty:   The indices of the non-penalized variables
## Output variables:
##   beta:      The estimated coefficient matrix
##   lambda:    The grid of lambdas
##   bic:        A sequence of cp values corresponding to each possible value of lambda
##   min_ind:   The index of the lambda which corresponds to the smallest PE
##   bic_min:   The minimum value of bic
##   lambda_op: The optimal lambda value corresponding to bic_min

adladlasso_bic <- function(Y=Y,X=X,beta0=beta0,lambda.min=lambda.min,lambda.max=lambda.max,len=50,penalty=c(1),sigma=3, print_progress = TRUE)
{
  require('MNM')
  
  lambda<-exp(seq(log(lambda.min),log(lambda.max),length=len))
  X <- cbind(1,X)
  p<-ncol(Y)
  q<-ncol(X)
  n<-nrow(Y) 
  npen<-length(penalty)
  
  value<-NULL
  h<-NULL
  bic<-NULL
  
  norms<-sqrt(diag(beta0%*%t(beta0)))
  gamma<-1/(norms+1/n)
  
  for (i in 1:len)
  {
    y1<-matrix(0,q-npen,p)
    y<-rbind(Y,y1)
    x1<-n*lambda[i]*diag(gamma)
    x1<-x1[-penalty,]
    x<-rbind(X,x1)
    mod1 <-mv.l1lm(y~-1+x,score="s",stand="o",maxiter = 10000, eps = 1e-4, eps.S = 1e-6)
    beta <- coef(mod1)
    E <- Y-X%*%beta
    E<-Y-X%*%beta   
    value[i]<-mean(sqrt(diag(E%*%t(E))))
    norms<-sqrt(diag(beta%*%t(beta)))
    h[i]<-sum(norms>1e-6)
    bic[i]<-value[i]/sigma+h[i]*log(n)/n
    if(print_progress) print(paste("i=",i,"value=",value[i],"h=",h[i],"bic=",bic[i]))
  }
  bic_min <- min(bic)
  lambda_op <- lambda[which.min(bic)]
  min_ind<-which.min(bic)
  
  y1<-matrix(0,q-npen,p)
  y<-rbind(Y,y1)
  x1<-n*lambda_op*diag(gamma)
  x1<-x1[-penalty,]
  x<-rbind(X,x1)
  mod1 <-mv.l1lm(y~-1+x,score="s",stand="o",maxiter = 10000, eps = 1e-4, eps.S = 1e-6)
  beta <- coef(mod1)
  
  plot(lambda,bic,type="l")
  abline(v=lambda_op)
  
  return(list(beta=beta,lambda=lambda, bic=bic, bic_min=bic_min, min_ind=min_ind, lambda_op = lambda_op))  
}

## Input variables:
##   Y:         
##   X:         
##   lambda.min: The minimum value in the grid
##   lambda.max: The maximum value in the grid
##   B:         The number of values in the grid
##   k:         The number of folds for cross validation, default is 5.
##   penalty:   The indices of the non-penalized variables
## Output variables:
##   beta:      The estimated coefficient matrix
##   qn:        
##   lambda:    The grid of lambdas
##   PE:        A sequence of prediction error corresponding to each possible value of lambda
##   min_ind:   The index of the lambda which corresponds to the smallest PE
##   PE_min:    The minimum value of prediction error
##   lambda_op: The optimal lambda value corresponding to PE_min

ladlasso_cv <- function(Y=Y300,X=X300,lambda.min=lambda.min,lambda.max=lambda.max,B=50,k=5,penalty=c(1), print_progress = TRUE)
{
  require('cvTools')
  require('MNM')
  
  lambda<-exp(seq(log(lambda.min),log(lambda.max),length=B))
  X <- cbind(1,X)
  p<-ncol(Y)
  q<-ncol(X)
  n<-nrow(Y) 
  npen<-length(penalty)
  C <- cvFolds(n, K = k,type = c("random"))
  qn <- rep(0,B)
  PE <- rep(0,B)
  
  for (j in 1:B)
  {
    PEi <- rep(0,k)
    qni <- rep(0,k)
    for (i in 1:k)
    {
      tra <- C$subsets[C$which!=i]
      te <- C$subsets[C$which==i]
      Yte <- Y[te,]
      Ytra <- Y[tra,]
      Xte <- X[te,]
      Xtra <- X[tra,]
      ntra <- dim(Ytra)[1] 
      y1<-matrix(0,q-npen,p)
      ytra<-rbind(Ytra,y1)
      x1<-ntra*lambda[j]*diag(q)
      x1<-x1[-penalty,]
      xtra<-rbind(Xtra,x1)
      mod1 <-mv.l1lm(ytra~-1+xtra,score="s",stand="o",maxiter = 10000, eps = 1e-4, eps.S = 1e-6)
      beta <- coef(mod1)
      E <- Yte-Xte%*%beta
      PEi[i] <- mean(sqrt(diag(E%*%t(E))))
      norms<-sqrt(diag(beta%*%t(beta)))
      qni[i]<-sum(norms>1e-6)-1
    }
    PE[j] <- mean(PEi)
    qn[j] <- mean(qni)
    if(print_progress) print(paste("j=",j,"PE=",PE[j], "lambda=",lambda[j]))
  }
  
  PE_min <- min(PE)
  lambda_op <- lambda[which.min(PE)]
  min_ind<-which.min(PE)
  
  y1<-matrix(0,q-npen,p)
  y<-rbind(Y,y1)
  x1<-n*lambda_op*diag(q)
  x1<-x1[-penalty,]
  x<-rbind(X,x1)
  mod1 <-mv.l1lm(y~-1+x,score="s",stand="o",maxiter = 10000, eps = 1e-4, eps.S = 1e-6)
  beta <- coef(mod1)
  norms<-sqrt(diag(beta%*%t(beta)))
  
  plot(lambda,PE,type="l",xlab="lambda",ylab="PE")
  abline(v=lambda_op,lty=2)
  
  #plot(norms[-1],type="h",xlab="marker",ylab=bquote(paste("Marker effect ||",beta[j],"||")))
  #points(c(50,75,100,150),c(-10,-10,-10,-10),type="h",col="red")
  
  return(list(beta=beta,qn=qn,lambda=lambda,PE=PE, min_ind=min_ind, PE_min=PE_min, 
              lambda_op = lambda_op))  
}



## Input variables:
##   Y:         
##   X:         
##   lambda.min: The minimum value in the grid
##   lambda.max: The maximum value in the grid
##   len:         The number of values in the grid
##   penalty:   The indices of the non-penalized variables
## Output variables:
##   beta:      The estimated coefficient matrix
##   lambda:    The grid of lambdas
##   bic:        A sequence of cp values corresponding to each possible value of lambda
##   min_ind:   The index of the lambda which corresponds to the smallest PE
##   bic_min:   The minimum value of bic
##   lambda_op: The optimal lambda value corresponding to bic_min

ladlasso_bic <- function(Y=Y,X=X,lambda.min=lambda.min,lambda.max=lambda.max,len=50,penalty=c(1),sigma=3, print_progress = TRUE)
{
  require('MNM')
  
  lambda<-exp(seq(log(lambda.min),log(lambda.max),length=len))
  X <- cbind(1,X)
  p<-ncol(Y)
  q<-ncol(X)
  n<-nrow(Y) 
  npen<-length(penalty)
  
  value<-NULL
  h<-NULL
  bic<-NULL
  
  for (i in 1:len)
  {
    y1<-matrix(0,q-npen,p)
    y<-rbind(Y,y1)
    x1<-n*lambda[i]*diag(q)
    x1<-x1[-penalty,]
    x<-rbind(X,x1)
    mod1 <-mv.l1lm(y~-1+x,score="s",stand="o",maxiter = 10000, eps = 1e-4, eps.S = 1e-6)
    beta <- coef(mod1)
    E <- Y-X%*%beta
    E<-Y-X%*%beta   
    value[i]<-mean(sqrt(diag(E%*%t(E))))
    norms<-sqrt(diag(beta%*%t(beta)))
    h[i]<-sum(norms>1e-6)
    bic[i]<-value[i]/sigma+h[i]*log(n)/n
    if(print_progress) print(paste("i=",i,"value=",value[i],"h=",h[i],"bic=",bic[i]))
  }
  bic_min <- min(bic)
  lambda_op <- lambda[which.min(bic)]
  min_ind<-which.min(bic)
  
  y1<-matrix(0,q-npen,p)
  y<-rbind(Y,y1)
  x1<-n*lambda_op*diag(q)
  x1<-x1[-penalty,]
  x<-rbind(X,x1)
  mod1 <-mv.l1lm(y~-1+x,score="s",stand="o",maxiter = 10000, eps = 1e-4, eps.S = 1e-6)
  beta <- coef(mod1)
  
  plot(lambda,bic,type="l")
  abline(v=lambda_op)
  
  return(list(beta=beta,lambda=lambda, bic=bic, bic_min=bic_min, min_ind=min_ind, lambda_op = lambda_op))  
}

