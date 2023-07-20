

################################################################################
# Function for simulating multivariate dataset
################################################################################

sim_mdata <-
  function(n = 100,
           p = 50,
           p0 = 10,
           q = 50,
           q0 = 10,
           c = 5,
           nrank = 3,
           s2n = 1,
           sigma = NULL,
           rho_X = 0.5,
           rho_E = 0,
           rho_C = 0) {
    CorrCS <- function(p, rho)
    {
      Sigma <- matrix(nrow = p, ncol = p, rho)
      diag(Sigma) <- 1
      Sigma
    }
    Sigma = CorrCS
    
    A1 <- matrix(ncol = nrank, nrow = q0, rnorm(q0 * nrank))
    A0 <- matrix(ncol = nrank, nrow = q - q0, 0)
    A <- rbind(A1, A0)
    B1 <- matrix(ncol = nrank, nrow = p0, rnorm(p0 * nrank))
    B0 <- matrix(ncol = nrank, nrow = p - p0, 0)
    B <- rbind(B1, B0)
    C <- B %*% t(A)
    
    Csigma <- Sigma(c, rho_C)
    CONF <- MASS::mvrnorm(n, rep(0, c), Csigma)
    Xsigma <- Sigma(p, rho_X)
    X <- MASS::mvrnorm(n, rep(0, p), Xsigma)
    A1_c <- matrix(ncol = nrank, nrow = c, rnorm(c * nrank))
    B1_c <- matrix(ncol = nrank, nrow = p, rnorm(p * nrank))
    C_c <- A1_c %*% t(B1_c)
    dim(C_c)
    X_C <- CONF %*% C_c + X
    X <- cbind(X[,-c(1:(0.2*p))],X_C[,c(1:(0.2*p))])
    UU <- MASS::mvrnorm(n, rep(0, q), Sigma(q, rho_E))
    
    #    ###Their definition is tr(CXC)/tr(E), which seems to be wrong
    #    sigma <- sqrt(sum(diag(t(C)%*%Sigma(p,rho_X)%*%C))/sum(diag(t(UU)%*%UU))/s2n)
    #    UU <- UU*sigma
    A1_cy <- matrix(ncol = nrank, nrow = c, rnorm(c * nrank))
    B1_cy <- matrix(ncol = nrank, nrow = q, rnorm(q * nrank))
    C_cy <- A1_cy %*% t(B1_cy)
    
    svdC_cy <- svd(C_cy)
    C3_cy <- svdC_cy$u[, nrank] %*% t(svdC_cy$v[, nrank]) * svdC_cy$d[nrank]
    
    svdC_xy <- svd(C)
    C3_xy <- svdC_xy$u[, nrank] %*% t(svdC_xy$v[, nrank]) * svdC_xy$d[nrank]
    Y3 <- X %*% C3_xy + CONF %*% C3_cy
    
    ##sigma <- sqrt(var(as.numeric(Y3))/var(as.numeric(UU))/s2n)
    ##the same
    if (is.null(sigma)) {
      sigma <- sqrt(sum(as.numeric(Y3) ^ 2) / sum(as.numeric(UU) ^ 2) / s2n)
    }
    UU <- UU * sigma
    
    Y <- matrix(nrow = n, ncol = q, NA)
    Y <- X %*% C3_xy + CONF %*% C3_cy + UU
    
    list(
      Y = Y,
      X = cbind(X,CONF),
      C = C,
      A = A,
      B = B,
      U = B,
      V = A,
      sigma = sigma,
      Xsigma = Xsigma,
      Csigma = Csigma
    )
    
  }


################################################################################
# Nesterov's accelerated gradient descent, for Graph-guided fused Lasso
################################################################################

# utils ------------------------------------------------------------------------

#' Standard soft-thresholding operator
#' @param v The vector to soft threshold
#' @param lambda The amount to soft-threshold by
#' @param covar Position of confounders in the X dataset.
#' @export
soft_threshold <- function(v, lambda, covar) {
  #----------changes
  res <- v 
  v[covar,] <- NA
  res[which(v > lambda)] <- v[which(v > lambda)] - lambda
  res[which(v < -lambda)] <- v[which(v < -lambda)] + lambda
  res[which(v > -lambda & v < lambda)] <- 0
  res
  #----------
}

#' Get the objective of the current estimates
#' @param X The data matrix.
#' @param B The J x K regression coefficient estimates.
#' @param C The L1 penalization matrix, returned by \code{get_C()}.
#' @param R The matrix of (thresholded) correlations between columns of Y
#' @return The error of the reconstruction according to the original (nonsmooth)
#' objective.
#' @export
objective <- function(X, B, Y, C, lambdaseq) {
  #(1 / 2) * sum( (Y - X %*% B) ^ 2) + lambda * sum(abs(B)) + sum(abs(B %*% t(C)))
  #----------changes
  (1 / 2) * sum( (Y - X %*% B) ^ 2) + sum(lambdaseq * abs(B)) + sum(abs(B %*% t(C)))
  #----------
}

# gradient-descent-funs --------------------------------------------------------

#' Get Optimal Alpha
#' S(C * W' / mu) where S projects into the interval [-1, 1].
#' @param C The L1 penalization matrix, returned by \code{get_C()}.
#' @param W The W matrix in Algorithm 1 of the reference
#' @param mu The smoothing parameter.
#' @return S(C * W' / mu)
#' @export
get_alpha_opt <- function(C, W, mu) {
  alpha <- C %*% t(W) / mu
  alpha[alpha > 1] <- 1
  alpha[alpha < -1] <- -1
  t(alpha)
}

#' Calculate the gradient for one step
#' @param Y The matrix of regression responses.
#' @param X The data matrix.
#' @param B The J x K regression coefficient estimates.
#' @param C The L1 penalization matrix, returned by \code{get_C()}.
#' @param mu The smoothing parameter.
#' @return The gradient descent direction for B, given the current estimates.
#' See equation (11) in the reference.
#' @references Smoothing Proximal Gradient Method for General Structured Sparse Regressoin
#' @export
get_grad_f <- function(X, Y, W, C, mu) {
  alpha <- get_alpha_opt(C, W, mu)
  t(X) %*% (X %*% W - Y) + alpha %*% C
}

#' Get next value of B in gradient descent
#' @param W The W matrix in Algorithm 1 of the reference
#' @param grad_f The gradient computed in Algorithm 1
#' @param L The lipshitz constant, which determines the step size
#' @param lambda The l1 regularization parameter.
#' @param covar Position of confounders in the X dataset.
#' @export
get_B_next <- function(W, grad_f, L, lambda, covar) {
  B_next <- W - (1 / L) * grad_f
  soft_threshold(B_next, lambda / L, covar)
}

#' Nesterov's Accelerated Gradient Descent
#' @param X The data matrix.
#' @param Y The matrix of regression responses.
#' @param H  The matrix H defined in the reference.
#' @param opts A list of gradient descent tuning parameters, \cr
#'   $iter_max: The maximum number of iterations.
#'   $mu: The smoothing parameter
#'   $lambda: The l1 regularization parameter
#'   $delta_conv: The convergence criterion for delta
#' @param covar Position of confounders in the X dataset.
#' @return A list containing the following elements \cr
#'   $B The final estimates of the multitask regression coefficients
#'   $obj The objective function across gradient descent iterations
#' @export
accgrad <- function(X, Y, C, opts, covar) {
  
  
  #----------changes
  if (!is.null(covar)) { 
    if (any(covar > ncol(X)) | any(covar < 1)) { 
      print("Error: Covariate index out of the dimensions of the input.") 
      break;
    }
  }
  penfactor <- rep(1,ncol(X))
  penfactor[covar] <- 0
  #penfactor=penfactor/sum(penfactor)*ncol(X)
  lambda_seq=opts$lambda*penfactor
  #----------
  
  
  # initialize results
  B <- opts$B0
  W <- B
  obj <- vector(length = opts$iter_max)
  
  theta <- 1
  if (opts$verbose) cat("\titer\t|\tobj\t|\t|B(t + 1) - B(t)| \n")
  for (iter in seq_len(opts$iter_max)) {
    # make a step
    grad_f <- get_grad_f(X, Y, W, C, opts$mu)
    B_next <- get_B_next(W, grad_f, opts$L, opts$lambda, covar)
    
    theta_next <- 2 / (iter + 1)
    W <- B_next + ((1 - theta) / theta) * (theta_next) * (B_next - B)
    
    # check convergence, and update counter
    delta <- sum(abs(B_next - B))
    B <- B_next
    obj[iter]  <- objective(X, B, Y, C, lambda_seq)
    if (iter %% 10 == 0 & opts$verbose) {
      cat(sprintf("%d \t | %f \t | %f \n", iter, obj[iter], delta))
    }
    if (delta < opts$delta_conv | iter > opts$iter_max) break
    theta <- theta_next
  }
  list(B = B, obj = obj[seq_len(iter)])
}


################################################################################
#  Cross-validation functions for Graph-guided fused Lasso.
################################################################################

## This can't just be done
## directly by caret since that package doesn't support multitask regression.
##
## author: Francisco Lima (https://github.com/monogenea), with revisions by
## sankaran.kris@gmail.com
## date: 12/26/2017

#' Compute the root mean squared error
rmse <- function(pred, y) {
  sqrt(mean( (pred - y) ^ 2 ))
}
#' k-fold Cross Validation for GFLasso
#'
#' @param Y The matrix of regression responses, scaled and centered as necessary.
#' @param X The data matrix, scaled and centered as necessary.
#' @param R The matrix of (thresholded) correlations between columns of Y.
#' @param additionalOpts Additional options to pass alongside lambda and gamma. See merge_proxgrad_opts().
#' @param k Number of folds.
#' @param times Number of repetitions. Total number of metric estimates = no. folds x no. times.
#' @param params The grid of lambda and gamma values to cross-validate.
#' @param nCores The number of CPU cores to be used.
#' @param seed Arbitrary number to ensure reproducibility. Defaults to 100.
#' @param err_fun A function that computes the metric (error/goodness-of-fit) between vectors of
#'   predicted and true responses. Defaults to rmse(pred, y) = sqrt(mean((pred - y) ^ 2)).
#' @param err_opt Specify whether do minimize ('min') or maximize ('max') `err_fun`.
#'   Default is 'min'.
#' @param covar Position of confounders in the X dataset.
#' @return cvMatrix A matrix of errors across a grid of lambda (row) and gamma
#'   (column) values.
#' @importFrom parallel mclapply detectCores
#' @importFrom pbapply pblapply
#' @importFrom caret createMultiFolds
#' @examples
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' u <- matrix(rnorm(10), 10, 1)
#' B <- u %*% t(u) + matrix(rnorm(10 * 10, 0, 0.1), 10, 10)
#' Y <- X %*% B + matrix(rnorm(100 * 10), 100, 10)
#' R <- ifelse(cor(Y) > .8, 1, 0)
#' system.time(testCV <- cv_gflasso(scale(X), scale(Y), R, nCores = 1))
#' system.time(testCV <- cv_gflasso(scale(X), scale(Y), R, nCores = 2))
#' cv_plot_gflasso(testCV)
#' @export
cv_gflasso <- function(X, Y, R, additionalOpts = list(), k = 5, times = 1,
                       params = seq(0, 1, by = 0.1), nCores = NULL, seed = 100,
                       err_fun = rmse, err_opt = 'min', covar) {
  
  additionalOpts <- merge_proxgrad_opts(additionalOpts, ncol(X), ncol(Y))
  if (is.null(nCores)) {
    nCores <- detectCores() - 1
  }
  
  cvFUN <- function(Y, X, R, opts, cvIndex) {
    sapply(seq_along(cvIndex), function(i){
      mod <- gflasso(Y = Y[-cvIndex[[i]],], X = X[-cvIndex[[i]],], R = R, opts = opts, covar=covar)
      pred <- X[cvIndex[[i]],] %*% mod$B
      err_fun(pred, Y[cvIndex[[i]], ])
    })
  }
  set.seed(seed)
  cvIndex <- caret::createMultiFolds(1:nrow(Y), k = k, times = times)
  cvArray <- array(NA, dim = c(rep(length(params), 2), k * times))
  dimnames(cvArray) <- list(params, params, names(cvIndex))
  
  grid <- expand.grid(lambda = params, gamma = params)
  allCV <- pbapply::pblapply(
    as.list(1:nrow(grid)),
    function(x) {
      if (additionalOpts$verbose && x %% 10 == 0) {
        cat(sprintf("CV grid %s/%s \n", x, nrow(grid)))
      }
      
      cvFUN(
        X = X, Y = Y, R = R,
        opts = list(lambda = grid[x, 1], gamma = grid[x, 2], additionalOpts),
        cvIndex = cvIndex
      )
    },
    cl = nCores
  )
  
  print(allCV[[1]])
  for(i in 1:nrow(grid)){
    cvArray[as.character(grid[i,1]),as.character(grid[i,2]),] <- allCV[[i]]
  }
  
  cvMean <- apply(cvArray, 1:2, mean)
  if(err_opt == 'min'){
    opt <- grid[which.min(cvMean), ]
  }else if(err_opt == 'max'){
    opt <- grid[which.max(cvMean), ]
  }
  list(
    "mean" = cvMean,
    "SE" = apply(cvArray, 1:2, sd) / sqrt(k * times),
    "optimal" = opt,
    "err_fun" = as.character(substitute(err_fun))
  )
}

#' Plot Results from Cross Validation
#'
#' @importFrom pheatmap pheatmap
#' @export
cv_plot_gflasso <- function(cv.gflasso){
  pheatmap(cv.gflasso$mean, cluster_rows = F, cluster_cols = F,
           main = paste("CV mean", cv.gflasso$err_fun, "\nOptimal pars:", "lambda =", cv.gflasso$optimal$lambda,
                        ",", "gamma =", cv.gflasso$optimal$gamma))
}

################################################################################
# Graph-guided fused Lasso
################################################################################

#' @title Merge default proximal gradient descent options
#' @param opts A potentially partially specified list of regularization and
#' gradient descent parameters. The currently supported options are,
#' $delta_conv How small does the change in B need to be to declare convergence?
#' $eps A tolerance used in calculating the degree of smoothing in the proximal
#' gradient objective.
#' $gamma The graph regularization parameter.
#' $iter_max What is the maximum number of iterations we should run?
#' $lambda The l1 regularization parameter.
#' $verbose Should the gradient descent print its progress?
#' @return A modified version of opts with defaults filled in.
#' @export
merge_proxgrad_opts <- function(opts, J, K) {
  default_opts <- list()
  default_opts$delta_conv <- 1e-2
  default_opts$eps <- 0.005
  default_opts$gamma <- 1
  default_opts$iter_max <- 1e3
  default_opts$lambda <- 1
  default_opts$verbose <- FALSE
  default_opts$B0 <- matrix(0, J, K)
  modifyList(default_opts, opts)
}

#' @title Graph-guided fused Lasso via Smoothed Proximal Gradient Descent
#' @param Y The matrix of regression responses, scaled and centered as necessary.
#' @param X The data matrix, scaled and centered as necessary.
#' @param R The matrix of (thresholded) correlations between columns of Y
#' @param opts A potentially partially specified list of regularization and
#' gradient descent parameters. See merge_proxgrad_opts().
#' @param covar Position of confounders in the X dataset.
#' @return A list containing the following quantities: \cr
#'   $B The estimated beta coefficient matrix. Its rows correspond columns of
#'    X, while its columns correspond to columns of Y. \cr
#'   $obj The graph fused lasso objective function, over iterations.
#'   $Lu The automatically calculated step size. \cr
#' reference.
#' @references Smoothing Proximal Gradient Method for General Structured Sparse Regression
#' @export
gflasso <- function(Y, X, R, opts = list(), covar) {
  
  # get opts
  opts <- merge_proxgrad_opts(opts, ncol(X), ncol(Y))
  
  # get L1 penalty matrix
  C <- opts$gamma * t(get_H(R))
  
  # calculate automatic step size
  D <- (1 / 2) * ncol(X) * (ncol(Y) + ncol(C) / 2)
  mu <- opts$eps / (2 * D)
  Lu <- get_Lu(X, C, opts$lambda, opts$gamma, mu)
  
  accgrad_opts <- list(lambda = opts$lambda, L = Lu, mu = mu,
                       iter_max = opts$iter_max, delta_conv = opts$delta_conv,
                       verbose = opts$verbose, B0 = opts$B0)
  optim_result <- accgrad(X, Y, C, accgrad_opts, covar)
  list(B = optim_result$B, obj = optim_result$obj, Lu = Lu)
}


#' @title GFLasso prediction
#'
#' @param model The model object from `gflasso`.
#' @param new.data The data from which to predict.
#'
#' @return A n x k matrix carrying all k predicted responses across all n samples
#' @export
predict_gflasso <- function(model, new.data){
  # Simple matrix multiplication
  return(new.data %*% model$B)
}


################################################################################
# Defining function for running MTL on different bootstrap subsamples
################################################################################
# This function is a wrapper to MTL main function, further implemmenting the estimation of P-values for the quantification of uncertainty in variable selection.
# for this demonstration, we will opt for running the standard MTL model in the whole dataset without the bootstrapping procedure (this is achieved through the argument "wholedataset=TRUE")
par_boostrapping_MTL <- function(x,y,covars,lam=0.1,wholedataset=TRUE) {
  
  if (!wholedataset) {
    indices <- sample(1:nrow(x), replace=TRUE)
  }else{ indices <- 1:nrow(x)}
  
  dataX <- x[indices,]
  dataY <- y[indices,]
  
  # Create a different object for confounders
  C <- dataX[,which(colnames(dataX) %in% covars)]
  
  # Remove confounders from the predictors
  dataX <- dataX[,-which(colnames(dataX) %in% covars)]
  
  xnamescat <- names(lapply(apply(dataX,2,table),length))
  whichxnamescat <- lapply(apply(dataX,2,table),length) < 10
  xnamescat <- xnamescat[whichxnamescat]
  
  # Residualization of Xs for confounders adjustment
  dataRES <- cbind(dataX,C)
  new_data_dataX <- dataX
  for (j in 1:ncol(dataX)) {
    f <- as.formula(paste(colnames(dataRES)[j]," ~ .",sep=""))
    if (j %in% which(colnames(dataX) %in% xnamescat)) {
    }else{
      resglm <- glm( f, family = gaussian, data=as.data.frame(dataRES[,c(j,(ncol(dataX)+1):ncol(dataRES))]))
      new_data_dataX[,j] <- resglm$residuals
    }
  }
  # Residualization of Ys for confounders adjustment
  dataRES <- cbind(dataY,C)
  new_data_dataY <- dataY
  for (j in 1:ncol(dataY)) {
    f <- as.formula(paste(colnames(dataRES)[j]," ~ .",sep=""))
    resglm <- glm( f, family = gaussian, data=as.data.frame(dataRES[,c(j,(ncol(dataY)+1):ncol(dataRES))]))
    new_data_dataY[,j] <- resglm$residuals
  }
  
  dataY <- new_data_dataY
  dataX_f <- do.call("list", replicate(ncol(dataY), new_data_dataX, simplify = FALSE))
  dataY_f <- list()                   # Create empty list
  for(k in 1:ncol(dataY)) {             # Using for-loop to add columns to list
    dataY_f[[k]] <- as.matrix(dataY[ , k])
  }
  names(dataY_f) <- colnames(dataY)
  
  return(MTL(dataX_f, dataY_f, type="Regression", Regularization="L21",
             Lam1=lam, Lam2=0, opts=list(init=0, tol=10^-6, maxIter=1500)))
  print(paste("model",as.character(i),"completed",sep=" "))
  
}
