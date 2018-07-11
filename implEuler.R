# Solution of Differential Equations
# Laureen Lake, 580722

# Function: implicit Euler

# Input: 
#   - f               Function of the ODE y'= f(t,y)
#   - Jf              Jacobi-matrix of f 
#   - t_interval      [t0,tn] tuple with boundaries of interval considered
#   - y0              Initial Value Y(t0)
#   - h               stepsize
#   - TOL             error tolerance (for Newton method)
#   - IterMax         maximum number of iterations (for Newton method)

# Output:
#   - t         column vector with considered timepoints
#   - y         matrix that contains the approximated solution of y(t(i)) in row i 
#   - num_iter  number of Newton iterations for each point in time

implEuler <- function(f, Jf, t_interval, y0, h=NA, TOL, IterMax) {
  
  # get boundaries of interval
  t0 <- t_interval[1]
  tn <- t_interval[2]
  
  # specify stepsize
  h <- ifelse(is.na(h), (tn-t0)/1000, h)
  # specify number of timesteps
  n <- ceiling((tn-t0)/h)
  
  # dimension of y
  m <- length(y0)
  
  # set initial values
  t        <- matrix(seq(t0,t0 + (n - 1) * h, by = h) , nrow = n, ncol = 1)
  y        <- matrix(rep(0,n*m), nrow = n, ncol = m)
  num_iter <- matrix(rep(0, n), nrow = n, ncol = 1)
  I        <- diag(m)
  
  y[1, ] <- y0
  
  # timeloop of implicit Euler
  for (i in 2:n) {
    y0     <- y[i-1, ]
    y[i, ] <- y0
    error  <- 2 * TOL
    
    # Newton method  
    while (num_iter[i] <= IterMax && error >= TOL) {
      # evaluation of the function and the Jacobi-matrix
      funy   <- y[i, ] - y0 - h * f(t[i], y[i, ])
      Jfuny  <- I - h * Jf(t[i], y[i, ])
      deltay <- optR(Jfuny, -funy, method= "gauss")$beta
      y[i, ] <- y[i, ] + t(deltay)
      error  <- norm(as.matrix(funy))
      
      if (num_iter[i] == IterMax && error >= TOL) {
        cat("Maximum number of iterations was reached. The residual equals: f(t,y) = ", funy)
        break
      }
      num_iter[i] <- num_iter[i] + 1
    }
      
  }
  
  return(list(t, y, num_iter))
}
