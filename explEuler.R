# Solution of Differential Equations
# Laureen Lake, 580722

# Function: explicit Euler

# Input: 
#   - f               Function of the ODE y'= f(t,y)
#   - t_interval      [t0,tn] tuple with boundaries of interval considered
#   - y0              Initial Value Y(t0)
#   - h               stepsize

# Output:
#   - t   column vector with considered timepoints
#   - y   matrix that contains the approximated solution of y(t(i)) in row i 

explEuler <- function(f, t_interval, y0, h=NA) {

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
  t <- matrix(rep(0, n), nrow = n, ncol = 1)
  y <- matrix(rep(0,n*m), nrow = n, ncol = m)
  
  t[1]   <- t0
  y[1, ] <- y0

  # timeloop of explicit Euler
  for (i in 1:(n-1)) {
    y[i+1, ] <- y[i, ] + h * f(t[i], y[i, ])
    t[i+1]   <- t[i] + h
  }
  
  return(list(t,y))
}
