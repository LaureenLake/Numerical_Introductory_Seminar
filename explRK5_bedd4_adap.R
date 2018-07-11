# Solution of Differential Equations
# Laureen Lake, 580722

# Function: Runge-Kutta method of order 5 with 
#           embedded method of order 4 with adaptive step-size

# Input: 
#   - f               Function of the ODE y'= f(t,y)
#   - t_interval      [t0,tn] tuple with boundaries of interval considered
#   - y0              Initial Value Y(t0)
#   - h0              initial stepsize
#   - TOL             tolerance for local error

# Output:
#   - t         column vector with considered timepoints
#   - y         matrix that contains the approximated solution of y(t(i)) in row i 

install.packages("deSolve")
library(deSolve)

explRK5_bedd4_adap <- function(f, t_interval, y0, h0=NA) {
  
  # get boundaries of interval
  t0 <- t_interval[1]
  tn <- t_interval[2]
  
  # dimension of y
  m <- length(y0)
  
  # define vectors and set initial values
  t   <- NULL
  h   <- NULL
  y   <- NULL
  yd  <- NULL
  err <- NULL
  
  t   <- rbind(t, t0)
  h   <- rbind(h, h0)
  y   <- rbind(y, y0)
  yd  <- rbind(yd, y0)
  err <- rbind(err, rep(0, m))
  
  # k1(new step) = k7(old step)
  k7 <- f(t[1], y[1, ])
  
  # timeloop of classical explicit Runge-Kutta method
  for (i in 1:(n-1)) {
    t <- rbind(t, t[i] + h[i])
    
    k1 <- k7
    k2 <- f(t[i] + h[i]/5, y[i, ] + h[i] * k1/5)
    k3 <- f(t[i] + h[i] * 3/10, y[i, ] + h[i] * (3/40 * k1 + 9/40 * k2))
    k4 <- f(t[i] + h[i] * 4/5, y[i, ] + h[i] * (44/45 * k1 - 56/15 * k2 + 32/9 * k3))
    k5 <- f(t[i] + h[i]* 8/9, y[i, ] + h[i] * (19372/6561 * k1 - 25360/2187 * k2 + 64448/6561 * k3 - 212/729 * k4))
    k6 <- f(t[i+1], y[i, ] + h[i] * (9017/3168 * k1 - 355/33 * k2 + 46732/5247 * k3 + 49/176 * k4 - 5103/18656 * k5))
    k7 <- f(t[i+1], y[i, ] + h[i] * (35/384 * k1 + 500/1113 * k3 + 125/192 * k4 - 2187/6784 * k5 + 11/84 * k6))
    
    y  <- rbind(y, y[i, ] + h[i] * (35/384 * k1 + 500/1113 * k3 + 125/192 * k4 - 2187/6784 * k5 + 11/84 * k6))
    yd <- rbind(y, y[i, ] + h[i] *(5179/57600 * k1 + 7571/16695 * k3 + 393/640 * k4 - 92097/339200 * k5 + 187/2100 * k6 + 1/40 * k7))
    
    err     <- rbind(err, abs(y[i+1, ] - yd[i+1, ]))
    loc_err <- norm(err[i+1, ])
    h       <- rbind(h, ifelse(loc_err >= TOL, 0.9 * h[i] * (TOL/loc_err)^(1/5), 0.9 * h[i] * (TOL/loc_err)^(1/4)))
  }
  
  return(list(t, y, err, h))
}