# /!\ Dont't forget to run the needed PART

# PART needed PART ----
# def. of all the stuff we need

# Trace of a matrix
tr <- function (m){
  total_sum <- 0
  if(is.matrix(m))
  {
    row_count <- nrow(m)
    col_count <- ncol(m)
    if(row_count == col_count)
    {
      total_sum <-sum(diag(m))
      total_sum
    }
    else
    {
      message ('Matrix is not square')
    }
  }
  else
  {
    message( 'Object is not a matrix')
    
  }
}

# function mu: return vector mu. Ex: mu(2, 4) return vector mu2 of dim 4
mu <- function (WhichMu, p){
  if (WhichMu == 0){
    mu <- rep(0, p)
  } else if (WhichMu == 1){
    mu <- rep(0.25, p)
  } else if (WhichMu == 2){
    mu <- c(rep(0, p %/% 3), rep(0.25, p %/% 3), rep(-0.25, p %/% 3))
  } else {
    print("ERROR: WhichMu must be 0, 1 or 2")
    next
  }
  return(mu)
}

# function Sigma: return matrix Sigma. Ex: mu(2, 4) return matrix Sigma2 of dim 4x4
Sigma <- function (WhichSigma, p){
  if (WhichSigma == 1){
    Sigma <- 'diag<-'(matrix(0.2, p, p), 1)
  } else if (WhichSigma == 2){
    Sigma <- matrix(0 , nrow = p, ncol = p)
    for (i in 1:p){
      for (j in 1:p){
        Sigma[i, j] <- 0.8 ** abs(i - j)
      }
    }
  } else if (WhichSigma == 3){
    D <- matrix(0 , nrow = p, ncol = p)
    for (i in 1:p){
      D[i, i] <- 2 + (p - i + 1) / p
    }
    R <- matrix(0 , nrow = p, ncol = p)
    for (i in 1:p){
      for (j in 1:p){
        R[i, j] <- (-1) ** (i + j) * (0.2) ** (abs(i - j) ** (0.1))
      }
    }
    Sigma <- D %*% R %*% D
  } else {
    print("ERROR: WhichSigma must be 1, 2 or 3")
    next
  }
  return(Sigma)
}

# Use P-value on each row of a matrix to return the power for each case. 
# Each row of Matrix_input corresponds to a case, for example n=20 and p=1000
# Null_Hypo is 0 as in the article
# alpha is the level of significance
Power_test <- function (Matrix_input, Null_Hypo, alpha){
  power_vector <- c()
  N_row <- nrow(Matrix_input)
  N_col <- ncol(Matrix_input)
  for (i in 1:N_row){
    std <- sd(Matrix_input[i,])
    power <- 0
    for (j in 1:N_col){
      Tn <- Matrix_input[i,j]
      p <- 2*(1-pnorm(Tn, Null_Hypo, std/sqrt(N_col)))
      if (p < alpha){
        power <- power + 1
      }
    }
    power <- power / N_col
    power_vector <- append(power_vector, power)
  }
  return(power_vector)
}

# Part Example 1 ----

#this function provides histo. of the power for 4 cases
#Ex: HistoPower_Ex1(2,3)
HistoPower_Ex1 <- function(WhichMu, WhichSigma){
  matrix_Tn <- matrix(0, nrow = 4, ncol = 1000)
  k <- 1
  # (n in c(10,12))
  for (n in c(20,50)) {
    # (p in c(21,30)
    for (p in c(1002,2001)){
      dat <- c()
      mu <- mu(WhichMu, p)
      Sigma <- Sigma(WhichSigma, p)
      print(paste0("n=",n," ", "p=",p))
      set.seed(17)              #Set seed for reproducibility 
      for (i in 1:1000){
        if (i %% 50 == 0){
          print(i)
        }
        Data_matrix <- matrix(0, nrow = n, ncol = p)
        for (i in 1:n){
          #random vector generated with cholesky decomposition
          Data_matrix[i,] <- t(mu + t(chol(Sigma)) %*% rnorm(p, 0, 1))
        }
        for (i in 1:n){
          #matrix with the vector Zi. Each row is a vector of p components
          Data_matrix[i,] <- Data_matrix[i,] / norm(Data_matrix[i,], type= '2')
        }
        Tn <- 0               
        for (i in 2:n){
          for (j in 1:(i-1)){
            Tn <- Tn + Data_matrix[i,] %*% Data_matrix[j,]
          }
        }  
        Z_star <- 0
        for (i in 1:n){
          Z_star <- Z_star + Data_matrix[i,]
        }
        Z_star <- Z_star / (n - 2)
        Matrix_Zj_ZjT <- matrix(0, nrow = p, ncol = p)
        for (i in 1:n){
          Matrix_Zj_ZjT <- Matrix_Zj_ZjT + Data_matrix[i,] %*% t(Data_matrix[i,])
        }
        Tr_B <- -n/(n - 2)**2 + ((n-1)/(n*(n-2)**2))* tr(Matrix_Zj_ZjT %*% Matrix_Zj_ZjT) +
          ((1-2*n)/(n*(n-1))) * (Z_star %*% Matrix_Zj_ZjT %*% Z_star) + 
          (2/n) * norm(Z_star, type = '2')**2 + (((n-2)**2)/(n*(n-1))) * norm(Z_star, type = '2')**4
          Tn_normalized <- Tn/((0.5*n*(n-1)*Tr_B)**0.5)
        dat <- append(dat, Tn_normalized)
      }
      matrix_Tn[k,] <- dat
      print(paste0("k=",k))
      k <- k + 1
    }
  }
  cases <- c("n=20,p=1000", "n=20,p=2000", "n=50,p=1000", "n=50,p=2000")
  power <- Power_test(matrix_Tn, 0, 0.05)
  print("histo shall be display:")
  print(power)
  d <- tibble::tibble(cases, power)
  p <- ggplot(d, aes(x = cases, y = power)) +
    geom_col(fill = 'blue', alpha = 0.3)
  title <- sprintf("Empirical power: n data generated from N_p(mu%d, Sigma%d)", WhichMu, WhichSigma)
  p + ggtitle(title) +
    xlab("Cases") + ylab("Power")
}