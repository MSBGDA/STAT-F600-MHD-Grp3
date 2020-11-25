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

A_matrix_Normal <- function(p, sigma){
  Id_matrix <- 'diag<-'(matrix(0, p, p), 1)
  A <- matrix(0 , nrow = p, ncol = p)
  N <- 1000
  for (i in 1:N){
    eps <- chol(sigma) %*% rnorm(p, 0, 1)
    eps_normed <- norm(eps, type= '2')
    a <- (1/eps_normed) * (Id_matrix - (1/eps_normed**2) * eps %*% t(eps))
    A <- A + a
  }
  A <- A/N
  return(A)
}

WhereIsTheBug <- function(matrix){
  print('Where are the NA values ?')
  N <- nrow(matrix)
  M <- ncol(matrix)
  for (i in 1:N){
    for (j in 1:M)
      if ( !(is.finite(matrix[i,j]))){
        a <- sprintf('row %d and col %d', i,j)
        print(a)
      }
  }
}

# Use P-value on each row of a matrix to return the power for each case. 
# Each row of Matrix_input corresponds to a case, for example n=20 and p=1000
# Null_Hypo is 0 as in the article
# alpha is the level of significance

####
# PV Power_test       /!\ alpha = 0.05
####
Power_test <- function (Matrix_input, Null_Hypo, alpha){
  power_vector <- c()
  N_row <- nrow(Matrix_input)
  N_col <- ncol(Matrix_input)
  for (i in 1:N_row){
    std <- sd(Matrix_input[i,], na.rm = T)
    power <- 0
    #print(std)
    for (j in 1:N_col){
      Tn <- Matrix_input[i,j]
      #print(Tn)
      p <- (1-pnorm(Tn, Null_Hypo, std/sqrt(N_col))) * 2
      #print(p)
      if (p < alpha){
        power <- power + 1
      }
      #print('ok')
    }
    power <- power / N_col
    power_vector <- append(power_vector, power)
  }
  return(power_vector)
  }


# Part Example 1 ----

#this function provides histo. of the power for 4 cases
#Ex: HistoPower_Ex1(2,3)
HistoPower_Ex1_Ex2 <- function(WhichMu, WhichSigma){
  matrix_Tn1 <- matrix(0, nrow = 4, ncol = 500)
  matrix_Tn2 <- matrix(0, nrow = 4, ncol = 500)
  k <- 1
  # (n in c(10,12))  (n in c(20,50))
  for (n in c(10,25)) {
    # (p in c(21,30))   (p in c(1002,2001))
    for (p in c(102,402)){
      dat1 <- c()
      dat2 <- c()
      mu <- mu(WhichMu, p)
      Sigma <- Sigma(WhichSigma, p)
      print(paste0("n=",n," ", "p=",p))
      set.seed(17)              #Set seed for reproducibility
      for (i in 1:500){
        Data_matrix1 <- matrix(0, nrow = n, ncol = p)
        Data_matrix2 <- matrix(0, nrow = n, ncol = p)
        if (i %% 50 == 0){
          print(i)
        }
        for (i in 1:n){
          #random vector generated with cholesky decomposition
          Data_matrix1[i,] <- t(mu + t(chol(Sigma)) %*% rnorm(p, 0, 1))
          Data_matrix2[i,] <- t(mu + t(chol(Sigma)) %*% rt(p, 3))
        }
        for (i in 1:n){
          #matrix with the vector Zi. Each row is a vector of p components
          Data_matrix1[i,] <- Data_matrix1[i,] / norm(Data_matrix1[i,], type= '2')
          Data_matrix2[i,] <- Data_matrix2[i,] / norm(Data_matrix2[i,], type= '2')
        }
        Tn1 <- 0
        Tn2 <- 0
        for (i in 2:n){
          for (j in 1:(i-1)){
            Tn1 <- Tn1 + Data_matrix1[i,] %*% Data_matrix1[j,]
            Tn2 <- Tn2 + Data_matrix2[i,] %*% Data_matrix2[j,]
          }
        }
        Z_star1 <- 0
        Z_star2 <- 0
        for (i in 1:n){
          Z_star1 <- Z_star1 + Data_matrix1[i,]
          Z_star2 <- Z_star2 + Data_matrix2[i,]
        }
        Z_star1 <- Z_star1 / (n - 2)
        Z_star2 <- Z_star2 / (n - 2)
        Matrix_Zj_ZjT1 <- matrix(0, nrow = p, ncol = p)
        Matrix_Zj_ZjT2 <- matrix(0, nrow = p, ncol = p)
        for (i in 1:n){
          Matrix_Zj_ZjT1 <- Matrix_Zj_ZjT1 + Data_matrix1[i,] %*% t(Data_matrix1[i,])
          Matrix_Zj_ZjT2 <- Matrix_Zj_ZjT2 + Data_matrix2[i,] %*% t(Data_matrix2[i,])
        }
        
        Tr_B1 <- -n/(n - 2)**2 + ((n-1)/(n*(n-2)**2))* tr(Matrix_Zj_ZjT1 %*% Matrix_Zj_ZjT1) +
          ((1-2*n)/(n*(n-1))) * (Z_star1 %*% Matrix_Zj_ZjT1 %*% Z_star1) + 
          (2/n) * norm(Z_star1, type = '2')**2 + (((n-2)**2)/(n*(n-1))) * norm(Z_star1, type = '2')**4
        
        Tr_B2 <- -n/(n - 2)**2 + ((n-1)/(n*(n-2)**2))* tr(Matrix_Zj_ZjT2 %*% Matrix_Zj_ZjT2) +
          ((1-2*n)/(n*(n-1))) * (Z_star2 %*% Matrix_Zj_ZjT2 %*% Z_star2) + 
          (2/n) * norm(Z_star2, type = '2')**2 + (((n-2)**2)/(n*(n-1))) * norm(Z_star2, type = '2')**4
        
        Tn_normalized1 <- Tn1/((0.5*n*(n-1)*Tr_B1)**0.5)
        Tn_normalized2 <- Tn2/((0.5*n*(n-1)*Tr_B2)**0.5)
        
        dat1 <- append(dat1, Tn_normalized1)
        dat2 <- append(dat2, Tn_normalized2)
        
      }
      matrix_Tn1[k,] <- dat1
      matrix_Tn2[k,] <- dat2
      
      print(paste0("k=",k))
      k <- k + 1
    }
  }
  library(ggplot2)
  title <- sprintf("n data generated from N_p(mu%d, Sigma%d)", WhichMu, WhichSigma)
  x <- rep(c("n=10 \np=100", "n=10 \np=400", "n=25 \np=100", "n=25 \np=400"),2)
  Ex1 <- Power_test(matrix_Tn1, 0, 0.05)
  Ex2 <- Power_test(matrix_Tn2, 0, 0.05)
  y <-c(Ex1, Ex2)
  type <-c(rep("Normal distri.", 4), rep("t distribution", 4))
  mydata <-data.frame(x, y)
  p <-ggplot(mydata, aes(x, y))
  p +geom_bar(stat = "identity", aes(fill = type), position = "dodge", alpha=0.7) +
    ylab("Power") + xlab("Cases") + 
    ggtitle(title) + geom_text(aes(label=y), vjust=1.5)
}

HistoPower2_Ex1_Ex2 <- function(WhichMu, WhichSigma){
  matrix_Tn1 <- matrix(0, nrow = 4, ncol = 500)
  matrix_Tn2 <- matrix(0, nrow = 4, ncol = 500)
  k <- 1
  #c(5)
  for (n in c(30)) {
    #c(30,90,210,510) c(6,9,12,18)
    for (p in c(30,90,210,510)){
      dat1 <- c()
      dat2 <- c()
      mu <- mu(WhichMu, p)
      Sigma <- Sigma(WhichSigma, p)
      print(paste0("n=",n," ", "p=",p))
      set.seed(17)              #Set seed for reproducibility
      for (i in 1:500){
        Data_matrix1 <- matrix(0, nrow = n, ncol = p)
        Data_matrix2 <- matrix(0, nrow = n, ncol = p)
        if (i %% 50 == 0){
          print(i)
        }
        for (i in 1:n){
          #random vector generated with cholesky decomposition
          Data_matrix1[i,] <- t(mu + t(chol(Sigma)) %*% rnorm(p, 0, 1))
          Data_matrix2[i,] <- t(mu + t(chol(Sigma)) %*% rt(p, 3))
        }
        for (i in 1:n){
          #matrix with the vector Zi. Each row is a vector of p components
          Data_matrix1[i,] <- Data_matrix1[i,] / norm(Data_matrix1[i,], type= '2')
          Data_matrix2[i,] <- Data_matrix2[i,] / norm(Data_matrix2[i,], type= '2')
        }
        Tn1 <- 0
        Tn2 <- 0
        for (i in 2:n){
          for (j in 1:(i-1)){
            Tn1 <- Tn1 + Data_matrix1[i,] %*% Data_matrix1[j,]
            Tn2 <- Tn2 + Data_matrix2[i,] %*% Data_matrix2[j,]
          }
        }
        Z_star1 <- 0
        Z_star2 <- 0
        for (i in 1:n){
          Z_star1 <- Z_star1 + Data_matrix1[i,]
          Z_star2 <- Z_star2 + Data_matrix2[i,]
        }
        Z_star1 <- Z_star1 / (n - 2)
        Z_star2 <- Z_star2 / (n - 2)
        Matrix_Zj_ZjT1 <- matrix(0, nrow = p, ncol = p)
        Matrix_Zj_ZjT2 <- matrix(0, nrow = p, ncol = p)
        for (i in 1:n){
          Matrix_Zj_ZjT1 <- Matrix_Zj_ZjT1 + Data_matrix1[i,] %*% t(Data_matrix1[i,])
          Matrix_Zj_ZjT2 <- Matrix_Zj_ZjT2 + Data_matrix2[i,] %*% t(Data_matrix2[i,])
        }
        
        Tr_B1 <- -n/(n - 2)**2 + ((n-1)/(n*(n-2)**2))* tr(Matrix_Zj_ZjT1 %*% Matrix_Zj_ZjT1) +
          ((1-2*n)/(n*(n-1))) * (Z_star1 %*% Matrix_Zj_ZjT1 %*% Z_star1) + 
          (2/n) * norm(Z_star1, type = '2')**2 + (((n-2)**2)/(n*(n-1))) * norm(Z_star1, type = '2')**4
        
        Tr_B2 <- -n/(n - 2)**2 + ((n-1)/(n*(n-2)**2))* tr(Matrix_Zj_ZjT2 %*% Matrix_Zj_ZjT2) +
          ((1-2*n)/(n*(n-1))) * (Z_star2 %*% Matrix_Zj_ZjT2 %*% Z_star2) + 
          (2/n) * norm(Z_star2, type = '2')**2 + (((n-2)**2)/(n*(n-1))) * norm(Z_star2, type = '2')**4
        
        Tn_normalized1 <- Tn1/((0.5*n*(n-1)*Tr_B1)**0.5)
        Tn_normalized2 <- Tn2/((0.5*n*(n-1)*Tr_B2)**0.5)

        dat1 <- append(dat1, Tn_normalized1)
        dat2 <- append(dat2, Tn_normalized2)
        
      }
      matrix_Tn1[k,] <- dat1
      matrix_Tn2[k,] <- dat2

      print(paste0("k=",k))
      k <- k + 1
    }
  }
  library(ggplot2)
  title <- sprintf("n=30 data generated from N_p(mu%d, Sigma%d)", 2, 3)
  x <- rep(c("p1=30", "p2=90", "p3=210", "p4=510"),2)
  Ex1 <- Power_test(matrix_Tn1, 0, 0.05)
  Ex2 <- Power_test(matrix_Tn2, 0, 0.05)
  y <-c(Ex1, Ex2)
  type <-c(rep("Normal distri.", 4), rep("t distribution", 4))
  mydata <-data.frame(x, y)
  p <-ggplot(mydata, aes(x, y))
  p +geom_bar(stat = "identity", aes(fill = type), position = "dodge", alpha=0.7) +
    ylab("Power") + xlab("Cases") + 
    ggtitle(title) + geom_text(aes(label=y), vjust=1.5)
}

# This one don't work because of NA values from Ex3
HistoPower2_Ex1_Ex2_Ex3 <- function(WhichMu, WhichSigma){
  matrix_Tn1 <- matrix(0, nrow = 4, ncol = 500)
  matrix_Tn2 <- matrix(0, nrow = 4, ncol = 500)
  matrix_Tn3 <- matrix(0, nrow = 4, ncol = 500)
  k <- 1
  #c(30)
  for (n in c(5)) {
    #c(30,90,210,510)
    for (p in c(6,9,12,18)){
      dat1 <- c()
      dat2 <- c()
      dat3 <- c()
      mu <- mu(WhichMu, p)
      Sigma <- Sigma(WhichSigma, p)
      print(paste0("n=",n," ", "p=",p))
      set.seed(17)              #Set seed for reproducibility 
      for (i in 1:500){
        Data_matrix1 <- matrix(0, nrow = n, ncol = p)
        Data_matrix2 <- matrix(0, nrow = n, ncol = p)
        Data_matrix3 <- matrix(0, nrow = n, ncol = p)
        if (i %% 50 == 0){
          print(i)
        }
        for (i in 1:n){
          #random vector generated with cholesky decomposition
          Data_matrix1[i,] <- t(mu + t(chol(Sigma)) %*% rnorm(p, 0, 1))
          Data_matrix2[i,] <- t(mu + t(chol(Sigma)) %*% rt(p, 3))
          Data_matrix3[i,] <- 0.9*t(mu + t(chol(Sigma)) %*% rnorm(p, 0, 1)) +
            0.1*t(mu + t(chol(9*Sigma)) %*% rnorm(p, 0, 1))
        }
        for (i in 1:n){
          #matrix with the vector Zi. Each row is a vector of p components
          Data_matrix1[i,] <- Data_matrix1[i,] / norm(Data_matrix1[i,], type= '2')
          Data_matrix2[i,] <- Data_matrix2[i,] / norm(Data_matrix2[i,], type= '2')
          Data_matrix3[i,] <- Data_matrix3[i,] / norm(Data_matrix3[i,], type= '2')
        }
        Tn1 <- 0
        Tn2 <- 0
        Tn3 <- 0
        for (i in 2:n){
          for (j in 1:(i-1)){
            Tn1 <- Tn1 + Data_matrix1[i,] %*% Data_matrix1[j,]
            Tn2 <- Tn2 + Data_matrix2[i,] %*% Data_matrix2[j,]
            Tn3 <- Tn3 + Data_matrix3[i,] %*% Data_matrix3[j,]
          }
        }
        Z_star1 <- 0
        Z_star2 <- 0
        Z_star3 <- 0
        for (i in 1:n){
          Z_star1 <- Z_star1 + Data_matrix1[i,]
          Z_star2 <- Z_star2 + Data_matrix2[i,]
          Z_star3 <- Z_star3 + Data_matrix3[i,]
        }
        Z_star1 <- Z_star1 / (n - 2)
        Z_star2 <- Z_star2 / (n - 2)
        Z_star3 <- Z_star3 / (n - 2)
        Matrix_Zj_ZjT1 <- matrix(0, nrow = p, ncol = p)
        Matrix_Zj_ZjT2 <- matrix(0, nrow = p, ncol = p)
        Matrix_Zj_ZjT3 <- matrix(0, nrow = p, ncol = p)
        
        for (i in 1:n){
          Matrix_Zj_ZjT1 <- Matrix_Zj_ZjT1 + Data_matrix1[i,] %*% t(Data_matrix1[i,])
          Matrix_Zj_ZjT2 <- Matrix_Zj_ZjT2 + Data_matrix2[i,] %*% t(Data_matrix2[i,])
          Matrix_Zj_ZjT3 <- Matrix_Zj_ZjT3 + Data_matrix3[i,] %*% t(Data_matrix3[i,])
        }
        
        Tr_B1 <- -n/(n - 2)**2 + ((n-1)/(n*(n-2)**2))* tr(Matrix_Zj_ZjT1 %*% Matrix_Zj_ZjT1) +
          ((1-2*n)/(n*(n-1))) * (Z_star1 %*% Matrix_Zj_ZjT1 %*% Z_star1) + 
          (2/n) * norm(Z_star1, type = '2')**2 + (((n-2)**2)/(n*(n-1))) * norm(Z_star1, type = '2')**4
        
        Tr_B2 <- -n/(n - 2)**2 + ((n-1)/(n*(n-2)**2))* tr(Matrix_Zj_ZjT2 %*% Matrix_Zj_ZjT2) +
          ((1-2*n)/(n*(n-1))) * (Z_star2 %*% Matrix_Zj_ZjT2 %*% Z_star2) + 
          (2/n) * norm(Z_star2, type = '2')**2 + (((n-2)**2)/(n*(n-1))) * norm(Z_star2, type = '2')**4
        
        Tr_B3 <- -n/(n - 2)**2 + ((n-1)/(n*(n-2)**2))* tr(Matrix_Zj_ZjT3 %*% Matrix_Zj_ZjT3) +
          ((1-2*n)/(n*(n-1))) * (Z_star1 %*% Matrix_Zj_ZjT3 %*% Z_star3) + 
          (2/n) * norm(Z_star3, type = '2')**2 + (((n-2)**2)/(n*(n-1))) * norm(Z_star3, type = '2')**4
        
        Tn_normalized1 <- Tn1/((0.5*n*(n-1)*Tr_B1)**0.5)
        Tn_normalized2 <- Tn2/((0.5*n*(n-1)*Tr_B2)**0.5)
        Tn_normalized3 <- Tn3/((0.5*n*(n-1)*Tr_B3)**0.5)
        dat1 <- append(dat1, Tn_normalized1)
        dat2 <- append(dat2, Tn_normalized2)
        dat3 <- append(dat3, Tn_normalized3)
      }
      matrix_Tn1[k,] <- dat1
      matrix_Tn2[k,] <- dat2
      matrix_Tn3[k,] <- dat3
      print(paste0("k=",k))
      k <- k + 1
    }
  }
  library(ggplot2)
  title <- sprintf("n data generated from N_p(mu%d, Sigma%d)", 2, 3)
  x <- rep(c("p1=30", "p2=90", "p3=210", "p4=510"),3)
  print('attention')
  Ex1 <- Power_test(matrix_Tn1, 0, 0.05)
  Ex2 <- Power_test(matrix_Tn2, 0, 0.05)
  WhereIsTheBug(matrix_Tn3)
  Ex3 <- Power_test(matrix_Tn3, 0, 0.05)
  y <-c(Ex1, Ex2, Ex3)
  type <-c(rep("type 1", 4), rep("type 2", 4), rep("type 3", 4))
  mydata <-data.frame(x, y)
  p <-ggplot(mydata, aes(x, y))
  p +geom_bar(stat = "identity", aes(fill = type), position = "dodge", alpha=0.7) +
    ylab("Power") + xlab("Cases") + 
    ggtitle(title)
}

#plot of the theoretical and empirical error for case 1
#power vs mu. So no need to specify it when you call the function
EmpiricalPowerPlot_Ex1 <- function(n, p, WhichSigma){
  #x <- seq(0,0.3,by=0.02)
  x <- seq(0.3,0.6,by=0.05)
  z_score <- c()
  eta <- c()
  matrix_Tn <- matrix(0, nrow = length(x), ncol = 500)
  k <- 1
  Sigma <- Sigma(WhichSigma, p)
  print('Computation of the matrix A takes time...  ¯|_(*_*)_|¯')
  A <- A_matrix_Normal(p, Sigma)
  print("OK done !")
  for (q in x){
      dat <- c()
      Tr_B_vect <- c()
      mu <- rep(q, p)
      print(q)
      set.seed(17)              #Set seed for reproducibility 
      for (i in 1:500){
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
        Tr_B_vect <- append(Tr_B_vect, Tr_B)
        dat <- append(dat, Tn_normalized)
        Mean_Tr_B <- mean(Tr_B_vect)
      }
      z <- -1.95 + sqrt((0.5 * n * (n-1))/ Mean_Tr_B) * (t(mu) %*% (A %*% A) %*% mu)
      eta <- append(eta, ((t(mu) %*% (A %*% A) %*% mu))/sqrt(Mean_Tr_B))
      matrix_Tn[k,] <- dat
      z_score <- append(z_score, z)
      k <- k + 1
    }
  library(ggplot2)
  
  dataPower <- data.frame(
    #x <- seq(0,0.30,by=0.02),
    x <- eta,
    power_theo <- pnorm(z_score),
    power_empi <- Power_test(matrix_Tn, 0, 0.05)
  )
  title <- sprintf("n=%d data generated from N_%d(mu, Sigma%d)", n, p, WhichSigma)
  
  ggplot()+
    geom_line(data=dataPower,aes(y=power_theo,x= x,colour="blue"),size=1 )+
    geom_line(data=dataPower,aes(y=power_empi,x= x,colour="red"),size=1) +
    scale_color_discrete(name = "power", labels = c("Theoretical", "Empirical"))+
    #labs(x = expression(paste(mu)))+
    labs(x = expression(paste(eta)))+
    labs(y = expression(paste(beta)))+
    ggtitle(title)
}

Test_Ex1 <- function(n, p, q, wichsigma, alpha){
  Sigma <- Sigma(WhichSigma, p)
  mu <- rep(q, p)
  set.seed(17)
  dat <- c()
  for (i in 1:500){
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
  N_col <- length(dat)
  std <- sd(dat)
  power <- 0
  for (j in dat){
      Tn <- j
      p <- (1-pnorm(Tn, 0, std/sqrt(N_col))) * 2
      if (p < alpha){
        power <- power + 1
      }
    }
  power <- power / N_col
  print('power is')
  return(power)
}

Test_Ex2 <- function(n, p, q, wichsigma, alpha){
  Sigma <- Sigma(WhichSigma, p)
  mu <- rep(q, p)
  set.seed(17)
  dat <- c()
  for (i in 1:500){
    if (i %% 50 == 0){
      print(i)
    }
    Data_matrix <- matrix(0, nrow = n, ncol = p)
    for (i in 1:n){
      #random vector generated with cholesky decomposition
      Data_matrix[i,] <- t(mu + t(chol(Sigma)) %*% rt(p, 3))
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
  N_col <- length(dat)
  std <- sd(dat)
  power <- 0
  for (j in dat){
    Tn <- j
    p <- (1-pnorm(Tn, 0, std/sqrt(N_col))) * 2
    if (p < alpha){
      power <- power + 1
    }
  }
  power <- power / N_col
  print('power is')
  return(power)
}
#
#PART testing stuff, don't pay attention ----

# DIRECT Power_test   /!\ alpha = 0.975
Power_test <- function (Matrix_input, Null_Hypo, alpha){
  power_vector <- c()
  N_row <- nrow(Matrix_input)
  N_col <- ncol(Matrix_input)
  for (i in 1:N_row){
    std <- sd(Matrix_input[i,])
    power <- 0
    for (j in 1:N_col){
      Tn <- mean(Matrix_input[i,])
      error <- qnorm(alpha)*(std/sqrt(N_col))
      left <- Null_Hypo - error
      right <- Null_Hypo + error
      Zleft <- (left - Tn)/(std/sqrt(N_col))
      Zright <-(right - Tn)/(std/sqrt(N_col))
      p <- pnorm(Zright)-pnorm(Zleft)
      power <- (1-p) +  power
    }
    power <- power / N_col
    power_vector <- append(power_vector, power)
  }
  return(power_vector)
}

HistoPower_Ex1 <- function(WhichMu, WhichSigma){
  matrix_Tn <- matrix(0, nrow = 4, ncol = 500)
  k <- 1
  # (n in c(10,12))  (n in c(20,50))
  for (n in c(10,12)) {
    # (p in c(21,30))   (p in c(1002,2001))
    for (p in c(21,30)){
      dat <- c()
      mu <- mu(WhichMu, p)
      Sigma <- Sigma(WhichSigma, p)
      print(paste0("n=",n," ", "p=",p))
      set.seed(17)              #Set seed for reproducibility 
      for (i in 1:500){
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
  library(ggplot2)
  cases <- c("n=10,p=100", "n=10,p=400", "n=25,p=100", "n=25,p=400")
  power <- Power_test(matrix_Tn, 0, 0.05)
  print("histo shall be display:")
  print(power)
  d <- tibble::tibble(cases, power)
  p <- ggplot(d, aes(x = cases, y = power)) +
    geom_col(fill = 'blue', alpha = 0.3) +
    geom_text(aes(label=power), vjust=1.5)
  title <- sprintf("n data generated from N_p(mu%d, Sigma%d)", WhichMu, WhichSigma)
  p + ggtitle(title) +
    xlab("Cases") + ylab("Power")
}

set.seed(1234)
wdata = data.frame(
  sex = factor(rep(c("F", "M"), each=200)),
  weight = c(rnorm(200, 55), rnorm(200, 58))
)
head(wdata, 4)
library(ggplot2)
theme_set(
  theme_minimal() +
    theme(legend.position = "right")
)
# Une autre option pour geom = "point"
ggplot(wdata, aes(x = weight)) +
  stat_ecdf(aes(color = sex,linetype = sex), 
            geom = "step", size = 1.5) +
  scale_color_manual(values = c("#00AFBB", "#E7B800"))+ 
  labs(y = "f(weight)")


data <- data.frame(
  cases <- c("n=10,p=100", "n=10,p=400", "n=25,p=100", "n=25,p=400"),
  power <- c(5,9,3,6),
  truc <- c(6,10,5,8)
)

library(ggplot2)
title <- sprintf("n data generated from N_p(mu%d, Sigma%d)", 2, 3)
x <-rep(c("n=10,p=100", "n=10,p=400", "n=25,p=100", "n=25,p=400"), 3)
Ex1 <- c(5,9,3,6)
Ex2 <- c(6,10,5,8)
Ex3 <- c(4,7,2,8)
y <-c(Ex1, Ex2, Ex3)
type <-c(rep("type 1", 4), rep("type 2", 4), rep("type 3", 4))
mydata <-data.frame(x, y)
p <-ggplot(mydata, aes(x, y))
p +geom_bar(stat = "identity", aes(fill = type), position = "dodge2", alpha=0.7) +
  ylab("Power") + xlab("Cases") + 
  ggtitle(title) + geom_text(aes(label = c(Ex1,Ex2,Ex3)), position = position_dodge(0.9), vjust = 1)
