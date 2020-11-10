# PART initialisation ----

n <- 18                                 #size vector
p <- 20                                 #dimension vector
mu <- rep(0, p)                         #mean vector equal 0
Gamma <- 'diag<-'(matrix(0, p, p), 1)   # Gamma matrix p x p: Identity for testing


# PART random error epsilon_i in p dim ----

# random vector uniformly distributed on the unit hypersphere in R^p
Ui <- function(x){
  #Muller , Marsaglia ('Normalised Gaussians') method
  #http://extremelearning.com.au/how-to-generate-uniformly-random-points-on-n-spheres-and-n-balls/
  
  #set.seed(0)          #Set seed for reproducibility
  v <- rnorm(p, 0, 1)
  norm <- norm(v, type = '2')
  U_i <- v / norm
  return(U_i)
}

# non negative random variable independent of Ui:
Ri <- function(x){
  #set.seed(0)          #Set seed for reproducibility
  R_i <- runif(1, 0, 100) # 100 is arbitrary
  return(R_i)
}

# special case Ri: Ri^2 ~ Chi-square distribution, p degree of freedom:
Ri_Chi <- function(x){
  # https://statisticsglobe.com/chi-square-distribution-in-r-dchisq-pchisq-qchisq-rchisq
  
  #set.seed(0)          #Set seed for reproducibility
  R_i <- rchisq(1, df = p) ** 0.5
  return(R_i)
}

epsilon_i <- function(x){
  epsi <- Gamma %*% (Ri_Chi() * Ui()) #nb: %*% is for matrix product
  #epsi <- Gamma %*% (Ri() * Ui())
  return(epsi)
}


# PART generate data ----

file.create("C:/Users/maria/OneDrive/Documents/Unif 2021/mult. high dim. stat/DataTest.txt")
FilePath <- "C:/Users/maria/OneDrive/Documents/Unif 2021/mult. high dim. stat/DataTest.txt"    #To modify

for (i in 1:n){
  Xi <- mu + epsilon_i()
  write.table(t(Xi) , FilePath, append = TRUE, sep = " ", dec = ".",
              row.names = FALSE, col.names = FALSE) #nb: t() is for transpose
}

# PART mean vector test Tn ----

#my_data <- read.delim(file.choose(), sep ="", header = FALSE, dec =".") #Choose the file 

my_data <- matrix(scan(file = FilePath),ncol=p,byrow=TRUE) #import data from file

Z <- matrix(0, nrow = n, ncol = p) #create empty matrix to fill with the vectors Zi

for (i in 1:n){
  #matrix with the vector Zi. Each row is a vector of p components
  Z[i,] <- my_data[i,] / norm(my_data[i,], type= '2')
}

Test_n <- function(x){
  Tn <- 0               
  for (i in 2:n){
    for (j in 1:(i-1)){
      Tn <- Tn + Z[i,] %*% Z[j,]
  return(Tn)
    }
  }
}

# PART Monte Carlo Simulation ----

p <- 30                                #dimension vector (multiple of 3)

mu0 <- rep(0, p)
mu1 <- rep(0.25, p)
mu2 <- c(rep(0, p %/% 3), rep(0.25, p %/% 3), rep(-0.25, p %/% 3))

Sigma1 <- 'diag<-'(matrix(0.2, p, p), 1)

Sigma2 <- matrix(0 , nrow = p, ncol = p)
for (i in 1:p){
  for (j in 1:p){
    Sigma2[i, j] <- 0.8 ** abs(i - j)
  }
}

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
Sigma3 <- D %*% R %*% D

# Example: Data_generator_Ex1(4, 6, mu2, Sigma3, 2, 3) 2 and 3 because we use mu2 ans Sigma3 in this examle
Data_generator_Ex1 <- function(n, p, mu, Sigma, WhichMu, WhichSig){
  Data_matrix <- matrix(0, nrow = n, ncol = p)
  for (i in 1:n){
    #set.seed(i)          #Set seed for reproducibility
    Data_matrix[i,] <- t(mu + t(chol(Sigma)) %*% rnorm(p, 0, 1))
  }
  name <- sprintf("Data_Ex1_n=%d_p=%d_mu%d_Sigma%d.txt", n, p, WhichMu, WhichSig)
  file.create(name)
  write.table(Data_matrix , name, append = TRUE, sep = " ", dec = ".",
              row.names = FALSE, col.names = FALSE)
}

Data_generator_Ex2 <- function(n, p, mu, Sigma, WhichMu, WhichSig){
  Data_matrix <- matrix(0, nrow = n, ncol = p)
  for (i in 1:n){
    #set.seed(i)          #Set seed for reproducibility
    Data_matrix[i,] <- t(mu + t(chol(Sigma)) %*% rt(p, 3))
  }
  name <- sprintf("Data_Ex2_n=%d_p=%d_mu%d_Sigma%d.txt", n, p, WhichMu, WhichSig)
  file.create(name)
  write.table(Data_matrix , name, append = TRUE, sep = " ", dec = ".",
              row.names = FALSE, col.names = FALSE)
}

Data_generator_Ex3 <- function(n, p, mu, Sigma, WhichMu, WhichSig){
  Data_matrix <- matrix(0, nrow = n, ncol = p)
  for (i in 1:n){
    #set.seed(i)          #Set seed for reproducibility
    Data_matrix[i,] <- t(0.9 * (mu + t(chol(Sigma)) %*% rnorm(p, 0, 1)) + 0.1*(mu + t(chol(9 * Sigma)) %*% rnorm(p, 0, 1)))
  }
  name <- sprintf("Data_Ex3_n=%d_p=%d_mu%d_Sigma%d.txt", n, p, WhichMu, WhichSig)
  file.create(name)
  write.table(Data_matrix , name, append = TRUE, sep = " ", dec = ".",
              row.names = FALSE, col.names = FALSE)
}



# PART generate data and Tn test it ----

#Because create 1000 (check article 3.1) data files for each cases will take more times
#this section is a first try to generate 1000 data set and create an histogram 
#of the Test function. To verify the disribution

p <- 30                                #dimension vector (multiple of 3)

mu0 <- rep(0, p)
mu1 <- rep(0.25, p)
mu2 <- c(rep(0, p %/% 3), rep(0.25, p %/% 3), rep(-0.25, p %/% 3))

Sigma1 <- 'diag<-'(matrix(0.2, p, p), 1)

Sigma2 <- matrix(0 , nrow = p, ncol = p)
for (i in 1:p){
  for (j in 1:p){
    Sigma2[i, j] <- 0.8 ** abs(i - j)
  }
}

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
Sigma3 <- D %*% R %*% D

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

HistoTest_Ex1 <- function(n, p, mu , Sigma, WhichMu, WhichSig){
  dat <- c()
  for (i in 1:1000){
    if (i %% 50 == 0){
      print(i)
    }
    Data_matrix <- matrix(0, nrow = n, ncol = p)
    for (i in 1:n){
      #set.seed(i)          #Set seed for reproducibility
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
      ((1-2*n)/(n*(n-1))) * Z_star %*% Matrix_Zj_ZjT %*% Z_star + 
      (2/n) * norm(Z_star, type = '2')**2 + (((n-2)**2)/(n*(n-1))) * norm(Z_star, type = '2')**4
    Tn_normalized <- Tn/((0.5*n*(n-1)*Tr_B)**0.5)
    dat <- append(dat, Tn_normalized)
  }
  title <- sprintf("Case 1. n=%d p=%d mu%d Sigma%d", n, p, WhichMu, WhichSig)
  h <- hist(dat,
        breaks = 20,
        prob = FALSE
    )
  h$counts <- h$counts / sum(h$counts)
  plot(h,main= title, freq=TRUE,xlab="Tn", ylab="Relative Frequency", ylim=c(0,0.5))
  #m<-mean(dat)
  #std<-sqrt(var(dat))
  curve(dnorm(x, mean=0, sd=1), 
        col="darkblue", lwd=2, add=TRUE, yaxt="n")
}
