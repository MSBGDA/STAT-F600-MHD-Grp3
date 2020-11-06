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

n <- 3                                 #size vector
p <- 6                                #dimension vector (multiple of 3)

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

Data_generator_Ex1 <- function(n, p, mu, Sigma, WhichMu, WhichSig){
  Data_matrix <- matrix(0, nrow = n, ncol = p)
  for (i in 1:n){
    #set.seed(i)          #Set seed for reproducibility
    Data_matrix[i,] <- rnorm(p, mu, Sigma)
  }
  name <- sprintf("Data_Ex1_n=%d_p=%d_mu%d_Sigma%d.txt", n, p, WhichMu, WhichSig)
  file.create(name)
  write.table(Data_matrix , name, append = TRUE, sep = " ", dec = ".",
              row.names = FALSE, col.names = FALSE)
}

