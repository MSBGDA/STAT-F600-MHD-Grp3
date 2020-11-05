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

my_data <- read.delim(file.choose(), sep ="", header = FALSE, dec =".") #Choose the file 


