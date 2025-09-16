## Name: Isty Rysava
## Date: 4/6/2024
## Code: Function to simulate transmission chains for any distribution

chain_alpha <- function(gens = Inf, I = 1, R0seas, alpha, pop){
  Z <- list() 
  Z[[1]] <- I
  j <- 1 
  while(sum(Z[[j]]) > 0 && j <= gens) { 
    lambda <- R0seas * (sum(Z[[j]])^alpha)/pop # calculate FOI 
    Z[[j+1]] <- rbinom(1, pop-sum(Z[[j]]), lambda)  
    j <- j+1 
  } 
  out <- c(length=length(Z)-2, total=sum(unlist(Z))-1)
}

chain_a <- function(gens = Inf, I = 1, a, bseas, cseas, mumseas, PDRseas, pop, rm){
  Z <- list() 
  Z[[1]] <- I
  j <- 1 
  while(sum(Z[[j]]) > 0 && j <= gens) { 
    lambdaM <-  (a * bseas * exp(-mumseas/(PDRseas+ec)))/(rm) * sum(Z[[j]])
    M <-  rpois(1, lambdaM)
    lambdaH <- a * cseas * 1/(mumseas) * M * ((pop-sum(Z[[j]]))/pop)
    Z[[j+1]] <- rpois(1, lambdaH)  
    j <- j+1 
  } 
  out <- c(length=length(Z)-2, total=sum(unlist(Z))-1)
}
 
chain_a_alpha <- function(gens = Inf, I = 1, a, bseas, cseas, mumseas, PDRseas, pop, rm, alpha){
  Z <- list() 
  Z[[1]] <- I
  j <- 1 
  while(sum(Z[[j]]) > 0 && j <= gens) { 
    seasbet <-  (a^2 * cseas * bseas * cseas * exp(-mumseas/(PDRseas+ec)))/(rm*mumseas)
    lambdaH <-  seasbet * sum(Z[[j]])^alpha * (pop-sum(Z[[j]]))/pop
    Z[[j+1]] <- rpois(1, lambdaH)  
    j <- j+1 
  } 
  out <- c(length=length(Z)-2, total=sum(unlist(Z))-1)
}
 
chain_a_k <- function(gens = Inf, I = 1, a, bseas, cseas, mumseas, PDRseas, pop, rm, k){
  Z <- list() 
  Z[[1]] <- I
  j <- 1 
  while(sum(Z[[j]]) > 0 && j <= gens) { 
    lambdaM <- (a * bseas * exp(-mumseas/(PDRseas+ec)))/(rm) * sum(Z[[j]])
    M <- rnbinom(1,  size=k, mu=lambdaM) 
    lambdaH <- a * cseas * 1/(mumseas) * M * (pop-sum(Z[[j]]))/pop
    Z[[j+1]] <- rnbinom(1,  size=k, mu=lambdaH)  
    j <- j+1 
  } 
  out <- c(length=length(Z)-2, total=sum(unlist(Z))-1)
}





