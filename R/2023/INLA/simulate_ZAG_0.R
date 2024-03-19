res_mesh = 150
res_pred = "station"
trans = NULL
myfamily = "ZAG"
minw = -10
maxw = 10
smoother = F
Rvar = "DEP"
year_to_report = year_to_report
trajet = 1
serie=""
name = paste0(Rvar, "_INLA_mesh", res_mesh,"km_pred_",res_pred,"_",myfamily,year_to_report,"_tajet",trajet, serie)

dir="results/INLA/"

load(paste0(dir,name,".Rdata"))

# Section 7: Can the ZAG cope with 85% of zeros? ------------

#' Can the ZAG cope with the 85% of zeros in the observed data?
#' To answer this question, we will:
#'  1. Simulate 1000 sets of regression parameters and random effects from the Bernoulli GAM
#'  2. Simulate 1000 sets of regression parameters and random effects from the Gamma GAM
#'  3. For each simulated set of regression parameters calculate Pi and mu, use these to simulate 
#'     Bernoulli and Gamma data, and use these to calculate ZAG expected values.
#'  4. Once we have 1000 simulated ZAG data sets, we determine the number of zeros
#'     for each simulated data set..and compare these to the number of zeros for the 
#'     observed data.



#* Section 7.1: Simulate ------------


#' When we executed the Bernoulli and GAM models, we used config = TRUE.
#' So, we do not have to do that again. To use similar code as in earlier
#' exercises, we use:
I1.sim <- Gbern
I2.sim <- Ggamma



#' Simulate regression parameters from their posterior distributions
NSim <- 1000 # Reduce computing time, we use only 50 simulations.
SimData01  <- inla.posterior.sample(n = NSim, result = I1.sim) #Simulation Bernoulli GAM
SimDataPos <- inla.posterior.sample(n = NSim, result = I2.sim) #Simulation Gamma GAM




#* Section 7.2: Where are the simulated betas? ------------

#' That was the easy part of the code. Now it gets ugly; we need to figure out 
#' where in the simulated output we can find the betas.
#' Option 1: Type SimData01[[1]]$latent and get the row numbers manually.
#' Option 2: Close your eyes and run the code below (it is ugly).

#' Following option 2:
#'  First we need to determine on which rows in SimData[[1]]$latent 
#'  the betas are.
#'  Problem: Some INLA versions put ':1' at the end of the
#'           names of the regression coefficients (the betas). 
#'           And that may cause some problems with the 'grep'
#'           function.
#'  Our solution:
#'    1. Determine the name of the last regression parameter (LastBeta)
#'    2. See whether its last 2 characters equal ':1'. (Last2Character)
#'    3. If so, then add ':1' to our beta names (BetasInModel), resulting in BetasInINLA
#'    4. Determine the row numbers in SimData[[1]]$latent where we can find our
#'       regression parameters (BetasInINLA): BetaRows

#' 1.
nd <- length(rownames(SimData01[[1]]$latent))
LastBeta <- rownames(SimData01[[1]]$latent)[nd]

#' 2.
substrRight <- function(x, n){ substr(x, nchar(x)-n+1, nchar(x)) }
Last2Character <- substrRight(LastBeta, 2)

#' 3. 
MyID2 <- function(x, SimulatedData){ which(rownames(SimulatedData[[1]]$latent) == x) }
Betas01InModel  <- rownames(I1.sim$summary.fixed) #<--Change for your model
BetasPosInModel <- rownames(I2.sim$summary.fixed) #<--Change for your model

Betas01InINLA  <- Betas01InModel
BetasPosInINLA <- BetasPosInModel
if (Last2Character == ":1") { 
  Betas01InINLA  <- paste(Betas01InModel, ":1", sep ="") 
  BetasPosInINLA  <- paste(BetasPosInModel, ":1", sep ="") 
}

#' 4.
Beta01Rows  <- lapply(Betas01InINLA, MyID2, SimulatedData = SimData01)
BetaPosRows <- lapply(BetasPosInINLA, MyID2, SimulatedData = SimDataPos)
Beta01Rows  <- as.numeric(Beta01Rows)
BetaPosRows <- as.numeric(BetaPosRows)
Beta01Rows
BetaPosRows





#* Section 7.3: Where are the simulated ws? ------------

#' Now the spatial correlated random effects.
MyGrep <- function(x, SimulatedData){ 
  # Small function to determine the row numbers of random effects 
  # in the output of inla.posterior.sample
  # This is partly code from Aaron Adamack (Aaron.Adamack@dfo-mpo.gc.ca)
  # SimulatedData is the name of the object containing the simulation results
  # x is an element of BetasInINLA
  names.SimData <- attributes(SimulatedData[[1]]$latent)$dimnames[[1]]
  names.SimData[grep(x, names.SimData)] 
}


#' Determine on which rows the spatial correlated random effects of the Bernoulli
#' GAM are.
names.SimData01  <- attributes(SimData01[[1]]$latent)$dimnames[[1]] 
wNames <- names.SimData01[grep("w:", names.SimData01)]  #<---Note the name of the w
head(wNames) #First 6 names. 
tail(wNames) #Last 6 names

#' Determine on which rows these random effects are
MyID  <- function(x){ which(rownames(SimData01[[1]]$latent) == x) }
RowsW01 <- lapply(wNames, MyID)
RowsW01 <- as.numeric(RowsW01)
SimData01[[1]]$latent[RowsW01]  #These are the simulated w_ks in the Bernoulli model

#' Do the same for the Gamma spatial random effects.
names.SimDataPos  <- attributes(SimDataPos[[1]]$latent)$dimnames[[1]] 
wNames <- names.SimDataPos[grep("wGamma:", names.SimDataPos)]  #<---Note the name of the w
MyID  <- function(x){ which(rownames(SimDataPos[[1]]$latent) == x) }
RowsWPos <- lapply(wNames, MyID)
RowsWPos <- as.numeric(RowsWPos)
SimDataPos[[1]]$latent[RowsWPos]  #These are the simulated w_ks in the Gamma model.



#* Section 7.4: Get the As ------------

#' Get the A matrix to convert the spatial correlated random
#' effects w into the u with u = A2 * w
ABern  <- as.matrix(A.Repl)
AGamma <- as.matrix(A.Repl)



#* Section 7.5: Get the X matrices ------------

#' We also need to get an X matrix (X01) containing all covariates for the
#' binary part, and an X matrix (XPos) containing all covariates for
#' the gamma part.

#' X matrix for the binary part
X01    <- model.matrix(~  -1 + station ,
                       data = dat)
X01   <- cbind("(Intercept)"=1,as.matrix(X01))





#' X matrix for the Gamma part (note that we are using the same covariates in both parts)
XPos    <- model.matrix(~ -1+ station,
                        data = dat)
XPos   <- cbind("(Intercept)"=1,as.matrix(XPos))



#* Section 7.6: Start a loop ------------

#' Start a loop to extract betas and random effects, calculate
#' the fitted values for the Bernoulli and Gamma models, and 
#' simulate data from the model.
N    <- nrow(dat)
Ysim <- matrix(nrow = N, ncol = NSim)
for (i in 1: NSim){
  Betas01    <- SimData01[[i]]$latent[Beta01Rows] #Betas for Bernoulli part
  w01        <- SimData01[[i]]$latent[RowsW01]       #w for Bernoulli part
  BetasPos   <- SimDataPos[[i]]$latent[BetaPosRows]  #Betas for gamma part
  wPos       <- SimDataPos[[i]]$latent[RowsWPos]     #w for gamma part
  
  #' Get the simulated parameter r from the gamma variance
  r <- SimDataPos[[i]]$hyperpar['Precision parameter for the Gamma observations']
  
  #' Calculate the Bernoulli and Gamma expected values
  Pi <- exp(X01 %*% Betas01 + ABern %*% w01) / (1 + exp(X01 %*% Betas01 + ABern %*% w01))
  mu <- exp(XPos %*% BetasPos + AGamma %*% wPos)
  
  #' Simulate a vector with zeros and ones 
  W <- rbinom(N, size = 1, prob = Pi)
  
  #' Simulate Gamma distributed data
  #' Help file rgamma:
  #'   a = shape, s = scale
  #'   The mean and variance of a Gamma distribution are 
  #'   E(X) = a*s and Var(X) = a*s^2.
  #' Here, use: shape = r   and  scale = mu / r
  Ygamma <- rgamma(n = N, shape = r, scale = mu/r )
  
  #' Calculate W * Ygamma
  Ysim[,i] <- W * Ygamma
}

#' The first 6 rows of the first 10 simulated data sets
head(Ysim[,1:10])



#* Section 7.7: Determine % of zeros ------------

#' Calculate the % of zeros in each of the simulated data sets, and also 
#' in the original data.
N <- nrow(dat)
ZerosInData <- 100 * sum(dat$DEP == 0) / N
ZerosInData

#' In the simulated data:
Zeros <- vector(length = NSim)
for(i in 1:NSim){
  Zeros[i] <- 100 * sum(Ysim[,i] == 0) / N
}



#* Section 7.8: Plot simulation results ------------
png("figures/INLA/models/simul0.png", width=5, height=5, units="in", res=300)
par(mar = c(5,5,2,2), cex.lab = 1.5)
hist(Zeros, 
     xlim = c(0, 100),
     xlab = "Percentage of zeros",
     ylab = "Frequency",
     main = "Simulation results")
points(x = ZerosInData, 
       y = 0, 
       pch = 16, 
       cex = 5, 
       col = 2)
dev.off()
#' The red dot is the percentage of zeros in the original data set.
#' The histogram shows the percentages of zeros for simulated ZAG data
#' from the model.
#' Conclusion: The model can cope with the large number of zeros.




