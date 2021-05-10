############################################
## STATISTICS LAB, ETH ZURICH (Spring 2018)
##
## Description:
## Power analysis with simulations to obtain
## sample size for different power values
##
## Authors:
##    - Beck, Elliot
##    - Fuentes, Eufemiano
##    - Park, Jun
############################################



rm(list = ls())

t <- Sys.time()

## Load libraries
suppressWarnings(suppressMessages({
library(bindata)
library(lmerTest)
library(foreach)
library(parallel)
library(doSNOW)}))

## Null hypothesis:
# Rack is not significant (none of rack levels are significant)
# Bypass is not significant
# Species is not significant
# Rack-Species interaction is not significant
# Bypass-Species interaction is not significant

## Function to compute column for interaction effects
interactions <- function(factor1, levels1, 
                         factor2, levels2, 
                         effects, df){
  
  y <- rep(0, times = nrow(df))
  
  # Explore every combination of factors 1 and 2
  for (i in levels1) {
    for (j in levels2) {
      
      # Store configurations (rows) where factor1 = i and factor2 = j
      pos <- which(df[, factor1] == i & df[, factor2] == j)
      
      # Set interaction effect for previous found positions configurations at rows 'pos'
      y[pos] <- effects[i, j]
    }
  }
  
  return(y)
}

## Function to simulate correlated binary responses
simulate <- function(probs, cormats, ngroups, col.species, nspecies = 5,
                     nfish.group = 3, nrack = 2, nbypass = 3) {
  
  nconf <- length(probs)
  season.rows <- nrack*nbypass*nspecies*nfish.group*ngroups
  y <- integer(length = season.rows * 2)
  
  # Generate correlated binary responses
  for (i in 1:nconf) {
    # For every configuration (every i in probs) compute nfish.group correlated binary
    # responses, using the correlation matrix of the fish species in that configuration
    
    # Get correlation matrix of fish species in configuration i
    cormat <- cormats[[col.species[i]]]
    
    # Replicate ngroups times, nfish.group correlated binary responses with probability
    # of configuration i (probs[i]) using correlation matrix of the corresponding fish
    # species
    reps <- rmvbin(n = ngroups,
                   margprob = rep(probs[i], nfish.group),
                   bincorr = cormat
    )
    
    # Store results of configuration i in y
    pos <- seq(from = 0,
               by = nfish.group*nspecies*nbypass*nrack,
               length.out = ngroups)
    pos <- rep(pos, each = nfish.group) + 1:nfish.group
    pos <- pos + ((i-1) %% (nconf/2))*nfish.group + ifelse(i > nconf/2, season.rows, 0)
    y[pos] <- as.vector(t(reps))
  }
  
  return(y)
}


# ---------------------------------------------------------------
# Assumptions about the model under the alternative hypothesis
# ---------------------------------------------------------------

## Settings
set.seed(123)
ngroups <- 4  # number of groups
nfish.group <- 3  # number of fish per group
mu <- 50
ncores <- detectCores()*3/4

## Levels of variables
season.levels  <- factor(c("Autumn", "Spring"))
rack.levels    <- factor(c("vertical", "horizontal"))
bypass.levels  <- factor(1:3)
species.levels <- factor(1:5)
levels(species.levels) <- 1:6
group.levels   <- factor(1:ngroups)
sample.levels  <- factor(1:nfish.group)

## Effect of variables on efficiency (already in %, will divide by 100 later)
season.effects  <- c(0, 0)
rack.effects    <- c(20, 0)
bypass.effects  <- c(10, 0, 0)
species.effects <- c(0, 0, 0, 0, 0, 0)  # 1:4 in both seasons, 5 in Autumn, 6 in Spring
cors            <- c(0.01, 0.25, 0.4, 0.66, .99, 0.8)

rack.species.effects   <- matrix(c(0, 0, 0, 0, 0, 0, 
                                   0, 0, 0, 0, 0, 0), nrow = 2, ncol = 6, byrow = TRUE)
bypass.species.effects <- matrix(c(0, 0, 0, 0, 0, 0, 
                                   0, 0, 0, 0, 0, 0, 
                                   0, 0, 0, 0, 0, 0), nrow = 3, ncol = 6, byrow = TRUE)

rownames(rack.species.effects) <- rack.levels
colnames(rack.species.effects) <- levels(species.levels)
rownames(bypass.species.effects) <- bypass.levels
colnames(bypass.species.effects) <- levels(species.levels)


# -----------------------------------------------------------------------------
# Compute correlation matrices, data frame for simulated data and efficiencies
# -----------------------------------------------------------------------------

# Compute correlation matrices for all fish species
cormats <- list()
for (i in 1:length(cors)) {
  corr <- cors[i]
  cormat <- matrix(rep(corr, nfish.group^2), nrow = nfish.group)
  diag(cormat) <- rep(1, nfish.group)
  cormats[[i]] <- cormat
}

## Create data frame for data from simulations (y's will be simulated inside the for loop)
data <- expand.grid(sample  = sample.levels,
                    species = species.levels,
                    bypass  = bypass.levels,
                    rack    = rack.levels,
                    season  = season.levels)
data <- data[, ncol(data):1]

## Change fish species level in Spring to 6
pos <- which(data[, "season"] == "Spring" & data[, "species"] == 5)
data[pos, "species"] <- 6

## Replicate for each season ngroups times
autumn <- data[rep(which(data$season == "Autumn"), times = ngroups), ]
spring <- data[rep(which(data$season == "Spring"), times = ngroups), ]
data <- rbind(autumn, spring)
rownames(data) <- 1:nrow(data)

## Add a column with days
ndays <- length(rack.levels)*length(bypass.levels)*ngroups
day <- rep(1:ndays, each = length(species.levels)*nfish.group)
day <- rep(day, times = length(season.levels))
data <- cbind(day, data)
data <- data[, c("season", "day", "rack", "bypass", "species", "sample")]

## Compute interaction effects for every interaction
df.nosamples <- expand.grid(species = species.levels,
                            bypass  = bypass.levels,
                            rack    = rack.levels, 
                            season  = season.levels)

## Change fish species level in Spring to 6
pos <- which(df.nosamples[, "season"] == "Spring" & df.nosamples[, "species"] == 5)
df.nosamples[pos, "species"] <- 6
df.nosamples <- df.nosamples[, ncol(df.nosamples):1]

inter.rack.species <- interactions(factor1 = "rack", levels1 = rack.levels, 
                                   factor2 = "species", levels2 = levels(species.levels),
                                   effects = rack.species.effects, df = df.nosamples)
inter.bypass.species <- interactions(factor1 = "bypass", levels1 = bypass.levels, 
                                     factor2 = "species", levels2 = levels(species.levels),
                                     effects = bypass.species.effects, df = df.nosamples)
# check <- cbind(df.nosamples, efficiencies)

## Compute expected efficiencies based on variables' effects
df.eff <- expand.grid(species = species.effects[1:5], 
                      bypass  = bypass.effects, 
                      rack    = rack.effects, 
                      season  = season.effects)

## Change fish species 6 effect
df.eff[pos, "species"] <- species.effects[6]

df.eff <- cbind(df.eff, inter.rack.species, inter.bypass.species)  # add interactions
efficiencies <- (mu + rowSums(df.eff)) / 100  # compute efficiency per configuration (row)


if (max(efficiencies) > 1 || min(efficiencies) < 0) {
  stop("ERROR: Some configuration has an efficiency > 1 or < 0")
} else {
  
  # -----------------------------------------------
  # Run simulations in parallel to obtain power
  # -----------------------------------------------
  
  ## Initialize number of simulations and alpha (1-confidence level)
  nsim <- 6L*5L
  alpha <- 0.05
  
  ## Set contrasts to sum to zero constraint
  options(contrasts = c(unordered = "contr.sum", ordered = "contr.sum"))
  
  ## Register the parallel backend
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)
  
  ## Setup for progress bar
  pb <- txtProgressBar(max = nsim, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  ## Run simulations in parallel
  results <- foreach(i=1:nsim, 
                     .packages = c("bindata", "lmerTest"), 
                     .combine = 'rbind', 
                     .options.snow = opts, 
                     .inorder = FALSE
  ) %dopar% {
    
    ## Simulate new response
    data$y <- simulate(probs = efficiencies, 
                       cormats  = cormats, 
                       ngroups  = ngroups, 
                       col.species = df.nosamples$species
    )
    
    ## Fit Generalized Linear Mixed-Effects model, only main effects
    fit <- glmer(y ~ season + (1 | season:day) + rack + bypass + species + 
                   (1 | season:day:species), data = data, family = binomial()
    )
    
    ## Get p-values for main effects except species using LRT (Chisq)
    pvalues <- drop1(fit, test = "Chisq")[[4]]
    
    ## Remove non-overlapping species (5 and 6)
    data.species <- data[-which(data$species == 5 | data$species == 6), ]
    
    ## Fit model again without non-overlapping fish species
    fit.species <- update(fit, data = data.species)

    ## Get p-value for species effect using LRT (Chisq) and filtered data
    pvalue.species <- drop1(fit.species, test = "Chisq")[[4]][5]
    
    ## Fit Generalized Linear Mixed-Effects models with interactions
    fit.rack.inter   <- update(fit, . ~ . + rack:species)
    fit.bypass.inter <- update(fit, . ~ . + bypass:species)

    ## Get p-values for interaction effects using LRT (Chisq)
    pvalue.rack.inter   <- anova(fit, fit.rack.inter, test = "Chisq")[[8]][2]
    pvalue.bypass.inter <- anova(fit, fit.bypass.inter, test = "Chisq")[[8]][2]

    ## Return results for main effects and interactions
    return(as.numeric(c(pvalues[3] < alpha,  # p-value of rack
                        pvalues[4] < alpha,  # p-value of bypass
                        pvalue.species < alpha,  # p-value of species (filtered)
                        pvalue.rack.inter < alpha,  # p-value of rack.species
                        pvalue.bypass.inter < alpha)))  # p-value of bypass.species
  }
  close(pb)
  stopCluster(cl)
  colnames(results) <- c("rack", "bypass", "species", "rack:species", "bypass:species")
}

Sys.time() - t

# ------------------------------------
# Results from simulation
# ------------------------------------

# Power = proportion of cases where we actually reject H_0
colMeans(results)

## Confidence intervals for power
# (power.rack.CI <- t.test(results[, "rack"])$conf)
# (power.bypass.CI <- t.test(results[, "bypass"])$conf)


