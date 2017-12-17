#--------------------------------------------------------------------------------------------#
# Various functions for calculating relevant effects sizes for the paper
# Authors: Written by A. Senior (@Sydney Uni), D Noble, and S. Nakagawa
# Thu Jan 26 11:56:30 2017
# Purpose: POLs sex-differences meta-analysis
#--------------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------------------#
## Code modified from A. Senior (@Sydney Uni). Calculate log response ratio 
## and sampling error variance
#---------------------------------------------------------------------------------------------#
	lnRR_es<-function(CMean, CN, CSD, EMean, EN, ESD){	
		lnRR<-log(EMean / CMean)
		EVar<-ESD^2
		CVar<-CSD^2
		v.lnRR<- (CVar / (CN * (CMean^2))) + (EVar / (EN * (EMean^2)))
		return(cbind(lnRR, v.lnRR))
	}

#-----------------------------------------------------------------------------------------------#
## Code from A. Senior (@Sydney Uni). Calculate log coefficient of variance ratio effect sizes. 
## This controls for changes in the mean.
#-----------------------------------------------------------------------------------------------#
	lnCVR<-function(CMean, CSD, CN, EMean, ESD, EN){	
		ES<-(log(ESD) - log(EMean) + 1 / (2*(EN - 1))) - (log(CSD) - log(CMean) + 1 / (2*(CN - 1)))
		return(ES)
	}

	v_lnCVR<-function(CMean, CSD, CN, EMean, ESD, EN, Equal.E.C.Corr=TRUE, repeated.control = FALSE, Control.IDs){
		
		if(Equal.E.C.Corr==TRUE){
				mvcorr<-cor.test(log(c(CMean, EMean)), log(c(CSD, ESD)))$estimate
				S2<- CSD^2 / (CN * (CMean^2)) + 1 / (2 * (CN - 1)) - 2 * mvcorr * sqrt((CSD^2 / (CN * (CMean^2))) * (1 / (2 * (CN - 1)))) + ESD^2 / (EN * (EMean^2)) + 1 / (2 * (EN - 1)) - 2 * mvcorr * sqrt((ESD^2 / (EN * (EMean^2))) * (1 / (2 * (EN - 1))))
			}
			else{
				Cmvcorr<-cor.test(log(CMean), log(CSD))$estimate
				Emvcorr<-cor.test(log(EMean), log(ESD))$estimate
		
				S2<- CSD^2 / (CN * (CMean^2)) + 1 / (2 * (CN - 1)) - 2 * Cmvcorr * sqrt((CSD^2 / (CN * (CMean^2))) * (1 / (2 * (CN - 1)))) + ESD^2 / (EN * (EMean^2)) + 1 / (2 * (EN - 1)) - 2 * Emvcorr * sqrt((ESD^2 / (EN * (EMean^2))) * (1 / (2 * (EN - 1))))
			}
		return(S2)
	}

#-----------------------------------------------------------------------------------------------#
## Code modified from  (@Sydney Uni). Calculating lnVR effect size and sampling error variance. 
## Note that this does not control for any changes in the mean across the sexes and so 
## it is useful to compare with lnCVR
#-----------------------------------------------------------------------------------------------#
	lnVR_es<-function(ESD, EN, CSD, CN){
		lnVR <- log(ESD / CSD) + (1 / (2 * (EN - 1))) - (1 / (2 * (CN - 1)))
		var.lnVR <- (1 / (2 * (CN - 1))) + (1 / (2 * (EN - 1)))
		return(cbind(var.lnVR, lnVR))
	}

#------------------------------------------------------------------------------------------------#
# Multi-model inference - glmulti and metafor
#------------------------------------------------------------------------------------------------#

  # Used for generating a model set.
	rma.mv_glmulti <- function(formula, data, vi = vi, random = random, ...){
		rma.mv(formula, vi, random = random, data = data, method = "ML", ...)
	}

#-------------------------------------------------------------------------------------------------#
# Useful for rounding data frames
#-------------------------------------------------------------------------------------------------#
	round_df <- function(x, digits) {
	    # round all numeric variables
	    # x: data frame 
	    # digits: number of digits to round
	    numeric_columns <- sapply(x, class) == 'numeric'
	    x[numeric_columns] <-  round(x[numeric_columns], digits)
	    x
	}


# Extract coefficients from metafor model
	extract <- function(model){
		coefs <- summary(model)[c("b", "se", "ci.lb", "ci.ub")]
			coefTab <- data.frame(matrix(nrow = length(model[[1]]), ncol = length(coefs)))
			for(i in 1:ncol(coefTab)){
				coefTab[[i]] <- coefs[[i]]
				colnames(coefTab) <- c("Est.", "se", "ci.lb", "ci.up")
				rownames(coefTab)<- rownames(coefs[[1]])
			}
			return(coefTab)
		}

# Will re-run the top models identified with glmulti and then extract all teh relevant stats needed for model averaging.


		TopModelEst <- function(modTable_glmulti){
			dataTab <-c()
			
			for(i in 1:nrow(mod_tab_lnCVR)){
				form <- formula(as.character(mod_tab_lnCVR[i,"model"]))
				modname <- rma.mv(form, V = VlnCVR, random = list(~1|study, ~1|species, ~1|obs), R = list(species = phylo_cor), data = data)
				extraction <- extract(modname)
				      aicc <- AICc(modname)
				      Estname <- as.character(rownames(extraction))
				dataTab <- rbind(dataTab, cbind(Estname, extraction, aicc, i))
			}

			dataTab$Estname <- as.character(dataTab$Estname)
			rownames(dataTab) <- 1:nrow(dataTab)
			dataTab$deltAICc <- dataTab$aicc - min(dataTab$aicc)
			dataTab$weight <- exp(-0.5*dataTab$deltAICc) / sum(exp(-0.5*unique(dataTab$deltAICc)))
			sum(unique(dataTab$weight))

			return(dataTab)
		}


# Model averaging function with top models within metafor. Using natural model averaging (i.e. only averages estimates over models they occur within.)
	mod_avg <- function(models, type = c("full", "natural")){
				#Grab unique variables in he 
				vars <- as.character(unique(models$Estname))
				modWeigts <- unique(models[,c("i", "deltAICc", "weight")])
				avergedEst <- c()
				
				type = match.arg(type)

				if(type == "natural"){
				#Extract the same estimates from big table, calculate model average parameter and average SEs
				for(i in 1:length(vars)){
					estimates <- models[which(models$Estname == vars[i]),]
					   avgEst <- sum(estimates$weight*estimates$Est.) / sum(estimates$weight)
					uncondSe <- sum((estimates$weight)*( sqrt((estimates$se)^2 + (estimates$Est. - mean(estimates$Est.))^2)))
					avergedEst <- rbind(avergedEst, c(avgEst, uncondSe))
					}
				}

				if(type == "full"){

					for(i in 1:length(vars)){
					estimates <- models[which(models$Estname == vars[i]),]
					estimates$Est. <- as.vector(estimates$Est.)
					     name <- unique(estimates$Estname)
					   
					   if(length(estimates$i) != length(modWeigts$i)){
					   estNot <- modWeigts[which(!modWeigts$i %in% estimates$i),]
					   estNotdata <- as.data.frame(cbind(name, matrix(data=0, nrow = nrow(estNot),ncol = ncol(estimates)-4), estNot))
					   colnames(estNotdata) <- colnames(estimates)
					   estDat <- rbind(estimates, estNotdata)
					   } else{estDat = estimates}

					   avgEst <- sum(estDat$weight*estDat$Est.)
					uncondSe <- sum((estDat$weight)*( sqrt((estDat$se)^2 + (estDat$Est. - mean(estDat$Est.))^2)))
					avergedEst <- rbind(avergedEst, c(avgEst, uncondSe))
					}
				}

				rownames(avergedEst) <- vars
				colnames(avergedEst) <- c("AvgEst.", "UncondSE")

			return(avergedEst)
		}


I2 <- function(model, v, sims = 1500, phylo = FALSE){
	
	if(class(model) != "MCMCglmm" && class(model) != "rma.mv" && class(model) != "rma"){
		stop("The model object is not of class 'MCMCglmm' or 'metafor'")
		}

		wi <- 1/v  #weight
		Vw <- sum((wi) * (length(wi) - 1))  / (((sum(wi)^2) - sum((wi)^2)))

	if("MCMCglmm" %in% class(model)){
		# Get posterior distribution 
		post <- model$VCV[,-match(c("sqrt(mev):sqrt(mev).meta"), colnames(model$VCV))]

		#Calculate total variance
		         VT <- rowSums(Matrix::cBind(post, Vw))
		          Vt <- rowSums(post)  # remove Vw
		
		# For each variance component divide by the total variance. Note this needs to be fixed for phylo, but does deal with variable random effects.
		  I2_re <- post / VT
		  I2_total  <- Vt / VT

		if(phylo == FALSE){
			tmpMatrix <- Matrix::cBind(I2_re, total = I2_total)
		}else{
		  	I2_phylo <- post[,match(phylo, colnames(sigma2))] / Vt
		  	tmpMatrix <- Matrix::cBind(I2_re, I2_phylo, total = I2_total)
		  	}

		   mode <- MCMCglmm::posterior.mode(coda::as.mcmc(tmpMatrix))
		   CI <- coda::HPDinterval(coda::as.mcmc(tmpMatrix))
		    colnames(CI) <- c("2.5% CI", "97.5% CI")

		    I2_Table <- as.data.frame(Matrix::cBind(I2_Est. = mode[-match("units", names(mode))], CI[-match("units", rownames(CI)),]))
		    class(I2_Table) <- c("metaAidR", "data.frame")

	return(round_df(I2_Table, digits = 4))
  	}

  	if("rma.mv" %in% class(model) | "rma" %in% class(model)){
  		#Monte Carlo Simulations
		# From metafor extract the important statistics
  		sigma2 <- matrix(model$sigma2, nrow = 1, ncol = length(model$sigma2))
  		colnames(sigma2) <- model$s.names
  		sigmaN <- model$s.nlevels

  		if("obs" %in% colnames(sigma2) == FALSE){
  			stop("The metafor object does not contain a residual variance estimate. Please include an observation-level random effect (~1|obs) when fitting model")
  		}

  		#For each variance estimate use Monte Carlo simulation of data
  		Sims <- data.frame(mapply(function(x,y) simMonteCarlo(x, y, sims = sims), sigma2, sigmaN))
		colnames(Sims) <- colnames(sigma2) 
		
		#Calculate total variance
		VT <- rowSums(Matrix::cBind(Sims, Vw))
		Vt <- rowSums(Sims)  # remove Vw
		
		# For each variance component divide by the total variance. Note this needs to be fixed for phylo, but does deal with variable random effects.
		 I2_re       <- Sims / VT
		 I2_total   <- Vt / VT

		  if(phylo == FALSE){
		  	## FIXED PROBLEM WITH COLUMN NAMES in I2. FIXED: Needs to be called data.frame.
		  	tmpMatrix <- as.data.frame(Matrix::cBind(I2_re[,-match("obs", colnames(I2_re))], I2_total))
		  		colnames(tmpMatrix) <- c(colnames(I2_re)[-match("obs", colnames(I2_re))], "total")
		   }else{
		  	I2_phylo <- Sims[, match(phylo, colnames(sigma2))] / Vt
		  	 tmpMatrix <- as.data.frame(Matrix::cBind(I2_re[,-match("obs", colnames(I2_re))], phylo = I2_phylo, total = I2_total))
		  }

		CI <- lapply(tmpMatrix, function(x) stats::quantile(x, c(0.025, 0.975), na.rm = TRUE)) ## Problem here!! is producing many CI's when it should match the cols.
		I_CI <- as.data.frame(do.call(rbind, CI))
		colnames(I_CI) <- c("2.5% CI", "97.5% CI")
		I2_table <- Matrix::cBind(I2_Est. = colMeans(tmpMatrix), I_CI )
		
		class(I2_table) <- c("metaAidR", "data.frame")

	return(round_df(I2_table, digits = 4))
  	}

}


#' @title Parametric simulation 
#' @description Function for calculating I2 estimates using parametric simulations of model estimates taken from metafor. Note that the effectiveness of these simulations depends on the accuracy of model variance estimates.
#' @param estimate The estimate (i.e. variance) from a metafor model
#' @param sims The number of simulations 
#' @param n The sample size used in estimating the variance 
#' @author Daniel Noble - daniel.noble@unsw.edu.au
#' @export
  simMonteCarlo <- function(estimate, n, sims){
  		set.seed(07)
  		tmp <- data.frame(num = rep(1:sims, each = n), y = stats::rnorm(n*sims, 0, sqrt(estimate)))
  		Var <- dplyr::summarise(dplyr::group_by(tmp, num), var = stats::var(y))
  		return(as.numeric(Var$var))
  	}


 # Functions for marginal estimates
 		
 	marginalize <- function(mod, vars){
 		modlist <- list()
		for(i in 1:length(vars)){
			modlist[[i]] <- effects::Effect(focal.predictors = vars[i], mod = mod, se = TRUE, confidence.level = 0.95)
		}
		return(modlist)
	}

	 # Takes a list, grab all estimates and then merges them together!
	 margTable <- function(modlist){
				as.data.frame(do.call(rbind, lapply(modlist, function(x) margEst(x))))

			}
	# Takes an "effects" object	
	 margEst <- function(margEsts){
				do.call(cbind, summary(margEsts)[c("effect", "lower", "upper")])
			}


	# Generates a within study covariance matrix
	covMatrix <- function(data, es_var, cor){
	 	tmp <- expand.grid(sqrt(data[, es_var]), sqrt(data[, es_var]))
	 	tmp$cor <- ifelse(tmp$Var1 == tmp$Var2, 1, 0.5)
		tmp$cov <- tmp$cor * tmp$Var1 * tmp$Var2
	  	corMat <- matrix(tmp$cov , nrow = nrow(data), ncol = nrow(data))
  	return(corMat)
	}

#'@title Function for calculating covariance between dependent effect sizes
#'@param data the data frame
#'@param es_var the character identifier effect size sampling error variance column in the data
#'@param depend the character identifier for the column in the data indicating which pairs of effect sizes are correlated
#'@param cor The correlation coefficient assumed for calculating the covariance between effect sizes. Shared control can be calculated directly, but is not yet implemented.
#' IMPORTANT: Note that the data frame SHOULD be arranged BEFORE using function on the dependence column. The reason is that sometimes you can get grouping spread out across a few rows. As of yet, I don't have checks for this and so making sure the data is arranged will ensure that the ORDER of the variances in the dataset match the cov matrix. 

VmCovMat <- function(data, es_var, depend, cor = 0.5){
    	data$dep<-paste(data[,"study"], data[,depend], sep="_")

    	tmp <- reshape::expand.grid.df(data.frame(sd1 = sqrt(data[, es_var]), stdy1 = data[,"dep"]), data.frame(sd2 = sqrt(data[, es_var]), stdy2 = data[,"dep"]))

    	tmp$cor <- ifelse((tmp$stdy1 == tmp$stdy2) & (tmp$sd1 != tmp$sd2), cor, 0)
		tmp$cov <- tmp$cor * tmp$sd1 * tmp$sd2
	  	corMat <- matrix(tmp$cov , nrow = nrow(data), ncol = nrow(data))
	  	diag(corMat) <- data[,es_var]
  	
  	return(corMat)
}
    
       