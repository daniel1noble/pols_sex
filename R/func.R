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
		
		if(repeated.control == TRUE){
				mean.control.for.cor<-CMean[match(unique(Control.IDs), Control.IDs)]
				sd.control.for.cor<-CSD[match(unique(Control.IDs), Control.IDs)]
			}
			else{
				mean.control.for.cor<-CMean
				sd.control.for.cor<-CSD
			}
		
		if(Equal.E.C.Corr==TRUE){
				mvcorr<-cor.test(log(c(mean.control.for.cor, EMean)), log(c(sd.control.for.cor, ESD)))$estimate
				S2<- CSD^2 / (CN * (CMean^2)) + 1 / (2 * (CN - 1)) - 2 * mvcorr * sqrt((CSD^2 / (CN * (CMean^2))) * (1 / (2 * (CN - 1)))) + ESD^2 / (EN * (EMean^2)) + 1 / (2 * (EN - 1)) - 2 * mvcorr * sqrt((ESD^2 / (EN * (EMean^2))) * (1 / (2 * (EN - 1))))
			}
			else{
				Cmvcorr<-cor.test(log(mean.control.for.cor), log(sd.control.for.cor))$estimate
				Emvcorr<-cor.test(log(EMean), (ESD))$estimate
		
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