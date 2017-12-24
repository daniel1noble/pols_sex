#----------------------------------------------------------------------------#
# Meta-Analysis of sex differences in behaviour, physiology to test sex-dependent POL
# Authors: Maja Tarka, Anja Günther, Petri Niemela, Shinichi Nakagawa and Daniel Noble
# 2017-04-03 15:53:29 AEST
#----------------------------------------------------------------------------#

	# clear working space
		rm(list=ls())

	# load functions & libraries
		source("./pols_sex/R/func.R")
		library(lattice)
		library(phytools)
		library(nlme)
		library(effects)
		library(dplyr)
		library(VIM) # Helps visualise missing data.
		library(metafor)
		library(Hmisc)
		library(rotl)
		library(devtools)
		library(MCMCglmm)
		library(MuMIn)
		library(corrplot)

	# load data
		data <- read.csv("./pols_sex/data/POLS_sex_db_20171222.csv", stringsAsFactors = FALSE)

# 1. Data Exploration & Processing
#----------------------------------------------------------------------------#
	# Exclude missingness.
		data <-data[data$include_2 == 1,] 
		
	#Calculate effect sizes and explore some more.
		#lnRR
		lnRR_estat <- lnRR_es(CMean = data$Mean_M, CN = data$M_n, CSD = data$SD_M,EMean = data$Mean_F, EN = data$F_n, ESD = data$SD_F )
		data <- cbind(data, lnRR_estat)

		#lnVR
		lnVR_estat <- lnVR_es(CN = data$M_n, CSD = data$SD_M, EN = data$F_n, ESD = data$SD_F)
		data <- cbind(data, lnVR_estat)

		#lnCVR
		lnCVR_es<- lnCVR(CMean = data$Mean_M, CN = data$M_n, CSD = data$SD_M, EMean = data$Mean_F, EN = data$F_n, ESD = data$SD_F)
		VlnCVR <- v_lnCVR(CMean = data$Mean_M, CN = data$M_n, CSD = data$SD_M, EMean = data$Mean_F, EN = data$F_n, ESD = data$SD_F)

		data <- cbind(data, cbind(lnCVR_es, VlnCVR))
 
	# Noticed some issues with categories labeled wrong; fix these
		      data$SSD <- ifelse(data$SSD == "no difference", "No difference", data$SSD)
		   data$mating <- ifelse(data$mating == "unkonwn", "Unknown", data$mating)
		   data$mating <- ifelse(data$mating == "polygyny", "Polygyny", data$mating)

		   data$mating <- ifelse(data$mating == "promiscuity", "Promiscuity", data$mating)
		   data$mating <- ifelse(data$mating == "social monogamy", "Monogamy", data$mating)
		   data$mating <- ifelse(data$mating == "monogamy", "Monogamy", data$mating)

		data$parenting <- ifelse(data$parenting == "both", "Both", data$parenting)
		data$parenting <- ifelse(data$parenting == "female", "Female", data$parenting)
		data$parenting <- ifelse(data$parenting == "Female/both", "Both", data$parenting)

		 data$breeding <- ifelse(data$breeding == "Iteroparous", "iteroparous", data$breeding)
		data$breeding <- ifelse(data$breeding == "Semelparous", "semelparous", data$breeding)

		data$climate_sp <- ifelse(data$climate_sp == "Arid", "arid", data$climate_sp)
		data$climate_sp <- ifelse(data$climate_sp == "Tropical", "tropical", data$climate_sp)
		data$climate_sp <- ifelse(data$climate_sp == "cold", "temperate", data$climate_sp)

		data$background1 <- ifelse(data$background == "semiwild", "wild", data$background)
		data$background1 <- ifelse(data$background1 == "domestic", "lab", data$background1)
		data$background1 <- ifelse(data$background1 == "captive", "lab", data$background1)


	# Check out species list
		#resolve species synonymn for Egernia whitii. Also, Sus scrofa subspecies can just be re-named to species so it is not dropped from Phylo construction.
		data$species <- ifelse(data$species == "Egernia whitii", "Liopholis whitii", data$species)
		data$species <- ifelse(data$species == "Sus scrofa domestica", "Sus scrofa", data$species)
		data$species <- ifelse(data$species == "Melospiza g. nigrescens", "Melospiza georgiana", data$species)
		data$species <- ifelse(data$species == "Melospiza g. georgiana", "Melospiza georgiana", data$species)
		data$species <- ifelse(data$species == "Lagopus muta hyperborea", "Lagopus muta", data$species)
		data$species <- ifelse(data$species == "Anolis carolensis", "Anolis carolinensis", data$species)
		data$species <- ifelse(data$species == "Anolis carolensis", "Anolis carolinensis", data$species)
	
	# Remove white space
		data$species <- trimws(data$species)
		 
	# Get some details on the number of species
		  spp <- unique(data$species)
		  spp <- gsub(" ", "_", spp)
		 nSpp <- length(spp)
		taxon <- unique(data$taxon)

	# Check directions
		direct <- gsub("high values =*", "", data$direction)
		direct <- gsub(" ", "", direct)
		data$direct <- direct
		data$direct <- replace(data$direct, which(is.na(data$direct)), "fast")

	# Reverse effect size values if they are slow. This is really only important for lnRR which makes specific predictions about directionality along a fast-slow continuum.
		data$lnRR_2 <- ifelse(data$direct == "slow", data$lnRR*(-1), data$lnRR)

	# Convert to factors
		data$category <- as.factor(data$category)
		data$background1 <- as.factor(data$background1)
		data$mating <- as.factor(data$mating)
		data$parenting <- as.factor(data$parenting)
		data$breeding <- as.factor(data$breeding)
		data$obs <- 1:nrow(data)

	# Get N for all predictors
		vars <- c("category" , "background1" , "mating" , "breeding")

		table(data$category)

		N <- c()
		for(i in 1:length(vars)){
			samp <- as.vector(t(table(data[,vars[i]])))
			   N <- c(N, samp)
		}

# 2. Exploratory plotting
#----------------------------------------------------------------------------#
		
		# Box plot function to plot the various variables quickly. Just change parameters
		#Parameters
		y <- "lnRR_2"
		x <- c("Trait type", "SSD", "Thermy", "Growth Strategy", "Mating system", "Parental Care", "Climate Zone", "Breeding style")
		var <- c("category", "SSD", "thermy", "growth", "mating", "parenting", "climate_sp", "breeding")

		# Function
		par(mfrow = c(0.5*length(var),2), mar = c(4,4,1,1))
		for(i in 1:length(var)){
			form <- as.formula(paste0(y, "~", var[i]))
			boxplot(form, ylab = as.character(y), xlab = as.character(x[i]), data = data, las =1, ylim = c(-1.8,2.2)) -> bp.out1
			text(paste("N = ", bp.out1$n), x = unique(bp.out1$group), y = 2)
		}
	
# 3. Create phylogenetic matrix
#----------------------------------------------------------------------------#
	# Phylogenetic correlation matrix

	# Build tree
		resolve_names <- tnrs_match_names(spp)
		  tree <- tol_induced_subtree(ott_ids = resolve_names$ott_id)
		  tree$tip.label <- gsub("_ott.+", "", tree$tip.label)

	  pdf(width = 5.982379, height = 6.237885, file = "tree_new_data2.pdf")
		par(mar=c(1,1,1,1))
		plot(tree, no.margin = TRUE, type = "phylogram", cex = 0.5)
	  dev.off()
		write.tree(tree, file="./pols_sex/output/tree_new_data2")

		# Now that we have matched species in ROTL, we need to convert the names of these species to the identified ID's in rotl to ensure the phylo matrix and species names match.
		data_names <- gsub("_", " ", firstup(resolve_names$search_string))
		rotl_names <- resolve_names$unique_name

		# Switch all species names within the data to the synonymns from rotl so that a phylogenetic correlation matrix can be included in the models. 
		spp_rotl <- data$species
		for(i in 1:length(data_names)){
		spp_rotl <- ifelse(data_names[i] == data$species, rotl_names[i], spp_rotl)
		}	

		data$spp_rotl <- spp_rotl

	#Remove ott labels on end to make sure to matches species in dataset
		phylo <- read.tree("./pols_sex/output/tree_new_data2")
		phylo$tip.label <- gsub("_ott.+", "", phylo$tip.label)
		phylo <- makeNodeLabel(phylo) 
		is.binary.tree(phylo)

	# Compute branch lengths
		phylo_BL <- compute.brlen(phylo, method = "Grafen", power = 0.5) # Note that larger values of power produce models with lower AICc estimates, but they are all pretty similar. Estimate of phylo variance increases with smaller power estimates. But, not greatly.

	# Create phylogenetic correlation matrix. For metafor. Make sure names match with species. 
		phylo_cor <- vcv(phylo_BL, corr = TRUE)
		names <- gsub("_", " ", rownames(phylo_cor))
		rownames(phylo_cor) <- colnames(phylo_cor) <- names

	# Check out correlation matrix
		corrplot(phylo_cor, tl.col = "black", tl.cex = 0.8)

		Ainv <- inverseA(phylo_BL, nodes = "ALL", scale = TRUE)$Ainv
		names <- gsub("_", " ", rownames(Ainv))
		rownames(Ainv) <- colnames(Ainv) <- names
	
# 4. Multi-level meta-analytic models (MLMA) - intercept only for hetero. 
#----------------------------------------------------------------------------#
	# Study and species level random effects are mostly confounded so study will probably capture most variation anyway, but worth attempting to estimate

	# lnRR
	   modRR_int <- rma.mv(lnRR_2 ~ 1, V = v.lnRR, random = list(~1|study, ~1|spp_rotl, ~1|obs), R = list(spp_rotl = phylo_cor), data = data)

	   # Sensitivity analysis 
		   # Within study correlation matrix between effects sharing the same individuals
		   VmatRR <- VmCorMat(data, "obs", "Dependency.individual.level")
		   rownames(VmatRR) <- colnames(VmatRR) <- data$obs
		   corrplot(as.matrix(VmatRR), tl.col = "black", tl.cex = 0.8, is.corr = FALSE, type = "lower", method = "color")	

			modRR_intSS <- rma.mv(lnRR_2 ~ 1, V = v.lnRR, random = list(~1|study, ~1|spp_rotl, ~1|obs), R = list(spp_rotl = phylo_cor, obs = VmatRR), data = data)	 
		  
	  #Generate heterogeneity measures and CI's. Note "species" is needed to do the correct calculations for phylogeny.
	   		I2(modRR_int, v = data$v.lnRR, phylo = "spp_rotl")

	   #I2 for various trait categories
		   Phys <- subset(data, data$category == "physiology")
		  modRR_intPhys <- rma.mv(lnRR_2 ~ 1, V = v.lnRR, random = list(~1|study, ~1|obs), data = Phys)
		 
	  #Generate heterogeneity measures and CI's. Note "species" is needed to do the correct calculations for phylogeny.
	   		PhysI2 <- I2(modRR_intPhys, v = Phys$v.lnRR)

	   #I2 for various trait categories
		   dev <- subset(data, data$category == "development")
		  modRR_intDev <- rma.mv(lnRR_2 ~ 1, V = v.lnRR, random = list(~1|study, ~1|obs), data = dev)
		 
	  #Generate heterogeneity measures and CI's. Note "species" is needed to do the correct calculations for phylogeny.
	   		DevI2 <- I2(modRR_intDev, v = dev$v.lnRR, phylo = FALSE)

	   #I2 for various trait categories
		   behav <- subset(data, data$category == "behavior")
		  modRR_intbehav <- rma.mv(lnRR_2 ~ 1, V = v.lnRR, random = list(~1|study, ~1|obs), data = behav)
		 
	  #Generate heterogeneity measures and CI's. Note "species" is needed to do the correct calculations for phylogeny.
	   		BehavI2 <- I2(modRR_intbehav, v = behav$v.lnRR, phylo = FALSE)

	    #I2 for various trait categories
		   LH <- subset(data, data$category == "life history")
		  modRR_intLH <- rma.mv(lnRR_2 ~ 1, V = v.lnRR, random = list(~1|study, ~1|obs), data = LH)
		 
	  #Generate heterogeneity measures and CI's. Note "species" is needed to do the correct calculations for phylogeny.
	   		LHI2 <- I2(modRR_intLH, v = LH$v.lnRR, phylo = FALSE)

	# lnCVR
	   modCVR_int <- rma.mv(lnCVR_es ~ 1, V = VlnCVR , random = list(~1|study, ~1|spp_rotl, ~1|obs), R = list(spp_rotl = phylo_cor), data = data)

	   # Sensitivity analysis 
	   	# Within study correlation matrix between effects sharing the same individuals
	   		VmatCVR <- VmCorMat(data, "obs", "Dependency.individual.level")
	   		rownames(VmatCVR) <- colnames(VmatCVR) <- data$obs
		    corrplot(as.matrix(VmatCVR), tl.col = "black", tl.cex = 0.8, is.corr = FALSE, type = "lower", method = "color")

	   #Generate heterogeneity measures and CI's. Note "species" is needed to do the correct calculations for phylogeny.
	   I2(modCVR_int, v = data$VlnCVR, phylo = "spp_rotl")

	   #I2 for various trait categories
	  modCVR_intPhys <- rma.mv(lnCVR_es ~ 1, V = VlnCVR, random = list(~1|study, ~1|obs), data = Phys)
	 
	  #Generate heterogeneity measures and CI's. Note "species" is needed to do the correct calculations for phylogeny.
	   PhysCVRI2 <- I2(modCVR_intPhys, v = Phys$VlnCVR)

	   #I2 for various trait categories
	  modCVR_intDev <- rma.mv(lnCVR_es ~ 1, V = VlnCVR, random = list(~1|study, ~1|obs), data = dev)
	 
	  #Generate heterogeneity measures and CI's. Note "species" is needed to do the correct calculations for phylogeny.
	   DevCVRI2 <- I2(modCVR_intDev, v = dev$VlnCVR)

	   #I2 for various trait categories
	  modCVR_intBehav <- rma.mv(lnCVR_es ~ 1, V = VlnCVR, random = list(~1|study, ~1|obs), data = behav)
	 
	  #Generate heterogeneity measures and CI's. Note "species" is needed to do the correct calculations for phylogeny.
	   BehavCVRI2 <- I2(modCVR_intBehav, v = behav$VlnCVR)

	    #I2 for various trait categories
	  modCVR_intLH <- rma.mv(lnCVR_es ~ 1, V = VlnCVR, random = list(~1|study, ~1|obs), data = LH)
	 
	  #Generate heterogeneity measures and CI's. Note "species" is needed to do the correct calculations for phylogeny.
	   LHCVRI2 <- I2(modCVR_intLH, v = LH$VlnCVR, phylo = FALSE)

	# Create a table of within trait I2 estimates.
		lnRR <- data.frame(rbind(PhysI2, DevI2, BehavI2, LHI2))
		lnCVR <- data.frame(rbind(PhysCVRI2, DevCVRI2, BehavCVRI2, LHCVRI2))

		write.csv(lnRR, file = "lnRR_het_22.12.17.csv")
		write.csv(lnCVR, file = "lnCVR_het_22.12.17.csv")
		
# 5. Marginal estimates / unconditional means for the groups. 
#----------------------------------------------------------------------------#


	# Note the below model, however, assumes now that the residual variance is fixed! So we must remove the estimation of the observation-level variance, or ADD it into the lme fit to get the same results between metafor and lme. But also problems including nested random effects in lme: https://biostatmatt.com/archives/2718 & here: https://stat.ethz.ch/pipermail/r-help/2002-September/025067.html.

	#lnRR
	    # Metafor
			
			modnameRR <- rma.mv(lnRR_2 ~ category + background1 + mating + breeding, V = v.lnRR, random = list(~1|study, ~1|spp_rotl, ~1|obs), method = "REML", R = list(spp_rotl = phylo_cor), data = data)
			AICc(modnameRR)
			coefRRTable1A <- round_df(data.frame(Est. = modnameRR$b, LCI = modnameRR$ci.lb, LCI = modnameRR$ci.ub), digits = 3)

			I2(modnameRR, v = data$v.lnRR, phylo = "spp_rotl")
			

			#Sensitivity Analysis. Covariance matrix.
				modnameRRDep <- rma.mv(lnRR_2 ~ category + background1 + mating + breeding, V = VmatRR, R = list(spp_rotl = phylo_cor, obs = VmatRR), random = list(~1|study, ~1|spp_rotl, ~1|obs), method = "REML", data = data)
				AICc(modnameRRDep)
				coefRRTable1B <- round_df(data.frame(Est. = modnameRRDep$b, LCI = modnameRRDep$ci.lb, LCI = modnameRRDep$ci.ub), digits =3)
			
				TableS1RR <- rbind(coefRRTable1A, coefRRTable1B)
				write.csv(TableS1RR, "./pols_sex/output/tables/TableS1RRBreeding_r1.csv")

			# Model predictions
				# Cerate new data we would like to predict for. 
					newDat.rmv <- expand.grid(list(unique(data$background1), unique(data$category),unique(data$mating), unique(data$breeding)), stringsAsFactors = TRUE)
					colnames(newDat.rmv) <- c("background1", "category","mating", "breeding")

				# Create model matrix
					X <- model.matrix(~category + background1 + mating + breeding, data = newDat.rmv)

				# Generate point estimate predictions and 95% confidence intervals
					newDat.rmv$pred <- X %*% coefficients(modnameRR)
					V = vcov(modnameRR)
					se2 <- rowSums((X %*% V) * X)

					alpha <- qt((1-0.95)/2, df = 304)
					newDat.rmv$CI_L <- newDat.rmv$pred + (-alpha*sqrt(se2))
					newDat.rmv$CI_U <- newDat.rmv$pred + (alpha*sqrt(se2))

					predictions <- newDat.rmv[newDat.rmv$background1 == "wild", ]


		# MCMCglmm

			       data$phylo <- gsub(" ", "_", data$species)
			       colnames(data)[match("trait", colnames(data))] <- "trait2"
			       data$esID <- 1:dim(data)[1]
			       Vmat <- as(solve(diag(data$v.lnRR)), "dgCMatrix")
			       colnames(Vmat) <- rownames(Vmat) <- data$esID
			
					prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002), G2 =  list(V = 1, nu = 0.002), G3 = list(V = 1, fix = 1))) 

				 mod_RR_mcmc <- MCMCglmm(lnRR_2 ~ category + background1 + mating + breeding, random = ~study + spp_rotl + esID, ginverse = list(spp_rotl = Ainv, esID = Vmat), data = data, prior = prior, nitt = 500000, thin = 100, pr = FALSE, family = "gaussian", verbose = FALSE)
				 summary(mod_RR_mcmc)

			# Lets get marginal mean estimates for each group 
			 	# extract fixed effects posterior distribution for each of the parameters estimated (contrasts)
			 	sol_lnRR <- mod_RR_mcmc$Sol

			 	# First we need the sample size for each of the categories as we'll produce a weighted average. We'll need all these throughout the process. 
			 		   n_category <- table(data$category)
			 		n_background1 <- table(data$background1)
			 			 n_mating <- table(data$mating)
			 		   n_breeding <- table(data$breeding)

			# Marginal estimates for each trait category, averaging across the other variables/parameters
			 		
			 	# Trait category
			 		behav_category <- as.mcmc(
			 			sol_lnRR[,"(Intercept)"] +
			 			sol_lnRR[,"background1wild"] * (n_background1[2] / sum(n_background1)) +
			 			sol_lnRR[,"matingPolygyny"] * (n_mating[2] / sum(n_mating)) +  
			 			sol_lnRR[,"matingPromiscuity"] * (n_mating[3] / sum(n_mating)) + 
			 			sol_lnRR[,"breedingsemelparous"] * (n_breeding[2] / sum(n_breeding))
			 			)

			 		# Now that we have an intercept that is essentially unconditioned / averaged across all other levels we can get the marginal estimates for all other levels 
			 			  develop_category <- behav_category + sol_lnRR[,"categorydevelopment"]
			 			life_hist_category <- behav_category + sol_lnRR[,"categorylife history"]
			 			     phys_category <- behav_category + sol_lnRR[,"categoryphysiology"]

			 			 category_lnRR <- rbind(cbind(mean(behav_category), HPDinterval(behav_category)), cbind(mean(develop_category), HPDinterval(develop_category)), cbind(mean(life_hist_category), HPDinterval(life_hist_category)), cbind(mean(phys_category), HPDinterval(phys_category)))
			 			 rownames(category_lnRR) <- c("behaviour", "development", "life history", "physiology")

			 	#Wild versus Lab
			 		lab_background <- as.mcmc(
			 			sol_lnRR[,"(Intercept)"] +
			 			sol_lnRR[,"matingPolygyny"] * (n_mating[2] / sum(n_mating)) +  
			 			sol_lnRR[,"matingPromiscuity"] * (n_mating[3] / sum(n_mating)) + 
			 			sol_lnRR[,"breedingsemelparous"] * (n_breeding[2] / sum(n_breeding)) + 
			 			sol_lnRR[,"categorydevelopment"] * (n_category[2] / sum(n_category)) + 
			 			sol_lnRR[,"categorylife history"] * (n_category[3] / sum(n_category)) +
			 			sol_lnRR[,"categoryphysiology"] * (n_category[4] / sum(n_category))
			 			)

			 		wild_background <- lab_background + sol_lnRR[,"background1wild"]

			 		 background_lnRR <- rbind(cbind(mean(lab_background), HPDinterval(lab_background)), cbind(mean(wild_background), HPDinterval(wild_background)))
			 		 rownames(background_lnRR) <- c("lab", "wild")

			 	# Mating system
			 		 monogamy_mating <- as.mcmc(
			 		 	sol_lnRR[,"(Intercept)"] +
			 			sol_lnRR[,"breedingsemelparous"] * (n_breeding[2] / sum(n_breeding)) + 
			 			sol_lnRR[,"categorydevelopment"] * (n_category[2] / sum(n_category)) + 
			 			sol_lnRR[,"categorylife history"] * (n_category[3] / sum(n_category)) +
			 			sol_lnRR[,"categoryphysiology"] * (n_category[4] / sum(n_category)) + 
			 			sol_lnRR[,"background1wild"] * (n_background1[2] / sum(n_background1))
			 			)

			 		 polygyny_mating <- monogamy_mating + sol_lnRR[,"matingPolygyny"]
			 		 promiscuity_mating <- monogamy_mating + sol_lnRR[,"matingPromiscuity"]

			 		 mating_lnRR <- rbind(cbind(mean(monogamy_mating), HPDinterval(monogamy_mating)), cbind(mean(polygyny_mating), HPDinterval(polygyny_mating)), cbind(mean(promiscuity_mating), HPDinterval(promiscuity_mating)))
			 		 rownames(mating_lnRR) <- c("monogamy", "polygyny", "promiscuity")

			 	# Breeding system
			 		 iteroparous_breeding <- as.mcmc(
			 		 	sol_lnRR[,"(Intercept)"] +
			 			sol_lnRR[,"categorydevelopment"] * (n_category[2] / sum(n_category)) + 
			 			sol_lnRR[,"categorylife history"] * (n_category[3] / sum(n_category)) +
			 			sol_lnRR[,"categoryphysiology"] * (n_category[4] / sum(n_category)) + 
			 			sol_lnRR[,"background1wild"] * (n_background1[2] / sum(n_background1)) +
			 			sol_lnRR[,"matingPolygyny"] * (n_mating[2] / sum(n_mating)) +  
			 			sol_lnRR[,"matingPromiscuity"] * (n_mating[3] / sum(n_mating))
			 			)

			 	     semelparous_breeding <- iteroparous_breeding + sol_lnRR[,"breedingsemelparous"]

			 	   breeding_lnRR <- rbind(cbind(mean(iteroparous_breeding), HPDinterval(iteroparous_breeding)), cbind(mean(semelparous_breeding), HPDinterval(semelparous_breeding)))
			 	   rownames(breeding_lnRR) <- c("iteroparous", "semelparous")

			 	  marg_estTab <- as.data.frame(rbind(category_lnRR, background_lnRR, mating_lnRR, breeding_lnRR))
			 	  marg_estTab$obs <- 1:nrow(marg_estTab)
			 	  marg_estTab$N <- N
			 	  colnames(marg_estTab)[1] <- c("effect")

			 	  # Calculate the percent increase of the sex with larger trait values
			 	  marg_estTab$per <- exp(marg_estTab$effect)
				  marg_estTab$per_inc <- ifelse(marg_estTab$per < 1, (1 - marg_estTab$per)*100, (marg_estTab$per-1)*100)
		
	#BEHAVIOURAL TRAITS
			
			behav <- subset(data, category == "behavior")
			behav$obs <- 1:dim(behav)[1]
			behav$behavioural.categories<- ifelse(behav$behavioural.categories == "dispersal", "activity", behav$behavioural.categories)
			behav$behavioural.categories <- as.factor(behav$behavioural.categories)
			behav$lnRR_2 <- ifelse(behav$direct == "high", behav$lnRR*(-1), behav$lnRR)

			VmatRR_behav <- VmCorMat(behav, "obs", "Dependency.behaviour")
			rownames(VmatRR_behav) <- colnames(VmatRR_behav) <- behav$obs
			corrplot(VmatRR_behav, tl.col = "black", tl.cex = 0.8, is.corr = FALSE, type = "lower", method = "color")

		# metafor
			# Model assuming independence
			mod_lnRR_behav <- rma.mv(lnRR_2 ~ behavioural.categories, V = v.lnRR, random = list(~1|study, ~1|spp_rotl, ~1|obs), data = behav) 
			
			mod_lnRR_behavSS <- rma.mv(lnRR_2 ~ behavioural.categories, V = v.lnRR, random = list(~1|study, ~1|spp_rotl, ~1|obs), R = list(obs = VmatRR_behav), data = behav) 

			AICc(mod_lnRR_behav)

		# MCMCglmm
			behav$esID <- 1:dim(behav)[1]
			Vmat_behav <- as(solve(diag(behav$v.lnRR)), "dgCMatrix")
			colnames(Vmat_behav) <- rownames(Vmat_behav) <- behav$esID

			prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002), G2 =  list(V = 1, nu = 0.002), G3 = list(V = 1, fix = 1))) 

			mod_lnRR_behav_mcmc <- MCMCglmm(lnRR ~ behavioural.categories, random = ~study + spp_rotl + esID, ginverse = list(spp_rotl = Ainv, esID = Vmat_behav), data = behav, prior = prior, nitt = 500000, thin = 100, pr = FALSE, family = "gaussian", verbose = FALSE)
			summary(mod_lnRR_behav_mcmc) 

			# Mean estimates in each level
				lnRR_behav_Sol <- mod_lnRR_behav_mcmc$Sol

				activity_behav <- lnRR_behav_Sol[,"(Intercept)"]
				aggression_behav <- activity_behav + lnRR_behav_Sol[,"behavioural.categoriesaggression"]
				bold_behav <- activity_behav + lnRR_behav_Sol[,"behavioural.categoriesboldness"]
				exlore_behav <- activity_behav + lnRR_behav_Sol[,"behavioural.categoriesexploration"]
				parent_behav <- activity_behav + lnRR_behav_Sol[,"behavioural.categoriesparenting"]
				stressCope_behav <- activity_behav + lnRR_behav_Sol[,"behavioural.categoriesstress-coping"]

			behav_lnRRTab <- as.data.frame(rbind(c(mean(activity_behav), HPDinterval(activity_behav)), c(mean(aggression_behav), HPDinterval(aggression_behav)), c(mean(bold_behav), HPDinterval(bold_behav)), c(mean(exlore_behav), HPDinterval(exlore_behav)), c(mean(parent_behav), HPDinterval(parent_behav)), c(mean(stressCope_behav), HPDinterval(stressCope_behav))))
			rownames(behav_lnRRTab) <- c("activity", "aggression", "boldness", "exploration", "parenting", "stress-coping")
			colnames(behav_lnRRTab) <- c("effect", "lower", "upper")
			behav_lnRRTab$obs <- 1:nrow(behav_lnRRTab)
			behav_lnRRTab$N <- table(behav$behavioural.categories)

	#PHYSIOLOGICAL TRAITS
			phys <- subset(data, category == "physiology")
			phys$obs <- 1:dim(phys)[1]
			phys$physiological.category <- as.factor(phys$physiological.category)

			# Model assuming independence
			mod_lnRR_phys <- rma.mv(lnRR ~ physiological.category, V = v.lnRR, random = list(~1|study, ~1|spp_rotl, ~1|obs), R = list(spp_rotl = phylo_cor), data = phys)

		#MCMCglmm
			phys$esID <- 1:dim(phys)[1]
			Vmat_phys <- as(solve(diag(phys$v.lnRR)), "dgCMatrix")
			colnames(Vmat_phys) <- rownames(Vmat_phys) <- phys$esID

			prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002), G2 =  list(V = 1, nu = 0.002), G3 = list(V = 1, fix = 1))) 

			mod_lnRR_phys_mcmc <- MCMCglmm(lnRR ~ physiological.category, random = ~study + spp_rotl + esID, ginverse = list(spp_rotl = Ainv, esID = Vmat_phys), data = phys, prior = prior, nitt = 500000, thin = 100, pr = FALSE, family = "gaussian", verbose = FALSE)
			summary(mod_lnRR_phys_mcmc) 

			# Mean estimates in each level
				lnRR_phys_Sol <- mod_lnRR_phys_mcmc$Sol

				baseline_phys <- lnRR_phys_Sol[,"(Intercept)"]
				immune_phys <- baseline_phys + lnRR_phys_Sol[,"physiological.categoryimmun"]
				other_phys <- baseline_phys + lnRR_phys_Sol[,"physiological.categoryother"]
				stressed_phys <- baseline_phys + lnRR_phys_Sol[,"physiological.categorystressed"]

			phys_lnRRTab <- as.data.frame(rbind(c(mean(baseline_phys), HPDinterval(baseline_phys)), c(mean(immune_phys), HPDinterval(immune_phys)), c(mean(other_phys), HPDinterval(other_phys)), c(mean(stressed_phys), HPDinterval(stressed_phys))))
			rownames(phys_lnRRTab) <- c("baseline", "immune", "other", "stressed")
			colnames(phys_lnRRTab) <- c("effect", "lower", "upper")
			phys_lnRRTab$obs <- 1:nrow(phys_lnRRTab)
			phys_lnRRTab$N <- table(phys$physiological.category)

	#lnCVR
	
		#metafor
			modnameCVR <- rma.mv(lnCVR_es ~ category + background1 + mating + breeding, V = VlnCVR, random = list(~1|study, ~1|spp_rotl, ~1|obs), method = "REML", R = list(spp_rotl = phylo_cor), data = data)
			AICc(modnameCVR)
			coefCVRTable1A <- round_df(data.frame(Est. = modnameCVR$b, LCI = modnameCVR$ci.lb, LCI = modnameCVR$ci.ub), digits = 3)

		#Sensitivity Analysis. Covariance matrix.
			modnameCVRDep <- rma.mv(lnCVR_es ~ category + background1 + mating + breeding, V = VlnCVR, R = list(spp_rotl = phylo_cor, obs = VmatCVR), random = list(~1|study, ~1|spp_rotl, ~1|obs), method = "REML", data = data)
			AICc(modnameCVRDep)
			coefCVRTable1B <- round_df(data.frame(Est. = modnameCVRDep$b, LCI = modnameCVRDep$ci.lb, LCI = modnameCVRDep$ci.ub), digits =3)
			
			TableS1CVR <- rbind(coefCVRTable1A, coefCVRTable1B)
			write.csv(TableS1CVR, "./pols_sex/output/tables/TableS1CVR_breed_r1.csv")

		#MCMCglmm

			prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002), G2 =  list(V = 1, nu = 0.002), G3 = list(V = 1, fix = 1))) 

			 mod_CVR_mcmc <- MCMCglmm(lnCVR_es ~ category + background1 + mating + breeding, random = ~study + spp_rotl + esID, ginverse = list(spp_rotl = Ainv, esID = Vmat), data = data, prior = prior, nitt = 500000, thin = 100, pr = FALSE, family = "gaussian", verbose = FALSE)
			 summary(mod_CVR_mcmc)

			# Lets get marginal mean estimates for each group 
			 	# extract fixed effects posterior distribution for each of the parameters estimated (contrasts)
			 	sol_lnCVR <- mod_CVR_mcmc$Sol

			 	# First we need the sample size for each of the categories as we'll produce a weighted average. We'll need all these throughout the process. 
			 		   n_category <- table(data$category)
			 		n_background1 <- table(data$background1)
			 			 n_mating <- table(data$mating)
			 		   n_breeding <- table(data$breeding)

			# Marginal estimates for each trait category, averaging across the other variables/parameters
			 		
			 	# Trait category
			 		behav_category_lnCVR <- as.mcmc(
			 			sol_lnCVR[,"(Intercept)"] +
			 			sol_lnCVR[,"background1wild"] * (n_background1[2] / sum(n_background1)) +
			 			sol_lnCVR[,"matingPolygyny"] * (n_mating[2] / sum(n_mating)) +  
			 			sol_lnCVR[,"matingPromiscuity"] * (n_mating[3] / sum(n_mating)) + 
			 			sol_lnCVR[,"breedingsemelparous"] * (n_breeding[2] / sum(n_breeding))
			 			)

			 		# Now that we have an intercept that is essentially unconditioned / averaged across all other levels we can get the marginal estimates for all other levels 
			 			  develop_category_lnCVR <- behav_category_lnCVR + sol_lnCVR[,"categorydevelopment"]
			 			life_hist_category_lnCVR <- behav_category_lnCVR + sol_lnCVR[,"categorylife history"]
			 			     phys_category_lnCVR <- behav_category_lnCVR + sol_lnCVR[,"categoryphysiology"]

			 			 category_lnCVR <- rbind(cbind(mean(behav_category_lnCVR), HPDinterval(behav_category_lnCVR)), cbind(mean(develop_category_lnCVR), HPDinterval(develop_category_lnCVR)), cbind(mean(life_hist_category_lnCVR), HPDinterval(life_hist_category_lnCVR)), cbind(mean(phys_category_lnCVR), HPDinterval(phys_category_lnCVR)))
			 			 rownames(category_lnCVR) <- c("behaviour", "development", "life history", "physiology")

			 	#Wild versus Lab
			 		lab_background_lnCVR <- as.mcmc(
			 			sol_lnCVR[,"(Intercept)"] +
			 			sol_lnCVR[,"matingPolygyny"] * (n_mating[2] / sum(n_mating)) +  
			 			sol_lnCVR[,"matingPromiscuity"] * (n_mating[3] / sum(n_mating)) + 
			 			sol_lnCVR[,"breedingsemelparous"] * (n_breeding[2] / sum(n_breeding)) + 
			 			sol_lnCVR[,"categorydevelopment"] * (n_category[2] / sum(n_category)) + 
			 			sol_lnCVR[,"categorylife history"] * (n_category[3] / sum(n_category)) +
			 			sol_lnCVR[,"categoryphysiology"] * (n_category[4] / sum(n_category))
			 			)

			 		wild_background_lnCVR <- lab_background_lnCVR + sol_lnCVR[,"background1wild"]

			 		 background_lnCVR <- rbind(cbind(mean(lab_background_lnCVR), HPDinterval(lab_background_lnCVR)), cbind(mean(wild_background_lnCVR), HPDinterval(wild_background_lnCVR)))
			 		 rownames(background_lnCVR) <- c("lab", "wild")

			 	# Mating system
			 		 monogamy_mating_lnCVR <- as.mcmc(
			 		 	sol_lnCVR[,"(Intercept)"] +
			 			sol_lnCVR[,"breedingsemelparous"] * (n_breeding[2] / sum(n_breeding)) + 
			 			sol_lnCVR[,"categorydevelopment"] * (n_category[2] / sum(n_category)) + 
			 			sol_lnCVR[,"categorylife history"] * (n_category[3] / sum(n_category)) +
			 			sol_lnCVR[,"categoryphysiology"] * (n_category[4] / sum(n_category)) + 
			 			sol_lnCVR[,"background1wild"] * (n_background1[2] / sum(n_background1))
			 			)

			 		 polygyny_mating_lnCVR <- monogamy_mating_lnCVR + sol_lnCVR[,"matingPolygyny"]
			 		 promiscuity_mating_lnCVR <- monogamy_mating_lnCVR + sol_lnCVR[,"matingPromiscuity"]

			 		 mating_lnCVR <- rbind(cbind(mean(monogamy_mating_lnCVR), HPDinterval(monogamy_mating_lnCVR)), cbind(mean(polygyny_mating_lnCVR), HPDinterval(polygyny_mating_lnCVR)), cbind(mean(promiscuity_mating_lnCVR), HPDinterval(promiscuity_mating_lnCVR)))
			 		 rownames(mating_lnCVR) <- c("monogamy", "polygyny", "promiscuity")

			 	# Breeding system
			 		 iteroparous_breeding_lnCVR <- as.mcmc(
			 		 	sol_lnCVR[,"(Intercept)"] +
			 			sol_lnCVR[,"categorydevelopment"] * (n_category[2] / sum(n_category)) + 
			 			sol_lnCVR[,"categorylife history"] * (n_category[3] / sum(n_category)) +
			 			sol_lnCVR[,"categoryphysiology"] * (n_category[4] / sum(n_category)) + 
			 			sol_lnCVR[,"background1wild"] * (n_background1[2] / sum(n_background1)) +
			 			sol_lnCVR[,"matingPolygyny"] * (n_mating[2] / sum(n_mating)) +  
			 			sol_lnCVR[,"matingPromiscuity"] * (n_mating[3] / sum(n_mating))
			 			)

			 	     semelparous_breeding_lnCVR <- iteroparous_breeding_lnCVR + sol_lnCVR[,"breedingsemelparous"]

			 	   breeding_lnCVR <- rbind(cbind(mean(iteroparous_breeding_lnCVR), HPDinterval(iteroparous_breeding_lnCVR)), cbind(mean(semelparous_breeding_lnCVR), HPDinterval(semelparous_breeding_lnCVR)))
			 	   rownames(breeding_lnCVR) <- c("iteroparous", "semelparous")

			 	  marg_estTab_CVR <- as.data.frame(rbind(category_lnCVR, background_lnCVR, mating_lnCVR, breeding_lnCVR))
			 	  marg_estTab_CVR$obs <- 1:nrow(marg_estTab_CVR)
			 	   marg_estTab_CVR$N <- N
			 	  colnames(marg_estTab_CVR)[1] <- c("effect")



# 6. Publication Bias
#----------------------------------------------------------------------------#
    # Eggers regression for lnRR. Modified version. Residuals should remove non-independence from multi-level model.
      
	       prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002), G2 =  list(V = 1, nu = 0.002), G3 = list(V = 1, fix = 1))) # Parameter expanded priors: V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000

       # Some mixing problems with phylogeny. Use species, which is pretty much the same. Actually, mixing problems with species too! Probably confound with study, but still different than metafor. Actually, I've figured this out. Turns out that when using mev argument this causes major mixing problems for species and phylogeny. This is WEIRD! So, added to ginverse, works pretty good and effective sample size for phylo goes up to normal. Not sure why mixing is compromised when using mev?? This problem above is NOT solved by using parameter expanded priors which should mix better than Inverse-Wishart.
	       modMCMCglmmRR <- MCMCglmm(lnRR_2 ~ 1, random = ~study + spp_rotl + esID, ginverse = list(spp_rotl = Ainv, esID = Vmat), data = data, prior = prior, nitt = 500000, thin = 100, pr = TRUE, family = "gaussian", verbose = FALSE)

	       I2(modMCMCglmmRR, v = data$v.lnRR, phylo = "spp_rotl")

       # Now that we have a multi-level model that estimates sources of non-independence and accounts for sampling error variance (esID), we can now extract residuals that are marginalized over the random effects. q
	       fitted <- predict(modMCMCglmmRR, marginal = ~esID) 
	       e <- data$lnRR_2 - fitted

       # Now grab the residuals from metafor to find out whether these residuals have also marginalised over the random effects in the model.
		   precisionlnRR <- 1/sqrt(data$v.lnRR)
	  
	   # Egger's regression
		   #WlnRR <- reslnRR*precisionlnRR
		   WlnRR <- e*precisionlnRR

		   eggerlnRR <- lm(WlnRR ~ precisionlnRR)
		   summary(eggerlnRR)

		   # Needed for funnel plots
		   metaResidRR <- rma(yi=e, vi = data$v.lnRR)
		   tfRR <-trimfill(metaResidRR)
	  
	# Eggers regression for lnCVR. Modified version. Residuals should remove non-independence from multi-level model.
		VmatCVR <- as(solve(diag(data$VlnCVR)), "dgCMatrix")
		colnames(VmatCVR) <- rownames(VmatCVR) <- data$esID
	 
	 # Re-fit the intercept only multi-level meta-analytic model. 
	 modMCMCglmmCVR <- MCMCglmm(lnCVR_es ~ 1, random = ~study + spp_rotl + esID, ginverse = list(spp_rotl = Ainv, esID = VmatCVR), data = data, prior = prior, nitt = 500000, thin = 100, pr = TRUE, family = "gaussian", verbose = FALSE)

	 # Extract residuals from model were predictions are conditioned on study and phylogenetic effects (BLUPS). Marginalise over the esID as this is fixed.
	 	fittedCVR <- predict(modMCMCglmmRR, marginal = ~esID)
	   reslnCVR <- data$lnCVR_es - fittedCVR

	 # Eggers regression
	   precisionlnCVR <- 1/sqrt(data$VlnCVR)
	   WlnCVR <- reslnCVR*precisionlnCVR

	   eggerlnCVR <- lm(WlnCVR ~ precisionlnCVR)
	   summary(eggerlnCVR)

	   # Needed for funnel plots. Take residuals and use as normal effect sizes (are independent)
		   metaResidCRR <- rma(yi=reslnCVR, vi = data$VlnCVR)
		   tfCVR <-trimfill(metaResidCRR)

# 7. Figures
#----------------------------------------------------------------------------#
	# Do some plotting. Funnel plots
	pdf(height = 4.519824, width = 9.938326, file = "./pols_sex/output/figures/figure3_r1.pdf")
			par(mfrow = c(1,2), mar = c(4, 5, 1, 1))
			
			#lnCVR
	   		funnel(metaResidRR, yaxis = "seinv", ylab = "Precision (1/SE)",xlab = "lnRR", pch = 21, digits = 0, las = 1, level = c(95, 99), back = "gray90")
			abline(v = 0, col = "red")
			mtext("A)", adj = -0.25, padj = 0.5)	
		
			#lnCVR
			funnel(metaResidCRR, yaxis = "seinv", ylab = "", xlab = "lnCVR", pch = 21, digits = 0, las = 1, level = c(95, 99), back = "gray90")
			abline(v = 0, col = "red")
			mtext("B)", adj = -0.25, padj = 0.5)	
	dev.off()		
	
	pdf(height = 7, width = 7, file = "./pols_sex/output/figures/figure1.pdf")
		  	par(mar = c(1,1,1,1))
		  	plot(phylo, cex = 0.85)  #type = "fan"
	dev.off()

	pdf(width = 12.678414, height = 5.506608, file = "./pols_sex/output/figures/Figure1_r1.pdf")
			par(mfrow = c(1,2),  bty = "n", mar = c(5,10,2,1))

			labels <- tolower(rownames(coefTabRR <- marg_estTab)) 
			yRef <- c(1:4, 7,8, 11:13, 16:17)
			
			#lnRR
			plot(obs~effect,  type = "n", xlim = c(-0.6, 0.6), ylim = c(0, max(yRef)+2), xlab = "lnRR", ylab = "", data = coefTabRR, yaxt='n', cex.lab = 1.5)
			abline(v = 0, lty = 2)
			
			points(yRef~coefTabRR[-30, "effect"], pch = 16) #-30 from other table
			arrows(x0=coefTabRR[-30,"effect"] , y0= yRef, x1= coefTabRR[-30,"lower"] , y1 = yRef, length = 0, angle = 90)
			arrows(x0=coefTabRR[-30,"effect"] , y0= yRef, x1= coefTabRR[-30,"upper"] , y1 = yRef, length = 0, angle = 90)
			mtext(side  = 2, labels, at = yRef, las = 1)
			mtext(side  = 2, expression(bold("A)")), at = max(yRef)+2, line = 6, las = 1, cex = 1.5, padj = -1.0)
			labRef <- c(5,9,14,18) #labRef <- c(5,9,13,20,25,31,37,42,46)
			titles <- c("Trait Type", "Study Envir.", "Mating Syst.", "Breeding")
			mtext(side  = 2, titles, font = 2, at = labRef, las = 1, cex = 1)
			text("Males 'faster'", x = -0.3, y = max(yRef)+2, cex = 1)
			text("Females 'faster'", x = +0.3, y = max(yRef)+2, cex = 1)

			#lnCVR
			coefTabCVR <- marg_estTab_CVR
			par(mar = c(5,1,2,10))
			plot(obs~effect,  type = "n", xlim = c(-0.6, 0.6), ylim = c(0, max(yRef)+2), xlab = "lnCVR", ylab = "", data = coefTabCVR, yaxt='n', cex.lab = 1.5)
			abline(v = 0, lty = 2)
						
			points(yRef~coefTabCVR[-30, "effect"], pch = 16)
			arrows(x0=coefTabCVR[-30,"effect"] , y0= yRef, x1= coefTabCVR[-30,"lower"] , y1 = yRef, length = 0, angle = 90)
			arrows(x0=coefTabCVR[-30,"effect"] , y0= yRef, x1= coefTabCVR[-30,"upper"] , y1 = yRef, length = 0, angle = 90)
			mtext(side  = 2, coefTabCVR$N, at = yRef, las = 1, adj = .5, line = 1)
			mtext(side = 2, expression(bold("N")), at = max(yRef)+2, las = 1, line = 0.5, cex = 1.2)
			mtext(side  = 2, expression(bold("B)")), at = max(yRef)+2, line = 1, las = 1, cex = 1.5, padj = -1.0)			
			text("Males high V", x = -0.3, y = max(yRef)+2, cex = 1)
			text("Females high V", x = +0.3, y = max(yRef)+2, cex = 1)
	dev.off()

	pdf(height = 7, width = 7, file = "./pols_sex/output/figures/FigureS1_r1.pdf")
		# Check out mean-variance relationships in each sex
			par(mar = c(5,5,1,1))
			plot(log(Mean_M) ~ log(SD_M), ylab = "log(Mean)", xlab = "log(SD)", data = data, col = "blue", las = 1, cex = 1.5, cex.axis = 1.5, cex.lab = 1.5)
			points(log(Mean_F) ~ log(SD_F),  data = data, col = "red", cex = 1.5)
			points(y = c(8, 7), x = c(-2, -2), col = c("blue", "red"), cex = 1.5)
			text(c("Males", "Females"), y = c(7.9, 6.9), x = c(-1.8, -1.8), adj = c(0,0))
			box()	
	dev.off()

	pdf(width=6.907489, height = 5.859031, file = "./pols_sex/output/figures/pred.Fig_r1.pdf")
		par(bty = "n", mar = c(5,10,2,1))
				predictions$yRef <- c(c(1:(0.5*(nrow(predictions)))), c(18:(nrow(predictions)+5)))
				Labels <- as.character(interaction(predictions$category, predictions$mating))
				#lnRR
				plot(yRef~pred,  type = "n", xlim = c(-0.5, 0.5), ylim = c(0, max(yRef)+2), xlab = "Predicted lnRR", ylab = "", data = predictions, yaxt='n', cex.lab = 1.5)
				abline(v = 0, lty = 2, col = "gray90")
				
				points(predictions$yRef~predictions[, "pred"], pch = 16) #-30 from other table
				arrows(x0=predictions[,"pred"] , y0= predictions$yRef, x1= predictions[,"CI_L"] , y1 = predictions$yRef, length = 0, angle = 90)
				arrows(x0=predictions[,"pred"] , y0= predictions$yRef, x1= predictions[,"CI_U"] , y1 = predictions$yRef, length = 0, angle = 90)

				text(x = 0, y = 30, "iteroparous", font = 2)
				text(x = 0, y = 14, "semelparous", font = 2)
				mtext(side  = 2, Labels, at = predictions$yRef, las = 1, cex = 0.8)
				abline(v=0, lty=2)
	dev.off()			

	pdf(width=13.736111, height = 7.111111, file = "./pols_sex/output/figures/behavPhys_fig2_r1.pdf")
			par(mfrow = c(1,2), bty = "n", mar = c(5,9,1,10))

			labels <- tolower(rownames(coefTabRR <- behav_lnRRTab)) 
			yRef <- c(1:6)
			
			#lnRR
			plot(obs~effect,  type = "n", xlim = c(-1, 1), ylim = c(0, max(yRef)+2), xlab = "lnRR", ylab = "", data = coefTabRR, yaxt='n', cex.lab = 1.5)
			abline(v = 0, lty = 2)
			
			points(yRef~coefTabRR[, "effect"], pch = 16) #-30 from other table
			arrows(x0=coefTabRR[,"effect"] , y0= yRef, x1= coefTabRR[,"lower"] , y1 = yRef, length = 0, angle = 90)
			arrows(x0=coefTabRR[,"effect"] , y0= yRef, x1= coefTabRR[,"upper"] , y1 = yRef, length = 0, angle = 90)
			mtext(side  = 2, labels, at = yRef, las = 1)
			
			labRef <- c(7.5) #labRef <- c(5,9,13,20,25,31,37,42,46)
			titles <- c("Behavioural Trait Type")
			mtext(side  = 2, titles, font = 2, at = labRef, las = 1, cex = 1)
			mtext(side  = 4, coefTabRR$N, at = yRef, las = 1, adj = .5, line = 1)
			mtext(side = 4, expression(bold("N")), at = max(yRef)+0.5, las = 1, line = 0.5, cex = 1.2)
			mtext(side  = 4, round(coefTabRR$effect, digits = 2), at = yRef, las = 1, adj = .5, line = 3)
			mtext(side  = 4, paste0("[", round(coefTabRR$lower, digits =2),",", round(coefTabRR$upper, digits = 2), "]"), at = yRef, las = 1, adj = .5, line = 6.5)
			mtext(side = 4, expression(bold("Est. [95% CI]")), at = labRef, las = 1, line = 3, cex = 1.2)
			text("Females Higher", font = 2, x = 0.5, y = labRef)
			text("Males Higher", font = 2, x = -0.5, y = labRef)
			mtext("A)", font=2, adj = -0.50, cex = 2, padj = 1)
			par(bty = "n", mar = c(5,9,1,10))

			labels <- tolower(rownames(coefTabRR <- phys_lnRRTab)) 
			yRef <- c(1:4)
			
			#lnRR
			plot(obs~effect,  type = "n", xlim = c(-1, 1), ylim = c(0, max(yRef)+1), xlab = "lnRR", ylab = "", data = coefTabRR, yaxt='n', cex.lab = 1.5)
			abline(v = 0, lty = 2)
			
			points(yRef~coefTabRR[, "effect"], pch = 16) #-30 from other table
			arrows(x0=coefTabRR[,"effect"] , y0= yRef, x1= coefTabRR[,"lower"] , y1 = yRef, length = 0, angle = 90)
			arrows(x0=coefTabRR[,"effect"] , y0= yRef, x1= coefTabRR[,"upper"] , y1 = yRef, length = 0, angle = 90)
			mtext(side  = 2, labels, at = yRef, las = 1)
			
			labRef <- c(4.7) #labRef <- c(5,9,13,20,25,31,37,42,46)
			titles <- c("Physiological Measure")
			mtext(side  = 2, titles, font = 2, at = labRef, las = 1, cex = 1)
			mtext(side  = 4, coefTabRR$N, at = yRef, las = 1, adj = .5, line = 1)
			mtext(side = 4, expression(bold("N")), at = max(yRef)+0.5, las = 1, line = 0.5, cex = 1.2)
			mtext(side  = 4, round(coefTabRR$effect, digits = 2), at = yRef, las = 1, adj = .5, line = 3)
			mtext(side  = 4, paste0("[", round(coefTabRR$lower, digits =2),",", round(coefTabRR$upper, digits = 2), "]"), at = yRef, las = 1, adj = .5, line = 6.5)
			mtext(side = 4, expression(bold("Est. [95% CI]")), at = labRef, las = 1, line = 3, cex = 1.2)
			text("Females Higher", font = 2, x = 0.5, y = labRef)
			text("Males Higher", font = 2, x = -0.5, y = labRef)
			mtext("B)", font=2, adj = -0.50, cex = 2, padj = 1)

	dev.off()

	