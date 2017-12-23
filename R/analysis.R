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
	
# 3. Create covariance matrix
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


	# Create within study dependency and test impacts with sensitivity analysis. Assume r = 0.5 to estimate the covariance between two effects.
		 #VmatRR <- VmCovMat(data, "v.lnRR", "Dependency.individual.level")
		#VmatCVR <- VmCovMat(data, "VlnCVR", "Dependency.individual.level")

	# Build a correlation matrix to test impact of dependency
		
		
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
	   		VmatCVR <- VmCorMat(data, "obs", "Dependency.individual.level")

		

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
	# Run some "full" models with the relevant variables discussed. Run in lme as we can then use the effects package
	# Run model in MCMCglmm
	       data$phylo <- gsub(" ", "_", data$species)
	       colnames(data)[match("trait", colnames(data))] <- "trait2"
	       data$esID <- 1:dim(data)[1]
	       Vmat <- as(solve(diag(data$v.lnRR)), "dgCMatrix")
	       colnames(Vmat) <- rownames(Vmat) <- data$esID

	# Note the below model, however, assumes now that the residual variance is fixed! So we must remove the estimation of the observation-level variance, or ADD it into the lme fit to get the same results between metafor and lme. But also problems including nested random effects in lme: https://biostatmatt.com/archives/2718 & here: https://stat.ethz.ch/pipermail/r-help/2002-September/025067.html.

	#lnRR
		modnameRR <- rma.mv(lnRR_2 ~ category + background1 + mating + breeding, V = v.lnRR, random = list(~1|study, ~1|spp_rotl, ~1|obs), method = "REML", R = list(spp_rotl = phylo_cor), data = data)
		AICc(modnameRR)
		coefRRTable1A <- round_df(data.frame(Est. = modnameRR$b, LCI = modnameRR$ci.lb, LCI = modnameRR$ci.ub), digits = 3)

		I2(modnameRR, v = data$v.lnRR, phylo = "spp_rotl")
		# Model was to test whether parental care was sufficiently distinct from other traits
		# modnameRR_parent <- rma.mv(lnRR_2 ~ category + background1 + mating + breeding + parenting, V = v.lnRR, random = list(~1|study, ~1|spp_rotl), method = "REML", R = list(spp_rotl = phylo_cor), data = data)

		#Sensitivity Analysis. Covariance matrix.
			modnameRRDep <- rma.mv(lnRR_2 ~ category + background1 + mating + breeding, V = VmatRR, R = list(spp_rotl = phylo_cor), random = list(~1|study, ~1|spp_rotl, ~1|obs), method = "REML", data = data)
			AICc(modnameRRDep)
			coefRRTable1B <- round_df(data.frame(Est. = modnameRRDep$b, LCI = modnameRRDep$ci.lb, LCI = modnameRRDep$ci.ub), digits =3)
		
			TableS1RR <- rbind(coefRRTable1A, coefRRTable1B)
			write.csv(TableS1RR, "./pols_sex/output/tables/TableS1RRBreeding_r1.csv")

		#PREDICTIONS for each trait category. 
			newDat <- expand.grid(list(0, unique(data$background1), unique(data$category),unique(data$mating), unique(data$breeding), 1, "M205"), stringsAsFactors = TRUE)
			colnames(newDat) <- c("lnRR_2", "background1", "category","mating", "breeding", "esID", "study")

			# Cerate new data we would like to predict for. 
				newDat.rmv <- expand.grid(list(unique(data$background1), unique(data$category),unique(data$mating), unique(data$breeding)), stringsAsFactors = TRUE)
				colnames(newDat.rmv) <- c("background1", "category","mating", "breeding")

			# Create model matrix
				X <- model.matrix(~category + background1 + mating + breeding, data = newDat.rmv)

			# Generate point estimate predictions
				newDat.rmv$pred <- X %*% coefficients(modnameRR)

				VarCorr(modnameRR)

			#modnameRR <- rma.mv(lnRR_2 ~ category + background1 + mating + breeding, V = v.lnRR, random = list(~1|study, ~1|spp_rotl), method = "REML", R = list(spp_rotl = phylo_cor), data = data)

			newDat.rmv_pred <- predict(modnameRR, newmods = newDat.rmv, addx = TRUE, tau2.levels = c("M205", "Sepsis cynipsea"))
			newDat.rmv_pred <- cbind(newDat.rmv, as.data.frame(pred))
			
			# The model for predictions. Much easier to predict with MCMCglmm objects and can integrate RE's
			prior = list(R = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000), G = list(G1 = list(V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000),  G2 = list(V = 1, fix = 1)))

			modnameRRMCMC <- MCMCglmm(lnRR_2 ~ category + background1 + mating + breeding, random = ~study + esID, ginverse = list(esID = Vmat), nitt = 1000000, thin = 1000, pr = TRUE, family = "gaussian", verbose = FALSE, data = data)
			summary(modnameRRMCMC)

		# Make predictions for trait categories in the various levels
			pred <- predict(modnameRRMCMC, newdata=newDat, marginal = ~esID, interval = "confidence")
			predictions <- cbind(newDat, pred)
			predictions <- predictions[predictions$background1 == "wild", ]

		# Re-fit in lme. Try a little trick: https://biostatmatt.com/archives/2718
			data2 <- data
			data2$Dummy <- factor(1)
			data2 <- groupedData(lnRR_2~1 |Dummy, data2)
			
			modRR <- lme(lnRR_2 ~ category + background1 + mating + breeding, random = pdBlocked(list(pdIdent(~study-1), pdIdent(~spp_rotl-1))), weights = varFixed(~v.lnRR), control=lmeControl(sigma = 1), method = "REML", data = data2)
			summary(modRR)

			# Try parenting. Try including phylogeny in the model using corBrownian and using the correlation argument in nlme. Can't quite get corBrownian to work. 
			modRR_parenting <- lme(lnRR_2 ~ category + background1 + parenting + breeding, random = pdBlocked(list(pdIdent(~study-1), pdIdent(~spp_rotl-1))), weights = varFixed(~v.lnRR), control=lmeControl(sigma = 1), method = "REML", data = data2)
			summary(modRR_parenting)
		

		# Get marginal / unconditional mean estimates from the model.
			marginalRR <- marginalize(mod = modRR, vars = c("category", "background1", "mating", "breeding")) 
			margTableRR <- margTable(marginalRR)
			margTableRR$N <- N
			margTableRR$obs <- 1:nrow(margTableRR)

		# Marginal estimates with parental care added in.
			marginalRR_parenting <- marginalize(mod = modRR_parenting, vars = c("category", "background1", "mating", "breeding", "parenting")) 
			margTableRR_parenting <- margTable(marginalRR_parenting)

		# Convert estimate back to percentage difference. Note that variances are extracted from the model excluding the sampling variance. 
			var <- sum(as.numeric(unique(VarCorr(modRR)[,"Variance"]))[1:3])
			margTableRR$per <- exp(margTableRR$effect)
			margTableRR$per_inc <- ifelse(margTableRR$per < 1, (1 - margTableRR$per)*100, (margTableRR$per-1)*100)

		# Models with parental care
			var <- sum(as.numeric(unique(VarCorr(modRR_parenting)[,"Variance"]))[1:3])
			margTableRR_parenting$per <- exp(margTableRR_parenting$effect)

		# Calculate the percent increase of the sex with larger trait values
			margTableRR_parenting$per_inc <- ifelse(margTableRR_parenting$per < 1, (1 - margTableRR_parenting$per)*100, (margTableRR_parenting$per-1)*100)

		# Adding subset analysis for behaviour and physiology to understand how specific traits within these categories respond.

		#BEHAVIOURAL TRAITS
			behav <- subset(data, category == "behavior")
			behav$obs <- 1:dim(behav)[1]
			behav$behavioural.categories <- as.factor(behav$behavioural.categories)
			behav$lnRR_2 <- ifelse(behav$direct == "high", behav$lnRR*(-1), behav$lnRR)

			VmatRR_behav <- VmCovMat(behav, "v.lnRR", "Dependency.behaviour")
			corrplot(VmatRR_behav, tl.col = "black", tl.cex = 0.8, is.corr = FALSE, type = "lower", method = "color")

		# Model assuming independence
			mod_lnRR_behav <- rma.mv(lnRR_2 ~ behavioural.categories, V = v.lnRR, random = list(~1|study, ~1|spp_rotl, ~1|obs), data = behav, method = )
			mod_lnRR_behav <- rma.mv(lnRR_2 ~ behavioural.categories, V = v.lnRR, random = list(~1|study, ~1|obs), data = behav)
			
			AICc(mod_lnRR_behav)

		# Fit in lmedata2 <- data
			data2_behav <- behav
			data2_behav$Dummy <- factor(1)
			data2_behav <- groupedData(lnRR_2 ~ 1 |Dummy, data2_behav)
			#control=lmeControl(sigma = 1)
			
			modRR_behav_lme <- lme(lnRR_2 ~ behavioural.categories, random = pdBlocked(list(pdIdent(~study-1), pdIdent(~spp_rotl-1), pdIdent(~obs-1))), control=lmeControl(sigma = 1), weights = varFixed(~v.lnRR), method = "REML", data = data2_behav)
			summary(modRR_behav_lme)

			modRR_behav_lme <- lme(lnRR_2 ~ behavioural.categories, random = ~1|study, weights = varFixed(~v.lnRR), method = "REML", data = data2_behav)
			summary(modRR_behav_lme)

			modRR_behav_lme4 <- lmer(lnRR_2 ~ behavioural.categories + (1|study) + (1|spp_rotl), data = data2_behav)
			summary(modRR_behav_lme4)

			# Marginal estimates for behavioural categories
			marginalRR_behav <- marginalize(mod = modRR_behav_lme, vars = c("behavioural.categories")) 
			marginalRR_behav$N <- rbind(table(behav$behavioural.categories, behav$background1)[,1], table(behav$behavioural.categories, behav$background1)[,2])
			marginalRR_behav <- margTable(marginalRR_behav)
			marginalRR_behav$obs <- 1:nrow(marginalRR_behav)
		
		# Model assessing the impact of dependency in that same individuals are run across assays
			mod_lnRR_behavSS <- rma.mv(lnRR ~ behavioural.categories, V = VmatRR_behav, random = list(~1|study, ~1|spp_rotl, ~1|obs), R = list(obs = ) data = behav)
			AICc(mod_lnRR_behavSS)

		#PHYSIOLOGICAL TRAITS
			phys <- subset(data, category == "physiology")
			phys$obs <- 1:dim(phys)[1]
			phys$physiological.category <- as.factor(phys$physiological.category)

			# Model assuming independence
			mod_lnRR_phys <- rma.mv(lnRR ~ physiological.category, V = v.lnRR, random = list(~1|study, ~1|spp_rotl, ~1|obs), R = list(spp_rotl = phylo_cor), data = phys)
		

	#lnCVR
	
		modnameCVR <- rma.mv(lnCVR_es ~ category + background1 + mating + breeding, V = VlnCVR, random = list(~1|study, ~1|spp_rotl), method = "REML", R = list(spp_rotl = phylo_cor), data = data)
		AICc(modnameCVR)
		coefCVRTable1A <- round_df(data.frame(Est. = modnameCVR$b, LCI = modnameCVR$ci.lb, LCI = modnameCVR$ci.ub), digits = 3)

		#Sensitivity Analysis. Covariance matrix.
			modnameCVRDep <- rma.mv(lnCVR_es ~ category + background1 + mating + breeding, V = VmatCVR, R = list(spp_rotl = phylo_cor), random = list(~1|study, ~1|spp_rotl), method = "REML", data = data)
			AICc(modnameCVRDep)
			coefCVRTable1B <- round_df(data.frame(Est. = modnameCVRDep$b, LCI = modnameCVRDep$ci.lb, LCI = modnameCVRDep$ci.ub), digits =3)
			
			TableS1CVR <- rbind(coefCVRTable1A, coefCVRTable1B)
			write.csv(TableS1CVR, "./pols_sex/output/tables/TableS1CVR_breed_r1.csv")

		# Re-fit in lme. Try a little trick: https://biostatmatt.com/archives/2718
			data2 <- data
			data2$Dummy <- factor(1)
			data2 <- groupedData(lnCVR_es~1 |Dummy, data2)
			modnameCVRMeta <- rma.mv(lnCVR_es ~ category + background1 + mating + breeding, V = VlnCVR, random = list(~1|study, ~1|spp_rotl), method = "REML", data = data2)
			modCVR <- lme(lnCVR_es ~ category + background1 + mating + breeding, random = pdBlocked(list(pdIdent(~study-1), pdIdent(~spp_rotl-1))), weights = varFixed(~VlnCVR), control=lmeControl(sigma = 1), method = "REML", data = data2)
			summary(modCVR)

		# Get marginal / unconditional estimates from the model.
			marginalCVR <- marginalize(mod = modCVR, vars = c("category", "background1", "mating", "breeding"))
			margTableCVR <- margTable(marginalCVR)
			margTableCVR$N <- N
			margTableCVR$obs <- 1:nrow(margTableCVR)

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
	
	pdf(height = 7, width = 7, file = "./output/figures/figure1.pdf")
		  	par(mar = c(1,1,1,1))
		  	plot(phylo, cex = 0.85)  #type = "fan"
	dev.off()

	pdf(width = 12.678414, height = 5.506608, file = "./pols_sex/output/figures/Figure1_r1.pdf")
			par(mfrow = c(1,2),  bty = "n", mar = c(5,10,2,1))

			labels <- tolower(rownames(coefTabRR <- margTableRR)) 
			yRef <- c(1:4, 7,8, 11:14, 17:18)
			
			#lnRR
			plot(obs~effect,  type = "n", xlim = c(-0.6, 0.6), ylim = c(0, max(yRef)+2), xlab = "lnRR", ylab = "", data = coefTabRR, yaxt='n', cex.lab = 1.5)
			abline(v = 0, lty = 2)
			
			points(yRef~coefTabRR[-30, "effect"], pch = 16) #-30 from other table
			arrows(x0=coefTabRR[-30,"effect"] , y0= yRef, x1= coefTabRR[-30,"lower"] , y1 = yRef, length = 0, angle = 90)
			arrows(x0=coefTabRR[-30,"effect"] , y0= yRef, x1= coefTabRR[-30,"upper"] , y1 = yRef, length = 0, angle = 90)
			mtext(side  = 2, labels, at = yRef, las = 1)
			mtext(side  = 2, expression(bold("A)")), at = max(yRef)+2, line = 6, las = 1, cex = 1.5, padj = -1.0)
			labRef <- c(5,9,15,19) #labRef <- c(5,9,13,20,25,31,37,42,46)
			titles <- c("Trait Type", "Study Envir.", "Mating Syst.", "Breeding")
			mtext(side  = 2, titles, font = 2, at = labRef, las = 1, cex = 1)
			text("Males 'faster'", x = -0.3, y = max(yRef)+2, cex = 1)
			text("Females 'faster'", x = +0.3, y = max(yRef)+2, cex = 1)

			#lnCVR
			coefTabCVR <- margTableCVR
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


	pdf(width=6.907489, height = 5.859031, file = "./output/figures/pred.Fig.pdf")
		par(bty = "n", mar = c(5,10,2,1))
				predictions$yRef <- c(c(1:(0.5*(nrow(predictions)))), c(18:(nrow(predictions)+1)))
				Labels <- as.character(interaction(predictions$category, predictions$mating))
				#lnRR
				plot(yRef~fit,  type = "n", xlim = c(-0.5, 0.5), ylim = c(0, max(yRef)+2), xlab = "Predicted lnRR", ylab = "", data = predictions, yaxt='n', cex.lab = 1.5)
				abline(v = 0, lty = 2, col = "gray90")
				
				points(predictions$yRef~predictions[, "fit"], pch = 16) #-30 from other table
				arrows(x0=predictions[,"fit"] , y0= predictions$yRef, x1= predictions[,"lwr"] , y1 = predictions$yRef, length = 0, angle = 90)
				arrows(x0=predictions[,"fit"] , y0= predictions$yRef, x1= predictions[,"upr"] , y1 = predictions$yRef, length = 0, angle = 90)

				text(x = 0, y = 34, "iteroparous", font = 2)
				text(x = 0, y = 17, "semelparous", font = 2)
				mtext(side  = 2, Labels, at = predictions$yRef, las = 1, cex = 0.8)
				abline(v=0, lty=2)
	dev.off()			

