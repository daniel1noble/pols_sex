#----------------------------------------------------------------------------#
# Meta-Analysis of sex differences in behaviour, physiology to test POLs
# Authors: Maja Tarka, Anja Günther, Petri Niemela, Shinichi Nakagawa and Daniel Noble
# 2017-04-03 15:53:29 AEST
#----------------------------------------------------------------------------#

	# clear working space
		rm(list=ls())

	# load functions & libraries
		source("./R/func.R")
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

	# load data
		data <- read.csv("./data/POLSsexdb_merged_20170315_recat.csv", stringsAsFactors = FALSE)

# 1. Data Exploration & Processing
#----------------------------------------------------------------------------#
	# Exclude missingness. 
		data <- data[complete.cases(data$SD_M, data$SD_F, data$M_n, data$F_n),]
		dim(data)

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
		   data$mating <- ifelse(data$mating == "unkonwn", "unknown", data$mating)
		   data$mating <- ifelse(data$mating == "polygyny", "Polygyny", data$mating)

		   data$mating <- ifelse(data$mating == "promiscuity", "Promiscuity", data$mating)

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

	# Check out species list
		#resolve species synonymn for Egernia whitii. Also, Sus scrofa subspecies can just be re-named to species so it is not dropped from Phylo construction.
		data$species <- ifelse(data$species == "Egernia whitii", "Liopholis whitii", data$species)
		data$species <- ifelse(data$species == "Sus scrofa domestica", "Sus scrofa", data$species)

		  spp <- unique(data$species)
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
	
		# Check out mean-variance relationships in each sex
		plot(log(Mean_M) ~ log(SD_M), ylab = "log(X)", xlab = "log(SD)", data = data, col = "blue", las = 1)
		points(log(Mean_F) ~ log(SD_F),  data = data, col = "red")
		points(y = c(11, 12), x = c(-5, -5), col = c("blue", "red"))
		text(c("Males", "Females"), y = c(10.5, 11.5), x = c(-4.5, -4.5), adj = c(0,0))
		abline(a = 0, b = 1)

# 3. Create covariance matrix
#----------------------------------------------------------------------------#
	# Phylogenetic correlation matrix
	# List of species
	#spp <- gsub(" ", "_", spp)
	#write.csv(spp, "./output/sppList.csv")

	# Pull species from the OTL database
	#tree <- tnrs_match_names(names = spp)

	#Cerate a tree based on IDs found. Problem here were "4124117" is just not found. Function does not work when included. This maybe dude to taxonomy change. Synonym of Egernia whitii is Liopholis whitti. Try changing this. 
	#tl <- tol_induced_subtree(ott_ids=tree$ott_id) 
	#write.tree(tl, "./output/tree")

	#Note that above code is only needed once to grab taxa and generate the tree. Now we can just import the tree file. 
	#Remove ott labels on end to make sure to matches species in dataset
	phylo <- read.tree("./output/tree")
	phylo$tip.label <- gsub("_ott.+", "", phylo$tip.label)
	phylo <- makeNodeLabel(phylo) 
	is.binary.tree(phylo)

	# Compute branch lengths
	phylo_BL <- compute.brlen(phylo, method = "Grafen", power = 0.5) # Note that larger values of power produce models with lower AICc estimates, but they are all pretty similar. Estimate of phylo variance increases with smaller power estimates. But, not greatly.

	# Create phylogenetic correlation matrix. For metafor. Make sure names match with species. 
	phylo_cor <- vcv(phylo_BL, corr = TRUE)
	names <- gsub("_", " ", rownames(phylo_cor))
	rownames(phylo_cor) <- colnames(phylo_cor) <- names

	Ainv <- inverseA(phylo_BL, nodes = "ALL", scale = TRUE)$Ainv

	# Create within study dependency and test impacts with sensitivity analysis. Assume r = 0.5 to estimate the covariance between two effects.
	 VmatRR <- VmCovMat(data, "v.lnRR", "dependence")
	VmatCVR <- VmCovMat(data, "VlnCVR", "dependence")
	
# 4. Multi-level meta-analytic models (MLMA) - intercept only for heterogeneity 
#----------------------------------------------------------------------------#
	# Study and species level random effects are mostly confounded so study will probably capture most variation anyway, but worth attempting to estimate
		data$obs <- 1:nrow(data)

	# lnRR
	   modRR_int <- rma.mv(lnRR_2 ~ 1, V = v.lnRR, random = list(~1|study, ~1|species, ~1|obs), R = list(species = phylo_cor), data = data)
	 
	  #Generate heterogeneity measures and CI's. Note "species" is needed to do the correct calculations for phylogeny.
	   I2(modRR_int, v = data$v.lnRR, phylo = "species")

	# lnCVR
	   modCVR_int <- rma.mv(lnCVR_es ~ 1, V = VlnCVR , random = list(~1|study, ~1|species, ~1|obs), R = list(species = phylo_cor), data = data)

	   #Generate heterogeneity measures and CI's. Note "species" is needed to do the correct calculations for phylogeny.
	   I2(modCVR_int, v = data$VlnCVR, phylo = "species")

# 5. Multi-level meta-regression models
#----------------------------------------------------------------------------#
	# First get N for all predictors
		vars <- c("category" , "background1" , "mating" , "parenting")

		N <- c()
		for(i in 1:length(vars)){
			samp <- as.vector(t(table(data[,vars[i]])))
			   N <- c(N, samp)
		}
	
	# lnRR
		# Fun univariate models for each predictor of interest. Get the mean in each of these groups.

		coefTabRR <- c()
		for(i in 1:length(vars)){
			form <- formula(paste0("lnRR_2", "~ -1 +", vars[i]))
			modname <- rma.mv(form, V = v.lnRR, random = list(~1|study, ~1|species, ~1|obs), R = list(species = phylo_cor), data = data)
			extraction <- extract(modname)
			extraction$transf <- exp(as.vector(extraction$Est.) + 0.5*sum(modname$sigma2))
			coefTabRR <- rbind(coefTabRR, extraction)
		}

		coefTabRR$N <- N
		coefTabRR$obs <- 1:nrow(coefTabRR)

	# lnCVR
		vars <- c("category" , "background1" , "thermy" , "climate_sp" , "SSD" , "ageclass" , "mating" , "parenting" , "breeding")

		coefTabCVR <- c()
		for(i in 1:length(vars)){
			form <- formula(paste0("lnCVR_es", "~ -1 +", vars[i]))
			modname <- rma.mv(form, V = VlnCVR, random = list(~1|study, ~1|species, ~1|obs), R = list(species = phylo_cor), data = data)
			extraction <- extract(modname)
			extraction$transf <- exp(as.vector(extraction$Est.) + 0.5*sum(modname$sigma2))
			coefTabCVR <- rbind(coefTabCVR, extraction)
		}

		coefTabCVR$N <- N
		coefTabCVR$obs <- 1:nrow(coefTabCVR)
		
# 6. Marginal estimates / unconditional means for the groups. 
#----------------------------------------------
	# Run some "full" models with the relevant variables discussed. Run in lme as we can then use the effects package

	# Note the below model, however, assumes now that the residual variance is fixed! So we must remove the estimation of the observation-level variance, or ADD it into the lme fit to get the same results between metafor and lme. But also problems including nested random effects in lme: https://biostatmatt.com/archives/2718 & here: https://stat.ethz.ch/pipermail/r-help/2002-September/025067.html.

	#lnRR
		modnameRR <- rma.mv(lnRR_2 ~ category + background1 + mating + parenting, V = v.lnRR, random = list(~1|study, ~1|species), method = "REML", R = list(species = phylo_cor), data = data)
		AICc(modnameRR)
		coefRRTable1A <- round_df(data.frame(Est. = modnameRR$b, LCI = modnameRR$ci.lb, LCI = modnameRR$ci.ub), digits = 3)

		#Sensitivity Analysis. Covariance matrix.
		modnameRRDep <- rma.mv(lnRR_2 ~ category + background1 + mating + parenting, V = VmatRR, R = list(species = phylo_cor), random = list(~1|study, ~1|species), method = "REML", data = data)
		AICc(modnameRRDep)
		coefRRTable1B <- round_df(data.frame(Est. = modnameRRDep$b, LCI = modnameRRDep$ci.lb, LCI = modnameRRDep$ci.ub), digits =3)
		
		TableS1RR <- rbind(coefRRTable1A, coefRRTable1B)
		write.csv(TableS1RR, "./output/tables/TableS1RR.csv")

		# Re-fit in lme. Try a little trick: https://biostatmatt.com/archives/2718
		data2 <- data
		data2$Dummy <- factor(1)
		data2 <- groupedData(lnRR_2~1 |Dummy, data2)
		modRR <- lme(lnRR_2 ~ category + background1 + mating + parenting, random = pdBlocked(list(pdIdent(~study-1), pdIdent(~species-1), pdIdent(~obs-1))), weights = varFixed(~v.lnRR), control=lmeControl(sigma = 1), method = "REML", data = data2)
		summary(modRR)

		# Get marginal / unconditional mean estimates from the model.
		marginalRR <- marginalize(mod = modRR, vars = c("category", "background1", "mating", "parenting")) 
		margTableRR <- margTable(marginalRR)
		margTableRR$N <- N
		margTableRR$obs <- 1:nrow(margTableRR)

		# Convert estimate back to percentage difference. Note that variances are extracted from the model excluding the sampling variance. 
		var <- sum(as.numeric(unique(VarCorr(modRR)[,"Variance"]))[1:3])
		margTableRR$per <- exp(margTableRR$effect + 0.5*(var))

	#lnCVR
	
		modnameCVR <- rma.mv(lnCVR_es ~ category + background1 + mating + parenting, V = VlnCVR, random = list(~1|study, ~1|species), method = "REML", R = list(species = phylo_cor), data = data)
		AICc(modnameCVR)
		coefCVRTable1A <- round_df(data.frame(Est. = modnameCVR$b, LCI = modnameCVR$ci.lb, LCI = modnameCVR$ci.ub), digits = 3)

		#Sensitivity Analysis. Covariance matrix.
		modnameCVRDep <- rma.mv(lnCVR_es ~ category + background1 + mating + parenting, V = VmatCVR, R = list(species = phylo_cor), random = list(~1|study, ~1|species), method = "REML", data = data)
		AICc(modnameCVRDep)
		coefCVRTable1B <- round_df(data.frame(Est. = modnameCVRDep$b, LCI = modnameCVRDep$ci.lb, LCI = modnameCVRDep$ci.ub), digits =3)
		
		TableS1CVR <- rbind(coefCVRTable1A, coefCVRTable1B)
		write.csv(TableS1CVR, "./output/tables/TableS1CVR.csv")

		# Re-fit in lme. Try a little trick: https://biostatmatt.com/archives/2718
		data2 <- data
		data2$Dummy <- factor(1)
		data2 <- groupedData(lnCVR_es~1 |Dummy, data2)
		modCVR <- lme(lnCVR_es ~ category + background1 + mating + parenting, random = pdBlocked(list(pdIdent(~study-1), pdIdent(~species-1))), weights = varFixed(~VlnCVR), control=lmeControl(sigma = 1), method = "REML", data = data2)
		summary(modCVR)

		# Get marginal / unconditional estimates from the model.
		marginalCVR <- marginalize(mod = modCVR, vars = c("category", "background1", "mating", "parenting"))
		margTableCVR <- margTable(marginalCVR)
		margTableCVR$N <- N
		margTableCVR$obs <- 1:nrow(margTableCVR)

# 7. Publication Bias
#----------------------------------------------------------------------------#
    # Eggers regression for lnRR. Modified version. Residuals should remove non-independence from multi-level model.
       # Run model in MCMCglmm
       data$phylo <- gsub(" ", "_", data$species)
       colnames(data)[match("trait", colnames(data))] <- "trait2"
       data$esID <- 1:dim(data)[1]
       Vmat <- as(solve(diag(data$v.lnRR)), "dgCMatrix")
       colnames(Vmat) <- rownames(Vmat) <- data$esID

       prior = list(R = list(V = 1, nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002), G2 =  list(V = 1, nu = 0.002), G3 = list(V = 1, fix = 1))) # Parameter expanded priors: V = 1, nu = 0.002, alpha.mu = 0, alpha.V = 1000

       # Some mixing problems with phylogeny. Use species, which is pretty much the same. Actually, mixing problems with species too! Probably confound with study, but still different than metafor.

       modMCMCglmmRR <- MCMCglmm(lnRR_2 ~ 1, random = ~study + species + esID, ginverse = list(esID = Vmat), data = data, nitt = 500000, thin = 100, pr = TRUE, family = "gaussian")
       #modMCMCglmmRR <- MCMCglmm(lnRR_2 ~ 1, mev = data$v.lnRR, random = ~study + phylo, ginverse = list(phylo = Ainv), data = data, nitt = 500000, thin = 100, pr = TRUE, family = "gaussian")
       summary(modMCMCglmmRR)

       fitted <- predict(modMCMCglmmRR)

	   reslnRR <- residuals(modRR_int)
	   precisionlnRR <- 1/sqrt(data$v.lnRR)
	   WlnRR <- reslnRR*precisionlnRR

	   eggerlnRR <- lm(WlnRR ~ precisionlnRR)
	   summary(eggerlnRR)
	  
	# Eggers regression for lnCVR. Modified version. Residuals should remove non-independence from multi-level model.
	   reslnCVR <- residuals(modCVR_int)
	   precisionlnCVR <- 1/sqrt(data$VlnCVR)
	   WlnCVR <- reslnCVR*precisionlnCVR

	   eggerlnCVR <- lm(WlnCVR ~ precisionlnCVR)
	   summary(eggerlnCVR)

# 8. Figures
#----------------------------------------------------------------------------#
	# Do some plotting. Funnel plots
	pdf(height = 4.519824, width = 9.938326, file = "./output/figures/funnels.pdf")
			par(mfrow = c(1,2), mar = c(4, 5, 1, 1))
			
			#lnCVR
	   		funnel(modRR_int, yaxis = "seinv", ylab = "Precision (1/SE)",xlab = "lnRR", pch = 21, digits = 0, las = 1, level = c(95, 99), back = "gray90")
			abline(v = 0, col = "red")
			mtext("A)", adj = -0.25, padj = 0.5)	
		
			#lnCVR
			funnel(modCVR_int, yaxis = "seinv", ylab = "", xlab = "lnCVR", pch = 21, digits = 0, las = 1, level = c(95, 99), back = "gray90")
			abline(v = 0, col = "red")
			mtext("B)", adj = -0.25, padj = 0.5)	
		dev.off()		
	
		pdf(height = 7, width = 7, file = "./output/figures/phylogeny.pdf")
		  	plot(phylo, cex = 0.85)
		dev.off()

	pdf(width = 12.678414, height = 5.506608, file = "./output/figures/Figure2.pdf")
			par(mfrow = c(1,2),  bty = "n", mar = c(5,10,2,1))

			labels <- tolower(rownames(coefTabRR <- margTableRR)) #toupper converts to caps, whereas tolower converts all to lower case. Useful function.
			yRef <- c(1:4, 7,8, 11:13, 16:18)
			#Labels <- c("Behaviour", "Development","Life History", "Physiology", "Lab", "Wild", "Ectotherm", "Endotherm", "Arid", "Artificial", "Global", "Temperate", "Tropical", "Females larger", "Males larger", "No difference", "Adults", "Juveniles", "Mixed", "Unknown", "Polygyny", "Promiscuity", "Monogamy", "Unknown", "Both", "Female", "None", "Iteroparous", "Semelparous")
			#lnRR
			plot(obs~effect,  type = "n", xlim = c(-0.6, 0.6), ylim = c(0, max(yRef)+2), xlab = "lnRR", ylab = "", data = coefTabRR, yaxt='n', cex.lab = 1.5)
			abline(v = 0, lty = 2)
			
			points(yRef~coefTabRR[-30, "effect"], pch = 16) #-30 from other table
			arrows(x0=coefTabRR[-30,"effect"] , y0= yRef, x1= coefTabRR[-30,"lower"] , y1 = yRef, length = 0, angle = 90)
			arrows(x0=coefTabRR[-30,"effect"] , y0= yRef, x1= coefTabRR[-30,"upper"] , y1 = yRef, length = 0, angle = 90)
			mtext(side  = 2, labels, at = yRef, las = 1)
			mtext(side  = 2, expression(bold("A)")), at = max(yRef)+2, line = 6, las = 1, cex = 1.5, padj = -1.0)
			labRef <- c(5,9,14,19) #labRef <- c(5,9,13,20,25,31,37,42,46)
			titles <- c("Trait Type", "Lab vs. Wild", "Mating Syst.", "Parental care")
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

