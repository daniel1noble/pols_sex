#----------------------------------------------------------------------------#
# Meta-Analysis of sex differences in behaviour, physiology to test POLs
# Authors: Maja Tarka, Anja Günther, Petri Niemela, Shinichi Nakagawa and Daniel Noble
# Thu Jan 26 11:56:30 2017
#----------------------------------------------------------------------------#

	# clear working space
	rm(list=ls())

	# load functions & libraries
	source("./R/func.R")
	library(lattice)
	library(dplyr)
	library(VIM) # Helps visualise missing data.
	library(metafor)
	library(Hmisc)
	library(rotl)

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

	# How N varies across categories
		table(data$taxon) # Small sample size in some categories of taxon. 
		table(data$thermy) # Good sample sizes in each group
		table(data$category) #Decent N in each group.
		table(data$SSD)  # Good N across levels
		table(data$mating) #unknown is pretty useless, but others are OK.
		table(data$breeding) # Most ok. But obviously unknown is useless.
		table(data$parenting) # good
		table(data$climate_sp) # Mostly ok.
		table(data$background1) # I think domestic and lab are equivalent and semi-wild and wild. Could merge

	# I don't think there is enough data to explore interactions, but just have a look at a few possible ones.
		table(data$mating, data$category) # Nope. Not enough data in many levels
		interaction.plot(data$mating, data$category, data$lnRR)

		table(data$SSD, data$category) # Nope. Not enough data in many levels
		interaction.plot(data$SSD, data$category, data$lnRR)

	# Check directions
		direct <- gsub("high values =*", "", data$direction)
		direct <- gsub(" ", "", direct)
		data$direct <- direct
		data$direct <- replace(data$direct, which(is.na(data$direct)), "fast")

	# Reverse effect size values if they are slow. This is really only important for lnRR which makes specific predictions about directionality along a fast-slow continuum.
		data$lnRR_2 <- ifelse(data$direct == "slow", data$lnRR*(-1), data$lnRR)

# 2. Exploratory plotting
#----------------------------------------------------------------------------#
		
		# Box plot function to plot the various variables quickly.
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
	spp <- gsub(" ", "_", spp)
	write.csv(spp, "./output/sppList.csv")

	# Pull species from the OTL database
	tree <- tnrs_match_names(names = spp)

	#Cerate a tree based on IDs found. Problem here were "4124117" is just not found. Function does not work when included. This maybe dude to taxonomy change. Synonym of Egernia whitii is Liopholis whitti. Try changing this. 
	tl <- tol_induced_subtree(ott_ids=tree$ott_id) 
	write.tree(tl, "./output/tree")

	#Remove ott labels on end to make sure to matches species in dataset
	phylo <- read.tree("./output/tree")
	phylo$tip.label <- gsub("_ott.+", "", phylo$tip.label)
	phylo <- makeNodeLabel(phylo) 
	is.binary.tree(phylo)

	# Compute branch lengths
	phylo_BL <- compute.brlen(phylo, method = "Grafen", power = 0.5)

	# Create phylogenetic correlation matrix. For metafor. Make sure names match with species. 
	phylo_cor <- vcv(phylo_BL, corr = TRUE)
	names <- gsub("_", " ", rownames(phylo_cor))
	rownames(phylo_cor) <- colnames(phylo_cor) <- names
	
	# Or inverse covariance matrix, which will be used for MCMCglmm, Makes sure names match with species.
	Ainv <- MCMCglmm::inverseA(phylo_BL, nodes = "ALL", scale = TRUE)$Ainv
	names <- gsub("_", " ", rownames(Ainv))
	rownames(Ainv) <- colnames(Ainv) <- names

	# Create covariance matrix between dependent effect sizes.
	data$depend<-0

    studies<-unique(data$study)
    for(i in 1:length(studies)){
    	locations<-which(data$study == studies[i])
    	data$depend[locations]<-seq(1, length(locations), 1)
    }

    data$depend<-paste(data$study, data$depend, sep="_")

    # This will grab all the rows that have values above zero, then paste the study from these rows with the actual value from dependence. 
    data$depend[which(data$dependence > 0)]<-paste(data$study[which(data$dependence > 0)], data$dependence[which(data$dependence > 0)], sep="_")

    dependency<-data$depend
    shared_cov <- which(dependency%in%dependency[duplicated(dependency)]==TRUE)
    combinations <- do.call("rbind", tapply(shared_cov, data[shared_cov,"depend"], function(x) t(combn(x,2))))

# 4. Multi-level meta-analytic models (MLMA) - intercept only for heterogeneity 
#----------------------------------------------------------------------------#
	# Study and species level random effects are mostly confounded so study will probably capture most variation anyway, but worth attempting to estimate
	data$obs <- 1:nrow(data)

	# lnRR
	   modRR_int <- rma.mv(lnRR_2 ~ 1, V = v.lnRR, random = list(~1|study, ~1|species, ~1|obs), R = list(species = phylo_cor), data = data)

	   #Generate heterogeneity measures and CI's
	   I2(modRR_int, v = data$v.lnRR, phylo = "species")


	# lnCVR
	   modCVR_int <- rma.mv(lnCVR_es ~ 1, V = VlnCVR , random = list(~1|study, ~1|species, ~1|obs), R = list(species = phylo_cor), data = data)

	   #Generate heterogeneity measures and CI's
	   I2(modCVR_int, v = data$VlnCVR, phylo = "species")

# 5. Multi-level meta-regression models
#----------------------------------------------------------------------------#
	# First get N for all predictors

		N <- c()
		for(i in 1:length(vars)){
			samp <- as.vector(t(table(data[,vars[i]])))
			   N <- c(N, samp)
		}
	
	# lnRR
		# Fun univariate models for each predictor of interest. Get the mean in each of these groups.
		vars <- c("category" , "background1" , "thermy" , "climate_sp" , "SSD" , "ageclass" , "mating" , "parenting" , "breeding")

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

		data$SSD_ord <- ifelse(data$SSD == "females larger", -1, ifelse(data$SSD == "No difference", 0, 1))

		modname2 <- rma.mv(lnRR_2 ~ -1 + category*SSD_ord + background1 + parenting, V = v.lnRR, random = list(~1|study, ~1|species, ~1|obs), R = list(species = phylo_cor), data = data)

	   # compute mean and CIs in each trait category in the wild and lab.
		   predRR <- predict(modRR_category, newmods = rbind(c(0,0,0), c(1,0,0), c(0,1,0), c(0,0,1)), tau2.levels = c(0,0,0), level = 0.95)
		   predWildRR <- predict(modRR_category, newmods = rbind(c(0,0,0,1), c(1,0,0,1), c(0,1,0,1), c(0,0,1,1)), tau2.levels = c(0,0,0), level = 0.95)
		   predLabRR <- predict(modRR_category, newmods = rbind(c(0,0,0,0), c(1,0,0,0), c(0,1,0,0), c(0,0,1,0)), tau2.levels = c(0,0,0), level = 0.95)
	unique(data$parenting)
	   # Back transform means account for Jensen's inequality. Use total variance (not including sampling)
		    wildTransRR <- 1 - exp(predWildRR$pred + 0.5*0.0683)
		    LabTransRR <- 1 - exp(predLabRR$pred + 0.5*0.0683)

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

	   # compute mean and CIs in each trait category in the wild and lab.
		   predCVR <- predict(modCVR_category, newmods = rbind(c(0,0,0), c(1,0,0), c(0,1,0), c(0,0,1)), tau2.levels = c(0,0,0), level = 0.95)
		   predWildCVR <- predict(modCVR_category, newmods = rbind(c(0,0,0,1), c(1,0,0,1), c(0,1,0,1), c(0,0,1,1)), tau2.levels = c(0,0,0), level = 
		   predLabCVR <- predict(modCVR_category, newmods = rbind(c(0,0,0,0), c(1,0,0,0), c(0,1,0,0), c(0,0,1,0)), tau2.levels = c(0,0,0), level = 0.95)

		# Run some "full" models with the relevant variables discussed. Run in lme as we can then use the effects package

		library(nlme)
		library(effects)
		data$category <- as.factor(data$category)
		data$background1 <- as.factor(data$background1)
		data$mating <- as.factor(data$mating)
		data$parenting <- as.factor(data$parenting)

		mod <- lme(lnCVR_es ~ category + background1 + mating + parenting, random = list(study = ~1), weights = varFixed(~VlnCVR), control=lmeControl(sigma = 1), method = "REML", data = data)

		mod <- lmer(lnCVR_es ~ category + background1 + mating + parenting + (1|study), weights = VlnCVR, control=lmerControl(sigma = 1), method = "REML", data = data)
		summary(mod)

		mod <- lme(lnCVR_es ~ category + background1 + mating + parenting, random = pdBlocked(list(pdSymm(~study-1), pdSymm(~species-1))), weights = varFixed(~VlnCVR), control=lmeControl(sigma = 1), method = "ML", data = data)
		summary(mod)
		

		margMeans <- Effect(focal.predictors = c("category"), mod = mod, se = TRUE, confidence.level = 0.95)
		summary(margMeans)

# 6. Figures
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

	pdf(width = 16.708333, height = 8.791667, file = "./output/figures/Figure2.pdf")
			par(mfrow = c(1,2),  bty = "n", mar = c(5,10,2,1))

			Labels <- c("Behaviour", "Development","Life History", "Physiology", "Lab", "Wild", "Ectotherm", "Endotherm", "Arid", "Artificial", "Global", "Temperate", "Tropical", "Females larger", "Males larger", "No difference", "Adults", "Juveniles", "Mixed", "Unknown", "Polygyny", "Promiscuity", "Monogamy", "Unknown", "Both", "Female", "None", "Iteroparous", "Semelparous")
			#lnRR
			plot(obs~Est.,  type = "n", xlim = c(-0.5, 0.5), ylim = c(0, 46), xlab = "lnRR", ylab = "", data = coefTabRR, yaxt='n', cex.lab = 1.5)
			abline(v = 0, lty = 2)
			yRef <- c(1:4, 7,8, 11,12, 15:19, 22:24, 27:30, 33:36, 39:41, 44:45)
			
			points(yRef~coefTabRR[-30, "Est."], pch = 16)
			arrows(x0=coefTabRR[-30,"Est."] , y0= yRef, x1= coefTabRR[-30,"ci.lb"] , y1 = yRef, length = 0, angle = 90)
			arrows(x0=coefTabRR[-30,"Est."] , y0= yRef, x1= coefTabRR[-30,"ci.up"] , y1 = yRef, length = 0, angle = 90)
			mtext(side  = 2, Labels, at = yRef, las = 1)
			mtext(side  = 2, expression(bold("A)")), at = 48, line = 6, las = 1, cex = 1.5, padj = -0.50)
			labRef <- c(5,9,13,20,25,31,37,42,46)
			titles <- c("Trait Type", "Lab vs. Wild", "Thermy", "Climate zone", "SSD", "Age class","Mating Syst.", "Parental care", "Growth")
			mtext(side  = 2, titles, font = 2, at = labRef, las = 1, cex = 1)
			text("Males 'faster'", x = -0.3, y = 47, cex = 1)
			text("Females 'faster'", x = +0.3, y = 47, cex = 1)

			#lnCVR
			par(mar = c(5,1,2,1))
			plot(obs~Est.,  type = "n", xlim = c(-0.5, 0.5), ylim = c(0, 46), xlab = "lnCVR", ylab = "", data = coefTabCVR, yaxt='n', cex.lab = 1.5)
			abline(v = 0, lty = 2)
						
			points(yRef~coefTabCVR[-30, "Est."], pch = 16)
			arrows(x0=coefTabCVR[-30,"Est."] , y0= yRef, x1= coefTabCVR[-30,"ci.lb"] , y1 = yRef, length = 0, angle = 90)
			arrows(x0=coefTabCVR[-30,"Est."] , y0= yRef, x1= coefTabCVR[-30,"ci.up"] , y1 = yRef, length = 0, angle = 90)
			mtext(side  = 2, coefTabCVR$N, at = yRef, las = 1, adj = .5, line = 1)
			mtext(side = 2, expression(bold("N")), at = 47, las = 1, line = 0.5, cex = 1.2)
			mtext(side  = 2, expression(bold("B)")), at = 48, line = 1, las = 1, cex = 1.5, padj = -0.50)			
			text("Males high V", x = -0.3, y = 47, cex = 1)
			text("Females high V", x = +0.3, y = 47, cex = 1)
	dev.off()

# 7. Publication Bias
#----------------------------------------------------------------------------#
  # Eggers regression
	   reslnRR <- residuals(modRR_int)
	   precisionlnRR <- 1/sqrt(data$v.lnRR)
	   WlnRR <- reslnRR*precisionlnRR

	   eggerlnRR <- lm(WlnRR ~ precisionlnRR)
	   summary(eggerlnRR)
	   
	   # Eggers regression
	   reslnCVR <- residuals(modCVR_int)
	   precisionlnCVR <- 1/sqrt(data$VlnCVR)
	   WlnCVR <- reslnCVR*precisionlnCVR

	   eggerlnCVR <- lm(WlnCVR ~ precisionlnCVR)
	   summary(eggerlnCVR)

# 8. Multi-model inference with metafor and glmulti. Can find some details here: http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti
#----------------------------------------------------------------------------#
	
	# Model averaging with MuMIn
		updated.rma.mv <- updateable(rma.mv)
		updated.rma.mv

		modname <- updated.rma.mv(lnRR_2 ~ 1 + category + mating + background1 + parenting, V = v.lnRR, random = list(~1|study, ~1|species, ~1|obs), R = list(species = phylo_cor), data = data)

		getCall(modname)
		
		# Function updates modname so it can be used by dredge etc.
		update(modname)

		test <- dredge(modname)

		test2 <- get.model(test, )

		mod_avg <- model.avg(test, subset = cumsum(weight) <= 0.95)
		summary(mod_avg)

	# Should be 2^9 = 512 models, only considering the main effects. Only needed initially.
	mods_lnRR <- glmulti(lnRR_2 ~ category + background1 + thermy + climate_sp + SSD + ageclass + mating + parenting + breeding, vi = data$v.lnRR, random = list(~1|study), R = list(species = phylo_cor), data = data, level = 1, fitfunction = rma.mv_glmulti, crit = "aicc", confsetsize=1000)
	saveRDS(mods_lnRR, "./output/mod_avg_res/mods_lnRR")

	# Extract model-averaged table, such that all models within 3 AICc units of top model
		   mods_lnRR <- readRDS("./output/mod_avg_res/mods_lnRR")
		mod_tab_lnRR <- weightable(mods_lnRR)
		mod_tab_lnRR <- mod_tab_lnRR[mod_tab_lnRR$aicc <= min(mod_tab_lnRR$aicc) + 3,]
		mod_tab_lnRR$deltAICc <- mod_tab_lnRR$aicc - min(mod_tab_lnRR$aicc)

		#Re-calc. weights with model set
		mod_tab_lnRR$weights <- (exp(-0.5*mod_tab_lnRR$deltAICc)) / (sum(exp(-0.5*mod_tab_lnRR$deltAICc)))
		# Problems model averaging and accounting for the random effects, so will do this post full subset model selection
		coefTabRR_modAvg <-c()
		for(i in 1:nrow(mod_tab_lnRR)){
			form <- formula(as.character(mod_tab_lnRR[i,"model"]))
			modname <- rma.mv(form, V = v.lnRR, random = list(~1|study, ~1|species, ~1|obs), R = list(species = phylo_cor), data = data)
			extraction <- extract(modname)
			aicc <- AICc(modname)
			coefTabRR_modAvg <- rbind(coefTabRR_modAvg, cbind(extraction, aicc, i))
		}

		coefTabRR_modAvg$deltAICc <- coefTabRR_modAvg$aicc - min(coefTabRR_modAvg$aicc)
		coefTabRR_modAvg <- arrange(coefTabRR_modAvg, deltAICc)

		parameters <- unique(rownames(coefTabRR_modAvg))
		
	# Should be 2^9 = 512 models, only considering the main effects. Only needed once and then can be re-loaded same as lnRR
	# mods_lnCVR <- glmulti(lnCVR_es ~ category + background1 + thermy + climate_sp + SSD + ageclass +mating + parenting + breeding, vi = data$VlnCVR, random = list(~1|study), data = data, level = 1, fitfunction = rma.mv_glmulti, crit = "aicc", confsetsize=1000)
	# saveRDS(mods_lnCVR, "./output/mod_avg_res/mods_lnCVR")

	# Extract model-averaged table, such that all models within 3 AICc units of top model
		   mods_lnCVR <- readRDS("./output/mod_avg_res/mods_lnCVR")
		mod_tab_lnCVR <- weightable(mods_lnCVR)
		mod_tab_lnCVR <- mod_tab_lnCVR[mod_tab_lnCVR$aicc <= min(mod_tab_lnCVR$aicc) + 3,]
		mod_tab_lnCVR$deltAICc <- mod_tab_lnCVR$aicc - min(mod_tab_lnCVR$aicc)

		#Re-calc. weights with model set
		mod_tab_lnCVR$weights <- (exp(-0.5*mod_tab_lnCVR$deltAICc)) / (sum(exp(-0.5*mod_tab_lnCVR$deltAICc)))

		coefTabCVR_modAvg <- TopModelEst(mod_tab_lnCVR)
		# Mod averaging
		mod.avgCVR <- mod_avg(coefTabCVR_modAvg)


# 9. Marginal estimates / unconditional means for the groups
#----------------------------------------------

# Modify functions from "effects" package

x1 <- rbinom(100, 1, 0.5)
x2 <- rbinom(100, 1, 0.5)
id <- rep(rnorm(10, 0, 1), each = 10)
e <- rnorm(100,0,1)

y <- 2 + 0.5*x1 + 1.5*x2 + id + e


data <- data.frame(y, x1, x2, id = as.factor(id))
library(lme4)

mod <- lme(y ~ x1 + x2, random= ~ 1|id, data = data)
summary(mod)



################# Extra Code that is no longer needed ################

	# Explore missing data and exclude data with no SD. Not relevant anymore. No missing data.
		#aggr(data2, cex.axis = 0.50)
		#marginplot(log(data2[,c("Mean_M", "SD_M")]))
		#vars <- c("Mean_M", "SD_M", "M_n")
		#males
		#marginmatrix(log(data2[,vars]),  col = c("blue", "red", "red4",
	       "orange", "orange4"))
		#females
		#vars <- c("Mean_F", "SD_F", "F_n")
		#marginmatrix(log(data2[,vars]),  col = c("blue", "red", "red4",
	       "orange", "orange4"))


	#Hmisc::subplot(
			#plot(1 / sqrt(v.lnRR) ~ lnRR, ylab = "", xlab = "", ylim = c(0,300), data = data, cex.lab = 1.5),  x = -1.5, y = 500, size = c(1, 1)
			#)

	# Check out odd ball points; checked a few and are correct.
	 #text(data$study, y = log(data$Mean_M)+0.30, x =  log(data$SD_M), cex = 0.5)

	 #modRR_thermy <- rma.mv(lnRR ~ thermy, V = v.lnRR, random = list(~1|study, ~1|species), data = data)
	#modRR_category <- rma.mv(lnRR_2 ~ 1 + category + background1 + thermy + climate_sp + growth, V = v.lnRR, random = list(~1|study, ~1|species), data = data)

	# Do wild and lab based studies differ in how the traits between sexes are impacted?
	#modRR_category <- rma.mv(lnRR_2 ~ -1+ category , V = v.lnRR, random = list(~1|study), data = data)

	#data$TraitBackground <- interaction(data$category, data$background1)
	#modRR_category <- rma.mv(lnRR_2 ~ -1+TraitBackground, V = v.lnRR, random = list(~1|study, ~1|species), data = data)
	
	# Subset analysis
	#modRR_category_wild <- rma.mv(lnRR_2 ~ -1 + category, V = v.lnRR, random = list(~1|study, ~1|species), data = subset(data, background1 == "wild"))
	#modRR_category_lab <- rma.mv(lnRR_2 ~ -1 + category, V = v.lnRR, random = list(~1|study, ~1|species), data = subset(data, background1 == "lab"))
	
	#modRR_taxon <- rma.mv(lnRR ~ relevel(as.factor(taxon), ref = "mammal"), V = v.lnRR, random = list(~1|study, ~1|species), data = data)
	#modRR_climate <- rma.mv(lnRR ~ climate_sp, V = v.lnRR, random = list(~1|study, ~1|species), data = data)

	#modCVR_thermy <- rma.mv(lnCVR_es ~ thermy, V = VlnCVR , random = list(~1|study, ~1|species), data = data)
	#modCVR_category <- rma.mv(lnCVR_es ~ category, V = VlnCVR , random = list(~1|study, ~1|species), data = data)
	#modCVR_taxon <- rma.mv(lnCVR_es ~ -1+ taxon, V = VlnCVR , random = list(~1|study, ~1|species), data = data)


	modRR_category <- rma.mv(lnRR_2 ~ -1+ category , V = v.lnRR, random = list(~1|study), data = data)

	modCVR_category <- rma.mv(lnCVR_es ~ -1 + category, V = VlnCVR , random = list(~1|study, ~1|species), data = data)

	#modVR_thermy <- rma.mv(lnVR ~ thermy, V = var.lnVR , random = list(~1|study, ~1|species), data = data)
	#modVR_category <- rma.mv(lnVR ~ category, V = var.lnVR , random = list(~1|study, ~1|species), data = data)
	#modVR_taxon <- rma.mv(lnVR ~ taxon, V = var.lnVR , random = list(~1|study, ~1|species), data = data)


	   # Compare with MCMCglmm. Looks fine assuming that a residual variance is correctly estimates. Otherwise, point estimate is affected. Presumably because study variance, which factors into weighting is impacted.
	   #modRR_intMCMC <- MCMCglmm(lnRR_2~1, mev = data$v.lnRR, random=~study + species, ginverse = list(species = Ainv), family = "gaussian", data = data)

    # Compare with MCMCglmm. Looks fine assuming that a residual variance is correctly estimates. Otherwise, point estimate is affected. Presumably because study variance, which factors into weighting is impacted.
	  # modCVR_intMCMC <- MCMCglmm(lnCVR_es~1, mev = data$VlnCVR, random=~study + species, ginverse = list(species = Ainv), family = "gaussian", data = data)

	    # Compare with MCMCglmm. Looks fine assuming that a residual variance is correctly estimates. Otherwise, point estimate is affected. Presumably because study variance, which factors into weighting is impacted.
	   #modVR_intMCMC <- MCMCglmm(lnVR~1, mev = data$var.lnVR, random=~study + species, ginverse = list(species = Ainv), family = "gaussian", data = data)