lnRR<-log(EMean / CMean)  # So when EMean is higher than ratio is above 0, effect size is positive and females have larger mean; when EMean is lower than ratio is below zero and effect size is negative and males have larger mean. Same is true of VR ratio.

ES<-(log(ESD) - log(EMean) + 1 / (2*(EN - 1))) - (log(CSD) - log(CMean) + 1 / (2*(CN - 1))); #When EMean is higher = females have increased variance and when EMean is smaller females have higher variance

# Marginal estimayes
# Modify functions from "effects" package

	x1 <- rbinom(100, 1, 0.5)
	x2 <- rbinom(100, 1, 0.5)
	id <- rep(rnorm(10, 0, 1), each = 10)
	e <- rnorm(100,0,1)

	y <- 2 + 0.5*x1 + 1.5*x2 + id + e
	
	summary(mod)
	mod2 <- lme(y ~ x1, random= ~ 1|id, method = "ML", data = data)
	summary(mod2)

	data <- data.frame(y, x1, x2, id = as.factor(id))
	library(lme4)

	mod <- lme(y ~ x1 + x2, random= ~ 1|id, method = "ML", data = data)
	summary(mod)
	mod2 <- lme(y ~ x1, random= ~ 1|id, method = "ML", data = data)
	summary(mod2)

	drop1(mod, test = "Chisq")


#Sexual selection, dimorphism and expect faster life in males in certain traits (most), expect corresponding variance in males. 

# Predictions
#Key findings
	#1) Males are faster developers (growth rate higher); not what you expect is females are faster in life history (caution: small sample size and generally lifespan so females); Physiology males have higher variance compared to females for physiological traits (mainly immune stuff) and somewhat development (but caution; CI's close to zero). Females have increased variance in life history 

	#2) Overall, in the wild we are getting males running faster POLs

	#3) In polygynous systems males have a faster pace of life and monogramy the opposite. Whereas when both sexes are promicuous we don't see a directional effect.

	#4) Both parental care then males exhibit.