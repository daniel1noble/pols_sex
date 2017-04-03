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



