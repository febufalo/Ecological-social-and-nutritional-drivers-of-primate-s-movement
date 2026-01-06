# FIT GLMM MODELS TO TEST IF CPs ARE LOCATED CLOSE TO POINTS OF INTERSECTION 
# ALONG BLT TRAVEL ROUTES. 
library(dplyr)
library(lme4)
library(ggbreak) #Breaks the y axis in the plot, suppressing values not needed
library(ggplot2)
library(effects)

load("data_intersection_2024_09_d05.rda")
# Scale distance categories
df_summary_cat$cat_dist_intersec_sl <- scale(df_summary_cat$cat_dist_intersec)
df_summary_cat$behav <- as.factor(df_summary_cat$behav)

######## Fit glmm models
glm1 <- glm(count ~ cat_dist_intersec_sl*area*behav, 
            data = df_summary_cat, family = "poisson")


glmm1 <- glmer(count ~ cat_dist_intersec_sl*behav + (1|area), 
               data = df_summary_cat, family = "poisson")

glmm2 <- glmer(count ~ cat_dist_intersec_sl*behav + 
                 (0+cat_dist_intersec_sl|area) + (1|area), 
               data = df_summary_cat, family = "poisson")

glmm3 <- glmer(count ~ cat_dist_intersec_sl*behav + 
                 (1+cat_dist_intersec_sl|area), 
               data = df_summary_cat, family = "poisson")

AIC(glmm1, glmm2, glmm3)  # glmm2 is best random effects structure
# use glmm2 random effects structure to explore fixed eff

#############
# using drop1() as reviewer suggested
drop1(glmm2, test = "Chisq")

############
# use set of alternative models:

glmm4 <- glmer(count ~ cat_dist_intersec_sl+behav + 
                 (0+cat_dist_intersec_sl|area) + (1|area), 
               data = df_summary_cat, family = "poisson")

glmm5 <- glmer(count ~ cat_dist_intersec_sl + 
                 (0+cat_dist_intersec_sl|area) + (1|area), 
               data = df_summary_cat, family = "poisson")

glmm6 <- glmer(count ~ 1 + (0+cat_dist_intersec_sl|area) + (1|area), 
               data = df_summary_cat, family = "poisson")

AIC(glmm2, glmm4, glmm5, glmm6)

# same result, best model, most supported by the data, is glmm2, 
# with an interaction dist*behaviour

summary(glmm2)

# btw, similar result when doing drop1(glmm2)
# best model is still with an interaction between dist*behav  

coef(glmm2)


# likellihood ratio test
anova(glmm2, glmm6, test="LRT")

### CREATING PLOTS
# Step 1: Generate predicted values from the model
df_summary_cat$predicted_count <- predict(glmm2, type = "response")

# Step 2: Create the plot with predicted values
#tiff("Figure4_2024_09_d26.tiff", width = 11000, height = 6500, res = 1000)
plot2 <- ggplot(df_summary_cat, aes(x = cat_dist_intersec_sl, y = predicted_count, color = behav)) +
  geom_point(aes(y = count), alpha = 0.5) +  # Overlay actual count points (optional)
  geom_smooth(se = FALSE, method = "loess") +  # Use loess for smoothing
  labs(title = "",
       x = "Scale(Distance to intersection)",
       y = "Predicted number of CPs") +
  scale_color_manual(values = c("FCPs" = "black", "LCPs" = "grey45", "OCPs" = "grey80")) +
  scale_fill_manual(values = c("FCPs" = "black", "LCPs" = "grey45", "OCPs" = "grey80")) +
  theme_bw() +
  #scale_y_break(c(10, 15.5)) +  # Suppress the section of y-axis between 10 and 15
  #scale_y_break(c(17, 27.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        plot.title = element_blank(), 
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_blank())
plot2
#dev.off()

#############################################################################
### Lines of code suggested by Reviewer 1 for reporting the results of the glmm:

# there seems to be an issue with the way that the output of the GLMM is 
# reported. I am adding below the lines of code that you need to run to report 
# the results of the GLMM (for a reference check Table 2 of Jang, H., Boesch, 
# C., Mundry, R., Ban, S. D., & Janmaat, K. R. (2019). Travel linearity and 
# speed of human foragers and chimpanzees during their daily search for food 
# in tropical rainforests. Scientific reports, 9(1), 11066.):

# when running a model with an interaction, you need to run the reduced version 
# of the model without the interaction when reporting the results associated 
# with predictors outside of the interactions (for details please consult Dr 
# Roger Mundry)

# First, you need to run the full model
summary(glmm2)
# Second, you calculate the Est and SE of the interaction
round(summary(glmm2)$coefficients,3)
# Third, you calculate the CI of the interaction
confint(glmm2,method="Wald")
# Fourth, you calculate the p value of the interaction
drop1(glmm2, test="Chisq")
# Then, you repeat the same procedure with the reduced model #(without the 
# interaction)#, and you report the same parameters but for the predictors 
# outside of the interaction
summary(glmm4)
round(summary(glmm4)$coefficients,3)
confint(glmm4,method="Wald")
drop1(glmm4, test="Chisq")

# Thus, in the table you will have three rows: Dist. to intersection, Kind of 
# CP and Dist. to intersection : Kind of CP (instead of having different rows 
# per category within predictors). If you would like to report differences 
# between specific categories, you can use post hoc tests.

#### Luca, here I am confused, because the code above (in lines 124 and 125) 
# still displays one line for each kind of CP different from FCPs (behavLCPs and 
# behavOCPs). The way I understand the comment from Reviewer #1, (s)he is 
# expecting us to report just one line for CPs in general. Do you have any 
# suggestion to deal with this? Does it make sense for you?

#################################################################################
#################################################################################


effects_glmm2 <- allEffects(glmm2)
summary(effects_glmm2)

# Extracting the effect estimates
effect_estimates <- effects_glmm2[[1]]$fit

# Extracting the SE values
se_values <- effects_glmm2[[1]]$se

# Extracting cat_dist_intersec_sl and behav for labeling
cat_dist_intersec_sl <- effects_glmm2[[1]]$x$cat_dist_intersec_sl
behav <- effects_glmm2[[1]]$x$behav

# Combining the data into a data frame
effects_with_se <- data.frame(
  cat_dist_intersec_sl = cat_dist_intersec_sl,
  behav = behav,
  effect_estimate = effect_estimates,
  SE = se_values
)

# Viewing the combined data
print(effects_with_se)

library(margins)
m <- margins(glmm2)
summary(m)
plot(m)

library(multcomp)
glht(glmm2)
posthoc_behav <- glht(glmm2, linfct = mcp(behav = "Tukey"))
summary(posthoc_behav) 
confint(posthoc_behav)

summary(glmm2)
confint_glmm <- confint(glmm2, method = "profile")
print(confint_glmm)

# Specify values for cat_dist_intersec_sl
values <- c(1, 5.5, 10, 14, 19)

# Get EMMs at those specific values
emmeans_interaction <- emmeans(glmm2, ~ behav | cat_dist_intersec_sl, at = list(cat_dist_intersec_sl = values))

# Define the contrasts for the interaction
contrast_interaction <- contrast(emmeans_interaction, interaction = "pairwise")

summary(contrast_interaction)


# Test model assumptions
library(DHARMa)
testDispersion(glmm2)
#tiff("FigureS4_2024_10_d22.tiff", width = 11000, height = 6500, res = 1000)
plotQQunif(simulationOutput)
#dev.off()