rm(list=ls());gc() #clear the memory

#' ## Load libraries and data
#+ echo=TRUE, message=FALSE, warning=FALSE
library(ResourceSelection)
library(glmmTMB)

setwd("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/RSF")

### load data ####
alldata <- read.csv("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/RSF/KRFP_UsedAvail_Covariates_210925_2.csv")

alldata$Sex <- as.character(alldata$Sex)
alldata$Sex[alldata$Sex == "f"] <- "F"

# create a warm/cool season covariate -- corresponds to May 1 to September 30
alldata$season <- ifelse(alldata$julian > 150 & alldata$julian < 274, "Warm", "Cool")

# fix the "period" column
alldata$Period <- ifelse(alldata$FisherYear < 2012, "Pre-Drought", (ifelse(alldata$FisherYear > 2012 & alldata$FisherYear < 2015, "Drought", "Tree Mortality")))

#' Scale and center variables
alldata$scale_Bare <- scale(alldata$Bare)
alldata$scale_forest <- scale(alldata$LiveForest)
alldata$scale_LST_SummerMean <- scale(alldata$LST_SummerMean)
alldata$scale_quke_ba <- scale(alldata$quke_ba)
alldata$scale_Shrub <- scale(alldata$Shrub)
alldata$scale_TreeMortality <- scale(alldata$TreeMortality)
alldata$scale_ogsi <- scale(alldata$ogsi)
alldata$scale_StreamDist <- scale(alldata$StreamDist.tif)
alldata$scale_TPI <- scale(alldata$TPI)

#' Use and available data by animal
#with(goats, prop.table(table(ID, STATUS), 1))
with(alldata, prop.table(table(ID_year, Use), 1))

# weight the available locations
alldata$weight <- 1000^(1-alldata$Use)

# subset the data for each temporal period
library(dplyr)
library(tidyr)
predroughtdata <- filter(alldata, Period == "Pre-Drought")
droughtdata <- filter(alldata, Period == "Drought")
treemortdata <- filter(alldata, Period == "Tree Mortality")

#### predrought resource selection model ####
# here, we will not have any tree mortality data
# random slopes for individual fishers for each covariate
# other than LST, which will have random slopes for season
fisher_predrought_temp_2 <- glmmTMB(Use ~ scale_Bare + scale_forest + scale_LST_SummerMean + scale_quke_ba + scale_Shrub + 
                                    scale_ogsi +  scale_StreamDist + scale_TPI + (1|ID_year) + (0+scale_forest|ID_year) +
                               (0+scale_Bare|ID_year) + (0+scale_LST_SummerMean|season) + (0+scale_quke_ba|ID_year) + 
                               (0+scale_Shrub|ID_year) + (0+scale_ogsi|ID_year) + (0+scale_StreamDist|ID_year) + (0+scale_TPI|ID_year), 
                             family=poisson(), weights = weight, data = predroughtdata)

summary(fisher_predrought_temp_2)
ranef(fisher_predrought_temp_2)

#### drought resource selection model ####
# here, we will not have any tree mortality data
# random slopes for individual fishers for each covariate
# other than LST, which will have random slopes for season
fisher_drought_temp_2 <- glmmTMB(Use ~ scale_Bare + scale_forest + scale_LST_SummerMean + scale_quke_ba + scale_Shrub + 
                                 scale_ogsi +  scale_StreamDist + scale_TPI + (1|ID_year) + 
                                 (0+scale_Bare|ID_year) + (0+scale_forest|ID_year) + (0+scale_LST_SummerMean|season) + (0+scale_quke_ba|ID_year) + 
                                 (0+scale_Shrub|ID_year) + (0+scale_ogsi|ID_year) + (0+scale_StreamDist|ID_year) + (0+scale_TPI|ID_year), 
                               family=poisson(), data = droughtdata, weights = weight)

summary(fisher_drought_temp_2)
ranef(fisher_drought_temp_2)

##' Then fix the standard deviation of the first random term, which is the `(1|id)` component  in the above model equation. We use $\sigma=10^3$, which corresponds to a variance of $10^6$:
##+ echo=TRUE, message=FALSE,cache=TRUE
#fisher_drought_temp$parameters$theta[1] <- log(1e3)
#
##' We need to tell `glmmTMB` not to change the first entry of the vector of variances, and give all other variances another indicator to make sure they can be freely estimated:
##+ echo=TRUE, message=FALSE,cache=TRUE
#fisher_drought_temp$mapArg <- list(theta=factor(c(NA, 1)))
#
##' Then fit the model and look at the results:
##+ echo=TRUE, message=FALSE, cache=TRUE 
#fisher_drought_rsf <- glmmTMB:::fitTMB(fisher_drought_temp)
#
#summary(fisher_drought_rsf)
#ranef(fisher_drought_rsf)




#### treemort resource selection model ####
# here, we will not have any tree mortality data
# random slopes for individual fishers for each covariate
# other than LST, which will have random slopes for season
fisher_treemort_2 <- glmmTMB(Use ~ scale_Bare + scale_forest + scale_LST_SummerMean + scale_quke_ba + scale_Shrub + 
                                       scale_ogsi +  scale_StreamDist + scale_TPI + scale_TreeMortality + (1|ID_year) + (0+scale_Bare|ID_year) + 
                                       (0+scale_forest|ID_year) + (0+scale_LST_SummerMean|season) + (0+scale_quke_ba|ID_year) + 
                                       (0+scale_Shrub|ID_year) + (0+scale_ogsi|ID_year) + (0+scale_StreamDist|ID_year) + (0+scale_TPI|ID_year) + (0+scale_TreeMortality|ID_year), 
                           family=poisson(), data = treemortdata, weights = weight)
summary(fisher_treemort_2)
ranef(fisher_treemort_2)

##### extract the model output to create graphs ####



#### create summary figures  #####
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(viridis)

## individual selection w/ population level
## pre drought selection
fixedeffects <- fixef(fisher_predrought_temp_2)
fixedeffects <- fixedeffects$cond
fixedeffects2 <- as.data.frame(fixedeffects)
fixedeffects2 <- as.data.frame(fixedeffects2[-c(1),])
names(fixedeffects2) <- c("value")
fixedeffects2$name <- c("Bare", "LiveForest", "LST", "Quke", "Shrub", "OGSI", "StreamDist", "TPI")
fixedeffects2
#1    0.1552340216       Bare
#2    0.4422218460 LiveForest
#3   -0.3087895197        LST
#4    0.0048482410       Quke
#5    0.0004397144        SDI
#6   -0.0161285404      Shrub
#7   -0.2074077441 StreamDist
#8   -0.3090479752        TPI

indeffs <- coef(fisher_predrought_temp_2)
library(data.table)
indeffs_fisher2 <- as.data.frame(indeffs$cond$ID_year)
setDT(indeffs_fisher2, keep.rownames = TRUE)[]
names(indeffs_fisher2) <- c("ID_year", "Intercept", "Bare", "LiveForest", "LST", "Quke", "Shrub", "OGSI", "StreamDist", "TPI")
indeffs_fisher2$Sex <- substr(indeffs_fisher2$ID_year,1,1)
indeffs_fisher2$Sex <- ifelse(indeffs_fisher2$Sex == "F", "Female", "Male")
indeffs_fisher2$Intercept <- NULL
indeffs_fisher2$LST <- NULL
droplevels(indeffs_fisher2)
indeffs_fisher2 <- indeffs_fisher2 %>% relocate(Sex, .after = ID_year)
indeffs_fisher_long2 <- indeffs_fisher2 %>%pivot_longer(
    cols = Bare:TPI,
    names_to = "name",
    values_to = "value",
    values_drop_na = TRUE)

indeffs_fisher2_season <- as.data.frame(indeffs$cond$season)
indeffs_fisher2_season <- setDT(indeffs_fisher2_season, keep.rownames = TRUE)[]
indeffs_fisher2_season <- indeffs_fisher2_season[,-c(6:10)]
indeffs_fisher2_season <- indeffs_fisher2_season[,-c(2)]
names(indeffs_fisher2_season) <- c("ID_year", "Sex", "name", "value")
indeffs_fisher2_season$Sex <- rep(NA, 2)
indeffs_fisher2_season$name <- rep("LST", 2)

indeffs_fisher2 <- rbind(indeffs_fisher_long2, indeffs_fisher2_season)

indselection_predrought<- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = indeffs_fisher_long2, aes(x=name, y=value),color = "black", size=3, alpha=0.7, position=position_jitter(w=0.1, h=0.01)) +
  #facet_wrap(~Sex)+
  geom_point(data = fixedeffects2, aes(x=name, y=value), color = "deepskyblue2", size=5) +
  #stat_summary(fun=median, geom="point", shape=19, color="black", aes(fill=name), size=5) +   
  ylim(-1.25,0.75)+
  labs(title = "Pre-Drought", x="", y = "")+
  theme(plot.title = element_text(hjust = 0.5, size =20),
        axis.text = element_text(size=16),
        axis.title= element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"))
ggsave("Fisher_PreDroughtSelection.jpg", plot = last_plot(), dpi=300)

## drought selection
drought_fixedeffects <- fixef(fisher_drought_temp_2)
drought_fixedeffects <- drought_fixedeffects$cond
drought_fixedeffects2 <- as.data.frame(drought_fixedeffects)
drought_fixedeffects2 <- as.data.frame(drought_fixedeffects2[-c(1),])
names(drought_fixedeffects2) <- c("value")
drought_fixedeffects2$name <- c("Bare", "LiveForest", "LST", "Quke", "Shrub", "OGSI", "StreamDist", "TPI")
drought_fixedeffects2
#1        0.17404254       Bare
#2        0.61218867 LiveForest
#3       -0.25655047        LST
#4        0.03689910       Quke
#5        0.10479461        SDI
#6        0.03811888      Shrub
#7       -0.27287790 StreamDist
#8       -0.24746625        TPI

drought_indeffs <- coef(fisher_drought_temp_2)
library(data.table)
drought_indeffs_fisher2 <- as.data.frame(drought_indeffs$cond$ID_year)
setDT(drought_indeffs_fisher2, keep.rownames = TRUE)[]
names(drought_indeffs_fisher2) <- c("ID_year", "Intercept", "Bare", "LiveForest", "LST", "Quke", "Shrub", "OGSI", "StreamDist", "TPI")
drought_indeffs_fisher2$Sex <- substr(drought_indeffs_fisher2$ID_year,1,1)
drought_indeffs_fisher2$Sex <- ifelse(drought_indeffs_fisher2$Sex == "F", "Female", "Male")
drought_indeffs_fisher2$Intercept <- NULL
drought_indeffs_fisher2$LST <- NULL
droplevels(drought_indeffs_fisher2)
drought_indeffs_fisher2 <- drought_indeffs_fisher2 %>% relocate(Sex, .after = ID_year)
drought_indeffs_fisher_long2 <- drought_indeffs_fisher2 %>%pivot_longer(
  cols = Bare:TPI,
  names_to = "name",
  values_to = "value",
  values_drop_na = TRUE)

drought_indeffs_fisher2_season <- as.data.frame(drought_indeffs$cond$season)
drought_indeffs_fisher2_season <- setDT(drought_indeffs_fisher2_season, keep.rownames = TRUE)[]
drought_indeffs_fisher2_season <- drought_indeffs_fisher2_season[,-c(6:10)]
drought_indeffs_fisher2_season <- drought_indeffs_fisher2_season[,-c(2)]
names(drought_indeffs_fisher2_season) <- c("ID_year", "Sex", "name", "value")
drought_indeffs_fisher2_season$Sex <- rep(NA, 2)
drought_indeffs_fisher2_season$name <- rep("LST", 2)

drought_indeffs_fisher2 <- rbind(drought_indeffs_fisher_long2, drought_indeffs_fisher2_season)

indselection_drought <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = drought_indeffs_fisher_long2, aes(x=name, y=value), color = "black", size=3, alpha=0.7, position=position_jitter(w=0.1, h=0.01)) +
  geom_point(data = drought_fixedeffects2, aes(x=name, y=value), color = "springgreen3", size=5) +
  ylim(-1.25,0.75)+
  labs(title = "Drought", x="", y = "Selection Coefficient")+
  theme(plot.title = element_text(hjust = 0.5, size =20),
        axis.text = element_text(size=16),
        axis.title= element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"))
ggsave("Fisher_DroughtSelection.jpg", plot = last_plot(), dpi=300)


## ## treemort selection
treemort_fixedeffects <- fixef(fisher_treemort_2)
treemort_fixedeffects <- treemort_fixedeffects$cond
treemort_fixedeffects2 <- as.data.frame(treemort_fixedeffects)
treemort_fixedeffects2 <- as.data.frame(treemort_fixedeffects2[-c(1),])
names(treemort_fixedeffects2) <- c("value")
treemort_fixedeffects2$name <- c("Bare", "LiveForest", "LST", "Quke", "Shrub", "OGSI", "StreamDist", "TPI", "TreeMort")
treemort_fixedeffects2
##1 -0.25334849       Bare
##2  0.22175409 LiveForest
##3 -0.31353859        LST
##4  0.07444690       Quke
##5 -0.26552167        SDI
##6 -0.05070582      Shrub
##7 -0.19334732 StreamDist
##8 -0.21647258        TPI
##9  0.16196473   TreeMort

treemort_indeffs <- coef(fisher_treemort_2)
library(data.table)
treemort_indeffs_fisher2 <- as.data.frame(treemort_indeffs$cond$ID_year)
setDT(treemort_indeffs_fisher2, keep.rownames = TRUE)[]
names(treemort_indeffs_fisher2) <- c("ID_year", "Intercept", "Bare", "LiveForest", "LST", "Quke", "Shrub", "OGSI", "StreamDist", "TPI", "TreeMort")
treemort_indeffs_fisher2$Sex <- substr(treemort_indeffs_fisher2$ID_year,1,1)
treemort_indeffs_fisher2$Sex <- ifelse(treemort_indeffs_fisher2$Sex == "F", "Female", "Male")
treemort_indeffs_fisher2$Intercept <- NULL
treemort_indeffs_fisher2$LST <- NULL
droplevels(treemort_indeffs_fisher2)
treemort_indeffs_fisher2 <- treemort_indeffs_fisher2 %>% relocate(Sex, .after = ID_year)
treemort_indeffs_fisher_long2 <- treemort_indeffs_fisher2 %>%pivot_longer(
  cols = Bare:TreeMort,
  names_to = "name",
  values_to = "value",
  values_drop_na = TRUE)

treemort_indeffs_fisher2_season <- as.data.frame(treemort_indeffs$cond$season)
treemort_indeffs_fisher2_season <- setDT(treemort_indeffs_fisher2_season, keep.rownames = TRUE)[]
treemort_indeffs_fisher2_season <- treemort_indeffs_fisher2_season[,-c(6:11)]
treemort_indeffs_fisher2_season <- treemort_indeffs_fisher2_season[,-c(2)]
names(treemort_indeffs_fisher2_season) <- c("ID_year", "Sex", "name", "value")
treemort_indeffs_fisher2_season$Sex <- rep(NA, 2)
treemort_indeffs_fisher2_season$name <- rep("LST", 2)

treemort_indeffs_fisher2 <- rbind(treemort_indeffs_fisher_long2, treemort_indeffs_fisher2_season)

indselection_treemort <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = treemort_indeffs_fisher_long2, aes(x=name, y=value), size=3, alpha=0.7, color = "black", position=position_jitter(w=0.1, h=0.01)) +
  ylim(-1.25,0.75)+
  geom_point(data = treemort_fixedeffects2, aes(x=name, y=value), color = "gold", size=5) +
  labs(title = "Tree Mortality" , x="Covariate", y = "")+
  theme(plot.title = element_text(hjust = 0.5, size =20),
        axis.text = element_text(size=16),
        axis.title= element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"))

png(filename = "KRFP_IndividualVariation_210926.png", width = 11, height = 15, units = "in", res = 300)
ggarrange(
  indselection_predrought,
  indselection_drought,
  indselection_treemort,
  ncol = 1,
  nrow = 3
)
dev.off()


#ggsave("Fisher_TreeMortSelection.jpg", plot = last_plot(), dpi=300)

#### combine the individual level effects ##
drought_indeffs_fisher_long2$Period <- rep("Drought", length(drought_indeffs_fisher_long2$ID_year))
indeffs_fisher_long2$Period <- rep("Pre-Drought", length(indeffs_fisher_long2$ID_year))
treemort_indeffs_fisher_long2$Period <- rep("Tree Mortality", length(treemort_indeffs_fisher_long2$ID_year))

combined_indeffs <- rbind(indeffs_fisher_long2, drought_indeffs_fisher_long2, treemort_indeffs_fisher_long2)
combined_indeffs_wide <- combined_indeffs %>% pivot_wider(names_from = name, values_from = value)


#### population level selection and variance for random effects #
# predrought
fixed_predrought <- as.data.frame(confint(fisher_predrought_temp_2, level = 0.90))
setDT(fixed_predrought, keep.rownames = TRUE)[]
fixed_predrought <- fixed_predrought[-c(1, 10:18), ]
names(fixed_predrought) <- c("Covariate", "Lower90", "Upper90", "Estimate")
fixed_predrought$Covariate<- c("Bare", "LiveForest", "LST", "Quke", "Shrub", "OGSI", "StreamDist", "TPI")

# predrought variance
predrought_indeffs <- combined_indeffs_wide[combined_indeffs_wide$Period == "Pre-Drought",]
predrought_indeffs$Barevar <- (predrought_indeffs$Bare-mean(predrought_indeffs$Bare))^2
predrought_indeffs$Forestvar <- (predrought_indeffs$LiveForest-mean(predrought_indeffs$LiveForest))^2
predrought_indeffs$Qukevar <- (predrought_indeffs$Quke-mean(predrought_indeffs$Quke))^2
predrought_indeffs$OGSIvar <- (predrought_indeffs$OGSI- mean(predrought_indeffs$OGSI))^2
predrought_indeffs$Shrubvar <- (predrought_indeffs$Shrub-mean(predrought_indeffs$Shrub))^2
predrought_indeffs$StreamDistvar <- (predrought_indeffs$StreamDist- mean(predrought_indeffs$StreamDist))^2
predrought_indeffs$TPIvar <- (predrought_indeffs$TPI-mean(predrought_indeffs$TPI))^2
indeffs_fisher2_season$LSTVar <- (indeffs_fisher2_season$value- mean(indeffs_fisher2_season$value))^2


BareVar <- as.data.frame(sum(predrought_indeffs$Barevar))
names(BareVar) <- c("variance")
ForestVar <- as.data.frame(sum(predrought_indeffs$Forestvar))
names(ForestVar) <- c("variance")
LSTVar <- as.data.frame(sum(indeffs_fisher2_season$LSTVar))
names(LSTVar) <- c("variance")
QukeVar <- as.data.frame(sum(predrought_indeffs$Qukevar))
names(QukeVar) <- c("variance")
OGSIVar <- as.data.frame(sum(predrought_indeffs$OGSIvar))
names(OGSIVar) <- c("variance")
ShrubVar <- as.data.frame(sum(predrought_indeffs$Shrubvar))
names(ShrubVar) <- c("variance")
StreamDistVar <- as.data.frame(sum(predrought_indeffs$StreamDistvar))
names(StreamDistVar) <- c("variance")
TPIVar <- as.data.frame(sum(predrought_indeffs$TPIvar))
names(TPIVar) <- c("variance")

predrought_var <- rbind(BareVar,
                        ForestVar,
                        LSTVar,
                        QukeVar,
                        OGSIVar,
                        ShrubVar,
                        StreamDistVar,
                        TPIVar)
predrought_var$Covariate <-  c("Bare", "LiveForest", "LST", "Quke", "OGSI", "Shrub", "StreamDist", "TPI")
View(predrought_var)


fixed_predrought$Covariate = forcats::fct_rev(factor(fixed_predrought$Covariate))
predrought_var$Covariate = forcats::fct_rev(factor(predrought_var$Covariate))
# fixed effect predrought plot
predrought_fixedplot <- ggplot(data = fixed_predrought, aes(x=Estimate, y=Covariate, fill = Covariate, color = Covariate)) +
  geom_point(size = 5)+
  geom_vline(xintercept = 0, linetype = "dashed") +
  #geom_point(, size=5) +
  geom_pointrange(data = fixed_predrought, aes(xmin=Lower90, xmax=Upper90), size =1)+
  xlim(-0.45,0.75)+
  scale_color_viridis(discrete = TRUE, option = "D") +
  scale_fill_viridis(discrete = TRUE) +
  labs(title = "Pre-Drought", x="", y = "")+
  theme(plot.title = element_text(hjust = 1.0, size =20),
        axis.text = element_text(size=16),
        axis.title= element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"))

# variance predrought plot
predrought_indvar <- ggplot(predrought_var, aes(x = variance, y = Covariate, fill = Covariate, color = Covariate))+
  scale_color_viridis(discrete = TRUE, option = "D") +
  scale_fill_viridis(discrete = TRUE) +
  geom_col(width = 0.5) +
  xlim(0,20)+
  labs(title = "", x = "", y = "") +
  theme(plot.title = element_text(hjust = 0.5, size =20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    axis.line = element_line(colour = "black"))

fixedeffects_indvariance_plot <- ggarrange(predrought_fixedplot,predrought_indvar, widths = c(1.5, 1))

# drought
fixed_drought <- as.data.frame(confint(fisher_drought_temp_2, level = 0.90))
setDT(fixed_drought, keep.rownames = TRUE)[]
fixed_drought <- fixed_drought[-c(1, 10:18), ]
names(fixed_drought) <- c("Covariate", "Lower90", "Upper90", "Estimate")
fixed_drought$Covariate<- c("Bare", "LiveForest", "LST", "Quke", "Shrub", "OGSI", "StreamDist", "TPI")

# drought variance
drought_indeffs <- combined_indeffs_wide[combined_indeffs_wide$Period == "Drought",]
drought_indeffs$Barevar <- (drought_indeffs$Bare-mean(drought_indeffs$Bare))^2
drought_indeffs$Forestvar <- (drought_indeffs$LiveForest-mean(drought_indeffs$LiveForest))^2
drought_indeffs$Qukevar <- (drought_indeffs$Quke-mean(drought_indeffs$Quke))^2
drought_indeffs$OGSIvar <- (drought_indeffs$OGSI-mean(drought_indeffs$OGSI))^2
drought_indeffs$Shrubvar <- (drought_indeffs$Shrub-mean(drought_indeffs$Shrub))^2
drought_indeffs$StreamDistvar <- (drought_indeffs$StreamDist-mean(drought_indeffs$StreamDist))^2
drought_indeffs$TPIvar <- (drought_indeffs$TPI-mean(drought_indeffs$TPI))^2
indeffs_fisher2_season$LSTVar <- (indeffs_fisher2_season$value-mean(indeffs_fisher2_season$value))^2


BareVar <- as.data.frame(sum(drought_indeffs$Barevar))
names(BareVar) <- c("variance")
ForestVar <- as.data.frame(sum(drought_indeffs$Forestvar))
names(ForestVar) <- c("variance")
LSTVar <- as.data.frame(sum(indeffs_fisher2_season$LSTVar))
names(LSTVar) <- c("variance")
QukeVar <- as.data.frame(sum(drought_indeffs$Qukevar))
names(QukeVar) <- c("variance")
OGSIVar <- as.data.frame(sum(drought_indeffs$OGSIvar))
names(OGSIVar) <- c("variance")
ShrubVar <- as.data.frame(sum(drought_indeffs$Shrubvar))
names(ShrubVar) <- c("variance")
StreamDistVar <- as.data.frame(sum(drought_indeffs$StreamDistvar))
names(StreamDistVar) <- c("variance")
TPIVar <- as.data.frame(sum(drought_indeffs$TPIvar))
names(TPIVar) <- c("variance")

drought_var <- rbind(BareVar,
                        ForestVar,
                        LSTVar,
                        QukeVar,
                        OGSIVar,
                        ShrubVar,
                        StreamDistVar,
                        TPIVar)
drought_var$Covariate <-  c("Bare", "LiveForest", "LST", "Quke", "OGSI", "Shrub", "StreamDist", "TPI")
View(drought_var)



fixed_drought$Covariate = forcats::fct_rev(factor(fixed_drought$Covariate))
drought_var$Covariate = forcats::fct_rev(factor(drought_var$Covariate))
# fixed effect drought plot
drought_fixedplot <- ggplot(data = fixed_drought, aes(x=Estimate, y=Covariate, fill = Covariate, color = Covariate)) +
  geom_point(size = 5)+
  geom_vline(xintercept = 0, linetype = "dashed") +
  #geom_point(, size=5) +
  geom_pointrange(data = fixed_drought, aes(xmin=Lower90, xmax=Upper90), size =1)+
  xlim(-0.45,0.75)+
  scale_color_viridis(discrete = TRUE, option = "D") +
  scale_fill_viridis(discrete = TRUE) +
  labs(title = "Drought", x="", y = "Covariate")+
  theme(plot.title = element_text(hjust = 1.0, size =20),
        axis.text = element_text(size=16),
        axis.title= element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"))

# variance drought plot
drought_indvar <- ggplot(drought_var, aes(x = variance, y = Covariate, fill = Covariate, color = Covariate))+
  scale_color_viridis(discrete = TRUE, option = "D") +
  scale_fill_viridis(discrete = TRUE) +
  geom_col(width = 0.5) +
  xlim(0,20)+
  labs(title = "",x = "", y = "") +
  theme(plot.title = element_text(hjust = 0.5, size =20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    axis.line = element_line(colour = "black"))

# tree mort
fixed_treemort <- as.data.frame(confint(fisher_treemort_2, level = 0.90))
setDT(fixed_treemort, keep.rownames = TRUE)[]
fixed_treemort <- fixed_treemort[-c(1, 11:20), ]
names(fixed_treemort) <- c("Covariate", "Lower90", "Upper90", "Estimate")
fixed_treemort$Covariate<- c("Bare", "LiveForest", "LST", "Quke", "Shrub", "OGSI", "StreamDist", "TPI", "TreeMort")

# treemort variance
treemort_indeffs <- combined_indeffs_wide[combined_indeffs_wide$Period == "Tree Mortality",]
treemort_indeffs$Barevar <- (treemort_indeffs$Bare- mean(treemort_indeffs$Bare))^2
treemort_indeffs$Forestvar <- (treemort_indeffs$LiveForest-mean(treemort_indeffs$LiveForest))^2
treemort_indeffs$Qukevar <- (treemort_indeffs$Quke-mean(treemort_indeffs$Quke))^2
treemort_indeffs$OGSIvar <- (treemort_indeffs$OGSI- mean(treemort_indeffs$OGSI))^2
treemort_indeffs$Shrubvar <- (treemort_indeffs$Shrub- mean(treemort_indeffs$Shrub))^2
treemort_indeffs$StreamDistvar <- (treemort_indeffs$StreamDist- mean(treemort_indeffs$StreamDist))^2
treemort_indeffs$TPIvar <- (treemort_indeffs$TPI- mean(treemort_indeffs$TPI))^2
treemort_indeffs$TreeMortvar <- (treemort_indeffs$TreeMort- mean(treemort_indeffs$TreeMort))^2
indeffs_fisher2_season$LSTVar <- (indeffs_fisher2_season$value- mean(indeffs_fisher2_season$LSTVar))^2


BareVar <- as.data.frame(sum(treemort_indeffs$Barevar))
names(BareVar) <- c("variance")
ForestVar <- as.data.frame(sum(treemort_indeffs$Forestvar))
names(ForestVar) <- c("variance")
LSTVar <- as.data.frame(sum(indeffs_fisher2_season$LSTVar))
names(LSTVar) <- c("variance")
QukeVar <- as.data.frame(sum(treemort_indeffs$Qukevar))
names(QukeVar) <- c("variance")
OGSIVar <- as.data.frame(sum(treemort_indeffs$OGSIvar))
names(OGSIVar) <- c("variance")
ShrubVar <- as.data.frame(sum(treemort_indeffs$Shrubvar))
names(ShrubVar) <- c("variance")
StreamDistVar <- as.data.frame(sum(treemort_indeffs$StreamDistvar))
names(StreamDistVar) <- c("variance")
TPIVar <- as.data.frame(sum(treemort_indeffs$TPIvar))
names(TPIVar) <- c("variance")
TreeMortVar <- as.data.frame(sum(treemort_indeffs$TreeMortvar))
names(TreeMortVar) <- c("variance")

treemort_var <- rbind(BareVar,
                     ForestVar,
                     LSTVar,
                     QukeVar,
                     OGSIVar,
                     ShrubVar,
                     StreamDistVar,
                     TPIVar,
                     TreeMortVar)
treemort_var$Covariate <-  c("Bare", "LiveForest", "LST", "Quke", "OGSI", "Shrub", "StreamDist", "TPI", "TreeMort")
View(treemort_var)



fixed_treemort$Covariate = forcats::fct_rev(factor(fixed_treemort$Covariate))
treemort_var$Covariate = forcats::fct_rev(factor(treemort_var$Covariate))
# fixed effect treemort plot
treemort_fixedplot <- ggplot(data = fixed_treemort, aes(x=Estimate, y=Covariate, fill = Covariate, color = Covariate)) +
  geom_point(size = 5)+
  geom_vline(xintercept = 0, linetype = "dashed") +
  #geom_point(, size=5) +
  geom_pointrange(data = fixed_treemort, aes(xmin=Lower90, xmax=Upper90), size =1)+
  xlim(-0.45,0.75)+
  scale_color_viridis(discrete = TRUE, option = "D") +
  scale_fill_viridis(discrete = TRUE) +
  labs(title = "Tree Mortality", x="Population-level Selection", y = "")+
  theme(plot.title = element_text(hjust = 1.0, size =20),
        axis.text = element_text(size=16),
        axis.title= element_text(size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        axis.line = element_line(colour = "black"))

# variance treemort plot
treemort_indvar <- ggplot(treemort_var, aes(x = variance, y = Covariate, fill = Covariate, color = Covariate))+
  scale_color_viridis(discrete = TRUE, option = "D") +
  scale_fill_viridis(discrete = TRUE) +
  xlim(0,20)+
  geom_col(width = 0.5) +
  labs(title = "", x = "Individual-level Variance", y = "") +
  theme(plot.title = element_text(hjust = 0.5, size =20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_blank(),
    axis.title = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    axis.line = element_line(colour = "black"))

png(filename = "KRFP_PopSelectionIndVar_210926.png", width = 8.5, height = 9, units = "in", res = 300)
ggarrange(
  predrought_fixedplot,
  predrought_indvar,
  drought_fixedplot,
  drought_indvar,
  treemort_fixedplot,
  treemort_indvar,
  widths = c(1.5, 1),
  ncol = 2,
  nrow = 3
)
dev.off()


#### resurce selection summary table ####
# summarize by sex and time-period
options(digits=2)
RSF_summary <- combined_indeffs_wide %>%
  group_by(Period) %>%
  summarize(Individuals = n_distinct(ID_year),
            Baremed = median(Bare, na.rm = TRUE),
            Baremin = min(Bare, na.rm = TRUE),
            Baremax = max(Bare, na.rm = TRUE),
            LiveForestmed = median(LiveForest, na.rm = TRUE),
            LiveForestmin = min(LiveForest, na.rm = TRUE),
            LiveForestmax = max(LiveForest, na.rm = TRUE),
            OGSImed = median(OGSI, na.rm = TRUE),
            OGSImin = min(OGSI, na.rm = TRUE),
            OGSImax = max(OGSI, na.rm = TRUE),
            Qukemed = median(Quke, na.rm = TRUE),
            Qukemin = min(Quke, na.rm = TRUE),
            Qukemax = max(Quke, na.rm = TRUE),
            Shrubmed = median(Shrub, na.rm = TRUE),
            Shrubmin = min(Shrub, na.rm = TRUE),
            Shrubmax = max(Shrub, na.rm = TRUE),
            StreamDistmed = median(StreamDist, na.rm = TRUE),
            StreamDistmin = min(StreamDist, na.rm = TRUE),
            StreamDistmax = max(StreamDist, na.rm = TRUE),
            TPImed = median(TPI, na.rm = TRUE),
            TPImin = min(TPI, na.rm = TRUE),
            TPImax = max(TPI, na.rm = TRUE),
            TreeMortmed = median(TreeMort, na.rm = TRUE),
            TreeMortmin = min(TreeMort, na.rm = TRUE),
            TreeMortmax = max(TreeMort, na.rm = TRUE))

RSF_summary$Bare <- paste0(round(RSF_summary$Baremed, digits = 2)," (",round(RSF_summary$Baremin, digits = 2),", ", round(RSF_summary$Baremax, digits = 2),")")
RSF_summary$LiveForest <- paste0(round(RSF_summary$LiveForestmed, digits = 2)," (",round(RSF_summary$LiveForestmin, digits = 2),", ", round(RSF_summary$LiveForestmax, digits = 2),")")
RSF_summary$OGSI <- paste0(round(RSF_summary$OGSImed, digits = 2)," (",round(RSF_summary$OGSImin, digits = 2),", ", round(RSF_summary$OGSImax, digits = 2),")")
RSF_summary$Quke <- paste0(round(RSF_summary$Qukemed, digits = 2)," (",round(RSF_summary$Qukemin, digits = 2),", ", round(RSF_summary$Qukemax, digits = 2),")")
RSF_summary$Shrub <- paste0(round(RSF_summary$Shrubmed, digits = 2)," (",round(RSF_summary$Shrubmin, digits = 2),", ", round(RSF_summary$Shrubmax, digits = 2),")")
RSF_summary$StreamDist <- paste0(round(RSF_summary$StreamDistmed, digits = 2)," (",round(RSF_summary$StreamDistmin, digits = 2),", ", round(RSF_summary$StreamDistmax, digits = 2),")")
RSF_summary$TPI <- paste0(round(RSF_summary$TPImed, digits = 2)," (",round(RSF_summary$TPImin, digits = 2),", ", round(RSF_summary$TPImax, digits = 2),")")
RSF_summary$TreeMort <- paste0(round(RSF_summary$TreeMortmed, digits = 2)," (",round(RSF_summary$TreeMortmin, digits = 2),", ", round(RSF_summary$TreeMortmax, digits = 2),")")

RSF_summary2 <- RSF_summary[,-c(4:27)]
#TPI = concatenate(median(TPI), "(",sd(TPI),")"),


write.csv(RSF_summary2, "AnnualHomeRange_RSFSummaries_210926.csv", row.names = FALSE)

##### save fixed effects
fixed_predrought$Period <- rep("Pre-Drought")
fixed_drought$Period <- rep("Drought")
fixed_treemort$Period <- rep("Tree Mortality")

fixed_summary <- rbind(fixed_predrought,fixed_drought,fixed_treemort)

write.csv(fixed_summary, "AnnualHomeRange_RSFSummaries_FixedEffects_210926.csv", row.names = FALSE)
