rm(list=ls());gc() #clear the memory

setwd("C:/Users/sean.matthews/Documents/MEM/OSU_INR_FacultyResearchAssistant/Projects/KRFP/")

pluginkdes <- read.csv("PluginOutput_All_210525.csv")

pluginkdes$Sex<-  substr(pluginkdes$ID_year,1,1)
pluginkdes$Year <-  sub('.*_',"", pluginkdes$ID_year)
pluginkdes$ID <-  sub('*_....',"", pluginkdes$ID_year)
##### summarize home range sizes #####
library(dplyr)
femaleareas <- pluginkdes %>%
  filter(Sex == "F")
# female summaries
# sample sizes
n <- unique(femaleareas$ID_year)
length(n)
# [1] 199
# unique  individuals
n <- unique(femaleareas$ID)
length(n)
# [1] 74


# annual summaries
options(digits=2)
annualfemale <- femaleareas %>%
  group_by(Year) %>%
  summarize(Individuals = n_distinct(ID),
            HR30_mean = mean(iso30areaKm),
            HR30_median = median(iso30areaKm),
            HR30_SD = sd(iso30areaKm),
            HR30_min = min(iso30areaKm),
            HR30_max = max(iso30areaKm),
            HR50_mean = mean(iso50areaKm),
            HR50_median = median(iso50areaKm),
            HR50_SD = sd(iso50areaKm),
            HR50_min = min(iso50areaKm),
            HR50_max = max(iso50areaKm),
            HR95_mean = mean(iso95areaKm),
            HR95_median = median(iso95areaKm),
            HR95_SD = sd(iso95areaKm),
            HR95_min = min(iso95areaKm),
            HR95_max = max(iso95areaKm))

write.csv(annualfemale, "AnnualHomeRange_FemaleSummaries_210525_2.csv", row.names = FALSE)

  
# mean HR size
mean(femaleareas$iso95areaKm)
#[1] 8.14
# median HR size
median(femaleareas$iso95areaKm)
# [1] 7.08
# min
min(femaleareas$iso95areaKm)
# [1] 1.41
# max
max(femaleareas$iso95areaKm)
# [1] 30.97
# sd
sd(femaleareas$iso95areaKm)
# [1] 4.52

# mean HR size
mean(femaleareas$iso90areaKm)
#[1] 6.523955
# median HR size
median(femaleareas$iso90areaKm)
# [1] 5.766389
# min
min(femaleareas$iso90areaKm)
# [1] 1.224589
# max
max(femaleareas$iso90areaKm)
# [1] 26.63327
# sd
sd(femaleareas$iso90areaKm)
# [1] 3.645277

#### 50% ####
# median HR size
median(femaleareas$iso50areaKm)
# [1] 1.44
# min
min(femaleareas$iso50areaKm)
# [1] 0.17
# max
max(femaleareas$iso50areaKm)
# [1] 9.59
# sd
sd(femaleareas$iso50areaKm)
# [1] 1.16


#### 30% ####
# median HR size
median(femaleareas$iso30areaKm)
# [1] 0.62
# min
min(femaleareas$iso30areaKm)
# [1] 0.04
# max
max(femaleareas$iso30areaKm)
# [1] 5.05
# sd
sd(femaleareas$iso30areaKm)
# [1] 0.58

#### male summaries ####
maleareas <- pluginkdes %>%
  filter(Sex == "M")
# male summaries
# sample sizes
n <- unique(maleareas$ID_year)
length(n)
# [1] 88
# unique individuals
n <- unique(maleareas$ID)
length(n)
# [1] 43
# median HR size
median(maleareas$iso90areaKm)
# [1] 19.20
# min
min(maleareas$iso90areaKm)
# [1] 5.38
# max
max(maleareas$iso90areaKm)
# [1] 113.26
# sd
sd(maleareas$iso90areaKm)
# [1] 18.33
# male summaries
# median HR size
median(maleareas$iso50areaKm)
# [1] 4.92
# min
min(maleareas$iso50areaKm)
# [1] 1.12
# max
max(maleareas$iso50areaKm)
# [1] 36.64
# sd
sd(maleareas$iso50areaKm)
# [1] 5.50
# male summaries
# median HR size
median(maleareas$iso30areaKm)
# [1] 1.92
# min
min(maleareas$iso30areaKm)
# [1] 0.43
# max
max(maleareas$iso30areaKm)
# [1] 16.86
# sd
sd(maleareas$iso30areaKm)
# [1] 2.65

options(digits=2)
annual_male <- maleareas %>%
  group_by(Year) %>%
  summarize(Individuals = n_distinct(ID),
            HR30_mean = mean(iso30areaKm),
            HR30_median = median(iso30areaKm),
            HR30_SD = sd(iso30areaKm),
            HR30_min = min(iso30areaKm),
            HR30_max = max(iso30areaKm),
            HR50_mean = mean(iso50areaKm),
            HR50_median = median(iso50areaKm),
            HR50_SD = sd(iso50areaKm),
            HR50_min = min(iso50areaKm),
            HR50_max = max(iso50areaKm),
            HR90_mean = mean(iso90areaKm),
            HR90_median = median(iso90areaKm),
            HR90_SD = sd(iso90areaKm),
            HR90_min = min(iso90areaKm),
            HR90_max = max(iso90areaKm),
            HR95_mean = mean(iso95areaKm),
            HR95_median = median(iso95areaKm),
            HR95_SD = sd(iso95areaKm),
            HR95_min = min(iso95areaKm),
            HR95_max = max(iso95areaKm))


write.csv(annual_male, "AnnualHomeRange_MaleSummaries_210525_3.csv", row.names = TRUE)

