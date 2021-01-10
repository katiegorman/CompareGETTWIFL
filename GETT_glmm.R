library(ggplot2)
library(lubridate)
library(dplyr)

GETT <- read.csv("GETT_fixdates.csv", header = TRUE)
head(GETT)

## get julian date
GETT$night_date.x <- as.Date(GETT$night_date.x)
GETT$jdate <- yday(GETT$night_date.x)

ggplot(GETT, aes(jdate, MYOSEP, color = site.x)) + geom_point()

## save jus t cols you need (myse calls)
t(names(GETT))
GETT <- GETT[,c(2:3,12,24,30)]

## bring in veg spreadsheet (this has total area for each habitat type,
## proporotion of area for each habitat type, and cols for open, forested, 
## and riparian %%) from https://www.nps.gov/im/vmi-gett-eise.htm
veg <- read.csv("GETT_prop_veg.csv", header = TRUE)       
head(veg)

t(names(veg))
veg <- veg[,c(1,19:21)]


## attach to call data
colnames(GETT)[names(GETT) == "site.x"] <- "site"

GETT.all <- left_join(GETT, veg, by = "site")
View(GETT.all)

GETT.all$timing <- ifelse(GETT.all$jdate < 150, "early", "late")

write.csv(GETT.all, "GETT_myse_prop_veg.csv", row.names = FALSE)


### glmms ###

library(glmmTMB)
library(pscl)
library(ggplot2)
library(corrplot)
library(dplyr)
library(egg) ## to plot two at once
library(extrafont) # to change font on plots

t(names(GETT.all))
cp <- cor(GETT.all[,c(3, 5:8)], use = "complete.obs")
cp

shapiro.test(GETT.all$MYOSEP)
## data not normally distributed, need to do nb

## null model
z1 <- zeroinfl(MYOSEP ~ 1, dist = "negbin", data = GETT.all)

## all vars alone
z2 <- zeroinfl(MYOSEP ~ timing, dist = "negbin", data = GETT.all)
z3 <- zeroinfl(MYOSEP ~ open, dist = "negbin", data = GETT.all)
z4 <- zeroinfl(MYOSEP ~ forested, dist = "negbin", data = GETT.all)
z5 <- zeroinfl(MYOSEP ~ riparian, dist = "negbin", data = GETT.all)

## global model
z5 <- zeroinfl(MYOSEP ~ timing + open + forested + riparian, dist = "negbin", data = GETT.all)

## models with one var removed
# z6 <- zeroinfl(MYOSEP ~ open + forested + riparian, dist = "negbin", data = GETT.all)
z7 <- zeroinfl(MYOSEP ~ timing + forested + riparian, dist = "negbin", data = GETT.all)
z8 <- zeroinfl(MYOSEP ~ timing + open * riparian, dist = "negbin", data = GETT.all)
z9 <- zeroinfl(MYOSEP ~ timing + open + forested, dist = "negbin", data = GETT.all)

## models with two or three vars
z10 <- zeroinfl(MYOSEP ~ forested + riparian, dist = "negbin", data = GETT.all)
z11 <- zeroinfl(MYOSEP ~ timing + open, dist = "negbin", data = GETT.all)
z12 <- zeroinfl(MYOSEP ~ timing + forested, dist = "negbin", data = GETT.all)
z13 <- zeroinfl(MYOSEP ~ timing + riparian, dist = "negbin", data = GETT.all)
z14 <- zeroinfl(MYOSEP ~ forested * riparian, dist = "negbin", data = GETT.all)


n=nrow(GETT.all)#or whatever the lenzth of your df is
table4 = AIC(z1,z2,z3,z4,z5,z9,z10,z11,z12,z13,z14)
table4$k<-c(z1$rank,z2$rank,z3$rank,z4$rank,z5$rank,
            z9$rank,z11$rank,z10$rank,z12$rank,
            z13$rank,z14$rank)
table4=table4[order(table4$AIC),]
#calculate delta AIC
table4$dAIC = table4$AIC - min(table4$AIC)
#you use the next two lines to zet weizhts
table4$edel<-exp(-0.5*table4$dAIC) 
table4$wt<-table4$edel/sum(table4$edel)
table4

summary(z10)


## get confidence intervals
confint(z11, level = 0.95)


GETT.all$predict <- predict(z10, type = "response")


f.plot <- ggplot(GETT.all, aes(forested, predict)) + stat_smooth(method = "glm") + 
  labs(title = "Gettysburg", x = "Proportion forested", y = "MYSE calls") + theme_classic() + 
  theme(text = element_text(family = "Times New Roman", face = "bold", size = 16)) + 
  theme(plot.title = element_text(hjust = 0.5))
f.plot 

r.plot <- ggplot(GETT.all, aes(riparian, predict)) + stat_smooth(method = "glm") + 
  labs(title = "Gettysburg", x = "Proportion riparian", y = "MYSE calls") + theme_classic() + 
  theme(text = element_text(family = "Times New Roman", face = "bold", size = 16)) + 
  theme(plot.title = element_text(hjust = 0.5))
r.plot

ggarrange(f.plot, r.plot, ncol = 2)

save.image("GETT_glmm.RDS")
