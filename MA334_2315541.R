# clearing environment console to avoid messy look on output
if(!is.null(dev.list()))
dev.off()
rm(list = ls())
cat("\014")

# importing required packages
library(tidyverse)
library(DescTools)

# set working directory to access the csv file and assigned to a variable

setwd("D:/DHANASEKAR/Essex/Semester 1/MA 334 DA & STATISTICS WIHT R/induvidual assignment")
BD11 <-  read.csv("proportional_species_richness_NAs_removed.csv")

# Creating column name for allocated 5 taxonomic groups

BD5_names <- c('Bees', 'Hoverflies', 'Isopods', 'Ladybirds', 'Grasshoppers_._Crickets')

# finding mean value for each rows in BD5

BD5_row_mean <- rowMeans(BD11[,BD5_names])

# separating allocated 5 taxonomic groups out of 11 and ignoring others

BD5 <- BD11 %>% 
  select("Location", all_of(BD5_names), "Easting", "Northing", "dominantLandClass","ecologicalStatus", "period") %>% 
  mutate(eco_status_5=BD5_row_mean)

# Converting these three variables into factor

BD5$period <- as.factor(BD5$period) 
BD5$dominantLandClass <- as.factor(BD5$dominantLandClass)
BD5$Location <- as.factor((BD5$Location))

# filter rows with the pattern 'w' from dominant land class for sample data

BD5<- BD5 %>%
  filter(grepl("w",dominantLandClass)) # Wales

# Function to find winsorized mean

win_mean <- function(x){
  sor <- sort(x)
  win_num <- floor(0.1 * length(sor)) 
  l <- length(sor)
  for (i in 1:win_num){
    sor[i] <- sor[win_num + 1]}
  for (i in seq(l - win_num + 1, l, 1)){
    sor[i] <- sor[l - win_num]}
  return(mean(sor))
}

# ---------------------------------------------------------------------------- #
# Summary and 20% winsorized mean for BD5

stat_7 <- data.frame()
for(i in c(2:6)){
  stat_7 <- rbind(stat_7,
                 c(names(BD5)[i],
                   summary(BD5[,i]),
                   win_mean(BD5[,i])
                   ))}
colnames(stat_7) <- c("Taxi_group", "Min", "1st Q", "Median", "Mean", "3rd Q", "Max", "Twenty_p_win_mean")
stat_7$Mean <- round(as.numeric(stat_7$Mean), digits = 2)
stat_7$Twenty_p_win_mean <- round(as.numeric(stat_7$Twenty_p_win_mean), digits = 2)
View(stat_7)

# Correlation for each taxi group 

bd5 <- BD5 %>% select(c(2:6))
correlation <- round(x = cor(bd5, use = 'pairwise.complete.obs'), digits = 2)
View(correlation)

# boxplot

par(mfrow = c(1,1))
boxplot(BD5$Ladybirds~BD5$period, main = 'Boxplot for Ladybird over two periods')

# ---------------------------------------------------------------------------- #
# Hypothesis test 1 
# Null hypothesis (H0) is the two samples come from the same distribution.

eco_status_11 <- BD5%>%pull(ecologicalStatus)
eco_status_5 <- BD5%>%pull(eco_status_5)
eco_period <- BD5%>%pull(period)

par(mfrow = c(1, 2))
boxplot(eco_status_11)
boxplot(eco_status_5)

qqplot(eco_status_11, eco_status_5)
abline(0,1, col = 'red')

BD5_cdf <- ecdf(eco_status_5)
BD11_cdf <- ecdf(eco_status_11)
plot(BD11_cdf, col = 'red', main = 'Plot of BD11 CDF and BD5 CDF')
lines(BD5_cdf, col = 'green')


ks.test(eco_status_11, eco_status_5)
# Alternate hypothesis(Ha) Got solid evidence that the two samples come from different distributions.

# Hypothesis test 2
# Null hypothesis (H0) is the mean of BD5_change is equal to mean of BD11_change.

par(mfrow = c(1, 2))
BD5_split <- BD5 %>% select(Location, period, eco_status_5) %>% 
  pivot_wider(names_from = period, values_from = eco_status_5) %>% 
  mutate(BD5_change = Y00 - Y70)
BD11_split <- BD5 %>% select(Location, period, ecologicalStatus) %>% 
  pivot_wider(names_from = period, values_from = ecologicalStatus) %>% 
  mutate(BD11_change = Y00 - Y70)

BD11_change <- BD11_split %>% pull(BD11_change)
BD5_change <- BD5_split %>% pull(BD5_change)
avg <- mean(BD11_change)


boxplot(BD11_change, main = 'Boxplot for BD11 change')
boxplot(BD5_change, main = 'Boxplot for BD5 change')

t.test(BD5_change, mu = avg)
# The mean of BD5_change is approximately near to mean of BD11_change

# ---------------------------------------------------------------------------- #
# Contingency table for BD11up against BD5up

eco_chng_BD5 <- BD5_split %>% select(Location, BD5_change)
eco_chng_BD11 <-BD11_split %>% select(Location, BD11_change)
both_eco_chng <- inner_join(eco_chng_BD5, eco_chng_BD11, by = 'Location')
both_eco_chng <- both_eco_chng %>% 
  mutate(BD11up=ifelse(both_eco_chng$BD11_change>0,'UP', 'DOWN'))%>%
  mutate(BD5up=ifelse(both_eco_chng$BD5_change>0,'UP', 'DOWN'))
BD5UP <- both_eco_chng %>% pull(BD5up)
BD11UP <- both_eco_chng %>% pull(BD11up)

table_up_down <- table(BD11UP, BD5UP)

# Contingency table for independent model

g_total <- addmargins(table_up_down)
a <- round(g_total[3, 1]*g_total[1, 3]/ g_total[3,3])
b <- round(g_total[3, 1]*g_total[2, 3]/ g_total[3,3])
c <- round(g_total[3, 2]*g_total[1, 3]/ g_total[3,3])
d <- round(g_total[3, 2]*g_total[2, 3]/ g_total[3,3])
independent_table <-as.table(cbind(c(a, c), c(b, d))) 
colnames(independent_table) <- c("DOWN","UP");rownames(independent_table) <- c("DOWN","UP")

independent_table

t4 <- write.csv(independent_table, file = 'independent table.csv')

# Likelihood-ratio test

GTest(table_up_down)

# Odds-ratio
p_BD11_up <- table_up_down[1, 1]/table_up_down[1, 2]
p_BD5_up <- table_up_down[2, 1]/table_up_down[2, 2]
odds_ratio <- p_BD11_up/p_BD5_up

# or 
OddsRatio(table_up_down)

# Sensitivity
sensitivity <- table_up_down[2,2]/colSums(table_up_down)[2]

# Specificity
specificity <- table_up_down[1,1]/colSums(table_up_down)[1]

# Youden's index
youden <- sensitivity + specificity - 1

# ---------------------------------------------------------------------------- #
# Simple linear regression

BD1 <- BD11 %>% filter(grepl('w', dominantLandClass))
nrow(BD1)

l_mod <- lm(BD1$Carabids ~ BD5$eco_status_5, y = T)
summary(l_mod)
cor(l_mod$fitted.values, l_mod$y)

d1 <- data.frame(eco_status_5 = BD5$eco_status_5, Carabids = BD1$Carabids)

ggplot(d1, aes(x = eco_status_5, y = Carabids)) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  geom_abline(color = 'red')+
  geom_point()

par(mfrow = c(1, 3))
plot(BD1$Carabids ~ BD5$eco_status_5, xlab="Eco Status 5", ylab="Carabids",main = 'Plot for Eco status 5 VS Carabids')
abline(0,1, col = 'red')
abline(l_mod, col = 'green')

plot(jitter(l_mod$fitted.values),l_mod$residuals, xlab="Fitted", ylab="Residuals", main= 'Plot for Residuals')
abline(h=0,col="blue")

qqnorm(l_mod$residuals)
qqline(l_mod$residuals, col="red")

# ---------------------------------------------------------------------------- #
# Multiple linear regression
# 1 -------------------------------------------------------------------------- #

mlr_mod <- lm(BD1$Carabids~., data = BD5[c(BD5_names)], y = T)

summary(mlr_mod)
AIC(mlr_mod); cor(mlr_mod$fitted.values, mlr_mod$y) 

# 2 -------------------------------------------------------------------------- #

mlr_mod_red <- lm(BD1$Carabids~., data =
                        BD5[c('Bees', 'Hoverflies', 'Isopods', 'Grasshoppers_._Crickets')], y = T)

summary(mlr_mod_red)
AIC(mlr_mod, mlr_mod_red); cor(mlr_mod_red$fitted.values, mlr_mod_red$y)

# 3 -------------------------------------------------------------------------- #

mlr_mod_int <- lm(BD1$Carabids ~ Bees+Hoverflies+Isopods+Grasshoppers_._Crickets
                        +Bees*Isopods, data=BD5,y = T)
summary(mlr_mod_int)
AIC(mlr_mod,mlr_mod_red,mlr_mod_int); cor(mlr_mod_int$fitted.values, mlr_mod_int$y)

par(mfrow = c(1, 3))
plot(BD1$Carabids ~ mlr_mod_int$fitted.values, xlab="Fitted", ylab="Carabids", main = 'Plot for fitted values')
abline(0,1, col = 'red')

plot(jitter(mlr_mod$fitted.values), mlr_mod$residuals, xlab="Fitted", ylab="Residuals", main = 'Plot for Residuals')
abline(h=0,col="blue")

qqnorm(mlr_mod$residuals)
qqline(mlr_mod$residuals, col="red")

# 4 -------------------------------------------------------------------------- #

training <- BD11 %>% filter(period == "Y70") # training set
test <- BD11 %>% filter(period == "Y00") # test set

nrow(training); nrow(test)

tr_mod <- lm(training$Carabids~., data=training[c(BD5_names)], y = T)
summary(tr_mod)

par(mfrow = c(1, 4))
plot(tr_mod$residuals~tr_mod$fitted.values, xlab="Fitted", ylab="Residuals", main = 'Plot for Residuals') # look for unwanted pattern in residuals
abline(0,0,col="red")

qqnorm(tr_mod$residuals)
qqline(tr_mod$residuals,col="red")

predict_Y00 <- predict(tr_mod, newdata=test, type='response')

te_mod <- lm(training$Carabids~predict_Y00)
summary(te_mod)

te_mod1 <- lm(test$Carabids~predict_Y00)

AIC(tr_mod, te_mod, te_mod1)
anova(tr_mod, te_mod, te_mod1)

plot(training$Carabids ~ tr_mod$fitted.values, xlab="Fitted", ylab="Carabids", main = 'Plot for trained vlaues')
abline(0,1,col="red")

plot(test$Carabids~predict_Y00, xlab="Predicted", ylab="Carabids", main = 'Plot for predicted values')
abline(0,1,col="red")

mean((training$Carabids - tr_mod$fitted.values)^2)  # MSE on train data set 
mean((test$Carabids-predict_Y00)^2)  # MSE on test data (higher)

# ---------------------------------------------------------------------------- #
# Open Analysis

BD2 <- BD11 %>% 
  mutate(eco_status_5=BD5_row_mean) %>% 
  filter(grepl('s', dominantLandClass)&(period == 'Y70')) # Scoltland 
nrow(BD2)

open_analysis <- lm(eco_status_5 ~ dominantLandClass, data=BD2, y = T)
summary(open_analysis)
cor(open_analysis$fitted.values, open_analysis$y)

par(mfrow = c(1,1))
qqnorm(open_analysis$residuals)
qqline(open_analysis$residuals, col = 'red')
