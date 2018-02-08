# install.packages("Lahman")
# install.packages("dplyr")
install.packages("fExtremes")
# install.packages("ggplot2")
install.packages("extRemes")
# install.packagse("Lahman", "dplyr", "fExtremes", "ggplot2", "extRemes", "evd")
library(ggplot2)
library(dplyr)
library(Lahman)
library(fExtremes)
library(extRemes)

#Load data
# data(package = 'Lahman')
data("Batting")

#Data Description
#G : Games Played
#AB : Bats
#R : Run
#H : Hits()
# X2B : Doulbs hits on which the batter reached second base safely
# X3B : Triples: hits on which the batter reached third base safely
# HR : Homeruns
# RBI : Runs Batted In
# SB : Stolen Bases
# CS : Caught Stealing
# BB : Base on Balls
# SO : Strikeouts
# IBB : Intentional walks
# HBP : Hit by pitch
# SH : Sacrifice Hits
# SF : Sacrifice files
# GIDP : Grounded into double plays

Lahman::Salaries
#Data Handling
batting <- battingStats()
#add salary
batting <- merge(batting, 
                 Lahman::Salaries[,c("playerID", "yearID", "teamID", "salary")], 
                 by=c("playerID", "yearID", "teamID"), all.x=TRUE)

#Add name, age and bat hand information:
masterInfo <- Lahman::Master[, c('playerID', 'birthYear', 'birthMonth',
                         'nameLast', 'nameFirst', 'bats')]
batting <- merge(batting, masterInfo, all.x = TRUE)
batting$age <- with(batting, yearID - birthYear -
                      ifelse(birthMonth < 10, 0, 1))
#sorting
batting <- arrange(batting, playerID, yearID, stint)

#filter for eligible hitters over PA > 450 after 1900.
eligibleHitters <- batting %>% dplyr::filter(yearID >= 1900 & PA > 450)

topHitters  <- eligibleHitters %>% group_by(yearID) %>%dplyr::filter(BA == max(BA)|BA > .400)

# Create a factor variable to distinguish the .400 hitters
topHitters$ba400 <- with(topHitters, BA >= 0.400)

# Sub-data frame for the .400 hitters plus the outliers after 1950
# (averages above .380) - used to produce labels in the plot below
bignames <- rbind(subset(topHitters, ba400),
                  subset(topHitters, yearID > 1950 & BA > 0.380))
# Cut to the relevant set of variables
bignames <- subset(bignames, select = c('playerID', 'yearID', 'nameLast',
                                        'nameFirst', 'BA'))

# Positional offsets to spread out certain labels
#                     NL TC JJ TC GS TC RH GS HH RH RH BT TW TW  RC GB TG 
bignames$xoffset <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, -2, 0, 0)
bignames$yoffset <- c(0, 0, -0.003, 0, 0, 0, 0, 0, -0.004, 0, 0, 0, 0, 0, -0.003, 0, 0)  +  0.002



ggplot(topHitters, aes(x = yearID, y = BA)) +
  geom_point(aes(colour = ba400), size = 2.5) +
  geom_hline(yintercept = 0.400, size = 1) +
  geom_text(data = bignames, aes(x = yearID + xoffset, y = BA + yoffset,
                                 label = nameLast), size = 3) +
  ylim(0.330, 0.430) +
  xlab('Year') +
  scale_y_continuous('Batting average',
                     breaks = seq(0.34, 0.42, by = 0.02),
                     labels = c('.340', '.360', '.380', '.400', '.420')) +
  geom_smooth(method = "loess")


#Pitcher
Lahman::Pitching
#Data Handling
batting <- battingStats()9
#add salary
batting <- merge(batting, 
                 Lahman::Salaries[,c("playerID", "yearID", "teamID", "salary")], 
                 by=c("playerID", "yearID", "teamID"), all.x=TRUE)

#Add name, age and bat hand information:
masterInfo <- Lahman::Master[, c('playerID', 'birthYear', 'birthMonth',
                                 'nameLast', 'nameFirst', 'bats')]
batting <- merge(batting, masterInfo, all.x = TRUE)
batting$age <- with(batting, yearID - birthYear -
                      ifelse(birthMonth < 10, 0, 1))
#sorting
batting <- arrange(batting, playerID, yearID, stint)

###########
# EVD, fExtreme, extRemes
gev_fit_model <- extRemes::fevd(topHitters$BA, type = "GEV", method = "MLE")
summary(gev_fit_model)
plot(gev_fit_model)
