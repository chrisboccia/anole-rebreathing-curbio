###This script takes the BORIS summary csv produced by boris_py.py
###and generates individual scores (best performance and )

###use line below if using RStudio to setwd to folder containing script####
if(Sys.getenv("RSTUDIO") == "1")
	setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##otherwise, setwd to suit your placement of the files##
#setwd("Your path")


####load packages####
library(geiger)
library(ape)
library(phytools)
library(ggplot2)
library(nlme)
library(dplyr)
library(car)

####get functions####

##loads 'matchstick' and 'var_endings' custom functions
source("functions.R")

####Rebreathing behavioural data#####

###RAW ethogram data is processed using a python script (see readme markdown for behaviour pipeline)###
###After processing, the resulting data is matched with individual metadata using the function below###
###Data are then filtered (juveniles and trials below 15s duration are removed)
###If you are not re-processing the raw data, you do not need to run the following matching script
###to generate the appropriate data structure for downstream analysis. Instead, 
###load in "matched_BORIS_behav_df_2021_02_all.csv" 
###and follow the script from that point (filtering)
###If you do not plan to change filtering, you also do not need to run these lines


#####load BORIS trial data, splice together with metadata######

#load BORIS summary csv (output)
fm1 <- read.csv("./boris_summary.csv", header=TRUE)

##need to remove .csv endings from VideoID for proper metadata matching##
library(dplyr)

#turn VideoID to characters
fm1$VideoID <- sapply(fm1$VideoID, as.character)

#substring character strings from start to end -4 (".csv")
fm1 <- fm1 %>%
  mutate(VideoID = substr(VideoID, 1, nchar(VideoID)-4))

#turn VideoID back into a string
fm1$VideoID <- sapply(fm1$VideoID, toString)

#load metadata for all trials
for_match <- read.csv("./behaviour_metadata.csv", header = TRUE)

#match BORIS output with trial/individual metadata (uses custom match function defined above)

#format: matches dataset 1 and 2 using a matching column; arguments 3 and 4
#identify these columns in the respective datasets. In this case...
#it's df1: fm1$VideoID; df2: column 1. The 5th argument tells the function
#to turn the selected columns of df2 to character (otherwise they would be factors)
tester <- matchstick(fm1, for_match,fm1$VideoID, 1, c(3,4,5,6,7,10,11))

#uncomment line below to print out matched data set#
#write.csv(tester, "matched_BORIS_behav_df_2021_02_all.csv", row.names = FALSE)

##################Data filtering##########################
#can skip to here if match boris summary already exists#
###remove juveniles, trials below 15 seconds in duration
###set time cut-off below; species removal based on low n
###is now done in the 'analyze' files

cutoff_1 <- 15 #removes trials less than 15s in duration

####remove juveniles####
unfiltered <- read.csv("matched_BORIS_behav_df_2021_02_all.csv", header = T)###read in file
filtered_1 <- unfiltered[!(unfiltered$Sex=="J"),] ###remove juveniles

#uncomment line below to output a file without juveniles
# write.csv(filtered_1,"matched_boris_nj.csv", row.names = FALSE) ##output to file

####remove trials shorter than 15s#####
filtered_2 <- filtered_1[filtered_1$Duration >= cutoff_1,] ##remove trials with durations below 15s

#uncomment line below to output a file with both filters applied
# write.csv(filtered_2,"matched_boris_nj_dur15filt.csv", row.names = FALSE) ##output to file


###summarize data for each individual, across all of its trials####

##set some columns as character##
#species,BocciaID,site,submersion type,submersion termination type sex, BORISID#
char_cols <- c(33:37, 40,41, 2)
for(t in 1:length(char_cols)){
  filtered_2[,char_cols[t]] <- as.character(filtered_2[,char_cols[t]])
}
c_hold <- colnames(filtered_2) ##get column names
lfactor <- unique(filtered_2$Boccia.ID)## get list of individuals

##from here, script splits into best performance vs average performance generation code

####################
##Best performance##
####################

###names for best performance trial variables, determined by number of reinhales####
concat_max <- var_endings(c(c_hold[c(5:20, 24:31)]), "bp")
all <- c(colnames(filtered_2[1,c(32:41)]),concat_max,"num trials")

###########combine metadata, best performance column names
indv <- as.data.frame(matrix(ncol = length(all), nrow=length(lfactor)))###set up a blank data frame for individual data
colnames(indv) <- all
######

####iterate through all individuals, selecting trial with best performance for each#
for(i in 1:length(lfactor)) {
  sub <- filtered_2[filtered_2$Boccia.ID ==lfactor[i],]##get all trials for one individual
  
  num_trials <- length(sub[,1])##record number of total trials for that individual
  
  max_rb.c <- max(sub$num_all_rb)##find max observed reinhales in one trial for that individual
  
  if(max_rb.c > 0){##if the lizard reinhaled at least once during all of its trials
    sub <- sub[sub$num_all_rb == max_rb.c,]##pare trials down to just the best reinhale trial by count

  }
  max_dur <- max(sub$Duration)##find max duration in remaining trials (all if no rebreathing was observed)
  sub <- sub[sub$Duration == max_dur,]##pare trials down to just the longest duration trial (if trials were tied)
  #####
  
  ##separate out metadata
  rest <- sub[1,c(32:41)]
  
  ###line below is a remnant from older code where best score per category was taken###
  ###current version selects best rb trial, uses it for all variables##
  ##however, this statement achieves the same thing, since sub contains
  ##only the 'best performance' trial

  max_vec <- sapply(c(sub[, c(5:20, 24:31)]), max)

  ##combine scores with metadata
  indv[i,] <- as.vector(c(rest, max_vec, num_trials))
}
write.csv(indv, "bp_individual.csv", row.names = FALSE)

####################
##Mean performance##
####################

char_cols <- c(33:37, 40,41, 2)
for(t in 1:length(char_cols)){
  filtered_2[,char_cols[t]] <- as.character(filtered_2[,char_cols[t]])
}
c_hold <- colnames(filtered_2) ##get column names
lfactor <- unique(filtered_2$Boccia.ID)## get list of individuals

###generate variable names for mean trial stats (summarized across all trials)####
concat_mean <- var_endings(c(c_hold[c(5:20, 24:31)]), "mean")
all <- c(colnames(filtered_2[1,c(33:41)]),concat_mean,"num trials")

###########combine metadata, mean/best performance column names
indv_m <- as.data.frame(matrix(ncol = length(all), nrow=length(lfactor)))###set up a blank data frame for individual data
colnames(indv_m) <- all
######

####iterate through all individuals, averaging across trials##
for(i in 1:length(lfactor)) {
  
  sub <- filtered_2[filtered_2$Boccia.ID ==lfactor[i],]##get all trials for one individual
 
  num_trials <- length(sub[,1])##record number of total trials for that individual
  
  max_rb.c <- max(sub$num_all_rb)##find max observed reinhales in one trial for that individual

  ##separate out metadata
  
  rest <- sub[1,c(33:41)]

  #generate designation for rebreathing 'type' based on best trial
  #all other fields will be averaged
  if(max_rb.c >4){
    rb_des <- "sustained"
  }else if(max_rb.c >0){
    rb_des <- "incidental"
  }else{
    rb_des <- "no"
  }
  rest$Rebreathing.type <- rb_des
  
  ###max mean and sum vectors (can take max performance of any/all variables,#### 
  ###sum of variables if desired; uncomment lines below)
  
  mean_vec <- sapply(c(sub[,c(5:20, 24:31)]), mean)##works in both trial averaged and best performance trials

  ##combine new scores
  indv_m[i,] <- as.vector(c(rest, mean_vec, num_trials))
}

#write averaged score file
write.csv(indv_m, "mean_individual.csv", row.names = FALSE)
