##Unified script for assessing rebreathing traces##

#set working directory to folder that script is in (adjust if necessary)#
if(Sys.getenv("RSTUDIO") == "1")
	setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#if not using RStudio, adjust to suit your system
#e.g., setwd("C:/your_directory")

library(lmtest)###load packages
library(ggplot2)
library(MuMIn)
library(lme4)

##read in filepaths for each oxygen trace .csv to be analyzed and plotted
lis <- read.csv("rebreathe_test_path_tl.csv", header = T)

##remove juveniles
lis <- lis[-c(43,53,54),]

###############################################################
####Generate curve stats, AICs, etc for all cleaned traces#####
###############################################################

##open pdf/png for printing (adjust classpath if needed)
pdf("all_traces.pdf", width = 20, height = 20)
par(mfrow=c(7,8))##adjust plotting presents; sets it up so R prints 7x8 graphs per page
par("mar" = c(3,3,3,3), mgp=c(2,1,0))##adjust plotting margins, etc.

##initialize dataframe to hold trial trace model stats
new_out <- data.frame(matrix(ncol =30 , nrow = length(lis[,1])))
colnames(new_out) <- c("plot", "species", "duration","reads","r2_lin","r2_cur", "f_read", "f_read_fit", "w_temp", "lin_aic", "quad_aic", "exp_aic", "asym_aic", "lin_lik", "quad_lik", "exp_lik", "asym_lik", "lin_int", "lin_slope", "quad_int", "quad_x", "quad_x2", "quad_dir", "exp_int", "exp_cur", "asym_int", "asym_cur", "asym_asym", "p_lin","p_exp")

#initialize lists, counters
countr <- 0
liklist <- list()
k_lis <- list()

for(i in 1:length(lis[,1])){##loops through plots
  #load string manipulation package
  require("stringi")
  tex <- ""
  path <- lis[i,1] ##Paste in path from file
  p2 <- paste(path, ".csv", sep="")##add extension (this var used to open file)
  nam <- substr(path, stri_locate_last_fixed(path, "/")+1, stri_length(path))#trial name, used for output df (removes path except for file name)
  
  start <- lis[i,2]##load in start of trial time from file
  end <- lis[i,3]##load in end of trial time from file
  duration <- end-start #calculate duration
  
  air <- lis[i,4]##load in average lab air value from file
  water <- lis[i, 6]##load in average water value from file

  oxy <- read.csv(p2, header = T)#read in oxygen data from csv
  colnames(oxy) <- c("time", "oxygen")#rename data columns
  
  freg <- oxy[(oxy$time > start & oxy$time <=end),]#downsize oxygen dataset to just region of dive
  num_rb <- length(freg$oxygen[!is.na(freg$oxygen)]) #get number of non-blank reads (bubble reads)
  freg$time= freg$time-freg$time[1]##standarize all dives to start at 0s (by subtracing first read from all read indices)
  
  ##compare/generate linear and quadratic fits
  ##linear fit
  fitF <- lm(oxygen ~ time, data=freg)
  
  ##quadratic fit
  fitR <- lm(oxygen ~ time + I(time^2), data=freg)
  
  ##AICs
  lin_a <- AICc(fitF)
  q_a <- AICc(fitR)
  
  sfitF <- summary(fitF)##can print these if interested
  sfitR <- summary(fitR)
  
  ##pass r2 values to holding variables
  r2F  <- sfitF$r.squared
  r2R <- sfitR$r.squared
  
  #pass p-values to holding variables
  pvalF <- summary(fitF)$coefficient[,"Pr(>|t|)"][2]
  pvalR <- summary(fitR)$coefficient[,"Pr(>|t|)"][2]
  
  #pass linear slope to holding variable
  l_slope <- sfitF$coefficients[,"Estimate"][2]
  
  #exponential decay model
  nls_t <- nls(oxygen~a*exp(-b*time),data=freg, start=list(a=180,b=0.003))
  nls_a <- AICc(nls_t)
  pvalexp <- summary(nls_t)$coefficients[2,4]
  
  ##asymptotic exponential decay
  nls_2 <- "NULL"#set variables to null first in case analysis fails
  nls_2aic <- "NULL"
  nls_2lik <- "NULL"
  nls_2int <- "NULL"
  nls_2cur <- "NULL"
  nls_2asym <- "NULL"
  #run in try statement to avoid crashing script (Some curve fits are singular and are not used); error messages silenced
  try({nls_2 <- nls(oxygen ~ (a - (a - b) * exp(- c * time)), data=freg, start=list(b=180, a=20, c=0.03 ))
  
  #pass asymptotic model fit scores and parameters to file
  nls_2aic <- AICc(nls_2)
  nls_2lik <- logLik(nls_2)
  nls_2int <- coef(nls_2)[1]
  nls_2cur <- coef(nls_2)[3]
  nls_2asym <- coef(nls_2)[2]
  
  }, silent=TRUE)
  
  #####last values, read length, final and extrapolated###
  rsnp <- freg[!is.na(freg[,2]),]
  readlen <- rsnp[length(rsnp[,2]),1] - freg[1,1]
  finv <- rsnp[length(rsnp[,2]),2]
  finv_fit <- predict(fitF, newdata = data.frame(time=c(duration)))
  
  ###homo/heteroscedasticity tests####
  require(lmtest)
  hetero_p <- bptest(fitF)$p.value
  
  ##print traces, including a species label
  ##tex and text lines below (commented out) add in additional metadata details to each plot
  ##metadata is set for each species; uncomment to print
   if(i<18){###separates out Anolis aquaticus trials
      if(countr==1){
        plot(freg[,2] ~ freg[,1], ylab = "Oxygen (hPa)", xlab="Time (s)", ylim=c(50, 165), xlim=c(0, 240))##prints plot with species label in x-axis and color
      }else{
        plot(freg[,2] ~ freg[,1], ylab = NA, xlab=expression(italic("aquaticus")), ylim=c(50, 165), xlim=c(0, 240))##prints plot with species label in x-axis and color
        
      }
      # tex <- paste("Temp:", lis$temp[i], sep=" ")##prints trial water temperature to graph
      # tex <- paste(tex,"\n Dur:", duration, "\n reads:", num_rb)
      # tex <- paste(tex,"\n readDur:", readlen, " \n FinVR:", round(finv,5), " \n Finfit:", round(finv_fit,5), " \n r2lin:", round(r2F,5)," \n p_lin:",round(pvalF,5),"\n slope:", round(l_slope,5), " \n r2c:", round(r2R,5), " \n lik", cur, " \n homhetp:", round(hetero_p, 5), "\n ", nam)
      # # text(x=200, y=140, labels=tex, cex=0.5)
      species_l <- "aquaticus"
    } else if(i>17 && i < 43){###separates out barkeri trials
      plot(freg[,2] ~ freg[,1], ylab = NA, xlab=expression(italic("barkeri")), ylim=c(60, 225), xlim=c(0, 580), col="blue")
      # tex <- paste("Temp:", lis$temp[i], sep=" ")
      # tex <- paste(tex,"\n Dur:", duration, "\n reads:", num_rb)
      # tex <- paste(tex,"\n readDur:", readlen, " \n FinVR:", round(finv,5), " \n Finfit:", round(finv_fit,5), " \n r2lin:", round(r2F,5)," \n p_lin:",round(pvalF,5),"\n slope:", round(l_slope,5), " \n r2c:", round(r2R,5), " \n lik", cur, " \n homhetp:", round(hetero_p, 5), "\n ", nam)
      # # text(x=423, y=180, labels=tex, cex=0.5)
      species_l <- "barkeri"
    } else if(i>42 && i<53){###separates out maculigula trials
      plot(freg[,2] ~ freg[,1], ylab = NA, xlab=expression(italic("maculigula")), ylim=c(95, 190), xlim=c(0, 475), col="purple")
      # tex <- paste("Temp:", lis$temp[i], sep=" ")
      # tex <- paste(tex,"\n Dur:", duration, "\n reads:", num_rb)
      # tex <- paste(tex,"\n readDur:", readlen, " \n FinVR:", round(finv,5), " \n Finfit:", round(finv_fit,5), " \n r2lin:", round(r2F,5)," \n p_lin:",round(pvalF,5),"\n slope:", round(l_slope,5), " \n r2c:", round(r2R,5), " \n lik", cur, " \n homhetp:", round(hetero_p, 5), "\n ", nam)
      # # text(x=390, y=123, labels=tex, cex=0.5)
      species_l <- "maculigula"
    } else if(i>52 && i<83){###separates out oxylophus (swierk) trials
      plot(freg[,2] ~ freg[,1], ylab = NA, xlab=expression(italic("oxylophus")), ylim=c(95, 210), xlim=c(0, 475), col="orange")
      # tex <- paste("Temp:", lis$temp[i], sep=" ")
      # tex <- paste(tex,"\n Dur:", duration, "\n reads:", num_rb)
      # tex <- paste(tex,"\n readDur:", readlen, " \n FinVR:", round(finv,5), " \n Finfit:", round(finv_fit,5), " \n r2lin:", round(r2F,5)," \n p_lin:",round(pvalF,5),"\n slope:", round(l_slope,5), " \n r2c:", round(r2R,5), " \n lik", cur, " \n homhetp:", round(hetero_p, 5),"\n ", nam)
      # # text(x=390, y=123, labels=tex, cex=0.5)
      species_l <- "oxylophus"
    } else if(i>=83){##mahler lynchi trials
      plot(freg[,2] ~ freg[,1], ylab = NA, xlab=expression(italic("lynchi")), ylim=c(95, 180), xlim=c(0, 475), col="red")
      # tex <- paste("Temp:", lis$temp[i], sep=" ")
      # tex <- paste(tex,"\n Dur:", duration, "\n reads:", num_rb)
      # tex <- paste(tex,"\n readDur:", readlen, " \n FinVR:", round(finv,5), " \n Finfit:", round(finv_fit,5), " \n r2lin:", round(r2F,5)," \n p_lin:",round(pvalF,5),"\n slope:", round(l_slope,5), " \n r2c:", round(r2R,5), " \n lik", cur, " \n homhetp:", round(hetero_p, 5), "\n ", nam)
      # # text(x=390, y=123, labels=tex, cex=0.5)
      species_l <- "lynchi"
    }
    #prints linear model slope (red)
    abline(fitF$coefficients[1], fitF$coefficients[2], col='red')##print linear fit
    
    #prints exponential decay curve fit (green)
    lines(c(freg[,1]), predict(nls_t,list(time=c(freg[,1]))),lwd=2, col = "green")
    
    #prints asymptotic curve fit line (blue)
    try(lines(c(freg[,1]), predict(nls_2,list(time=c(freg[,1]))),lwd=2, col = "blue"), silent=TRUE)
    
    #prints quadratic fit line (pink)
    predicted.intervals <- predict(fitR,newdata=data.frame(time=c(freg[,1])),interval='confidence',
                                   level=0.99)##generate quadratic fit numbers for plotting
    lines(freg[,1],predicted.intervals[,1],col='pink',lwd=3)#plot quadratic fit

  
  #determine quadratic curve behaviour (opens 'up' or 'down')
  quad_dir <- "up"
  if(coef(fitR)[3] < 0){
    quad_dir <- "down"
  }
  
  new_out[i,] <- c(as.character(lis[i,1]), species_l, duration,num_rb,r2F,r2R,finv, finv_fit, lis$temp[i], AICc(fitF), AICc(fitR), AICc(nls_t), nls_2aic, logLik(fitF)[1], logLik(fitR)[1], logLik(nls_t)[1], nls_2lik, coef(fitF)[1], coef(fitF)[2], coef(fitR)[1],coef(fitR)[2], coef(fitR)[3], quad_dir,coef(nls_t)[1], coef(nls_t)[2], nls_2int, nls_2cur,nls_2asym, pvalF,pvalexp)

}

dev.off()  ##turn off printing device 

#output trace model stats to .csv 
write.csv(new_out, file="4model_trace_stats_final.csv")

###################
##Trace filtering##
###################

#load data
data <- read.csv("4model_trace_stats_final.csv",header = TRUE)

meta <- read.csv("trace_indv_corres.csv", header=TRUE)

#load matching function
source("./functions.R")

##match data with metadata##
mat <- matchstick(data,meta,data$plot,1,c(2,4))

data <- mat

# filter out shortest trials (shorter than 100s) / also add other quality filters
#duration filter (>=100 s)
data <- subset(data, data$duration >= 100)
#reads filter ( >=10)
data <- subset(data, data$reads >= 10)
#quality filter (r2 (linear) > 0.75)
data <- subset(data, data$r2_lin >= 0.75)

#ensure data/species are character vectors
data$species <- as.character(data$species)
data$plot <- as.character(data$plot)

#get list of individuals
indv_l <- unique(data$individual)

#generate storage df, assign colnames
ndf <- data.frame(matrix(ncol = length(data[1,]), nrow = length(indv_l)))
colnames(ndf) <- colnames(data)

#get longest trial for each individual; retain only this trial
for(t in 1:length(indv_l)){
  get_bp <- subset(data, data$individual==indv_l[t])
  ndf[t,] <- get_bp[get_bp$duration==max(get_bp$duration),]
}
ndf$individual <- indv_l
ndf$d_aic <- ndf$lin_aic - ndf$exp_aic

lin_bf <- subset(ndf, ndf$d_aic < 2)
lin_bf$bf <- rep("lin", length(lin_bf[,1]))
exp_bf <- subset(ndf, ndf$d_aic >= 2)
exp_bf$bf <- rep("exp", length(exp_bf[,1]))

combo <- rbind(lin_bf, exp_bf)
##print out R2 averages for table
mean(combo$r2_lin[combo$species=="aquaticus"])
mean(combo$r2_lin[combo$species=="barkeri"])
mean(combo$r2_lin[combo$species=="lynchi"])
mean(combo$r2_lin[combo$species=="maculigula"])
mean(combo$r2_lin[combo$species=="oxylophus"])

##print out p-value averages for table
mean(combo$p_lin[combo$species=="aquaticus"])
mean(combo$p_lin[combo$species=="barkeri"])
mean(combo$p_lin[combo$species=="lynchi"])
mean(combo$p_lin[combo$species=="maculigula"])
mean(combo$p_lin[combo$species=="oxylophus"])

##check if model fit is correlated with trial duration/etc
mod_f <- glm(formula = as.factor(bf)~w_temp*log(mass)*duration, data=combo, family="binomial")
summary(mod_f)

#check for relationship with speecies
mod_f <- glm(formula = as.factor(bf)~species, data=combo, family="binomial")
summary(mod_f) 

######################
##linear fit analyses##
#######################
lin <- lin_bf
#duration
plot(lin_slope~duration, data=lin)#plot linear slope vs duration
mod_exp <- lm(lin_slope~duration, data=lin)#regression
abline(mod_exp)#add slope to plot
summary(mod_exp)#assess lm
mod_exp$coefficients[2]#get slope for table

#species
mod_exp <- lm(lin_slope~as.factor(species), data=lin)
summary(mod_exp)
boxplot(lin_slope~species, data=lin, ylim=c(-0.5,0.35))#print boxplot
tk <- aov(lin_slope~as.factor(species), data=lin)#input for tukeyhsd
TukeyHSD(tk)#run multiple comparisons

#water temp
plot(lin_slope~w_temp, data=lin)#plot linear slope vs water temp
mod_exp <- lm(lin_slope~w_temp, data=lin)#regression
abline(mod_exp)#add slope to plot
summary(mod_exp)#assess lm
mod_exp$coefficients[2]#get slope

#log(mass)
plot(lin_slope~log(mass), data=lin)#plot linear slope vs ln(mass)
mod_exp <- lm(lin_slope~log(mass), data=lin)#regression
abline(mod_exp)#add slope to plot
summary(mod_exp)#assess lm
options(scipen=10)#adjust scientific notation sensitivity
mod_exp$coefficients[2]#slope

###########
#run analyses on half-life
###########

##see notes on code from above analyses
#convert curve parameter to half-life, calculate half-life change/s
hlf <- exp_bf
hlf$t_2 <- log(2)/hlf$exp_cur

#duration
plot(t_2~duration, data=hlf)
mod_exp <- lm(t_2~duration, data=hlf)
abline(mod_exp)
summary(mod_exp)

#species
boxplot(t_2~species, data=hlf)
mod_exp <- lm(t_2~as.factor(species), data=hlf)
summary(mod_exp)
tk <- aov(t_2~as.factor(species), data=hlf)
TukeyHSD(tk)

#water temp
plot(t_2~w_temp, data=hlf)
mod_exp <- lm(t_2~w_temp, data=hlf)
abline(mod_exp)
summary(mod_exp)
mod_exp$coefficients[2]#slope

#log(mass)
plot(t_2~log(mass), data=hlf)
mod_exp <- lm(t_2~log(mass), data=hlf)
abline(mod_exp)
summary(mod_exp)
options(scipen=10)
mod_exp$coefficients[2]#slope

##################################################
#Plot all trials that meet filtering requirements#
##################################################

combo <- combo[order(combo$species),]

####plot only best trials 2020-07-18####
lis <- read.csv("rebreathe_test_path_tl.csv", header = T)

#image printing mode; uncomment, and comment pdf to use
#png("filtered_traces.png", width = 20, height = 20, units = "in", res=72) ###image version

pdf("filtered_traces.pdf", width = 20, height = 20)
par(mfrow=c(5,5), cex.axis=1.5,font=2)##adjust plotting presents; sets it up so R prints 7x8 graphs per page
par("mar" = c(6,6,6,6), mgp=c(2,1,0))##adjust plotting margins, etc.


for(i in 1:length(combo$plot)){
  require("stringi")
  path <- combo[i,2] ##Paste in path from df
  p2 <- paste(path, ".csv", sep="")##add extension
  
  i_lis <- lis[lis$path==path,]#get appropriate trial line from metadata df
  start <- as.integer(i_lis[2])##load in start of trial time from file
  end <- as.integer(i_lis[3])##load in end of trial time from file
  duration <- end-start#calculate duration
  
  air <- i_lis[4]##load in average lab air value from file
  water <- i_lis[6]##load in average water value from file
  
  oxy <- read.csv(p2, header = T)#read in oxygen data from csv
  colnames(oxy) <- c("time", "oxygen")#rename data columns
  
  freg <- oxy[(oxy$time > start & oxy$time <=end),]#downsize oxygen dataset to just region of dive
  num_rb <- length(freg$oxygen[!is.na(freg$oxygen)])
  freg$time= freg$time-freg$time[1]##standarize all dives to start at 0s
  fitF <- lm(oxygen ~ time, data=freg)
  nls_t <- nls(oxygen~a*exp(-b*time),data=freg, start=list(a=180,b=0.003))
  
  if(combo$species[i]=="aquaticus"){
    if(i==1){
      plot(freg[,2] ~ freg[,1], ylim=c(50, 165), xlim=c(0, 240), xlab=NA, ylab=NA)##prints plot with species label in x-axis and color
      mtext(side=2, line=3, expression("O"[2]*" Partial Pressure (hPa)"), cex=2, font=1)
      mtext(side=1, line=4, "Time (s)", cex=2, font=1)
    }else{
      plot(freg[,2] ~ freg[,1], ylab = NA, xlab=NA, ylim=c(50, 165), xlim=c(0, 240))##prints plot with species label in x-axis and color
      mtext(side=1, line=4, bquote(italic("A. aquaticus")), cex=2)
    }
  }else if(combo$species[i]=="barkeri"){#checks if species label should be barkeri, selects appropriate scale, point colour
    plot(freg[,2] ~ freg[,1], ylab = NA, xlab=NA, ylim=c(60, 225), xlim=c(0, 580), col="blue")
    mtext(side=1, line=4, bquote(italic("A. barkeri")), cex=2)
  }else if(combo$species[i]=="lynchi"){#checks if species label should be lynchi, selects appropriate scale, point colour
    plot(freg[,2] ~ freg[,1], ylab = NA, xlab=NA, ylim=c(95, 180), xlim=c(0, 475), col="red")
    mtext(side=1, line=4, bquote(italic("A. lynchi")), cex=2)
  }else if(combo$species[i]=="maculigula"){#checks if species label should be maculigula, selects approriate scale, point colour
    plot(freg[,2] ~ freg[,1], ylab = NA, xlab=NA, ylim=c(95, 190), xlim=c(0, 475), col="purple")
    mtext(side=1, line=4, bquote(italic("A. maculigula")), cex=2)
  }else if(combo$species[i]=="oxylophus"){#checks if species label should be oxylophus, selects approriate scale, point colour
    plot(freg[,2] ~ freg[,1], ylab = NA, xlab=NA, ylim=c(95, 210), xlim=c(0, 475), col="orange")
    mtext(side=1, line=4, bquote(italic("A. oxylophus")), cex=2)
  }
  if(combo$bf[i] =="lin"){#if best fit is linear, add red linear slope
    abline(fitF$coefficients[1], fitF$coefficients[2], col='red')##print linear fit
  }else{#if best fit is exponential, add green line
    lines(c(freg[,1]), predict(nls_t,list(time=c(freg[,1]))),lwd=2, col = "green")
  }
}
dev.off()#turn off printing device

#######################################
##Plot representative trials for paper#
##Figure 4#########################
#######################################

###print representative trace plots for each semi-aquatic species###

if(Sys.getenv("RSTUDIO") == "1")
	setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##load required packages##
library(stringr)
library(stringi)
library(ggplot2)
#load version of propagate package with only taylor series estimation
#CRAN version uses both taylor and MC methods; MC takes a long time
#and yields identical results
source("propagate_taylor_only.R")
lis <- read.csv("rebreathe_test_path_tl.csv", header = T)

##extract 'best' traces for each species####
nl <- lis[c(60,88,36,55,14),]

spec <- c("oxylophus", "lynchi", "barkeri", "maculigula","aquaticus", "sham")
pl_r <- list()##create list for plots

####generate plots, put in list obj###
for(i in 1:5){
  path <- nl[i,1] ##Paste in path from file
  p2 <- paste(path, ".csv", sep="")##add extension
  nam <- substr(path, stri_locate_last_fixed(path, "/")+1, stri_length(path))
  start <- nl[i,2]##load in start of trial time from file
  end <- nl[i,3]##load in end of trial time from file
  duration <- end-start
  air <- nl[i,4]##load in average lab air value from file
  
  
  oxy <- read.csv(p2, header = T)#read in oxygen data from csv
  colnames(oxy) <- c("time", "oxygen")#rename data columns
  freg <- oxy[(oxy$time > start & oxy$time <=end),]#downsize oxygen dataset to just region of dive
  num_rb <- length(freg$oxygen[!is.na(freg$oxygen)])
  freg$time= freg$time-freg$time[1]##standarize all dives to start at 0s
  
  fitF <- lm(oxygen ~ time, data=freg)
  sfitF <- summary(fitF)##can print these if interested; commented out for plot only output
  r2F  <- sfitF$r.squared
  
  model_nls <- nls(oxygen ~ a* exp(-b * time), data=freg, start=list(a=180, b=0.03))
  time <- seq(from=1, to=length(freg$time), by = 1)
  
  #currently, predict.nls does not do confidence/prediction intervals
  # pred_nls <- predict(model_nls, list(time=time), interval = "confidence")
  
  pred_nls <- predictNLS(model_nls,data.frame(time), interval = "confidence", nsim=10)
  pred_nls_df <- data.frame(pred_nls)
  pred_nls <- cbind(freg, pred_nls_df[, c("X2.5.","X97.5.","mu.1")])
  pred_nls$lwr<- as.numeric(pred_nls$`X2.5.`)
  pred_nls$upr<- as.numeric(pred_nls$`X97.5.`)
  l_slope <- sfitF$coefficients[,"Estimate"][2]
  
  #plot oxygen data along with exponential fit (best model), and exponential fit CIs
  p <-ggplot(data=freg, aes(x=time,y=oxygen)) + 
    geom_point(size=7, fill="#4f66f0", color="#233175", pch=21) + 
    #geom_abline(slope=l_slope, intercept=sfitF$coefficients[1], col="#4251ad", size=1.25)+
    geom_abline(slope = 0, intercept = air, col="red", size=2)+ 
    #geom_smooth(method=lm, color='#4251ad', size=2) +
    #geom_smooth(method="nls",formula=y~a*exp(-b*x), method.args = list(start=c(a=180,b=0.003)), color="#7FFF00") +
    geom_line(data=pred_nls, aes(y=mu.1), color="#00c8ff", size=2)+
    geom_line(data=pred_nls, aes(y=upr), color = "black", linetype = "dashed",size=1)  +
    geom_line(data=pred_nls,aes(y=lwr), color = "black", linetype = "dashed",size=1)  +
    labs(x="Time (s)", y=expression("O"[2]*" Partial Pressure (hPa)"))
  
  r <-p +theme_light() + ylim(80, 220) + xlim(0, 425)
  
  
  r <-r + scale_color_manual(values = c("blue", "red")) + theme(axis.title = element_text(size=40,face="bold"),
                                                                axis.line = element_line(colour = "black", size = 2, linetype = "solid"),
                                                                axis.text = element_text(face="bold", color="black", size=25),
                                                                legend.text = element_text(face="bold", color="black", size=20),
                                                                axis.ticks = element_line(size=2),
                                                                plot.margin = unit(c(1,1,1,1), "cm"))
  
  
  
  pl_r[[i]] <- r + ggtitle(spec[i])+theme(plot.title = element_text(face="italic", size=30, hjust=0.5))
}

####add sham trial to paper plot###
##load trial, metadata##
ske <- read.csv("./sham trials/sham se key.csv", header=TRUE)
sh <-read.csv("./sham trials/sham trial 1.csv", header=TRUE)

red <-sh[(sh$time >= ske[1,2] & sh$time <=ske[1,3]),]#downsize oxygen dataset to just region of dive
red$time= red$time-red$time[1]##standarize all dives to start at 0s
# summary(lm(oxygen~time,data=red))

sham_mod <- lm(oxygen~time,data=red)
tsq <- seq(from=1, to=length(red$time), by=1)
ts <- na.omit(red)$time

pred_sham <- data.frame(time=ts,predict(sham_mod, interval="confidence"))

airsh <- ske[1,4]##get air value

##plot representative sham###
p <-ggplot(data=red, aes(x=time,y=oxygen)) +
  geom_point(size=7, fill="#4f66f0", color="#233175", pch=21) +
  geom_abline(slope = 0, intercept = airsh, col="red", size=2)+
  geom_smooth(method=lm, color='#00c8ff', size=2)+
  geom_line(data=pred_sham, aes(y=fit), color="#00c8ff", size=2)+
  geom_line(data=pred_sham, aes(y=upr), color = "black", linetype = "dashed",size=1)  +
  geom_line(data=pred_sham, aes(y=lwr), color = "black", linetype = "dashed",size=1)+
  labs(x="Time (s)", y=expression("O"[2]*" Partial Pressure (hPa)"))

r <-p +theme_light() + ylim(80, 220) + xlim(0, 425)


r <-r + scale_color_manual(values = c("blue", "red")) + theme(axis.title = element_text(size=40,face="bold"),
                                                              axis.line = element_line(colour = "black", size = 2, linetype = "solid"),
                                                              axis.text = element_text(face="bold", color="black", size=25),
                                                              legend.text = element_text(face="bold", color="black", size=20),
                                                              axis.ticks = element_line(size=2),
                                                              plot.margin = unit(c(1,1,1,1), "cm"))

pl_r[[6]] <- r+ggtitle("sham")+theme(plot.title = element_text(size=30, hjust=0.5)) ###add sham plot to list

#add patchwork library for stitching plots together
library(patchwork)

##load seymour plot function
source("oxygen_consumption.R")

#get oxygen consumption plot with seymour data
pcon_t <- get_seymour_plot()

###edit axis labels###
pl_r[[1]] <- pl_r[[1]] + theme(axis.title.x=element_blank())
pl_r[[2]] <- pl_r[[2]] + theme(axis.text.y=element_blank(),
                               axis.title=element_blank(), axis.ticks.y = element_blank())
pl_r[[3]] <- pl_r[[3]] + theme(axis.text.y=element_blank(),
                               axis.title=element_blank(), axis.ticks.y = element_blank())
pl_r[[5]] <- pl_r[[5]]
pl_r[[4]] <- pl_r[[4]] + theme(axis.text.y=element_blank(),
                               axis.title.y=element_blank(), axis.ticks.y = element_blank(),axis.title.x=element_blank())
pl_r[[6]] <- pl_r[[6]] + theme(axis.text.y=element_blank(),
                               axis.title=element_blank(), axis.ticks.y = element_blank())
#adjust title fonts, margins
pcon_t <- pcon_t + theme_light()+theme(axis.title = element_text(size=40,face="bold"),
                                       axis.line = element_line(colour = "black", size = 2, linetype = "solid"),
                                       axis.text = element_text(face="bold", color="black", size=25),
                                       legend.text = element_text(face="bold", color="black", size=20),
                                       axis.ticks = element_line(size=2),
                                       plot.margin = unit(c(1,1,1,1), "cm"))


#plot spacing design
design <- "
AABB#CCDD
EEFF#GGGG
"

pdf("fig4.pdf", width = 40, height = 21)
pl_r[[1]]+pl_r[[2]]+pl_r[[3]]+pl_r[[4]]+pl_r[[5]]+pl_r[[6]]+pcon_t+plot_layout(design = design, guides = "collect")
dev.off()

#################################################
###sham bubble analysis and plotting, for SOM####
###Figure S3#####################################

fna <- c("sham trial 1.csv", "sham trial 2.csv", "sham trial 3.csv", 
         "sham trial 4.csv", "sham trial 5.csv")

ske <- read.csv("./sham trials/sham se key.csv", header=TRUE) # load metadata for sham trials

sham_l <- list()

for(j in 1:5){#run through
  airsh <- ske[j,4]
  pa <- paste("./sham trials/",fna[j],sep="")#adjust trial name to directory
  sham <- read.csv(pa, header=TRUE)#read in sham trial
  sre <-sham[(sham$time >= ske[j,2] & sham$time <=ske[j,3]),]#downsize oxygen dataset to just region of sham dunk
  sre$time= sre$time-sre$time[1]##standarize all sham dunks to start at 0s
  ls <- summary(lm(sre$oxygen ~ sre$time))
  
  #print out stats for S3 (shown within plots in SOM)
  print(paste("for trial ", fna[j], " Fstat= ",ls$fstatistic[1]," df1= ",ls$fstatistic[2], " df2= ", ls$fstatistic[3]," slope= ", ls$coefficients[2]," r2= ",ls$r.squared, sep="") )
  
  #plot sham trace
  p <-ggplot(data=sre, aes(x=time,y=oxygen)) + 
    geom_point(size=5, col="#4f66f0") + 
    geom_abline(slope = 0, intercept = airsh, col="red", size=2)+
    labs(x="Time (s)", y="Partial Pressure of Oxygen (hPa)")
  
  #adjust x, y axes
  r <-p + theme_classic() + ylim(80, 220) + xlim(0, 600)
  
  #adjust line colours, axis titles
  r <-r + scale_color_manual(values = c("blue", "red")) + theme(axis.title = element_text(size=25,face="bold"),
                                                                axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
                                                                axis.text = element_text(face="bold", color="black", size=15),
                                                                legend.text = element_text(face="bold", color="black", size=12),
                                                                axis.ticks = element_line(size=2),
                                                                plot.margin = unit(c(1,1,1,1), "cm"))
  sham_l[[j]] <- r

}

#design for patchwork plot stitching
design <- "
AABB##
CCDDEE
"

##print figure S3 to pdf (stats in SOM were added manually using Word text boxes)
pdf("s3.pdf", width = 30, height = 21)
sham_l[[1]]+sham_l[[2]]+sham_l[[3]]+sham_l[[4]]+sham_l[[5]]+plot_layout(design = design, guides = "collect")
dev.off()

