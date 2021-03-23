#####subsample 3 version#####
####subsample all species down to 3###

##this code will generate the components of figure S1 (best-performance subsampled results)
##the latter half of this script will plot trial-averaged subsampled results (not shown in paper)

#####################
#load packages, data#
#####################

###use line below if using RStudio to setwd to folder with script####
##otherwise, setwd to suit your placement of the files##
if(Sys.getenv("RSTUDIO") == "1")
	setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

####load packages####
library(geiger)
library(ape)
library(phytools)
library(ggplot2)
library(nlme)
library(dplyr)
library(car)
library(caper)

####define functions--run these lines into memory####
source("functions.R")
source("rb_assess.R")
##end of functions

##To recreate plots used in paper SOM, skip to readRDS statement on line 242
##otherwise, this code will generate a novel subsampled dataset which will not exactly resemble
##the data shown in the paper

aq_subsample <- function(indv, spec, ht, reps, type='bperformance'){
  # ##debug##
  # indv <- indv
  # spec <-list_keep
  # ht <-ht
  # reps <- 1000
  # type='bperformance'
  # ##
  if(type=='bperformance'){
    pt<-1
  }else if(type=='mean'){
    pt<-2
  }else{
    print("Error: incorrect analysis type (must be 'bperformance' or 'mean'")
    break;
  }
  
  
  #generate list
  df <- list()
  
  #get species vector (drop any species that have been filtered out from factor)
  if(is.factor(spec))
    spec <- as.character(droplevels(spec))
  
  #get number of rows needed for output dataframe (reps*species)
  row_n1 <-length(spec)*reps
  
  #create output df for subsample sun results
  pdfr <- data.frame(matrix(nrow=row_n1,ncol = 13))
  
  #create output df for print F values, significance scores
  prfsc <- data.frame(matrix(nrow=reps,ncol = 12))

  
  for(v in 1:reps) {
    #print row header so that each run can be distinguished
    cradd <- NA
    for (i in 1:length(spec)) {
      #concatenate run + run# + species
      cradd[i] <- paste("run",v,spec[i])
      
      #grab all data in table for one species
      sel <- data.frame(indv[indv$Species == spec[i], ])
      
      #
      sel <- sel[sample(1:nrow(sel), 3, replace = FALSE), ]
      if(i==1){
        new_in <- sel
      }else{
        new_in <- rbind(new_in, sel)
      }
    }
    
    if(pt==1){
        for_ag <- new_in[, c("Species",
                         "num_all_rb.bp",
                         "rb_full_rate.bp",
                         "rate_gul.bp")]
    }else if(pt==2){
      for_ag <- new_in[, c("Species",
                           "num_all_rb.mean",
                           "rb_full_rate.mean",
                           "rate_gul.mean")]
      
    }
    agg <- aggregate(. ~ Species, for_ag, mean)
    nn <- aggregate(. ~ Species, for_ag, length)
    agg$n <- nn[,2]
    rownames(agg) <- agg$Species
    
    agg <- agg[spec,]
    poetree <- read.tree("Poe2017timetree_sis_ed.tre")
    nc <- name.check(poetree, agg, data.names = agg$species)
    npt <- drop.tip(poetree, tip=nc$tree_not_data)
    
    ##generate table of rebreathing types
    rb_type <- table(new_in$Species, new_in$Rebreathing.type)
    
    ###convert table to dataframe
    rbt <- data.frame(species=rownames(rb_type),
                      incidental=as.vector(rb_type[,1]), 
                      no=as.vector(rb_type[,2]), 
                      sustained=as.vector(rb_type[,3]))
    
    #get sample size
    rbt$n <- rbt$incidental+rbt$no+rbt$sustained 
    
    ##drop species with fewer than 3 individuals
    rbt <- rbt[rbt$n >2,]
    
    ##match up habitats
    rbt <- matchstick(rbt, ht, rbt$species, 1,2)
    
    ###fix up column names, make habitat a factor
    rn <- colnames(rbt)
    rn[6] <- "habitat"
    colnames(rbt) <- rn
    rbt[,6] <- as.factor(rbt[,6])
    ####
    
    ##turn count data into proportion data###
    rbt$rb_any <- (rbt$sustained+rbt$incidental) / rbt$n
    rbt$incidental <- rbt$incidental/rbt$n
    rbt$no <- rbt$no/rbt$n
    rbt$sustained <- rbt$sustained/rbt$n
    
    ##logit transform props that will be analyzed
    ##adjusts 0 to 0.025 and 1 to 0.975 to avoid 0/1
    rbt$logit_sr <- logit(rbt$sustained)
    rbt$logit_rb_any <- logit(rbt$rb_any)
    rownames(rbt) <- rbt$species
    rbt <- rbt[spec,]
    
    comb <- cbind(agg,rbt[,-1])
    
    ###analyze logit proportions using caper####
    compphy <- comparative.data(npt, comb, names.col = 'Species')
    pglsModel <- pgls(logit_sr~habitat,
                      data = compphy, lambda=1)
    pglsModel2 <- pgls(logit_rb_any~habitat,
                       data = compphy, lambda=1)
    if(pt==1){
      pglsModel3 <- pgls(num_all_rb.bp~habitat,
                       data = compphy, lambda=1)
      pglsModel4 <- pgls(rb_full_rate.bp~habitat,
                       data = compphy, lambda=1)
      pglsModel5 <- pgls(rate_gul.bp ~habitat,
                       data = compphy, lambda=1)
    }else if (pt==2){
      pglsModel3 <- pgls(num_all_rb.mean~habitat,
                         data = compphy, lambda=1)
      pglsModel4 <- pgls(rb_full_rate.mean~habitat,
                         data = compphy, lambda=1)
      pglsModel5 <- pgls(rate_gul.mean ~habitat,
                         data = compphy, lambda=1)
    }
    ##output rebreathing type PGLS results
    prfsc[v,]<- c(anova(pglsModel)$`F value`[1],
                  anova(pglsModel)$`Pr(>F)`[1],
                  anova(pglsModel2)$`F value`[1],
                  anova(pglsModel2)$`Pr(>F)`[1],
                  anova(pglsModel3)$`F value`[1],
                  anova(pglsModel3)$`Pr(>F)`[1], 
                  anova(pglsModel4)$`F value`[1],
                  anova(pglsModel4)$`Pr(>F)`[1],
                  anova(pglsModel5)$`F value`[1],
                  anova(pglsModel5)$`Pr(>F)`[1],
                  anova(pglsModel)$Df[1],
                  anova(pglsModel)$Df[2]
    )
    rownames(comb) <- comb$Species
    tpm <- comb[spec,]
    rownames(tpm) <- cradd
    if(v==1){
      # iterate=v+4
      pdfr <- tpm
      # iterate <- v+5
    }else{
      # i2 <- iterate+4
      pdfr <- rbind(pdfr,tpm)
      # iterate <- iterate+5
    }
    
    
  }
  colnames(pdfr) <- colnames(comb)
  df[[1]] <- pdfr
  colnames(prfsc) <- c("f_lsr", "pr_lsr", "f_lar", "pr_lar", "f_rbc", "pr_rbc", "f_rbr", "pr_rbr", "f_gul", "pr_gul", "df1", "df2")
  df[[2]] <- prfsc
  return(df)
}


#################################
###aquatic subsample analysis####
#################################


#best performance#
##load data from file (if using a custom filtered file, change name here)##
indv <- read.csv("bp_individual.csv", header=TRUE)

###load habitat for each individual###
ht <- read.csv("habitats.csv", header = TRUE)##load habitat associations

##match up habitats w individuals
te <- matchstick(indv, ht, indv$Species, 1,2)##match up habitats

###fix up column names, make habitat a factor
ncn <- colnames(te)
ncn[36] <- "habitat"
colnames(te) <- ncn
te[,36] <- as.factor(te[,36])

###generate aggregated max reinhales observed for each species####
ag_n <- aggregate(indv$num_all_rb.bp, list(species=indv$Species), FUN=function(x) c(max=max(x), sd=sd(x),n=length(x), se=sd(x)/sqrt(length(x))))
tem <- data.frame(ag_n$x)
ag_n <- cbind(ag_n$species, tem)
colnames(ag_n) <- c("species", colnames(tem))
#####

##make list of species with at least 3 individuals
list_keep <- ag_n$species[ag_n$n >2]

###test out aquatic subsample###
df2 <-aq_subsample(indv, spec=list_keep, ht=ht, reps=10000)

#subsampled dataframe takes awhile to generate; uncomment line below to save
# saveRDS(df2, "aq_subsample_bp10000.rds")

#can load previously generated subsampled df from file; uncomment line below
df2 <- readRDS("aq_subsample_bp10000.rds")

df10 <- df2[[1]]
df3<- data.frame(df2[[2]])

##boxplots of subsampled species means for variables of interest#

#################
##gular pumping##
#################

#############
##boxplot for sustained rebreathing props
nv <- df10[,c("Species", "sustained", "incidental","no", "habitat")]
vio_sust2<-ggplot(data=nv, aes(x=Species,y=sustained,fill=habitat)) +
  geom_boxplot(width=0.5)+
  labs(x = "Species", y = as.character("Avg. Prop. sustained rebreathing"), fill = "Habitat")
vio_sust2 <- vio_sust2+ scale_fill_manual(labels = c("Semi-aquatic", "Non-aquatic"),values=c("#5DADE2", "#E59866"))+theme_bw()


vio_sust2 <- vio_sust2 + theme(text = element_text(size=20), 
                        axis.text.x = element_text(angle=90, hjust=1))+
                        ylim(0,1)
png("s1a.png", width = 10, height=8, units="in", res=72 )
vio_sust2
dev.off()

##threshold lsr
thresh <- qf(.95, df1=1, df2=18)
mean(df3$pr_lsr)
mean(df3$f_lsr)
mod <- df3$f_lsr[df3$f_lsr<thresh]
perc <- length(mod)/10000*100
perc
countr <- 0

##checks if aquatic mean is larger than non-aquatic for all 10000 runs
countr <- 0
for(i in 1:10000){
  tg1 <- (i-1)*20 + 1
  tg2 <- tg1 + 19
  selrun <- df10[tg1:tg2,]
  if(mean(selrun[selrun$habitat=="aquatic",]$sustained) > 
     mean(selrun[selrun$habitat=="non-aquatic",]$sustained)){
    countr <- countr + 1
  }
}
countr

png("s1b.png", width = 10, height=8, units="in", res=72 )
hist(df3$f_lsr, xlab="F-statistic", breaks=35, xlim=c(0,100), main=NULL)
abline(v=thresh, col="red", lty=2, lwd=3)
dev.off()

########
###subsampled incidental rebreathing
nv <- df10[,c("Species", "incidental","sustained", "habitat")]
nv$any <- nv$incidental+nv$sustained

vio_sust2<-ggplot(data=nv, aes(x=Species,y=any,fill=habitat)) +
  geom_boxplot(width=0.5)+
  labs(x = "Species", y = as.character("Avg. Prop. any rebreathing"), fill = "Habitat")
vio_sust2 <- vio_sust2+ scale_fill_manual(labels = c("Semi-aquatic", "Non-aquatic"),values=c("#5DADE2", "#E59866"))+theme_bw()


vio_sust2 <- vio_sust2 + theme(text = element_text(size=20), 
                               axis.text.x = element_text(angle=90, hjust=1))+
  ylim(0,1)
png("s1c.png", width = 10, height=8, units="in", res=72 )
vio_sust2
dev.off()


##threshold lar
thresh <- qf(.95, df1=1, df2=18)
mean(df3$pr_lar)
mean(df3$f_lar)
mod <- df3$f_lar[df3$f_lar<thresh]
perc <- length(mod)/10000*100
perc

countr <- 0
for(i in 1:10000){
  tg1 <- (i-1)*20 + 1
  tg2 <- tg1 + 19
  selrun <- df10[tg1:tg2,]
  if(sum(mean(selrun[selrun$habitat=="aquatic",]$incidental),mean(selrun[selrun$habitat=="aquatic",]$sustained)) > 
     sum(mean(selrun[selrun$habitat=="non-aquatic",]$sustained),mean(selrun[selrun$habitat=="non-aquatic",]$sustained))){
    countr <- countr + 1
  }
}
countr

png("s1d.png", width = 10, height=8, units="in", res=72 )
hist(df3$f_lar, xlab="F-statistic", breaks=35, xlim=c(0,60), main=NULL)
abline(v=thresh, col="red", lty=2, lwd=3)
dev.off()

#####
##boxplot, rb count
nv <- df10[,c("Species","num_all_rb.bp", "habitat")]
vio_sust2<-ggplot(data=nv, aes(x=Species,y=num_all_rb.bp,fill=habitat)) +
  geom_boxplot(width=0.5)+
  labs(x = "Species", y = as.character("Avg. reinspires/trial"), fill = "Habitat")
vio_sust2 <- vio_sust2+ scale_fill_manual(labels = c("Semi-aquatic", "Non-aquatic"),values=c("#5DADE2", "#E59866"))+theme_bw()


vio_sust2 <- vio_sust2 + theme(text = element_text(size=20), 
                               axis.text.x = element_text(angle=90, hjust=1))+
  ylim(0,20)
png("s1e.png", width = 10, height=8, units="in", res=72 )
vio_sust2
dev.off()

##threshold rb count
thresh <- qf(.95, df1=1, df2=18)
mean(df3$pr_rbc)
mean(df3$f_rbc)
mod <- df3$f_rbc[df3$f_rbc<thresh]
perc <- length(mod)/10000*100
perc

png("s1f.png", width = 10, height=8, units="in", res=72 )
hist(df3$f_rbc, xlab="F-statistic", breaks=35, xlim=c(0,200), main=NULL)
abline(v=thresh, col="red", lty=2, lwd=3)
dev.off()

##get number of trials in which aquatic means exceeds non-aquatic
countr <- 0
for(i in 1:10000){
  tg1 <- (i-1)*20 + 1
  tg2 <- tg1 + 19
  selrun <- df10[tg1:tg2,]
  if(mean(selrun[selrun$habitat=="aquatic",]$num_all_rb.bp) > 
     mean(selrun[selrun$habitat=="non-aquatic",]$num_all_rb.bp)){
    countr <- countr + 1
  }
}
countr

#####
##boxplot, rb rate
nv <- df10[,c("Species","rb_full_rate.bp", "habitat")]
vio_sust2<-ggplot(data=nv, aes(x=Species,y=rb_full_rate.bp,fill=habitat)) +
  geom_boxplot(width=0.5)+
  labs(x = "Species", y = as.character("Avg. rebreathing rate"), fill = "Habitat")
vio_sust2 <- vio_sust2+ scale_fill_manual(labels = c("Semi-aquatic", "Non-aquatic"),values=c("#5DADE2", "#E59866"))+theme_bw()


vio_sust2 <- vio_sust2 + theme(text = element_text(size=20), 
                               axis.text.x = element_text(angle=90, hjust=1))+
  ylim(0,0.20)
png("s1g.png", width = 10, height=8, units="in", res=72 )
vio_sust2
dev.off()

##threshold rate
thresh <- qf(.95, df1=1, df2=18)
mean(df3$pr_rbr)
mean(df3$f_rbr)
mod <- df3$f_rbr[df3$f_rbr<thresh]
perc <- length(mod)/10000*100
perc


png("s1h.png", width = 10, height=8, units="in", res=72 )
hist(df3$f_rbr, xlab="F-statistic", breaks=35, xlim=c(0,100), main=NULL)
abline(v=thresh, col="red", lty=2, lwd=3)
dev.off()

#get number of trials in which aquatic mean exceeded non-aquatic
countr <- 0
for(i in 1:10000){
  tg1 <- (i-1)*20 + 1
  tg2 <- tg1 + 19
  selrun <- df10[tg1:tg2,]
  if(mean(selrun[selrun$habitat=="aquatic",]$rb_full_rate.bp) > 
     mean(selrun[selrun$habitat=="non-aquatic",]$rb_full_rate.bp)){
    countr <- countr + 1
  }
}
countr

###########
##boxplot for gular pumping
nv <- df10[,c("Species", "rate_gul.bp", "habitat")]
#init ggplot
vio_sust<-ggplot(data=nv, aes(x=Species,y=rate_gul.bp,fill=habitat)) +
  geom_boxplot(width=0.5)+
  labs(x = "Species", y = as.character("Avg. gular pumping rate"), fill = "Habitat")
vio_sust <- vio_sust+ scale_fill_manual(labels = c("Semi-aquatic", "Non-aquatic"),values=c("#5DADE2", "#E59866"))+theme_bw()

#adjust text size, angle
vio_sust <- vio_sust + theme(text = element_text(size=20), 
                             axis.text.x = element_text(angle=90, hjust=1))+
  ylim(0,0.2)

png("s1i.png", width = 10, height=8, units="in", res=72 )
vio_sust
dev.off()

#generate F stat cutoff
thresh <- qf(.95, df1=1, df2=18)

#get average f-stat, p-value
mean(df3$f_gul)
mean(df3$pr_gul)

#get all subsamples with p < thresh
mod <- df3$f_gul[df3$f_gul<thresh]

#see what percentage are nonsignificant
perc <- length(mod)/10000*100
perc

##checks if aquatic mean is larger than non-aquatic for all 10000 runs
countr <- 0
for(i in 1:10000){
  #generate start/end indices to grab appropriate rows from df for one run
  tg1 <- (i-1)*20 + 1
  tg2 <- tg1 + 19
  selrun <- df10[tg1:tg2,]
  #if aquatic mean is greater, add to counter
  if(mean(selrun[selrun$habitat=="aquatic",]$rate_gul.bp) > 
     mean(selrun[selrun$habitat=="non-aquatic",]$rate_gul.bp)){
    countr <- countr + 1
  }
}
countr


#print F-stat distribution
png("s1j.png", width = 10, height=8, units="in", res=72 )
hist(df3$f_gul, xlab="F-statistic", breaks=20, xlim=c(0,5), main=NULL)
abline(v=thresh, col="red", lty=2, lwd=3)
dev.off()


#############################
#Mean/trial-averaged version#
#############################

## The results of these analyses are not included in the paper or SOM ##

indv <- read.csv("mean_individual.csv", header=TRUE)

##load habitat associations for each individual
ht <- read.csv("habitats.csv", header = TRUE)

#match up habitat for each indidual###
te <- matchstick(indv, ht, indv$Species, 1,2)

###fix up column names, make habitat a factor
ncn <- colnames(te)
ncn[35] <- "habitat"
colnames(te) <- ncn
te[,35] <- as.factor(te[,35])
###

###generate aggregated max reinhales observed for each species
ag_n <- aggregate(indv$num_all_rb.mean, list(species=indv$Species), FUN=function(x) c(max=max(x), sd=sd(x),n=length(x), se=sd(x)/sqrt(length(x))))

##dataframe formatting
tem <- data.frame(ag_n$x)
ag_n <- cbind(ag_n$species, tem)
colnames(ag_n) <- c("species", colnames(tem))

##make list of species with at least 3 individuals
list_keep <- ag_n$species[ag_n$n >2]

##subsample mean individual data set
###test out aquatic subsample###
df2 <-aq_subsample(indv, spec=list_keep, ht=ht, reps=10000, type="mean")

#subsampled dataframe takes awhile to generate; uncomment line below to save
#saveRDS(df2, "aq_subsample_mean10000.rds")

#can load previously generated subsampled df from file; uncomment line below
#df2 <- readRDS("aq_subsample_mean10000.rds")

#split list from aq_subsample into dfs
df10 <- df2[[1]]
df3<- data.frame(df2[[2]])

##boxplots of subsampled species means for variables of interest#
##sustained/any rebreathing plots will be roughly the same (with some simulation differences
##So these are not included for the mean dataset. They can be recomputed
##using the best performance version

#############

#####
##boxplot, rb count (trial-averaged version)

nv <- df10[,c("Species","num_all_rb.mean", "habitat")]
vio_sust2<-ggplot(data=nv, aes(x=Species,y=num_all_rb.mean,fill=habitat)) +
  geom_boxplot(width=0.5)+
  labs(x = "Species", y = as.character("Avg. reinspires/trial"), fill = "Habitat")
vio_sust2 <- vio_sust2+ scale_fill_manual(labels = c("Semi-aquatic", "Non-aquatic"),values=c("#5DADE2", "#E59866"))+theme_bw()


vio_sust2 <- vio_sust2 + theme(text = element_text(size=20), 
                               axis.text.x = element_text(angle=90, hjust=1))+
  ylim(0,20)

png("s1e_mean.png", width = 10, height=8, units="in", res=72 )
vio_sust2
dev.off()

##threshold rb count
thresh <- qf(.95, df1=1, df2=18)
mean(df3$pr_rbc)
mean(df3$f_rbc)
mod <- df3$f_rbc[df3$f_rbc<thresh]
perc <- length(mod)/10000*100
perc

png("s1f_mean.png", width = 10, height=8, units="in", res=72 )
hist(df3$f_rbc, xlab="F-statistic", breaks=35, xlim=c(0,200), main=NULL)
abline(v=thresh, col="red", lty=2, lwd=3)
dev.off()

##get number of trials in which aquatic means exceeds non-aquatic
countr <- 0
for(i in 1:10000){
  tg1 <- (i-1)*20 + 1
  tg2 <- tg1 + 19
  selrun <- df10[tg1:tg2,]
  if(mean(selrun[selrun$habitat=="aquatic",]$num_all_rb.mean) > 
     mean(selrun[selrun$habitat=="non-aquatic",]$num_all_rb.mean)){
    countr <- countr + 1
  }
}


#####
##boxplot, rb rate
nv <- df10[,c("Species","rb_full_rate.mean", "habitat")]
vio_sust2<-ggplot(data=nv, aes(x=Species,y=rb_full_rate.mean,fill=habitat)) +
  geom_boxplot(width=0.5)+
  labs(x = "Species", y = as.character("Avg. rebreathing rate"), fill = "Habitat")
vio_sust2 <- vio_sust2+ scale_fill_manual(labels = c("Semi-aquatic", "Non-aquatic"),values=c("#5DADE2", "#E59866"))+theme_bw()


vio_sust2 <- vio_sust2 + theme(text = element_text(size=20), 
                               axis.text.x = element_text(angle=90, hjust=1))+
  ylim(0,0.20)
png("s1g_mean.png", width = 10, height=8, units="in", res=72 )
vio_sust2
dev.off()

##threshold rate
thresh <- qf(.95, df1=1, df2=18)
mean(df3$pr_rbr)
mean(df3$f_rbr)
mod <- df3$f_rbr[df3$f_rbr<thresh]
perc <- length(mod)/10000*100
perc


png("s1h_mean.png", width = 10, height=8, units="in", res=72 )
hist(df3$f_rbr, xlab="F-statistic", breaks=35, xlim=c(0,100), main=NULL)
abline(v=thresh, col="red", lty=2, lwd=3)
dev.off()

#get number of trials in which aquatic mean exceeded non-aquatic
countr <- 0
for(i in 1:10000){
  tg1 <- (i-1)*20 + 1
  tg2 <- tg1 + 19
  selrun <- df10[tg1:tg2,]
  if(mean(selrun[selrun$habitat=="aquatic",]$rb_full_rate.mean) > 
     mean(selrun[selrun$habitat=="non-aquatic",]$rb_full_rate.mean)){
    countr <- countr + 1
  }
}

###########
##boxplot for gular pumping
nv <- df10[,c("Species", "rate_gul.mean", "habitat")]
#init ggplot
vio_sust<-ggplot(data=nv, aes(x=Species,y=rate_gul.mean,fill=habitat)) +
  geom_boxplot(width=0.5)+
  labs(x = "Species", y = as.character("Avg. gular pumping rate"), fill = "Habitat")
vio_sust <- vio_sust+ scale_fill_manual(labels = c("Semi-aquatic", "Non-aquatic"),values=c("#5DADE2", "#E59866"))+theme_bw()

#adjust text size, angle
vio_sust <- vio_sust + theme(text = element_text(size=20), 
                             axis.text.x = element_text(angle=90, hjust=1))+
  ylim(0,0.2)

png("s1i_mean.png", width = 10, height=8, units="in", res=72 )
vio_sust
dev.off()

#generate F stat cutoff
thresh <- qf(.95, df1=1, df2=18)

#get average f-stat, p-value
mean(df3$f_gul)
mean(df3$pr_gul)

#get all subsamples with p < thresh
mod <- df3$f_gul[df3$f_gul<thresh]

#see what percentage are nonsignificant
perc <- length(mod)/10000*100
perc

##checks if aquatic mean is larger than non-aquatic for all 10000 runs
countr <- 0
for(i in 1:10000){
  #generate start/end indices to grab appropriate rows from df for one run
  tg1 <- (i-1)*20 + 1
  tg2 <- tg1 + 19
  selrun <- df10[tg1:tg2,]
  #if aquatic mean is greater, add to counter
  if(mean(selrun[selrun$habitat=="aquatic",]$rate_gul.mean) > 
     mean(selrun[selrun$habitat=="non-aquatic",]$rate_gul.mean)){
    countr <- countr + 1
  }
}

png("s1j_mean.png", width = 10, height=8, units="in", res=72 )
hist(df3$f_gul, xlab="F-statistic", breaks=20, xlim=c(0,5), main=NULL)
abline(v=thresh, col="red", lty=2, lwd=3)
dev.off()