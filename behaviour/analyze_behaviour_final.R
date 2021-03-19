####This script analyzes and plots the behavioural data#####
####based on best individual performances              #####

#####################
#load packages, data#
#####################

###use line below if using RStudio to setwd to folder with script####
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
library(caper)

####define functions####
source("functions.R")
source("./rb_assess.R")
#####End of functions

##this script analyzes a 'best performance' filtered data set

##load data from file (if using a custom filtered file, change name here)##
indv <- read.csv("bp_individual.csv", header=TRUE)

###load habitat for each individual###
ht <- read.csv("habitats.csv", header = TRUE)##load habitat associations

###generate aggregated max reinhales observed for each species####
ag_n <- aggregate(indv$num_all_rb.bp, list(species=indv$Species), FUN=function(x) c(max=max(x), sd=sd(x),n=length(x), se=sd(x)/sqrt(length(x))))
tem <- data.frame(ag_n$x)
ag_n <- cbind(ag_n$species, tem)
colnames(ag_n) <- c("species", colnames(tem))
#####

##make list of species with at least 3 individuals
list_keep <- ag_n$species[ag_n$n >2]
###this species list will be carried through all subsequent analyses###

##associate rows w species name##
rownames(ag_n) <- ag_n$species

##cut down data to just species that make n>2 filter
ag_n <- ag_n[list_keep,]

###add habitat
###create habitat vector for species###
hab_vec <- c("Non-aquatic", rep("Aquatic", 2), rep("Non-aquatic",5), "Aquatic", "Non-aquatic", "Aquatic", 
             rep("Non-aquatic", 2), "Aquatic", rep("Non-aquatic", 6))
ag_n$habitat <- as.factor(hab_vec)

###initial setup for combined data frame for plotting###
##add species, habitat vectors
sp_df <- cbind(as.character(ag_n$species), as.character(ag_n$habitat),as.numeric(ag_n$n))
colnames(sp_df) <- c("species", "habitat", "n")
tcol <- colnames(ag_n)
tcol <- tcol[c(2)]

##add reinspire max to primary df
to_add <- ag_n[,c(2)]
tcol<- var_endings(tcol, "rb_max")
sp_df <- data.frame(sp_df, rb_max=to_add)


###load phylo tree to conduct phylogenetic ANOVA
poetree <- read.tree("Poe2017timetree_sis_ed.tre")

###trim tree to focal taxa####
nc <- name.check(poetree, ag_n, data.names = ag_n$species)
npt <- drop.tip(poetree, tip=nc$tree_not_data)

##recheck tree, names
name.check(npt, ag_n)

###Brownian motion pgls, anova for max reinhales
compphy <- comparative.data(npt, ag_n, names.col = 'species')
pglsModel <- pgls(formula=max~habitat, 
                  data = compphy, lambda = 1)

# #####################
# ###nlme + ape mode###
# #####################
# 
# pglsModel<-gls(max~habitat, data=ag_n, correlation=
#                  corPagel(value=0.8, phy=npt, fixed=FALSE, 
#                          form = ~species), method="ML")
# #####################

##output caper PGLS results
anova(pglsModel)
#uncomment line below for additional info
# summary(pglsModel)

##Due to the tree used, lambda estimates were unstable when using ML estimation
##We set lambda at 1
##output caper lambda (and lambda ML plot)
#lm.lk<-pgls.profile(pglsModel, which="lambda")
#plot(lm.lk)
###############

###process/analyze additional variables, add to dataframe for plotting
###uses 'rb_assess' function
##uses indv data frame, selects column (2nd argument)
###to see scores, plots, run one by one
sp_df <- rb_assess(indv,12,"rb_count_avg", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,13,"rb_rate_avg", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,24,"gular_rate_avg", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,17,"rb_period_avg", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,22,"br_rate_avg", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,30,"ex_dur_avg", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,31,"soh_avg", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,32,"toh_avg", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,33,"hbn_avg", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,28,"exhale_pct", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,11,"duration", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,8,"svl", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,7,"mass", sp_df, list_keep, hab_vec)

#######################################################
####dive duration v habitat (best performance version)#
#best performance version of supmat plot S2g###########
#not included in paper#################################
#######################################################

##ggplot of duration v habitat
plex<-ggplot(data=sp_df, aes(x=habitat,y=`mean duration`,fill=habitat)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(fill="white", color="black", width=0.3)+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=3, fill="white", binwidth = 5)+
  labs(x = "Habitat", y = as.character("Average trial duration (s)"), fill = "Habitat")

##create legend, alter labels
plex <- plex + scale_fill_manual(labels = c("Semi-aquatic", "Non-aquatic"),values=c("#5DADE2", "#E59866"))+theme_bw()+
  scale_x_discrete(labels=c("Semi-aquatic", "Non-aquatic")) +
  ylim(c(0, 400))

##uncomment png, dev.off to output as an image
# png("duration_plot_2020-07-16.png", width = 8, height = 15, units = "in", res=72)
#pdf("S2g_bp.pdf",width = 8, height = 15)
plex
# dev.off()

####End figure#####


###create rebreathing type table for analysis, plotting###

##generate table of rebreathing types
rb_type <- table(indv$Species, indv$Rebreathing.type)

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

###analyze logit proportions using caper####
compphy <- comparative.data(npt, rbt, names.col = 'species')
pglsModel <- pgls(logit_sr~habitat,
                  data = compphy, lambda = 1)
#pglsModel<-gls(logit_rb_any~habitat, data=rbt, correlation=
                  #   corPagel(value=1, phy=npt, fixed=FALSE,
                    #          form = ~species), method="ML")

##output rebreathing type PGLS results
anova(pglsModel)
#summary(pglsModel)

##remove??##
# write.csv(rbt, "updated_rb_type_prop_table.csv")

####assess species average duration, other factors in a model
###that includes svl and mass (using PGLS)####
dur_syn <- data.frame(sp_df$`mean rb_rate_avg`,sp_df$`mean rb_count_avg`,sp_df$`rb_max`,sp_df$`mean gular_rate_avg`,sp_df$species,sp_df$habitat, sp_df$`mean duration`, sp_df$`mean svl`,sp_df$`mean mass`)
colnames(dur_syn) <- c("rbr","rbc","rbmax","gulr","species", "habitat", "duration", "svl", "mass")
rownames(dur_syn) <- dur_syn$species
compphy <- comparative.data(npt, dur_syn, names.col = 'species')

#rebreathing rate
#pglsModel <- pgls(formula=rbr~habitat*svl*log(mass), data = compphy, lambda=1) #covariate mode, results equivalent
pglsModel <- pgls(formula=rbr~habitat, data = compphy, lambda=1)
anova(pglsModel)
#summary(pglsModel)

#rebreathing count
#pglsModel <- pgls(formula=rbc~habitat*svl*log(mass), data = compphy, lambda=1)#covariate mode, results equivalent
pglsModel <- pgls(formula=rbc~habitat, data = compphy, lambda=1)

anova(pglsModel)

#rebreathing max
#pglsModel <- pgls(formula=rbmax~habitat*svl*log(mass), data = compphy, lambda=1)#covariates included, results equivalent
pglsModel <- pgls(formula=rbmax~habitat, data = compphy, lambda=1)
anova(pglsModel)
#summary(pglsModel)

#gular pumping
#pglsModel <- pgls(formula=gulr~habitat*svl*log(mass), data = compphy, lambda=1)#covariates included, results equivalent
pglsModel <- pgls(formula=gulr~habitat, data = compphy, lambda=1)
anova(pglsModel)
#summary(pglsModel)

#duration
pglsModel <- pgls(formula=duration~habitat*log(mass), data = compphy, lambda=1)
anova(pglsModel)
#summary(pglsModel)

##check svl/mass impact on rebreathing prop models

#add rebreathing type proportions to df
dur_syn$rb_any <- rbt$logit_rb_any
dur_syn$rb_sr <- rbt$logit_sr
compphy <- comparative.data(npt, dur_syn, names.col = 'species')

#test sustained rebreathing
pglsModel <- pgls(formula=rb_sr~habitat*svl*mass, data = compphy, lambda=1)
anova(pglsModel)

#test sustained rebreathing
pglsModel <- pgls(formula=rb_any~habitat*svl*mass, data = compphy, lambda=1)
anova(pglsModel)

#######Figure 3###########
####generate pie chart plot####

##collect bubble location columns
locs <- sp_df[,c(1,2,3,23,26,29)]

##print pie charts of locations for each species; uncomment to print to file##
#png("fig3.png", width = 28, height = 21, units = "in", res=72)
pdf("fig3.pdf", width = 28, height = 21)

##adjust plotting presents; sets it up so R prints 4x5 graphs per page
par(mfrow=c(4,5))
par(mai=c(0.15,0.15,0.15,0.15))
locs <- locs[order(locs$habitat,locs$species),]

##iterate through species
for(x in 1:20){
  if(sum(locs[x,4], locs[x,5],locs[x,6])>0){
    #print pie chart
    pie(as.numeric(locs[x,c(4,5,6)]), labels = NA, col = c("red", "purple", "blue"))
    
    ##print species name (italicized) and n
    ##bquote takes anything identified with the '.' operator and interprets it as a variable
    ##~ is a spacing operator
    mtext(bquote(italic("A.")~italic(.(as.character(locs$species[x]))) ~ (  .(as.character(locs$n[x])) ) ), side=1, cex=3)
  }
}

#uncomment when printing to file#
dev.off()
###for pie chart labels, sub in: labels = c("narial high", "top of head", "side of head")

##########
#Figure 2#
##########

#####load packages####
require(ggplot2)
require(ggtree)

###create phylogenetic tree plot for panel, tip order will be x labels for other plots#####
treeg <- ggtree(npt) +geom_tiplab(align = TRUE)
ord_c <- data.frame(treeg$data[7], treeg$data[4])
ord_c <- ord_c[1:20,]
ord_c <- ord_c[order(ord_c$y, decreasing=TRUE),]
order_p <- ord_c$label

##make species names row identifiers
rownames(sp_df) <- sp_df$species

##order df so that names will match up w phylo
sp_df <- sp_df[order_p,]

##necessary so that species will be in tree order
sp_df$species <- factor(sp_df$species, levels=order_p)

####max observed reinhale species plot (no error bars, single point)####
##left side, third plot from bottom##
pl4<-ggplot(data=sp_df, aes(x=species,y=rb_max, color=habitat)) +
  geom_point(stat="identity", size=4)+
  labs(x = "Species", y = as.character("Max reinspires observed"), fill = "Habitat")

pl4 <- pl4 + scale_color_manual(labels = c("Semi-aquatic", "Non-aquatic"),values=c("#5DADE2", "#E59866"))+theme_bw()


pl4 <- pl4 + theme(text = element_text(size=20), 
                   axis.text.x = element_text(angle=90, hjust=1))+
  ylim(0,72)##set ylims to equal

##remove axis labels for plots high up in the grid
pl4 <- pl4 +theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  plot.margin = unit(c(0.25,0,0.25,0), "cm"))
##uncomment to remove legend###
pl4 <- pl4 +theme(legend.position = "none")

###uncomment to remove y-axis labels####
# pl4 <- pl4 + theme(axis.text.y=element_blank(),
# axis.title.y=element_blank(),)
#####

####max observed reinhale habitat plot (violin)
##right side, third from bottom##
pl4h<-ggplot(data=sp_df, aes(x=habitat,y=rb_max,fill=habitat)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1)+
  labs(x = "Species", y = as.character("Max reinspires observed"), fill = "Habitat")


pl4h <- pl4h + scale_fill_manual(labels = c("Semi-aquatic", "Non-aquatic"),values=c("#5DADE2", "#E59866"))+theme_bw()


pl4h <- pl4h + theme(text = element_text(size=20), 
                     axis.text.x = element_text(angle=90, hjust=1))+
  ylim(0,72)
pl4h <- pl4h +theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.title.y=element_blank(),
                    plot.margin = unit(c(0.25,0,0.25,0), "cm")) #top, right, bottom, left
#uncomment to remove legend
# pl4h <- pl4h +theme(legend.position = "none")

####avg gular species plot####
###left side, bottom plot###
pl3<-ggplot(data=sp_df, aes(x=species,y=`mean gular_rate_avg`, color=habitat)) +
  geom_point(stat="identity", size=4)+
  geom_errorbar(aes(ymin=`mean gular_rate_avg`-(2*`se gular_rate_avg`), ymax=`mean gular_rate_avg`+(2*`se gular_rate_avg`)), width=.5, size=1)+
  labs(x = "Species", y = as.character("Avg. gular pumps/s"), fill = "Habitat")

pl3 <- pl3 + scale_color_manual(labels = c("Semi-aquatic", "Non-aquatic"),values=c("#5DADE2", "#E59866"))+theme_bw()


pl3 <- pl3 + theme(text = element_text(size=20), 
                   axis.text.x = element_text(angle=90, hjust=1),
                   plot.margin = unit(c(0.25,0,0.25,0), "cm"))+
  ylim(-0.05,0.35)##set ylims to equal

##uncomment to remove axis labels
# pl3 <- pl3 +theme(
#                   # axis.title.x=element_blank(),
#                   # axis.text.x=element_blank(),
#                   axis.text.y=element_blank(),
#                   axis.title.y=element_blank())

##comment to remove legend
pl3 <- pl3 +theme(legend.position = "none")
pl3 <- pl3 + coord_cartesian(ylim = c(0,0.35))
#####

####avg gular rate habitat plot (violin)
###right side, bottom plot###
pl3h<-ggplot(data=sp_df, aes(x=habitat,y=`mean gular_rate_avg`,fill=habitat)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1)+
  labs(x = "Species", y = as.character("Avg. gular pumps/s"), fill = "Habitat")


pl3h <- pl3h + scale_fill_manual(labels = c("Semi-aquatic", "Non-aquatic"),values=c("#5DADE2", "#E59866"))+theme_bw()


pl3h <- pl3h + theme(text = element_text(size=20), 
                     axis.text.x = element_text(angle=90, hjust=1))+
  ylim(-0.05,0.35)
pl3h <- pl3h +theme(
  ###uncomment to remove x-axis labels, legend
  # axis.title.x=element_blank(),
  # axis.text.x=element_blank(),
  # legend.position = "none",
  axis.text.y=element_blank(),
  axis.title.y=element_blank(),
  plot.margin = unit(c(0.25,0,0.25,0), "cm"))
pl3h <- pl3h + coord_cartesian(ylim = c(0,0.35))

####avg rebreathing rate species plot####
##left side, second from bottom
pl2<-ggplot(data=sp_df, aes(x=species,y=`mean rb_rate_avg`, color=habitat)) +
  geom_point(stat="identity", size=4)+
  geom_errorbar(aes(ymin=`mean rb_rate_avg`-(2*`se rb_rate_avg`), ymax=`mean rb_rate_avg`+(2*`se rb_rate_avg`)), width=.5, size=1)+
  labs(x = "Species", y = as.character("Avg. reinspires/s"), fill = "Habitat")

pl2 <- pl2 + scale_color_manual(labels = c("Semi-aquatic", "Non-aquatic"),values=c("#5DADE2", "#E59866"))+theme_bw()


pl2 <- pl2 + theme(text = element_text(size=20), 
                   axis.text.x = element_text(angle=90, hjust=1))+
  ylim(-0.05,0.3)##set ylims to equal

##remove axis labels for plots high up in the grid
pl2 <- pl2 +theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  plot.margin = unit(c(0.25,0,0.25,0), "cm"))

##uncomment to remove legend#
pl2 <- pl2 +theme(legend.position = "none")

##run to remove y-axis labels##
# pl2 <- pl2 +theme(axis.text.y=element_blank(),
#                   axis.title.y=element_blank())
##

pl2 <- pl2 + coord_cartesian(ylim = c(0,0.3))
#####

####avg rebreathing rate habitat plot (violin)
###right side, second from bottom
pl2h<-ggplot(data=sp_df, aes(x=habitat,y=`mean rb_rate_avg`,fill=habitat)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1)+
  labs(x = "Species", y = as.character("Avg. reinspires/s"), fill = "Habitat")


pl2h <- pl2h + scale_fill_manual(labels = c("Semi-aquatic", "Non-aquatic"),values=c("#5DADE2", "#E59866"))+theme_bw()


pl2h <- pl2h + theme(text = element_text(size=20), 
                     axis.text.x = element_text(angle=90, hjust=1))+
  ylim(-0.05,0.3)
pl2h <- pl2h +theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.title.y=element_blank(),
                    plot.margin = unit(c(0.25,0,0.25,0), "cm"))
##uncomment to remove legend
# pl2h <- pl2h +theme(legend.position = "none")

pl2h <- pl2h + coord_cartesian(ylim = c(0,0.25))


####avg rebreathing count species plot####
##left side, 4th from bottom
pl1<-ggplot(data=sp_df, aes(x=species,y=`mean rb_count_avg`, color=habitat)) +
  geom_point(stat="identity", size=4)+
  geom_errorbar(aes(ymin=`mean rb_count_avg`-(2*`se rb_count_avg`), ymax=`mean rb_count_avg`+(2*`se rb_count_avg`)), width=.5, size=1)+
  labs(x = "Species", y = as.character("Avg. reinspires/trial"), fill = "Habitat")

pl1 <- pl1 + scale_color_manual(labels = c("Semi-aquatic", "Non-aquatic"),values=c("#5DADE2", "#E59866"))+theme_bw()


pl1 <- pl1 + theme(text = element_text(size=20), 
                   axis.text.x = element_text(angle=90, hjust=1))+
  ylim(-5,35)##set ylims to equal

##remove axis labels for plots high up in the grid
pl1 <- pl1 +theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  plot.margin = unit(c(0.25,0,0.25,0), "cm"))

##uncomment to remove legend##
pl1 <- pl1 +theme(legend.position = "none")

###uncomment to remove y-axis labels
# pl1 <- pl1 +theme(axis.text.y=element_blank(),
#                   axis.title.y=element_blank())

pl1 <- pl1 + coord_cartesian(ylim = c(0,35))
#####

####avg rebreathing count habitat plot (violin)
##right side, 4th from bottom
pl1h<-ggplot(data=sp_df, aes(x=habitat,y=`mean rb_count_avg`,fill=habitat)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1)+
  labs(x = "Species", y = as.character("Avg. reinspires/trial"), fill = "Habitat")


pl1h <- pl1h + scale_fill_manual(labels = c("Semi-aquatic", "Non-aquatic"),values=c("#5DADE2", "#E59866"))+theme_bw()


pl1h <- pl1h + theme(text = element_text(size=20), 
                     axis.text.x = element_text(angle=90, hjust=1))+
  ylim(-5,35)
pl1h <- pl1h +theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.text.y=element_blank(),
                    axis.title.y=element_blank(),
                    plot.margin = unit(c(0.25,0,0.25,0), "cm"))
##uncomment to remove legend##
# pl1h <- pl1h +theme(legend.position = "none")

pl1h <- pl1h + coord_cartesian(ylim = c(0,35))

########habitat/rebreathing type plots code--needs cleaning######
###uses data from rbt, rebreathing type table
###reorder rbt###
rownames(rbt) <- rbt$species ##make species names row identifiers
rbt <- rbt[order_p,] ##order df so that names will match up
rbt$sr_se <- sqrt((rbt$sustained*(1-rbt$sustained))/rbt$n)
rbt$ir_se <- sqrt((rbt$incidental*(1-rbt$incidental))/rbt$n)

###generate dataframe for plotting habitat/rb type pairs
stp_d <- data.frame(species=rep(rbt$species,2), habitat=rep(rbt$habitat,2), type=c(rep("sr", 20), rep("ir", 20)),prop=c(rbt$sustained,rbt$incidental), se=c(rbt$sr_se, rbt$ir_se))
stp_d$ccat <- paste(stp_d$habitat,stp_d$type)
stp_d$species <- factor(stp_d$species, levels = order_p)

#####rebreathing extent plot (SR,IR)
##top left plot
pl<-ggplot(data=stp_d, aes(x=species,y=prop,fill=type)) +
  geom_bar(position="stack",stat="identity")+
  labs(x = "Species", y = "Prop. of indiv.", fill = "Rebreathing frequency")

pl <- pl + scale_fill_manual(labels = c("Rebreathing (any)", "Sustained rebreathing"),values=c("gray","blue"))+theme_bw()

pl <- pl+ theme(text = element_text(size=20), 
                axis.text.x = element_text(angle=90, hjust=1, margin=margin(15,0,0,0)))+
  ylim(0,1)
pl <- pl +theme(axis.title.x=element_blank(),
                axis.text.x=element_blank(),
                plot.margin = unit(c(0.25,0,0.25,0), "cm"))
####uncomment to remove legend
# pl <- pl +theme(legend.position = "none")

# ###uncomment to remove y-axis labels##
# pl <- pl +theme(axis.text.y=element_blank(),
#                 axis.title.y=element_blank())

####hab vioplot for pl
##top right plot
plh<-ggplot(data=stp_d, aes(x=ccat,y=prop,fill=type)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1)+
  labs(x = "Habitat", y = "Prop. of indiv.", fill = "Rebreathing frequency")

plh <- plh + scale_fill_manual(labels = c("Rebreathing (any)", "Sustained rebreathing"),values=c("gray","blue"))+theme_bw()
plh <- plh + scale_x_discrete(labels=c("Aquatic IR", "Aquatic SR","Non-aquatic IR", "Non-aquatic SR" ))
plh <- plh+ theme(text = element_text(size=20), 
                  axis.text.x = element_text(angle=90, hjust=1))+
  ylim(0,1)
plh <- plh +theme(axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.text.y=element_blank(),
                  axis.title.y=element_blank(),
                  plot.margin = unit(c(0.25,0,0.25,0), "cm"))
###uncomment to remove legend###
plh <- plh +theme(legend.position = "none")

###################
##rotated ggtree plot###
treeg <- ggtree(npt) + coord_flip()

###create list of plots####
plts <- list(pl, plh,pl1,pl1h,pl4,pl4h,pl2,pl2h,pl3, pl3h)

library(patchwork)
###edit plot text sizes, etc; general template###
addo <- theme(axis.title.y = element_text(face="bold", size=25),
              axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
              axis.text.y = element_text(face="bold", color="black", size=25),
              axis.text.x = element_text(face="italic", size=30),
              legend.text = element_text(face="bold", color="black", size=20),
              axis.ticks = element_line(size=2))

plt_e <- lapply(plts, FUN = "+", addo)### add general template to all plots

##separate out plots with and without y axis labels, apply label removals
no_lab<- plt_e[c(2,4,6,8)]
##define theme for plots with no labels
add2 <- theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              plot.margin = unit(c(0.25,0,0.25,0), "cm"))

no_lab <- lapply(no_lab, FUN = "+", add2)

#separate out x label only plots
x_only <- plt_e[c(10)]

##define theme for plots with x labels only
add2 <- theme(axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              axis.title.x=element_blank(),
              axis.text.x = element_text(margin=margin(15,0,0,0)),
              plot.margin = unit(c(0.25,0,0.25,0), "cm"))

##apply x only theme
x_only <- lapply(x_only, FUN = "+", add2)

#separate out y and x label plots
y_ax <- plt_e[c(9)]

##define y and x label plot theme
add2 <- theme(axis.title.x=element_blank(),
              axis.text.x = element_text(margin=margin(15,0,0,0)),
              plot.margin = unit(c(0.25,0,0.25,0), "cm"))
##apply y and x label  plot theme
y_ax <- lapply(y_ax, FUN = "+", add2)

#separate out y label only plots
y_only <- plt_e[c(1,3,5,7)]

##define y label only theme
add2 <- theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              plot.margin = unit(c(0.25,0,0.25,0), "cm"))

##apply y label only theme
y_only <- lapply(y_only, FUN = "+", add2)

##best performance version##
##print using patchwork
#uncomment lines below to output figure 2 to file
#png("fig2.png", width = 21, height = 28, units = "in", res=72) ###image version
#pdf("fig2.pdf", width = 21, height = 28) ###pdf version
y_only[[1]]+no_lab[[1]]+y_only[[2]]+no_lab[[2]]+y_only[[3]]+no_lab[[3]]+y_only[[4]]+no_lab[[4]]+
  y_ax[[1]]+x_only[[1]]+treeg+ plot_layout(ncol=2, guides = "collect")
##uncomment to output as pdf
#dev.off()
