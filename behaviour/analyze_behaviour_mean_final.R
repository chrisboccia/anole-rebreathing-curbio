####This script analyzes and plots the behavioural data#####
####based on averaged individual performances across trials#
####see 'analyze_behaviour.R' for the in-text figures#######


#############################
#load functions and packages#
#############################

###run line below if using RStudio to setwd to folder with script####
if(Sys.getenv("RSTUDIO") == "1")
	setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##otherwise, setwd to suit your placement of the files##
#setwd("Your path")

####load packageVs####
library(geiger)
library(ape)
library(phytools)
library(ggplot2)
library(nlme)
library(dplyr)
library(car)
library(caper)
####

####define functions####
source("functions.R")
source("rb_assess.R")
#####End of functions

############
#load data##
############

##this script analyzes a trial_averaged filtered data set

#load mean data (pre-filtered)
#can run 'generate_indv' to recreate this file with different filters
indv <- read.csv("mean_individual.csv", header=TRUE)

##load habitat associations for each individual
ht <- read.csv("habitats.csv", header = TRUE)

################
#start analysis#
################

####trial mean version###
##the rb_assess function does the tasks laid out below;
##this block of code initializes the multipurpose dataframe
##that will contain all variables and metrics of interest

###generate aggregated max reinhales observed for each species
ag_n <- aggregate(indv$num_all_rb.mean, list(species=indv$Species), FUN=function(x) c(max=max(x), sd=sd(x),n=length(x), se=sd(x)/sqrt(length(x))))

##dataframe formatting
tem <- data.frame(ag_n$x)
ag_n <- cbind(ag_n$species, tem)
colnames(ag_n) <- c("species", colnames(tem))

##make list of species with at least 3 individuals
list_keep <- ag_n$species[ag_n$n >2]

###make habitat vector to match species

rownames(ag_n) <- ag_n$species
###this species list will be carried through all subsequent analyses###
ag_n <- ag_n[list_keep,]##cut down data to just species that make n>2 filter

###add habitat
hab_vec <- c("Non-aquatic", rep("Aquatic", 2), rep("Non-aquatic",5), "Aquatic", "Non-aquatic", "Aquatic",
             rep("Non-aquatic", 2), "Aquatic", rep("Non-aquatic", 6))
ag_n$habitat <- as.factor(hab_vec)

###initial setup for combined data frame for plotting###
##add species, habitat vectors
sp_df <- cbind(as.character(ag_n$species), as.character(ag_n$habitat),as.numeric(ag_n$n))
colnames(sp_df) <- c("species", "habitat", "n")
tcol <- colnames(ag_n)
tcol <- tcol[c(2)]

#get max rb (for reference; this analysis was already done in the best performance version)
to_add <- ag_n[,c(2)]
tcol<- var_endings(tcol, "rb_max")
sp_df <- data.frame(sp_df, rb_max=to_add)


###load phylo tree to conduct phylogenetic ANOVA
poetree <- read.tree("Poe2017timetree_sis_ed.tre")

###trim tree to focal taxa###
nc <- name.check(poetree, ag_n, data.names = ag_n$species)
npt <- drop.tip(poetree, tip=nc$tree_not_data)

##recheck tree, names
name.check(npt, ag_n)


###process/analyze additional variables, add to dataframe for plotting
###summarizes individual results by species, calculates se,sd,n, etc
###also runs caper pgls (see lines in rb_assess function) and outputs results
###variable names in quotes; can add additional variables at the end if desired
###argument format is: df_w_data, column_of_interest, colname_for_new_variable
sp_df <- rb_assess(indv,11,"rb_count_avg", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,12,"rb_rate_avg", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,23,"gular_rate_avg", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,16,"rb_period_avg", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,21,"br_rate_avg", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,29,"ex_dur_avg", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,30,"soh_avg", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,31,"toh_avg", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,32,"hbn_avg", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,27,"exhale_pct", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,10,"duration", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,6,"svl", sp_df, list_keep, hab_vec)
sp_df <- rb_assess(indv,7,"mass", sp_df, list_keep, hab_vec)

############
#figure S4g#
############
###initial ggplot
plex<-ggplot(data=sp_df, aes(x=habitat,y=`mean duration`,fill=habitat)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(fill="white", color="black", width=0.15)+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=3, fill="white", binwidth = 5)+
  labs(x = "Habitat", y = as.character("Average trial duration (s)"), fill = "Habitat")

##add labels, adjust ylim
plex <- plex + scale_fill_manual(labels = c("Semi-aquatic", "Non-aquatic"),values=c("#5DADE2", "#E59866"))+theme_bw()+
  scale_x_discrete(labels=c("aquatic", "non-aquatic")) +
  ylim(c(0, 400))

###uncomment to output plot image as png
#pdf("S4g.pdf", width = 8, height = 8)
plex
#dev.off()

####End figure S4g mean#####


####assess species average duration/habitat using PGLS####
dur_syn <- data.frame(sp_df$species,sp_df$habitat, sp_df$`mean duration`, sp_df$`mean svl`,sp_df$`mean mass`,sp_df$`mean rb_count_avg`, sp_df$`mean rb_rate_avg`, sp_df$`mean gular_rate_avg`)
colnames(dur_syn) <- c("species", "habitat", "duration", "svl", "mass","rbc","rbr","gulr")
rownames(dur_syn) <- dur_syn$species

compphy <- comparative.data(npt, dur_syn, names.col = 'species')

#rebreathing rate
# pglsModel <- pgls(formula=rbr~habitat*svl*log(mass), data = compphy, lambda=1)#covariate version (results equivalent)
pglsModel <- pgls(formula=rbr~habitat, data = compphy, lambda=1)
anova(pglsModel)

#rebreathing count
#pglsModel <- pgls(formula=rbc~habitat*svl*log(mass), data = compphy, lambda=1)#covariate version (results equivalent)
pglsModel <- pgls(formula=rbc~habitat, data = compphy, lambda=1)

anova(pglsModel)

#gular pumping
#pglsModel <- pgls(formula=gulr~habitat*svl*log(mass), data = compphy, lambda=1)#covariate version (results equivalent)
pglsModel <- pgls(formula=gulr~habitat, data = compphy, lambda=1)
anova(pglsModel)

#duration (stats reported in caption for fig. S4g)
#pglsModel <- pgls(formula=duration~habitat*svl*log(mass), data = compphy, lambda=1)#covariate version (results equivalent)
pglsModel <- pgls(formula=duration~habitat*log(mass), data = compphy, lambda=1)
anova(pglsModel)

###########
#Figure S4#
#(a-f)#####
###########

####load packages####
require(ggplot2)
require(ggtree)

###Assemble all panels for panel plot###

###create phylogenetic tree plot for panel, get x labels for other plots#####
treeg <- ggtree(npt) +geom_tiplab(align = TRUE)
ord_c <- data.frame(treeg$data[7], treeg$data[4])
#get phylo tip labels (so that plots can be put in the same order)
ord_c <- ord_c[1:20,]
ord_c <- ord_c[order(ord_c$y, decreasing=TRUE),]
order_p <- ord_c$label

##make species names row identifiers
rownames(sp_df) <- sp_df$species

##order df so that names will match up
sp_df <- sp_df[order_p,]

##necessary so that species will be in tree order
sp_df$species <- factor(sp_df$species, levels=order_p)

##left side, bottom plot above phylogeny##
####avg gular species plot####
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

##uncomment to remove legend
pl3 <- pl3 +theme(legend.position = "none")

pl3 <- pl3 + coord_cartesian(ylim = c(0,0.35))
#####

##bottom right plot##
####avg gular rate habitat plot (violin)
pl3h<-ggplot(data=sp_df, aes(x=habitat,y=`mean gular_rate_avg`,fill=habitat)) +
  geom_violin(trim=FALSE)+
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

##top left plot###
####avg rebreathing rate species plot####
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

##comment to remove legend#
pl2 <- pl2 +theme(legend.position = "none")

##run to remove y-axis labels##
# pl2 <- pl2 +theme(axis.text.y=element_blank(),
#                   axis.title.y=element_blank())
##

pl2 <- pl2 + coord_cartesian(ylim = c(0,0.3))
#####

##top right plot##
####avg rebreathing rate habitat plot (violin)
pl2h<-ggplot(data=sp_df, aes(x=habitat,y=`mean rb_rate_avg`,fill=habitat)) +
  geom_violin(trim=FALSE)+
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

##middle left plot##
####avg rebreathing count species plot####
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

##middle right plot###
####avg rebreathing count habitat plot (violin)
pl1h<-ggplot(data=sp_df, aes(x=habitat,y=`mean rb_count_avg`,fill=habitat)) +
  geom_violin(trim=FALSE)+
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

###################
##rotated ggtree plot###
treeg <- ggtree(npt) + coord_flip()

###create list of plots####
plts <- list(pl1,pl1h,pl2,pl2h,pl3,pl3h,plex)

library(patchwork)
###edit plot text sizes, etc; general template###
addo <- theme(axis.title.y = element_text(face="bold", size=25),
              axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
              axis.text.y = element_text(face="bold", color="black", size=25),
              axis.text.x = element_text(face="italic", size=30),
              legend.text = element_text(face="bold", color="black", size=20),
              axis.ticks = element_line(size=2))

### add general template to all plots using lapply
plt_e <- lapply(plts, FUN = "+", addo)

#remove legend from bottom right corner plot
plt_e[[7]] <- plt_e[[7]] + theme(legend.position = "none")
##separate out plots with and without y axis labels, apply label removals

##identify plots with no labels
no_lab<- plt_e[c(2,4,6)]

##generate theme with no labels
add2 <- theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.title.y=element_blank(),
              plot.margin = unit(c(0.25,0,0.25,0), "cm"))

##apply no label theme
no_lab <- lapply(no_lab, FUN = "+", add2)

#apply x-axis label only theme
# x_only <- plt_e[c(10)]
# add2 <- theme(axis.text.y=element_blank(),
#               axis.title.y=element_blank(),
#               axis.title.x=element_blank(),
#               axis.text.x = element_text(margin=margin(15,0,0,0)),
#               plot.margin = unit(c(0.25,0,0.25,0), "cm"))
# x_only <- lapply(x_only, FUN = "+", add2)

#apply y and x axis labels theme
y_ax <- plt_e[c(5,7)]
add2 <- theme(axis.title.x=element_blank(),
              axis.text.x = element_text(margin=margin(15,0,0,0)),
              plot.margin = unit(c(0.25,0,0.25,0), "cm"))
y_ax <- lapply(y_ax, FUN = "+", add2)

#apply y-axis only label theme
y_only <- plt_e[c(1,3)]
add2 <- theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              plot.margin = unit(c(0.25,0,0.25,0), "cm"))

y_only <- lapply(y_only, FUN = "+", add2)

##print figure S2##
##mean version##

#design for plot output
design <- "
AABB
CCDD
EEFF
GGHH"

#uncomment lines below to print to file
#pdf("figS4_stitched.pdf", width = 21, height = 28)
y_only[[1]]+no_lab[[1]]+y_only[[2]]+no_lab[[2]]+
  y_ax[[1]]+no_lab[[3]]+treeg+y_ax[[2]]+ plot_layout(design=design, guides = "collect")
#
#uncomment below to print to file
#dev.off()
