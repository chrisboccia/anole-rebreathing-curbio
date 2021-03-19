####run this block of code to define a function that is used below
###to summarize a behavioural variable of interest, run a phylo ANOVA
###and add it to a cumulative data frame for later plotting###
rb_assess <- function(indv,var_num,var_name,sp_df,spec,hab){
  library(nlme)
  library(caper)
  library(phytools)
  
  #aggregate data by species for the specified variable of interest
  ag_n <- aggregate(indv[,var_num], list(species=indv$Species), FUN=function(x) c(mean=mean(x, na.rm = TRUE), sd=sd(x, na.rm=TRUE),n=sum(!is.na(x)), se=sd(x, na.rm=TRUE)/sqrt(sum(!is.na(x)))))
  
  #convert to dataframe
  tem <- data.frame(ag_n$x)
  
  #combine with species names, add colnames
  ag_n <- cbind(ag_n$species, tem)
  colnames(ag_n) <- c("species", colnames(tem))
  
  ##make list of species with at least 3 individuals (argument passed by user)
  list_keep <- spec 
  hab_vec <- hab
  rownames(ag_n) <- ag_n$species
  
  ##cut down data to just species that make n>2 filter
  ag_n <- ag_n[list_keep,]
  
  ###add habitat
  ag_n$habitat <- as.factor(hab_vec)
  fphy <- ag_n[!is.na(ag_n$mean),]
  fphy <- fphy[fphy$n >2,]
  print(fphy)
  
  ###load tree for phylo anova, trim to data###
  poetree <- read.tree("Poe2017timetree_sis_ed.tre")
  nc <- name.check(poetree, fphy, data.names = fphy$species)
  npt <- drop.tip(poetree, tip=nc$tree_not_data)
  
  ###Brownian motion pgls, anova
  ##generate comparative data object for caper analysis
  compphy <- comparative.data(npt, fphy, names.col = 'species')
  
  ##run phylo anova on variable passed to function; set lambda using ML
  pglsModel <- pgls(formula=mean~habitat, data = compphy, lambda= 1)
  
  #lines below run an alternate PGLS formulation that uses nlme and ape
  # pglsModel<-gls(mean~habitat, data=fphy, correlation=
  #                      corPagel(value=1, phy=npt, fixed=FALSE, 
  #                               form = ~species), method="ML")
  
  ##output anova results
  print(anova(pglsModel))
  #print(summary(pglsModel))
  
  ##output lambda (and lambda ML plot)
  #lm.lk<-pgls.profile(pglsModel, which="lambda")
  #plot(lm.lk)
  plot(mean~habitat, data = ag_n)
  
  ##select columns, colnames to add to df
  tcol <- colnames(ag_n)
  tcol <- tcol[c(2,3,5)]
  
  ##add columns to df (passed to function)
  to_add <- ag_n[,c(2,3,5)]
  tcol<- var_endings(tcol, var_name)
  colnames(to_add) <- tcol 
  sp_df <- cbind(sp_df, to_add)
  
  ##return updated df
  return(sp_df)
}
