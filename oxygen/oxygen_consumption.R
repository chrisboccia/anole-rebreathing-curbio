#####add Seymour data/plastron supply/oxygen demand lines####
get_seymour_plot <- function(){
  #setwd to directory script is in
  if(Sys.getenv("RSTUDIO") == "1")
	setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  
  ##read in csv of oxygen consumption data
  ##see rbmath sheet for calculation cells
  ox <- read.csv("oxygen_consumption.csv", header = TRUE)
  
  #split data by source
  ox_s <- ox[ox$src=="seymour",]
  
  ox_b <- ox[ox$src=="boccia",]
  
  ###take individual average of nmol (masses are the same at all measurements but are included for code simplicity)
  ox_ib <- aggregate(cbind(mass,nmol,cham_upr_nmol,gecko_upr_nmol,tubin_lwr_nmol)~indv, ox_b, FUN = "mean")
  ox_ib$species <- c(rep("oxylophus",6), rep("lynchi",1), 
                     rep("aquaticus",5), rep("maculigula", 2),
                     rep("barkeri", 7),"maculigula",rep("barkeri", 2))
  
  ##aggregate seymour data
  ox_is <- aggregate(cbind(mass,nmol)~indv, ox_s, FUN = "mean")
  ox_is$species <- c("riolus_cupreus","elmis_maugei","phytobius_velatus", rep("macroplea_mutica",2),
                     "agraptocorixa_eurynome","aphelocheirus_aestivalis","dytiscus_marginalis","hydrous_piceus")
  
  ##take species average of mass, nmol
  ox_is <- aggregate(cbind(mass,nmol)~species, ox_is, FUN = "mean")
  
  ox_ib <- aggregate(cbind(mass,nmol,cham_upr_nmol,gecko_upr_nmol,tubin_lwr_nmol)~species, ox_ib, FUN = "mean")
  
  ##convert seymour mass values to g
  ox_is$mass <- ox_is$mass/1000
  
  #log original vars
  ox_ib$log_nmol <- log(ox_ib$nmol,10)
  ox_ib$log_mass <- log(ox_ib$mass,10)
  
  ox_is$log_nmol <- log(ox_is$nmol,10)
  ox_is$log_mass <- log(ox_is$mass,10)
  
  ##generate linear models##
  model_sey <- lm(log_nmol~log_mass , data = ox_is)
  model_boc <- lm(log_nmol~log_mass , data = ox_ib)
  
  ##can be used for model line comparisons##
  # mod_met <- lm(log(nmol)~log(mass), data = dfmeta)
  # 
  # ##use to check if any slopes of interest are significantly different
  # fit1 <- model_met
  # s1 <- summary(fit1)$coefficients
  # fit2 <- model_boc
  # s2 <- summary(fit2)$coefficients
  # db <- (s2[2,1]-s1[2,1])
  # sd <- sqrt(s2[2,2]^2+s1[2,2]^2)
  # df <- (fit1$df.residual+fit2$df.residual)
  # td <- db/sd
  # 2*pt(-abs(td), df)
  
  ##Generate prediction intervals (not plotted) 
  pred_sey <- predict(model_sey, interval = "prediction")
  pred_d_sey <- cbind(ox_is, pred_sey)
  
  pred_boc <- predict(model_boc, interval = "prediction")
  pred_d_boc <- cbind(ox_ib, pred_boc)
  
  #back transform prediction intervals (not plotted in current version)
  pred_d_sey$lwr <- 10^pred_d_sey$lwr
  pred_d_sey$upr <- 10^pred_d_sey$upr
  
  pred_d_boc$lwr <- 10^pred_d_boc$lwr
  pred_d_boc$upr <- 10^pred_d_boc$upr
  
  ##create respiration/provisioning lines from Seymour and Matthews 2013##
  
  #seymour actual data (mass in g, O2 consump. in nmol/s)
  
  #demand line
  dfmeta <- data.frame(mass=c(0.0001,10),nmol=c(0.00205959408648516,25.9287533345715))
  
  #plastron supply, 20um boundary layer thickness
  df20 <- data.frame(mass=c(0.0001,10),nmol=c(0.014658748434923,31.5813161406487))
  
  #plastron supply, 80 um boundary layer thickness
  df80 <- data.frame(mass=c(0.0001,10),nmol=c(0.00366468710873076,7.89532903516219))
  
  #plastron supply, 200 um boundary layer thickness
  df200 <- data.frame(mass=c(0.0001,10),nmol=c(0.0014658748434923,3.15813161406487))
  
  #plastron supply, 1000 um boundary layer thickness
  df1000 <- data.frame(mass=c(0.0001,10),nmol=c(0.000293174968698461,0.631626322812975))
  
  #plot using ggplot#
  
  #load package#
  library("ggplot2")
  
  #commented lines below remove all plastron lines except 20um
  #prediction intervals for seymour data, linear model fit/CI for aquatic anole data are also not plotted
  options(scipen=1000)
  
  #main plot with points, CIs
  pcon <- ggplot(pred_d_sey , aes(mass, nmol)) +
    geom_point(size=15) +
    #stat_smooth(method = lm, color="black") + 
    geom_point(size=13,data=pred_d_boc,fill="#4f66f0", pch=21,color="#233175")
    #stat_smooth(data=pred_d_boc, method = lm, color="red")
    #geom_point(data=ox_ib, aes(x=mass,y=cham_upr_nmol), color="gray") +
    #geom_point(data=ox_ib, aes(x=mass,y=gecko_upr_nmol), color="yellow") +
    #geom_point(data=ox_ib, aes(x=mass,y=tubin_lwr_nmol), color="purple")
    
    # Add demand, plastron lines intervals
  pcon <- pcon + 
    #stat_smooth(data=pred_d_sey, aes(y=lwr,x=mass), method=lm, se=FALSE, color = "black", linetype = "dashed", fullrange=T)+
    #stat_smooth(data=pred_d_sey, aes(y=upr,x=mass), method=lm, se=FALSE, color = "black", linetype = "dashed", fullrange=T)  +
    #geom_line(data=pred_d_boc, aes(y=upr), color = "red", linetype = "dashed")  +
    #geom_line(data=pred_d_boc,aes(y=lwr), color = "red", linetype = "dashed")  +
    stat_smooth(data=df20,aes(x=mass,y=nmol), method=lm, se=FALSE, color="blue", linetype="dotted",fullrange=T,size=3)+
    #stat_smooth(data=df80,aes(x=mass,y=nmol), method=lm, se=FALSE, color="blue", linetype="dotted",fullrange=T)+
    #stat_smooth(data=df200,aes(x=mass,y=nmol), method=lm, se=FALSE, color="blue", linetype="dotted",fullrange=T)+
    #stat_smooth(data=df1000,aes(x=mass,y=nmol), method=lm, se=FALSE, color="blue", linetype="dotted",fullrange=T)+
    stat_smooth(data=dfmeta,aes(x=mass,y=nmol), method=lm, se=FALSE, color="green", linetype="solid", fullrange=T,size=3)+
    scale_x_continuous(trans='log10', limits = c(0.0001,100)) + 
    scale_y_continuous(trans='log10')+
    labs(x="Mass (g)", y=expression("Oxygen consumption (nmol/s)"))
  
  #uncomment these lines to plot this graph only, output to file#
  #pdf("seymour_boccia_oxygen_con.pdf")
  #print(pcon)
  #dev.off()
  return(pcon)
}

