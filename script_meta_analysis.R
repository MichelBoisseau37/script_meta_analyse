rm(list=ls())
options (max.print = 5000)
library(leaflet)
require("ggrepel")
require(reshape)
library(cowplot)
require(dplyr)
require(vegan)
require(Hmisc)
require(pvclust)
require(ade4)
require(nlme)
require(gplots)
require(lme4)
require(plyr)
require(lattice)
require(lsmeans)
require(car)
require(reshape2)
library(forcats)
require(FactoMineR)
require(factoextra)
require(viridis)
library(RColorBrewer)
library(tidyr)
library(tidygraph)
library(ggraph)
library(ecoCopula)
require(funrar)
require(cooccur)
require(mvabund)
require(lme4)
require(fitdistrplus)
setwd("/home/mboisseau/Bureau/metaanalyse")
# sessionInfo()
# #loaded via a namespace (and not attached):
#   #   [1] bitops_1.0-6         tools_3.6.3          backports_1.1.4      R6_2.4.0             rpart_4.1-15         KernSmooth_2.23-17   lazyeval_0.2.2      
#   # [8] mgcv_1.8-31          colorspace_1.4-1     nnet_7.3-14          withr_2.1.2          tidyselect_1.1.0     gridExtra_2.3        curl_4.2            
#   # [15] compiler_3.6.3       htmlTable_1.13.2     flashClust_1.01-2    sandwich_2.5-1       caTools_1.17.1.2     scales_1.0.0         checkmate_1.9.4     
#   # [22] mvtnorm_1.0-11       stringr_1.4.0        digest_0.6.21        foreign_0.8-76       minqa_1.2.4          rio_0.5.16           base64enc_0.1-3     
#   # [29] pkgconfig_2.0.3      htmltools_0.4.0      htmlwidgets_1.3      rlang_0.4.6          readxl_1.3.1         rstudioapi_0.10      zoo_1.8-6           
#   # [36] gtools_3.8.1         acepack_1.4.1        zip_2.0.4            magrittr_1.5         leaps_3.0            Rcpp_1.0.2           munsell_0.5.0       
#   # [43] abind_1.4-5          lifecycle_0.2.0      scatterplot3d_0.3-41 stringi_1.4.3        multcomp_1.4-10      MASS_7.3-51.6        grid_3.6.3          
#   # [50] parallel_3.6.3       gdata_2.18.0         crayon_1.3.4         haven_2.3.1          splines_3.6.3        hms_0.5.3            knitr_1.25          
#   # [57] pillar_1.4.4         boot_1.3-25          estimability_1.3     codetools_0.2-16     glue_1.3.1           latticeExtra_0.6-28  data.table_1.12.2   
#   # [64] vctrs_0.3.1          nloptr_1.2.1         cellranger_1.1.0     gtable_0.3.0         purrr_0.3.4          assertthat_0.2.1     xfun_0.9            
#   # [71] openxlsx_4.1.0.1     xtable_1.8-4         coda_0.19-3          tibble_3.0.1         cluster_2.1.0        TH.data_1.0-10       ellipsis_0.3.0  

rm(list=ls())

# we create the dataset with abundance and prevalence tables

abu=read.csv("dataabu_brut0607.csv", header=TRUE,sep=";", dec=",")
head(abu)
# author yearsPubli     continent      country                   area Nb_Horse types_of_horses
# 2              Anjos       2006 South_america       brazil             Seropedica       33        domestic
# 
# breed       age         condition    method    CLR   CLA ASY HYB    LON BID  GOL    CAL    MIN   POC
# 2                      unk       unk               unk  necropsy    9.6  18.9   0   0  342.7   7   51   48.8   80.5   0.1
# .0
# ALV    TET  CAT  COR   PAT MON   LBI    LBR    NAS    LEP BRE   INS ELO TRO RAD AUR   ASH   ULT PZI PMB    CBI  PTM
# 2    0  371.4    0    0  48.2   0  18.9    9.6  447.5  171.5   0   0.5   0   0 1.2   0   5.2   0.2   0   0    4.7  0.2
# 
# CAP MNR TZI NIP    CCR   PEU  TGB    SER   EDE ACU   EQU TBV   TEN   VUL SKR SAG ADE PKR CLO MUC ORB
# 2    2.2 0.9   0   0    2.7   0.8    0    0.4   6.0 0.1   0.8 0.0   0.2   0.2   0   0   0   0   0   0   0

prev=read.csv("data_preva.csv",header=TRUE,sep=";", dec=",")
prev=prev[which(prev$author !="traversa"),]
prev=prev[which(prev$author !="stancampiano"),]
prev=prev[which(prev$author !="barus"),]
head(prev)
# study_ID  author yearsPubli     continent country       area Nb_Horse types_of_horses             breed  age condition
# 1       33   Silva       1999 South_america  brazil     brazil       36             unk               unk  unk       unk
# 
# method   CLA   CLR ASY HYB   LON BID SKR   GOL   CAL   MIN   POC  ALV TET    CAT   COR   PAT MON  LBI LBR SAG TRO
# 1  necropsy 41.67 19.44   0   0 100.0   0   0 91.67 80.56 88.89 13.89 2.78   0  94.44 22.22 66.67   0  0.0   0   0   0
# 
# NAS   LEP   BRE   INS  ELO  RAD AUR   ASH ADE   ULT PKR   PZI   PMB   CBI PTM   CAP CLO  EDE  EQU   VUL ACU MUC  SER
# 1  97.22 69.44 19.44 30.55  0.0 25.0   0 52.78   0 11.11   0  8.33  2.78 27.78   0 44.44   0  0.0  0.0   0.0   0   0  0.0
# 
# TEN MNR TZI  TBV NIP   PEU TGB ORB CCR
# 1   0   0   0  0.0 0.0 27.78   0   0   0

# we removed studies present in abu and prev tables
abu$authorY=paste(abu$author, abu$year, sep="_")
abu=abu[which(abu$authorY !="Silva_1999"),]
abu=abu[which(abu$authorY !="foster_1936"),]
abu=abu[which(abu$authorY !="Kornas_2011"),]
abu=abu[which(abu$authorY !="chapman_2003"),]
abu=abu[which(abu$authorY !="Torbert_1986"),]
abu=abu[which(abu$authorY !="Collobert_Laugier_2002"),]
abu=abu[which(abu$authorY !="Reinemeyer_1984"),]
abu=abu[which(abu$authorY !="kuzmina_2005"),]
abu=abu[which(abu$authorY !="Slivinska_2009"),]
abu=abu[which(abu$authorY !="Slivinska_2006"),]
abu=abu[which(abu$authorY !="Schankova_2015"),]
abu=abu[which(abu$authorY !="Bucknell_1995"),]
abu=abu[which(abu$authorY !="Morariu_2016"),]
abu=abu[which(abu$author !="stancampiano"),]
abu$authorY=NULL
abu$season=NULL
abu$commentaire=NULL

PaB=bind_rows(prev,abu)
rm(abu, prev)

PaB$CAS=NULL ### remove this species, bad name, or not strongyl
PaB$CLO=NULL
PaB$TZI=NULL
PaB$study_ID=seq(1:dim(PaB)[1])
dd=read.table('coordo_geo3.csv',sep=';',header=T)
df=merge(dd,PaB, by="area")

rm(dd)

dd=df[,c('area','lat','lon')]
doublons <- which(duplicated(dd[,c('area')]))
dd <-dd[-doublons,]
head(dd)
# area        lat       lon
# 1  Australia -35.288681 149.13183
# 2     brazil -14.235004 -51.92528
# 3   BTA_univ  49.794265  30.11198
# 4 canal_zone   8.956376 -79.53795
# 5    Chiapas  16.215600  93.58310
# 6   cracovie  50.065638  19.95818
clim=read.table(file='koppen_1901-2010.tsv',header=T)
#  longitude latitude p1901_2010
# 1   -179.75    71.25         ET
# 2   -179.75    68.75         ET
# 3   -179.75    68.25         ET
# 4   -179.75    67.75         ET
# 5   -179.75    67.25         ET
# 6   -179.75    66.75         ET
colnames(clim)=c("lon","lat","kop")
head(clim)
# lon   lat kop
# 1 -179.75 71.25  ET
# 2 -179.75 68.75  ET
# 3 -179.75 68.25  ET
# 4 -179.75 67.75  ET
# 5 -179.75 67.25  ET
# 6 -179.75 66.75  ET
clim$kop=as.character(clim$kop)


dd$kop=0
v=seq(1,dim(dd)[1])
temp_lat =0
temp_lon=0

for(i in v){
  temp_lat = dd$lat[i]
  temp_lon = dd$lon[i]
  temp = clim[which(clim$lon> temp_lon-2 & clim$lon< temp_lon+2 &
                      clim$lat> temp_lat-2 & clim$lat< temp_lat+2),]
  tt=array(0,dim(temp)[1])
  t=seq(1,dim(temp)[1])
  for(k in t){
    dlat=abs(temp_lat - temp$lat[k])
    dlon=abs(temp_lon - temp$lon[k])
    d=dlat+dlon
    tt[k] = d
  }
  dd$kop[i]=temp$kop[which.min(tt)]
}

## Add cimate category
df_eq=data.frame(kop = unique(substr(dd$kop,1,1)),
                 clim=c('Temperate','Tropical','Continental','Arid'))
dd$climate = df_eq$clim[match(substr(dd$kop,1,1),df_eq$kop)]
head(unique(dd[,c('area','kop', 'climate')]))
# area kop     climate
# 1  Australia Cfb   Temperate
# 2     brazil  Aw    Tropical
# 3   BTA_univ Dfb Continental
# 4 canal_zone  Am    Tropical
# 5    Chiapas  Am    Tropical
# 6   cracovie Dfb Continental
rm(d,dlat,dlon,i,k,t,temp_lat,temp_lon,tt,v,temp, df_eq,df, doublons)
dd$lat=NULL
dd$lon=NULL
df_abu= merge(x = PaB, y=dd, by=c('area')) ### merge by climatic data to add climatic and kop data 

rm(clim, PaB,dd)

df_abu$continent <- fct_recode(df_abu$continent, "Oceania"="Oceanie", "Africa"="South_africa","Africa"="South_Africa","Asia"="Asie", "America"="South_america",
                               "America"="North_america","America"="North_America")
df_abu$country <- fct_recode(df_abu$country, "Brazil"="brazil", "France"="france", "Italy"="italy", "Mexico"="mexique", "Panama"="panama",
                             "Poland"="poland", "Czech_republic"="republique_tchèque", "Romania"="romania", "South_africa"="south_africa",
                             "United_Kingdom"="uk", "Ukraine"="ukraine","Germany"="germany", "Kazakhstan"="kazakhstan","Lithuania"="lithuania",
                             "Panama"="panama", "Czech_republic"="republiquetcheque", "Romania"="romania","Russia"="russie","Ukraine"="ukraine")

##### Add clim2, some levels of kop didn t have enough data (2 studies or less)
lgt <- dim(df_abu)[1]
dtaNorm= array(NA,lgt)
for (i in 1:lgt){
  if(df_abu$kop[i]=="Am"|
     df_abu$kop[i]=="Aw"|
     df_abu$kop[i]=="BSh"|
     df_abu$kop[i]=="Csa"|
     df_abu$kop[i]=="As"|
     df_abu$kop[i]=="Csb") {dtaNorm[i]="Others"}
  if(df_abu$kop[i]=="Cfa") {dtaNorm[i]="Cfa"}
  if(df_abu$kop[i]=="Dfb") {dtaNorm[i]="Dfb"}
  if(df_abu$kop[i]=="Cfb") {dtaNorm[i]="Cfb"}
}
abuClim2=cbind(df_abu,dtaNorm) 
colnames(abuClim2)=c( "area","study_ID","author","yearsPubli","continent","country",
                      "Nb_Horse","types_of_horses" ,"breed","age","condition","method",
                      "CLA","CLR","ASY","HYB","LON","BID",
                      "SKR","GOL","CAL","MIN","POC","ALV",
                      "TET","CAT","COR","PAT","MON","LAB",
                      "LBR","SAG","NAS","LEP","BRE","INS",
                      "ELO","RAD","AUR","ASH","ADE","ULT",
                      "PKR","RAT","IMP","BIC","MET","CAP",
                      "EDE","EQU","VUL","ACU","MUC","SER",
                      "TEN","MNR","TBV","NIP","EUP","GOB",
                      "ORB","CCR","kop","climate","kop2")
rm(dtaNorm, i, lgt, df_abu)

# we removed useless variables
abuClim2$area=NULL
abuClim2$country=NULL
abuClim2$types_of_horses=NULL
abuClim2$breed=NULL
abuClim2$kop=NULL
abuClim2$condition=NULL
abuClim2$age=NULL
abuClim2$lat=NULL
abuClim2$lon=NULL

# CLA and CLR are the same species LBI and LBR corresponding to the names given actually
abuClim2$LBI <- rowSums(abuClim2[,c("LBI","CLA")],na.rm=TRUE)
abuClim2$LBR <- rowSums(abuClim2[,c("LBR","CLR")],na.rm=TRUE)
abuClim2$CLR=NULL
abuClim2$CLA=NULL

# refomat table with all species in only one column
abuClim22 <- melt(abuClim2, id.vars=c("study_ID","author","yearsPubli","continent","Nb_Horse",          
                                      "method","climate","kop2"),
                  variable.name="species", value.name="val")
rm(abuClim2)

# add presence/absence and remove quantitatives values
abuClim22$pab=ifelse(abuClim22$val > 0,1,0)
abuClim22$val=NULL

write.table(abuClim22, "PaB4.csv", sep=";", row.names=FALSE)


####################### plot studies on map #########################
rm(list=ls())
require(ggmap)
require(ggrepel)

co=read.csv(file = 'dtaForCGEO2.csv', header=T, sep=';')
head(co)
# area N study_ID                 val        country   method kop     climate       lat        lon
# 1  Australia 3       49      Australia(n=3)      Australia necropsy Cfb   Temperate -35.28868 149.131835
# 2     brazil 3        1         Brazil(n=3)         brazil necropsy  Aw    Tropical -14.23500 -51.925280
# 3     prague 1       25 Czech-republic(n=1) Czech-republic necropsy Cfb   Temperate  50.07237  14.446392
# 4  normandie 1       15         France(n=1)         france necropsy Cfb   Temperate  49.24262  -0.350946
# 5    germany 2        6        Germany(n=2)        germany necropsy Cfb   Temperate  52.52013  13.421293
# 6 kazakhstan 1        7     Kazakhstan(n=1)     kazakhstan necropsy Dfb Continental  51.16509  71.475951
clcol=viridis_pal(option='D')(50)
mp <- NULL

mapWorld <- borders("world",ylim = c(-50, 80),xlim=c(-180,180), colour="gray50", fill="white") # create a layer of borders

mp <- ggplot() + 
  mapWorld +
  theme_bw()

#Now Layer the cities on top

mp <- mp + geom_point(data=co,aes(x=lon, y=lat ,color=method), size=3) +
  xlab("Longitude") + 
  ylab("Latitude") +
  ggtitle("")+
  geom_label_repel(data=co,segment.alpha = 0,size=4,  
                   aes(x=lon, y=lat, col=method,label=val),show.legend = FALSE,max.overlaps = Inf) +
  scale_color_brewer(palette="Dark2")+
  labs(color = "Recovery method")+
  theme(legend.position=c(.1,.2),text=element_text(size=10))+ 
  theme_bw()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14), 
        legend.position="none")

mp
############################################################
##############  Plot prev abd average  #####################
############################################################

# we use two tables corresponding to the initial abundance and prevalence tables, 
# where the species are all put together in a single column (melt function)
# and the abundances are transformed into relative abundance

abu=read.csv("abdFinal.csv", header=TRUE,sep=";")
Prev=read.csv("prevalences_Finales.csv", header=TRUE,sep=";")
Prev$species <- fct_recode(Prev$species, "LBI"="CLA")
abu$species <- fct_recode(abu$species, "LBI"="CLA")
abu$species <- fct_recode(abu$species, "LBR"="CLR")
par(mfrow=c(1,1))


c=aggregate(prev ~ species , FUN=mean,data=Prev)
a=aggregate(prev ~ species , FUN=sd,data=Prev)
c=merge(c,a, by=c("species"))
colnames(c)=c("species", "prev", "sdP")
d=aggregate(abdRel ~ species, FUN=mean,data=abu)
f=aggregate(abdRel ~ species, FUN=sd,data=abu)
d=merge(d,f, by=c("species"))
colnames(d)=c("species","abd", "sda")
e=merge(c,d,by=c("species"))
e$species <- fct_recode(e$species, "BIC"="CBI","LAB"="LBI","EUP"="PEU","IMP"="PMB",
                        "MET"="PTM","RAT"="PZI","GOB"="TGB" )

nm=read.csv("annexe.csv", header=TRUE,sep=";")
nm$full=paste(nm$L,nm$Species, sep=".")
colnames(nm)=c("species", "genus","l", "sp", 'speciesF')
e=merge(e,nm, by="species")
e$l=NULL
rm(a,c,d,f, abu, Prev, nm)

lgt <- dim(e)[1]
dtaNorm= array(NA,lgt)
for (i in 1:lgt){
  if(e$prev[i] > 80) {dtaNorm[i]="Dominant"}
  if(e$prev[i] < 80 & e$prev[i] > 25) {dtaNorm[i]="Subdominant"}
  #if(e$prev[i] < 50 & e$prev[i] > 20) {dtaNorm[i]="Background"}
  if(e$prev[i] < 25) {dtaNorm[i]="Rare"}
}
e=cbind(e,dtaNorm)
colnames(e)=c("species","prev","sdP","abd","sda","genus","sp","speciesF", "grp" )
rm(dtaNorm,i,lgt)
revalue(e[e$speciesF=="C.insisgnis",], "C.insigne")

e$speciesF <- fct_recode(e$speciesF, "C.insigne"="C.insisgnis")
e$sdP=ifelse(e$prev<25,0,e$sdP)
e$sda=ifelse(e$abd<	0.035,0,e$sda)


# We make a figure by calculating the prevalences of species using the presence/absence of our dataset. 
PAB=read.csv("PaB4.csv", header=TRUE,sep=";") 
PAB$GrpYear=ifelse(PAB$yearsPubli > 1960, "After", "Before")
PAB$env = factor(paste(substr(PAB$continent,1,2),substr(PAB$clim,1,4),substr(PAB$meth,1,3),sep='-')) ## Factor to be tested too
colnames(PAB)=c( "study_ID","author","yearsPubli", "continent",  "nh",   "method","climate","kop2","species",
                 "pab","GrpYear","env" )

df= dcast(PAB, study_ID + author + yearsPubli  + continent + nh + method + climate + kop2 +      
            GrpYear + env ~ species, value.var="pab")

df=df[,-c(1:10)]
df2=as.data.frame(t(df))
df2$richness=rowSums(df2)
df2$sp=rownames(df2)

df2=df2[,-c(1:69)]
nstu=length(df[,1])

df2$prev2=df2$richness/nstu*100
colnames(df2)=c("eff" , "species" , "prev2")
rm(PAB, df, nstu)

e=merge(e,df2, by="species")
e=e[e$species!="MON",]
a1=ggplot(e, aes(x=prev, y=abd,label = speciesF,ymin=abd-sda, ymax=abd, xmin=prev-sdP, xmax=prev))+
  geom_point(aes(color=genus, size=15))+
  guides(colour = guide_legend(override.aes = list(size=7)))+
  #scale_color_manual(values=c("red2","royalblue3","chartreuse4"))+
  scale_color_viridis(discrete=TRUE,option="plasma", direction = -1) +
  geom_text_repel(size=5)+
  #  scale_y_continuous(limits=c(0, 0.6))+
  #scale_x_continuous(limits=c(0, 100))+
  xlab('Average of reported prevalences')+
  ylab('Average relative abundance')+
  labs(shape = "Type of species", color="Genus") +
  scale_shape_manual(values=c(3, 16, 17))+
  guides(size = "none")+
  geom_vline(xintercept = c(25,80), linetype="dashed", color = "red")+
  theme_bw()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14))

a2=ggplot(e, aes(x=prev, y=prev2,label = speciesF))+
  geom_point(aes(color=genus, size=15))+
  #scale_color_manual(values=c("red2","royalblue3","chartreuse4"))+
  scale_color_viridis(discrete=TRUE,option="plasma", direction = -1) +
  geom_text_repel(size=5)+
  guides(size = "none")+
  #  scale_y_continuous(limits=c(0, 0.6))+
  #scale_x_continuous(limits=c(0, 100))+
  xlab("Average of reported prevalences")+
  ylab("Observed worlwide prevalence")+
  labs(shape = "Type of species", color="Genus") +
  scale_shape_manual(values=c(3, 16, 17))+
  geom_vline(xintercept = c(25,80), linetype="dashed", color = "red")+
  theme_bw()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.position="none")

library(ggpubr)
ggarrange(mp,                # Première ligne contenant le line plot
          # Deuxième ligne avec box plots et dot plots
          ggarrange(a1, a2, ncol = 2, labels = c("B", "C")), 
          nrow = 2, 
          labels = "A" )      # Étiquette du line plot) 


mean(e[e$grp=="Rare",]$abd) #  0.003446516
mean(e[e$grp=="Dominant",]$abd) #  0.178484
mean(e[e$grp=="Subdominant",]$abd) # 0.04119527
sd(e[e$grp=="Rare",]$abd) # 0.002714834
sd(e[e$grp=="Dominant",]$abd) # 0.01481067
sd(e[e$grp=="Subdominant",]$abd) # 0.01924178

mean(e[e$grp=="Rare",]$prev) # 8.937745
mean(e[e$grp=="Subdominant",]$prev) # 54.965
mean(e[e$grp=="Dominant",]$prev) # 83.97327
sd(e[e$grp=="Rare",]$prev) # 6.701526
sd(e[e$grp=="Subdominant",]$prev) # 12.68996
sd(e[e$grp=="Dominant",]$prev) #  3.386756

mean(e[e$grp=="Rare",]$prev2) # 38.53999
mean(e[e$grp=="Subdominant",]$prev2) # 88.22464
mean(e[e$grp=="Dominant",]$prev2) #  97.10145
sd(e[e$grp=="Rare",]$prev2) # 22.30501
sd(e[e$grp=="Subdominant",]$prev2) # 15.16925
sd(e[e$grp=="Dominant",]$prev2) #  2.898551


e$grp=factor(e$grp)
kruskal.test(e$prev ~ e$grp)
# Kruskal-Wallis rank sum test
# 
# data:  e$prev by e$grp
# Kruskal-Wallis chi-squared = 23.381, df = 2, p-value = 8.375e-06
pairwise.wilcox.test(e$prev , e$grp,p.adj = "bonf")
# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  e$prev and e$grp 
# 
# Dominant Rare   
# Rare        0.0015   -      
#    Subdominant 0.0364   2.5e-07
# 
# P value adjustment method: bonferroni 
kruskal.test(e$abd ~ e$grp)

# Kruskal-Wallis rank sum test
# 
# data:  e$abd by e$grp
# Kruskal-Wallis chi-squared = 23.381, df = 2, p-value = 8.375e-06

pairwise.wilcox.test(e$abd , e$grp, p.adj = "bonf")

# Pairwise comparisons using Wilcoxon rank sum test 
# 
# data:  e$abd and e$grp 
# 
# Dominant Rare   
# Rare        0.0015   -      
#    Subdominant 0.0364   2.5e-07
# 
# P value adjustment method: bonferroni 

e$sdP=NULL
e$sda=NULL
e$grp=NULL
e$eff=NULL
e$sp=NULL
e$genus=NULL
colnames(e)=c("species_code", "prevalence1", "abundance", "species", "prevalence2")
write.table(e, "SupplementaryTab2.csv", sep=";", row.names=FALSE)

####################################################################################################################
par(mfrow=c(1,1)) 
rm(list=ls())
######---- let s see the construction of the dataset

ma = read.csv(file = './PaB4.csv', header=T, sep=';')
colnames(ma)=c("study_ID","author","year", "continent" , "nh"  , "meth","clim","kop2","species" ,"pab")
dt = ma %>% dcast(study_ID + year + clim + continent + nh + kop2 + meth ~ species, value.var = 'pab')
dt[is.na(dt)] = 0
dt$meth = factor(dt$meth)
dt = data.frame(dt)

### Contingency table
table(dt$continent,dt$meth)

# deworming necropsy
# Africa          0        3
# America         0       15
# Asia            0        1
# Europe         31       16
# Oceania         0        3

table(dt$continent,dt$kop)
# Cfa Cfb Dfb Others
# Africa    0   2   0      1
# America   6   2   0      7
# Asia      0   0   1      0
# Europe    1  13  33      0
# Oceania   0   2   0      1

table(dt$continent,dt$clim)
# Arid Continental Temperate Tropical
# Africa     0           0         3        0
# America    0           0         8        7
# Asia       0           1         0        0
# Europe     0          33        14        0
# Oceania    1           0         2        0

table(dt$clim, dt$meth)
# deworming necropsy
# Arid                0        1
# Continental        28        6
# Temperate           3       24
# Tropical            0        7

table(dt$continent[dt$meth=='necropsy'],
      dt$clim[dt$meth=='necropsy'])
# Arid Continental Temperate Tropical
# Africa     0           0         3        0
# America    0           0         8        7
# Asia       0           1         0        0
# Europe     0           5        11        0
# Oceania    1           0         2        0
dt$env = factor(paste(substr(dt$continent,1,2),substr(dt$clim,1,4),substr(dt$meth,1,3),sep='-')) ## Factor to be tested too
env5=names(table(dt$env))[table(dt$env)>=5]
env3=names(table(dt$env))[table(dt$env)>=3]

### Other option based on contclim
#table(dt$contclim[dt$meth!='necropsy']) 
# Africa Cfb  Africa Others    America Cfa    America Cfb America Others       Asia Dfb 
#          2              1              6              2              5              2 
# Europe Cfa     Europe Cfb     Europe Dfb    Oceania Cfb Oceania Others 
#          1              8             33              2              1 
#dt = dt[dt$meth=='necropsy' & 
#          dt$contclim %in% c('America Cfa','America Others','Europe Cfb','Europe Dfb'),]

## Filter conditions
dt = dt[dt$env %in% env5,] ## keeping Africa destabilizes NMDS and PCA

### Keep species with enough variation
colnames(dt)[2]='year'
#sp5 = names(which(colSums(dt[,9:(dim(dt)[2]-1)])>=10 & colSums(dt[,9:(dim(dt)[2]-1)])<50))
sp5 = names(which(colSums(dt[,9:(dim(dt)[2]-1)])>=2 & colSums(dt[,9:(dim(dt)[2]-1)])<53))

#sp5 = colnames(dt[,9:(dim(dt)[2]-1)])
env = dt[,c('study_ID','year','nh','env')]
sp = dt[,c(which(colnames(dt) %in% sp5))]
sp = data.frame(lapply(sp,as.integer))
dt2 = cbind(env,sp)
dt2 = data.frame(dt2)
dt2$env = factor(dt2$env)
ma2 = melt(dt2,1:4) ## for LMER
dt2$study_ID = NULL ## for PCA

### PCA
pc = dudi.mix(dt2,nf = 2,scannf=F)
s.corcircle(pc$co) ## ASH vs. year

####-- Some species are not found elsewhere
## Check that species names did not evolve
yeartrend = reshape2::melt(dt2[,c('env','year','ASH')],1:2)
ggplot(yeartrend,aes(y = year,x = factor(value),col=env,shape = variable, group = variable)) +
  geom_point(size = 3) #+ facet_wrap(~ variable, ncol=1)

### NMDS
sp = dt2[,-c(1:3)]
env = dt2[,1:3]


ord <- metaMDS(sp,try = 20,k=3, distance = 'jaccard')
fit <- envfit(ord, env, perm = 999)

scores(fit, "vectors")
#        NMDS1   NMDS2
# year  0.6110 0.08513
# nh   -0.5847 0.03168

### First two axes
#par(mfrow=c(1,1)) 
### Ellipsoid hulls : Method
plot(ord,type = 't')
with(dt, ordiellipse(ord, meth, kind = "ehull", label = TRUE))

plot(ord,type = 't')
with(dt, ordiellipse(ord, continent, kind = "ehull", label = TRUE))

### Ellipsoid hulls : CLIMATE
plot(ord,type = 't')
with(dt, ordiellipse(ord, clim, kind = "ehull", label = TRUE))

### Ellipsoid hulls : CONTINENT
plot(ord,type = 't')
with(dt, ordiellipse(ord, env, kind = "ehull", label = TRUE))

### Position year and nh
plot(ord,type = 't')
plot(fit)
with(dt, ordiellipse(ord, meth, kind = "ehull", label = TRUE))

### Axes 1,3
### First two axes
### Ellipsoid hulls : Method
plot(ord,type = 't',choices = c(3,2))
with(dt, ordiellipse(ord, meth, kind = "ehull", label = TRUE))

### Ellipsoid hulls : CLIMATE
plot(ord,type = 't',choices = c(3,2))
with(dt, ordiellipse(ord, clim, kind = "ehull", label = TRUE))

### Ellipsoid hulls : CONTINENT
plot(ord,type = 't',choices = c(3,2))
with(dt, ordiellipse(ord, continent, kind = "ehull", label = TRUE))


####-------======== Species test
sp5 = names(which(colSums(dt[,9:(dim(dt)[2]-1)])>=10 & colSums(dt[,9:(dim(dt)[2]-1)])<50))
env = dt[,c('study_ID','year','nh','env')]
sp = dt[,c(which(colnames(dt) %in% sp5))]
sp = data.frame(lapply(sp,as.integer))
dt3 = cbind(env,sp)
dt3 = data.frame(dt3)
dt3$env = factor(dt3$env)
ma3 = melt(dt3,1:4) ## for LMER
ma3$study_ID =factor(match(ma3$study_ID,levels(factor(ma3$study_ID))))

### Mixed model
m = glmer(value ~ nh + year + env*variable + (1|study_ID),
          family = binomial(link = 'logit'),
          control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
          data = ma3)
Anova(m,'III')
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: value
# Chisq Df Pr(>Chisq)    
# (Intercept)    0.0006  1    0.98094    
# nh             4.6430  1    0.03118 *  
#    year           0.0022  1    0.96277    
# env           12.9010  4    0.01177 *  
#    variable      21.4830 24    0.61009    
# env:variable 178.3020 96  6.925e-07 ***
#    ---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
s=data.frame(summary(m)$coefficients)
s[s$Pr...z..<0.05,]


rm(list=ls())
##### let s try on Europe x continental subset

PAB=read.csv("PaB4.csv", header=TRUE,sep=";") ## Big data set with pres/abs from prev and abd
PAB$GrpYear=ifelse(PAB$yearsPubli > 1960, "After", "Before")
PAB$env = factor(paste(substr(PAB$continent,1,2),substr(PAB$clim,1,4),substr(PAB$meth,1,3),sep='-')) ## Factor to be tested too
colnames(PAB)=c( "study_ID","author","yearsPubli", "continent",  "nh",   "method","climate","kop2","species",
                 "pab","GrpYear","env" )
head(PAB)
#  study_ID  author yearsPubli continent        nh  method     climate   kop2 species pab GrpYear         env
# 1       49  boxell       2004   Oceania       29  necropsy   Temperate    Cfb     CLA   0   After Oc-Temp-nec
# 2        1   Silva       1999   America       36  necropsy    Tropical Others     CLA   1   After Am-Trop-nec
# 3        2 kuzmina       2008    Europe       12 deworming Continental    Dfb     CLA   1   After Eu-Cont-dew
# 4        3  foster       1936   America       17  necropsy    Tropical Others     CLA   0  Before Am-Trop-nec
# 5       57  guiris       2010   America        2  necropsy    Tropical Others     CLA   1   After Am-Trop-nec
# 6        4  Kornas       2011    Europe       14 deworming Continental    Dfb     CLA   1   After Eu-Cont-dew


#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
### we selected species with some variation in order to facilitate the convergence of the model 
df=as.data.frame(table(PAB$species, PAB$pab))
colnames(df)=c("species", "pab", "freq")
df=df[df$pab==1,]
df=df[df$freq > 14,] # 10/80 % de prev? 62 et 7 55/14 
df=df[df$freq < 55,]

df$species=factor(df$species)
spe=df$species
# ACU ASH CAP CBI EDE ELO LBR PEU PMB POC PTM PZI RAD SER TBV ULT VUL
PAB2=subset(PAB, species %in% spe)
PAB2$species=factor(PAB2$species)
levels(PAB2$species)
rm(df, spe)
PAB2Obs=subset(PAB2, env %in% c("Af-Temp-nec", "Am-Temp-nec", "Am-Trop-nec" ,"As-Cont-nec" ,"Eu-Cont-dew" ,
                                "Eu-Cont-nec", "Eu-Temp-dew", "Eu-Temp-nec" ,"Oc-Arid-nec" ,"Oc-Temp-nec"))
PAB2analysis=subset(PAB2, env %in% c("Am-Temp-nec", "Am-Trop-nec" ,"Eu-Cont-dew" ,
                                     "Eu-Cont-nec", "Eu-Temp-dew", "Eu-Temp-nec" ))
mean(PAB[PAB$meth=="necropsy",]$nh)  # 36.84211
mean(PAB[PAB$meth=="deworming",]$nh) # 18.90323

# unique(PAB[, c("meth", "author", "nh")])
ggplot(PAB)+
  geom_point(aes(x=method, y=nh))+
  geom_boxplot(aes(x=method, y=nh))+
  theme_bw()

wilcox.test(nh~method, data=PAB)
# Wilcoxon rank sum test with continuity correction
# 
# data:  nh by method
# W = 645869, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0
s=data.frame(summary(glmer(pab ~ env*species+log(nh)+GrpYear+ (1|study_ID), family=binomial,data=PAB2analysis))$coefficients)
s[s$Pr...z..<0.05,]

# Estimate Std..Error   z.value     Pr...z..
# envEu-Temp-nec            -3.5270528  1.4860597 -2.373426 0.0176239254
# speciesBRE                -3.2343076  1.4667653 -2.205061 0.0274497828
# log(nh)                    0.8065656  0.2268852  3.554951 0.0003780500
# envEu-Cont-dew:speciesASH  6.6730034  1.7358643  3.844196 0.0001209481
# envAm-Trop-nec:speciesBRE  4.5996614  1.8893164  2.434564 0.0149097484
# envEu-Temp-nec:speciesBRE  4.6507139  1.9662804  2.365234 0.0180186657
# envEu-Temp-nec:speciesCBI  4.4230348  2.0036756  2.207460 0.0272819050
# envEu-Temp-nec:speciesELO  5.8783352  1.9582380  3.001849 0.0026834504
# envEu-Cont-dew:speciesSER  2.8189107  1.3439850  2.097427 0.0359557754
# envEu-Temp-nec:speciesULT  3.9152547  1.7290563  2.264388 0.0235502468
# envEu-Cont-dew:speciesVUL  3.0821013  1.3530183  2.277945 0.0227298521

df2=as.data.frame(table(PAB2Obs$species,PAB2Obs$pab, PAB2Obs$env))
colnames(df2)=c("species","abs" ,"env", "prop")
df2=subset(df2, species %in% c("LBR", "ASH", "CBI", "ELO", "SER", "ULT", "VUL"))
ggplot(data=df2, aes(x=env, y=prop, fill=abs)) +
  geom_bar(stat="identity", position="fill")+
  xlab("Presence or absence")+
  ylab("Proportion of presence and absence")+
  facet_wrap(~species)
rm(PAB2analysis,PAB2Obs, df2, s)


# we have to remove species without variation between deworming and necropsy

#### effet de la methode sur les donnees europeennes
DfEur=subset(PAB, env %in% c( "Eu-Cont-dew", "Eu-Cont-nec", "Eu-Temp-dew","Eu-Temp-nec" ))
DfEur$env=factor(DfEur$env)
DfEur$continent=factor(DfEur$continent)
DfEur$clim=factor(DfEur$clim)

ggplot(DfEur,aes(x=method, y=as.factor(pab)))+
  geom_point()+
  geom_jitter( height = 0.02)+
  theme_bw()+
  facet_wrap(~species)+
  xlab("Method")+
  ylab('Presence or absence')+
  ggtitle("")


DfEur2=subset(DfEur, clim=="Continental")
DfEur2$species=factor(DfEur2$species)
DfEur2$continent=factor(DfEur2$continent)
DfEur2$env=factor(DfEur2$env)
DfEur2$clim=factor(DfEur2$clim)
SpMeth = data.frame(table(DfEur2$species,DfEur2$pab,DfEur2$method))
ggplot(SpMeth,aes(x = Var1, y = Freq, group = paste0(Var3,Var1), fill = Var2))+
  geom_bar(stat='identity') +
  facet_wrap(~ Var3,ncol=1) + theme_bw()
manu=c("ACU","BRE", "EDE", "EQU", "HYB" , "NIP" , "PAT" , "PEU" , "SER" , "TBV" , "ULT" , "VUL")
dt = DfEur2 %>% dcast(study_ID + yearsPubli + climate + continent + nh + kop2 + method ~ species, value.var = 'pab')
sp5 = names(which(colSums(dt[,8:(dim(dt)[2]-1)])>=5 & colSums(dt[,8:(dim(dt)[2]-1)])<30))
DfEur3=subset(DfEur2, species %in% manu)

s=data.frame(summary(glmer(pab ~ method*species+log(nh)+(1|study_ID), family=binomial,data=DfEur3))$coefficients)
s[s$Pr...z..<0.05,]
# Estimate Std..Error   z.value     Pr...z..
# (Intercept)             -4.173056  0.8813005 -4.735111 2.189355e-06
# speciesASH               5.075813  1.2368119  4.103949 4.061571e-05
# speciesBRE              -1.502335  0.7403787 -2.029143 4.244368e-02
# speciesCLA               2.568006  0.7045219  3.645034 2.673568e-04
# speciesELO               1.591366  0.6421544  2.478167 1.320593e-02
# speciesLBR              -1.500195  0.7400950 -2.027030 4.265933e-02
# speciesPMB               1.383767  0.6346138  2.180487 2.922138e-02
# speciesPOC               1.595901  0.6423401  2.484510 1.297297e-02
# speciesPTM               1.383392  0.6346019  2.179937 2.926211e-02
# speciesPZI              -1.886266  0.7996874 -2.358755 1.833638e-02
# speciesSER               2.045764  0.6653866  3.074550 2.108203e-03
# speciesVUL               2.297478  0.6824845  3.366344 7.617163e-04
# log(nh)                  1.270419  0.2739997  4.636571 3.542359e-06
# methnecropsy:speciesLBI  4.389700  1.8288421  2.400262 1.638334e-02
# methnecropsy:speciesLBR  4.707842  1.8434226  2.553859 1.065364e-02
rm(s)

df2=as.data.frame(table(DfEur2$species,DfEur2$pab, DfEur2$meth))
colnames(df2)=c("species","abs" ,"meth", "prop")
df2=subset(df2, species %in% c("LBI", "LBR"))
ggplot(data=df2, aes(x=meth, y=prop, fill=abs)) +
  geom_bar(stat="identity", position="fill")+
  xlab("Presence or absence")+
  ylab("Proportion of presence and absence")+
  ggtitle("Proportion of presence and absence between sampling method")+
  facet_wrap(~species)
rm(DfEur, DfEur2, curWarnings, df2)

################################################################################################################################################################################
########################################################## Data for graphic with variable with low levels ######################################################################
################################################################################################################################################################################

# selection of necropsy only 
DfNecObs=subset(PAB2, env %in% c("Af-Temp-nec", "Am-Temp-nec", "Am-Trop-nec", "As-Cont-nec","Eu-Cont-nec","Eu-Temp-nec", "Oc-Arid-nec", "Oc-Temp-nec"))
DfNecObs$env=factor(DfNecObs$env)
DfNecObs$clim=factor(DfNecObs$clim)
DfNecObs$continent=factor(DfNecObs$continent)
DfNecObs$species=factor(DfNecObs$species)

# subset for clim diff 
DfNeclimObs=subset(DfNecObs, clim %in% c("Continental", "Temperate","Tropical"))
DfNeclimObs$clim=factor(DfNeclimObs$clim)
DfNeclimObs$continent=factor(DfNeclimObs$continent)
DfNeclimObs$env=factor(DfNeclimObs$env)
DfNeclimObs$species=factor(DfNeclimObs$species)

# subset for continent diff 
DfNecontObs=subset(DfNecObs, continent %in% c("America", "Europe"))
DfNecontObs$continent=factor(DfNecontObs$continent)
DfNecontObs$species=factor(DfNecontObs$species)
DfNecontObs$env=factor(DfNecontObs$env)
DfNecontObs$clim=factor(DfNecontObs$clim)


# contclim creation
DfNecObs$contclim2=factor(paste(DfNecObs$continent, DfNecObs$clim, sep=""))
DfNecontclim2Obs=subset(DfNecObs, contclim2 %in% c("AmericaTropical", "EuropeContinental", "EuropeTemperate", "AmericaTemperate"))
DfNecontclim2Obs$contclim2=factor(DfNecontclim2Obs$contclim2)
table(DfNecontclim2Obs$contclim2, DfNecontclim2Obs$species)

################################################################################################################################################################################
##########################################################   Data for analysis wiThout var with low lvl   ######################################################################
################################################################################################################################################################################

# selection of necropsy only 
DfNecAnalysis=subset(PAB2, env %in% c("Am-Temp-nec", "Am-Trop-nec","Eu-Cont-nec","Eu-Temp-nec"))
DfNecAnalysis$env=factor(DfNecAnalysis$env)
DfNecAnalysis$clim=factor(DfNecAnalysis$clim)
DfNecAnalysis$continent=factor(DfNecAnalysis$continent)
DfNecAnalysis$species=factor(DfNecAnalysis$species)

# subset for clim diff 
DfNeclimAnalysis=subset(DfNecAnalysis, clim %in% c("Continental", "Temperate","Tropical"))
DfNeclimAnalysis$clim=factor(DfNeclimAnalysis$clim)
DfNeclimAnalysis$continent=factor(DfNeclimAnalysis$continent)
DfNeclimAnalysis$env=factor(DfNeclimAnalysis$env)
DfNeclimAnalysis$species=factor(DfNeclimAnalysis$species)

# subset for continent diff 
DfNecontAnalysis=subset(DfNecAnalysis, continent %in% c("America", "Europe"))
DfNecontAnalysis$continent=factor(DfNecontAnalysis$continent)
DfNecontAnalysis$species=factor(DfNecontAnalysis$species)
DfNecontAnalysis$env=factor(DfNecontAnalysis$env)
DfNecontAnalysis$clim=factor(DfNecontAnalysis$clim)


# contclim creation
DfNecAnalysis$contclim2=factor(paste(DfNecAnalysis$continent, DfNecAnalysis$clim, sep=""))
DfNecontclim2Analysis=subset(DfNecAnalysis, contclim2 %in% c("AmericaTropical", "EuropeContinental", "EuropeTemperate", "AmericaTemperate"))
DfNecontclim2Analysis$contclim2=factor(DfNecontclim2Analysis$contclim2)
DfNecontclim2Analysis$species=factor(DfNecontclim2Analysis$species)
table(DfNecontclim2Analysis$contclim2, DfNecontclim2Analysis$species)

# global effect? 
s=data.frame(summary(glmer(pab ~ env*species+GrpYear+nh+ (1|study_ID), family=binomial,data=DfNecAnalysis))$coefficients)
s[s$Pr...z..<0.05,]
# envEu-Temp-nec            -3.083757   1.437553 -2.145144 0.031941368
# envEu-Temp-nec:speciesCBI  4.233539   1.972499  2.146282 0.031850511
# envEu-Temp-nec:speciesELO  5.541861   1.929836  2.871674 0.004083035
# envEu-Temp-nec:speciesLBR  3.438940   1.737315  1.979457 0.047764561
# envEu-Temp-nec:speciesULT  3.553181   1.709201  2.078856 0.037630622
DfNecontclim2Analysis$contclim2=relevel(DfNecontclim2Analysis$contclim2, ref="AmericaTropical")
# effet geoclimatique sur les especes?
s=data.frame(summary(glmer(pab ~ contclim2*species+GrpYear+nh+ (1|study_ID), family=binomial,data=DfNecontclim2Analysis))$coefficients)
s[s$Pr...z..<0.05,]
# Estimate Std..Error   z.value    Pr...z..
# contclim2EuropeTemperate            -3.226820   1.455805 -2.216519 0.026655968
# contclim2EuropeTemperate:speciesCBI  4.295002   1.984201  2.164600 0.030418300
# contclim2EuropeTemperate:speciesELO  5.613583   1.941682  2.891092 0.003839054
# contclim2EuropeTemperate:speciesLBR  3.493980   1.747493  1.999425 0.045562422
# contclim2EuropeTemperate:speciesULT  3.605502   1.719716  2.096568 0.036031814
rm(s)

df2=as.data.frame(table(DfNecontclim2Obs$species,DfNecontclim2Obs$pab, DfNecontclim2Obs$contclim2))
colnames(df2)=c("species","abs" ,"contclim2", "prop")
df2=subset(df2, species %in% c("CBI", "LBI", "ELO", "LBR", "RAD"))
ggplot(data=df2, aes(x=contclim2, y=prop, fill=abs)) +
  geom_bar(stat="identity", position="fill")+
  xlab("Presence or absence")+
  ylab("Proportion of presence and absence")+
  ggtitle("Proportion of presence and absence between sampling method")+
  facet_wrap(~species)
rm(df2, DfNecontclim2Obs,DfNecontclim2Analysis)


# effet continental?
s=data.frame(summary(glmer(pab ~ continent*species+GrpYear+nh+ (1|study_ID), family=binomial,data=DfNecontAnalysis))$coefficients)
s[s$Pr...z..<0.05,]
rm(s)
# Estimate Std..Error   z.value    Pr...z..
# speciesHYB                 -2.570392   1.283026 -2.003383 0.045136221
# continentEurope:speciesELO  4.605862   1.783468  2.582531 0.009807844
# continentEurope:speciesLBI  5.405788   1.924392  2.809088 0.004968201
# continentEurope:speciesLBR  3.644684   1.542667  2.362586 0.018147935


df2=as.data.frame(table(DfNecontObs$species,DfNecontObs$pab, DfNecontObs$continent,DfNecontObs$env ))
colnames(df2)=c("species","abs" ,"continent", "env","prop")
df2=subset(df2, species %in% c("CBI", "ELO","LBI", "LBR", "PMB", "POC", "PZI"))
ggplot(data=df2, aes(x=continent, y=prop, fill=abs)) +
  geom_bar(stat="identity", position="fill")+
  xlab("Presence or absence")+
  ylab("Proportion of presence and absence")+
  facet_wrap(~species)

rm(s, df2, DfNecontObs,DfNecontAnalysis)
# effet climatique?

s=data.frame(summary(glmer(pab ~ clim*species+GrpYear+nh+ (1|study_ID), family=binomial,data=DfNeclimAnalysis))$coefficients)
s[s$Pr...z..<0.05,]
# [1] Estimate   Std..Error z.value    Pr...z..  
# <0 lignes> (ou 'row.names' de longueur nulle)
# > 

rm(DfNec, curWarnings,DfNecontclim2,DfNeclim,DfNecont,PAB2, DfNecAnalysis,DfNecObs,s,DfNeclimAnalysis, DfNeclimObs)

# why some species are affected by the continent and not by env? 
df2=subset(DfNeclimAnalysis, species %in% c("CBI", "ELO","LBI", "LBR", "PMB", "POC", "PZI"))
df2$species=factor(df2$species)
s=data.frame(summary(glmer(pab ~ clim*species+GrpYear+nh+ (1|study_ID), family=binomial,data=df2))$coefficients)
s[s$Pr...z..<0.05,]


#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################
#####################################################################################################################################################################

rm(list=ls())
######---- let s see the construction of the dataset

ma = read.csv(file = './PaB4.csv', header=T, sep=';')
colnames(ma)=c("study_ID","author","year", "continent" , "nh"  , "meth","clim","kop2","species" ,"pab")
ma$env = factor(paste(substr(ma$continent,1,2),substr(ma$clim,1,4),substr(ma$meth,1,3),sep='-'))
dt = ma %>% dcast(study_ID + year + clim + continent + nh + kop2 + meth + env ~ species, value.var = 'pab')
dt[is.na(dt)] = 0
dt$meth = factor(dt$meth)
dt = data.frame(dt)


####-------======== Species test - ENV ======--------------

## 10% prevalence
sp5 = names(which(colSums(dt[,9:(dim(dt)[2]-1)])>=0.1*dim(dt)[1] & colSums(dt[,9:(dim(dt)[2]-1)])<0.9*dim(dt)[1]))

env = dt[,c('study_ID','year','nh','env')]
sp = dt[,c(which(colnames(dt) %in% sp5))]
sp = data.frame(lapply(sp,as.integer))
dt3 = cbind(env,sp)
dt3 = data.frame(dt3)

dt3$env = factor(dt3$env)

ma3 = melt(dt3,1:4) ## for LMER

ma3$study_ID =factor(match(ma3$study_ID,levels(factor(ma3$study_ID))))



###--- Consider method effect

Eu = dt[dt$continent=='Europe' & dt$clim=='Continental',]

dim(Eu)

# [1]  33 56
table(Eu$meth)

# deworming  necropsy

# 28         5



## Focus on species with 10% prevalence within this subset for model convergence

sp5 = names(which(colSums(Eu[,9:(dim(Eu)[2]-1)])>=0.1*dim(Eu)[1] & colSums(Eu[,9:(dim(Eu)[2]-1)])<0.9*dim(Eu)[1]))

sp5

# "ACU" "ASH" "BID" "BRE" "CAP" "CBI" "COR" "EDE" "ELO" "EQU" "HYB" "LBR" "NIP" "PAT" "PEU" "PMB" "POC" "PTM" "PZI" "RAD"
#  "SER" "TBV" "TEN" "ULT"

length(sp5)

#[1] 24

env = Eu[,c('study_ID','year','nh','env')]

sp = Eu[,c(which(colnames(Eu) %in% sp5))]

sp = data.frame(lapply(sp,as.integer))

Eu3 = cbind(env,sp)

Eu3 = data.frame(Eu3)

Eu3$env = factor(Eu3$env)

ma3 = melt(Eu3,1:4) ## for LMER

ma3$study_ID =factor(match(ma3$study_ID,levels(factor(ma3$study_ID))))

ma3$env = factor(gsub('Eu-Cont-','',ma3$env))



### Try model implementation on this dataset



mMeth1 = glmer(value ~ log(nh) + env*variable + (1|study_ID),
               
               family = binomial(link = 'logit'),
               
               control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
               
               data = ma3)

# Warning messages:

#   1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :

#                     unable to evaluate scaled gradient

#                   2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :

#                                     Model failed to converge: degenerate  Hessian with 13 negative eigenvalues

summary(mMeth1)
# Fixed effects:
#    Estimate Std. Error z value Pr(>|z|)    
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)        -4.269e+00  8.987e-01  -4.750 2.04e-06 ***
#    log(nh)             1.307e+00  2.811e-01   4.651 3.30e-06 ***
#    envnec             -2.319e+00  1.431e+00  -1.621  0.10509    
# variableASH         5.075e+00  1.231e+00   4.123 3.75e-05 ***
#    variableBID        -2.100e-01  6.399e-01  -0.328  0.74282    
# variableBRE        -1.521e+00  7.449e-01  -2.042  0.04119 *  
#    variableCAP         8.019e-06  6.330e-01   0.000  0.99999    
# variableCBI         7.616e-06  6.330e-01   0.000  0.99999    
# variableCOR         2.056e+00  6.671e-01   3.082  0.00206 ** 
#    variableEDE         4.012e-01  6.256e-01   0.641  0.52130    
# variableELO         1.602e+00  6.441e-01   2.487  0.01289 *  
#    variableEQU        -1.199e+00  7.070e-01  -1.695  0.09002 .  
# variableHYB         7.861e-06  6.330e-01   0.000  0.99999    
# variableLBR         2.888e+00  7.359e-01   3.924 8.71e-05 ***
#    variableNIP        -6.640e-01  6.629e-01  -1.002  0.31649    
# variablePAT         3.667e+00  8.423e-01   4.354 1.34e-05 ***
#    variablePEU        -9.175e-01  6.812e-01  -1.347  0.17804    
# variablePMB         1.391e+00  6.365e-01   2.185  0.02888 *  
#    variablePOC         1.602e+00  6.441e-01   2.487  0.01289 *  
#    variablePTM         1.391e+00  6.365e-01   2.185  0.02888 *  
#    variablePZI        -1.909e+00  8.049e-01  -2.371  0.01773 *  
#    variableRAD         7.920e-01  6.251e-01   1.267  0.20517    
# variableSER         2.056e+00  6.671e-01   3.082  0.00206 ** 
#    variableTBV         5.970e-01  6.246e-01   0.956  0.33917    
# variableTEN        -6.640e-01  6.629e-01  -1.002  0.31649    
# variableULT        -6.640e-01  6.629e-01  -1.002  0.31649    
# envnec:variableASH -5.075e+00  2.132e+00  -2.381  0.01728 *  
#    envnec:variableBID -1.614e+01  2.924e+03  -0.006  0.99560    
# envnec:variableBRE  1.521e+00  1.893e+00   0.803  0.42179    
# envnec:variableCAP  1.931e+01  2.621e+03   0.007  0.99412    
# envnec:variableCBI  1.933e+01  2.649e+03   0.007  0.99418    
# envnec:variableCOR  1.727e+01  2.639e+03   0.007  0.99478    
# envnec:variableEDE -4.013e-01  1.849e+00  -0.217  0.82821    
# envnec:variableELO  1.779e+01  2.721e+03   0.007  0.99478    
# envnec:variableEQU  1.199e+00  1.878e+00   0.638  0.52345    
# envnec:variableHYB  2.226e+00  1.705e+00   1.306  0.19172    
# envnec:variableLBR  1.651e+01  2.740e+03   0.006  0.99519    
# envnec:variableNIP  6.639e-01  1.862e+00   0.357  0.72145    
# envnec:variablePAT -2.803e-01  1.929e+00  -0.145  0.88448    
# envnec:variablePEU  3.144e+00  1.724e+00   1.824  0.06819 .  
# envnec:variablePMB  1.797e+01  2.686e+03   0.007  0.99466    
# envnec:variablePOC  1.778e+01  2.709e+03   0.007  0.99476    
# envnec:variablePTM  1.799e+01  2.709e+03   0.007  0.99470    
# envnec:variablePZI  2.115e+01  2.527e+03   0.008  0.99332    
# envnec:variableRAD  1.854e+01  2.649e+03   0.007  0.99442    
# envnec:variableSER -2.056e+00  1.864e+00  -1.103  0.27000    
# envnec:variableTBV -5.970e-01  1.849e+00  -0.323  0.74678    
# envnec:variableTEN -1.562e+01  2.831e+03  -0.006  0.99560    
# envnec:variableULT  2.890e+00  1.717e+00   1.684  0.09224 ..  
# ---
anova(mMeth1)
# Analysis of Variance Table
# npar  Sum Sq Mean Sq F value
# log(nh)         1   3.811  3.8106  3.8106
# env             1   2.636  2.6362  2.6362
# variable       23 128.525  5.5881  5.5881
# env:variable   23  17.067  0.7421  0.7421

### Check for complete separation with env

tabsp = data.frame(table(ma3$variable,ma3$value,ma3$env))
colnames(tabsp)=c("species","Presence", "Recovery_method", "Frequency")
tabsp$Recovery_method=ifelse(tabsp$Recovery_method=="dew", "Deworming", "Necropsy")
ggplot(tabsp,aes(x = species, y = Frequency, group = paste0(Recovery_method,species), fill = Presence))+
  
  geom_bar(stat='identity') + theme_bw() +
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text = element_text(face="bold", size=14),
        strip.background = element_rect(fill="white", colour="white",size=5))+
  
  
  facet_wrap(~ Recovery_method,ncol=1)

tabsp$Recovery_method=ifelse(tabsp$Recovery_method=="Deworming", "dew", "nec")
colnames(tabsp)=c("Var1","Var2", "Var3","Frequency")


## Complete separation observed for a few species

filt = tabsp[tabsp$Var2==1,]



nstu = data.frame(table(Eu3$env))

colnames(nstu)=c('Var3','tot')

nstu$Var3 = gsub('Eu-Cont-','',nstu$Var3)



filt = merge(nstu,filt,by = 'Var3')

filt$CompSep = filt$Freq/filt$tot

## Remove species with at least 1 separation

elim_sp = names(table(filt$Var1[filt$CompSep %in% c(0,1)])[table(filt$Var1[filt$CompSep %in% c(0,1)])>=1])

elim_sp

#[1] "BID" "CAP" "CBI" "COR" "ELO" "LBR" "PMB" "POC" "PTM" "PZI" "RAD" "TEN"



ma4 = ma3[!(ma3$variable %in% elim_sp),]

ma4$variable = factor(ma4$variable)

table(ma4$variable)

# ACU ASH BRE EDE EQU HYB NIP PAT PEU SER TBV ULT 
# 33  33  33  33  33  33  33  33  33  33  33  33



tabstu = data.frame(table(ma4$variable,ma4$value,ma4$study_ID))

a=aggregate(Freq ~ Var2+ Var3, FUN=sum, data=tabstu) ## check complete sep by stud

a[a$Freq<3 |a$Freq>sum(a$Freq[c(1,2)])-3,]

# Var2 Var3 Freq
# 5     0    3   10
# 6     1    3    2
# 13    0    7   11
# 14    1    7    1
# 23    0   12   10
# 24    1   12    2
# 35    0   18    1
# 36    1   18   11
# 49    0   25   10
# 50    1   25    2
# 59    0   30   11
# 60    1   30    1
# 61    0   31   10
# 62    1   31    2
# 65    0   33   10
# 66    1   33    2



stu_elim = unique(a[a$Freq<3 |a$Freq>sum(a$Freq[c(1,2)])-3,])

#

ma4 = ma4[!(ma4$study_ID %in% stu_elim),]

ma4$study_ID = factor(ma4$study_ID)

table(ma4$study_ID)

### Mixed model for species

m = glmer(value ~ log(nh) + year + env*variable + (1|study_ID),
          
          family = binomial(link = 'logit'),
          
          control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
          
          data = ma4)

car::Anova(m)

# Response: value

# Chisq Df Pr(>Chisq)    
# log(nh)      15.5448  1  8.057e-05 ***
#    year         14.0948  1  0.0001738 ***
#    env           2.2186  1  0.1363594    
# variable     46.0953 11  3.110e-06 ***
#    env:variable 25.0262 11  0.0090364 ** 
#    ---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 
# 


s=data.frame(summary(m)$coefficients)

s[s$Pr...z..<0.05,]
# 
# Intercept)        -23.077466028 5.116324071 -4.510556 6.465794e-06
# log(nh)              1.237474911 0.313865675  3.942690 8.057293e-05
# year                 0.009464842 0.002521064  3.754305 1.738232e-04
# variableASH          4.646695998 1.130958743  4.108634 3.980056e-05
# variableBRE         -1.465297939 0.739971745 -1.980208 4.768017e-02
# variablePAT          3.359366242 0.819635311  4.098611 4.156370e-05
# variableSER          1.917766859 0.668636795  2.868174 4.128480e-03
# envnec:variableASH  -4.646728152 2.051299376 -2.265261 2.349668e-02


# Removal of year effect does not change the conclusions

m2 = glmer(value ~ log(nh) + env*variable + (1|study_ID),
           
           family = binomial(link = 'logit'),
           control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
           
           data = ma4)

car::Anova(m2)

# Response: value
#Chisq Df Pr(>Chisq)    
# log(nh)      16.8347  1  4.078e-05 ***
#    env           4.3457  1   0.037102 *  
#    variable     46.0949 11  3.111e-06 ***
#    env:variable 26.6394 11   0.005208 **


s2=data.frame(summary(m2)$coefficients)

s2[s2$Pr...z..<0.05,]
# Intercept)        -4.314496  1.0146732 -4.252104 2.117719e-05
# log(nh)             1.334328  0.3252073  4.103006 4.078163e-05
# variableASH         4.965403  1.2218822  4.063733 4.829406e-05
# variableBRE        -1.550929  0.7548505 -2.054617 3.991601e-02
# variablePAT         3.606966  0.8584971  4.201489 2.651655e-05
# variableSER         2.057294  0.6822200  3.015588 2.564819e-03
# envnec:variableASH -4.965394  2.1143867 -2.348385 1.885502e-02
rm(a,dt3,env,Eu, Eu3, filt, m, m2, ma3, ma4, mMeth1, nstu,sp,stu_elim, tabsp, tabstu, elim_sp, s, dt, s2)

########################################################################################
## test on contclim 
ma$contclim=paste(ma$continent, ma$clim, sep="")
dt = ma %>% dcast(study_ID + year + clim + continent + nh + kop2 + meth + env + contclim ~ species, value.var = 'pab')
dt[is.na(dt)] = 0
dt$meth = factor(dt$meth)
dt = data.frame(dt)

###--- Consider method effect
Cclim=subset(dt, contclim %in% c("AmericaTropical", "EuropeContinental", "EuropeTemperate", "AmericaTemperate"))
dim(Cclim)

# [1] 61 58

table(Cclim$contclim)
# AmericaTemperate   AmericaTropical EuropeContinental   EuropeTemperate 
# 8                 7                33                14 


## Focus on species with 10% prevalence within this subset for model convergence
sp5 = names(which(colSums(Cclim[,10:(dim(Cclim)[2]-1)])>=0.1*dim(Cclim)[1] & colSums(Cclim[,10:(dim(Cclim)[2]-1)])<0.9*dim(Cclim)[1]))
sp5
# "ACU" "ASH" "ASY" "BID" "BRE" "CAP" "CBI" "COR" "EDE" "ELO" "EQU" "HYB" "LBR" "NIP" "PEU" "PMB" "POC" "PTM" "PZI" "RAD"
# [21] "SER" "TBV" "TEN" "TET" "ULT"
# > 

length(sp5)

#[1] 25

env = Cclim[,c('study_ID','year','nh','contclim')]

sp = Cclim[,c(which(colnames(Cclim) %in% sp5))]

sp = data.frame(lapply(sp,as.integer))

Cclim2 = cbind(env,sp)

Cclim2 = data.frame(Cclim2)

Cclim2$contclim = factor(Cclim2$contclim)

ma3 = melt(Cclim2,1:4) ## for LMER

ma3$study_ID =factor(match(ma3$study_ID,levels(factor(ma3$study_ID))))

rm(env, sp)
### Try model implementation on this dataset

mMeth1 = glmer(value ~ log(nh) + contclim*variable + (1|study_ID),
               
               family = binomial(link = 'logit'),
               
               control = glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),
               
               data = ma3)

# Warning messages:

#   1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :

#                     unable to evaluate scaled gradient

#                   2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :

#                                     Model failed to converge: degenerate  Hessian with 13 negative eigenvalues

summary(mMeth1)
s=data.frame(summary(mMeth1)$coefficients)

s[s$Pr...z..<0.05,]
# log(nh)                                0.7487403  0.1672455  4.476894 0.0000075737
# contclimEuropeTemperate               -3.7269078  1.4194443 -2.625611 0.0086493728
# variableBRE                           -3.0829570  1.4324435 -2.152236 0.0313787438
# contclimEuropeContinental:variableASH  4.4043739  1.3392939  3.288579 0.0010069452
# contclimAmericaTropical:variableBRE    4.4231477  1.8534048  2.386498 0.0170096788
# contclimEuropeTemperate:variableBRE    4.4496826  1.9262793  2.309988 0.0208888012
# contclimEuropeContinental:variableELO  2.6109025  1.2755612  2.046866 0.0406712751
# contclimEuropeTemperate:variableELO    4.6037688  1.6963183  2.713977 0.0066480797
# contclimEuropeTemperate:variableULT    3.5217317  1.6714357  2.107010 0.0351167206


tabsp = data.frame(table(ma3$variable,ma3$value,ma3$contclim))
colnames(tabsp)=c("species","Presence", "Contclim", "Frequency")
ggplot(tabsp,aes(x = species, y = Frequency, group = paste0(Contclim,species), fill = Presence))+
  
  geom_bar(stat='identity') + theme_bw() +
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14), 
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        strip.text = element_text(face="bold", size=14),
        strip.background = element_rect(fill="white", colour="white",size=5))+
  facet_wrap(~ Contclim,ncol=1)

#####################################################################################################################################################################
###############                           Diversity Analysis                                      ###################################################################
#####################################################################################################################################################################
PAB=read.csv("PaB4.csv", header=TRUE,sep=";") ## Big data set with pres/abs from prev and abd
PAB$GrpYear=ifelse(PAB$yearsPubli > 1960, "After", "Before")
PAB$env = factor(paste(substr(PAB$continent,1,2),substr(PAB$clim,1,4),substr(PAB$meth,1,3),sep='-')) ## Factor to be tested too
colnames(PAB)=c( "study_ID","author","yearsPubli", "continent",  "nh",   "method","climate","kop2","species",
                 "pab","GrpYear","env" )
dfRichness= dcast(PAB, study_ID + author + yearsPubli  + continent + nh + method + climate + kop2 +      
                    GrpYear + env ~ species, value.var="pab")
dfRichness$richness=rowSums(dfRichness[,-c(1:10)])
dfRichness$contclim=paste(dfRichness$continent,dfRichness$clim, sep="")

table(dfRichness$continent)
# Africa America    Asia  Europe Oceania 
# 3      15       1      47       3 


######################################################################################
############################## analysis on richness index ############################
######################################################################################

dfRichnessAnalysis=subset(dfRichness, env %in% c("Am-Temp-nec", "Am-Trop-nec" ,"Eu-Cont-dew" ,
                                                 "Eu-Cont-nec", "Eu-Temp-dew", "Eu-Temp-nec" ))
dfRichnessAnalysis$climate=factor(dfRichnessAnalysis$climate)
ggplot(dfRichnessAnalysis)+
  geom_boxplot(aes(x=env, y=richness, fill=env))+
  xlab("")+
  ylab("Richness")+
  scale_fill_brewer(palette="Accent")+
  theme_bw()
summary(glm(richness ~ env+nh+GrpYear, 
            data=dfRichnessAnalysis, family="gaussian"))


dfRichnessEur=subset(dfRichness, env %in% c("Eu-Cont-dew" ,"Eu-Cont-nec", "Eu-Temp-dew", "Eu-Temp-nec" ))
dfRichnessEur$climate=factor(dfRichnessEur$climate)
table(dfRichnessEur$clim, dfRichnessEur$meth)

# deworming necropsy
#   Continental        28        5
#   Temperate           3       11
dfRichnessEur=subset(dfRichnessEur, climate %in% c("Continental"))
dfRichnessEur$climate=factor(dfRichnessEur$climate)
dfRichnessEur$continent=factor(dfRichnessEur$continent)
dfRichnessEur$contclim=factor(dfRichnessEur$contclim)
dfRichnessEur$env=factor(dfRichnessEur$env)
shapiro.test(dfRichnessEur$richness) # W = 0.95475, p-value = 0.1828, normality ok.

# effect of the method on the richness index? 
ggplot(dfRichnessEur)+
  geom_boxplot(aes(x=method, y=richness, fill=method))+
  labs(fill="Recovery method")+
  xlab("")+
  ylab("Richness index")+
  ggtitle("Repartition of richness index between recovery method")+
  scale_fill_brewer(palette="Set1")+
  theme_bw()
summary(glm(richness ~ method+nh+GrpYear, 
            data=dfRichnessEur, family="gaussian")) ### no effect of the method
# Coefficients:
#    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)   16.82180    1.32601  12.686 2.32e-13 ***
#    methnecropsy  -1.32818    3.27019  -0.406  0.68762    
# nh             0.18040    0.05307   3.399  0.00198 ** 
#    GrpYearbefore  2.60210    4.33434   0.600  0.55294      

rm(dfRichnessEur)

###################"" effect of climate (diff between continental and temperate)
dfRichnessClimO=subset(dfRichness, climate %in% c("Continental", "Temperate", "Tropical"))
dfRichnessClim=subset(dfRichnessClimO, continent %in% c("Europe", "America"))
dfRichnessClim$climate=factor(dfRichnessClim$climate)
dfRichnessClim$continent=factor(dfRichnessClim$continent)
dfRichnessClim$contclim=factor(dfRichnessClim$contclim)
dfRichnessClim$env=factor(dfRichnessClim$env)

# diff between climate?
ggplot(dfRichnessClimO,aes(x=climate, y=richness, fill=climate))+
  geom_dotplot(binaxis='y',stackdir='center', binwidth = 0.45)+
  labs(color = "Continent-climat")+
  theme_bw()+
  scale_fill_viridis(discrete=TRUE,option="magma")+
  theme(axis.text.x = element_blank())+
  xlab("Climat")+
  ylab('Richness index')+
  ggtitle("Repartition of richness index between climate")

summary(glm(richness ~ climate*method+GrpYear+nh, data=dfRichnessClim, family="gaussian"))

# diff between contclim?
ggplot(dfRichnessClim,aes(x=contclim, y=richness, fill=contclim))+
  geom_dotplot(binaxis='y',stackdir='center', binwidth = 0.45)+
  labs(fill = "Continent-climat")+
  theme_bw()+
  scale_fill_viridis(discrete=TRUE,option="magma")+
  theme(axis.text.x = element_blank())+
  xlab("")+
  ylab('Richness index')+
  ggtitle("Repartition of richness index between climate")

dfRichnessClim$contclim=relevel(dfRichnessClim$contclim, ref="AmericaTropical")
summary(glm(richness ~ contclim+GrpYear+nh, data=dfRichnessClim, family="gaussian"))
# Call:
#    glm(formula = richness ~ contclim + GrpYear + nh, family = "gaussian", 
#        data = dfRichnessClim)
# 
# Deviance Residuals: 
#    Min        1Q    Median        3Q       Max  
# -10.5728   -2.7366    0.4272    3.2067   10.1946  
# 
# Coefficients:
#    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               21.18461    2.04005  10.384 1.42e-14 ***
#    contclimAmericaTropical    0.36232    2.66819   0.136  0.89248    
# contclimEuropeContinental -1.89626    1.99292  -0.951  0.34551    
# contclimEuropeTemperate   -6.15918    2.27531  -2.707  0.00903 ** 
#    GrpYearBefore              0.20788    2.32316   0.089  0.92902    
# nh                         0.09482    0.03578   2.650  0.01050 *  
#    ---
rm(dfRichnessClim, dfRichnessClimO)

###################"" effect of continent (within temperate data)
dfRichnessContO=subset(dfRichness, climate %in% c("Temperate"))
dfRichnessCont=subset(dfRichnessContO, method %in% c("necropsy"))
dfRichnessCont=subset(dfRichnessCont, continent %in% c("Europe", "America"))
dfRichnessCont$clim=factor(dfRichnessCont$climate)
dfRichnessCont$continent=factor(dfRichnessCont$continent)
dfRichnessCont$contclim=factor(dfRichnessCont$contclim)
dfRichnessCont$env=factor(dfRichnessCont$env)

ggplot(dfRichnessContO,aes(x=continent, y=richness, fill=continent))+
  geom_dotplot(binaxis='y',stackdir='center', binwidth = 0.45)+
  labs(color = "Continent-climat")+
  theme_bw()+
  scale_fill_viridis(discrete=TRUE,option="magma")+
  theme(axis.text.x = element_blank())+
  xlab("Continent")+
  ylab('Richness index')+
  ggtitle("Repartition of richness index between Continent")+
  facet_wrap(~climate)


# dfRichnessCont$continent=relevel(dfRichnessCont$continent, ref="Europe")
summary(glm(richness ~ continent+GrpYear+nh, data=dfRichnessCont,family="gaussian"))
# Coefficients:
#    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     24.50694    2.44132  10.038 8.91e-08 ***
#    continentEurope -3.66092    2.36434  -1.548    0.144    
# GrpYearBefore   -2.34383    3.80122  -0.617    0.547    
# nh              -0.02227    0.06301  -0.353    0.729  
rm(dfRichness,dfRichnessCont, dfRichnessContO,dfRichnessAnalysis)
#==========================================================================================================================================================#
#------------------------------------------------------------------   	 beta diversity      	---------------------------------------------------------#
#==========================================================================================================================================================#

df2=dcast(PAB, study_ID+author+yearsPubli+continent+nh+method+climate+kop2+GrpYear+env ~ species, value.var="pab")
df2[is.na(df2)] <- 0
df2=subset(df2, env %in% c("Am-Temp-nec", "Am-Trop-nec" ,"Eu-Cont-dew" ,
                           "Eu-Cont-nec", "Eu-Temp-dew", "Eu-Temp-nec" ))
df2$env=factor(df2$env)
df2$climate=factor(df2$climate)
df2$continent=factor(df2$continent)

# global effect? 

dfEnv <- df2[,-c(11:58)] #### we need to keep only meta data 
head(dfEnv)
# study_ID  author yearsPubli continent nh    method     climate   kop2 GrpYear         env
# 2        1   Silva       1999   America 36  necropsy    Tropical Others   After Am-Trop-nec
# 3        2 kuzmina       2008    Europe 12 deworming Continental    Dfb   After Eu-Cont-dew
# 4        3  foster       1936   America 17  necropsy    Tropical Others  Before Am-Trop-nec
# 5       57  guiris       2010   America  2  necropsy    Tropical Others   After Am-Trop-nec
# 6        4  Kornas       2011    Europe 14 deworming Continental    Dfb   After Eu-Cont-dew
# 7       65   Salle       2018    Europe 14 deworming Continental    Dfb   After Eu-Cont-dew

dfEnv$author=NULL
dfEnv$climate=NULL
dfEnv$kop2=NULL
dfEnv$continent=NULL
dfEnv$yearsPubli=NULL
dfEnv$GrpYear=NULL
dfEnv$method=NULL

dfsp <- df2[,c(11:58)]
head(dfsp)
# ACU ADE ALV ASH ASY AUR BID BRE CAL CAP CAT CBI CCR COR EDE ELO EQU GOL HYB INS LBI LBR LEP LON MIN MNR MON MUC NAS NIP ORB
# 2   0   0   1   1   0   0   0   1   1   1   1   1   0   1   0   0   0   1   0   1   1   1   1   1   1   0   0   0   1   0   0
# 3   0   0   0   1   0   0   0   1   1   1   1   0   0   1   1   1   1   1   0   1   1   1   1   1   1   0   0   0   1   1   0
# 4   0   0   0   0   0   0   0   1   1   1   1   0   0   1   1   1   1   1   1   1   1   1   1   1   1   1   0   1   1   0   0
# 5   0   0   0   0   1   0   1   0   0   0   0   1   1   0   1   0   1   0   0   1   1   1   1   1   1   0   0   0   0   0   0
# 6   0   0   0   1   0   0   0   1   1   0   1   0   0   1   0   1   0   1   0   1   1   1   1   1   1   0   0   0   1   0   0
# 7   0   0   0   1   0   0   0   0   1   0   1   1   0   1   1   1   0   1   0   1   1   1   1   1   1   0   0   0   1   0   0
# PAT PEU PKR PMB POC PTM PZI RAD SAG SER SKR TBV TEN TET TGB TRO ULT VUL
# 2   1   1   0   1   1   0   1   1   0   0   0   0   0   0   0   0   1   0
# 3   1   0   0   0   1   0   0   1   0   1   0   1   0   0   0   0   0   1
# 4   1   1   0   1   1   0   1   1   0   1   0   1   1   0   0   0   1   1
# 5   0   0   0   1   1   0   0   0   0   1   0   0   1   1   1   0   0   1
# 6   1   0   0   1   1   1   0   1   0   0   0   0   0   0   0   0   0   0
# 7   1   0   0   1   1   1   0   1   0   1   0   1   0   0   0   0   0   1
b=vegdist(dfsp, method="jacc", na.rm=TRUE)
adonis(b ~  ., data=dfEnv[,-1], perm=9999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# nh         1    0.5568 0.55681  6.9939 0.08702  1e-04 ***
#    env        5    1.5424 0.30848  3.8747 0.24106  1e-04 ***
#    Residuals 54    4.2991 0.07961         0.67191           
# Total     60    6.3983                 1.00000           
c=metaMDS(b)
mds_data <- as.data.frame(c$points)
mds_data$study_ID <- dfEnv$study_ID
mds_data <- merge(mds_data, dfEnv, by="study_ID")

ggplot(mds_data, aes(x = MDS1, y = MDS2, color=env)) +
  #scale_color_manual(values=c("darkred", "royalblue2"))+
  #scale_color_brewer(palette="Dark2")+
  geom_point(size = 5 , alpha = 0.6)+
  theme_bw()
rm(b,c,dfEnv,mds_data,dfsp)

# european subset for method effect

df=subset(df2, env %in% c("Eu-Cont-dew" ,
                          "Eu-Cont-nec", "Eu-Temp-dew", "Eu-Temp-nec" ))
df=subset(df, climate %in% c("Continental" ))
df$env=factor(df$env)
df$climate=factor(df$climate)
df$continent=factor(df$continent)

dfEnv <- df[,-c(11:58)] #### we need to keep only meta data 
dfEnv$author=NULL
dfEnv$climate=NULL
dfEnv$kop2=NULL
dfEnv$continent=NULL
dfEnv$yearsPubli=NULL
dfEnv$GrpYear=NULL
dfEnv$env=NULL
dfsp <- df[,c(11:58)]
dis = vegdist(dfsp,method="jacc")
mod = betadisper(dis, dfEnv$method)
par(mfrow=c(1,2))
plot(mod)
anova(mod)
b=vegdist(dfsp, method="jacc", na.rm=TRUE)
adonis(b ~  ., data=dfEnv[,-1], perm=9999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# nh         1   0.48829 0.48829  6.8438 0.16527 0.0003 ***
#    method     1   0.32579 0.32579  4.5663 0.11027 0.0012 ** 
#    Residuals 30   2.14042 0.07135         0.72446           
# Total     32   2.95450                 1.00000           
c=metaMDS(b)
mds_data <- as.data.frame(c$points)
mds_data$study_ID <- dfEnv$study_ID
mds_data <- merge(mds_data, dfEnv, by="study_ID")

ggplot(mds_data, aes(x = MDS1, y = MDS2, color=method)) +
  #scale_color_manual(values=c("darkred", "royalblue2"))+
  scale_color_brewer(palette="Dark2")+
  geom_point(size = 5 , alpha = 0.6)+
  theme_bw()+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.position = "none",
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14))

rm(c,df,dfEnv, dfsp, mds_data,b, mod)


#################################################### Necropsy only 
# effet du continent?
df=subset(df2, env %in% c("Am-Temp-nec", "Am-Trop-nec" ,
                          "Eu-Cont-nec", "Eu-Temp-nec" ))


df=subset(df, climate %in% c("Temperate" ))
df$env=factor(df$env)
df$climate=factor(df$climate)
df$continent=factor(df$continent)


dfEnv <- df[,-c(11:58)] #### we need to keep only meta data 
dfEnv$author=NULL
dfEnv$climate=NULL
dfEnv$kop2=NULL
dfEnv$yearsPubli=NULL
dfEnv$GrpYear=NULL
#dfEnv$env=NULL
dfEnv$method=NULL

dfsp <- df[,c(11:58)]
dis = vegdist(dfsp,method="jacc")
mod = betadisper(dis, dfEnv$env)
plot(mod)
anova(mod)
b=vegdist(dfsp, method="jacc", na.rm=TRUE)

adonis(b ~  ., data=dfEnv[,-1], perm=9999)
# continent  1   0.25084 0.250840  3.3061 0.15238 0.0094 **
#    nh         1   0.25728 0.257278  3.3909 0.15629 0.0076 **
#    Residuals 15   1.13808 0.075872         0.69134          
# Total     17   1.64620                  1.00000         

c=metaMDS(b)
mds_data <- as.data.frame(c$points)
mds_data$study_ID <- dfEnv$study_ID
mds_data <- merge(mds_data, dfEnv, by="study_ID")

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = continent)) +
  geom_point(size = 5 , alpha = 0.6)+
  scale_color_brewer(palette="Accent")+
  # scale_color_manual(values=c("red", "orange", "blue"))+
  theme_bw()+ labs(color = "Continent")

rm(a,c, mds_data,b,dfsp,dfEnv, df)


# effet du climat?
df=subset(df2, env %in% c("Am-Temp-nec", "Am-Trop-nec" ,
                          "Eu-Cont-nec", "Eu-Temp-nec" ))

df$env=factor(df$env)
df$climate=factor(df$climate)
df$continent=factor(df$continent)


dfEnv <- df[,-c(11:58)] #### we need to keep only meta data 
dfEnv$author=NULL
#dfEnv$continent=NULL
dfEnv$kop2=NULL
dfEnv$yearsPubli=NULL
dfEnv$GrpYear=NULL
dfEnv$env=NULL
dfEnv$method=NULL

dfsp <- df[,c(11:58)]

b=vegdist(dfsp, method="jacc", na.rm=TRUE)

adonis(b ~  ., data=dfEnv[,-1], perm=9999)
# continent  1   0.25084 0.250840  3.3061 0.15238 0.0094 **
#    nh         1   0.25728 0.257278  3.3909 0.15629 0.0076 **
#    Residuals 15   1.13808 0.075872         0.69134          
# Total     17   1.64620                  1.00000         

c=metaMDS(b)
mds_data <- as.data.frame(c$points)
mds_data$study_ID <- dfEnv$study_ID
mds_data <- merge(mds_data, dfEnv, by="study_ID")

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = climate, shape=continent)) +
  geom_point(size = 5 , alpha = 0.6)+
  scale_color_brewer(palette="Accent")+
  # scale_color_manual(values=c("red", "orange", "blue"))+
  theme_bw()+ labs(color = "Continent")

rm(c,dfEnv,dfsp, mds_data, b, df)

# Effet contclim?

df=subset(df2, env %in% c("Am-Temp-nec", "Am-Trop-nec" ,
                          "Eu-Cont-nec", "Eu-Temp-nec" ))

df$env=factor(df$env)
df$climate=factor(df$climate)
df$continent=factor(df$continent)
df$contclim=paste(df$continent,df$climate,sep="")

dfEnv <- df[,-c(11:58)] #### we need to keep only meta data 
dfEnv$author=NULL
dfEnv$continent=NULL
dfEnv$kop2=NULL
dfEnv$yearsPubli=NULL
dfEnv$GrpYear=NULL
dfEnv$env=NULL
dfEnv$method=NULL
dfEnv$climate=NULL
dfsp <- df[,c(11:58)]

dis = vegdist(dfsp,method="jacc")
mod = betadisper(dis, dfEnv$contclim)
plot(mod)
anova(mod)

b=vegdist(dfsp, method="jacc", na.rm=TRUE)
adonis(b ~  ., data=dfEnv[,-1], perm=9999)
# continent  1   0.25084 0.250840  3.3061 0.15238 0.0094 **
#    nh         1   0.25728 0.257278  3.3909 0.15629 0.0076 **
#    Residuals 15   1.13808 0.075872         0.69134          
# Total     17   1.64620                  1.00000         

c=metaMDS(b)
mds_data <- as.data.frame(c$points)
mds_data$study_ID <- dfEnv$study_ID
mds_data <- merge(mds_data, dfEnv, by="study_ID")

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = contclim)) +
  geom_point(size = 5 , alpha = 0.6)+
  scale_color_brewer(palette="Accent")+
  # scale_color_manual(values=c("red", "orange", "blue"))+
  theme_bw()+ labs(color = "Continent")


rm(list=ls())
###################### Relation between FEC and Gini simpson/species richness from individual data ####################
gs = read.csv(file = './dataGuillaume2018.csv',sep=';',header=T)
gsp = gs[,9:dim(gs)[2]]
fec = gs$Pre.treatment.Faecal.Egg.Count

## nsp
nsp = rowSums(ifelse(gsp>0,1,0))
require(vegetarian)
## Compute Gini-Simpson
gi = NULL
for(i in 1:dim(gsp)[1]){
  gi[i] = H(gsp[i,], lev = "alpha", wts = FALSE, q = 2, HCDT = FALSE, gini = TRUE, 
            boot = FALSE, boot.arg = list(s.sizes = NULL, num.iter = 100))
}


### Relationship with FEC
df = data.frame(N = nsp, FEC=fec, G = gi)
df$Fcat = 1
df$Fcat[df$FEC>=250 & df$FEC<500] = 2
df$Fcat[df$FEC>=500 & df$FEC<1000] = 3
df$Fcat[df$FEC>=1000] = 4

df$Fcat = factor(df$Fcat)
table(df$Fcat)
#  1  2  3  4 
# 14 14 10 10 

## Plot N et G sur le meême graphe (meme echelle)
ggplot(df,aes(x = N/max(N), y = G, col = Fcat)) + 
  geom_point(size = 4) + geom_smooth(method = 'lm')

ggplot(df,aes(x = Fcat, y = G, col = Fcat)) + 
  geom_boxplot()

ggplot(df,aes(x = Fcat, y = N, col = Fcat)) + 
  geom_boxplot()

summary(glm(N ~ Fcat, family = poisson(link='log'),data= df)) ## no diff with Fcat4


A=ggplot(df, aes(x=N, y=log(FEC+50)))+
  geom_point(color="royalblue2")+
  geom_smooth(method=lm,color="royalblue2")+
  xlab('Log(FEC)')+
  ylab("Species richness")+
  ggtitle("")+
  theme_bw()


B=ggplot(df, aes(x=G, y=log(FEC+50)))+
  geom_point(color="royalblue2")+
  geom_smooth(method=lm,color="royalblue2")+
  xlab('Log(FEC)')+
  ylab("Gini-Simpson index")+
  ggtitle("")+
  theme_bw()

rcorr(as.matrix(df),type = 'spearman')
#         N  FEC    G Fcat
# N    1.00 0.37 0.67 0.31
# FEC  0.37 1.00 0.24 0.97
# G    0.67 0.24 1.00 0.15
# Fcat 0.31 0.97 0.15 1.00
# 
# n= 48 
# 
# 
# P
#        N      FEC    G      Fcat  
# N           0.0090 0.0000 0.0324
# FEC  0.0090        0.1063 0.0000
# G    0 .0000 0.1063        0.3042
# Fcat 0.0324 0.0000 0.3042 


## Kuzmina data
kuz = read.csv(file='./data3_kuzmina2016.csv',sep=';',header=T)

kuzsp = kuz[,6:dim(kuz)[2]]
kuzfec = kuz$EPG
## nsp
nsp = rowSums(ifelse(kuzsp>0,1,0))

## Compute Gini-Simpson
giK = NULL
for(i in 1:dim(kuzsp)[1]){
  giK[i] = H(kuzsp[i,], lev = "alpha", wts = FALSE, q = 2, HCDT = FALSE, gini = TRUE, 
             boot = FALSE, boot.arg = list(s.sizes = NULL, num.iter = 100))
}


df = reshape2::melt(data.frame(N = nsp, FEC=kuzfec, G = giK),2)

ggplot(df,aes(x = FEC, y = value,col=variable)) +
  geom_point()

rcorr(df$N,df$FEC,type='pearson')

rcorr(df$FEC,df$G,type='spearman')


### Relationship with FEC
df = data.frame(N = nsp, FEC=kuzfec, G = giK)
df$Fcat = 1
df$Fcat[df$FEC>=250 & df$FEC<500] = 2
df$Fcat[df$FEC>=500 & df$FEC<1000] = 3
df$Fcat[df$FEC>=1000] = 4

df$Fcat = factor(df$Fcat)
table(df$Fcat)


## Plot N et G sur le meême graphe (meme echelle)
ggplot(df,aes(x = N/max(N), y = G, col = Fcat)) + 
  geom_point(size = 4) + geom_smooth(method = 'lm')

ggplot(df,aes(x = Fcat, y = G, col = Fcat)) + 
  geom_boxplot()

ggplot(df,aes(x = Fcat, y = N, col = Fcat)) + 
  geom_boxplot()

summary(glm(N ~ Fcat, family = poisson(link='log'),data= df)) ## no diff with Fcat4

plot(df$N,log(df$FEC+50))


C=ggplot(df, aes(x=N, y=log(FEC+50)))+
  geom_point(color="orange2")+
  geom_smooth(method=lm,color="orange2")+
  xlab('Log(FEC)')+
  ylab("Species richness")+
  ggtitle("")+
  theme_bw()

plot(df$G,log(df$FEC+50))
D=ggplot(df, aes(x=G, y=log(FEC+50)))+
  geom_point(color="orange2")+
  geom_smooth(method=lm,color="orange2")+
  xlab('Log(FEC)')+
  ylab("Gini-Simpson index")+
  ggtitle("")+
  theme_bw()


plot_grid(A,B,C,D, labels=c("A", "B","C","D"), ncol = 2, nrow = 2)
rcorr(as.matrix(df),type = 'spearman')


#==========================================================================================================================================================#
#------------------------------------------------------------------  	Rarecurves on individual data   	---------------------------------------------------------#
#==========================================================================================================================================================#
rm(list=ls())

#####################################################------------- donnees individuelles de l etude de laugier et al -------------############################################## 
par(mfrow=c(3,2)) 

AbCae=read.csv("caecum.csv", header=TRUE,sep=",")

AbCoV=read.csv("ColonV.csv", header=TRUE,sep=",")

AbCoD=read.csv("ColonD.csv", header=TRUE,sep=",")


library(dplyr)
T2 <-bind_rows(AbCae, AbCoV) %>%
  group_by(chevaux) %>%
  summarise_all(funs(sum(., na.rm = TRUE)))

T3 <-bind_rows(T2, AbCoD) %>%
  group_by(chevaux) %>%
  summarise_all(funs(sum(., na.rm = TRUE)))

rm(AbCae,AbCoV,AbCoD, T2)

T3 <- T3[,-1]
raremax <- min(rowSums(T3))
rarecurve(T3, step = 27, xlim=c(1,2000), cex.main=0.95,cex.lab=1.2,xlab="Number of worms collected", ylab="Richness index", col="violetred1", 
          main="A", label=FALSE) 
T3=as.data.frame(T3)
rareslope(T3, sample=1000)
min(rareslope(T3, sample=1000)) # 0
max(rareslope(T3, sample=1000)) # 0.08645256
sp1 <- specaccum(T3)
plot(sp1, col="blue", lwd=2, ci.lty=0, cex.main=0.95,cex.lab=1.2,ci.col="lightblue", xlab="Number of horses sampled", 
     ylab="Richness index", main="B")
rm(T3, sp1, raremax)


#####################################################------------- donnees individuelles de l etude de kuzmina 2005 -------------##############################################
Kuz=read.csv("data_kuzmina2016.csv", header=TRUE,sep=",")
Kuz<- Kuz[,6:38]
raremax <- min(rowSums(Kuz))
rarecurve(Kuz, step = 27, xlim=c(1,2000),cex.main=0.95,cex.lab=1.2, xlab="Number of worms collected", ylab="Richness index", col="violetred1", 
          main="C", label=FALSE) 
rareslope(Kuz, sample=197)
min(rareslope(Kuz, sample=197)) # 0
max(rareslope(Kuz, sample=197)) # 0.03195337

sp1 <- specaccum(Kuz)
plot(sp1, col="blue", lwd=2, ci.lty=0,  cex.main=0.95,cex.lab=1.2,ci.col="lightblue", xlab="Number of horses sampled", 
     ylab="Richness index", main="D")

rm(Kuz, sp1, raremax)

#####################################################------------- donnees individuelles de l etude de Salle 2018 -------------##############################################
dta2018 <- read.csv("data2018.csv", header=TRUE,sep=";")
dta2018<- dta2018[,9:32]
raremax <- min(rowSums(dta2018))
rarecurve(dta2018, step = 27, xlim=c(1,2000), cex.main=0.95,cex.lab=1.2,xlab="Number of worms collected", ylab="Richness index", col="violetred1",
          main="E", label=FALSE) 
rareslope(dta2018,sample=raremax)
sp1 <- specaccum(dta2018)
plot(sp1, col="blue", lwd=2, ci.lty=0, cex.main=0.95,cex.lab=1.2, ci.col="lightblue",xlab="Number of horses sampled", 
     ylab="Richness index", main="F")
rareslope(dta2018, sample=48)
min(rareslope(dta2018, sample=48)) # 0.005037152
max(rareslope(dta2018, sample=48)) # 0.07712526
(rm(dta2018, raremax, sp1))

