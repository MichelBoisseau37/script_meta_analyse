
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#####========= HSMC analysis
setwd('~/Documents/CYATHOMIX/WP3/MetaAnalysis/')
require(lme4)
require(fitdistrplus)

### rarecurves
require(vegetarian)
require(vegan)
require(indicspecies)
###------============= NORMANDY data for niche partitioning + network analysis =========---------------
cae = read.csv(file = './caecum.csv',header=T)
cae = cae[order(cae$chevaux),]
cv =  read.csv(file = './ColonV.csv',header=T)
cv = cv[order(cae$chevaux),]
cd =  read.csv(file = './ColonD.csv',header=T)
cd = cd[order(cae$chevaux),]

cl = cd[,-1] + cv[,-1] + cae[,-1]
cl$G.strongles = NULL
cl$oxyures = NULL
cl$chevaux = cae$chevaux

cl = cl[which(rowSums(cl[,-dim(cl)[2]])!=0),]
cl$Cyathos. = NULL
cl$Cylicodon. = NULL
cl$Anoploce. = NULL
cl$P.strongles = NULL
cl$Poteriosto. = NULL
cl$Triodonto. = NULL

cl[is.na(cl)] <- 0 ## replace missing species

clenv = read.csv(file = './Age_samplingCL.csv',header=T,sep=';')
clenv$chevaux = clenv$Nom

dim(cld)
#[1] 41 25

cld = merge(clenv,cl,by='chevaux')
clenv = cld[,c(1,3)]
cldsp = cld[,-c(1:5)]

### Retain horses with parasites
clenv = clenv[which(rowSums(cldsp)!=0),]
cldsp = cldsp[which(rowSums(cldsp)!=0),]
dim(cldsp)
#[1] 37 20

### Retain species with 5% prevalence in the retained set (parasite positive horses)
ny = dim(cldsp)[1]
ny
#[1] 37
prev = colSums(1*(cldsp>0))
# C.nassatus   C.coronatum   C.leptosto.    C.labratum    C.labiatum     C.insigne       C.goldi   C.tetracan.   C.calicatus 
#         24            26             4             7            14            19            16             9             7 
# C.minutus   C.longibur.   C.catinatum   C.ultrajec.   C.elongatus   C.pateratum   C.bicorona.  P.impariden.      P.ratzii 
#         9             5            23             9             2             3             3            13             2 
# T.tenuicollis    T.serratus 
#             1             8
selSP = prev>=ny*0.1
names(which(selSP))
# [1] "C.nassatus"   "C.coronatum"  "C.leptosto."  "C.labratum"   "C.labiatum"   "C.insigne"    "C.goldi"      "C.tetracan." 
# [9] "C.calicatus"  "C.minutus"    "C.longibur."  "C.catinatum"  "C.ultrajec."  "P.impariden." "T.serratus" 

ns = sum(selSP)
ns
#[1] 15

#####-----------------################## HMSC - Prev 10% - Test if cooccurrence is affected by niche or age
require(Hmsc) ## 
require(corrplot)
require(parallel)
### Hmsc R code is reused from Mol Ecol paper available under the following links:
#Article:https://doi.org/10.1111/mec.15516 
#https://datadryad.org/stash/dataset/doi:10.5061/dryad.9kd51c5dp

# BUILD DATASETS ##################################################################

# XData
env = data.frame(niche = factor(rep(c('cd','cv','cae'), each = dim(clenv)[1])), 
                 age = rep(clenv$âge.en.mois,3))
env$yr = 1
env$yr[env$age/12>2 & env$age/12<10]=2
env$yr[env$age/12>=10]=3

which(is.na(env$age))
#[1]  3 40 77

### Species count matrix
ni = rbind(cd,cv,cae) 
ni$niche = factor(rep(c('cd','cv','cae'), each =dim(cd)[1]))

ni$G.strongles = NULL
ni$oxyures = NULL
ni$Cyathos. = NULL
ni$Cylicodon. = NULL
ni$Anoploce. = NULL
ni$P.strongles = NULL
ni$Poteriosto. = NULL
ni$Triodonto. = NULL

### Keep horses with postivie burden only
dim(ni)
#[1] 138  22
ni = ni[ni$chevaux %in% clenv$chevaux,]
dim(ni)
#[1] 111  22

### Add infection level for each horse
sumByniche=data.frame(chevaux = as.character(ni$chevaux),
                      wbni = as.integer(rowSums(ni[,-c(1,dim(ni)[2])])))
wb = aggregate(wbni ~ chevaux, FUN=sum,data=sumByniche)
env$wb = rep(wb$wbni,3)

### Species count matrix
nisp = ni[,-c(1,dim(ni)[2])]

### Filter for 10% prevalence species
nisp = nisp[,selSP]
dim(nisp)
#[1] 111  15
nisp = nisp[which(!is.na(env$age)),]
dim(nisp)
#[1] 108  15
env = na.omit(env) ## 1 horse with missing age

env$age = NULL
dim(env)
#[1] 108   3


Nisp.abu = nisp
Nisp.abu[nisp==0] = NA
Nisp.abu = log(Nisp.abu)
for (i in 1:ns){
  Nisp.abu[,i] = Nisp.abu[,i] - mean(Nisp.abu[,i],na.rm=TRUE)
  Nisp.abu[,i] = Nisp.abu[,i] / sqrt(var(Nisp.abu[,i],na.rm=TRUE))
}

### Random effect and study design
studyDesign = data.frame(ID = as.factor(as.character(1:dim(nisp)[1])), Horse = as.factor(rep(seq(1:(dim(nisp)[1]/3)),3)))
rL.Horse = HmscRandomLevel(units = studyDesign$Horse)
rL.ID = HmscRandomLevel(units = studyDesign$ID)

# SET DIRECTORIES AND PARAMETERS ##################################################################

localDir = "."
dataDir = file.path(localDir, ".")
ModelDir = file.path(localDir, "models")
MixingDir = file.path(localDir, "mixing")

samples = 1000
nChains = 4
dataset = 1

# RUN MODELS ##################################################################
#
##### Record *samples* posterior samples while skipping *thin* MCMC step between samples
## from *nChains* chains after discarding the first *transient* MCMC steps

for (thin in c(1, 10, 100)){
  for (modeltype in 1:2){ ## 1 = Pres/Abs vs. 2 = Abundance
    for (model in 1:2){
      if(model==1){ ## Raw data
        XFormula =~ wb ## Infection level
      }
      if(model==2){
        XFormula =~ wb + niche + yr ## Age and niche correction
      } 
      print(paste0("thin = ",as.character(thin),", modeltype = ",c("pa","abu")[modeltype],
                   ", model = ",as.character(model)))
      
      set.seed(1)
      m = Hmsc(Y = if(modeltype==1){1*(nisp>0)} else {Nisp.abu},
               XData = env,  XFormula = XFormula,
               distr=if(modeltype==1){"probit"} else {"normal"},
               studyDesign = studyDesign, 
               ranLevels = list(ID=rL.ID, Horse = rL.Horse))
      
      ptm = proc.time()
      m = sampleMcmc(m, samples = samples, thin = thin, verbose = samples,
                     adaptNf = rep(ceiling(0.4*samples*thin),2),
                     transient = ceiling(0.5*samples*thin),
                     nChains = nChains, nParallel = nChains)
      computational.time =  proc.time() - ptm
      
      print(computeWAIC(m))
      
      filename = file.path(ModelDir, paste("dataset_",as.character(dataset),"_",c("pa","abundance")[modeltype],
                                           "_model",as.character(model),
                                           "_thin_", as.character(thin),"_samples_", as.character(samples),
                                           "_chains_",as.character(nChains),
                                           ".Rdata",sep = ""))
      save(m,file=filename,computational.time)
    }
  }
}

# MIXING STATISTICS ####################################################################

samples = 1000
nChains = 4
thin_max = 100

for (modeltype in 1:2){
  for (model in 1:2){
    for(thin in c(1, 10, 100)){ #= thin_max
      cat(modeltype,dataset,model,thin)
      cont = TRUE
      while(cont){
        filename = file.path(ModelDir, paste("dataset_",as.character(dataset),
                                             "_",c("pa","abundance")[modeltype],
                                             "_model",as.character(model),
                                             "_thin_", as.character(thin),
                                             "_samples_", as.character(samples),
                                             "_chains_",as.character(nChains),
                                             ".Rdata",sep = ""))
        if(file.exists(filename)){
          load(filename)
          
          mpost = convertToCodaObject(m)
          
          es.beta = effectiveSize(mpost$Beta)
          ge.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
          
          es.gamma = effectiveSize(mpost$Gamma)
          ge.gamma = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
          
          es.V = effectiveSize(mpost$V)
          ge.V = gelman.diag(mpost$V,multivariate=FALSE)$psrf
          
          mpost$temp = mpost$Omega[[1]]
          for(i in 1:length(mpost$temp)){
            mpost$temp[[i]] = mpost$temp[[i]][,1:ns^2]
          }
          es.omega = effectiveSize(mpost$temp)
          ge.omega = gelman.diag(mpost$temp,multivariate=FALSE)$psrf
          
          
          mixing = list(es.beta=es.beta, ge.beta=ge.beta,
                        es.gamma=es.gamma, ge.gamma=ge.gamma,
                        es.V=es.V, ge.V=ge.V,
                        es.omega=es.omega, ge.omega=ge.omega)
          
          filename = file.path(MixingDir, 
                               paste("dataset_",as.character(dataset),"_",
                                     c("pa","abundance")[modeltype],"_thin",as.character(thin),
                                     "_model",as.character(model),".Rdata",sep = ""))
          
          save(file=filename, mixing)
          
          cont=FALSE
        } else {
          thin = thin/10}
        if (thin<1){
          print("not found")
          cont=FALSE
        }
      }
    }
  }
}

# PLOT MIXING STATISTICS (START) ####################################################################
thin = thin_max
for (zzz in 1:2){
  # tiff(filename=c(paste0("./effective_sample_size.tif"),
  #                 paste0("./psrf.tif"))[zzz], 
  #      res=300, unit="cm", height=17, width=17)
  pdf(file = c(paste0("./effective_sample_size.pdf"),
                                paste0("./psrf.pdf"))[zzz])
  par(mfrow=c(2,3))
  for (modeltype in 1:2){
    for (model in 1:2){
      for(thin in c(1,10,100)){
        filename = file.path(MixingDir, paste("dataset_",as.character(dataset),"_",
                                              c("pa","abundance")[modeltype],
                                              "_thin",as.character(thin),
                                              "_model",as.character(model),".Rdata",sep = ""))
        load(filename)
        resa=list()
        for (param in 1:2){
          if(param == 1){
            if(zzz==1){
              resa[[param]]=as.vector(mixing$es.beta)
            } else {
              resa[[param]]=as.vector(mixing$ge.beta)
            }
          } else {
            if(zzz==1){
              resa[[param]]=as.vector(mixing$es.omega)
            } else {
              resa[[param]]=as.vector(mixing$ge.omega)
            } 
          }
        }
        boxplot(resa, main = paste0(c("P-A","Abu")[[modeltype]],
                                    " (model ",as.character(model),")",
                                    " - thin", as.character(thin)),
                names=c("beta","omega"),
                ylim=if(zzz==1){c(0,6000)} else {c(0,2)})
      }
    }
  }
  dev.off()
}

# PARAMETERS ESTIMATES ####################################################################
sp = names(which(selSP))

plotVP =
  function (hM, VP, cols=NULL, main = "Variance Partitioning", ...)
  {
    ng = dim(VP$vals)[1]
    if(is.null(cols)){
      cols = heat.colors(ng, alpha = 1)
    }
    leg = VP$groupnames
    for (r in 1:hM$nr) {
      leg = c(leg, paste("Random: ", hM$rLNames[r], sep = ""))
    }
    means = round(100 * rowMeans(VP$vals), 1)
    for (i in 1:ng) {
      leg[i] = paste(leg[i], " (mean = ", toString(means[i]),
                     ")", sep = "")
    }
    VP2 = reshape2::melt(VP$vals)
    p=ggplot(VP2,aes(x = Var2, y = value,fill=Var1)) +
      geom_bar(stat='identity') + theme_bw() +
      scale_fill_manual(values = viridis::viridis_pal(option='D')(length(levels(VP2$Var1))),
                        labels = leg) +
      xlab('Species') + ylab('Variance proportion') +
      theme(legend.position = 'bottom',legend.title = element_blank(),
            axis.text.x = element_text(angle=45,vjust = 0.5))
    print(p)
  }

for (modeltype in 1:2){
  for (model in c(1,2)){
    thin = thin_max
    cont = TRUE
    while(cont){
      filename = file.path(ModelDir, paste("dataset_",as.character(dataset),"_",c("pa","abundance")[modeltype],
                                           "_model",as.character(model),
                                           "_thin_", as.character(thin),"_samples_", as.character(samples),
                                           "_chains_",as.character(nChains),
                                           ".Rdata",sep = ""))
      if(file.exists(filename)){
        load(filename)
        cont=FALSE
      } else {
        thin=thin/10}
      if (thin<1){
        print("not found")
        cont=FALSE
      }
    }
    
    
    OmegaCor = computeAssociations(m)
    supportLevel = 0.95
    toPlot = ((OmegaCor[[1]]$support > supportLevel) + (OmegaCor[[1]]$support < (1-supportLevel)) > 0) * OmegaCor[[1]]$mean
    
    pos = (sum((OmegaCor[[1]]$support > supportLevel))-m$ns)/(m$ns*(m$ns-1))
    neg = mean((OmegaCor[[1]]$support < (1-supportLevel)))
    
    cat("modeltype = ",modeltype,"dataset = ",dataset, "model = ", model,"pos = ",pos, "neg = ",neg,"\n")
    #    if(dataset==4){
    
    if(TRUE){
      predY = computePredictedValues(m, expected=TRUE)
      MF = evaluateModelFit(hM=m, predY=predY)
      if(modeltype==1){
        ta = cbind(MF$AUC,MF$TjurR2)
        colnames(ta) = c("AUC","TjurR2")
      } else {
        ta = cbind(MF$R2)
        colnames(ta) = c("R2")
      }
      filnam = paste0("res/MF_",c("pa","abu")[modeltype],"_coc_",c("raw_","res_")[model])
      row.names(ta) = m$spNames
      write.csv(ta,file = paste0(filnam,".csv"))
    }
    
    plotOrder = sp[corrMatOrder(toPlot[sp,sp], order = "AOE")]
    
    filnam = paste0("res/",c("pa","abu")[modeltype],
                    "_coc_",c("raw_","res_")[model],"order_groups")
    tiff(filename=paste0(filnam,".tiff"),res=300, unit="cm", height=17, width=17)
    corrplot(toPlot[sp,sp],method='circle',
             col=viridis::viridis_pal(option='D')(200),tl.col = 'darkblue',type='upper')
    dev.off()
    namelist=as.data.frame(m$spNames)
    write.csv(namelist,file = paste0(filnam,".csv"))
    
    ### Save corr matrix
    assign(paste0('corr_modelT',as.character(modeltype),
                  '_model',as.character(model)),toPlot[sp,sp])
    
    # head(m$X)
    #   (Intercept)   wb nichecd nichecv yr
    # 1           1  450       1       0  3
    # 2           1  950       1       0  3
    # 4           1   40       1       0  2
    # 5           1 2420       1       0  2
    # 6           1  390       1       0  1
    # 7           1 1500       1       0  2
    if (model==2){
      groupnames = c("WB","Niche","Age")
      group = c(1,1,2,2,3)
      VP = computeVariancePartitioning(m,group = group,groupnames = groupnames)
      filnam = paste0("res/",c("pa","abu")[modeltype],"_VP_",c("raw","res")[model])
      #tiff(filename=paste0(filnam,".tiff"),res=300, unit="cm", height=17, width=17)
      pdf(file = paste0(filnam,".pdf"),width = 14,height =8)
      plotVP(m,VP)
      dev.off()
      
      # VP = computeVariancePartitioning(m)
      # filnam = paste0("res/",c("pa","abu")[modeltype],"_VP_not_grouped_",c("raw_","res_")[model])
      # tiff(filename=paste0(filnam,".tiff"),res=300, unit="cm", height=17, width=17)
      # plotVariancePartitioning(m,VP)
      # dev.off()
    }
    
    p = (sum(OmegaCor[[1]]$support > supportLevel)-ns)/(ns*(ns-1))
    
    res = mean((OmegaCor[[1]]$support < (1-supportLevel)))
    cat(p,res,'\n')
    # 
    # cat("p11,p22,p33,p12,p13,p23 = ",res,"\n")
    # n11 = mean((OmegaCor[[1]]$support[sp1,sp1] < (1-supportLevel)))
    # n22 = mean((OmegaCor[[1]]$support[sp2,sp2] < (1-supportLevel)))
    # n33 = mean((OmegaCor[[1]]$support[sp3,sp3] < (1-supportLevel)))
    # n12 = mean((OmegaCor[[1]]$support[sp1,sp2] < (1-supportLevel)))
    # n13 = mean((OmegaCor[[1]]$support[sp1,sp3] < (1-supportLevel)))
    # n23 = mean((OmegaCor[[1]]$support[sp2,sp3] < (1-supportLevel)))
    # res = c(n11,n22,n33,n12,n13,n23)
    # cat("n11,n22,n33,n12,n13,n23 = ",res,"\n")
  }
}

###--------======= Combine corr matrix for pa/abu models

### Residual vs. raw - Presence/Absence
#corr_modelT1_model1 raw
#corr_modelT1_model2 residual
corr_pa = corr_modelT1_model1
# replace lower values by residual cooccurences
corr_pa[lower.tri(corr_pa)] <- corr_modelT1_model2[lower.tri(corr_modelT1_model2)]


### Residual vs. raw - Abundance
#corr_modelT2_model1 raw
#corr_modelT2_model1 residual
corr_abu = corr_modelT2_model1
# replace lower values by residual cooccurences
corr_abu[lower.tri(corr_abu)] <- corr_modelT2_model2[lower.tri(corr_modelT2_model2)]

### #Give species their full names
rownames(corr_pa) = gsub('P.impariden.','P.imparidentatum', rownames(corr_pa))
rownames(corr_pa) = gsub('C.ultrajec.','C.ultrajectinus', rownames(corr_pa))
rownames(corr_pa) = gsub('C.longibur.','C.longibursatus', rownames(corr_pa))
rownames(corr_pa) = gsub('C.tetracan.','C.tetracanthum', rownames(corr_pa))
rownames(corr_pa) = gsub('C.leptosto.','C.leptostomum', rownames(corr_pa))
colnames(corr_pa) = gsub('P.impariden.','P.imparidentatum', colnames(corr_pa))
colnames(corr_pa) = gsub('C.ultrajec.','C.ultrajectinus', colnames(corr_pa))
colnames(corr_pa) = gsub('C.longibur.','C.longibursatus', colnames(corr_pa))
colnames(corr_pa) = gsub('C.tetracan.','C.tetracanthum', colnames(corr_pa))
colnames(corr_pa) = gsub('C.leptosto.','C.leptostomum', colnames(corr_pa))

rownames(corr_abu) = gsub('P.impariden.','P.imparidentatum', rownames(corr_abu))
rownames(corr_abu) = gsub('C.ultrajec.','C.ultrajectinus', rownames(corr_abu))
rownames(corr_abu) = gsub('C.longibur.','C.longibursatus', rownames(corr_abu))
rownames(corr_abu) = gsub('C.tetracan.','C.tetracanthum', rownames(corr_abu))
rownames(corr_abu) = gsub('C.leptosto.','C.leptostomum', rownames(corr_abu))
colnames(corr_abu) = gsub('P.impariden.','P.imparidentatum', colnames(corr_abu))
colnames(corr_abu) = gsub('C.ultrajec.','C.ultrajectinus', colnames(corr_abu))
colnames(corr_abu) = gsub('C.longibur.','C.longibursatus', colnames(corr_abu))
colnames(corr_abu) = gsub('C.tetracan.','C.tetracanthum', colnames(corr_abu))
colnames(corr_abu) = gsub('C.leptosto.','C.leptostomum', colnames(corr_abu))

### Save to figure
dfpa = reshape2::melt(corr_pa) ##Var1 = x, Var2 = y
dfpa$value[dfpa$value==0]=NA

p1 = ggplot(dfpa,aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = viridis_pal(option='D')(100),
                       na.value = 'white') +
  xlab('') + ylab('') +
  theme(legend.position = 'bottom',text = element_text(size = 12),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2),
        legend.text = element_text(size = 8, angle=45, vjust = .5),
        legend.title = element_blank())+
  labs('Co-occurrence')

dfabu = reshape2::melt(corr_abu) ##Var1 = x, Var2 = y
dfabu$value[dfabu$value==0]=NA

p2 = ggplot(dfabu,aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = viridis_pal(option='C')(100),
                       na.value = 'white') +
  xlab('') + ylab('') +
  theme(legend.position = 'bottom',text = element_text(size = 12),
        axis.text.x = element_text(angle=90,hjust=0.95,vjust=0.2),
        legend.text = element_text(size = 8, angle=45, vjust =.5),
        legend.title = element_blank()) +
  labs(colour = 'Co-occurrence')

pdf(file = './Figure4.pdf', width = 14, height= 8)
multiplot(p1,p2,cols=2)
dev.off()
