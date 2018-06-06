#### boardgames ####

THISFOLDER = "" # define the path to the folder where this file and the scripts file are

setwd(THISFOLDER)
source("advection_functions.R")  # source the functions and packages

#### define parameters ####

mincount=5       # minimal occurrence threshold to include element in advection calculation
topiclength=30   # topic vector length
nfree=0           # how many cores to leave free (0 uses all cores for parallel computing)
smoothing=c(-1,0) # change to 0 to disable topic smoothing

eras = list(1950:1980, 2010:2013) # define the eras (lists of sequences of years)

PATH = "/games" # define folder where the encoded xml boardgame files are (120167-to-120247.xml etc). a handy dump may be found here: https://github.com/inaimathi/all-boardgames-ever
MAINFOLDER = getwd()  # or define folder where to store the temporary files


#### after defining all of the above, run the code starting here ####

#### folders
PERIODFOLDER = file.path(MAINFOLDER,"periods") 
TCMFOLDER = file.path(MAINFOLDER,"tcm")
RELEFOLDER = file.path(MAINFOLDER, "relvectors")
sapply(c(MAINFOLDER, PERIODFOLDER, TCMFOLDER, RELEFOLDER), dir.create)


#### data prep ####

n=164340
years = rep(NA, n)
i = 1; mechs = list();length(mechs)=n

Sys.time()
for(p in paths){
  try({
    bg = xmlToList(xmlParse(p), simplify = F)
    for(b in 1:length(bg)){
      try({
        if(length(bg[[b]]) > 1) { # skip termsofuse, last item; and errors (item not found)
          years[i] = bg[[b]]$yearpublished
          x = unlist(bg[[b]][grep("boardgamemechanic",names(bg[[b]]))])
          mechs[[i]] = as.character(x[names(x)=="boardgamemechanic.text"])
          names(mechs)[i] = gsub(" ","_", bg[[b]]$name$text)
          i=i+1
          if(i %% 1000 == 0) print(i)
        }
      })
    }
  })
}
Sys.time()
save("mechs", "years", file=file.path(MAINFOLDER, "bgtest.RData"))
#system("pmset sleepnow")


#### operate data ####

n2 = which(is.na(years))[1]-1
ok = which(lapply(mechs, length)>0 ) # 58771
mechs2 = mechs[ok]
years2 = years[ok]

plot(sort(table(years2[years2>1900])),xlim=c(1900,2013))
length(years2[years2>1949])

allyears = table(years2)
percounts = numeric()
for(p in eras){
  percounts[as.character(p[1])] = sum(allyears[names(allyears)%in% p ])
}

allmechs = unique(unlist(mechs2[years2>1949]))
bgcountmat = matrix(0, ncol=length(eras), nrow=length(allmechs));rownames(bgcountmat)=allmechs
for(i in 1:ncol(bgcountmat)){
  x = table(unlist(mechs2[years2 %in% eras[[i]] ]))
  bgcountmat[,i] = x[allmechs]/sum(x)*10000 # normalize to per10k
}
bgcountmat[is.na(bgcountmat)]=0
barplot(t(bgcountmat), horiz = T, las=1)
colSums(bgcountmat) # is normalized

# log freq change:
bgfreqdifmat = matrix(NA,nrow=length(allmechs), ncol=ncol(bgcountmat)); rownames(bgfreqdifmat)=allmechs
for(i in 2:ncol(bgfreqdifmat)){
  fr1 = bgcountmat[allmechs,i-1] + 1
  fr2 = bgcountmat[allmechs,i]   + 1
  bgfreqdifmat[,i] = log(fr2/fr1)    #freqdif1
}
bgfreqdifmat[bgcountmat==0] = NA

# save era data
period = mechs2[years2 %in% eras[[1]]]
save(period, file=file.path(PERIODFOLDER, "1.RData"))
period = mechs2[years2 %in% eras[[2]]]
save(period, file=file.path(PERIODFOLDER, "2.RData"))


perioddatas(rdats = list.files(PERIODFOLDER, full.names = T), comlex=NA, minc=mincount, tcmwindow=200, relevancy_threshold=topiclength, nfree=nfree,wordclassesr="", tcmfolder=TCMFOLDER,relefolder=RELEFOLDER, interpol=smoothing, yearnames=c(1,2), windowweights="uniform")
allrelevants = aggregaterelevants(RELEFOLDER)
advection = dotopicsloop(periods=allrelevants, difmat = bgfreqdifmat, countmat=bgcountmat, wordclassesr="", threshold=topiclength, estimate0 = F, estimate2 = F)
# 30 (ie allow max length contexts (23 is max) is slightly better than 10 for both)


cor(advection$inertiamat[,2], bgfreqdifmat[rownames(advection$inertiamat),2], use="pairw")^2
# R2 = 0.3418684


#### tartu plots ####

# boardgames
titl= "Boardgame mechanics, 1950-1980 vs 2010-2013"
ann = "R<sup>2</sup>=0.21 (no smooth)\nR<sup>2</sup>=0.34 (smoothing)"
tartuplots(a=bg_a, f=bg_f, isnew=NULL, isgone=NULL, titl=titl, ann=ann, s=9)



