#### IMDb ####

THISFOLDER = "" # define the path to the folder where this file and the scripts file are

setwd(THISFOLDER)
source("advection_functions.R")  # source the functions and packages

#### define parameters ####

mincount=25       # minimal occurrence threshold to include element in advection calculation
topiclength=100   # topic vector length
nfree=0           # how many cores to leave free (0 uses all cores for parallel computing)
smoothing=c(-1,0) # change to 0 to disable topic smoothing

MINVOTES=100
EXCLUDE = c("Documentary", "Reality-TV","Short","Talk-Show","Musical","News", "Music", "Game-Show", "Animation", "Adult", "Reality-tv", "Lifestyle")
INCLUDE = "Drama"
COUNTRY = "USA"
ERAS = list(c(1960:1979), c(2010:2018)) # define the eras (lists of sequences of years)

PATH = "/imdb_dump13.11.17" # important, define folder where the imdb dump files are (eg from ftp://ftp.funet.fi/pub/mirrors/ftp.imdb.com/pub/frozendata/)
MAINFOLDER = getwd()  # or define folder where to store the temporary files


#### after defining all of the above, run the code starting here ####

#### folders
PERIODFOLDER = file.path(MAINFOLDER,"periods") 
TCMFOLDER = file.path(MAINFOLDER,"tcm")
RELEFOLDER = file.path(MAINFOLDER, "relvectors")
sapply(c(MAINFOLDER, PERIODFOLDER, TCMFOLDER, RELEFOLDER), dir.create)


#### filter countries ####

Sys.setlocale('LC_ALL','C') 
countries = readLines(file.path(PATH, "countries.txt"))
countries = matrix(unlist(strsplit(countries, "\t{1,}")), ncol = 2, byrow = TRUE)
countries = countries[-grep("\\(T{0,1}VG{0,1}\\)|([{}])", countries[,1], value = F),]
nrow(countries)
usa = countries[countries[,2]==COUNTRY,1] # get only usa movies

#### ratings ####
ratings = readLines(file.path(PATH, "ratings.txt"))
rat2 = matrix(unlist(strsplit(ratings, " {2,}")[-1]), ncol = 5, byrow = TRUE)
# rerun from here if test sample before:
rat3 = data.frame(n=as.numeric(rat2[,3]),score=as.numeric(rat2[,4]), names = rat2[,5], stringsAsFactors = F)
rat3$year=as.numeric(gsub("[^(]{4,}\\(([0-9]{4})\\).*", "\\1", rat3$names))
rat3 = rat3[!is.na(rat3$year),]
rat3 = rat3[-grep("\\(T{0,1}VG{0,1}\\)|([{}])", rat3$names, value = F),]

rat3 = rat3[rat3$names %in% usa, ] # get country
rat3 = rat3[rat3$n>=MINVOTES,]

#head(rat3)
#hist(rat3$year,100);summary(rat3$n)

# drop unneeded movies
rat3 = rat3[rat3$year %in% unlist(ERAS), ]
#rat3 = rat3[rat3$year %in% c(1960,2016), ]  # test subset, comment out ----
#rat3 = rat3[rat3[,"n"]>5000,]               # test subset
nrow(rat3) 
# ---

rat3 = rat3[,c("names", "year")]   # drop other unneeded columns

#### genres ####
genres = readLines(file.path(PATH, "genres.txt"))
genres = matrix(unlist(strsplit(genres, "\t{1,}")), ncol = 2, byrow = TRUE)
genres = genres[-grep("\\(T{0,1}VG{0,1}\\)|([{}])", genres[,1], value = F),]
genrelist = split(genres[,1], genres[,2])
#unique(genres[,2]) # 31 genres after {} removal

# filter for certain genre, and exclude TV-specific stuff, pr0n, cartoons and whatnot.
uExclude = unique(unlist(genrelist[EXCLUDE]))
uInclude = unique(unlist(genrelist[INCLUDE]))
okgenre = setdiff(uInclude, uExclude)

rat4 = rat3[which(rat3$names %in% okgenre),]  # this will serve as basis for getting eras and keywords
nrow(rat4)



#### keywords ####

kw = readLines(file.path(PATH, "keywords.txt"))
kwm = matrix(unlist(strsplit(kw, "\t{1,}")), ncol = 2, byrow = TRUE)
kwm = kwm[-grep('\\(T{0,1}VG{0,1}\\)|([{}])|^\\"', kwm[,1], value = F),]
kwm = kwm[kwm[,1] %in% rat4$names,]
nrow(kwm)
kwlist = split(kwm[,2], kwm[,1]) # named list, each vector is tags for one movie
length(kwlist)

##### eras ####
# quoted tv series; current code should leave out those even when slipped in from previous files, since no match with keyword file, where they are quoted.
alleras=vector("list", length(ERAS))
for(era in 1:length(ERAS)){
  pmovies = rat4$names[rat4$year %in% ERAS[[era]] ]
  period = kwlist[pmovies]
  period = period[sapply(period, length)>=5]
  alleras[[era]] = period
  save(period,file=file.path(PERIODFOLDER, paste0(era, ".RData")) )
}


# counts
t1 = table(unlist(alleras[[1]]))
t2 = table(unlist(alleras[[2]]))
trms = union(names(t1), names(t2))
countmat=matrix(NA,ncol=2, nrow=length(trms), dimnames = list(trms, c(1,2)) )
countmat[,1] = t1[rownames(countmat)]
countmat[,2] = t2[rownames(countmat)]
countmat[is.na(countmat)]=0
erasms = colSums(countmat)
# normalize to 100k:
countmat[,1] = (countmat[,1]/erasms[1])*100000   # normalize to per-100000
countmat[,2] = (countmat[,2]/erasms[2])*100000
countmat = countmat[nchar(rownames(countmat))>0, ] # doublecheck to avoid empty strings
# </hacky counting


perioddatas(rdats = list.files(PERIODFOLDER, full.names = T), comlex=NA, minc=mincount, tcmwindow=1000, relevancy_threshold=topiclength, nfree=nfree,wordclassesr="", tcmfolder=TCMFOLDER,relefolder=RELEFOLDER, interpol=smoothing, yearnames=c(1,2), windowweights="uniform")

allrelevants = aggregaterelevants(RELEFOLDER)
freqdifmat = dofreqdif(countmat, smooth=1, lexicons=rownames(countmat))
advection = dotopicsloop(periods=allrelevants, difmat = freqdifmat, countmat=countmat, wordclassesr="", threshold=topiclength, estimate0 = F, estimate2 = F)


#### results&plots ####

cor(advection$inertiamat[,2], freqdifmat[rownames(advection$inertiamat),2], use="pairw")^2

i_a=advection$inertiamat[!is.na(advection$inertiamat[,2]),2]
i_f=freqdifmat[names(i_a),2]

i_isnew=countmat[names(i_a),1]==0
i_isgone=countmat[names(i_a),2]==0

# imdb
titl= "US drama genre movies keywords, 1960-1979 vs 2010-2018"
ann = "R<sup>2</sup>=0.18 (no smooth)\nR<sup>2</sup>=0.7 (smoothing)"
tartuplots(a=i_a, f=i_f, isnew=i_isnew, isgone=i_isgone, titl=titl, ann=ann)
# nosmooth 0.176951 (0.1127138 before weighting fix)
# smooth 0.7011762  (0.6729361 before weighting fix)