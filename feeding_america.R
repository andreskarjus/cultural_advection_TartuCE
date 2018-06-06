#### Feeding America ####

THISFOLDER = "" # define the path to the folder where this file and the scripts file are

setwd(THISFOLDER)
source("advection_functions.R")  # source the functions and packages

#### define parameters ####

mincount=25       # minimal occurrence threshold to include element in advection calculation
topiclength=100   # topic vector length
nfree=0           # how many cores to leave free (0 uses all cores for parallel computing)
smoothing=c(-1,0) # change to 0 to disable topic smoothing


FOODCLASSES =  c("fruitvegbeans","meatfishgame","eggscheesedairy","breadsweets","soups","accompaniments","beverages") # which recepies to consider for recepies, from http://archive.lib.msu.edu/dinfo/feedingamerica/cookbook_dtd.txt
ERAS = list(c(1803:1860), c(1890:1913)) # define the eras (lists of sequences of years)

PATH = "/cookbook_textencoded" # important, define folder where the encoded xml cookbooks from https://www.lib.msu.edu/feedingamericadata/ are
MAINFOLDER = getwd()  # or define folder where to store the temporary files


#### after defining all of the above, run the code starting here ####

#### folders
PERIODFOLDER = file.path(MAINFOLDER,"periods") 
TCMFOLDER = file.path(MAINFOLDER,"tcm")
RELEFOLDER = file.path(MAINFOLDER, "relvectors")
sapply(c(PERIODFOLDER, TCMFOLDER, RELEFOLDER), dir.create)


#### data preparation ####

# load and preparse all books
files = list.files(PATH, full.names = T) %>% grep("amem\\.xml|fore\\.xml",., invert = T, value = T) # files, excluding two empty xml files that occur in the corpus
booklist = vector("list", length(files)); years=numeric()
for(i in 1:length(files)){
  booklist[[i]] = xmlRoot(xmlParse(files[i]))
  # get year metadata - inconsistent, will attempt different tags:
  bookid = xmlAttrs(booklist[[i]])["bookID"]
  if(grepl("[0-9]{4}", bookid)){
    idyear = sub("[^0-9]*([0-9]{4})[^0-9]*.*", "\\1", bookid)
    years[i] = as.numeric(idyear)
  } else {
    bookmeta = xmlToList(booklist[[i]][["meta"]]) 
    if(!is.null(bookmeta)){
      if(bookmeta[[1]] != "\n"){
        years[i] = as.numeric(bookmeta[["dcDate"]][[1]])
      }
    }
  }
}

# extract books in chosen eras, extract recepies from books, vectorize ingredients, save as lists:
#ERAS=list(c(1800:1810)) #test
alleras=vector("list", length(ERAS))
for(era in 1:length(ERAS)){
  ibook = which(years %in% ERAS[[era]]) # find which books fall within era
  period = list()
  for(i in ibook){   # for each book in the era, do:
    recipes = xmlElementsByTagName(booklist[[i]], "recipe", T) # fetch all recipes from a book
    # [["body"]][["chapter"]] doesnt work for all books
    # vectorize recepies:
    ingredientlist = vector(mode="list", length=length(recipes))
    for(ir in 1:length(recipes)){
      try({
        if(xmlAttrs(recipes[[ir]])["class1"] %in% FOODCLASSES){
          tmprec = xmlElementsByTagName(recipes[[ir]], "ingredient", T)
          if(length(unique(tmprec))>1){
            try({
              # do basic cleaning, stemming, lowercasing and such:
              unl = sapply(1:length(tmprec), function(x) xmlToList(tmprec[[x]][[1]])) %>%
                unique() %>%
                gsub("[[:punct:][:space:][:digit:]]", "", .) %>% 
                tolower()
              # simple stemming:
              unl = ifelse(grepl("ies$", unl), sub("ies$", "y", unl),
                           ifelse(grepl("oes$", unl), sub("oes$", "o", unl),
                                  ifelse(grepl("leaves$", unl), sub("leaves$", "leaf", unl),
                                         ifelse(grepl("s$", unl), sub("s$", "", unl), unl)))
              )
              unl = unl[nchar(unl)>1]  # remove abnormally short and zero-length ingredients
              ingredientlist[[ir]] = unl
            })
          }
        }
      })
    }
    ingredientlist = ingredientlist[which(sapply(ingredientlist, length)>1)]
    period = c(period, ingredientlist)
  }
  # save period list into rdata for later processing
  alleras[[era]] = period
  save(period,file=file.path(PERIODFOLDER, paste0(era, ".RData")) )
  # .... file.path(PERIODFOLDER, paste0(era, ".RData"))
}

# count matrices
t1 = table(unlist(alleras[[1]]))
t2 = table(unlist(alleras[[2]])) # t2=t2[-1]
trms = union(names(t1), names(t2))
countmat=matrix(NA,ncol=2, nrow=length(trms), dimnames = list(trms, c(1,2)) )
countmat[,1] = t1[rownames(countmat)]
countmat[,2] = t2[rownames(countmat)]
countmat[is.na(countmat)]=0
erasms = colSums(countmat)
# normalize to 100k:
countmat[,1] = (countmat[,1]/erasms[1])*100000  # normalize counts to per-100000
countmat[,2] = (countmat[,2]/erasms[2])*100000
countmat = countmat[nchar(rownames(countmat))>0, ] # doublecheck to avoid empty strings
# </hacky counting

perioddatas(rdats = list.files(PERIODFOLDER, full.names = T), comlex=NA, minc=mincount, tcmwindow=200, relevancy_threshold=topiclength, nfree=3,wordclassesr="", tcmfolder=TCMFOLDER,relefolder=RELEFOLDER, interpol=smoothing, yearnames=c(1,2), windowweights="uniform")
# change interpol=0 to use the no-smoothing version

allrelevants = aggregaterelevants(RELEFOLDER)
freqdifmat = dofreqdif(countmat, smooth=1, lexicons=rownames(countmat))
advection = dotopicsloop(periods=allrelevants, difmat = freqdifmat, countmat=countmat, wordclassesr="", threshold=topiclength, estimate0 = F, estimate2 = F)

print(cor(advection$inertiamat[,2], freqdifmat[rownames(advection$inertiamat),2], use="pairw")^2)

fa_a=advection$inertiamat[!is.na(advection$inertiamat[,2]),2]
fa_f=freqdifmat[names(fa_a),2]
fa_r = residuals(lm(fa_f~fa_a))
fa_isnew=countmat[names(fa_a),1]==0
fa_isgone=countmat[names(fa_a),2]==0

titl= "Cookbook ingredients 1803-1860 vs 1890-1913"
ann = "R<sup>2</sup>=0.23 (no smooth)\nR<sup>2</sup>=0.79 (smoothing)"
tartuplots(a=fa_a, f=fa_f, isnew=fa_isnew, isgone=fa_isgone, titl=titl, ann=ann)


