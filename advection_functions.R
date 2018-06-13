#
# Functions for the scripts to calculate the topical-cultural advection values on various datasets.
# Andres Karjus, University of Edinburgh
#

# The following packages are required; this bit will attempt to install the missing ones.
p = setdiff(c("visNetwork", "igraph", "XML", "magrittr", "plotly", "crosstalk", "parallel", "fastmatch", "text2vec", "grr", "Matrix", "compiler"), rownames(installed.packages()))
if(length(p)>0) {install.packages(p)}

# Note: text2vec is in beta, the scripts here used v0.5.0.

# load packages
library(visNetwork)
library(igraph)
library(XML)
library(magrittr)
library(plotly)
library(crosstalk)
library(parallel)
library(fastmatch)
library(text2vec)
library(grr)
library(Matrix)
library(compiler)


# functions

perioddatas = function(rdats = list.files(paste0(FOLDER,"cohaperiods"), full.names = T), comlex, minc, tcmwindow, windowweights="weighted", relevancy_threshold, nfree=0,wordclassesr="", tcmfolder=paste0(FOLDER,"cohaperiodtcm/"),relefolder=paste0(FOLDER,"cohaperiodrelevants/"), interpol=0, yearnames=seq(1810,2000,10)){
  library(parallel)
  nc = detectCores() - nfree
  print(c(paste(Sys.time(), "start tcm phase")))
  
  # TCMs
  if(length(rdats)>0){
    cl = makeCluster(nc)
    clusterExport(cl, c("rdats","comlex", "minc", "tcmwindow","wordclassesr", "interpol", "yearnames", "tcmfolder", "windowweights"),envir = environment())  # from local
    clusterExport(cl, c("make_ppmi_tcm", "makeppmi")) # from global
    clusterEvalQ(cl, c(library(fastmatch), library(igraph), library(text2vec), library(grr), library(Matrix)))
    
    tryCatch({
      parLapply(cl, 1:length(rdats), function(x){
        # load(x) # load period data -> in function instead; send just period index
        # do and save tcm
        # NB: assumes ordered period files!
        tcm = make_ppmi_tcm(nperiod=x,rdats=rdats, comlex=comlex, minc=minc, gramwin = tcmwindow, wordclassesr=wordclassesr, wordclassesc="", interpol=interpol, windowweights=windowweights)
        save(tcm, file=file.path(tcmfolder, paste0(yearnames[x],".RData" )))
        rm(period);gc() # remove period data
      }
      )
    }, error=function(e){print(e)}, finally = stopCluster(cl))  
    gc()
  }
  print(c(paste(Sys.time(), "TCMs done")))
  
  # get relevance sets
  if(length(tcmfolder)>0){
    tcmdats = list.files(tcmfolder, full.names = T)
    cl = makeCluster(nc)
    clusterEvalQ(cl,  c(library(grr), library(Matrix)))
    clusterExport(cl, c("dorelevants"))
    clusterExport(cl, c("relefolder","relevancy_threshold"),envir = environment())  # from local
    tryCatch({
      parLapply(cl, tcmdats, function(x){
        load(x)
        relevants = dorelevants(tcm,relevancy_threshold=relevancy_threshold)
        names(relevants) = rownames(tcm)
        save(relevants, file=file.path(relefolder, paste0(gsub("[^0-9]","",x),".RData" )))
        rm(tcm);gc()
      })
    }, error=function(e){print(e)}, finally = stopCluster(cl)) 
    gc()
  }
  print(c(paste(Sys.time(), "relevance vectors done")))
  
}


make_ppmi_tcm = function(nperiod, rdats, minc, gramwin, comlex, wordclassesr, wordclassesc= ".", interpol=c(0), 
                         simuversion=F, maxperiod=NA,segment_indices=NA, 
                         ppmi=T, windowweights ){   
  #  3 for simu only; simuversion is only for compatibility with ac->spok change simulation
  # NB, period files are up to ~400mb, so concatenating will eat RAM fast...
  library(text2vec)
  library(Matrix)
  if(!simuversion){
    periods=list()
    for(ni in (nperiod+interpol) ){ # 0 is current, -1 previous, 1 future, etc
      if(ni > 0 & !(ni > length(rdats)) ){ # avoid indices beyond list
        load(rdats[ni]) # load period
        periods = append(periods, period)
        #print(c(length(period), length(periods)) ) # debug
      }
    }
    rm(period) # clean lastly loaded text data object
  }
  
  if(simuversion){
    periods=list()
    for(ni in (nperiod+interpol) ){ # 0 is current, -1 previous, 1 future, etc
      if(ni > 0 & !(ni > maxperiod) ){ # avoid indices beyond list
        period = c(rdats$academic[segment_indices[[1]][[ni]] ], rdats$spoken[segment_indices[[2]][[ni]] ])
        periods = append(periods, period)
        #print(c(length(period), length(periods)) ) # debug
      }
    }
    rm(period,rdats);gc()
  }
  
  # do iter,vocab, tcm
  it = itoken(periods, progressbar = FALSE) # period list from load()
  rm(periods); gc() # clean periods - iterator envir now contains all data
  vocab2 <- create_vocabulary(it, ngram = c(ngram_min = 1L, ngram_max = 1L))
  #docmin = vocab2$vocab$doc_counts/vocab2[["document_count"]]
  vocab2 <- prune_vocabulary(vocab2, term_count_min = minc, doc_proportion_max=1, doc_proportion_min=0)# 2/vocab2$document_count) # avoid single-doc words
  
  vectorizer <- vocab_vectorizer(vocab2)
  if(windowweights == "uniform"){
    tcm = create_tcm(it, vectorizer, skip_grams_window = gramwin, 
                     weights = 1 / rep(1,gramwin) )
  } else {
    tcm = create_tcm(it, vectorizer, skip_grams_window = gramwin,
                     weights =  1 / seq_len(gramwin))
    
  }
  rm(it, vocab2, vectorizer);gc() # free some memory
  tcm = tcm + Matrix::t(triu(tcm))  # originally triangular, copy to make it work
  # NB: vocab_vectorizer and create_tcm are changed in text2vec 0.5.0, 
  # now params in latter; also has symmetric window and weights options!
 
  Matrix::diag(tcm) = 0
  
  # keep only persistent lexicon as context (and pos if set)
  if(!is.na(comlex[1])){
    if(nchar(wordclassesc) > 1){ comlex = grep(wordclassesc, comlex, value=T)} # if filtered comlex
    # NB: currently filtering wordclass in columns works only with comlex, tweak ifelse if needed
    if(nchar(wordclassesr)>1) {
      tcm = tcm[grep(wordclassesr,rownames(tcm)), comlex ]  }
    else { tcm = tcm[, comlex ] }
  } else { 
    # if NO COMLEX
    if(nchar(wordclassesr)>1) {
      tcm = tcm[grep(wordclassesr,rownames(tcm)),  ]  
    }   # NO else, use full tcm, no filtering for wordclass or comlex
    
  }  
  gc()
  if(ppmi){
    tcm = makeppmi(tcm)               # run ppmi weigthing
  }
  #print(paste(Sys.time(), "tcm done")) # in cluster, nowhere to print
  return(tcm)
}



makeppmi = function(pmat, positive=T){
  library(Matrix, quietly = T)
  pmat = Matrix::t(pmat)
  #pmat = (tcmfornews)
  #set.seed(1)
  #pmat = matrix(sample(c(0,0,0,0,0,0,1,10),5*10,T), 5,10, byrow=T)
  #pmat = Matrix::t(Matrix(pmat, sparse=T))
  
  tcmrs = Matrix::rowSums(pmat)
  tcmcs = Matrix::colSums(pmat)
  N = sum(tcmrs)
  colp = tcmcs/N
  rowp = tcmrs/N
  
  pp = pmat@p+1
  ip = pmat@i+1
  tmpx = rep(0,length(pmat@x))
  for(i in 1:(length(pmat@p)-1) ){ 
    #for(i in 1:100 ){ 
    #not0 = which(pmat[, i] > 0)
    ind = pp[i]:(pp[i+1]-1)
    not0 = ip[ind]
    icol = pmat@x[ind]
    #print(icol)
    #tmp = log( (pmat[not0,i]/N) / (rowp[not0] * colp[i] ))
    tmp = log2( (icol/N) / (rowp[not0] * colp[i] ))
    tmpx[ind] = tmp
    #print(tmp)
    # tmp = ifelse(tmp < 0, 0, tmp)
    #pmat2[not0, i] = tmp
  }
  if(positive){
    tmpx[tmpx<0] = 0 # see jätab nullid sisse, aga samas ei tõsta suurust ka
  }
  pmat2 = pmat
  pmat2@x = tmpx
  pmat2 = Matrix::t(pmat2) # tagasi sõnad ridadeks
  pmat2 = drop0(pmat2,is.Csparse=T) # nullid: võtab oma 100 mega maha, väga hea
  return(pmat2)
}


aggregaterelevants = function(foldr){
  # aggregate relevants, 20-30mb per period
  allrelevants = list()
  for(x in list.files(foldr, full.names = T)){
    load(x)
    allrelevants[[gsub("[^0-9]","",x) ]] = relevants 
  }
  gc()
  return(allrelevants)
}


dofreqdif = function(countmat, smooth, lexicons){
  #lexicons = getlexicons(countmat[[1]],countmat[[2]], minraw) # 1 comlex for columns, 2 lexicon for dataframes
  Sys.time()
  freqdifmat = matrix(nrow=length(lexicons), ncol=ncol(countmat)); rownames(freqdifmat)=lexicons
  for(i in 2:ncol(countmat)){
    fr1 = countmat[lexicons,i-1] +smooth
    fr2 = countmat[lexicons,i]   +smooth
    #freqdif1 = (abs(fr1 - fr2))/fr1*100; freqdif1 = ifelse(fr1 > fr2, freqdif1*-1, freqdif1)
    freqdifmat[,i] = log(fr2/fr1)    #freqdif1
  }
  gc()
  return(freqdifmat)
}
library(compiler)
dofreqdif = cmpfun(dofreqdif, options=list(optimize=3))

dotopicsloop = function(periods, difmat, countmat=NULL, wordclassesr, threshold,
                        estimate0=F, estimate2=F){ 
  # wt = list(1:relevancy_threshold, 
  #           log(rev(1:relevancy_threshold)), 
  #           rep(1,relevancy_threshold), 
  #           (-log(1:relevancy_threshold)+5),  
  #           (-log(1:relevancy_threshold)+5)^2, 
  #           (-log(1:relevancy_threshold)+5)^5) [[wtype]]
  if(length(wordclassesr) == 1){
    words = grep(wordclassesr, rownames(difmat), value = T)  # by regex; eg do only noun words
  } else {
    words = wordclassesr
  } # if len>1 then assume wordlist, use that instead
  
  inertiamat=matrix(NA, nrow=length(words),ncol=ncol(difmat)); rownames(inertiamat) = words
  #inertiamat2 = inertiamat
  estimat = inertiamat
  for(x in 2:length(periods)){
    #print(x)
    relevants=periods[[x]]
    inperiod = which(words %in% names(relevants))
    
    # determine if estimating based on previous or future relevant vectors
    # disable if using interpolation at the TCM level
    canestimate0=NULL
    if(estimate0){
      canestimate0 = setdiff(which(words %in% names(periods[[x-1]])), inperiod)
    }
    canestimate2=NULL
    if(estimate2){
      if(x < length(periods)){  # if not last period
        canestimate2 = setdiff(which(words %in% names(periods[[x+1]])),
                               union(canestimate0, inperiod))
      }
    }
    #estimat[inperiod,x] = F; estimat[union(canestimate0, canestimate2), x] = T # others NA
    
    # DO ADVECTION
    # advection from current values of present relevants
    for(i in inperiod){ 
      nrelevants = 1:(min(threshold, length(relevants[[words[i] ]])))
      inertiamat[i,x] = weighted.mean(difmat   [ names(relevants[[words[i] ]][nrelevants] ), x],
                                      w = relevants[[ words[i] ]][nrelevants] , na.rm = T) 
      #inertiamat2[i,x] = weighted.mean(countmat[ names(relevants[[words[i] ]] ), x],
      #                                w = relevants[[ words[i] ]], na.rm = T) 
    }
    
    # interpolate if missing:
    # advection from current values of relevants from t-1:
    if(length(canestimate0)>0){
      relevants=periods[[x-1]]
      for(i in canestimate0){ 
        nrelevants = 1:(min(threshold, length(relevants[[words[i] ]])))
        inertiamat[i,x] = weighted.mean(difmat   [ names(relevants[[words[i] ]][nrelevants] ), x],
                                        w = relevants[[ words[i] ]][nrelevants] , na.rm = T)
        # inertiamat2[i,x] = weighted.mean(countmat[ names(relevants[[words[i] ]]), x],
        #                                 w = relevants[[ words[i] ]], na.rm = T) 
      }
    }
    # advection from current values of relevants from t+1:
    if(length(canestimate2)>0){ # does not go into if last period (length null)
      relevants=periods[[x+1]]
      for(i in canestimate2){
        nrelevants = 1:(min(threshold, length(relevants[[words[i] ]])))
        inertiamat[i,x] = weighted.mean(difmat[ names(relevants[[words[i] ]][nrelevants] ), x],
                                        w = relevants[[ words[i] ]][nrelevants], na.rm = T)
        # inertiamat2[i,x] = weighted.mean(countmat[ names(relevants[[words[i] ]]), x],
        #                              w = relevants[[ words[i] ]], na.rm = T) 
      }
    }
    rm(relevants);gc()
  }
  #attr(inertiamat, "estimat") = estimat
  return( list( inertiamat=inertiamat, estimat=estimat))
}
library(compiler)
dotopicsloop = cmpfun(dotopicsloop, options=list(optimize=3))

dorelevants = function(ppmimat,relevancy_threshold){
  library(grr)
  relevants = list()
  print(paste(Sys.time(), "do relevance vectors"))

  relevants = list(); length(relevants) = nrow(ppmimat)
  ppmimat = Matrix::as.matrix(ppmimat)
  for(x in 1:nrow(ppmimat)){
    tmp = ppmimat[x,-which(ppmimat[x,]==0)]
    y=tmp[rev(order2(tmp))] # sort2 is decreasing=F BUT MESSES UP NAMES, USE order2
    y=y[1:min(length(y),relevancy_threshold)] # but this need top ones, so rev
    relevants[[x]] = y; names(relevants[[x]]) = names(y)
  }
  print(paste(Sys.time(),"relevance vectors done"))
  return(relevants)
}
dorelevants = cmpfun(dorelevants, options=list(optimize=3))



tartuplots  = function(a, f, isnew, isgone, titl, ann, s=5, rs=F, alph=0.3, rsline=T){
  if(is.null(isnew)){
    cols=rgb(0,0,0,alph)
  } else {
    cols=ifelse(isnew, rgb(0.7,0,0,0.7), ifelse(isgone, rgb(0,0,0.7,0.7), rgb(0,0,0,alph)))
  }
  ok = !is.na(a)
  a=a[ok]
  f=f[ok]
  fit=fitted(lm(a~f))
  sx = SharedData$new(data.frame(a,f), group=titl)
  p =plot_ly(type="scatter", mode="markers", hoverinfo = 'text', text=~gsub("S:","",sx$key()), 
             hoverlabel = list(bgcolor="white",  font=list(size=20, color=cols)) ) %>% 
    add_trace(x = ~f, y = ~a, data=sx, marker=list(color= cols, size=s )) %>% 
    add_lines(x = ~f, y = ~fit, line = list(color = rgb(0,0,0)) ,hoverinfo="none", visible=rsline) %>% 
    layout(showlegend = FALSE,autosize=T,
           #title = titl, 
           #titlefont=list(size=13),
           paper_bgcolor=rgb(0,0,0,0.04), plot_bgcolor=rgb(0,0,0,0),
           annotations = list(text = ann,  x = max(f,na.rm=T)*0.7, y = min(a, na.rm=T)*0.95,showarrow=FALSE ),
           margin = list(b = 100) #, hovermode = 'compare'
    )
  
  if(rs){
    p = add_segments(p, x = f, y = 0, xend = f, yend = a, line = list(color = rgb(0,0,0,0.2),width = 0.5), hoverlabel = list(bgcolor="white",  font=list(size=20, color="black")))
    p = layout(p, xaxis = list(title="fitted (frequency change ~ advection)", color="gray", gridcolor="white", zerolinecolor="white",zerolinewidth=4), 
               yaxis = list(title="residuals", color="gray", gridcolor="white",
                            zerolinecolor="black",zerolinewidth=2))
  } else {
    p = layout(p, 
               xaxis = list(title="log frequency change", color="gray", gridcolor="white", 
                            zerolinecolor="white",zerolinewidth=4), 
               yaxis = list(title="advection (log context change)", color="gray", gridcolor="white",
                            zerolinecolor="white",zerolinewidth=4) )
  }
  
  
  p = highlight(p, on = "plotly_click", off="plotly_doubleclick", persistent = F, selectize = TRUE)
  #htmltools::div(p, align="center", style="{height: 10%; width: 10%}") 
  p$sizingPolicy$defaultHeight = 1000
  p$sizingPolicy$browser$defaultHeight = 1000
  p$sizingPolicy$knitr$defaultHeight = 1000
  #p$sizingPolicy$browser$defaultHeight = "20%"
  p
  
}
