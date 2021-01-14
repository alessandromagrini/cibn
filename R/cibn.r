require(graph)
require(gRain)


# cut code at a symbol (internal use only)
strcut <- function(x,sym,keep=F) {
  for(i in 1:length(x)) {
    if(grepl(sym,x[i])) {
      iaux <- strsplit(x[i],sym)[[1]]
      if(keep==T) {
        for(j in 1:(length(iaux)-1)) {
          iaux[j] <- paste(iaux[j],gsub("\\\\","",sym),sep="")
          }
        }
      if(i==1) {
        x <- c(iaux,x[(i+1):length(x)])
        } else if(i<length(x)) {
        x <- c(x[1:(i-1)],iaux,x[(i+1):length(x)])
        } else {
        x <- c(x[1:(i-1)],iaux)
        }
      }
    }
  x
  }

# get string between two symbols (internal use only)
gbetw <- function(x,sym1,sym2) {
  if(grepl(sym1,x) & grepl(sym2,x)) {
    strsplit(strsplit(x,sym1)[[1]][-1],sym2)[[1]][1]
    } else {
    x
    }
  }
  
# scan code sections (internal use only)
scanbl <- function(x) {
  x <- strcut(x,"\\{",keep=T)
  x <- strcut(x,"\\}",keep=T)
  for(i in 1:length(x)) {
    grepl("\\{",x[i])
    }
  aux1 <- grep("\\{",x)
  aux2 <- grep("\\}",x)
  for(i in 1:length(aux1)) {
    if(length(aux2)<i) {
      stop("Unexpected left brace in model code at line ",aux1[i],call.=F)
      } else if(aux1[i]>aux2[i]) {
      stop("Unexpected right brace in model code at line ",aux2[i],call.=F)
      } else {
      if(i==length(aux1)&length(aux2)>i) stop("Unexpected right brace in model code at line ",aux2[i],call.=F)
      }
    }
  cbind(aux1,aux2)
  }

# trim a string (internal use only)
trimstr <- function(x) {
  x <- gsub("^\\s+|\\s+$","",x)
  x <- gsub("\\( {1,}","\\(",x)
  x <- gsub(" {1,}\\)","\\)",x)
  x <- gsub(" {1,},",",",x)
  x <- gsub(", {1,}",",",x)
  x <- gsub(" {1,}=","=",x)
  x <- gsub("= {1,}","=",x)
  x <- gsub("\\(,","\\(",x)
  x <- gsub(",\\)","\\)",x)
  x
  }

# clean code section (internal use only)
cleanbl <- function(x) {                                                              
  #auxcod <- strsplit(strsplit(paste(x,collapse=";"),"\\{")[[1]][-1],"\\}")[[1]][1]
  auxcod <- gbetw(paste(x,collapse=";"),"\\{","\\}")
  auxcod <- trimstr(strsplit(auxcod,";")[[1]])                                       
  auxcod[which(sapply(auxcod,nchar)>0)]
  }

# split an argument (internal use only)
argsplit <- function(x) {
  res <- strsplit(x,"=")[[1]]
  if(substr(x,nchar(x),nchar(x))=="=") {
    c(res,"")
    } else {
    res
    }
  }

# check node name (internal use only)
#   (must begin with a capital letter, no special characters excepting '_')
checkNodeName <- function(x) {
  n <- nchar(x)
  ifelse(n<1 || n>20 || (substr(x,1,1) %in% toupper(letters))==F ||
    grepl("^LAMBDA",toupper(x)) || grepl("^AUX",toupper(x)), ok <- F, ok <- T)
  if(ok==T & n>1) {
    for(i in 2:n) {
      if(length(intersect(substr(x,i,i),c(0:9,letters,toupper(letters),"_")))==0) {
        ok <- F
        break()
        }
      }
    }  
  ok
  }

# check state name (internal use only)
#   (no special characters excepting '_', cannot begin with 'LAMBDA' or 'AUX')
checkStateName <- function(x) {
  n <- nchar(x)
  ifelse(n<1 || n>50 || #(substr(x,1,1) %in% c(toupper(letters),letters))==F ||
    grepl("^LAMBDA",toupper(x)) || grepl("^AUX",toupper(x)), ok <- F, ok <- T)
  if(ok==T & n>1) {
    for(i in 2:n) {
      if(length(intersect(substr(x,i,i),c(0:9,letters,toupper(letters),"_")))==0) {
        ok <- F
        break()
        }
      }
    }  
  ok
  }
 
# find a topological order for a graph (internal use only)
topolog <- function(parSet) {
  nomi <- names(parSet)
  L <- c()
  S <- nomi[which(sapply(parSet,length)==0)]
  while(length(S)>0) {
    xaux <- S[1]
    S <- setdiff(S,xaux)
    L <- c(L,xaux)
    sch <- c()
    for(j in 1:length(parSet)) {
      if(xaux %in% parSet[[j]]) sch <- c(sch,nomi[j])
      }
    if(length(sch)>0) {
      for(j in 1:length(sch)) {
        parSet[[sch[j]]] <- setdiff(parSet[[sch[j]]],xaux)
        if(length(parSet[[sch[j]]])==0) S <- c(S,sch[j])  
        }
      }
    }
  if(sum(sapply(parSet,length))==0) L else NULL
  }

# create a new network
new.cibn <- function(model.code=NULL,path=NULL,maximal=TRUE) {
  if(maximal) nmax <- 2 else nmax <- Inf
  if(!is.null(model.code)) {
    model.code <- strsplit(model.code,"\n")[[1]]
    } else {
    if(!is.null(path)) {
      model.code <- readLines(path,warn=F)
      } else {
      stop("One among arguments 'model.code' and 'path' must be provided",call.=F)
      }
    }
  if(length(model.code)==0) stop("The model code is empty",call.=F)
  for(i in 1:length(model.code)) {
    if(grepl("#",model.code[i])) model.code[i] <- strsplit(model.code[i],"#")[[1]][1]
    model.code[i] <- gsub("^\\s+|\\s+$","",model.code[i])
    }
  blockMat <- scanbl(model.code)
  if(nrow(blockMat)==0) stop("The model code contains no command",call.=F)  
  ### process commands 'variable' and 'interaction'
  varNames <- varType <- auxline <- intNames <- c()
  omega <- parSets <- descr <- list()
  for(i in 1:nrow(blockMat)) {     
    iblock <- model.code[blockMat[i,1]:blockMat[i,2]]  
    auxnam <- strsplit(iblock[1]," ")[[1]]
    auxcomm <- auxnam[1]
    if(auxcomm=="variable") {
      if(length(auxnam)==3) {
        auxvar <- auxnam[2]
        if(checkNodeName(auxvar)==F) stop("Invalid variable name '",auxvar,"' in command 'variable' at line ",blockMat[i,1],call.=F)
        if(auxvar %in% varNames) stop("Duplicated command 'variable' at line ",blockMat[i,1],call.=F)
        if(auxvar %in% intNames) stop("Non-unique name in command 'variable' at line ",blockMat[i,1],call.=F)
        varNames <- c(varNames,auxvar)
        } else {
        stop("Malformed command 'variable' at line ",blockMat[i,1],call.=F)
        }
      auxfeat <- cleanbl(iblock)
      ijstr <- lapply(auxfeat,argsplit)
      auxopt <- unname(sapply(ijstr,function(x){x[1]}))                                 
      auxdupl <- auxopt[duplicated(auxopt)]
      if(length(auxdupl)>0) stop("Duplicated argument '",auxdupl[1],"' in command 'variable' at line ",blockMat[i,1],call.=F)
      auxchk1 <- setdiff(c("type","states","parents"),auxopt)
      if(length(auxchk1)>0) stop("Missing argument '",auxchk1[1],"' in command 'variable' at line ",blockMat[i,1],call.=F)
      auxchk2 <- setdiff(auxopt,c("type","states","parents","description"))
      if(length(auxchk2)>0) stop("Unknown argument '",auxchk2[1],"' in command 'variable' at line ",blockMat[i,1],call.=F)
      auxline[auxvar] <- blockMat[i,1]
      auxij.typ <- which(auxopt=="type")
      if(length(ijstr[[auxij.typ]])!=2) stop("Malformed argument 'type' in command 'variable' at line ",blockMat[i,1],call.=F) 
      if((ijstr[[auxij.typ]][2] %in% c("GRAD","DGRAD","NOM")==F)) stop("Invalid value for argument 'type' in command 'variable' at line ",blockMat[i,1],call.=F) ##### 
      varType[auxvar] <- ijstr[[auxij.typ]][2]
      auxij.lev <- ijstr[[which(auxopt=="states")]]
      if(length(auxij.lev)!=2 || substr(auxij.lev[2],1,1)!="(" || substr(auxij.lev[2],nchar(auxij.lev[2]),nchar(auxij.lev[2]))!=")") stop("Malformed argument 'states' in command 'variable' at line ",blockMat[i,1],call.=F) 
      auxst <- strsplit(gbetw(auxij.lev[2],"\\(","\\)"),",")[[1]]
      auxstdup <- auxst[duplicated(auxst)]
      auxchk3 <- which(sapply(auxst, checkStateName)==F)
      if(length(auxchk3)>0) stop("Invalid state name '",auxst[auxchk3[1]],"' in command 'variable' at line ",blockMat[i,1],call.=F)
      if(length(auxstdup)>0) stop("Duplicated state '",auxstdup[1],"' in command 'variable' at line ",blockMat[i,1],call.=F)
      if(varType[auxvar]=="DGRAD" & (length(auxst) %% 2)==0) stop("Invalid number of states in command 'variable' at line ",blockMat[i,1],call.=F) 
      omega[[auxvar]] <- auxst                       
      auxij.par <- ijstr[[which(auxopt=="parents")]]
      if(length(auxij.par)!=2 || substr(auxij.par[2],1,1)!="(" || substr(auxij.par[2],nchar(auxij.par[2]),nchar(auxij.par[2]))!=")") stop("Malformed argument 'parents' in command 'variable' at line ",blockMat[i,1],call.=F) 
      auxpar <- strsplit(gbetw(auxij.par[2],"\\(","\\)"),",")[[1]] 
      parSets[[auxvar]] <- auxpar
      auxij.descr <- which(auxopt=="description")
      if(length(auxij.descr)>0) {
        idescr <- ijstr[[auxij.descr]][2]
        if(substr(idescr,1,1)!='<'||substr(idescr,nchar(idescr),nchar(idescr))!='>') stop("Malformed argument 'description' in command 'variable' at line ",blockMat[i,1],call.=F)
        descr[auxvar] <- substr(idescr,2,nchar(idescr)-1)
        } else {
        descr[auxvar] <- ""
        }
      } else if(auxcomm=="interaction") {
      if(length(auxnam)==3) {
        auxvar <- auxnam[2]
        if(checkNodeName(auxvar)==F) stop("Invalid variable name in command 'interaction' at line ",blockMat[i,1],call.=F)
        if(auxvar %in% intNames) stop("Duplicated command 'interaction' at line ",blockMat[i,1],call.=F)
        if(auxvar %in% varNames) stop("Non-unique name in command 'interaction' at line ",blockMat[i,1],call.=F)
        } else {
        stop("Malformed command 'interaction' at line ",blockMat[i,1],call.=F)
        }  
      auxfeat <- cleanbl(iblock)
      if(length(auxfeat)==0) stop("Empty command 'function' at line ",blockMat[i,1],call.=F) 
      ifeamat <- list()                                
      for(j in 1:length(auxfeat)) {                    
        ifeamat[[j]] <- argsplit(auxfeat[j])
        }
      auxline[auxvar] <- blockMat[i,1]
      auxarg <- sapply(ifeamat,function(x){x[1]})
      if(("from" %in% auxarg)==F) stop("Missing argument 'from' in command 'interaction' at line ",blockMat[i,1],call.=F)  
      auxcheck <- auxarg[which(auxarg %in% c("from","description")==F)]
      if(length(auxcheck)>0) stop("Unknown argument '",auxcheck[1],"' in command 'interaction' at line ",blockMat[i,1],call.=F)
      auxdupl <- auxarg[duplicated(auxarg)]
      if(length(auxdupl)>0) stop("Duplicated argument '",auxdupl[1],"' in command 'interaction' at line ",blockMat[i,1],call.=F)      
      iauxfr <- ifeamat[[which(auxarg=="from")]]
      if(length(iauxfr)!=2) stop("Malformed argument 'from' in command 'interaction' at line ",blockMat[i,1],call.=F)
      iauxpar <- gsub("\\(","",gsub("\\)","",strsplit(iauxfr[-1],",")[[1]]))
      ifrdup <- iauxpar[duplicated(iauxpar)]
      if(length(ifrdup)>0) stop("Duplicated variable '",ifrdup[1],"' provided to argument 'from' in command 'interaction' at line ",blockMat[i,1],call.=F)
      auxij.descr <- which(auxarg=="description")
      if(length(auxij.descr)>0) {
        idescr <- ifeamat[[which(auxarg=="description")]][-1]
        if(substr(idescr,1,1)!='<'||substr(idescr,nchar(idescr),nchar(idescr))!='>') stop("Malformed argument 'description' in command 'interaction' at line ",blockMat[i,1],call.=F)
        descr[auxvar] <- substr(idescr,2,nchar(idescr)-1)
        } else {
        descr[auxvar] <- ""
        }
      intNames <- c(intNames,auxvar)
      parSets[[auxvar]] <- iauxpar
      } else {
      if(auxcomm!="model") stop("Unknown command '",auxcomm,"' at line ",blockMat[i,1],call.=F)
      }
    }
  ### check parent sets
  for(i in 1:length(varNames)) {
    ipset <- parSets[[varNames[i]]]
    auxchk <- which((ipset %in% c(varNames,intNames))==F)
    if(length(auxchk)>0) stop("Unknown variable '",ipset[[auxchk[1]]],"' provided to argument 'parents' in command 'variable' at line ",auxline[varNames[i]],call.=F)
    auxdupl <- ipset[duplicated(ipset)]
    if(length(auxdupl)>0) stop("Duplicated variable '",auxdupl[1],"' provided to argument 'parents' in command 'variable' at line ",blockMat[i,1],call.=F)      
    }
  if(length(intNames)>0) {
    for(i in 1:length(intNames)) {
      ipset <- parSets[[intNames[i]]]
      auxchk1 <- which((ipset %in% varNames)==F)
      if(length(auxchk1)>0) stop("Unknown variable '",ipset[[auxchk1[1]]],"' provided to argument 'from' in command 'interaction' at line ",auxline[intNames[i]],call.=F)
      omega[[intNames[i]]] <- apply(cartesP(omega[ipset]),1,paste,collapse="+")
      varType[[intNames[i]]] <- "INTER"
      }
    }
  if(is.null(topolog(parSets))) stop("The model contains directed cycles",call.=F)
  ### process commands 'model'
  thetaList <- list()
  for(i in 1:nrow(blockMat)) {     
    iblock <- model.code[blockMat[i,1]:blockMat[i,2]]  
    auxnam <- strsplit(iblock[1]," ")[[1]]
    auxcomm <- auxnam[1]
    if(auxcomm=="model") {
      if(length(auxnam)==3) {
        auxvar <- auxnam[2]
        if((auxvar %in% varNames)==F) stop("Unknown variable '",auxvar,"' in command 'model' at line ",blockMat[i,1],call.=F)
        if(auxvar %in% names(thetaList)) stop("Duplicated command 'model' at line ",blockMat[i,1],call.=F)
        } else {
        stop("Malformed command 'model' at line ",blockMat[i,1],call.=F)
        }
      auxfeat <- cleanbl(iblock)
      ifeamat <- c()
      if(length(auxfeat)==0) stop("Empty command 'model' at line ",blockMat[i,1],call.=F) 
      for(j in 1:length(auxfeat)) {
        jthstr0 <- argsplit(auxfeat[j])
        if(length(jthstr0)!=2) stop("Malformed command 'model' at line ",blockMat[i,1],call.=F) 
        jthstr1 <- strsplit(jthstr0[1],":")[[1]]
        #jthstr2 <- strsplit(jthstr0[2],"\\*")[[1]]
        jthstr2 <- jthstr0[2]
        if(jthstr1[1]!="omitted") {
          if((jthstr1[1] %in% parSets[[auxvar]])==F) stop("Unknown parent '",jthstr1[1],"' in command 'model' at line ",blockMat[i,1],call.=F)
          if(is.na(jthstr1[2])) stop("Missing state for parent '",jthstr1[1],"' in command 'model' at line ",blockMat[i,1],call.=F)
          if((jthstr1[2] %in% omega[[jthstr1[1]]])==F) stop("Unknown state for parent '",jthstr1[1],"' in command 'model' at line ",blockMat[i,1],call.=F) 
          if(varType[[jthstr1[1]]]=="DGRAD") {
            jthneu <- omega[[jthstr1[1]]][(1+length(omega[[jthstr1[1]]]))/2]
            } else {
            jthneu <- omega[[jthstr1[1]]][1]            
            }
          if(identical(jthstr1[2],jthneu)) stop("The reference state for parent '",jthstr1[1],"' appears in command 'model' at line ",blockMat[i,1],call.=F)
          }
        if(substr(jthstr2,1,1)!="("||substr(jthstr2,nchar(jthstr2),nchar(jthstr2))!=")") stop("Malformed parameter value for '",jthstr1[1],"' in command 'model' at line ",blockMat[i,1],call.=F)
        #
        #jthval <- gbetw(jthstr2[1],"\\(","\\)")
        #options(warn=-1)
        #auxcheck <- as.numeric(strsplit(jthval,",")[[1]])
        #
        jthval <- paste("c(",gbetw(jthstr2[1],"\\(","\\)"),")",sep="")
        auxcheck <- try(eval(parse(text=jthval)),silent=T)
        #
        if(!identical(class(auxcheck),"numeric")) stop("Non-numerical parameter value for '",jthstr1[1],"' in command 'model' at line ",blockMat[i,1],call.=F)
        #
        ### <--------------------------------- GESTIRE VALORI NON NUMERICI
        #
        if(length(auxcheck)!=length(omega[[auxvar]])) stop("Uncorrect parameter length for '",jthstr1[1],"' in command 'model' at line ",blockMat[i,1],call.=F)
        if(sum(is.na(auxcheck))>0) stop("Non-numeric parameter value for '",jthstr1[1],"' in command 'model' at line ",blockMat[i,1],call.=F)
        if(sum(is.na(auxcheck))>0 || sum(auxcheck<0)>0) stop("Invalid parameter value for '",jthstr1[1],"' in command 'model' at line ",blockMat[i,1],call.=F)
        options(warn=0)
        ifeamat <- rbind(ifeamat,c(jthstr1[1],jthstr1[2],jthval))
        }
      colnames(ifeamat) <- c("Parent","Parent_state","Param")
      rownames(ifeamat) <- NULL
      auxcheck <- setdiff(parSets[[auxvar]],ifeamat[,"Parent"])
      if(length(auxcheck)>0) stop("Missing parameter for parent '",auxcheck[1],"' in command 'model' at line ",blockMat[i,1],call.=F)
      auxdupl <- ifeamat[duplicated(ifeamat[,"Parent"]),"Parent"]
      if(length(auxdupl)>0) {
        for(j in 1:length(auxdupl)) {
          auxjdu <- duplicated(ifeamat[which(ifeamat[,"Parent"]==auxdupl[j]),"Parent_state"])
          if(sum(auxjdu)>0) stop("Duplicated parameter for parent '",auxdupl[j],"' in command 'model' at line ",blockMat[i,1],call.=F)
          }
        }
      auxuni <- unique(setdiff(ifeamat[,"Parent"],"omitted"))
      if(length(auxuni)>0) {
        for(j in 1:length(auxuni)) {
          if(identical(sort(omega[[auxuni[j]]]),sort(ifeamat[which(ifeamat[,"Parent"]==auxuni[j]),"Parent_state"]))) {
            stop("Too much parameters for parent '",auxuni[j],"' in command 'model' at line ",blockMat[i,1],call.=F)
            }
          }
        }
      thetaList[[auxvar]] <- ifeamat     
      }
    }
  modchk <- setdiff(varNames,names(thetaList))
  if(length(modchk)>0) stop("Missing command 'model' for variable'",modchk[1],"'",call.=F)
  params <- list()
  nomi <- c()
  for(i in 1:length(varNames)) {
    inam <- varNames[i]
    ipar <- parSets[[inam]]
    ist <- omega[[inam]]
    ithet <- thetaList[[inam]]
    idis0 <- rep(0,length(ist))
    if(varType[[inam]]=="DGRAD") {
      idis0[(1+length(ist))/2] <- 1
      } else {
      idis0[1] <- 1
      }
    if("omitted" %in% ithet[,1]) {
      ip0 <- eval(parse(text=ithet[which(ithet[,1]=="omitted"),3]))
      } else {
      ip0 <- idis0
      }
    if(varType[[inam]]=="GRAD") {           
      ### create the CPT for a GRAD node
      igamnam <- paste(paste("Lambda1.",inam,"..",sep=""),c("omitted",ipar),sep="")
      if("omitted" %in% ithet[,1]) {
        ip0 <- eval(parse(text=ithet[which(ithet[,1]=="omitted"),3]))
        } else {
        ip0 <- idis0
        }
      itab <- array(ip0,dim=length(ist))
      idimnam <- list(ist)
      names(idimnam) <- igamnam[1]
      dimnames(itab) <- idimnam
      class(itab) <- "table"
      params <- c(params,list(itab))
      nomi <- c(nomi,igamnam[1])
      if(length(ipar)>0) {
        for(w in 1:length(ipar)) {
          iwst <- omega[[ipar[w]]]
          itab <- array(dim=c(length(ist),length(iwst)))
          idimnam <- list(ist,iwst)
          names(idimnam) <- c(igamnam[w+1],ipar[w])
          dimnames(itab) <- idimnam
          class(itab) <- "table"
          for(j in 1:length(iwst)) {
            if(ipar[w] %in% ithet[,1]) {
              if(iwst[j] %in% ithet[which(ithet[,1]==ipar[w]),2]) {
                ipj <- eval(parse(text=ithet[which(ithet[,1]==ipar[w]&ithet[,2]==iwst[j]),3]))
                } else {
                ipj <- idis0
                }
              } else {
              ipj <- idis0
              }
            itab[,j] <- ipj
            }
          params <- c(params,list(itab))
          nomi <- c(nomi,igamnam[w+1])
          }
        }
      #if(is.null(max.parents)) {
      #  nmax <- length(ipar)+1
      #  } else {
      #  nmax <- max.parents
      #  }
      ifun <- decompCpt(inam,ist,ipar,id=1,start=1,nmax=nmax,type="grad")
      names(ifun)[length(ifun)] <- inam
      params <- c(params,ifun)
      nomi <- c(nomi,names(ifun))
      } else if(varType[[inam]]=="DGRAD") {
      ### create the CPT for a DGRAD node
      iaux1 <- 1:((1+length(ist))/2)
      iaux2 <- ((1+length(ist))/2):length(ist)
      ist1 <- ist[iaux1]
      ist2 <- ist[iaux2]
      idis0_1 <- idis0[iaux1]
      idis0_2 <- idis0[iaux2]
      igamnam1 <- paste(paste("Lambda1.",inam,"..",sep=""),c("omitted",ipar),sep="")
      igamnam2 <- paste(paste("Lambda2.",inam,"..",sep=""),c("omitted",ipar),sep="")
      if("omitted" %in% ithet[,1]) {
        ip0 <- eval(parse(text=ithet[which(ithet[,1]=="omitted"),3]))
        ip0_1 <- ip0[iaux1]
        ip0_2 <- ip0[iaux2]
        } else {
        ip0_1 <- idis0_1
        ip0_2 <- idis0_2
        }
      itab1 <- array(rev(ip0_1),dim=length(ist1))
      idimnam1 <- list(rev(ist1))
      names(idimnam1) <- igamnam1[1]
      dimnames(itab1) <- idimnam1
      class(itab1) <- "table"
      itab2 <- array(ip0_2,dim=length(ist2))
      idimnam2 <- list(ist2)
      names(idimnam2) <- igamnam2[1]
      dimnames(itab2) <- idimnam2
      class(itab2) <- "table"
      params <- c(params,list(itab1,itab2))
      nomi <- c(nomi,igamnam1[1],igamnam2[1])
      if(length(ipar)>0) {
        for(w in 1:length(ipar)) {
          iwst <- omega[[ipar[w]]]
          itab1 <- array(dim=c(length(ist1),length(iwst)))
          idimnam1 <- list(rev(ist1),iwst)
          names(idimnam1) <- c(igamnam1[w+1],ipar[w])
          dimnames(itab1) <- idimnam1
          class(itab1) <- "table"
          itab2 <- array(dim=c(length(ist2),length(iwst)))
          idimnam2 <- list(ist2,iwst)
          names(idimnam2) <- c(igamnam2[w+1],ipar[w])
          dimnames(itab2) <- idimnam2
          class(itab2) <- "table"
          for(j in 1:length(iwst)) {
            if(ipar[w] %in% ithet[,1]) {
              if(iwst[j] %in% ithet[which(ithet[,1]==ipar[w]),2]) {
                ipj <- eval(parse(text=ithet[which(ithet[,1]==ipar[w]&ithet[,2]==iwst[j]),3]))
                } else {
                ipj <- idis0
                }
              } else {
              ipj <- idis0
              }
            itab1[,j] <- rev(ipj[iaux1])
            itab2[,j] <- ipj[iaux2]
            }
          params <- c(params,list(itab1,itab2))
          nomi <- c(nomi,igamnam1[w+1],igamnam2[w+1])
          }
        }
      #if(is.null(max.parents)) {
      #  nmax <- length(ipar)+1
      #  } else {
      #  nmax <- max.parents
      #  }
      ifun1 <- decompCpt(inam,rev(ist[iaux1]),ipar,id=1,start=1,nmax=nmax,type="dgrad")
      end1 <- max(getID(names(ifun1)))
      ifun2 <- decompCpt(inam,ist[iaux2],ipar,id=2,start=end1+1,nmax=nmax,type="dgrad")          
      end2 <- max(getID(names(ifun2)))
      params <- c(params,ifun1,ifun2,list(dgradCpt(inam,ist,end1,end2)))
      nomi <- c(nomi,names(ifun1),names(ifun2),inam)
      } else {
      ### create the CPT for a NOM node
      if("omitted" %in% ithet[,1]) {
        ip0 <- eval(parse(text=ithet[which(ithet[,1]=="omitted"),3]))
        } else {
        ip0 <- idis0
        }
      igamnam <- list()
      for(z in 2:length(ist)) {
        igamnam[[z-1]] <- paste(paste("Lambda",z-1,".",inam,"..",sep=""),c("omitted",ipar),sep="")
        }
      itab <- list()
      for(z in 2:length(ist)) {
        itab[[z-1]] <- array(c(sum(ip0)-ip0[z],ip0[z]),dim=2)
        ijdnam <- list(c(0,1))
        names(ijdnam) <- igamnam[[z-1]][1]
        dimnames(itab[[z-1]]) <- ijdnam
        class(itab[[z-1]]) <- "table"
        }
      params <- c(params,itab)
      nomi <- c(nomi,sapply(igamnam,function(v){v[1]}))
      if(length(ipar)>0) {
        for(w in 1:length(ipar)) {
          iwst <- omega[[ipar[w]]]
          itab <- list()
          for(z in 2:length(ist)) {
            itab[[z-1]] <- array(dim=c(2,length(iwst)))
            ijdnam <- list(c(0,1),iwst)
            names(ijdnam) <- c(igamnam[[z-1]][w+1],ipar[w])
            dimnames(itab[[z-1]]) <- ijdnam
            class(itab[[z-1]]) <- "table"
            }
          for(j in 1:length(iwst)) {
            if(ipar[w] %in% ithet[,1]) {
              if(iwst[j] %in% ithet[which(ithet[,1]==ipar[w]),2]) {
                ipj <- eval(parse(text=ithet[which(ithet[,1]==ipar[w]&ithet[,2]==iwst[j]),3]))
                } else {
                ipj <- idis0
                }
              } else {
              ipj <- idis0
              }
            for(z in 2:length(ist)) {
              itab[[z-1]][,j] <- c(sum(ipj)-ipj[z],ipj[z])
              }
            }
          params <- c(params,itab)
          nomi <- c(nomi,sapply(igamnam,function(z){z[w+1]}))
          }
        }
      #if(is.null(max.parents)) {
      #  nmax <- length(ipar)+1
      #  } else {
      #  nmax <- max.parents
      #  }
      ifun <- c()
      end0 <- 0
      for(j in 2:length(ist)) {
        ijf <- decompCpt(inam,c(0,1),ipar,id=j-1,start=end0[j-1]+1,nmax=nmax,type="nom")
        ifun <- c(ifun,ijf)
        end0 <- c(end0,max(getID(names(ijf))))
        nomi <- c(nomi,names(ijf))
        }
      params <- c(params,ifun,list(nomCpt(inam,ist,end0[-1])))
      nomi <- c(nomi,inam)
      }
    }
  if(length(intNames)>0) {
    for(i in 1:length(intNames)) {
      ipar <- parSets[[intNames[i]]]
      ipst <- omega[ipar]      
      ist <- apply(cartesP(ipst),1,paste,collapse="+")
      icpt <- array(0,dim=c(length(ist),sapply(ipst,length)))
      kstr <- paste(paste("k",1:length(ipar),sep=""),collapse=",")
      mystr <- "auxcou <- 0; "
      for(j in 1:length(ipar)) {
        mystr <- paste(mystr,"for(k",j," in 1:",length(omega[[ipar[j]]]),") {; ",sep="")
        }
      mystr <- paste(mystr,"auxcou <- auxcou+1; icpt[auxcou,",kstr,"] <- 1",sep="")  
      mystr <- paste(mystr,paste(rep("}",length(ipar)),collapse="; "),sep="")
      eval(parse(text=mystr))
      idnam <- c(list(ist),ipst)
      names(idnam) <- c(intNames[i],ipar)
      dimnames(icpt) <- idnam
      class(icpt) <- "table"
      params <- c(params,list(icpt))
      nomi <- c(nomi,intNames[i])
      }
    }
  names(params) <- nomi
  RES <- list(CPTs=normalise(params),prior=params,posterior=NULL,description=descr)
  class(RES) <- "cibn"
  RES
  }

# get id of auxiliary nodes (internal use only)
getID <- function(x) {
  res <- c()
  for(i in 1:length(x)) {
    if(grepl("^AUX[1-9]{1,}.",x[i])) {
      res[i] <- as.numeric(strsplit(strsplit(x[i],"AUX")[[1]][2],"\\.")[[1]][1])
      }
    }
  res
  }

# create a DGRAD CPT (internal use only)
dgradCpt <- function(yname,ystates,id1,id2) {
  yaux1 <- ((1+length(ystates))/2):1
  yaux2 <- ((1+length(ystates))/2):length(ystates)
  yval <- c(1:length(ystates))-yaux2[1]
  yval1 <- yval[yaux1]
  yval2 <- yval[yaux2]
  res <- array(dim=c(length(ystates),length(yaux1),length(yaux2)))
  for(i in 1:length(yaux1)) {
    for(j in 1:length(yaux2)) {    
      ijp <- rep(0,length(ystates))
      ijval <- yval1[i]+yval2[j]
      ijp[ijval+(1+length(ystates))/2] <- 1
      res[,i,j] <- ijp
      }
    }
  dimnam <- list(ystates,ystates[yaux1],ystates[yaux2])
  names(dimnam) <- c(yname,paste("AUX",c(id1,id2),".",yname,sep=""))
  dimnames(res) <- dimnam
  class(res) <- "table"
  res
  }

# create a NOM CPT (internal use only)
nomCpt <- function(yname,ystates,idvet) {
  m <- length(ystates)-1
  res <- array(0,dim=c(m+1,rep(2,m)))
  kstr <- paste(paste("k",1:m,sep=""),collapse=",")
  mystr <- ""
  for(i in 1:m) {
    mystr <- paste(mystr,"for(k",i," in 1:2) {; ",sep="")
    }
  mystr <- paste(mystr,"auxind <- which(c(",kstr,")==2); if(length(auxind)==1) { res[auxind+1,",kstr,"] <- 1 } else { res[1,",kstr,"] <- 1 }",sep="")
  mystr <- paste(mystr,paste(rep("}",m),collapse="; "),sep="")
  eval(parse(text=mystr))
  dimnam <- vector("list",length=m+1)
  dimnam[[1]] <- ystates
  for(i in 1:m) {
    dimnam[[i+1]] <- c(0,1)
    }
  names(dimnam) <- c(yname,paste("AUX",idvet,".",yname,sep=""))
  dimnames(res) <- dimnam
  class(res) <- "table"
  res
  }

# create a MAX CPT (internal use only)
maxCpt <- function(yname,ystates,xnames,fun="max") {
  m <- length(ystates)
  p <- length(xnames)
  res <- array(0,dim=rep(m,p+1))
  kstr <- paste(paste("k",1:p,sep=""),collapse=",")
  mystr <- ""
  for(i in 1:p) {
    mystr <- paste(mystr,"for(k",i," in 1:",m,") {; ",sep="")
    }
  mystr <- paste(mystr,"res[",fun,"(c(",kstr,")),",kstr,"] <- 1",sep="")  
  mystr <- paste(mystr,paste(rep("}",p),collapse="; "),sep="")
  eval(parse(text=mystr))
  dimnam <- vector("list",length=p+1)
  for(i in 1:length(dimnam)) {
    dimnam[[i]] <- ystates
    }  
  names(dimnam) <- c(yname,xnames)
  dimnames(res) <- dimnam
  class(res) <- "table"
  res
  }

# decompose a CPT (internal use only)
decompCpt <- function(yname,ystates,xnames,id,start,nmax,type) {
  myfun <- function(z,s0) {
    res <- list()
    k <- floor(length(z)/nmax)
    if(length(z)<=nmax) {
      auxadd <- list(z)
      names(auxadd) <- yname
      res <- c(res,auxadd)
      list(res,s0)
      } else {
      aux <- 1
      for(i in 1:k) {                        
        auxadd <- list(z[aux:(aux+nmax-1)])
        names(auxadd) <- paste("AUX",s0+i,".",yname,sep="")
        res <- c(res,auxadd)
        aux <- aux+nmax
        }
      zadd <- names(res)
      if((length(z) %% nmax)>0) zadd <- c(zadd,z[(k*nmax+1):length(z)])
      zadd <- list(zadd)
      names(zadd) <- yname      
      list(c(res,zadd),s0+k)
      }
    }
  if(length(xnames)<2 & type=="grad") {
    itab <- maxCpt(yname,ystates,paste("Lambda1.",yname,"..",c("omitted",xnames),sep=""))
    RES <- list(itab)
    names(RES) <- yname
    } else {
    parnam <- paste(paste("Lambda",id,".",yname,"..",sep=""),c("omitted",xnames),sep="")
    newcalc <- myfun(parnam,start-1)
    newdag <- newcalc[[1]]
    newstart <- newcalc[[2]]
    current <- newdag[[yname]]
    OUT <- newdag[setdiff(names(newdag),yname)]   
    while(length(current)>nmax) {
      auxnam <- grep("AUX",current)
      newcalc <- myfun(current,newstart)
      newdag <- newcalc[[1]]
      newstart <- newcalc[[2]]
      OUT <- c(OUT,newdag[setdiff(names(newdag),yname)])
      current <- newdag[[yname]]
      }
    OUT <- c(OUT,list(current))
    names(OUT)[length(OUT)] <- paste("AUX",newstart+1,".",yname,paste="",sep="")
    RES <- list()
    for(i in 1:length(OUT)) {
      RES <- c(RES,list(maxCpt(names(OUT)[i],ystates,OUT[[i]])))
      }
    names(RES) <- names(OUT)
    }
  RES
  }

# find children of a node (internal use only)
nodchil <- function(counts,nodeName) {
  res <- c()
  nomi <- names(counts)
  for(i in 1:length(counts)) {
    inam <- nomi[i]
    idnam <- setdiff(names(dimnames(counts[[i]])),inam)
    if(nodeName %in% idnam) res <- c(res,inam)
    }
  if(length(res)>0) {
    sort(res)
    } else {
    character(0)
    }
  }

# cartesian product (internal use only)
cartesP <- function(omegList) {
  M <- expand.grid(rev(omegList))
  res <- as.matrix(M[,ncol(M):1])
  colnames(res) <- colnames(M)
  res
  }

# compute CPTs from counts (internal use only)
normalise <- function(counts) {
  theta <- counts
  for(i in 1:length(theta)) {
    idim <- dim(counts[[i]])                         
    if(length(idim)>1) {
      auxcfg <- list()
      for(w in 2:length(idim)) {
        auxcfg[[w-1]] <- 1:idim[w]
        }
      icfg <- cartesP(auxcfg)
      for(j in 1:nrow(icfg)) {
        cfgstr <- paste("for(w in 1:idim[1]) {; if(sum(unlist(counts[[i]][,",paste(icfg[j,],collapse=","),"]))==0) {; ",
          "theta[[i]][,",paste(icfg[j,],collapse=","),"][w] <- 1/idim[1];",
          "} else { ; theta[[i]][,",paste(icfg[j,],collapse=","),"][w] <- counts[[i]][,",
          paste(icfg[j,],collapse=","),"][w]/sum(counts[[i]][,",paste(icfg[j,],collapse=","),"])}; }",sep="")    
        eval(parse(text=cfgstr))
        }
      } else {
      for(j in 1:idim) {           
        if(sum(unlist(counts[[i]]))==0) {
          theta[[i]][j] <- 1/idim 
          } else {
          theta[[i]][j] <- counts[[i]][j]/sum(counts[[i]])          
          }
        }
      }
    }  
  theta
  } 

# get states (internal use only)
fullOmega <- function(counts) {
  omega <- list()
  for(i in 1:length(counts)) {
    idim <- dim(counts[[i]])
    if(length(idim)<2) {
      omega[[i]] <- names(counts[[i]])
      } else {
      auxst <- dimnames(counts[[i]])[1]
      attr(auxst,"names") <- NULL      
      omega[[i]] <- unlist(auxst)
      }
    }
  names(omega) <- names(counts)
  omega[sort(names(omega))]
  }

# get the full parent sets (internal use only)
fullPSets <- function(counts) {
  parSet <- list()
  for(i in 1:length(counts)) {
    if(length(dim(counts[[i]]))<2) {
      parSet[[i]] <- character(0)
      } else {
      parSet[[i]] <- sort(names(dimnames(counts[[i]]))[-1])      
      }
    }
  names(parSet) <- names(counts)   
  parSet[sort(names(parSet))]
  }

# find activator nodes (internal use only)
findLambda <- function(counts) {
  parSet <- fullPSets(counts)
  nomi <- names(parSet) 
  nomiLam <- nomi[which(grepl("^Lambda[0-9]{1,}.",nomi))]
  if(length(nomiLam)>0) {
    inam <- c()
    for(i in 1:length(nomiLam)) {
      inam[i] <- strsplit(strsplit(nomiLam[i],"\\.\\.")[[1]][1],"Lambda[0-9]{1,}.")[[1]][2]
      }
    res <- vector("list",length=length(unique(inam)))
    names(res) <- sort(unique(inam))
    for(i in 1:length(res)) {
      res[[i]] <- nomiLam[which(inam==names(res)[i])]
      }    
    res
    } else {
    NULL
    }
  }

# find auxiliary nodes (internal use only)
findAux <- function(counts) {
  parSet <- fullPSets(counts)
  nomi <- names(parSet) 
  nomiAux <- nomi[which(grepl("^AUX[0-9]{1,}.",nomi))]
  if(length(nomiAux)>0) {
    inam <- c()
    for(i in 1:length(nomiAux)) {
      inam[i] <- strsplit(nomiAux[i],"AUX[0-9]{1,}.")[[1]][2]
      }
    res <- vector("list",length=length(unique(inam)))
    names(res) <- sort(unique(inam))
    for(i in 1:length(res)) {
      res[[i]] <- nomiAux[which(inam==names(res)[i])]
      }    
    res
    } else {
    NULL
    }
  }

# find interaction nodes (internal use only)
findInter <- function(counts,resp) {
  parSet <- fullPSets(counts)
  nomi <- names(parSet)
  nomiLam <- nomi[which(grepl("^Lambda[0-9]{1,}.",nomi))]
  nomiAux <- nomi[which(grepl("^AUX[0-9]{1,}.",nomi))]
  nomiOK <- setdiff(names(parSet),c(nomiLam,nomiAux))
  res <- list()
  nomiInt <- c()
  for(i in 1:length(nomiOK)) {
    ipset <- parSet[[nomiOK[i]]]
    ipsOK <- setdiff(ipset,c(nomiLam,nomiAux))
    if(length(ipsOK)>0) {
      if(resp==T) {
        ichld <- nodchil(counts,nomiOK[i])
        iy <- c()
        if(length(ichld)>0) {
           for(j in 1:length(ichld)) {
            iy <- c(iy,strsplit(strsplit(ichld[j],"\\.\\.")[[1]][1],"Lambda[0-9]{1,}.")[[1]][2])
            }
          } else {
          iy <- character(0)  
          }
        res[[nomiOK[i]]] <- list(from=ipsOK,to=unique(iy))
        } else {
        res[[nomiOK[i]]] <- list(from=ipsOK)
        }
      }
    }
  res
  }

# find the relevant nodes for all variables (internal use only)
findRelevant <- function(counts) {
  lambSet <- findLambda(counts)
  nomi <- names(lambSet)
  auxSet <- findAux(counts)
  intSet <- findInter(counts,resp=T)
  intTo <- lapply(intSet,function(z){z$to})
  res <- list()
  for(i in 1:length(nomi)) {
    inam <- nomi[i]
    ipset <- lambSet[[inam]]
    ipar <- c()
    for(j in 1:length(ipset)) {
      ipar <- c(ipar,strsplit(ipset[j],"\\.\\.")[[1]][2])
      }
    ito <- names(which(sapply(intTo,function(z){inam %in% z})==T))
    iint <- intersect(names(intSet),ipar)
    if(length(iint)>0) {
      for(j in 1:length(iint)) {
        ipar <- c(ipar,intSet[[iint[j]]]$from)
        }
      ipar <- setdiff(ipar,iint)
      }
    ires1 <- setdiff(ipar,"omitted")
    ires2 <- c(lambSet[[inam]],auxSet[[inam]],ito,inam)
    res[[i]] <- list(par=sort(unique(ires1)),aux=sort(unique(ires2)))
    }
  names(res) <- nomi
  res
  }

# get the parent set after absorbing auxiliary nodes (internal use only)
absPSets <- function(counts) {
  rel <- findRelevant(counts)
  lapply(rel,function(z){z$par})
  }

# check d-separation
dSepCheck <- function(x,var1,var2,given=NULL) {
  if(!identical(class(x),"cibn")) stop("The first argument must be an object of class 'cibn'",call.=F)
  if(length(var1)!=1||length(var2)!=1) stop("Argument 'var1' and 'var2' must both have length one",call.=F)
  nomi <- names(x$CPTs)
  auxcheck <- setdiff(c(var1,var2,given),nomi)
  if(length(auxcheck)>0) {
    stop("Unknown variables: ",paste(auxcheck,collapse=", "),call.=F)
    }
  Gm <- moralize(ancestralGraph(c(var1,var2,given),as.graphNEL(x,full=F)))
  auxedg <- unlist(edgeList(Gm))
  if(!is.null(auxedg)) {
    xedg <- matrix(auxedg,ncol=2,byrow=T)     
    xedg <- rbind(xedg,xedg[,2:1])
    nomi <- sort(unique(xedg))
    borde <- var1
    reached <- c()
    neighb <- function(x,node) {
      Ne <- c(x[which(x[,1]==node),2],x[which(x[,2]==node),1])
      sort(unique(Ne))
      }
    while(length(borde)>0) {
      reached <- c(reached,borde)
      fan_borde <- c()
      for(i in 1:length(borde)) {
        fan_borde <- c(fan_borde,neighb(xedg,borde[i]))
        }
      borde <- setdiff(fan_borde,c(reached,given))
      if(length(intersect(borde,var2))>0) break()
      }
    ifelse(length(borde)>0, res <- F, res <- T)    
    } else {
    res <- T
    }
  res
  }

# compute a CPT after absorbing auxiliary nodes (internal use only)
cptCalc <- function(counts,resp,xpar,xrelev) {
  omega <- fullOmega(counts)
  idimnam <- omega[c(resp,xpar)]
  icpt <- array(dim=sapply(idimnam,length))
  icpts <- counts[xrelev]
  if(length(xpar)>0) {
    mystr <- ""
    for(j in 1:length(xpar)) {
      ijomeg <- omega[[xpar[j]]]
      ijcpt <- array(rep(1/length(ijomeg),length(ijomeg)),dim=length(ijomeg))
      ijdnam <- list(ijomeg)
      names(ijdnam) <- xpar[j]
      dimnames(ijcpt) <- ijdnam
      class(ijcpt) <- "table"
      icpts[[xpar[j]]] <- ijcpt
      mystr <- paste(mystr,"for(k",j," in 1:",length(ijomeg),") {; ",sep="")
      }
    mystr <- paste(mystr,"ijevid <- list(",paste("omega[[xpar[",1:length(xpar),"]]][k",1:length(xpar),"]",sep="",collapse=","),"); names(ijevid) <- xpar; ",
      "icpt[,",paste("k",1:length(xpar),sep="",collapse=","),"] <- propFun(ibn,resp,ijevid,'marginal')[[1]]; ",sep="")
    mystr <- paste(mystr,paste(rep("}",length(xpar)),collapse="; "),sep="")
    G <- list(CPTs=icpts,prior=icpts,posterior=NULL,description=NULL)
    class(G) <- "cibn"
    ibn <- as.grain(G)
    eval(parse(text=mystr))
    } else {
    G <- list(CPTs=icpts,prior=icpts,posterior=NULL,description=NULL)
    class(G) <- "cibn"
    ibn <- as.grain(G)
    icpt[] <- propFun(ibn,resp,NULL,"marginal")[[1]]
    }
  dimnames(icpt) <- idimnam
  class(icpt) <- "table"
  icpt
  }

# print method for class 'cibn'
print.cibn <- function(x,...) {
  auxpset <- absPSets(x$CPTs)
  nvar <- length(auxpset)
  narc <- sum(sapply(auxpset,length))
  if(is.null(x$posterior)) {
    isUp <- "FALSE"
    } else {
    isUp <- "TRUE"
    }
  #cat("A causal independence Bayesian network with ",nvar," variables and ",narc," edges. Updated = ",isUp,sep="","\n") 
  cat("A causal independence Bayesian network with ",nvar," variables and ",narc," edges",sep="","\n")
  }

# summary method for class 'cibn'
summary.cibn <- function(object,...) {
  auxomeg <- fullOmega(object$CPTs)
  auxpset <- absPSets(object$CPTs)
  omegOK <- auxomeg[sort(names(auxpset))]
  auxlen1 <- sapply(auxpset,length)
  auxlen2 <- sapply(omegOK,length)
  nvar <- length(auxpset)
  narc <- sum(sapply(auxpset,length))
  if(is.null(object$posterior)) {
    isUp <- "FALSE"
    } else {
    isUp <- "TRUE"
    }
  cat(" Name: ",deparse(substitute(object)),sep="","\n")
  cat(" Variables: ",nvar,sep="","\n")
  cat(" Edges: ",narc,sep="","\n")
  cat(" Average size of parent sets: ",round(mean(auxlen1),2),sep="","\n")
  cat(" Average dimension of sample spaces: ",round(mean(auxlen2),2),sep="","\n")
  #cat(" Updated: ",isUp,sep="")
  }

# plot method for class 'cibn'
plot.cibn <- function(x,full=FALSE,...) {
  G <- as.graphNEL(x, full=full)
  plot(G,...)
  }

# get variable names
getVariables <- function(x) {
  if(!identical(class(x),"cibn")) stop("The argument must be an object of class 'cibn'",call.=F)
  auxpset <- absPSets(x$CPTs)
  sort(names(auxpset))
  }

# get description
getDescription <- function(x) {
  if(!identical(class(x),"cibn")) stop("The argument must be an object of class 'cibn'",call.=F)
  x$description
  }

# get sample spaces
getStates <- function(x) {
  if(!identical(class(x),"cibn")) stop("The argument must be an object of class 'cibn'",call.=F)
  auxomeg <- fullOmega(x$CPTs)
  auxpset <- absPSets(x$CPTs)
  auxomeg[sort(names(auxpset))]
  }

# get variable types
getTypes <- function(x) {
  if(!identical(class(x),"cibn")) stop("The argument must be an object of class 'cibn'",call.=F)
  omega <- fullOmega(x$CPTs)
  nomi <- getVariables(x)
  xbin <- intersect(names(omega)[which(sapply(omega,length)==2)],nomi)
  res <- rep("GRAD",length(nomi))
  names(res) <- nomi
  #if(length(xbin)>0) res[xbin] <- "NOM"
  auxNod <- findAux(x$CPTs)
  if(length(auxNod)>0) {
    for(i in 1:length(auxNod)) {
      inam <- names(auxNod)[i]
      iomeg <- omega[c(inam,auxNod[[i]])]
      ilen <- sapply(iomeg,length)
      iind <- which(ilen==ilen[1])
      if(length(iind)!=length(ilen)) {
        ibin <- which(ilen[-1]==2)
        if(length(ibin)>0 && length(ibin)==length(ilen)-1) {
          res[inam] <- "NOM"
          } else {
          res[inam] <- "DGRAD"
          }
        }
      }
    }
  res  
  }

# get the parent sets
getParSets <- function(x, full=FALSE) {
  if(!identical(class(x),"cibn")) stop("The argument must be an object of class 'cibn'",call.=F)
  if(full) fullPSets(x$CPTs) else absPSets(x$CPTs)
  }

# find the child sets (internal use only)
chldsets <- function(x) {
  pset <- inEdges(x)
  findchld <- function(xname,ps) {names(which(sapply(ps,function(z){xname %in% z})==T))}
  nomi <- names(pset)
  res <- lapply(nomi,findchld,ps=pset)
  names(res) <- nomi
  res
  }

# get edges (internal use only)
getEdges <- function(x) {
  pset <- getParSets(x)
  nomi <- names(pset)
  edmat <- matrix(NA,nrow=0,ncol=2)
  for(i in 1:length(pset)) {
    inam <- nomi[i]
    if(length(pset[[inam]] )>0) edmat <- rbind(edmat,cbind(pset[[inam]],rep(inam,length(pset[[inam]]))))
    }
  colnames(edmat) <- c("from","to")
  edmat
  }

# get CPT
getCPT <- function(x, variables=NULL) {
  if(!identical(class(x),"cibn")) stop("The first argument must be an object of class 'cibn'",call.=F)
  nomi <- getVariables(x)
  if(is.null(variables)) variables <- nomi
  variables <- intersect(nomi,variables)
  if(length(variables)==0) stop("All the variables provided to argument 'variables' are unknown",call.=F)
  rel <- findRelevant(x$CPTs)
  res <- list()
  for(i in 1:length(variables)) {
    irel <- rel[[variables[i]]]
    res[[i]] <- cptCalc(x$CPTs,variables[i],irel$par,irel$aux)
    }
  names(res) <- variables
  res
  }

# draw a sample
sample.cibn <- function(x,nsam,seed=NULL) {
  if(!identical(class(x),"cibn")) stop("The first argument must be an object of class 'cibn'",call.=F)
  if(nsam<1 || round(nsam)!=nsam) {
    stop("Argument 'nsam' must be a positive integer number",call.=F)
    }
  gnet <- as.grain(x)
  mysam <- simulate(gnet,nsim=nsam,seed=seed)
  omega <- fullOmega(x$CPTs)
  for(i in 1:ncol(mysam)) {
    yst <- omega[[colnames(mysam)[i]]]
    mysam[,i] <- factor(yst[mysam[,i]],levels=yst)
    }
  nomiOK <- getVariables(x)
  mysam[,nomiOK]
  }

# convert from class 'cibn' to class 'graphNEL'
as.graphNEL <- function(x, full=FALSE) {
  if(!identical(class(x),"cibn")) stop("The argument must be an object of class 'cibn'",call.=F)
  if(full) {
    parSet <- fullPSets(x$CPTs)
    } else {
    parSet <- absPSets(x$CPTs)
    }
  nomi <- names(parSet)
  if(sum(sapply(parSet,length))==0) {
    graphNEL(nodes=nomi)
    } else {
    from <- to <- c()
    for(i in 1:length(parSet)) {
      inam <- nomi[i]
      if(length(parSet[[inam]])>0) {
        to <- c(to,rep(inam,length(parSet[[inam]])))
        from <- c(from,parSet[[inam]])
        }
      }
    ftM2graphNEL(cbind(from,to),V=nomi)
    }
  }

# convert from class 'cibn' to class 'grain'
as.grain <- function(x) {
  if(!identical(class(x),"cibn")) stop("The argument must be an object of class 'cibn'",call.=F)
  estim <- x$CPTs                       
  nomi <- names(estim)            
  parSet <- fullPSets(x$CPTs)
  omega <- fullOmega(x$CPTs) 
  theta <- list()             
  for(i in 1:length(estim)) {
    inam <- nomi[i]
    ipar <- parSet[[inam]]
    if(length(ipar)==0) {
      ithetastr <- paste("theta[[i]] <- cptable(~",inam,
        ",values=c(",paste(estim[[i]],collapse=",",sep=""),
        "),levels=omega[[inam]])",sep="")
      } else {
      ithetastr <- paste("theta[[i]] <- cptable(~",inam,
        "|",paste(ipar,collapse="+",sep=""),
        ",values=c(",paste(estim[[i]],collapse=",",sep=""),
        "),levels=omega[[inam]])",sep="")
      }                                
    eval(parse(text=ithetastr)) 
    }                           
  plist <- compileCPT(theta)
  grain(plist)
  }

# propagation function (internal use only)
propFun <- function(G,target,evidence,type) {
  pE <- c()
  if(length(evidence)>0) {
    G <- setEvidence(G,names(evidence),evidence)
    pE <- attr(G$equipot,"pEvidence")
    }
  #compile(G)
  out <- querygrain(G,nodes=target,type=type,exclude=F)
  if(type=="marginal") {
    for(i in 1:length(out)) {
      out[[i]] <- c(out[[i]])
      }
    res <- out[target]
    } else {
    res <- aperm(out,target)
    }
  attr(res,"pEvidence") <- pE
  res
  }

# evidence propagation
query.cibn <- function(x,target=NULL,evidence=NULL,type="marginal") {
  if(!identical(class(x),"cibn")) stop("The first argument must be an object of class 'cibn'",call.=F)
  nomi <- getVariables(x)
  if((type %in% c("marginal","joint"))==F) stop("Argument 'type' must be either 'marginal' or 'joint'",call.=F)
  if(length(target)==0) target <- nomi
  unknVar1 <- setdiff(target,names(x$CPTs))
  if(length(unknVar1)>0) {
    warning("Unknown variables in argument 'target' have been ignored: ",paste(unknVar1,collapse=", "),call.=F)
    target <- setdiff(target,unknVar1)
    }
  if(!is.null(evidence) && !is.list(evidence)) stop("Argument 'evidence' must be a list",call.=F)
  unknVar2 <- setdiff(names(evidence),names(x$CPTs))
  if(length(unknVar2)>0) {
    warning("Unknown variable names in argument 'evidence' have been ignored: ",paste(unknVar2,collapse=", "),call.=F)
    evidence <- evidence[setdiff(names(evidence),unknVar2)]
    }
  if(length(evidence)>0) {
    auxdupl <- names(evidence)[duplicated(names(evidence))]
    if(length(auxdupl)>0) stop("Duplicated variable name in argument 'evidence': ",auxdupl[1],call.=F)
    omega <- fullOmega(x$CPTs)
    for(i in 1:length(evidence)) {
      inam <- names(evidence)[i]
      ichk <- setdiff(evidence[[i]],omega[[inam]])
      if(length(ichk)>0) {
        warning("Unknown states in argument 'evidence' for variable ",inam," have been ignored: ",paste(ichk,collapse=", "),call.=F)
        if(length(evidence)==1) evidence <- NULL
        }
      }
    }
  gnet <- as.grain(x)
  propFun(gnet,target,evidence,type)
  }
