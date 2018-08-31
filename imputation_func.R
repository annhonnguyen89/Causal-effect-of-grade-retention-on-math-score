intertreat.cat_1 <- function(){
  addedinterdat <- NULL
  addedintermeth <- NULL
  n <- nrow(widedfnoid)
  
  createinter <- function(fixstr1, numlevtreat, fixstr2, numlevcat){
    # create interaction (treat_Y1, 3, sex, 2) --> leftside: [treat_Y1.1.sex.1, treat_Y1.2.sex.1,...], rightside: [~I(treat_Y1.1*sex.1), ~I(treat_Y1.2*sex.1),...]
    numlev1 <- numlevtreat - 1
    numlev2 <- numlevcat - 1
    leftside <- paste0(paste0(fixstr1, ".", c(1:numlev1)), ".", paste0(fixstr2,"." ,c(1:numlev2)))
    rightside <- paste0("~I(",paste0(fixstr1, ".", c(1:numlev1)), "*", paste0(fixstr2, ".", c(1:numlev2)), ")")
    return(rbind(leftside, rightside))
  }
  
  treat <- "treat_Y1"
  numlevtreat = 2
  
  addedintermeth <- NULL
  catlist = c("Sex", "Hlang", "Etn", "Net_Y1", "Net_Y2", "Net_Y3", "Net_Y4", "Net_Y5", "Net_Y6")
  # interaction treat_Y1:Sex
  for ( i in catlist){
    numlevcat <- length(levels(widedfnoid[,i]))
    t <- createinter(treat, numlevtreat, i, numlevcat)
    
    if(is.null(addedintermeth)){
      addedintermeth <- t
    } else{
      addedintermeth <- cbind(addedintermeth, t)  
    }
    
    # Extend data matrix (add interaction columns with NA value) 
    newvar <- t["leftside",]
    for(name in newvar){
      c <- paste0(name, " <- rep(NA, n)")
      v <- eval(parse(text=c))
      addedinterdat <- cbind(addedinterdat, v)
      m = ncol(addedinterdat)
      colnames(addedinterdat)[m] = name
    }
  }
  return(list(addedintermeth, addedinterdat))
}

modify_pred_inter_2 <- function(pred){
  treat <- "treat_Y1"
  numlevtreat = 2
  
  addedintermeth <- NULL
  catlist = c("Sex", "Hlang", "Etn", "Net_Y1", "Net_Y2", "Net_Y3", "Net_Y4", "Net_Y5", "Net_Y6")
  
  numlev1 = 2-1 # number of levels of treat_Y1
  for ( i in catlist){
    numlev2 = length(levels(widedfnoid[,i])) - 1 # number of levels of cat variable
    leftside <- paste0(paste0("treat_Y1.", c(1:numlev1)), ".", paste0(i,"." ,c(1:numlev2)))
    
    pred["treat_Y1", leftside] = 0
    pred[i, leftside] = 0
    
    pred[leftside, "treat_Y1"] = 1
    pred[leftside, i] = 1
  }
  return(pred)
}

modify_pred_inter_3 <- function(pred){
  
  #inter between treat_Yi and fixed cat
  # inter between treat_Yi and educat_Yi
  for (y in c(1:6)){
    treat <- paste0("treat_Y",y)
    numlevtreat = 2
    
    addedintermeth <- NULL
    #catlist = c("Sex", "Hlang", "Etn", "Net_Y1", "Net_Y2", "Net_Y3", "Net_Y4", "Net_Y5", "Net_Y6")
    
    numlev1 = 2-1 # number of levels of treat_Y1
    for ( i in fixedcatvar){
      numlev2 = length(levels(widedfnoid[,i])) - 1 # number of levels of cat variable
      leftside <- paste0(paste0("treat_Y1.", c(1:numlev1)), ".", paste0(i,"." ,c(1:numlev2)))
      
      pred[treat, leftside] = 0
      pred[i, leftside] = 0
      
      pred[leftside, treat] = 1
      pred[leftside, i] = 1
    }
    
    edulist <- paste0(catedu, "_Y", y)
    
    for ( i in edulist){
      numlev2 = length(levels(widedfnoid[,i])) - 1 # number of levels of cat variable
      leftside <- paste0(paste0(treat, ".", c(1:numlev1)), ".", paste0(i,"." ,c(1:numlev2)))
      
      pred[treat, leftside] = 0
      pred[i, leftside] = 0
      
      pred[leftside, treat] = 1
      pred[leftside, i] = 1
    }
  }
  
  return(pred)
}


interTreat.cat <- function(y){
  # interaction between treat_Yy and fixed categorical variables
  
  addedinterdat <- NULL
  addedintermeth <- NULL
  n <- nrow(widedfnoid)
  
  createinter <- function(fixstr1, numlevtreat, fixstr2, numlevcat){
    numlev1 <- numlevtreat - 1
    numlev2 <- numlevcat - 1
    leftside <- paste0(paste0(fixstr1, ".", c(1:numlev1)), ".", paste0(fixstr2,"." ,c(1:numlev2)))
    rightside <- paste0("~I(",paste0(fixstr1, ".", c(1:numlev1)), "*", paste0(fixstr2, ".", c(1:numlev2)), ")")
    return(rbind(leftside, rightside))
  }
  
  treat <- paste0("treat_Y", y)
  
  ############################### create interaction between treat and fixed categorical variable
  for (name in fixedcatvar){
    catvar <- name
    
    numlevtreat <- length(unique(widedfnoid[,treat][!is.na(widedfnoid[,treat])]))
    numlevcat <- length(unique(widedfnoid[,catvar][!is.na(widedfnoid[,catvar])]))
    
    t <- createinter(treat, numlevtreat, catvar, numlevcat)
    
    if(is.null(addedintermeth)){
      addedintermeth <- t
    } else{
      addedintermeth <- cbind(addedintermeth, t)  
    }
    
    # Extend data matrix (add interaction columns with NA value) 
    newvar <- t["leftside",]
    
    for(name in newvar){
      c <- paste0(name, " <- rep(NA, n)")
      v <- eval(parse(text=c))
      addedinterdat <- cbind(addedinterdat, v)
      m = ncol(addedinterdat)
      colnames(addedinterdat)[m] = name
    }
  }
  
  ####################### interaction between treatment and timevarying variables Koep_y1, Koep_y2, ... , PM_Y1, PM_Y2, ... 
  for(name in catedu){
    catvar <- paste0(name, "_Y", y)
    
    numlevtreat <- length(unique(widedfnoid[,treat][!is.na(widedfnoid[,treat])]))
    numlevcat <- length(unique(widedfnoid[,catvar][!is.na(widedfnoid[,catvar])]))
    
    t <- createinter(treat, numlevtreat, catvar, numlevcat)
    
    if(is.null(addedintermeth)){
      addedintermeth <- cbind(addedintermeth, t)
    } else{
      addedintermeth <- cbind(addedintermeth, t)  
    }
    
    # Extend data matrix (add interaction columns with NA value) 
    newvar <- t["leftside",]
    
    for(name in newvar){
      c <- paste0(name, " <- rep(NA, n)")
      v <- eval(parse(text=c))
      addedinterdat <- cbind(addedinterdat, v)
      m = ncol(addedinterdat)
      colnames(addedinterdat)[m] = name
    }
  }
  
  return(list(addedintermeth, addedinterdat))
}


createinter <- function(fixstr1, num, fixstr2, avg){
  # Create formula for interaction between each level of fixstr1 and fixstr2
  # ex: fixstr1 = "Koep_Y0", num = 7 (Koep_Y0 has 7 levels, fixstr2: "Math_Y0", avg = mean(Math_Y0))
  # result: 
  #    leftside: "Koep_Y0.1.Math_Y0", "Koep_Y0.2.Math_Y0", ...
  #    rightside: "~I(Koep_Y0.7*(Math_Y0-51.51465)", ...  
  num <- num - 1
  leftside <- paste0(paste0(fixstr1, ".", c(1:num)), ".", fixstr2)
  rightside <- paste0("~I(",paste0(fixstr1, ".",c(1:num)), "*(", fixstr2, "-", avg, "))")
  return(rbind(leftside, rightside))
}

interKoep.score <- function(y){
  addedinterdat <- NULL
  addedintermeth <- NULL
  n <- nrow(widedfnoid)
  
  for (name in scorename){
    koep <- paste0("Koep_Y",y)
    math <- paste0(name, "_Y", y)
    meanMath <- mean(widedfnoid[,math], na.rm = TRUE)
    
    level <- length(unique(widedfnoid[,koep][!is.na(widedfnoid[,koep])]))
    
    t <- createinter(koep, level, math, meanMath)
    z <- matrix(nrow = 1, ncol = level)
    
    if(is.null(addedintermeth)){
      addedintermeth <- t
    } else{
      addedintermeth <- cbind(addedintermeth, t)  
    }
    
    # Extend data matrix (add interaction columns with NA value) 
    newvar <- t["leftside",]
    
    for(name in newvar){
      
      c <- paste0(name, " <- rep(NA, n)")
      v <- eval(parse(text=c))
      
      if(is.null(addedinterdat)){
        addedinterdat <- v
      } else {
        addedinterdat <- cbind(addedinterdat, v)  
      }
    }
  }
  return(list(addedintermeth, addedinterdat))
}




modifypredmat <- function(fixstr1, num ,fixstr2list, pred){
  # pred[c("Koep_Y0","Math_Y0"), c("Koep_Y0.1.Math_Y0", "Koep_Y0.2.Math_Y0", "Koep_Y0.3.Math_Y0", "Koep_Y0.4.Math_Y0", "Koep_Y0.5.Math_Y0", "Koep_Y0.6.Math_Y0", "Koep_Y0.7.Math_Y0")] = 0
  
  num <- num - 1
  for (i in fixstr2list){
      leftside <- paste0(paste0(fixstr1, ".", c(1:num)), ".", i)
      pred[c(fixstr1, i), leftside] = 0

      # may be commented out
      pred[leftside, c(fixstr1, i)] = 1
  }
  return (pred)
}

# if x1 and x2 are correlated --> x1 is also correlated with t*x2
# 

modifypred_inter <- function(pred){

  rcount <- ncol(pred)
  predcol <- colnames(pred)
  usablecol <- colnames(usablepro)
  pvalcol <- colnames(cat_cont_pval)
  contcorrcol <- colnames(contcorr)
 
  # check correlation cont-cont
  for( i in contcorrcol){
    for (j in contcorrcol){
      tryCatch(
        if( i!=j && is.na(mutcorr[i, j]) == FALSE && is.na(contcorr[i, j]) == FALSE && abs(mutcorr[i, j]) > 0.4 && abs(contcorr[i, j]) > 0.4 ){
          pred[i,j] = 1
          pred[j,i] = 1
        } else {
          pred[i,j] = 0
          pred[j,i] = 0
        }, error = function(e) {print(paste(e,"-", i, "-", j))}
      )
    }
  }
  
  
  # check pvalue of relation cat_cont
  for (cat in catvars){
    for (cont in contvars){
      if(is.na(cat_cont_pval[cat, cont]) == FALSE && cat_cont_pval[cat, cont] < 0.05){
        pred[cont, cat] = 1
        pred[cat, cont] = 1
      } else {
        pred[cont, cat] = 0
        pred[cat, cont] = 0
      }
    }
  }
  
  # check pvalue of relation cat_cat
  for (cat1 in catvars){
    for (cat2 in catvars){
      if(cat1 != cat2 && is.na(cat_cat_pval[cat1, cat2]) == FALSE && cat_cat_pval[cat1, cat2] < 0.05){
        pred[cat1, cat2] = 1
        pred[cat2, cat1] = 1
      } else {
        pred[cat1, cat2] = 0
        pred[cat2, cat1] = 0
      }
    }
  }
      
  # check usable
  for(i in usablecol){
    for (j in usablecol){
      if(i!=j && (abs(usablepro[i, j]) < 0.30 || is.na(usablepro[i, j]) == TRUE)){
        pred[i, j] = 0
      } else {
        pred[i, j] = 1
      }
    }
  }
  
  
  # alltermpred <- colnames(pred)
  alltermpred <- colnames(pred)
  for (i in alltermpred){
    if(grepl(".", i, fixed=TRUE)){
      l <- unlist(strsplit(i, ".", fixed = TRUE))
      for(j in c(catvars, contvars)){
        tryCatch(
          if ( l[1] != j && l[3] != j && pred[l[3], j] == 0 && pred[l[1], j] == 0 ) {
            pred[i, j] == 0
          } else {pred[i, j] == 1}, error = function(e){print (paste(c(l[1],l[2],l[3]), collapse = "-"))}
        )
      }
    }
  }
  return (pred)

}

remove_inter_frompred <- function(pred){
  library(stringr)
  rnames <- rownames(pred)
  cnames <- colnames(pred)
  
  for (i in rnames){
    for (j in cnames){
      if (str_count(j,"\\.") == 2){
        loca <- str_locate_all(pattern ='\\.', j)
        end1 = loca[[1]][1] - 1
        var1 = substr(j, 1, end1)
        start2 = loca[[1]][2] + 1
        end2 = nchar(j)
        var2 = substr(j, start2, end2)
        if(pred[i,var1] == 1 && pred[i,var2] == 1){
          pred[i,j] = 1  
        } else{
          pred[i,j] = 0
        }
      } else if(str_count(j,"\\.") == 3){
        loca <- str_locate_all(pattern ='\\.', j)
        end1 = loca[[1]][1] - 1
        var1 = substr(j, 1, end1)
        start2 = loca[[1]][2] + 1
        end2 = loca[[1]][3] - 1
        var2 = substr(j, start2, end2)
        if(pred[i,var1] == 1 && pred[i,var2] == 1){
          pred[i,j] = 1  
        } else {
          pred[i,j] = 0
        }
      }
    }
  }
  return(pred)
}

modifypredmat1 <- function(pred){
  
  #pred <- modify_pred_inter_2(pred)
  #pred <- modify_pred_inter_3(pred)
  
  rcount <- ncol(pred)
  predcol <- colnames(pred)
  usablecol <- colnames(usablepro)
  pvalcol <- colnames(cat_cont_pval)
  contcorrcol <- colnames(contcorr)
  
  # check correlation cont-cont
  for( i in contcorrcol){
    for (j in contcorrcol){
      tryCatch(
        #if( i!=j && is.na(mutcorr[i, j]) == FALSE && is.na(contcorr[i, j]) == FALSE && abs(mutcorr[i, j]) > 0.1 && abs(contcorr[i, j]) > 0.1 ){
        if( i!=j && is.na(contcorr[i, j]) == FALSE && abs(contcorr[i, j]) > 0.3 ){
          pred[i,j] = 1
          pred[j,i] = 1
        } else {
          pred[i,j] = 0
          pred[j,i] = 0
        }, error = function(e) print(paste(e,"-", i, "-", j))
      )
    }
  }
  
  # check pvalue of relation cat_cont
  for (cat in catvars){
    for (cont in contvars){
      if(is.na(cat_cont_pval[cat, cont]) == FALSE && cat_cont_pval[cat, cont] < 0.05){
        pred[cont, cat] = 1
        pred[cat, cont] = 1
      } else {
        pred[cont, cat] = 0
        pred[cat, cont] = 0
      }
    }
  }
  
  # check pvalue of relation cat_cat
  for (cat1 in catvars){
    for (cat2 in catvars){
      if(cat1 != cat2 && is.na(cat_cat_pval[cat1, cat2]) == FALSE && cat_cat_pval[cat1, cat2] < 0.05){
        pred[cat1, cat2] = 1
        pred[cat2, cat1] = 1
      } else {
        pred[cat1, cat2] = 0
        pred[cat2, cat1] = 0
      }
    }
  }
  
  # check usable
  for(i in usablecol){
    for (j in usablecol){
      if(i!=j && (abs(usablepro[i, j]) < 0.3 || is.na(usablepro[i, j]) == TRUE)){
        pred[i, j] = 0
      }
    }
  }
  
  #pred <- remove_inter_frompred(pred)
  return(pred)
}



squeeze <- function(ini, varnames){
  post <- ini$post
  dat <- ini$data
  #post["chl"] <- "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(100,300))"
  for (v in varnames){
    if(class(dat[,v]) == 'numeric'){
      mi <- min(dat[,v], na.rm = TRUE)
      ma <- max(dat[,v], na.rm = TRUE)
      
      statement <- paste0("imp[[j]][,i] <- squeeze(imp[[j]][,i],c(",mi,",",ma,"))")
      
      #post[,v] <- "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(100,300))"
      post[v] <- statement
    }
  }
  return(post)
}





