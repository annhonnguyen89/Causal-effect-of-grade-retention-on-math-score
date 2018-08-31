# g-estimate function

library(nnet)
library(DataCombine)
library(DTRreg)
library(mice)
library(logistf)
library(caret)
library(ROCR)
library(pROC)
library(heatmapFit)

check_perfect_sep <- function(cat1, cat2){
  tab = table(cat1, cat2)
}

detransform <- function(imp){
  imp$Math_Y0 <- (imp$Math_Y0)^(1/2)
  imp$Math_Y1 <- (imp$Math_Y1)^(1/2)
  imp$Math_Y2 <- (imp$Math_Y2)^(1/2.5)
  imp$Math_Y3 <- (imp$Math_Y3)^(1/2.5)
  imp$Math_Y4 <- (imp$Math_Y4)^(1/2.5)
  imp$Math_Y5 <- (imp$Math_Y5)^(1/2.5)
  imp$Math_Y6 <- (imp$Math_Y6)^(1/2.5)
  return(imp)
}


select_model <- function(null.mod, full.mod){
  data <- full.mod$data
  data <- na.omit(data)
  mod <- step(null.mod, scope = list(lower = null.mod, upper = full.mod), direction = "both", trace = 1)
  
  diagnose <- function(mod){
    # http://data.princeton.edu/wws509/notes/c3s8.html
    # diagnose residual of logistic model:
    # - Deviance residual: distribution of dev residual coresponding to y = 1 is similar to y = 0
    # - Working residual: 
    # - pearson residual: 
    
    fit <- mod$fitted.value
    
    #working
    workres<- resid(mod,type="working")
    plot(workres, col=c("blue","red")[mod$y])
    plot(fit, workres, col=c("blue","red")[mod$y])
    
    #deviance
    res_dev<- resid(mod, type="dev")
    
    plot(fit, res_dev, col=c("blue","red")[mod$y])
    plot(res_dev, col=c("blue","red")[mod$y])
    
    #response
    resres <- resid(mod, type="resp")
    plot(fit, resres, col=c("blue","red")[mod$y])
    plot(resres, col=c("blue","red")[mod$y])
    
    #pearson
    pearres <- resid(mod, type="pear")
    plot(fit, pearres, col = c("blue","red")[mod$y])
    plot(pearres, col = c("blue","red")[mod$y])
  }
  
  #diagnose(mod)
  return(mod)
}

diagnose_logistic <- function(mod){
  # http://data.princeton.edu/wws509/notes/c3s8.html
  # diagnose residual of logistic model:
  # - Deviance residual: distribution of dev residual coresponding to y = 1 is similar to y = 0
  # - Working residual: 
  # - pearson residual: 
  
  fit <- mod$fitted.value
  
  Hosmer.Lemeshow(mod)
  heatmap.fit(imp$treat_Y1, fit, reps = 1000, compress.obs = TRUE, 
              init.grid = 2000, span.l = "aicc") # error ??
  
  #working
  workres<- resid(mod,type="working")
  plot(workres, col=c("blue","red")[mod$y])
  plot(fit, workres, col=c("blue","red")[mod$y])
  
  #deviance
  res_dev<- resid(mod, type="dev")
  
  plot(fit, res_dev, col=c("blue","red")[mod$y])
  plot(res_dev, col=c("blue","red")[mod$y])
  
  #response
  resres <- resid(mod, type="resp")
  plot(fit, resres, col=c("blue","red")[mod$y])
  plot(resres, col=c("blue","red")[mod$y])
  
  #pearson
  pearres <- resid(mod, type="pear")
  plot(fit, pearres, col = c("blue","red")[mod$y])
  plot(pearres, col = c("blue","red")[mod$y])
  
}

creat_trt <- function(time, removed=NULL, simple = TRUE){
  # school providers in year time, same year with the retention
  allvar <- fixedvar
  t = time - 1
  #names <- setdiff(nametimevarying, c("T"))
  acacore1 <- c("TC_Lang",  "TC_Math", "Math")
  names <- c(charscore, acacore1)
  
  if(simple == TRUE){
    name_y <- sapply(names, function(x){paste0(x,"_Y",t)})
    allvar <- c(allvar, name_y)
  } else{
    for(i in c(0:t)){
      name_y <- sapply(names, function(x){paste0(x,"_Y",i)})
      allvar <- c(allvar, name_y)
    }  
  }
  
  schoolvar <- sapply(c("Koep", "PM", "Schooltype", "Net", "PercGOK"), function(x){paste0(x, "_Y", time)})
  allvar <- c(allvar, schoolvar)
  allvar <- allvar[!allvar %in% removed]
  allvar = setdiff(allvar, "treat_Y0")
  formular <- paste0(paste0("treat_Y",time, " ~ "), paste0(allvar, collapse = "+"))
  return(formular)  
}

creat_trt1 <- function(time, removed=NULL, simple = TRUE){
  # the school providers are the one in year time, 1 year before the year of retention
  
  allvar <- fixedvar
  t = time - 1
  #names <- setdiff(nametimevarying, c("T"))
  acacore1 <- c("TC_Lang",  "TC_Math", "Math")
  names <- c(charscore, acacore1)
  
  if(simple == TRUE){
    name_y <- sapply(names, function(x){paste0(x,"_Y",t)})
    allvar <- c(allvar, name_y)
  } else{
    for(i in c(0:t)){
      name_y <- sapply(names, function(x){paste0(x,"_Y",i)})
      allvar <- c(allvar, name_y)
    }  
  }
  
  schoolvar <- sapply(c("Koep", "PM", "Schooltype", "Net", "PercGOK"), function(x){paste0(x, "_Y", t)})
  allvar <- c(allvar, schoolvar)
  allvar <- allvar[!allvar %in% removed]
  allvar = setdiff(allvar, "treat_Y0")
  formular <- paste0(paste0("treat_Y",time, " ~ "), paste0(allvar, collapse = "+"))
  return(formular)  
}

creat_blip <- function(time_treat, time_math, removed = NULL){
  time <- time_treat
  treat <- paste0("treat_Y", time)
  math <- paste0("Math_Y", time_math)
  ps <- paste0("ps_Y", time_treat)
  
  #treat_fixvar <- sapply(fixedvar, function(x){paste0(treat, ":", x)})
  t = time - 1
  names <- setdiff(nametimevarying, c("T"))
  
  name_y = NULL
  for (i in c(0:t)){
    if(is.null(name_y)){
      name_y <- sapply(names, function(x){paste0(x,"_Y",i)})
    } else {
      temp <- sapply(names, function(x){paste0(x,"_Y",i)})
      name_y <- c(name_y, temp)
    }
  }
  
  name_y <- setdiff(name_y, "treat_Y0")
  name_y <- setdiff(name_y, c("Anx_Y1", "Gap_Y1", "Anx_Y3", "Gap_Y4", "Gap_Y6"))
  
  #schoolvar <- sapply(c("Koep", "PM", "Schooltype", "Net", "PercGOK"), function(x){paste0(treat,":",x, "_Y", time)})
  schoolvar <- sapply(c("Koep", "PM", "Schooltype", "Net", "PercGOK"), function(x){paste0(x, "_Y", time)})
  
  allvar <- c(fixedvar, treat , name_y, schoolvar, ps)
  allvar <- allvar[! allvar %in% removed]
  
  formular <- paste0(math, "~ ", paste0(allvar, collapse = "+"))
  return(formular)
}


creat_blip_1 <- function(time_treat, time_math, removed = NULL){
  time <- time_treat
  treat <- paste0("treat_Y", time)
  math <- paste0("Math_Y", time_math)
  
  treat_fixvar <- sapply(fixedvar, function(x){paste0(treat, "*", x)})
  t = time - 1
  names <- setdiff(nametimevarying, c("T"))
  
  name_y = NULL
  for (i in c(0:t)){
    if(is.null(name_y)){
      name_y <- sapply(names, function(x){paste0(treat,"*",x,"_Y",i)})
    } else {
      temp <- sapply(names, function(x){paste0(treat,"*",x,"_Y",i)})
      name_y <- c(name_y, temp)
    }
  }
  
  name_y <- setdiff(name_y, paste0(treat,"*treat_Y0"))
  name_y <- setdiff(name_y, paste0(treat, "*", c("Anx_Y1", "Gap_Y1", "Anx_Y3", "Gap_Y4", "Gap_Y6")))
  
  schoolvar <- sapply(c("Koep", "PM", "Schooltype", "Net", "PercGOK"), function(x){paste0(treat,"*",x, "_Y", time)})
  allvar <- c(treat_fixvar, treat , name_y, schoolvar)
  allvar <- allvar[! allvar %in% removed]
  
  formular <- paste0(math, "~ ", paste0(allvar, collapse = "+"))
  return(formular)
}


creat_blip_2 <- function(time_treat, time_math, removed = NULL){
  time <- time_treat
  treat <- paste0("treat_Y", time)
  math <- paste0("Math_Y", time_math)
  
  treat_fixvar <- fixedvar
  t = time - 1
  names <- setdiff(nametimevarying, c("T"))
  
  name_y = NULL
  for (i in c(0:t)){
    if(is.null(name_y)){
      name_y <- sapply(names, function(x){paste0(x,"_Y",i)})
    } else {
      temp <- sapply(names, function(x){paste0(x,"_Y",i)})
      name_y <- c(name_y, temp)
    }
  }
  
  name_y <- setdiff(name_y, "treat_Y0")
  name_y <- setdiff(name_y, c("Anx_Y1", "Gap_Y1", "Anx_Y3", "Gap_Y4", "Gap_Y6"))
  
  schoolvar <- sapply(c("Koep", "PM", "Schooltype", "Net", "PercGOK"), function(x){paste0(x, "_Y", time)})
  allvar <- c(treat_fixvar, treat , name_y, schoolvar)
  allvar <- allvar[! allvar %in% removed]
  
  formular <- paste0(math, "~ ", paste0(allvar, collapse = "+"))
  res <- c(formular, allvar)
  return(res)
}


creat_treatfree <- function(time_treat, time_math, removed = NULL){
  # the treatment-free model of treat_Y2 on Math_Y4:
  #      baseline + all measures in Y0, Y1 + Koep_Y2 + PM_Y2 + schooltype_Y2 + Net_Y2 + PercGOK_Y2
  
  time <- time_treat
  math <- paste0("Math_Y", time_math)
  names <- setdiff(nametimevarying, c("T"))
  t = time_treat - 1
  
  
  treat_var = paste0("treat_", time_treat)
  
  removed = c(removed, treat_var)
  
  name_y = NULL
  if(is.null(name_y)){
    name_y <- sapply(names, function(x){paste0(x, "_Y", t)})
  } else {
    temp <- sapply(names, function(x){paste0(x, "_Y", t)})
    name_y <- c(name_y, temp)
  }
  
  name_y <- setdiff(name_y, "treat_Y0")
  schoolvar <- sapply(c("Koep", "PM", "Schooltype", "Net", "PercGOK"), function(x){paste0(x, "_Y", time_treat)})
  allvar <- c(fixedvar, name_y, schoolvar)
  #allvar <- c(fixedvar, name_y)
  allvar <- allvar[! allvar %in% removed]
  formular <- paste0(math, "~ ", paste0(allvar, collapse = "+"))
  reslist = c(formular, allvar)
  return(reslist)
}


select_model <- function(null.mod, full.mod){
  data <- full.mod$data
  data <- na.omit(data)
  mod <- step(null.mod, scope = list(lower = null.mod, upper = full.mod), direction = "both", trace = 1)
  
  diagnose <- function(mod){
    # http://data.princeton.edu/wws509/notes/c3s8.html
    # diagnose residual of logistic model:
    # - Deviance residual: distribution of dev residual coresponding to y = 1 is similar to y = 0
    # - Working residual: 
    # - pearson residual: 
    
    fit <- mod$fitted.value
    
    #working
    workres<- resid(mod,type="working")
    plot(workres, col=c("blue","red")[mod$y])
    plot(fit, workres, col=c("blue","red")[mod$y])
    
    #deviance
    res_dev<- resid(mod, type="dev")
    
    plot(fit, res_dev, col=c("blue","red")[mod$y])
    plot(res_dev, col=c("blue","red")[mod$y])
    
    #response
    resres <- resid(mod, type="resp")
    plot(fit, resres, col=c("blue","red")[mod$y])
    plot(resres, col=c("blue","red")[mod$y])
    
    #pearson
    pearres <- resid(mod, type="pear")
    plot(fit, pearres, col = c("blue","red")[mod$y])
    plot(pearres, col = c("blue","red")[mod$y])
  }
  
  diagnose(mod)
  return(mod)
}

blip_remove_na_refit <- function(time_treat, time_math, removed){
  null_formu <- paste0(paste0("Math_Y", time_math),"~ ",paste0("treat_Y", time_treat))
  blip_treaty1_mathy2.null.mod <- glm(null_formu, data = imputed1)
  
  blip_treaty1_mathy2.full <- creat_blip_1(time_treat, time_math, removed = removed)
  blip_treaty1_mathy2.full.mod <- glm(formula = blip_treaty1_mathy2.full, data = imputed1)
  
  summary(blip_treaty1_mathy2.full.mod)
  coeff <- blip_treaty1_mathy2.full.mod$coefficients
  name_coeff <- names(coeff)
  removed <- name_coeff[is.na(coeff)]
  removed
  
  blip_treaty1_mathy2 <- select_model(blip_treaty1_mathy2.null.mod, blip_treaty1_mathy2.full.mod)
  return(blip_treaty1_mathy2)
}

extract_inter <- function(x){
  pos1 <- gregexpr('_', x )[[1]]
  pos2 <- regexpr(':', x)[[1]]
  
  part1 = NULL
  part2 = NULL
  
  if(length(pos1) > 1){
    part1 = substr(x, 1, pos1[1]+2)
    part2 = substr(x, pos2[1]+1, pos1[2]+2)
  } else {
    part1 = substr(x, 1, pos1[1]+2)
    n = nchar(x)
    part2 = substr(x, pos2[1]+1, n)
  }
  inter <- paste0(part1,":",part2)
  return(inter)
}

na_variables <- function(dat){
  # create list of variables and interactions having na values
  varnames <- colnames(dat)
  na_list = NULL
  
  for (v in varnames){
    if(sum(is.na(dat[,v])) > 0){
      pos <- gregexpr('\\.', v )[[1]]
      if(pos[1] == -1){
        na_list = c(na_list, v)
      } else {
        if(length(pos) == 2){
          inter = paste0(substr(v, 1, pos[1]-1), ":", substr(v, pos[2]+1, length(v)))
          na_list <- c(na_list, inter)
        } else {
          inter = paste0(substr(v, 1, pos[1]-1), ":", substr(v, pos[2]+1, pos[3]-1))
          na_list <- c(na_list, inter)
        }
      }
    } 
  }
  return(na_list)
}

na_singvar <- function(dat){
  # find single variables having na value
  varnames <- colnames(dat)
  na_list = NULL
  
  for (v in varnames){
    if(sum(is.na(dat[,v])) > 0){
      pos <- gregexpr('\\.', v )[[1]]
      if(pos[1] == -1){
        na_list = c(na_list, v)
      } 
    } 
  }
  return(na_list)
}


remove_na_fromfull <- function(time_treat, time_math){
  blip_treaty1_mathy2.full <- creat_blip_1(time_treat, time_math, removed = na_var_list)
  blip_treaty1_mathy2.full.mod <- glm(formula = blip_treaty1_mathy2.full, data = imputed1)
  coeff <- blip_treaty1_mathy2.full.mod$coefficients
  name_coeff <- names(coeff)
  removed <- name_coeff[is.na(coeff)]
  
  return(removed)
}

treat_free_build_mod <- function(treat_time, math_time, removed_var, dat){
  
  treatfree_treaty1_mathy2.full <- creat_treatfree(treat_time, math_time, removed = removed_var)
  treatfree_treaty1_mathy2.full.mod <- glm(formula = treatfree_treaty1_mathy2.full, data = dat)
  
  coeff <- treatfree_treaty1_mathy2.full.mod$coefficients
  name_coeff <- names(coeff)
  removed <- name_coeff[is.na(coeff)]
  
  null.formu <- as.formula(paste0("Math_Y", math_time,"~1"))
  null.mod <- glm(null.formu, data = dat)
  
  treatfree_treaty1_mathy2 <- select_model(null.mod, treatfree_treaty1_mathy2.full.mod)
  plot(treatfree_treaty1_mathy2$fitted.values, treatfree_treaty1_mathy2$residuals)
  return(treatfree_treaty1_mathy2)
  
}


build_treatfree <- function(varlist1, varlist2, bliplist, trtlist){
  n1 = length(varlist1)
  n2 = length(varlist2)
  
  flag1 <- rep(0, n1)
  flag2 <- rep(0, n2)
  
  getse <- function(formu1, formu2){
    # can use for both stage 1 and stage 2
    tflist <- list(formu1, formu2)
    mod <- DTRreg( Math_Y2, bliplist, trtlist, tflist, data = imp, method = "gest", var.estim = "sandwich" )
    se1 = sqrt(mod$covmat[[1]][1])
    se2 = sqrt(mod$covmat[[2]][1])
    return(list(se1, se2))
  }
  
  
  firstformu <- function(varlist, otherformu, one_or_two){
    # can use for both stage 1 and stage 2
    formu <- as.formula(paste0("~ ", varlist[1]))
    if(one_or_two == 1){
      minse = getse(formu, otherformu)[1][[1]][1]  
    } else{
      minse = getse(otherformu, formu)[1][[1]][1]  
    }
    
    n = length(varlist)
    flag <- rep(0, n)
    
    for (i in c(1:n)){
      tmp <- as.formula(paste0("~ ", varlist[i]))
      if(one_or_two == 1){
        se <- getse(tmp, otherformu)  
      } else {
        se <- getse(otherformu, tmp)
      }
      
      if(se[1][[1]][1] <= minse){
        minse = se[1][[1]][1]
        formu <- tmp
        min_i = i
      }
    }
    flag[min_i] = 1
    return(list(formu, flag))
  }
  
  creat_formu_byflag <- function(varlist, varflag){
    # can use for stage 1 and stage 2
    formu = "~"
    n = length(varlist)
    
    for (i in c(1:n)){
      if(varflag[i] == 1){
        if (formu != "~"){
          formu <- paste0(formu, "+", varlist1[i])
        } else {
          formu <- paste0(formu, varlist[i])
        }  
      }
    }
    
    formu <- as.formula(formu)
    return(formu)
  }
  
  formu2 <- as.formula(paste0("~ ", varlist2[1]))
  flag2[1] = 1
  
  init1 <- firstformu(varlist1, formu2, 1)
  formu1 <- init1[[1]]
  flag1 <- init1[[2]]
  
  delse1_thresh = 0.05 # if max(delta_se) < 0.05 then stop
  maxdel1 = 0
  
  build <- function(varlist, otherformu, one_or_two){
    # DO IT TONIGHT
    delse1_thresh = 0.05 # if max(delta_se) < 0.05 then stop
    maxdel1 = 0
    while(maxdel1 == 0 || maxdel1 > delse1_thresh){
      formu1 <- creat_formu_byflag(varlist1, flag1)
      maxvar = 0
      maxse1 = 0
      se = getse(formu1, formu2)
      for (i in c(1:n1)){
        if(flag1[i] == 0){
          tmp <- as.character(formu1)
          tmp <- paste(tmp, collapse = " ")
          tmp <- as.formula(paste0(tmp, "+", varlist1[i]))
          tmpse <- getse(tmp, formu2)
          del = se[1][[1]][1] - tmpse[1][[1]][1]
          if(del > maxdel1){
            maxdel1 = del
            maxvar = i
          }
        }
      }
      if(maxvar > 0){
        flag1[maxvar] = 1
        maxvar = 0  
      }
    }
    return(flag1)
  }
  formu1 <- creat_formu_byflag(varlist1, flag1)
  return(formu1)
}



Hosmer.Lemeshow <- function(fit, NG = 10) {
  # number of groups is ten by default
  
  predict <- predict(fit, type="response")
  if(length(unique(predict)) < 10) {
    stop("number of unique patterns < 10")
  }
  
  # group the predicted probabilities in 10 groups
  groups <- cut(predict,
                breaks=quantile(x=predict,seq(0,1,1/NG)),
                include.lowest=T )
  
  # sum of predicted probabilities by group
  egroups <- tapply(predict, groups, sum)
  
  # sum of observed 1/0 responses by group
  ogroups <- tapply(fit$y, groups, sum) #$
  
  # number of observations per group
  ngroups <- table(groups)
  
  # Hosmer Lemeshow statistic
  statistic <- sum( (ogroups - egroups)^2 /
                      ( egroups*(1-egroups/ngroups)) )
  
  # p-value
  pvalue <- 1-pchisq(statistic, df=(NG-2))
  list(statistic=statistic, df=(NG-2), pvalue=pvalue)
}

logistic_acuracy <- function(mod, imp, year){
  
  prob <- predict(mod, data = imp, type='response')
  imp$prob = prob
  actual_trt <- imp[,paste0("treat_Y",year)]
  pred_trt <- as.factor(sapply(prob, function(x) if (x>0.5) return(1) else return(0)))
  
  confmat <- confusionMatrix(pred_trt, actual_trt)
  precision <- posPredValue(pred_trt, actual_trt, positive="1")
  print("precision")
  print(precision)
  recall <- sensitivity(pred_trt, actual_trt, positive="1")
  print("recall")
  print(recall)
  pred <- prediction(prob, pred_trt)
  F1 <- (2 * precision * recall) / (precision + recall)
  F1
  # ROC area under the curve
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  print("auc")
  print(auc)
}

logistic_envolope <- function(fit, imp){
  
  plothfquant <- function(fit, imp){
    k<-length(predict(fit))
    hnplot <- matrix(0, nrow=k, ncol=19)
    predicted <- predict(fit)
    
    for(i in 1:19){
      set.seed(i)
      new_trt <- rbinom(k, 1, exp(predicted)/(1 + exp(predicted)))
      imp$new_trt <- new_trt
      fitnew <- glm(new_trt ~ TC_Math_Y0 + Math_Y0 + Monthb + Koep_Y1 + Ind_Y0 + 
                      Etn + Anx_Y0 + Rel_Y0 + Schooltype_Y1 + Math_Y0_squared + 
                      Math_Y0_squared:Koep_Y1, data = imp, family=binomial)
      
      hnplot[ ,i] <- sort(abs(rstandard(fitnew)), na.last = TRUE)
    }
    
    hfmin<-apply(hnplot,1,min)
    hfmax<-apply(hnplot,1,max)
    hfquant<-qnorm((k+seq(1,k,by=1))/(2*k+1))
    
    s <- sum(hfquant < hfmax && hfquant > hfmin)
    print(s)
    
    plot(hfquant,sort(abs(rstandard(fit)), na.last = TRUE), xlab="Half-normal quantiles", ylab="Absolute Value of Dev")
    lines(hfquant,hfmin,lty=2)
    lines(hfquant,hfmax,lty=2)
    return(hnplot)
  }
  
  mat <- plothfquant(fit, imp)
  res <-rstandard(fit)
  imp$res <- res
  
  res <-rstandard(fit)
  sum(is.na(res))
  logistic_acuracy(fit, imp, 1)
}

check_ass_ind_dep <- function(inds, dep, imp){
  for (v in inds){
      plot(imp[,dep] ~ imp[,v], main=paste0(c(v, "vs", dep), collapse = " "))
      lines(lowess(imp[,dep] ~ as.numeric(imp[,v])), col=2, lty=2)  
  }
}

