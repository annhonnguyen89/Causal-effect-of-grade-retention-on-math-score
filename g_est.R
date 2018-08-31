# G-estimation
# this file is run after file propensity_wald_test
library(nnet)
library(DataCombine)
library(DTRreg)
library(mice)
library(logistf)
library(heatmapFit)



load("widedf1")
source("allvar.R")
source("g_estimate_func.R")
load("imp_inter_treat_01_02")


############# test preprocess
navars = c("Anx_Y1", "Gap_Y1", "Gap_Y4", "Gap_Y6")
preprocess <- function(imp){
  allvars = colnames(imp)
  retainedvars = setdiff(allvars, navars)
  imp <- imp[,retainedvars]
  
  imp <- within(imp, {
    n1 = as.numeric(treat_Y1)
    if(min(n1) > 0) treat_Y1 = as.factor(n1-1)
    else treat_Y1 = as.factor(n1)
    
    n1 = as.numeric(treat_Y2)
    if(min(n1) > 0) treat_Y2 = as.factor(n1-1)
    else treat_Y2 = as.factor(n1)
    
    n1 = as.numeric(treat_Y3)
    if(min(n1) > 0) treat_Y3 = as.factor(n1-1)
    else treat_Y3 = as.factor(n1)
    
    n1 = as.numeric(treat_Y4)
    if(min(n1) > 0) treat_Y4 = as.factor(n1-1)
    else treat_Y4 = as.factor(n1)
    
    n1 = as.numeric(treat_Y5)
    if(min(n1) > 0) treat_Y5 = as.factor(n1-1)
    else treat_Y5 = as.factor(n1)
    
    n1 = as.numeric(treat_Y6)
    if(min(n1) > 0) treat_Y6 = as.factor(n1-1)
    else treat_Y6 = as.factor(n1)
    
    squared_Math_Y0 <- Math_Y0
    squared_Math_Y1 <- Math_Y1   
    squared_Math_Y2 <- Math_Y2
    squared_Math_Y3 <- Math_Y3
    squared_Math_Y4 <- Math_Y4
    squared_Math_Y5 <- Math_Y5
    squared_Math_Y6 <- Math_Y6
    
    Math_Y0 <- (Math_Y0)^(1/2)
    Math_Y1 <- (Math_Y1)^(1/2)
    Math_Y2 <- (Math_Y2)^(1/2.5)
    Math_Y3 <- (Math_Y3)^(1/2.5)
    Math_Y4 <- (Math_Y4)^(1/2.5)
    Math_Y5 <- (Math_Y5)^(1/2.5)
    Math_Y6 <- (Math_Y6)^(1/2.5)
  })
  
  # remove category Niet gekend in
  imp$Net_Y1 <- as.character(imp$Net_Y1)
  imp <- imp[imp$Net_Y1 != "Niet gekend",]
  imp$Net_Y1 <- as.factor(imp$Net_Y1)
  
  
  imp$Net_Y2 <- as.character(imp$Net_Y2)
  imp <- imp[imp$Net_Y2 != "Niet gekend",]
  imp$Net_Y2 <- as.factor(imp$Net_Y2)
  
  ##########
  
  imp <- imp
  imp$Koep_Y0 <- as.character(imp$Koep_Y0)
  imp[imp$Koep_Y0 == "Ander" | imp$Koep_Y0=="IPCO",]$Koep_Y0 = "Ander"
  imp$Koep_Y0 <- as.factor(imp$Koep_Y0)
  
  imp <- imp
  imp$Koep_Y1 <- as.character(imp$Koep_Y1)
  imp[imp$Koep_Y1 == "Ander" | imp$Koep_Y1=="IPCO" | imp$Koep_Y1 == "VOOP" | imp$Koep_Y1 == "FRSV",]$Koep_Y1 = "Ander"
  imp$Koep_Y1 <- as.factor(imp$Koep_Y1)
  
  imp <- imp
  imp$Koep_Y2 <- as.character(imp$Koep_Y2)
  imp[imp$Koep_Y2 == "Ander" | imp$Koep_Y2=="IPCO" | imp$Koep_Y2 == "VOOP" | imp$Koep_Y2 == "POV" | imp$Koep_Y2 == "FRSV",]$Koep_Y2 = "Ander"
  imp$Koep_Y2 <- as.factor(imp$Koep_Y2)
  
  return(imp)
}

imp1 <- complete(ini,1)
imp2 <- complete(ini,2)
imp3 <- complete(ini,3)

imp <- imp1
imp <- preprocess(imp)
imp2 <- preprocess(imp2)
imp3 <- preprocess(imp3)


################# treatfree model

treatfree_treaty1_mathy2.full <- creat_treatfree(1, 1, removed = navars)
formu <- treatfree_treaty1_mathy2.full[1]
treatfree_treaty1_mathy2.full.mod <- glm(formula = formu, data = imp)
null.mod <- glm(Math_Y1 ~ 1, data = imp)
treatfree_treaty1_mathy2 <- select_model(null.mod, treatfree_treaty1_mathy2.full.mod)
treatfree_treaty1_mathy2$formula
summary(treatfree_treaty1_mathy2)

treatfree_treaty2_mathy2.full <- creat_treatfree(2, 2, removed = navars)
treatfree_treaty2_mathy2.full.mod <- glm(formula = treatfree_treaty2_mathy2.full, data = imp)
null.mod <- glm(Math_Y2 ~ 1, data = imp)
treatfree_treaty2_mathy2 <- select_model(null.mod, treatfree_treaty2_mathy2.full.mod)
summary(treatfree_treaty2_mathy2)

##################### marginal effect

imp$treat_Y1 <- as.numeric(imp$treat_Y1)
imp$treat_Y2 <- as.numeric(imp$treat_Y2)

imp2$treat_Y1 <- as.numeric(imp2$treat_Y1)
imp2$treat_Y2 <- as.numeric(imp2$treat_Y2)

imp3$treat_Y1 <- as.numeric(imp3$treat_Y1)
imp3$treat_Y2 <- as.numeric(imp3$treat_Y2)

trtlist <- list(treat_Y1 ~ Math_Y0 + Rel_Y0 + Monthb + Koep_Y0 + 
                  TC_Lang_Y0 + Etn + Ind_Y0 + PM_Y0 + PercGOK_Y0 + squared_Math_Y0 + 
                  SES + Hyp_Y0 + Math_Y0:Koep_Y0+ Math_Y0:Etn + Rel_Y0:Etn + PercGOK_Y0:Etn + squared_Math_Y0:Etn,
                treat_Y2 ~ Math_Y1 + Math_Y0 + treat_Y1 + Rel_Y1 + 
                  TC_Lang_Y1 + Att_Y1 + squared_Math_Y1 + PSE_Y1 + Hyp_Y1 + 
                  Wel_Y1 + Schooltype_Y1 + Ind_Y1 + Net_Y1 + Peer_Y1 + Etn + 
                  Coop_Y1 + Sex + Math_Y1:Net_Y1 + Math_Y0:Net_Y1 + Rel_Y1:Schooltype_Y1 + 
                  Att_Y1:Net_Y1+squared_Math_Y1:Net_Y1 +
                  Schooltype_Y1:Ind_Y1 + Net_Y1:Peer_Y1)

treatfreelist <- list(~ Math_Y0 + TC_Math_Y0 + PSE_Y0 + PM_Y0 + Sex + PercGOK_Y1 + 
                  Ind_Y0 + Gap_Y0 + Hyp_Y0 + Koep_Y1 + GOK + Aso_Y0 + Att_Y0 + 
                  Peer_Y0 + Self_Y0 + Rel_Y0 + squared_Math_Y0,
                ~ Math_Y1 + TC_Math_Y1 + Rel_Y1 + PercGOK_Y1 + Sex + 
                  PSE_Y1 + treat_Y1 + Ind_Y1 + Koep_Y2 + Hlang + GOK + Hyp_Y1 + 
                  Schooltype_Y1 + Etn + Self_Y1 + Att_Y1 + Pro_Y1 + Agr_Y1 + 
                  TC_Lang_Y1 + squared_Math_Y1)

bliplist <- list(~1, ~1)
mod1 <- DTRreg( Math_Y2, bliplist, trtlist, treatfreelist, data = imp, method = "gest", var.estim = "sandwich" )
#Error in solve.default(dT.dalpha) : 
#system is computationally singular: reciprocal condition number = 1.37602e-16
#--> simplify treat model by remove Net_Y1:Peer_Y1 from treatment-free model 2

trtlist <- list(treat_Y1 ~ Math_Y0 + Rel_Y0 + Monthb + Koep_Y0 + 
                  TC_Lang_Y0 + Etn + Ind_Y0 + PM_Y0 + PercGOK_Y0 + squared_Math_Y0 + 
                  SES + Hyp_Y0+ Math_Y0:Koep_Y0+ Math_Y0:Etn + Rel_Y0:Etn + PercGOK_Y0:Etn + squared_Math_Y0:Etn,
                treat_Y2 ~ Math_Y1 + Math_Y0 + treat_Y1 + Rel_Y1 + 
                  TC_Lang_Y1 + Att_Y1 + squared_Math_Y1 + PSE_Y1 + Hyp_Y1 + 
                  Wel_Y1 + Schooltype_Y1 + Ind_Y1 + Net_Y1 + Peer_Y1 + Etn + 
                  Coop_Y1 + Sex + Math_Y1:Net_Y1 + Math_Y0:Net_Y1 + Rel_Y1:Schooltype_Y1 + 
                  Att_Y1:Net_Y1+squared_Math_Y1:Net_Y1 +
                  Schooltype_Y1:Ind_Y1)

mod1 <- DTRreg( Math_Y2, bliplist, trtlist, treatfreelist, data = imp, method = "gest", var.estim = "sandwich" )
mod1

mod2 <- DTRreg( Math_Y2, bliplist, trtlist, treatfreelist, data = imp2, method = "gest", var.estim = "sandwich" )
mod2

mod3 <- DTRreg( Math_Y2, bliplist, trtlist, treatfreelist, data = imp3, method = "gest", var.estim = "sandwich" )
mod3

### pooling
pool_analysis <- function(mod1,  mod2, mod3, varname){
  psi1 <- (mod1$psi[[1]] + mod2$psi[[1]] + mod3$psi[[1]])/3
  psi2 <- (mod1$psi[[2]] + mod2$psi[[2]] + mod3$psi[[2]])/3
  n = length(psi1)
  var1 = 0*c(1:n)
  var2 = 0*c(1:n)
  se1 = 0*c(1:n)
  se2 = 0*c(1:n)
  low1 = 0*c(1:n)
  low2 = 0*c(1:n)
  
  up1 = 0*c(1:n)
  up2 = 0*c(1:n)
  names(psi1) = names(psi2) = names(se1) = names(se2) = names(low1) = names(up1) = names(low2) = names(up2)= varname
  
  for (i in c(1:n)){
    var1 = mean(c(c(mod1$covmat[[1]][i,i], mod2$covmat[[1]][i,i], mod3$covmat[[1]][i,i]))) +
      (1+1/3)/2*((mod1$psi[[1]][i] - psi1[i])^2 + (mod2$psi[[1]][i] - psi1[i])^2 + (mod3$psi[[1]][i] - psi1[i])^2 )
    
    var2 = mean(c(c(mod1$covmat[[2]][i,i], mod2$covmat[[2]][i,i], mod3$covmat[[2]][i,i]))) +
      (1+1/3)/2*((mod1$psi[[2]][i] - psi2[i])^2 + (mod2$psi[[2]][i] - psi2[i])^2 + (mod3$psi[[2]][i] - psi2[i])^2 )
    
    se1[i] = sqrt(var1)
    se2[i] = sqrt(var2)
    
  }
  low1 = psi1 - 1.96*se1/sqrt(3)
  up1 = psi1 + 1.96*se1/sqrt(3)
  
  low2 = psi2 - 1.96*se2/sqrt(3)
  up2 = psi2 + 1.96*se2/sqrt(3)
  
  df = data.frame(year1 = psi1, year2 = psi2, se1 = se1, se2 = se2, low1 = low1, up1 = up1, low2 = low2, up2 = up2)
  
  return(df)
}

margin_result = pool_analysis(mod1, mod2, mod3, c("int"))
margin_result
mod1

####################### causal effect conditioning on previous math score
center_math <- function(dat){
  dat <- within(dat,{
    Math_Y0 <- Math_Y0 - mean(Math_Y0)
    Math_Y1 <- Math_Y1 - mean(Math_Y1)
    Math_Y2 <- Math_Y2 - mean(Math_Y2)
  })
  return(dat)
}

imp <- center_math(imp)
imp2 <- center_math(imp2)
imp3 <- center_math(imp3)

blip_math_list <- list(~ 1 + Math_Y0, ~ 1 + Math_Y1)
mod1 <- DTRreg( Math_Y2, blip_math_list, trtlist, treatfreelist, data = imp, method = "gest", var.estim = "sandwich" )
mod2 <- DTRreg( Math_Y2, blip_math_list, trtlist, treatfreelist, data = imp2, method = "gest", var.estim = "sandwich" )
mod3 <- DTRreg( Math_Y2, blip_math_list, trtlist, treatfreelist, data = imp3, method = "gest", var.estim = "sandwich" )
mod1

## error mod3 --> simplify trt model by removing squared_Math_Y1:Net_Y1 + Schooltype_Y1:Ind_Y1 in second tr
trtlist <- list(treat_Y1 ~ Math_Y0 + Rel_Y0 + Monthb + Koep_Y0 + 
                  TC_Lang_Y0 + Etn + Ind_Y0 + PM_Y0 + PercGOK_Y0 + squared_Math_Y0 + 
                  SES + Hyp_Y0 + Math_Y0:Koep_Y0+ Math_Y0:Etn + Rel_Y0:Etn + PercGOK_Y0:Etn + squared_Math_Y0:Etn,
                treat_Y2 ~ Math_Y1 + Math_Y0 + treat_Y1 + Rel_Y1 + 
                  TC_Lang_Y1 + Att_Y1 + squared_Math_Y1 + PSE_Y1 + Hyp_Y1 + 
                  Wel_Y1 + Schooltype_Y1 + Ind_Y1 + Net_Y1 + Peer_Y1 + Etn + 
                  Coop_Y1 + Sex + Math_Y1:Net_Y1 + Math_Y0:Net_Y1 + Rel_Y1:Schooltype_Y1 + 
                  Att_Y1:Net_Y1
                  )
mod1 <- DTRreg( Math_Y2, blip_math_list, trtlist, treatfreelist, data = imp, method = "gest", var.estim = "sandwich" )
mod2 <- DTRreg( Math_Y2, blip_math_list, trtlist, treatfreelist, data = imp2, method = "gest", var.estim = "sandwich" )
mod3 <- DTRreg( Math_Y2, blip_math_list, trtlist, treatfreelist, data = imp3, method = "gest", var.estim = "sandwich" )
mod1
mod2
mod3
pool_analysis(mod1, mod2, mod3, c("int", "prev_math"))

######################## effect conditioning on prev umbrella
preprocess_koep <- function(imp){
  levels(imp$Koep_Y0)
  imp$FOPEM0 = 0
  imp$FRSV0 = 0
  imp$GO0 = 0
  imp$OVSG0 = 0
  imp$VSKO0 = 0
  imp$otherKoep0 = 0
  
  imp[imp$Koep_Y0 == "Ander",]$otherKoep0 = 1
  
  imp[imp$Koep_Y0 == "FOPEM",]$FOPEM0 = 1
  imp[imp$Koep_Y0 == "FRSV",]$FRSV0 = 1
  imp[imp$Koep_Y0 == "GO!",]$GO0 = 1
  imp[imp$Koep_Y0 == "OVSG",]$OVSG0 = 1
  imp[imp$Koep_Y0 == "VSKO",]$VSKO0 = 1
  
  levels(imp$Koep_Y1)
  
  table(imp$Koep_Y1, imp$treat_Y1)
  
  imp$FOPEM1 = 0
  imp$GO1 = 0
  imp$OVSG1 = 0
  imp$VSKO1 = 0
  imp$otherKoep1 = 0
  
  imp[imp$Koep_Y1 == "FOPEM",]$FOPEM1 = 1
  imp[imp$Koep_Y1 == "GO!",]$GO1 = 1
  imp[imp$Koep_Y1 == "OVSG",]$OVSG1 = 1
  imp[imp$Koep_Y1 == "VSKO",]$VSKO1 = 1
  imp[imp$Koep_Y1 == "Ander",]$otherKoep1 = 1
  
  return(imp)
}
imp <- preprocess_koep(imp)
imp2 <- preprocess_koep(imp2)
imp3 <- preprocess_koep(imp3)

blip_koep_list <- list(~ FOPEM0 + FRSV0 + GO0 + OVSG0 + VSKO0 + 1, 
                            ~ FOPEM1 + GO1 + OVSG1 + VSKO1 + 1  )

mod1 <- DTRreg( Math_Y2, blip_koep_list, trtlist, treatfreelist, data = imp, method = "gest", var.estim = "sandwich" )
mod2 <- DTRreg( Math_Y2, blip_koep_list, trtlist, treatfreelist, data = imp2, method = "gest", var.estim = "sandwich" )
mod3 <- DTRreg( Math_Y2, blip_koep_list, trtlist, treatfreelist, data = imp3, method = "gest", var.estim = "sandwich" )


pool_koep_effect <- function(mod1, mod2, mod3, varnames1, varnames2){
  psi1 <- (mod1$psi[[1]] + mod2$psi[[1]] + mod3$psi[[1]])/3
  psi2 <- (mod1$psi[[2]] + mod2$psi[[2]] + mod3$psi[[2]])/3
  n1 = length(psi1)
  n2 = length(psi2)
  
  m1 = psi1
  m2 = psi2
  #browser()
  for (i in c(2:n1)){
    m1[i] = psi1[1] + psi1[i]
  }
  
  for (i in c(2:n2)){
    m2[i] = psi2[1] + psi2[i]
  }
  
  var1 = 0*c(1:n1)
  var2 = 0*c(1:n2)
  se1 = 0*c(1:n1)
  se2 = 0*c(1:n2)
  low1 = 0*c(1:n1)
  low2 = 0*c(1:n2)
  
  up1 = 0*c(1:n1)
  up2 = 0*c(1:n2)
  names(m1) = names(se1) = names(low1) = names(up1) = varnames1
  names(m2) = names(se2) = names(low2) = names(up2) = varnames2
  
  for (i in c(1:n1)){
    var1[i] = mean(c(c(mod1$covmat[[1]][i,i], mod2$covmat[[1]][i,i], mod3$covmat[[1]][i,i]))) +
      (1+1/3)/2*((mod1$psi[[1]][i] - psi1[i])^2 + (mod2$psi[[1]][i] - psi1[i])^2 + (mod3$psi[[1]][i] - psi1[i])^2 )
  }
  for (i in c(1:n2)){
    var2[i] = mean(c(c(mod1$covmat[[2]][i,i], mod2$covmat[[2]][i,i], mod3$covmat[[2]][i,i]))) +
      (1+1/3)/2*((mod1$psi[[2]][i] - psi2[i])^2 + (mod2$psi[[2]][i] - psi2[i])^2 + (mod3$psi[[2]][i] - psi2[i])^2 )
  }
  
  se1[1] = sqrt(var1[1])
  se2[1] = sqrt(var2[1])
  
  low1[1] = m1[1] - 1.96*se1[1]/sqrt(3)
  up1[1] = m1[1] + 1.96*se1[1]/sqrt(3)
  
  low2[1] = m2[1] - 1.96*se2[1]/sqrt(3)
  up2[1] = m2[1] + 1.96*se2[1]/sqrt(3)
  

  k = n1
  for (i in c(2:k)){
    m = i
    se1[m] = sqrt((var1[1] + var1[i])/2) #+ 2*(mod1$covmat[[1]][i,1] + mod2$covmat[[1]][i,1] + mod3$covmat[[1]][i,1])/3)
    
    low1[m] = m1[m] - 1.96*se1[m]/sqrt(3)
    up1[m] = m1[m]  + 1.96*se1[m]/sqrt(3)
  }
  
  k = n2
  for (i in c(2:k)){
    m = i
    se2[m] = sqrt((var2[1] + var2[i])/2)  #+ 2*(mod1$covmat[[2]][i,1] + mod2$covmat[[2]][i,1] + mod3$covmat[[2]][i,1])/3)
    
    low2[m] = m2[m] - 1.96*se2[m]/sqrt(3)
    up2[m] = m2[m] + 1.96*se2[m]/sqrt(3)   
  }
  #browser()
  df1 = data.frame(m1 = m1, se1 = se1, low1 = low1, up1 = up1)
  df2 = data.frame(m2 = m2, se2 = se2, low2 = low2, up2 = up2)
  return(list(df1, df2))
}
mod3
names1 = c("FOPEM0", "FRSV0", "GO0", "OVSG0", "VSKO0", "other0")
names2= c("FOPEM1", "GO1", "OVSG1", "VSKO1", "other1")
pool_koep_effect(mod1, mod2, mod3, names1, names2)



