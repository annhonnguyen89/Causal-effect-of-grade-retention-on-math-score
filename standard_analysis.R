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

navars = c("Anx_Y1", "Gap_Y1", "Gap_Y4", "Gap_Y6")
imp1 <- complete(ini,1)
imp2 <- complete(ini,2)
imp3 <- complete(ini,3)

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
  imp[imp$Koep_Y1 == "Ander" | imp$Koep_Y1=="IPCO" | imp$Koep_Y1 == "VOOP",]$Koep_Y1 = "Ander"
  imp$Koep_Y1 <- as.factor(imp$Koep_Y1)
  
  imp <- imp
  imp$Koep_Y2 <- as.character(imp$Koep_Y2)
  imp[imp$Koep_Y2 == "Ander" | imp$Koep_Y2=="IPCO" | imp$Koep_Y2 == "VOOP",]$Koep_Y2 = "Ander"
  imp$Koep_Y2 <- as.factor(imp$Koep_Y2)
  
  imp$Math_Y0 <- imp$Math_Y0 - mean(imp$Math_Y0)
  imp$Math_Y1 <- imp$Math_Y1 - mean(imp$Math_Y1)
  imp$Math_Y2 <- imp$Math_Y2 - mean(imp$Math_Y2)  
  
  return(imp)
}

imp <- imp1
imp <- preprocess(imp)
imp2 <- preprocess(imp2)
imp3 <- preprocess(imp3)

##################################### adjust for covariates before retention in year 1

full.formu = "Math_Y2 ~ Sex+Monthb+Hlang+SES+GOK+Etn+treat_Y1+treat_Y2 + 
PM_Y0+Koep_Y0+Net_Y0+Schooltype_Y0+PercGOK_Y0+Wel_Y0+Att_Y0+Aso_Y0+Coop_Y0+
Ind_Y0+Peer_Y0+Pro_Y0+Hyp_Y0+Agr_Y0+Self_Y0+SuppHE_Y0+
Rel_Y0+PSE_Y0+TC_Lang_Y0+Koep_Y0+PM_Y0+Schooltype_Y0+Net_Y0+PercGOK_Y0 + Math_Y0 + squared_Math_Y0"



full.mod <- lm(as.formula(full.formu), data = imp)
null.mod <- lm(Math_Y2 ~ treat_Y1 + treat_Y2 + Math_Y0, data = imp)
mod <- step(null.mod, scope = list(lower = null.mod, upper = full.mod), direction = "both", trace = 1)


formu = "Math_Y2 ~ treat_Y1 + treat_Y2 + Math_Y0 + PSE_Y0 + 
    Sex + Att_Y0 + PM_Y0 + PercGOK_Y0 + squared_Math_Y0 + Rel_Y0 + 
    Monthb + GOK + Etn + Aso_Y0 + Self_Y0 + Hyp_Y0 + Schooltype_Y0 + 
    Coop_Y0"

cont_list = c("Math_Y0","PSE_Y0", "Att_Y0","PercGOK_Y0", "squared_Math_Y0", 
              "Rel_Y0", "Monthb", "Aso_Y0", "Self_Y0", "Hyp_Y0", "Coop_Y0")

cat_list = c("Sex", "GOK", "Schooltype_Y0", "Etn", "PM_Y0")


##### creat model
removedvar = "PM_Y0"

  create_inter <- function(cont_list, cat_list){
    interlist = c()
    for (cont in cont_list)
      for (cat in cat_list){
        term = paste0(cont, ":", cat)
        interlist = append(interlist, term)
      }
    return(interlist)
  }
  
  cat_list = setdiff(cat_list, removedvar)
  inter_list = create_inter(cont_list, cat_list)
  #inter_list = setdiff(inter_list, removedvar)
  
  full_formu = paste0(c(cont_list, cat_list, inter_list), collapse = "+")
  full_formu = paste0("Math_Y2 ~ ", full_formu)
  
  full_formu = as.formula(full_formu)
  null.mod <- lm(Math_Y2 ~ treat_Y1 + treat_Y2 + Math_Y0, data = imp)
  full <- lm(Math_Y2 ~ treat_Y1 + treat_Y2 + Math_Y0 + PSE_Y0 + Att_Y0 + PercGOK_Y0 + squared_Math_Y0 + 
                   Rel_Y0 + Monthb + Aso_Y0 + Self_Y0 + Hyp_Y0 + Coop_Y0 + Sex + 
                   GOK + Schooltype_Y0 + Etn + Math_Y0:Sex + Math_Y0:GOK + Math_Y0:Schooltype_Y0 + 
                   Math_Y0:Etn + PSE_Y0:Sex + PSE_Y0:GOK + PSE_Y0:Schooltype_Y0 + 
                   PSE_Y0:Etn + Att_Y0:Sex + Att_Y0:GOK + Att_Y0:Schooltype_Y0 + 
                   Att_Y0:Etn + PercGOK_Y0:Sex + PercGOK_Y0:GOK + PercGOK_Y0:Schooltype_Y0 + 
                   PercGOK_Y0:Etn + squared_Math_Y0:Sex + squared_Math_Y0:GOK + 
                   squared_Math_Y0:Schooltype_Y0 + squared_Math_Y0:Etn + Rel_Y0:Sex + 
                   Rel_Y0:GOK + Rel_Y0:Schooltype_Y0 + Rel_Y0:Etn + Monthb:Sex + 
                   Monthb:GOK + Monthb:Schooltype_Y0 + Monthb:Etn + Aso_Y0:Sex + 
                   Aso_Y0:GOK + Aso_Y0:Schooltype_Y0 + Aso_Y0:Etn + Self_Y0:Sex + 
                   Self_Y0:GOK + Self_Y0:Schooltype_Y0 + Self_Y0:Etn + Hyp_Y0:Sex + 
                   Hyp_Y0:GOK + Hyp_Y0:Schooltype_Y0 + Hyp_Y0:Etn + Coop_Y0:Sex + 
                   Coop_Y0:GOK + Coop_Y0:Schooltype_Y0 + Coop_Y0:Etn, data = imp)
  
  selected_mod <- step(null.mod, scope = list(lower = null.mod, upper = full), direction = "both", trace = 1)

selected_mod$call
selected_mod_trt1 <- lm(Math_Y2 ~ treat_Y1 + treat_Y2 + Math_Y0 + PSE_Y0 + 
                          Sex + Att_Y0 + PercGOK_Y0 + squared_Math_Y0 + Rel_Y0 + Monthb + 
                          Etn + GOK + Aso_Y0 + Self_Y0 + Hyp_Y0 + Coop_Y0 + Math_Y0:Sex + 
                          Sex:squared_Math_Y0 + Etn:Aso_Y0 + Sex:Aso_Y0 + Monthb:GOK + 
                          squared_Math_Y0:Etn + Math_Y0:Etn + PercGOK_Y0:GOK, data = imp)
# because data of (treat_Y1, treat_Y2): (1,0) = 341 rows, (1,1) = 0, (0,0) = 4787 rows, (0,1) = 486 rows ==> do not add interaction treat_Y1*treat_Y2

summary(selected_mod_trt1)

#y1 + y2 + y1*y2

#1,0 y1 
#0,0 0
#0,1 y2
# 1, 1 y1 + y2 + y1*y2 ==> cong huong neu ca 2 cung xay ra

# G-est: 1,0; 0,1 = 1,1 = 
selected_mod_trt1$call
beforetrt1_mod1 <- lm(Math_Y2 ~ treat_Y1 + treat_Y2 + Math_Y0 + PSE_Y0 + 
                        Sex + Att_Y0 + PercGOK_Y0 + squared_Math_Y0 + Rel_Y0 + Monthb + 
                        Etn + GOK + Aso_Y0 + Self_Y0 + Hyp_Y0 + Coop_Y0 + Math_Y0:Sex + 
                        Sex:squared_Math_Y0 + Etn:Aso_Y0 + Sex:Aso_Y0 + Monthb:GOK + 
                        squared_Math_Y0:Etn + Math_Y0:Etn + PercGOK_Y0:GOK, data = imp)
beforetrt1_mod2 <- lm(Math_Y2 ~ treat_Y1 + treat_Y2 + Math_Y0 + PSE_Y0 + 
                        Sex + Att_Y0 + PercGOK_Y0 + squared_Math_Y0 + Rel_Y0 + Monthb + 
                        Etn + GOK + Aso_Y0 + Self_Y0 + Hyp_Y0 + Coop_Y0 + Math_Y0:Sex + 
                        Sex:squared_Math_Y0 + Etn:Aso_Y0 + Sex:Aso_Y0 + Monthb:GOK + 
                        squared_Math_Y0:Etn + Math_Y0:Etn + PercGOK_Y0:GOK, data = imp2)

beforetrt1_mod3 <- lm(Math_Y2 ~ treat_Y1 + treat_Y2 + Math_Y0 + PSE_Y0 + 
                        Sex + Att_Y0 + PercGOK_Y0 + squared_Math_Y0 + Rel_Y0 + Monthb + 
                        Etn + GOK + Aso_Y0 + Self_Y0 + Hyp_Y0 + Coop_Y0 + Math_Y0:Sex + 
                        Sex:squared_Math_Y0 + Etn:Aso_Y0 + Sex:Aso_Y0 + Monthb:GOK + 
                        squared_Math_Y0:Etn + Math_Y0:Etn + PercGOK_Y0:GOK, data = imp3)
# pool analysis
pool_analysis <- function(mod1, mod2, mod3){
  varnames = names(mod1$coefficients)
  
  m = (mod1$coefficients + mod2$coefficients + mod3$coefficients)/3
  n  = length(varnames)
  se = 0*c(1:n)
  low = 0*c(1:n)
  high = 0*c(1:n)
  names(m) = varnames
  names(se) = varnames
  names(high) = varnames
  names(low) = varnames
  for (i in varnames){
    se[i] = sqrt((vcov(mod1)[i,i] + vcov(mod2)[i,i] + vcov(mod3)[i,i])/3)
    low[i] = m[i]- 1.96*se[i]/sqrt(3)
    high[i] = m[i] + 1.96*se[i]/sqrt(3)  
  }
  df = data.frame(m = m, se = se,low = low, high = high)
  return(df)
}
pool_beforetrt1 <- pool_analysis(beforetrt1_mod1, beforetrt1_mod2, beforetrt1_mod3)
pool_beforetrt1
##################################### adjust for all covariates before retention Y2

full.formu = "Math_Y2 ~ Sex+Monthb+Hlang+SES+GOK+Etn+treat_Y1+treat_Y2 + 
PM_Y0+Koep_Y0+Net_Y0+Schooltype_Y0+PercGOK_Y0+Wel_Y0+Att_Y0+Aso_Y0+Coop_Y0+
Ind_Y0+Peer_Y0+Pro_Y0+Hyp_Y0+Agr_Y0+Self_Y0+SuppHE_Y0+
Rel_Y0+PSE_Y0+TC_Lang_Y0+Koep_Y0+PM_Y0+Schooltype_Y0+Net_Y0+PercGOK_Y0 + Math_Y0 + squared_Math_Y0 +
PM_Y1+Koep_Y1+Net_Y1+Schooltype_Y1+PercGOK_Y1+Wel_Y1+Att_Y1+Aso_Y1+Coop_Y1+
Ind_Y1+Peer_Y1+Pro_Y1+Hyp_Y1+Agr_Y1+Self_Y1+SuppHE_Y1+
Rel_Y1+PSE_Y1+TC_Lang_Y1+Koep_Y1+PM_Y1+Schooltype_Y1+Net_Y1+PercGOK_Y1 + Math_Y1 + squared_Math_Y1"

full.mod <- lm(as.formula(full.formu), data = imp)
null.mod <- lm(Math_Y2 ~ treat_Y1 + treat_Y2 + Math_Y0, data = imp)
mod <- step(null.mod, scope = list(lower = null.mod, upper = full.mod), direction = "both", trace = 1)
mod$call
cat_list = c("PM_Y0", "Koep_Y1", "Etn", "GOK")
cont_list = c("Math_Y0", "squared_Math_Y1", "Rel_Y1", "PSE_Y0","PercGOK_Y1","Ind_Y1", "Att_Y1", 
              "Monthb","Att_Y0","PSE_Y1", "Pro_Y1", "Agr_Y1", "Agr_Y0", "Self_Y0", "SuppHE_Y0", "TC_Lang_Y1",
              "Hyp_Y1", "Aso_Y0", "Self_Y1", "Math_Y1")

create_inter <- function(cont_list, cat_list){
  interlist = c()
  for (cont in cont_list)
    for (cat in cat_list){
      term = paste0(cont, ":", cat)
      interlist = append(interlist, term)
    }
  return(interlist)
}


inter_list = create_inter(cont_list, cat_list)
#inter_list = setdiff(inter_list, removedvar)

full_formu = paste0(c(cont_list, cat_list, inter_list), collapse = "+")
full_formu = paste0("Math_Y2 ~ ", full_formu)

null.mod <- lm(Math_Y2 ~ treat_Y1 + treat_Y2 + Math_Y0 + Math_Y1, data = imp)
full <- lm(Math_Y2 ~ treat_Y1 + treat_Y2 + Math_Y0+squared_Math_Y1+Rel_Y1+PSE_Y0+PercGOK_Y1+
             Ind_Y1+Att_Y1+Monthb+Att_Y0+PSE_Y1+Pro_Y1+Agr_Y1+Agr_Y0+Self_Y0+SuppHE_Y0+
             TC_Lang_Y1+Hyp_Y1+Aso_Y0+Self_Y1+Math_Y1+PM_Y0+Koep_Y1+Etn+GOK+Math_Y0:PM_Y0+
             Math_Y0:Koep_Y1+Math_Y0:Etn+Math_Y0:GOK+squared_Math_Y1:PM_Y0+squared_Math_Y1:Koep_Y1+
             squared_Math_Y1:Etn+squared_Math_Y1:GOK+Rel_Y1:PM_Y0+Rel_Y1:Koep_Y1+Rel_Y1:Etn+Rel_Y1:GOK+
             PSE_Y0:PM_Y0+PSE_Y0:Koep_Y1+PSE_Y0:Etn+PSE_Y0:GOK+PercGOK_Y1:PM_Y0+PercGOK_Y1:Koep_Y1+
             PercGOK_Y1:Etn+PercGOK_Y1:GOK+Ind_Y1:PM_Y0+Ind_Y1:Koep_Y1+Ind_Y1:Etn+Ind_Y1:GOK+Att_Y1:PM_Y0+
             Att_Y1:Koep_Y1+Att_Y1:Etn+Att_Y1:GOK+Monthb:PM_Y0+Monthb:Koep_Y1+Monthb:Etn+Monthb:GOK+Att_Y0:PM_Y0+
             Att_Y0:Koep_Y1+Att_Y0:Etn+Att_Y0:GOK+PSE_Y1:PM_Y0+PSE_Y1:Koep_Y1+PSE_Y1:Etn+PSE_Y1:GOK+Pro_Y1:PM_Y0+
             Pro_Y1:Koep_Y1+Pro_Y1:Etn+Pro_Y1:GOK+Agr_Y1:PM_Y0+Agr_Y1:Koep_Y1+Agr_Y1:Etn+Agr_Y1:GOK+Agr_Y0:PM_Y0+
             Agr_Y0:Koep_Y1+Agr_Y0:Etn+Agr_Y0:GOK+Self_Y0:PM_Y0+Self_Y0:Koep_Y1+Self_Y0:Etn+Self_Y0:GOK+SuppHE_Y0:PM_Y0+
             SuppHE_Y0:Koep_Y1+SuppHE_Y0:Etn+SuppHE_Y0:GOK+TC_Lang_Y1:PM_Y0+TC_Lang_Y1:Koep_Y1+TC_Lang_Y1:Etn+TC_Lang_Y1:GOK+
             Hyp_Y1:PM_Y0+Hyp_Y1:Koep_Y1+Hyp_Y1:Etn+Hyp_Y1:GOK+Aso_Y0:PM_Y0+Aso_Y0:Koep_Y1+Aso_Y0:Etn+Aso_Y0:GOK+Self_Y1:PM_Y0+
             Self_Y1:Koep_Y1+Self_Y1:Etn+Self_Y1:GOK+Math_Y1:PM_Y0+Math_Y1:Koep_Y1+Math_Y1:Etn+Math_Y1:GOK, data = imp)

selected_mod <- step(null.mod, scope = list(lower = null.mod, upper = full), direction = "both", trace = 1)
selected_mod$call

beforetrt2_mod1 <- lm(Math_Y2 ~ treat_Y1 + treat_Y2 + Math_Y0 + Math_Y1 + 
                        Rel_Y1 + Ind_Y1 + PercGOK_Y1 + squared_Math_Y1 + PSE_Y0 + 
                        PM_Y0 + Pro_Y1 + Agr_Y0 + Monthb + Att_Y1 + Etn + GOK + Self_Y0 + 
                        Agr_Y1 + PSE_Y1 + Koep_Y1 + Att_Y0 + Self_Y1 + PM_Y0:Monthb + 
                        Pro_Y1:GOK + Pro_Y1:Koep_Y1 + PSE_Y0:Koep_Y1 + Etn:Self_Y0 + 
                        PM_Y0:Att_Y1, data = imp)
beforetrt2_mod2 <- lm(Math_Y2 ~ treat_Y1 + treat_Y2 + Math_Y0 + Math_Y1 + 
                        Rel_Y1 + Ind_Y1 + PercGOK_Y1 + squared_Math_Y1 + PSE_Y0 + 
                        PM_Y0 + Pro_Y1 + Agr_Y0 + Monthb + Att_Y1 + Etn + GOK + Self_Y0 + 
                        Agr_Y1 + PSE_Y1 + Koep_Y1 + Att_Y0 + Self_Y1 + PM_Y0:Monthb + 
                        Pro_Y1:GOK + Pro_Y1:Koep_Y1 + PSE_Y0:Koep_Y1 + Etn:Self_Y0 + 
                        PM_Y0:Att_Y1, data = imp2)
beforetrt2_mod3 <- lm(Math_Y2 ~ treat_Y1 + treat_Y2 + Math_Y0 + Math_Y1 + 
                        Rel_Y1 + Ind_Y1 + PercGOK_Y1 + squared_Math_Y1 + PSE_Y0 + 
                        PM_Y0 + Pro_Y1 + Agr_Y0 + Monthb + Att_Y1 + Etn + GOK + Self_Y0 + 
                        Agr_Y1 + PSE_Y1 + Koep_Y1 + Att_Y0 + Self_Y1 + PM_Y0:Monthb + 
                        Pro_Y1:GOK + Pro_Y1:Koep_Y1 + PSE_Y0:Koep_Y1 + Etn:Self_Y0 + 
                        PM_Y0:Att_Y1, data = imp3)

pool_beforetrt2 <- pool_analysis(beforetrt2_mod1, beforetrt2_mod2, beforetrt2_mod3)
pool_beforetrt2





