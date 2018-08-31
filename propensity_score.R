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

########### build logistic regression treat_Y1

null.mod <- glm(treat_Y1 ~ Math_Y0, data = imp, family = binomial)
full.mod <- glm(treat_Y1 ~ Sex+Monthb+Hlang+SES+GOK+Etn+
                  PM_Y0+Koep_Y0+Net_Y0+Schooltype_Y0+PercGOK_Y0+Wel_Y0+Att_Y0+Aso_Y0+Coop_Y0+
                  Ind_Y0+Peer_Y0+Pro_Y0+Hyp_Y0+Agr_Y0+Self_Y0+SuppHE_Y0+
                  Rel_Y0+PSE_Y0+TC_Lang_Y0+Koep_Y0+PM_Y0+Schooltype_Y0+Net_Y0+
                  PercGOK_Y0 + Math_Y0 + squared_Math_Y0, data = imp, family = binomial)
mod <- step(null.mod, scope = list(lower = null.mod, upper = full.mod), direction = "both", trace = 1)
mod$call
treat_Y1 ~ Math_Y0 + Rel_Y0 + Monthb + Koep_Y0 + 
  TC_Lang_Y0 + Etn + Ind_Y0 + PM_Y0 + PercGOK_Y0 + squared_Math_Y0 + 
  SES + Hyp_Y0

cont_list = c("Math_Y0", "Rel_Y0", "Monthb", "TC_Lang_Y0", "Ind_Y0", "PercGOK_Y0", "squared_Math_Y0", "Hyp_Y0", "SES")
cat_list = c("Koep_Y0", "Etn")

create_inter <- function(cont_list, cat_list){
  interlist = c()
  for (cont in cont_list)
    for (cat in cat_list){
      #browser()
      term = paste0(cont, ":", cat)
      interlist = append(interlist, term)
    }
  return(interlist)
}

inter_list = create_inter(cont_list, cat_list)
paste0(inter_list, collapse = "+")
newmod <- glm(treat_Y1 ~ Math_Y0 + Rel_Y0 + Monthb + Koep_Y0 + 
               TC_Lang_Y0 + Etn + Ind_Y0 + PM_Y0 + PercGOK_Y0 + squared_Math_Y0 + 
               SES + Hyp_Y0 + Math_Y0:PM_Y0+
                Math_Y0:Koep_Y0+Math_Y0:Etn+Rel_Y0:Koep_Y0+Rel_Y0:Etn+Monthb:Koep_Y0+Monthb:Etn+TC_Lang_Y0:Koep_Y0+
                TC_Lang_Y0:Etn+Ind_Y0:Koep_Y0+Ind_Y0:Etn+PercGOK_Y0:Koep_Y0+PercGOK_Y0:Etn+squared_Math_Y0:Koep_Y0+
                squared_Math_Y0:Etn+Hyp_Y0:Koep_Y0+Hyp_Y0:Etn+SES:Koep_Y0+SES:Etn, data = imp, family = binomial)

summary(newmod)
library(car)
Anova(newmod, type = "II")

newmod <- glm(treat_Y1 ~ Math_Y0 + Rel_Y0 + Monthb + Koep_Y0 + 
                TC_Lang_Y0 + Etn + Ind_Y0 + PM_Y0 + PercGOK_Y0 + squared_Math_Y0 + 
                SES + Hyp_Y0 + Math_Y0:Koep_Y0+ Math_Y0:Etn + Rel_Y0:Etn + PercGOK_Y0:Etn + squared_Math_Y0:Etn, data = imp, family = binomial)
Anova(newmod, type = "II")
summary(newmod)

########### check balance
library(MatchIt)
library(dplyr)
library(ggplot2)
mod_match <- matchit(treat_Y1 ~ Math_Y0 + Rel_Y0 + Monthb + Koep_Y0 + 
                       TC_Lang_Y0 + Etn + Ind_Y0 + PM_Y0 + PercGOK_Y0 + squared_Math_Y0 + 
                       SES + Hyp_Y0 + Math_Y0:Koep_Y0+ Math_Y0:Etn + Rel_Y0:Etn + PercGOK_Y0:Etn + squared_Math_Y0:Etn,
                     method = "nearest", data = imp)
dta_m <- match.data(mod_match)

dta_m$reten_Y1 <- "promoted"
dta_m[dta_m$treat_Y1 == 1,]$reten_Y1 <- "retained"


fn_bal <- function(dta, variable, newvar) {
  dta$variable <- dta[, variable]
  dta$variable <- as.numeric(dta$variable)
  
  dta$treat_Y1 <- as.factor(dta$treat_Y1)
  support <- c(min(dta$variable), max(dta$variable))
  ggplot(dta, aes(x = distance, y = variable, color = reten_Y1)) +
    geom_point(alpha = 0.2, size = 1.3) +
    geom_smooth(method = "loess", se = F) +
    xlab("Propensity score") +
    ylab(newvar) +
    theme_bw() +
    ylim(support)
}

library(gridExtra)
library(ggplot2)
grid.arrange(
  fn_bal(dta_m, "Math_Y0", "Math Y0") + theme(legend.position = "none"),
  fn_bal(dta_m, "Koep_Y0", "Umbrella Y0"),
  fn_bal(dta_m, "Ind_Y0", "Independent Y0") + theme(legend.position = "none"),
  fn_bal(dta_m, "Etn", "Etn"),
  nrow = 2, widths = c(1, 0.8)
)

## density of matched propensity score
# 682 rows
dens_retained <- dta_m[dta_m$treat_Y1 == 0,]$distance
dens_prom <- dta_m[dta_m$treat_Y1 == 1,]$distance

plot(density(dens_retained), xlab = "Propensity score",
     main= "Propensity score of the ovelapped \n between retained and promoted in year 1", col = "blue")
lines(density(dens_prom) ,col = "red")

legend("topright",legend = c("retained", "promoted"), fill = c("blue", "red"))

sum(dta_m$treat_Y1 == 1)

######### check balance for other covariates

time_varying_list = c("PercGOK_Y0","Wel_Y0","Att_Y0","Aso_Y0","Coop_Y0",
                      "Ind_Y0","Peer_Y0","Pro_Y0","Hyp_Y0","Agr_Y0","Self_Y0","SuppHE_Y0",
                      "Rel_Y0","PSE_Y0","TC_Lang_Y0","PercGOK_Y0","Math_Y0","squared_Math_Y0")

pvalue_list <- c()
for (v in time_varying_list){
  t <- t.test(dta_m[, v] ~ dta_m$treat_Y1)
  p = t$p.value
  pvalue_list = append(pvalue_list, p)
}


names(pvalue_list) <- time_varying_list
df = data.frame(pvalue_list)
df$name <- time_varying_list

########### propensity score of retention year 2
null.mod <- glm(treat_Y2 ~ Math_Y1 + Math_Y0 + treat_Y1, data = imp, family = binomial)

full.mod <- glm(treat_Y2 ~ Sex+Monthb+Hlang+SES+GOK+Etn+
                  PM_Y1+Koep_Y1+Net_Y1+Schooltype_Y1+PercGOK_Y1+Wel_Y1+Att_Y1+Aso_Y1+Coop_Y1+
                  Ind_Y1+Peer_Y1+Pro_Y1+Hyp_Y1+Agr_Y1+Self_Y1+SuppHE_Y1+
                  Rel_Y1+PSE_Y1+TC_Lang_Y1+Koep_Y1+PM_Y1+Schooltype_Y1+Net_Y1+
                  PercGOK_Y1 + Math_Y1 + squared_Math_Y1 + Math_Y0 + treat_Y1, data = imp, family = binomial)

mod <- step(null.mod, scope = list(lower = null.mod, upper = full.mod), direction = "both", trace = 1)
mod$call

cont_list <- c("Math_Y1", "Math_Y0", "Rel_Y1", "TC_Lang_Y1","Att_Y1",
               "squared_Math_Y1","PSE_Y1","Hyp_Y1","Wel_Y1",
               "Ind_Y1", "Peer_Y1", "Coop_Y1")

cat_list <- c("treat_Y1", "Schooltype_Y1", "Net_Y1", "Etn", "Sex")
inter_list = create_inter(cont_list, cat_list)

inter_formu <- paste0(inter_list, collapse = "+")

full.formu <- glm(treat_Y2 ~ Math_Y1 + Math_Y0 + treat_Y1 + Rel_Y1 + 
                    TC_Lang_Y1 + Att_Y1 + squared_Math_Y1 + PSE_Y1 + Hyp_Y1 + 
                    Wel_Y1 + Schooltype_Y1 + Ind_Y1 + Net_Y1 + Peer_Y1 + Etn + 
                    Coop_Y1 + Sex + Math_Y1:treat_Y1+Math_Y1:Schooltype_Y1+Math_Y1:Net_Y1+
                    Math_Y1:Etn+Math_Y1:Sex+Math_Y0:treat_Y1+Math_Y0:Schooltype_Y1+Math_Y0:Net_Y1+
                    Math_Y0:Etn+Math_Y0:Sex+Rel_Y1:treat_Y1+Rel_Y1:Schooltype_Y1+Rel_Y1:Net_Y1+
                    Rel_Y1:Etn+Rel_Y1:Sex+TC_Lang_Y1:treat_Y1+TC_Lang_Y1:Schooltype_Y1+TC_Lang_Y1:Net_Y1+
                    TC_Lang_Y1:Etn+TC_Lang_Y1:Sex+Att_Y1:treat_Y1+Att_Y1:Schooltype_Y1+Att_Y1:Net_Y1+Att_Y1:Etn+
                    Att_Y1:Sex+squared_Math_Y1:treat_Y1+squared_Math_Y1:Schooltype_Y1+squared_Math_Y1:Net_Y1+
                    squared_Math_Y1:Etn+squared_Math_Y1:Sex+PSE_Y1:treat_Y1+PSE_Y1:Schooltype_Y1+PSE_Y1:Net_Y1+
                    PSE_Y1:Etn+PSE_Y1:Sex+Hyp_Y1:treat_Y1+Hyp_Y1:Schooltype_Y1+Hyp_Y1:Net_Y1+Hyp_Y1:Etn+Hyp_Y1:Sex+
                    Wel_Y1:treat_Y1+Wel_Y1:Schooltype_Y1+Wel_Y1:Net_Y1+Wel_Y1:Etn+Wel_Y1:Sex+Ind_Y1:treat_Y1+Ind_Y1:Schooltype_Y1+
                    Ind_Y1:Net_Y1+Ind_Y1:Etn+Ind_Y1:Sex+Peer_Y1:treat_Y1+Peer_Y1:Schooltype_Y1+Peer_Y1:Net_Y1+Peer_Y1:Etn+Peer_Y1:Sex+
                    Coop_Y1:treat_Y1+Coop_Y1:Schooltype_Y1+Coop_Y1:Net_Y1+Coop_Y1:Etn+Coop_Y1:Sex, family = binomial, data = imp)
Anova(full.formu, type = "II")

full.formu <- glm(treat_Y2 ~ Math_Y1 + Math_Y0 + treat_Y1 + Rel_Y1 + 
                    TC_Lang_Y1 + Att_Y1 + squared_Math_Y1 + PSE_Y1 + Hyp_Y1 + 
                    Wel_Y1 + Schooltype_Y1 + Ind_Y1 + Net_Y1 + Peer_Y1 + Etn + 
                    Coop_Y1 + Sex + Math_Y1:Net_Y1 + Math_Y0:Net_Y1 + Rel_Y1:Schooltype_Y1 + 
                    Att_Y1:Net_Y1+squared_Math_Y1:Net_Y1 +
                    Schooltype_Y1:Ind_Y1 + Net_Y1:Peer_Y1, family = binomial, data = imp)
Anova(full.formu, type = "II")
summary(full.formu)
mod_match <- matchit(treat_Y2 ~ Math_Y1 + Math_Y0 + treat_Y1 + Rel_Y1 + 
                       TC_Lang_Y1 + Att_Y1 + squared_Math_Y1 + PSE_Y1 + Hyp_Y1 + 
                       Wel_Y1 + Schooltype_Y1 + Ind_Y1 + Net_Y1 + Peer_Y1 + Etn + 
                       Coop_Y1 + Sex + Math_Y1:Net_Y1 + Math_Y0:Net_Y1 + Rel_Y1:Schooltype_Y1 + 
                       Att_Y1:Net_Y1+squared_Math_Y1:Net_Y1 +
                       Schooltype_Y1:Ind_Y1 + Net_Y1:Peer_Y1 + squared_Math_Y0,
                     method = "nearest", data = imp)
dta_m <- match.data(mod_match)

dta_m$retention_Y2 <- "promoted"
dta_m[dta_m$treat_Y2 == 1,]$retention_Y2 <- "retained"

dens_retained <- dta_m[dta_m$treat_Y2 == 0,]$distance
dens_prom <- dta_m[dta_m$treat_Y2 == 1,]$distance

plot(density(dens_retained), xlab = "Propensity score",
     main= "Propensity score of the ovelapped \n between retained and promoted in year 2", col = "blue")
lines(density(dens_prom) ,col = "red")

legend("topleft",legend = c("retained", "promoted"), fill = c("blue", "red"))

sum(dta_m$treat_Y2 == 1)

dta_m <- match.data(mod_match)

dta_m$retention_Y2 <- "promoted"
dta_m[dta_m$treat_Y2 == 1,]$retention_Y2 <- "retained"


fn_bal <- function(dta, variable, newvar) {
  dta$variable <- dta[, variable]
  dta$variable <- as.numeric(dta$variable)
  
  dta$treat_Y2 <- as.factor(dta$treat_Y2)
  support <- c(min(dta$variable), max(dta$variable))
  ggplot(dta, aes(x = distance, y = variable, color = retention_Y2)) +
    #geom_point(alpha = 0.2, size = 1.3) +
    geom_smooth(method = "loess", se = F) +
    xlab("Propensity score") +
    ylab(newvar) +
    theme_bw() +
    ylim(support)
}

library(gridExtra)
library(ggplot2)
grid.arrange(
  fn_bal(dta_m, "Math_Y0", "Math Y0") + theme(legend.position = "none"),
  fn_bal(dta_m, "Koep_Y1", "Umbrella Y1"),
  fn_bal(dta_m, "Math_Y1", "Math Y1") + theme(legend.position = "none"),
  fn_bal(dta_m, "Etn", "Etn"),
  nrow = 2, widths = c(1, 0.8)
)

### test balance
cont_list <- c("Math_Y1", "Math_Y0", "Rel_Y1", "TC_Lang_Y1","Att_Y1",
               "squared_Math_Y1","PSE_Y1","Hyp_Y1","Wel_Y1",
               "Ind_Y1", "Peer_Y1", "Coop_Y1")

cat_list <- c("treat_Y1", "Schooltype_Y1", "Net_Y1", "Etn", "Sex")

pvalue_list <- c()
for (v in cont_list){
  t <- t.test(dta_m[, v] ~ dta_m$treat_Y2)
  p = t$p.value
  pvalue_list = append(pvalue_list, p)
}

df <- data.frame(pvalue_list)
df$name <- cont_list

ps <- dta_m$distance
t <- lm(Math_Y1 ~ distance*treat_Y2, data = dta_m)

summary(t)
Anova(t, type = "II")
