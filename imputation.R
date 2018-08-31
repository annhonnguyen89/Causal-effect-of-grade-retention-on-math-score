library(mice)
library(lattice)
library(VIM)


load("widedf1")
load("Stud4b.Rdata")
widedf[,c("treat_Y0")] <- 0
widedf[,c("treat_Y1", "treat_Y2", "treat_Y3", "treat_Y4", "treat_Y5", "treat_Y6")] <- Stud4b[,c("treat_Y1", "treat_Y2", "treat_Y3", "treat_Y4", "treat_Y5", "treat_Y6")]
widedfnoid[,c("treat_Y1", "treat_Y2", "treat_Y3", "treat_Y4", "treat_Y5", "treat_Y6")] <- Stud4b[,c("treat_Y1", "treat_Y2", "treat_Y3", "treat_Y4", "treat_Y5", "treat_Y6")] 
ncol = ncol(widedf)
widedfnoid <- as.data.frame(widedf[,2:ncol])
source("allvar.R")

widedfnoid$treat_Y1 <- NA
widedfnoid$treat_Y2 <- NA
widedfnoid$treat_Y3 <- NA
widedfnoid$treat_Y4 <- NA
widedfnoid$treat_Y5 <- NA
widedfnoid$treat_Y6 <- NA

beforeimpute <- function(){
  t <- md.pattern(widedfnoid)
  tmp <- widedfnoid[imputeorder]
  pat <- md.pattern(tmp)
  
  sum(is.na(widedfnoid[,"Sex"]))
  sum(is.na(widedfnoid[,"Monthb"]))
  sum(is.na(widedfnoid[,"Hlang"]))
  sum(is.na(widedfnoid[,"SES"]))
  sum(is.na(widedfnoid[,"GOK"]))
  sum(is.na(widedfnoid[,"Etn"]))
  
}

modify_treat <- function(dat){
  n = nrow(dat)
  # treat_Y1
  for (i in c(1:n)){
    print(i)
    if(is.na(as.character(dat[i,"T_Y1"])) == FALSE && as.character(dat[i,"T_Y1"])!="BuLO" && is.na(as.character(dat[i,"T_Y0"])) == FALSE && as.character(dat[i,"T_Y0"])!="BuLO")
    {  
      if(as.character(dat[i,"T_Y1"]) == as.character(dat[i, "T_Y0"]))
        dat[i,"treat_Y1"] = 1
      else
        dat[i,"treat_Y1"] = 0
    } else dat[i,"treat_Y1"] = NA
    
    if(is.na(as.character(dat[i,"T_Y2"])) == FALSE && as.character(dat[i,"T_Y2"])!="BuLO" && is.na(as.character(dat[i,"T_Y1"])) == FALSE && as.character(dat[i,"T_Y1"])!="BuLO")
    {
      if(as.character(dat[i,"T_Y2"]) == as.character(dat[i, "T_Y1"]))
        dat[i,"treat_Y2"] = 1
      else
        dat[i,"treat_Y2"] = 0
    } else dat[i,"treat_Y2"] = NA
      
    if(is.na(as.character(dat[i,"T_Y3"])) == FALSE && as.character(dat[i,"T_Y3"])!="BuLO" && is.na(as.character(dat[i,"T_Y2"])) == FALSE && as.character(dat[i,"T_Y2"])!="BuLO")
    {
      if(as.character(dat[i,"T_Y3"]) == as.character(dat[i, "T_Y2"]))
        dat[i,"treat_Y3"] = 1
      else
        dat[i,"treat_Y3"] = 0
    } else dat[i,"treat_Y3"] = NA
    
    if(is.na(as.character(dat[i,"T_Y4"])) == FALSE && as.character(dat[i,"T_Y4"])!="BuLO" && is.na(as.character(dat[i,"T_Y3"])) == FALSE && as.character(dat[i,"T_Y3"])!="BuLO")
    {
      if(as.character(dat[i,"T_Y4"]) == as.character(dat[i, "T_Y3"]))
        dat[i,"treat_Y4"] = 1
      else
        dat[i,"treat_Y4"] = 0
    } else dat[i,"treat_Y4"] = NA
      
    if(is.na(as.character(dat[i,"T_Y5"])) == FALSE && as.character(dat[i,"T_Y5"])!="BuLO" && is.na(as.character(dat[i,"T_Y4"])) == FALSE && as.character(dat[i,"T_Y4"])!="BuLO")
    {
      if(as.character(dat[i,"T_Y5"]) == as.character(dat[i, "T_Y4"]))
        dat[i,"treat_Y5"] = 1
      else
        dat[i,"treat_Y5"] = 0
    } else dat[i,"treat_Y5"] = NA
      
    if(is.na(as.character(dat[i,"T_Y6"])) == FALSE && as.character(dat[i,"T_Y6"])!="BuLO" && is.na(as.character(dat[i,"T_Y5"])) == FALSE && as.character(dat[i,"T_Y5"])!="BuLO")
    {
      if(as.character(dat[i,"T_Y6"]) == as.character(dat[i, "T_Y5"]))
        dat[i,"treat_Y6"] = 1
      else
        dat[i,"treat_Y6"] = 0
    } else dat[i,"treat_Y6"] = NA
  }
  return(dat)
}


widedfnoid <- modify_treat(widedfnoid)
widedfnoid <- within(widedfnoid,{
  treat_Y1 <- as.factor(treat_Y1)
  treat_Y2 <- as.factor(treat_Y2)
  treat_Y3 <- as.factor(treat_Y3)
  treat_Y4 <- as.factor(treat_Y4)
  treat_Y5 <- as.factor(treat_Y5)
  treat_Y6 <- as.factor(treat_Y6)
})

source("D:\\Gent_Thesis\\code\\main\\allvar.R")

########################################################################################
transf <- function(df){
  df$Math_Y0 <- (df$Math_Y0)^2
  df$Math_Y1 <- (df$Math_Y1)^2
  df$Math_Y2 <- (df$Math_Y2)^2.5
  df$Math_Y3 <- (df$Math_Y3)^2.5
  df$Math_Y4 <- (df$Math_Y4)^2.5
  df$Math_Y5 <- (df$Math_Y5)^2.5
  df$Math_Y6 <- (df$Math_Y6)^2.5
  return(df)
}

widedf <- transf(widedf)
widedfnoid <- transf(widedfnoid)


source("allvar.R")  
source("imputation_func.R")
widedfnoid <- widedfnoid[imputeorder]


ini <- mice(widedfnoid, max = 0, print = TRUE, nnet.MaxNWts = 5000)
meth <- ini$method
pred <- ini$pred
misspat <- md.pattern(widedfnoid)
nrow <- nrow(misspat)
names <- names(widedfnoid)

# run again to get new form of prediction matrix which include the interaction
ini <- mice(widedfnoid, meth = meth, max = 0, print = TRUE, nnet.MaxNWts = 5000)
meth <- ini$method
pred <- ini$pred
varnames <- names(widedfnoid)
post <- squeeze(ini, varnames)

pred1 <- modifypredmat1(pred)

ini <- mice(widedfnoid, method = 'cart', predictorMatrix = pred1, post = post, print = TRUE, max = 20, m = 3, seed = 9212, ridge = 0.001, nnet.MaxNWts = 10000)
meth <- ini$method
pred <- ini$pred
pred0 <- pred

################### diagnotics
library(mice)
library(lattice)
library(VIM)

plot(ini, c("Math_Y0", "Math_Y1", "Math_Y2"))
plot(ini, c("Math_Y3", "Math_Y4", "Math_Y5", "Math_Y6"))
plot(ini, c("Aso_Y0", "Aso_Y1", "Aso_Y2", "Aso_Y3"))

library(ggplot2)
library(mice)
plot(ini, c("Math_Y0", "Math_Y1", "Math_Y2", "Math_Y3", "Math_Y4", "Math_Y5", "Math_Y6"))
plot(ini, "Math_Y6")

assess_imputed_treat <- function(trt_year, stream){
  trt_year = 1
  stream = 1
  
  tmp <- complete(ini,stream)
  tmp$trtna <- "observed"
  targetvar <- paste0("Math_Y", trt_year)
  tmp[is.na(widedf[,targetvar]) == TRUE, ]$trtna <- "missing"
  tmp$trtna <- as.factor(tmp$trtna)
  
  t = as.matrix(table(tmp[,targetvar], tmp$trtna))
  barplot(t, ylim = c(0, 5000))
}

assess_dens <- function(math_year, stream){
  library(ggplot2)
  tmp <- complete(ini,stream)
  tmp$Mathna <- "observed"
  targetvar <- paste0("Math_Y", math_year)
  tmp[is.na(widedf[,targetvar]) == TRUE, ]$Mathna <- "missing"
  tmp$Mathna <- as.factor(tmp$Mathna)
  
  
  p <- ggplot(tmp, aes(x=tmp[,targetvar])) + 
    geom_density(aes(group=tmp$Mathna, colour=tmp$Mathna)) +
    theme(legend.title=element_blank(), legend.position="none") +
    labs(x = targetvar)
  return(p)
  
}

p0_1 <- assess_dens(math_year = 0, stream = 1)
p0_2 <- assess_dens(math_year = 0, stream = 2)
p0_3 <- assess_dens(math_year = 0, stream = 3)

p1_1 <- assess_dens(math_year = 1, stream = 1)
p1_2 <- assess_dens(math_year = 1, stream = 2)
p1_3 <- assess_dens(math_year = 1, stream = 3)

p2_1 <- assess_dens(math_year = 2, stream = 1)
p2_2 <- assess_dens(math_year = 2, stream = 2)
p2_3 <- assess_dens(math_year = 2, stream = 3)

multiplot(p0_1, p1_1, p2_1, p0_2, p1_2, p2_2, p0_3, p1_3, p2_3, cols=3)

p3_1 <- assess_dens(math_year = 3, stream = 1)
p3_2 <- assess_dens(math_year = 3, stream = 2)
p3_3 <- assess_dens(math_year = 3, stream = 3)

p4_1 <- assess_dens(math_year = 4, stream = 1)
p4_2 <- assess_dens(math_year = 4, stream = 2)
p4_3 <- assess_dens(math_year = 4, stream = 3)

p5_1 <- assess_dens(math_year = 5, stream = 1)
p5_2 <- assess_dens(math_year = 5, stream = 2)
p5_3 <- assess_dens(math_year = 5, stream = 3)

p6_1 <- assess_dens(math_year = 6, stream = 1)
p6_2 <- assess_dens(math_year = 6, stream = 2)
p6_3 <- assess_dens(math_year = 6, stream = 3)

multiplot(p3_1, p4_1, p6_1, p3_2, p4_2, p6_2, p3_3, p4_3, p6_3, cols=3)

multiplot(p6_1, p6_2, p6_3, cols=3)

imp1 <- complete(ini,1)
md.pattern(imp1)

imp1 <- complete(ini,2)
md.pattern(imp1)

imp1 <- complete(ini,3)
md.pattern(imp1)


