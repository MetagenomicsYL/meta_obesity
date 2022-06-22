###########################function library#####################################
############# some of these functions come from MicrobeR ########
################################################################

reqpkg <- c("tidyverse", "broom", "ggthemes", "UpSetR", "RColorBrewer")
for(i in reqpkg)
{
  library(i, quietly = TRUE, verbose = FALSE, warn.conflicts = FALSE, character.only = TRUE)
}

TidyConvert.WhatAmI<-function(TABLE){
  class(TABLE)[1]
}

TidyConvert.ToMatrix<-function(TABLE, KEY){
  if(TidyConvert.WhatAmI(TABLE)=="matrix"){return(TABLE)}
  else if(TidyConvert.WhatAmI(TABLE)=="data.frame"){return(TABLE %>% as.matrix())}
  else if(TidyConvert.WhatAmI(TABLE)=="data.table"){return(TABLE %>% as.data.frame() %>% column_to_rownames(KEY) %>% as.matrix())}
  else if(TidyConvert.WhatAmI(TABLE)=="tibble"){return(TABLE %>% as.data.frame() %>% column_to_rownames(KEY) %>% as.matrix())}
}


Confidence.Filter<-function(FEATURES,MINSAMPS,MINREADS,VERBOSE){
  FEATURES<-TidyConvert.ToMatrix(FEATURES, colnames(FEATURES)[1])
  total_observations<-sum(FEATURES)  #get total number of reads
  
  if(missing(VERBOSE)){VERBOSE=T}
  if(VERBOSE==T){
    message(paste("Filtering features such that they are present in at least", MINSAMPS, "samples with a total of at least", MINREADS, "reads."))
    message(paste("...There are", total_observations, "reads and", nrow(FEATURES), " features"))
  }
  
  filtered<-FEATURES[apply(FEATURES,1,function(x){length(grep("TRUE",x!=0))}>=MINSAMPS),]
  filtered<-filtered[((rowSums(filtered))>=MINREADS),]
  
  if(VERBOSE==T){
    message(paste("...After filtering there are", sum(filtered), "reads and", nrow(filtered), "features"))
  }
  return(filtered)
}

REMfun <- function(tableupmeta, raceinfo){
  tableplot <- gather(tableupmeta,
                      key = "genus", value = "abundance", -Run, -StudyID, -BMI, -group, -race)
  BMIresult <- tibble()
  for (i in unique(tableplot$genus)){
    datalmm <- subset(tableplot, genus==i)
    fit<-lmerTest:::lmer(abundance~BMI+(1|StudyID), data=datalmm)
    beta <- summary(fit)$coefficients["BMI", "Estimate"]
    Pval <- summary(fit)$coefficients["BMI", "Pr(>|t|)"]
    OR <- exp(summary(fit)$coefficients["BMI", "Estimate"])
    cf<-lme4::confint.merMod(object = fit,method = "Wald")
    ORCI <- exp(cf)
    lowOR <- ORCI["BMI",1]
    highOR <- ORCI["BMI",2]
    low <- cf["BMI",1]
    high <- cf["BMI",2]
    so <- summary(fit)
    corr <- so$vcov@factors$correlation["BMI",1]
    BMIresult <- bind_rows(BMIresult, tibble(
      Study=raceinfo,Genus=i,OR=OR,low_OR=lowOR,high_OR=highOR,beta=beta,low_CI=low, high_CI=high, Pval= Pval,correlation=corr
    ))
  }
  return(BMIresult)
}

cal_logfun <- function(metagroup, fungather, acolname){
  funlog <- list()
  for (i in unique(metagroup$StudyID)){
    funsty <- fungather%>%
      filter(StudyID==i)
    fun_smr <- funsty %>%
      group_by(eval(parse(text=acolname)), group) %>%
      summarize(mean=mean(log2(reads))) %>%
      spread(key=group, value=mean) %>%
      rename(mean_log2_Healthy=Healthy, mean_log2_Obese=Obese)
    names(fun_smr)[1] <- acolname
    fun_joinlog <- left_join(funsty, fun_smr)
    fun_joinlog$log2FC <- log2(fun_joinlog$reads)-(fun_joinlog$mean_log2_Healthy)
    funlog[[i]] <- select(fun_joinlog, StudyID, everything())
  }
  return(funlog)
}

cal_log2fc <- function(metagroup, funlog, acolname){
  eachst <- list()
  st <- list()
  for (i in unique(metagroup$StudyID)){
    st$rawlog <- funlog[[i]]
    st$logstats <-  st$rawlog %>%
      group_by(eval(parse(text=acolname))) %>%  
      do(
        broom::tidy(t.test(log2FC~group, data=., conf.int=TRUE, conf.level=0.95)) 
      ) %>%
      mutate(StudyID=i)
    names(st$logstats)[1] <- acolname
    st$logstats <- st$logstats%>%
      select(StudyID, all_of(acolname), log2FC=estimate, Pvalue=p.value, mean_Healthy=estimate1, mean_Obese=estimate2, CI_low=conf.high, CI_high=conf.low)
    st$logstats$log2FC <- -st$logstats$log2FC
    st$logstats$CI_low <- -st$logstats$CI_low
    st$logstats$CI_high <- -st$logstats$CI_high
    eachst[[i]] <- st
  }
  return(eachst)
}


combine_log <- function(eachst){
  totaldf<-lapply(eachst, function(x) x$rawlog)%>%
    do.call(bind_rows, .) %>%
    mutate(group=factor(group, levels=c("Healthy","Obese")))
  return(totaldf)
}

combine_logfc <- function(totaldf, acolname){
  funCombined<-tibble()
  for(i in unique(totaldf[, acolname])){
    fit<-lmerTest:::lmer(log2FC~group+(1|StudyID), data=subset(totaldf, eval(parse(text=acolname))==i))
    cf<-confint.merMod(object = fit, method = "Wald")
    funCombined<-bind_rows(funCombined, tibble(
      StudyID="Combined", 
      acolname=i, 
      log2FC=summary(fit)$coefficients["groupObese", "Estimate"], 
      Pvalue=anova(fit)$`Pr(>F)`, 
      mean_Healthy=NA, 
      mean_Obese=NA, 
      CI_low=cf["groupObese",1], 
      CI_high=cf["groupObese",2]
    ))
  }
  names(funCombined)[2] <-acolname
  return(funCombined)
}


