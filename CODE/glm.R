#D4, D5 worst feeling for patient
#D4, D5, D6, D7 can vary
#D2, D8 could look more similar 
#MSUrine Metabolite variation (biggest change in metabolite number) Emphasis on D5, D6
##################################################################################

library(tidyverse)

raw<-read_csv("MS-Urine.csv")
raw$ID<-str_split_i(raw$ID_and_DAY, "-D",1)
raw<-raw%>% relocate(ID, .before=ID_and_DAY)

# lets compare all five days compared to d2 
metabolite_df<-raw%>% janitor::clean_names()
metabolite_df<-metabolite_df%>% select(id:age, s_3_hydroxyisobutyric_acid_hmdb0000023:xanthosine_hmdb0000299)

#process metabolite_df 
metabolite_df<-metabolite_df%>% mutate(sex=case_when(sex=="F"~1, sex=="M"~0))
covariates<-c("sex", "age")
metabolites<-names(metabolite_df)[6:ncol(metabolite_df)]

#remove outliers and scale
# Function to remove outliers from a single column
remove_outliers_column<- function(column) {
  median <- median(column,na.rm = TRUE)
  IQR <- IQR(column, na.rm=TRUE)
  lower_bound <- median - 5 * IQR
  upper_bound <- median + 5 * IQR
  column[column < lower_bound | column > upper_bound] <- NA
  return(column)
}

met_outliers<-metabolite_df
met_outliers[metabolites] <- lapply(met_outliers[metabolites], remove_outliers_column)
met_out_scale<-met_outliers
met_out_scale[,metabolites]<-scale(met_out_scale[,metabolites])


# Association test ------------------------------------------------------

perform_mwas<- function(input, covariates,sensitivity=FALSE){
  
  res<-tibble()
  sense<-tibble()
  
  for (i in 1:length(metabolites)){
    
    biomarker_name<- metabolites[i]
    f<-formula(paste("label ~ ", paste(c(biomarker_name, covariates), collapse=" + ")))

    fit<-glm(formula = f, family=binomial(link="logit"), data=input)
    
    coef<-summary(fit)$coef
    temp<-tibble(metabolite=biomarker_name, predictor=c("intercept", "metabolite", covariates), as_tibble(coef))
    names(temp)<-c("metabolite", "predictor","estimate", "se", "z", "p")
    temp$auc<-NA
    res<-bind_rows(res, temp)
    
    test_prob = predict(fit, newdata = input, type = "response")
    pROC_obj <- pROC::roc(input$label,test_prob,
                    plot=FALSE)
    auc<-tibble(metabolite=biomarker_name,
           predictor="auc",
           estimate=NA, 
           se=NA, 
           z=NA, 
           p=NA, 
           auc=as.numeric(pROC_obj$auc))
    res<-bind_rows(res, auc)
    
    if (sensitivity==TRUE){
      ci<-confint(fit)
      ev<-evalues.OR(est = exp(as.numeric(coef(fit)[2])),
                     lo=exp(ci)[2,1], hi=exp(ci)[2,2], rare=TRUE)
      ev<-apply(ev, 2, as.numeric) %>% as_tibble()
      
      ev<-ev%>% mutate(metabolite=biomarker_name)%>% relocate(metabolite, .before="point")
      ev<-ev%>% mutate(type=c("OR", "Eval")) %>% relocate(type, .after="metabolite")
      sense<-bind_rows(sense, ev)
    }
  }
  res<-res %>%
    pivot_wider(names_from = predictor,
                values_from = c(estimate, se, z, p, auc),
                names_glue = "{predictor}_{.value}")
  auc<-res%>%select(auc_auc)%>% rename(auc="auc_auc")
  res<-res%>% select(-contains("_auc"), -contains("auc_"))%>% bind_cols(auc)
  
  res<-res[ , order(names(res))]
  res<-res%>% select(metabolite,starts_with("intercept"), starts_with("metabolite"), auc, everything())

  if(sensitivity==TRUE){
    return(list(res=res, sense=sense))
  }else{
    return(res)
  }
}

df<-metabolite_df
df<-met_out_scale

days<-unique(raw$Day)[-1]
list_res<-list() 

for (i in 1:length(days)) {
  
  input<-df%>% filter(day %in% c("D2", days[i]))
  input<-input%>% mutate(label=case_when(day=="D2"~0, 
                                       day==days[i]~1))
  #input[,metabolites]<-scale(input[,metabolites])
  
  
  res<- perform_mwas(input, covariates)
  
  list_res<-append(list_res, list(res))
  names(list_res)[i]<-days[i]
  
}

# Extract betas-----
extract_betas<- function(list_res, col, odds.ratio=FALSE){
  
  output<-list_res[[1]]%>% select(metabolite)
  
  for (i in 1:length(list_res)){
    df<-as.data.frame(list_res[i])
    suffix<-names(list_res)[i]
    
    df<-df%>% select(contains(paste0(col, "_estimate")), contains(paste0(col, "_p")))
    if (odds.ratio==TRUE){
      or<-df%>%select(contains("_estimate"))%>% exp()
      names(or)<-paste0(suffix,"_", col, "_or")
      df<-bind_cols(or, df)%>% select(-contains("_estimate"))
    }
    output<-bind_cols(output, df)
  }
  return (output)
}

biomarker_betas<-extract_betas(list_res, "metabolite")
write_csv(biomarker_betas, "res_urine/metabolite_betas_outlier_removed.csv") 
biomarker_or<-extract_betas(list_res, "metabolite", odds.ratio = TRUE)
write_csv(biomarker_or, "res_urine/metabolite_or_outlier_removed.csv") #use for finding sig biomarkerz

#Prepare df for forest plots
extract_summary<- function(list_res, col){
  
  for (i in 1:length(list_res)){
    df<-as.data.frame(list_res[i])
    suffix<-names(list_res)[i]
    
    trait_df<-df %>% rename_with(~str_remove(., paste0(suffix, "."))) %>% mutate(label=suffix)
    if (col=="metabolite"){
      trait_df<-trait_df %>% select(metabolite, label, starts_with("metabolite"))
    } else{
      trait_df<- trait_df %>% select(metabolite, label, starts_with(col))
    }
    
    if (i==1){
      output<-trait_df
    } else {
      output<-bind_rows(output, trait_df)
    }
  }
  return (output)
}

summary_betas<-extract_summary(list_res, "metabolite")
write_csv(summary_betas, "res_urine/summary_betas_outlier_removed.csv") #use this for ggforest
