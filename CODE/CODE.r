setwd("C:/Users/ah437/OneDrive/Documents/Work")

library(tidyverse)

raw_data<-read_csv("MS Blood.csv")
raw_data<-janitor:::clean_names(raw_data)

#change col names
names(raw_data)[1:4]<-c("id", "class", "day", "participant_id")


#order rows
raw_data<-raw_data%>% arrange(participant_id)

##### prepare datasets #####

#extract label
label<-raw_data%>% select(id:participant_id, severity)

#get sample size
sample<-label%>% select(-id, -day)%>% distinct()
sample%>% count(class) # ME n=61, HC n=19

#extract demographics (non-metaboltie data)
demo<-raw_data%>% select(id:order)

#extract metabolite data
demo_names<-names(demo)[!names(demo)%in%names(label)]
met_df<-raw_data%>%select(-any_of(demo_names))

# get d4-d1 (deltas) for case-control analysis 
d1<-met_df%>% filter(day==1)
d4<-met_df%>% filter(day==4)

#check if the participants in both d1 and d4 are irmn the same order
identical(d1$participant_id, d4$participant_id)
identical(d1$participant_id, sample$participant_id)

#calculate the metabolite change from d4-d1 
delta4_1<-d4[, 6:ncol(d4)]-d1[, 6:ncol(d1)]
delta4_1<-bind_cols(sample, delta4_1)

#cfs with all the days 
cfs_met<-met_df%>% filter(class=="ME")


#### data analysis ##### 
#### pca (unsupervised clustering) #####
pca_score<-function(data, label, colour, shape=NULL, pc1=1, pc2=2, title, loadings=FALSE){
  
  temp<-data%>% select(-any_of(c("id", "day", "class", "participant_id", "severity")))%>%mutate(across(where(is.numeric), ~replace_na(., median(., na.rm=TRUE))))
  res_pca<-prcomp(temp, scale=TRUE)
  
  eigs<-res_pca$sdev^2
  var_df<-tibble(
    PC=seq(1:length(eigs)),
    SD = sqrt(eigs),
    Proportion = eigs/sum(eigs),
    Cumulative = cumsum(eigs)/sum(eigs))                          
  
  score_df<-bind_cols(label, as_tibble(res_pca$x))
  
  if(is.null(shape)==TRUE){
    
    score_df2<-score_df%>% select(any_of(c(colour, paste0("PC",1), paste0("PC", 2))))
    names(score_df2)<-c("colour", "pc1", "pc2")
    
    s<-ggplot(score_df2, aes(x=pc1, y=pc2, color=colour))+
      geom_point()+theme_minimal()+labs(color=colour)+ggtitle(title)+
      theme(panel.border = element_rect(color="black", fill=NA),
            panel.grid = element_blank())+ylab(paste0("PC2 \u2013", round(var_df$Proportion[pc2],3)*100, "%"))+
      xlab(paste0("PC1 \u2013", round(var_df$Proportion[pc1],3)*100, "%"))
    
  }else{
    score_df2<-score_df%>% select(any_of(c(colour, shape, paste0("PC",1), paste0("PC", 2))))
    names(score_df2)<-c("colour", "shape", "pc1", "pc2")
    score_df2$colour<-factor(score_df2$colour)
    score_df2$shape<-factor(score_df2$shape)
    
    s<-ggplot(score_df2, aes(x=pc1, y=pc2, color=colour))+
      geom_point(aes(shape=shape))+theme_minimal()+labs(color=colour, shape=shape)+ggtitle(title)+
      theme(panel.border = element_rect(color="black", fill=NA),
            panel.grid = element_blank())+
      ylab(paste0("PC2 \u2013", round(var_df$Proportion[pc2],3)*100, "%"))+
      xlab(paste0("PC1 \u2013", round(var_df$Proportion[pc1],3)*100, "%"))
  }
  
  if(loadings==TRUE){
    loadings<-tibble(metabolite=names(temp), loadings=abs(res_pca$rotation[,pc1]))
    loadings<-loadings%>% arrange(desc(loadings))
    
    l<-loadings%>% slice(1:10)%>%
      ggplot(aes(x=metabolite, y=loadings))+
      geom_bar(stat = "identity")+
      theme(panel.border = element_rect(color="black", fill=NA),
            panel.grid = element_blank())+
      xlab("")+ylab(" Loadings")+coord_flip()
    return(list(score=s, loadings=l))
  }else{
    list(score=s)
  }
}

patchwork:::wrap_plots(c(pca_score(met_df, met_df%>% select(any_of(c(names(label)))), "class", "day", title="All datapoints", loadings = TRUE),
pca_score(delta4_1, delta4_1%>% select(any_of(c(names(label)))), "class", shape=NULL, title="Delta D4-D1", loadings = TRUE),
pca_score(cfs_met, cfs_met%>% select(any_of(c(names(label)))), "severity", "day", title="ME/CFS only", loadings = TRUE)), ncol = 2)

ggsave("MS_blood_pca.svg", width=8.7, height=7, units = "in")

#### plsda supervised clustering) ##### 

plsda_score<-function(data, colour,title){
  
  temp<-data%>% select(-any_of(c("id", "day", "class", "participant_id", "severity")))%>%mutate(across(where(is.numeric), ~replace_na(., median(., na.rm=TRUE))))
  
  y<-data%>% select(all_of(colour))
  names(y)<-"temp"
  groups<-y%>% count(temp)%>% mutate(group=c("class1","class2"))
  y<-y$temp
  
  x<-temp
  res_plsda<-mixOmics:::plsda(x,y, ncomp=2, scale=TRUE)
  var_explained<-mixOmics:::explained_variance(res_plsda$X, res_plsda$variates$X, ncomp=2)
  
  score_df<-bind_cols(data%>% select(all_of(colour)), as_tibble(res_plsda$variates$X))
  names(score_df)[1]<-"colour"
  
  s<-ggplot(score_df, aes(x=comp1, y=comp2, color=colour))+
      geom_point()+theme_minimal()+ggtitle(title)+labs(colour=colour)+
      theme(panel.border = element_rect(color="black", fill=NA),
            panel.grid = element_blank())+
    ylab(paste0("PC2 \u2013", round(as.numeric(var_explained[2]),3)*100, "%"))+
      xlab(paste0("PC1 \u2013", round(as.numeric(var_explained[1]),3)*100, "%"))

  vip_scores<-mixOmics:::vip(res_plsda)
  vip_df<-tibble(metabolites=names(temp), comp1=vip_scores[,1])
  vip_df<-vip_df%>% arrange(desc(comp1))%>% slice(1:16)
  
  met_mean<-data%>% select(all_of(c(colour, vip_df$metabolites)))
  names(met_mean)[1]<-"class"
  met_mean<-met_mean%>% group_by(class)%>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE), .names = "{col}"))
  test<-t(met_mean)
  met_mean<-tibble(metabolite=rownames(test),class1=test[,1], class2=test[,2])
  names(met_mean)<-met_mean[1,]
  met_mean<-met_mean%>% slice(-1)
  names(met_mean)<-c("metabolites", groups$group)
  met_mean$class1<-as.numeric(met_mean$class1)
  met_mean$class2<-as.numeric(met_mean$class2)
  met_mean<-met_mean%>% mutate(in_class2=case_when(class2>class1~"Increase",
                                            TRUE~"Decrease"))
  vip_df<-left_join(vip_df, met_mean)
  
  vip_names<-str_split_i(vip_df$metabolites, "_hmdb",1)
  vip_df$metabolites<-vip_names
  vip_df$metabolites<-factor(vip_df$metabolites, levels = rev(vip_df$metabolites))
  vip_df$in_class2<-factor(vip_df$in_class2, levels = c("Increase", "Decrease"))
  
  lab_name<-groups%>% filter(group=="class2")
  v<-ggplot(vip_df, aes(x=metabolites, y=comp1, fill=in_class2))+
    geom_bar(stat = "identity")+theme_minimal()+
    theme(panel.border = element_rect(color="black", fill=NA),
          panel.grid = element_blank())+labs(fill=lab_name$temp)+
    xlab("")+ylab(" VIP score")+coord_flip()+scale_fill_manual(values=c("#9E0142","#053061"))
  
  patchwork:::wrap_plots(s,v, ncol=2)
}

patchwork::wrap_plots(plsda_score(delta4_1, "class", "Delta D4-D1"), 
                      plsda_score(cfs_met, "day", "ME/CFS only"), nrow=2)
ggsave("MS_blood_plsda.svg", width=9, height=6.5, units = "in")

#### fold change #### 

#calculate fold change
fc_df<-cfs_met%>% select(-any_of(c("id", "class", "participant_id", "severity")))
fc_df<-met_df%>%filter(day==1) %>% select(-any_of(c("id", "day", "participant_id", "severity")))

fc_df$day<-factor(fc_df$day)
median_df<-fc_df%>% group_by(day)%>% 
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE), .names = "{col}"))
median_df<-t(median_df)
median_df<-tibble(metabolites=rownames(median_df), d1=median_df[,1], d4=median_df[,2])%>% slice(-1)
median_df$HC<-as.numeric(median_df$HC)
median_df$ME<-as.numeric(median_df$ME)
median_df$d1<-as.numeric(median_df$d1)
median_df$d4<-as.numeric(median_df$d4)
median_df<-median_df%>% mutate(fc=d4/d1)%>%mutate(log2fc=log2(fc))

#calculate p-value 
wilcox_pvalue <- function(column, condition) {
  wilcox.test(column ~ condition)$p.value
}

# Apply the function to each metabolite column
p_values <- apply(fc_df[, -1], 2, wilcox_pvalue, condition = fc_df$day)

# Convert p-values to a data frame
p_values_df <- tibble(
  metabolites = names(p_values),
  pval = p_values
)

median_df<-left_join(median_df, p_values_df)
median_df<-median_df%>% mutate(neg_log10p=-log10(pval))


# P value adjustment here look at the p values up or down with out fold change 
sig_fc<-0.6
sig_pval<-0.05

median_df$diffexpressed<-"Not sig"
median_df$diffexpressed[median_df$log2fc > sig_fc & median_df$pval < sig_pval] <- "Up"
median_df$diffexpressed[median_df$log2fc < -sig_fc & median_df$pval < sig_pval] <- "Down"
median_df$diffexpressed<-factor(median_df$diffexpressed, levels=c("Up", "Not sig", "Down"))
median_df$metabolites<-str_split_i(test$metabolites, "_hmdb",1)
median_df<-median_df%>% mutate(label=case_when(diffexpressed=="Not sig"~"",
                                               TRUE~as.character(metabolites)))

#plot

ggplot(median_df, aes(x=log2fc, y=-log10(pval), col=diffexpressed, label=label))+
  geom_point() + 
  theme_minimal() + theme(panel.border = element_rect(color="black", fill=NA),
                          panel.grid = element_blank())+
  ggrepel:::geom_text_repel() + ggtitle("ME/CFS only; D4:D1")+
  #scale_color_manual(values=c("black", "blue")) +
  scale_color_manual(values=c("red", "black", "blue")) +
  geom_vline(xintercept=c(-sig_fc, sig_fc), col="grey", linetype="dashed") +
  geom_hline(yintercept=-log10(sig_pval), col="grey", linetype="dashed")+xlim(c(-1.5,1.5))+labs(col="")

ggsave("MS_blood_cfs_only_volcano_cfs.svg", width=8.7, height=5, units = "in")

