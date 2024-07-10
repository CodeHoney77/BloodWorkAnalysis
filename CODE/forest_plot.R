library(tidyverse)
library(ggforestplot)
library(patchwork)

beta_summary <- read_csv("res_urine/summary_betas_outlier_removed.csv")
df<-beta_summary
df$metabolite<-str_split_i(df$metabolite, "_hmdb", 1)


#----select subgroups for each col----

metabolites<-unique(df$metabolite)
sg1<-metabolites[1:36]
sg2<- metabolites[37:73]
sg3<- metabolites[74:110]
sg4<- metabolites[111:147]
sg5<- metabolites[148:184]
sg6<- metabolites[185:length(metabolites)]


df2<-df %>% mutate(facet=case_when(metabolite %in% sg1~"col_1",
                                   metabolite %in% sg2~"col_2",
                                   metabolite %in% sg3~"col_3",
                                   metabolite %in% sg4~"col_4",
                                   metabolite %in% sg5~"col_5",
                                   metabolite %in% sg6~"col_6"))
 
sig<-0.05

#df3<-df2%>% filter(label%in% c("ME/CFS", "ME/CFS only"))
#----PLot----

plot_col<-function(df, col_as_string, legend=FALSE){
  g<-ggforestplot::forestplot(
    df = df %>% filter(facet==col_as_string),
    name = metabolite,
    estimate = metabolite_estimate,
    se = metabolite_se,
    colour = label,
    pvalue=metabolite_p,
    psignif=sig, 
    xlab = "Odds Ratio (95% CI), per 1-SD",
    logodds = TRUE) + 
    labs(colour="")
    #scale_color_manual(values = rev(colour$col[1:3]))
  if (legend==TRUE){
    g<-g+
      theme(legend.title = element_text(size=8),
                   legend.position = "top",
                   legend.justification = "top",
                   legend.text = element_text(size=8),
            axis.title=element_text(size=8),
            axis.text = element_text(size=8),
            strip.text = element_text(size=8)) +
      guides(colour=guide_legend(ncol=3,
                                 reverse = TRUE))
  }else{
    g<-g+
      theme(legend.position="none",
            axis.title=element_text(size=8),
            axis.text = element_text(size=8),
            strip.text = element_text(size=8))
  }
return(g)
}
plot_page<-function(df, col1, col2, col3){
  g1<-plot_col(df, col1, legend = TRUE)
  g2<-plot_col(df, col2, legend=TRUE)
  #g3<-plot_col(df, col3, legend=FALSE)
  p<-wrap_plots(g1,g2, ncol=2)
  return(p)
}

p1<-plot_page(df2, "col_1","col_2")
p2<-plot_page(df2, "col_3", "col_4")
p3<-plot_page(df2, "col_5", "col_6")


svg(file="res_urine/forest_page1_outlier_removed.svg", 
    width = 8.3, height=11.7)
p1
dev.off()

svg(file="res_urine/forest_page2_outlier_removed.svg", 
    width = 8.3, height=11.7)
p2
dev.off()

svg(file="res_urine/forest_page3_outlier_removed.svg", 
    width = 8.3, height=11.7)
p3
dev.off()
