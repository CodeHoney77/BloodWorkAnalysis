library(tidyverse)

raw<-read_csv("MS-Urine.csv")
raw$ID<-str_split_i(raw$ID_and_DAY, "-D",1)
raw<-raw%>% relocate(ID, .before=ID_and_DAY)
raw<-raw%>% janitor::clean_names()

outlier_met<-read_csv("outlier_count_metabolite.csv")
# lets compare all five days compared to d2 

met_df<-raw%>% select(id:age, s_3_hydroxyisobutyric_acid_hmdb0000023:xanthosine_hmdb0000299)
met_select<-met_df%>% select(all_of(c("id_and_day", "id", "day", outlier_met$metabolite)))

met_scale<-met_select
met_scale[,outlier_met$metabolite]<-scale(met_scale[,outlier_met$metabolite])

#complex heatmap 
library(ComplexHeatmap)
library(circlize)


mat<-t(as.matrix(met_scale%>%select(-id, -id_and_day, -day)))
col_order<-met_scale%>% select(id, id_and_day, day)

col_fun<-colorRamp2(c(-2.5, 0, 8.5), c("blue", "white", "red"))

rownames(mat)<-str_split_i(rownames(mat), "_hmdb", 1)

column_ha = HeatmapAnnotation(ID = col_order$id, Day = col_order$day,
                              annotation_name_gp= gpar(fontsize = 8), height = unit(2, "mm"))

Heatmap(mat,name="Z-score", col=col_fun,
             show_row_dend = TRUE, cluster_columns=TRUE,
             row_names_gp = gpar(fontsize = 4),
        top_annotation = column_ha)
