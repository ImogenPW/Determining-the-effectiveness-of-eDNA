
#title: "Chapter 4 -eDNA vs Tradtional morthological sampling"
#author: "Imogen Poyntz-Wright"
#date: "2024-02-29"


#Libaries
library(ggplot2)
library(gridExtra)
library(grid)
library(readr)
library(dplyr)
library(ggplot2)
library(maps)
library(ggplot2)
library(sf)
library(cowplot)
library(stringr)
library(MASS)  
library(lmtest)
library(DHARMa)
library(MuMIn)
library(olsrr)
library(car)
library(lme4)
library(brglm)
library(ggeffects)
library(foreign)
library(logistf)  
library(rms)
library(sandwich)
library(lmtest)
library(writexl)
library(openxlsx)


#datasets

###datasets for counting/percentage occurence of taxa across all samples/studies
df_ocurr_co1<- read.csv('C:/Users/imoge/OneDrive/phd/Research/Objective 7/Analysis/across all sites dataset/Corrected taxa missed by eDNA but detected by trad co1.csv')
df_ocurr_co1_terr<- read.csv('C:/Users/imoge/OneDrive/phd/Research/Objective 7/Analysis/across all sites dataset/Corrected taxa missed by eDNA but detected by trad co1 terrestrial.csv')
df_ocurr_16s<- read.csv('C:/Users/imoge/OneDrive/phd/Research/Objective 7/Analysis/across all sites dataset/Corrected taxa missed by eDNA but detected by trad 16s.csv')
df_ocurr_16s_terr<- read.csv('C:/Users/imoge/OneDrive/phd/Research/Objective 7/Analysis/across all sites dataset/Corrected taxa missed by eDNA but detected by trad 16s terrestrial.csv')
df_ocurr_18s<- read.csv('C:/Users/imoge/OneDrive/phd/Research/Objective 7/Analysis/across all sites dataset/Corrected taxa missed by eDNA but detected by trad 18s.csv')
percentage_not_recorded <- read.csv('C:/Users/imoge/OneDrive/phd/Research/Objective 7/Analysis/across all sites dataset/percentage and count all markers.csv')
percentage_not_recorded_terrestrial <- read.csv('C:/Users/imoge/OneDrive/phd/Research/Objective 7/Analysis/across all sites dataset/percentage and count all markers terrestrial.csv')

###dataset for number of taxa only detected by eDNA as a percentage of all taxa recorded present
percentage_not_recorded_trad_absent_present <- read.csv('C:/Users/imoge/OneDrive/phd/Research/Objective 7/Analysis/across all sites dataset/percentage and count all markers_ TAXA missed or detected by trad technique.csv')

###datasets for the number of unique taxa detected
df_all_counts_pre_abs <- read.csv('C:/Users/imoge/OneDrive/phd/Research/Objective 7/Analysis/unique taxa counts rivers.csv')
df_all_counts_pre_abs_TERR <- read.csv('C:/Users/imoge/OneDrive/phd/Research/Objective 7/Analysis/unique taxa counts terrestrial.csv')

###dataset for the numebr of studies which use which primer and marker 
summary_of_sites <-read.csv("C:/Users/imoge/OneDrive/phd/Research/Objective 7/Analysis/Summary of sites.csv")

###datasets for modelling (also the same datasets used to to produce occurence dataframes)
data_16_model <- read.csv('C:/Users/imoge/OneDrive/phd/Research/Objective 7/Analysis/final datasets/all_16s_taxa_for_modelling.csv')
data_18_model <- read.csv('C:/Users/imoge/OneDrive/phd/Research/Objective 7/Analysis/final datasets/all_18s_taxa_for_modelling.csv')
data_co1_model <- read.csv('C:/Users/imoge/OneDrive/phd/Research/Objective 7/Analysis/final datasets/all_co1_taxa_for_modelling.csv')





################################################################################

#SECTION 1; OCCURENCE OF MISSED TAXA OVER ALL STUDIES PER MARKER (supplementary figures)

#CO1
##family
df_ocurr_co1 <- df_ocurr_co1[df_ocurr_co1$Family != "" & !is.na(df_ocurr_co1$Family), ]
value_counts_Family <- table(df_ocurr_co1$Family)
value_counts_Family <- table(df_ocurr_co1$Family)
value_counts_Family_df <- as.data.frame(value_counts_Family)
value_counts_Family_df$Var1 <- factor(value_counts_Family_df$Var1, levels = value_counts_Family_df$Var1[order(-value_counts_Family_df$Freq)])
value_counts_Family_df$FacetGroup <- cut(as.numeric(value_counts_Family_df$Var1), breaks = 4, labels = c("Group1", "Group2", "Group3", "Group4"))
ggplot(value_counts_Family_df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Family", y = "Number of taxa") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 20),  # Increase x-axis label font size
        axis.text.y = element_text(size = 20),  # Increase y-axis label font size
        axis.title.x = element_text(size = 20),  # Increase x-axis title font size
        axis.title.y = element_text(size = 20),  # Increase y-axis title font size
        strip.text = element_text(size = 20),  # Increase facet label font size
        plot.title = element_text(size = 20)) +
  facet_wrap(~ FacetGroup, scales = "free_x", ncol = 1, nrow = 4)
#ggsave("CO1_family_taxa_occurrence.tiff", width = 30, height = 20, dpi = 330)

##phylum
selected_columnsco1 <- percentage_not_recorded[, c("count_not_detected_co1", "Phylum")]
selected_columnsco1 <- na.omit(selected_columnsco1)
plot_1<- ggplot(selected_columnsco1, aes(x = reorder(Phylum, -count_not_detected_co1), y = count_not_detected_co1)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Phyla", y = "Number of taxa") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle =45, hjust = 1,size = 20),  # Increase x-axis label font size
        axis.text.y = element_text(size = 20),  # Increase y-axis label font size
        axis.title.x = element_text(size = 20),  # Increase x-axis title font size
        axis.title.y = element_text(size = 20),  # Increase y-axis title font size
        strip.text = element_text(size = 20),  # Increase facet label font size
        plot.title = element_text(size = 20)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
#ggsave("CO1 phylum taxa occurence.tiff", dpi = 330)


##16s
#family
df_ocurr_16s <- df_ocurr_16s[df_ocurr_16s$Family != "" & !is.na(df_ocurr_16s$Family), ]
value_counts_Family <- table(df_ocurr_16s$Family)
value_counts_Family_df <- as.data.frame(value_counts_Family)
value_counts_Family_df$Var1 <- factor(value_counts_Family_df$Var1, levels = value_counts_Family_df$Var1[order(-value_counts_Family_df$Freq)])
ggplot(value_counts_Family_df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Family", y = "Number of taxa") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 20),  # Increase x-axis label font size
        axis.text.y = element_text(size = 20),  # Increase y-axis label font size
        axis.title.x = element_text(size = 20),  # Increase x-axis title font size
        axis.title.y = element_text(size = 20),  # Increase y-axis title font size
        strip.text = element_text(size = 20),  # Increase facet label font size
        plot.title = element_text(size = 20)) 
#ggsave("16S family taxa occurence.tiff", width = 15, height = 6, dpi = 330)

#phylum
selected_columnsco1 <- percentage_not_recorded[, c("count_not_detected_16", "Phylum")]
selected_columnsco1 <- na.omit(selected_columnsco1)# Create a bar plot with faceting into four graphs in one column
plot_2 <-ggplot(selected_columnsco1, aes(x = reorder(Phylum, -count_not_detected_16), y = count_not_detected_16)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Phyla", y = "Number of taxa") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 20),  # Increase x-axis label font size
        axis.text.y = element_text(size = 20),  # Increase y-axis label font size
        axis.title.x = element_text(size = 20),  # Increase x-axis title font size
        axis.title.y = element_text(size = 20),  # Increase y-axis title font size
        strip.text = element_text(size = 20),  # Increase facet label font size
        plot.title = element_text(size = 20)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
#ggsave("16s phylum taxa occurence.tiff",dpi = 330)


##18s
#family
value_counts_Family <- table(df_ocurr_18s$Family)
value_counts_Family_df <- as.data.frame(value_counts_Family)
value_counts_Family_df$Var1 <- factor(value_counts_Family_df$Var1, levels = value_counts_Family_df$Var1[order(-value_counts_Family_df$Freq)])
ggplot(value_counts_Family_df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Family", y = "Number of taxa") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 20),  # Increase x-axis label font size
        axis.text.y = element_text(size = 20),  # Increase y-axis label font size
        axis.title.x = element_text(size = 20),  # Increase x-axis title font size
        axis.title.y = element_text(size = 20),  # Increase y-axis title font size
        strip.text = element_text(size = 20),  # Increase facet label font size
        plot.title = element_text(size = 20)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
#ggsave("18s family taxa occurence.tiff",  dpi = 330)

#phylum
selected_columnsco1 <- percentage_not_recorded[, c("count_not_detected_18", "Phylum")]
selected_columnsco1 <- na.omit(selected_columnsco1)# 
plot_3<- ggplot(selected_columnsco1, aes(x = reorder(Phylum, -count_not_detected_18), y = count_not_detected_18)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Phyla", y = "Number of taxa") +
  theme_bw() + 
  theme(axis.text.x = element_text( angle = 45,hjust = 1,size = 20),  # Increase x-axis label font size
        axis.text.y = element_text(size = 20),  # Increase y-axis label font size
        axis.title.x = element_text(size = 20),  # Increase x-axis title font size
        axis.title.y = element_text(size = 20),  # Increase y-axis title font size
        strip.text = element_text(size = 20),  # Increase facet label font size
        plot.title = element_text(size = 20)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
#ggsave("18s phylum taxa occurence.tiff",dpi = 330)




#TERRESTRIAL

##co1
#family
df_ocurr_co1_terr <- df_ocurr_co1_terr[df_ocurr_co1_terr$family != "" & !is.na(df_ocurr_co1_terr$family), ]
value_counts_Family <- table(df_ocurr_co1_terr$family)
value_counts_Family <- table(df_ocurr_co1_terr$family)
value_counts_Family_df <- as.data.frame(value_counts_Family)
value_counts_Family_df$Var1 <- factor(value_counts_Family_df$Var1, levels = value_counts_Family_df$Var1[order(-value_counts_Family_df$Freq)])
value_counts_Family_df$FacetGroup <- cut(as.numeric(value_counts_Family_df$Var1), breaks = 4, labels = c("Group1", "Group2", "Group3", "Group4"))
ggplot(value_counts_Family_df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Family", y = "Number of taxa") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 20),  # Increase x-axis label font size
        axis.text.y = element_text(size = 20),  # Increase y-axis label font size
        axis.title.x = element_text(size = 20),  # Increase x-axis title font size
        axis.title.y = element_text(size = 20),  # Increase y-axis title font size
        strip.text = element_text(size = 20),  # Increase facet label font size
        plot.title = element_text(size = 20)) +
  facet_wrap(~ FacetGroup, scales = "free_x", ncol = 1, nrow = 4)
#ggsave("CO1_family_taxa_occurrence_terrestrial.tiff", width = 30, height = 20, dpi = 330)

#phylum

selected_columnsco1 <- percentage_not_recorded_terrestrial[, c("count_not_detected_co1", "Phylum")]
selected_columnsco1 <- na.omit(selected_columnsco1)# 
plot_1<- ggplot(selected_columnsco1, aes(x = reorder(Phylum, -count_not_detected_co1), y = count_not_detected_co1)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.5) +
  labs(x = "Phyla", y = "Number of taxa") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle =45, hjust = 1,size = 20),  # Increase x-axis label font size
        axis.text.y = element_text(size = 20),  # Increase y-axis label font size
        axis.title.x = element_text(size = 20),  # Increase x-axis title font size
        axis.title.y = element_text(size = 20),  # Increase y-axis title font size
        strip.text = element_text(size = 20),  # Increase facet label font size
        plot.title = element_text(size = 20)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
#ggsave("CO1 phylum taxa occurence_terrestrial.tiff", dpi = 330)


##16s
#family
value_counts_Family <- table(df_ocurr_16s_terr$family)
value_counts_Family_df <- as.data.frame(value_counts_Family)
value_counts_Family_df$Var1 <- factor(value_counts_Family_df$Var1, levels = value_counts_Family_df$Var1[order(-value_counts_Family_df$Freq)])
ggplot(value_counts_Family_df, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "skyblue", width =0.5) +
  labs(x = "Family", y = "Number of taxa") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 20),  # Increase x-axis label font size
        axis.text.y = element_text(size = 20),  # Increase y-axis label font size
        axis.title.x = element_text(size = 20),  # Increase x-axis title font size
        axis.title.y = element_text(size = 20),  # Increase y-axis title font size
        strip.text = element_text(size = 20),  # Increase facet label font size
        plot.title = element_text(size = 20)) 
#ggsave("16S family taxa occurence_terrestrial.tiff", dpi = 330)

#phylum
selected_columnsco1 <- percentage_not_recorded_terrestrial[, c("count_not_detected_16", "Phylum")]
selected_columnsco1 <- na.omit(selected_columnsco1)# 
plot_2 <-ggplot(selected_columnsco1, aes(x =  reorder(Phylum, -count_not_detected_16), y = count_not_detected_16)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.5) +
  labs(x = "Phyla", y = "Number of taxa") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,size = 20),  # Increase x-axis label font size
        axis.text.y = element_text(size = 20),  # Increase y-axis label font size
        axis.title.x = element_text(size = 20),  # Increase x-axis title font size
        axis.title.y = element_text(size = 20),  # Increase y-axis title font size
        strip.text = element_text(size = 20),  # Increase facet label font size
        plot.title = element_text(size = 20)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
#ggsave("16s phylum taxa occurence_terrestrial.tiff",dpi = 330)


################################################################################


#SECTION 2: PRESENCE AND ABSENCE OF UNIQUE TAXA (figure 4a)


#aquatic - all
selected_columnsco1 <- df_all_counts_pre_abs[, c("Phylum", "Count", "Fill")]
selected_columnsco1 <- na.omit(selected_columnsco1)
ggplot(selected_columnsco1, aes(x = reorder(Phylum, -Count), y = Count, fill = Fill)) +
  geom_bar(stat = "identity") +
 # scale_y_log10() +
  labs(x = "Phylum", y = "Number of unique taxa", fill = "Dataframe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = c("lightblue", "grey"), name = " ")  # Set bar colors 
#ggsave("Figure 4A.tiff",dpi = 330)


#terrestrial - all
selected_columnsco1 <- df_all_counts_pre_abs_TERR[, c("Phylum", "Count", "Fill")]
selected_columnsco1 <- na.omit(selected_columnsco1)
ggplot(selected_columnsco1, aes(x = reorder(Phylum, -Count), y = Count, fill = Fill)) +
  geom_bar(stat = "identity") +
  labs(x = "Phylum", y = "Number of unique taxa", fill = "Dataframe") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = c("lightblue", "grey"), name = " ")   # Set bar colors
  #scale_y_log10()  # Set y-axis to log scale
ggsave("Figure UNIQUE terrestrial TAXA.tiff",dpi = 330)


################################################################################

#SECTION 3 - PERCENTAGE OF TAXA WHICH IS NOT RECORDED BY eDNA but by trad (figure 4b - d)
selected_columnsco1 <- percentage_not_recorded[, c("percentage_not_detected_co1", "Phylum")]
selected_columnsco1 <- na.omit(selected_columnsco1)

P1<-ggplot(selected_columnsco1, aes(x = reorder(Phylum, -percentage_not_detected_co1 ), y = percentage_not_detected_co1)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Phyla", y = "Percentage of taxa") +
  theme_bw() + 
  theme(axis.text.x = element_text( angle = 45,hjust = 1,size = 20),  # Increase x-axis label font size
        axis.text.y = element_text(size = 20),  # Increase y-axis label font size
        axis.title.x = element_text(size = 20),  # Increase x-axis title font size
        axis.title.y = element_text(size = 20),  # Increase y-axis title font size
        strip.text = element_text(size = 20),  # Increase facet label font size
        plot.title = element_text(size = 20)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
#ggsave("CO1 percentage taxa not recorded.tiff", dpi = 330)


selected_columns16s <- percentage_not_recorded[, c("percentage_not_detected_16", "Phylum")]
selected_columns16s <- na.omit(selected_columns16s)

P2<-ggplot(selected_columns16s, aes(x = reorder(Phylum, -percentage_not_detected_16 ), y = percentage_not_detected_16)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Phyla", y = "Percentage of taxa") +
  theme_bw() + 
  theme(axis.text.x = element_text( angle = 45,hjust = 1,size = 20),  # Increase x-axis label font size
        axis.text.y = element_text(size = 20),  # Increase y-axis label font size
        axis.title.x = element_text(size = 20),  # Increase x-axis title font size
        axis.title.y = element_text(size = 20),  # Increase y-axis title font size
        strip.text = element_text(size = 20),  # Increase facet label font size
        plot.title = element_text(size = 20)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
#ggsave("16s percentage taxa not recorded.tiff", dpi = 330)


selected_columns18s <- percentage_not_recorded[, c("percentage_not_detected_18", "Phylum")]
selected_columns18s <- na.omit(selected_columns18s)

P3<-ggplot(selected_columns18s, aes(x = reorder(Phylum, -percentage_not_detected_18 ), y = percentage_not_detected_18)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Phyla", y = "Percentage of taxa") +
  theme_bw() + 
  theme(axis.text.x = element_text( angle = 45,hjust = 1,size = 20),  # Increase x-axis label font size
        axis.text.y = element_text(size = 20),  # Increase y-axis label font size
        axis.title.x = element_text(size = 20),  # Increase x-axis title font size
        axis.title.y = element_text(size = 20),  # Increase y-axis title font size
        strip.text = element_text(size = 20),  # Increase facet label font size
        plot.title = element_text(size = 20)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
#ggsave("18s percentage taxa not recorded.tiff", dpi = 330)


#TERRESTRIAL (supplementary figures)

selected_columnsco1 <- percentage_not_recorded_terrestrial[, c("percentage_not_detected_co1", "Phylum")]
selected_columnsco1 <- na.omit(selected_columnsco1)

P1<-ggplot(selected_columnsco1, aes(x = reorder(Phylum, -percentage_not_detected_co1 ), y = percentage_not_detected_co1)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Phyla", y = "Percentage of taxa") +
  theme_bw() + 
  theme(axis.text.x = element_text( angle = 45,hjust = 1,size = 20),  # Increase x-axis label font size
        axis.text.y = element_text(size = 20),  # Increase y-axis label font size
        axis.title.x = element_text(size = 20),  # Increase x-axis title font size
        axis.title.y = element_text(size = 20),  # Increase y-axis title font size
        strip.text = element_text(size = 20),  # Increase facet label font size
        plot.title = element_text(size = 20)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
#ggsave("CO1 percentage taxa not recorded terrestrial.tiff", dpi = 330)


selected_columns16s <- percentage_not_recorded_terrestrial[, c("percentage_not_detected_16", "Phylum")]
selected_columns16s <- na.omit(selected_columns16s)

P2<-ggplot(selected_columns16s, aes(x = reorder(Phylum, -percentage_not_detected_16 ), y = percentage_not_detected_16)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Phyla", y = "Percentage of taxa") +
  theme_bw() + 
  theme(axis.text.x = element_text( angle = 45,hjust = 1,size = 20),  # Increase x-axis label font size
        axis.text.y = element_text(size = 20),  # Increase y-axis label font size
        axis.title.x = element_text(size = 20),  # Increase x-axis title font size
        axis.title.y = element_text(size = 20),  # Increase y-axis title font size
        strip.text = element_text(size = 20),  # Increase facet label font size
        plot.title = element_text(size = 20)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
#ggsave("16s percentage taxa not recorded terrestrial.tiff", dpi = 330)



################################################################################
#SECTION 3B - Number of riverine taxa detected by eDNA and missed by traditional methods (as percentage of total taxa recorded by both trad and eDNA) (figure 5)

#Absent = the number of taxa detected by only eDNA
#Total = the total number of all taxa detected by both eDNA and traditional methods

trad_missed_16s <- percentage_not_recorded_trad_absent_present[percentage_not_recorded_trad_absent_present$Fill_trad_per_16  == "Absent", ]
selected_columns16s <- trad_missed_16s [, c("Percentage_not_detected_Trad_16", "Phylum" )]
selected_columns16s <- na.omit(selected_columns16s)

P1<-ggplot(selected_columns16s, aes(x = reorder(Phylum, -Percentage_not_detected_Trad_16   ), y = Percentage_not_detected_Trad_16  )) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Phyla", y = "Percentage of taxa") +
  theme_bw() + 
  theme(axis.text.x = element_text( angle = 45,hjust = 1,size = 20),  # Increase x-axis label font size
        axis.text.y = element_text(size = 20),  # Increase y-axis label font size
        axis.title.x = element_text(size = 20),  # Increase x-axis title font size
        axis.title.y = element_text(size = 20),  # Increase y-axis title font size
        strip.text = element_text(size = 20),  # Increase facet label font size
        plot.title = element_text(size = 20)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
#ggsave("Percentage of taxa missed only detected by eDNA 16s.tiff", dpi = 330)


trad_missed_18s <- percentage_not_recorded_trad_absent_present[percentage_not_recorded_trad_absent_present$Fill_trad_per_18  == "Absent", ]
selected_columns18s <- trad_missed_18s [, c("Percentage_not_detected_Trad_18", "Phylum" )]
selected_columns18s <- na.omit(selected_columns18s)

P2<-ggplot(selected_columns18s, aes(x = reorder(Phylum, -Percentage_not_detected_Trad_18   ), y = Percentage_not_detected_Trad_18  )) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Phyla", y = "Percentage of taxa") +
  theme_bw() + 
  theme(axis.text.x = element_text( angle = 45,hjust = 1,size = 20),  # Increase x-axis label font size
        axis.text.y = element_text(size = 20),  # Increase y-axis label font size
        axis.title.x = element_text(size = 20),  # Increase x-axis title font size
        axis.title.y = element_text(size = 20),  # Increase y-axis title font size
        strip.text = element_text(size = 20),  # Increase facet label font size
        plot.title = element_text(size = 20)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
#ggsave( "Percentage of taxa only detected by eDNA 18s.tiff", dpi = 330)


trad_missed_co1 <- percentage_not_recorded_trad_absent_present[percentage_not_recorded_trad_absent_present$Fill_per_trad_co1  == "Absent", ]
selected_columnsco1 <- trad_missed_co1 [, c("Percentage_not_detected_Trad_co1", "Phylum" )]
selected_columnsco1 <- na.omit(selected_columnsco1)

P3<-ggplot(selected_columnsco1, aes(x = reorder(Phylum, -Percentage_not_detected_Trad_co1 ), y = Percentage_not_detected_Trad_co1  )) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(x = "Phyla", y = "Percentage of taxa") +
  theme_bw() + 
  theme(axis.text.x = element_text( angle = 45,hjust = 1,size = 20),  # Increase x-axis label font size
        axis.text.y = element_text(size = 20),  # Increase y-axis label font size
        axis.title.x = element_text(size = 20),  # Increase x-axis title font size
        axis.title.y = element_text(size = 20),  # Increase y-axis title font size
        strip.text = element_text(size = 20),  # Increase facet label font size
        plot.title = element_text(size = 20)) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))
#ggsave( "Percentage of taxa only detected by eDNA co1.tiff", dpi = 330)




################################################################################
#SECTION 4 - MAP (figure 1)

#make dataset
data <- data.frame(
  city = c("Wales", "Switzerland", "Japan", "Spain", "Belarus", "Greece", "South Africa", "Canada", "Germany", "New Zealand", "China", "England"),
  lat = c(52.30, 46.8182, 36.2048, 39.3999, 53.7098, 38.9954, -29.0469, 56.1304, 51.0647, -36.848461, 39.916668, 51.509865),
  lon = c( -3.70, 8.2275, 138.2529, -3.2264, 27.9534, 23.8052, 25.0979, -106.3468, 12.0128, 174.8860, 116.383331, -0.118092),
  count = c(1, 3, 1,2,1,1,1,1, 1, 2, 1, 1), 
  sampling_type = c("Aquatic", "Aquatic", "Aquatic", "Aquatic", "Aquatic", "Aquatic", "Aquatic", "Aquatic", "Terrestrial", "Terrestrial", "Terrestrial", "Terrestrial"))
#no.compar= c(14, )

# Plot a world map
world_map <- ggplot() +
  geom_polygon(data = map_data("world"), aes(x = long, y = lat, group = group), 
               fill = "lightblue", color = "white") +
  geom_point(data = data, aes(x = lon, y = lat,  shape = sampling_type, size = count,), color = "red") +
  scale_size_continuous(range = c(1, 3), breaks = c(1, 2, 3)) +  # Adjust size range as needed
  scale_shape_manual(values = c("Terrestrial" = 17, "Aquatic" = 19)) +  # Set shapes for specific groups
  theme_void()
print(world_map)
ggsave("world_map.tiff", world_map,  dpi=330)


#PLOTS FOR MAP - NUMBER OF STUDIES PER MARKER/PRIMER

##MARKER 
selected_columns <- c('Marker', 'Comparisons_markersite', 'Sample_type1')

# Use subset function to select specific columns
selected_data <- subset(summary_of_sites, select = selected_columns)
selected_data <- na.omit(selected_data)

p2<- ggplot(selected_data, aes(x  = Comparisons_markersite, y = reorder(Marker, -Comparisons_markersite), fill =Sample_type1)) +
  geom_bar(stat = "identity",  alpha = 0.7, position ='dodge') +
  labs(x = "Number of studies", y = "eDNA Marker") + theme_classic() +
  scale_fill_manual(
    values = c("lightgreen", "grey"),
    name = "Renamed Groups",
    labels = c("Riverine", "Terrestrial")
  ) +  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),  # Increase x-axis text size
             axis.text.y = element_text(size = 20),  # Increase y-axis text size
             axis.title = element_text(size = 20),  # Increase axis title size
             plot.title = element_text(size = 20)) 
#ggsave("Marker_studies.tiff", p2, width = 10, height = 5, dpi=330)

##PRIMER

selected_columns <- c('Primer', 'Comparisons_primer', 'Sample_type2')

# Use subset function to select specific columns
selected_data2 <- subset(summary_of_sites, select = selected_columns)
selected_data2 <- na.omit(selected_data2)
p1 <- ggplot(selected_data2, aes(x  = Comparisons_primer, y = reorder(Primer, -Comparisons_primer), fill =Sample_type2)) +
  geom_bar(stat = "identity",  alpha = 0.7, position ='dodge') +
  labs(x = "Number of studies", y = "eDNA Primer") + theme_classic() +
  scale_fill_manual(
    values = c("lightgreen", "grey"),
    name = "Renamed Groups",
    labels = c("Riverine", "Terrestrial")
  ) +  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20),  # Increase x-axis text size
             axis.text.y = element_text(size = 20),  # Increase y-axis text size
             axis.title = element_text(size = 20),  # Increase axis title size
             plot.title = element_text(size = 20)) 
#scale_y_discrete(labels = function(labels) gsub("/", "\n", labels))
ggsave("primer_studies.tiff", p1, width = 10, height = 5, dpi=330)


############################################################################################
#Section 5 - binomial mixed model for testing taxon-specific bias

#datasets
head(data_18_model) 
head(data_co1_model)
head(data_16_model) 


#fit  binomial models

#model - CO1
model <- glm(eDNA ~  title , data = data_co1_model, family = binomial) #model for impact of phylum and study on presence/absence of taxa
model <- glm(eDNA ~ Phylum , data = data_co1_model, family = binomial) #model for impact of phylum and study on presence/absence of taxa

# Assuming mydf has a "Phylum" column
mydf <- ggpredict(model, terms = c("Phylum"))
# Create a data frame with all unique combinations of Phylum - Compute heteroscedasticity-robust standard errors
robust_se <- sqrt(diag(vcovHC(model, type = "HC1")))
# Repeat the robust standard errors for each row
mydf$robust_se <- rep(robust_se, length.out = nrow(mydf))
# Display the modified data frame
print(mydf)
#MODEL PHYLUM AND TITLE SEPERATELY

#write_xlsx(mydf, "CO1_model_phylum.xlsx")
p1 <- plot(mydf)+   theme(axis.title = element_text(size = 14)) +   labs(y = "Predicted likelihood of taxa detected")
#ggsave(p1, file = "co1-model output phylum.tiff", height = 6.89, width = 8, dpi = 330)


#title
mydf<- ggpredict(model, terms = "title")
# Create a data frame with all unique combinations of Phylum - Compute heteroscedasticity-robust standard errors
robust_se <- sqrt(diag(vcovHC(model, type = "HC1")))
# Repeat the robust standard errors for each row
mydf$robust_se <- rep(robust_se, length.out = nrow(mydf))
# Display the modified data frame
print(mydf)
#write_xlsx(mydf, "CO1_model_title.xlsx")

p2 <- plot(mydf) +   theme(axis.title = element_text(size = 14)) +   labs(y = "Predicted likelihood of taxa detected")
#ggsave(p2, file = "co1-model output title.tiff", height = 6.89, width = 8, dpi = 330)

# clustering
      # Checking the clustering of phylum in the plot
residuals_df <- data.frame(Residuals = residuals(model, type = "pearson"), Phylum = data_co1_model$Phylum)
boxplot(Residuals ~ Phylum, data = residuals_df, xlab = "Phylum", ylab = "Residuals", main = "Residuals by Phylum")

      # Checking the clustering of study in the plot
residuals_df <- data.frame(Residuals = residuals(model, type = "pearson"), title_df1 = data_co1_model$title_df1)
boxplot(Residuals ~ title_df1, data = residuals_df, xlab = "title_df1", ylab = "Residuals", main = "Residuals by title_df1")

# Check for heteroscedasticity
# Extract residuals using resid() function
fitted_values <- model$fitted.values
residuals <- residuals(model)
# Check the lengths
length(fitted_values)
length(residuals)
# Now, create the plot
plot(fitted_values, residuals, xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs. Fitted Values")
abline(h = 0, col = "red", lty = 2)
lines(lowess(fitted_values, residuals), col = "green", lty = 2)
      
      #Scale-Location Plot (Square Root of Standardized Residuals vs. Fitted Values):
sqrt_abs_std_resid <- sqrt(abs(rstandard(model)))
plot(model$fitted.values, sqrt_abs_std_resid, xlab = "Fitted Values", ylab = "Square Root of Standardized Residuals", main = "Scale-Location Plot")

      #leverage plot
plot(hatvalues(model), residuals(model)^2, xlab = "Leverage", ylab = "Residuals Squared", main = "Leverage vs. Residuals Squared")
bptest(model)

#test overdispersion and normality of data
testDispersion(model)
model_test<- simulateResiduals(model, refit=T)
testOverdispersion(model_test)
plot(model_test) #NOT OVERDISPERSED + normality of data tested

#multicolinaerity
vif(model) #not multicolinear

#or
chi_square_test <- chisq.test(table(combined_df$Phylum, combined_df$marker))
print(chi_square_test)

#outliers
plot(cooks.distance(model), type = "p", pch = 1, main = "Cook's distance")
abline(h = 4/(nobs(model)-ncol(model)), col = "red", lty = 2)
influence_obj <- influence.measures(model)
cook_distance <- influence_obj$cooks.distance
# Identify observations with high Cook's distance
outliers <- which(cook_distance > 4 / nrow(data_16_model3))  # Using a common threshold of 4/n
# Print or inspect the indices of the identified outliers
print(outliers) #NO outliers

#AICc selection
options(na.action = "na.fail")
Sys.setlocale("LC_ALL", "en_US.UTF-8")
dredge1 <- dredge(model, rank = "AICc") #AIC values all within value of 6
dredge1 #null model is much worse
write_xlsx(dredge1, "AICc_CO1_phyla.xlsx")



###model - 16s - insufficent data
subset_data <- data_16_model_corr[, c("title", "eDNA", "Phylum")]

#check for separation
# Check for Perfect Prediction
freq_table <- table(subset_data$eDNA)
print(freq_table) #there is seperation - low values

# Create Logistic Regression Model (using logistf for penalized likelihood)
model_firth <- logistf(eDNA ~ Phylum, data = subset_data, family = binomial)
summary(model_firth)
colnames(subset_data)

# Display the robust standard errors
model <-glm(eDNA ~ Phylum, data = subset_data, family = binomial)
robust_se <- sqrt(diag(vcovHC(model, type = "HC1")))
# Display the robust standard errors
cat("Robust Standard Errors:", robust_se, "\n")
# Conduct hypothesis tests using the robust standard errors
summary(model, robust = TRUE)

# Assuming "Phylum" is the only predictor variable in your subset_data
predictor_matrix <- as.matrix(subset_data[, "Phylum", drop = FALSE])
response_vector <- as.vector(subset_data$eDNA)

# Create a dummy matrix with a second column of zeros
dummy_matrix <- cbind(predictor_matrix, 0)

# Use cv.glmnet with predefined lambda values
lambda_values <- 10^seq(10, -2, length = 100)
model_lasso <- cv.glmnet(dummy_matrix, response_vector, family = "binomial", alpha = 1, lambda = lambda_values)

mydf<- ggpredict(model, terms = "Phylum")
# Create a data frame with all unique combinations of Phylum - Compute heteroscedasticity-robust standard errors
  robust_se <- sqrt(diag(vcovHC(model, type = "HC1")))
# Repeat the robust standard errors for each row
  mydf$robust_se <- rep(robust_se, length.out = nrow(mydf))
# Display the modified data frame
print(mydf)

p3<- plot(mydf)

#ggsave(p3, file = "16s-model output phylum.tiff", dpi = 330)

mydf<- ggpredict(model, terms = "title_df1")

p3<- plot(mydf)

# clustering
# Checking the clustering of phylum in the plot
residuals_df <- data.frame(Residuals = residuals(model, type = "pearson"), Phylum = data_16_model_corr$Phylum)
boxplot(Residuals ~ Phylum, data = residuals_df, xlab = "Phylum", ylab = "Residuals", main = "Residuals by Phylum")

# Checking the clustering of study in the plot
residuals_df <- data.frame(Residuals = residuals(model, type = "pearson"), title_df1 = data_16_model_corr$title_df1)
boxplot(Residuals ~ title_df1, data = residuals_df, xlab = "title_df1", ylab = "Residuals", main = "Residuals by title_df1")

# Check for heteroscedasticity
# Extract residuals using resid() function
fitted_values <- model$fitted.values
residuals <- residuals(model)
# Check the lengths
length(fitted_values)
length(residuals)
# Now, create the plot
plot(fitted_values, residuals, xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs. Fitted Values")
abline(h = 0, col = "red", lty = 2)
lines(lowess(fitted_values, residuals), col = "green", lty = 2)

#Scale-Location Plot (Square Root of Standardized Residuals vs. Fitted Values):
sqrt_abs_std_resid <- sqrt(abs(rstandard(model)))
plot(model$fitted.values, sqrt_abs_std_resid, xlab = "Fitted Values", ylab = "Square Root of Standardized Residuals", main = "Scale-Location Plot")

#leverage plot
plot(hatvalues(model), residuals(model)^2, xlab = "Leverage", ylab = "Residuals Squared", main = "Leverage vs. Residuals Squared")
bptest(model)

#test overdispersion and normality of data
testDispersion(model)
model_test<- simulateResiduals(model, refit=T)
testOverdispersion(model_test)
plot(model_test) #NOT OVERDISPERSED + normality of data tested

#multicolinaerity
vif(model) #not multicolinear

#outliers
plot(cooks.distance(model), type = "p", pch = 1, main = "Cook's distance")
abline(h = 4/(nobs(model)-ncol(model)), col = "red", lty = 2)
influence_obj <- influence.measures(model)
cook_distance <- influence_obj$cooks.distance
# Identify observations with high Cook's distance
outliers <- which(cook_distance > 4 / nrow(data_16_model3))  # Using a common threshold of 4/n
# Print or inspect the indices of the identified outliers
print(outliers) #NO outliers

#AICc selection
options(na.action = "na.fail")
dredge1 <- dredge(model_bayesian, rank = "AICc") #AIC values all within value of 6
dredge1 #null model is much worse




###model - 18s - insufficent data
model <- glm(eDNA ~ Phylum , data = data_18_model, family = binomial) #model for impact of phylum and study on presence/absence of taxa
# Compute heteroscedasticity-robust standard errors
robust_se <- sqrt(diag(vcovHC(model, type = "HC1")))
# Display the robust standard errors
cat("Robust Standard Errors:", robust_se, "\n")
# Conduct hypothesis tests using the robust standard errors
summary(model, robust = TRUE)

mydf<- ggpredict(model, terms = "Phylum")

plot(mydf)

# clustering
# Checking the clustering of phylum in the plot
residuals_df <- data.frame(Residuals = residuals(model, type = "pearson"), Phylum = data_18_model$Phylum)
boxplot(Residuals ~ Phylum, data = residuals_df, xlab = "Phylum", ylab = "Residuals", main = "Residuals by Phylum")

# Check for heteroscedasticity
# Extract residuals using resid() function
fitted_values <- model$fitted.values
residuals <- residuals(model)
# Check the lengths
length(fitted_values)
length(residuals)
# Now, create the plot
plot(fitted_values, residuals, xlab = "Fitted Values", ylab = "Residuals", main = "Residuals vs. Fitted Values")
abline(h = 0, col = "red", lty = 2)
lines(lowess(fitted_values, residuals), col = "green", lty = 2)

#Scale-Location Plot (Square Root of Standardized Residuals vs. Fitted Values):
sqrt_abs_std_resid <- sqrt(abs(rstandard(model)))
plot(model$fitted.values, sqrt_abs_std_resid, xlab = "Fitted Values", ylab = "Square Root of Standardized Residuals", main = "Scale-Location Plot")

#leverage plot
plot(hatvalues(model), residuals(model)^2, xlab = "Leverage", ylab = "Residuals Squared", main = "Leverage vs. Residuals Squared")
bptest(model)

#test overdispersion and normality of data
testDispersion(model)
model_test<- simulateResiduals(model, refit=T)
testOverdispersion(model_test)
plot(model_test) #NOT OVERDISPERSED + normality of data tested

#multicolinaerity
vif(model) #not multicolinear

#outliers
plot(cooks.distance(model), type = "p", pch = 1, main = "Cook's distance")
abline(h = 4/(nobs(model)-ncol(model)), col = "red", lty = 2)
influence_obj <- influence.measures(model)
cook_distance <- influence_obj$cooks.distance
# Identify observations with high Cook's distance
outliers <- which(cook_distance > 4 / nrow(data_16_model3))  # Using a common threshold of 4/n
# Print or inspect the indices of the identified outliers
print(outliers) #NO outliers

#AICc selection
options(na.action = "na.fail")
dredge1 <- dredge(model, rank = "AICc") #AIC values all within value of 6
dredge1 #null model is much worse


