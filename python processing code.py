#title: "Chapter 4 -eDNA vs Tradtional morthological sampling"
#author: "Imogen Poyntz-Wright"
#date: "2024-02-29"

#Libraries
 library(reticulate)
 py_config()
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
from scipy.stats import shapiro
from scipy.stats import anderson
from numpy import log as ln
import math
import matplotlib.pyplot as plt
import scipy.stats as stats



#Load dataset

##calculate log-ratios
df_1 = pd.read_excel('dataset - analysis 1.xlsx')
df_16s = pd.read_excel('dataset - analysis 1 16s.xlsx')
df_18s = pd.read_excel('dataset - analysis 1 18s.xlsx')
df_terrestrial= pd.read_excel('dataset - analysis 1 terrestrial.xlsx')

##occurence of riverine taxa
df_2 = pd.read_excel('dataset - analysis 2.xlsx')
df_2_16s = pd.read_excel('dataset - analysis 2_16s.xlsx')
df_2_18s = pd.read_excel('dataset - analysis 2_18ss.xlsx')
df_2_all = pd.read_excel('dataset - analysis 2_all.xlsx')

##counts of the number of taxa absent per study but recored by trad
df_spainCO1 = pd.read_excel('counts_absent_spainCO1.xlsx')
df_spain18S = pd.read_excel('counts_absent_spain18S.xlsx')
df_canada = pd.read_excel('counts_absent_canada.xlsx')
df_wales = pd.read_excel('counts_absent_wales.xlsx')
df_switz4 = pd.read_excel('counts_absent_switz4.xlsx')
df_us = pd.read_excel('counts_absent_us.xlsx')
df_switz3 = pd.read_excel('counts_absent_switz3.xlsx')
df_japan = pd.read_excel('counts_absent_japan.xlsx')
df_switz2 = pd.read_excel('counts_absent_switz2.xlsx')
df_belarus = pd.read_excel('counts_absent_bareus.xlsx')
df_greece = pd.read_excel('counts_absent_greece.xlsx')
df_southafrica = pd.read_excel('counts_absent_south_africa.xlsx')
df_spain2 = pd.read_excel('counts_absent_spain2.xlsx')
df_demark_missing = pd.read_excel('missing denmark 1 terr edna.xlsx')
df_china1_missing = pd.read_excel('missing china 1 terr edna.xlsx')
df_newzeal2_missing = pd.read_excel('missing newzealand 2 terr edna.xlsx')
df_engalnd1_missing = pd.read_excel('missing england 1 terr edna.xlsx')
df_newzeal1_missing = pd.read_excel('missing newzealand 1 terr edna.xlsx')

#for log-ratio calculation of CO1 markers
BF2_BR2_co1_for_log = pd.read_excel('BF2 BR2_co1_for_log.xlsx')
fwhF2_EPTDr2n_co1_for_log = pd.read_excel('fwhF2 EPTDr2n_co1_for_log.xlsx')
m1COIintF_jgHCO2198_co1_for_log = pd.read_excel('m1COIintF jgHCO2198_co1_for_log.xlsx')
LCO1490HCO2198_co1_for_log = pd.read_excel('LCO1490HCO2198_co1_for_log.xlsx')

#occurence of terrestrial taxa
CO1_terrestrial = pd.read_excel('CO1 terrestrial data.xlsx')
S16_terrestrial = pd.read_excel('16s terrestrial taxa.xlsx')

#outputs of log-ratio
merged_df_co1_corr = pd.read_excel('merged_df_co1.xlsx')
merged_18S_df_corr = pd.read_excel('merged_18S_df.xlsx')
merged_16S_df_corr = pd.read_excel('merged_16S_df.xlsx')
merged_16S_df_terr = pd.read_excel('merged_16S_terr_df.xlsx')
merged_co1_df_terr = pd.read_excel('merged_df_terr_CO1.xlsx')

#For all taxa with phyla linked
all_taxa_present_corr = pd.read_excel('all_taxa_present.xlsx')



################################################################################
#SECTION 1 - log-ratios 
##Calculate the log-ratios per marker (CO1/18S/16S). This will determine which eDNA marker detects more taxa relative to the number of taxa detected by traditional methods. 

#CALCULATE: Ln(A/B) 

# Calculate Ln(eDNA_no_taxa_marker_CO1 / trad_no_taxa) and create a new column
df_1['eDNA_no_taxa_marker_CO1'] = df_1['eDNA_no_taxa_marker_CO1'].replace(0, 0.001) #otherwise log will not work
CO1_df = df_1[df_1['eDNA_no_taxa_marker_CO1'].notna()].copy()
CO1_df['Ln_CO1'] = CO1_df.apply(lambda row: ln(row['eDNA_no_taxa_marker_CO1'] / row['trad_no_taxa']) if not math.isnan(row['eDNA_no_taxa_marker_CO1']) and not math.isnan(row['trad_no_taxa']) else math.nan, axis=1)
#Each row in the dataframe has both a value for eDNA and traditional method

# Calculate Ln(eDNA_no_taxa_marker_18S / trad_no_taxa) and create a new column
df_18S = df_1[df_1['eDNA_no_taxa_ maker_18S'].notna()].copy()
df_18S['Ln_18S'] = df_18S.apply(lambda row: ln(row['eDNA_no_taxa_ maker_18S'] / row['trad_no_taxa']) if not math.isnan(row['eDNA_no_taxa_ maker_18S']) and not math.isnan(row['trad_no_taxa']) else math.nan, axis=1)

# Calculate Ln(eDNA_no_taxa_marker_16S / trad_no_taxa) and create a new column
df_16S = df_16s[df_16s['eDNA no. taxa_marker 16S'].notna()].copy()
df_16S['Ln_16S'] = df_16S.apply(lambda row: ln(row['eDNA no. taxa_marker 16S'] / row['Trad no. taxa']) if not math.isnan(row['eDNA no. taxa_marker 16S']) and not math.isnan(row['Trad no. taxa']) else math.nan, axis=1)
average_16S = df_16S['Ln_16S'].mean(skipna=True)
print(f"The average is: {average_16S}")

    
#PLOT: Density graphs

##CO1
sns.displot(CO1_df, x="Ln_CO1", kind="kde", fill = 'black', alpha=.5)
plt.title("density plot co1", loc = 'left')
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.tight_layout()
#plt.savefig('Density plot CO1.tiff',dpi = 330) 
plt.show()

###AVERGAE value is:
CO1_average = CO1_df.loc[np.isfinite(CO1_df['Ln_CO1']), 'Ln_CO1'].mean()
print(f"The average is: {CO1_average}") 

###Make plot for CO1 - Per family
sns.displot(CO1_df, x="Ln_CO1", kind="kde", fill='black', alpha=.5, hue="taxa_level") #per taxa group
plt.title("Density Plot CO1", loc='left')
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.tight_layout()
#plt.savefig('Density plot CO1 and coloured by taxa.tiff',dpi = 330)
plt.show()


##18S 
sns.displot(df_18S, x="Ln_18S", kind="kde", fill = 'black', alpha=.5)
plt.title("density plot 18s", loc = 'left')
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.tight_layout()
#plt.savefig('Density plot 18s.tiff',dpi = 330)
plt.show()
average_18S = df_18S['Ln_18S'].mean(skipna=True)
print(f"The average is: {average_18S}")


##16s
sns.displot(df_16S, x="Ln_16S", kind="kde", fill = 'black', alpha=.5)
plt.title("density plot 16s", loc = 'left')
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.tight_layout()
#plt.savefig('Density plot 16s.tiff',dpi = 330)  
average_16S = df_16S['Ln_16S'].mean(skipna=True)
print(f"The average is: {average_16S}")    


################################################################################


#SECTION 2 -identifying missing taxa 
##Identifying the number of taxa missed per sample comparison (n = 170)
### comparison = the comparison number out of 170. site.no = the study site. Sample_no = the sample taken from the site. Really only need to filter by comparison.
### The number of unique taxa were counted directly from dataframe using EXCEL

#Extracted the needed columns for eDNA dataframes and traditional dataframes
CO1_taxa = df_2[['title','comparison', 'site_no', 'sample_no', 'CO1_marker_ eDNA_taxa']] #select columns
CO1_taxa = CO1_taxa.dropna(subset=['CO1_marker_ eDNA_taxa'])
CO1_taxa1 = CO1_taxa.copy()

S18_taxa = df_2_18s[['title','comparison', 'site_no', 'sample_no', '18S_marker_Edna_taxa']]
S18_taxa = S18_taxa.dropna(subset=['18S_marker_Edna_taxa'])
S18_taxa1 = S18_taxa.copy()
S18_taxa1.loc[:, 'method'] = '18S'

S16_taxa = df_2_16s[['title','comparison', 'site_no', 'sample_no', '16S_marker_eDNA_taxa']]
S16_taxa = S16_taxa.dropna(subset=['16S_marker_eDNA_taxa'])
S16_taxa1 = S16_taxa.copy()
S16_taxa1.loc[:, 'method'] = '16S'

trad_taxa_CO1 = df_2[['title','comparison', 'site_no', 'sample_no', 'tradno_taxa']]
trad_taxa_CO1 = trad_taxa_CO1.dropna(subset=['tradno_taxa'])
trad_taxa_CO2 = trad_taxa_CO1.copy()
#trad_taxa_CO2.loc[:, 'method'] = 'Traditional'

trad_taxa_18S1 = df_2_18s[['title','comparison', 'site_no', 'sample_no', 'trad_no_taxa']]
trad_taxa_18S1 = trad_taxa_18S1.dropna(subset=['trad_no_taxa'])
trad_taxa_18S2 = trad_taxa_18S1.copy()
trad_taxa_18S2.loc[:, 'method'] = 'Traditional'

trad_taxa_16S1 = df_2_16s[['title','comparison', 'site_no', 'sample_no', 'trad_no_taxa']]
trad_taxa_16S1 = trad_taxa_16S1.dropna(subset=['trad_no_taxa'])
trad_taxa_16S2 = trad_taxa_16S1.copy()
trad_taxa_16S2.loc[:, 'method'] = 'Traditional'

#Edit dataframes so they are ready to merge (add marker and rename columns)

###CO1 marker
CO1_taxa1_diff_cleaned = CO1_taxa1[['CO1_marker_ eDNA_taxa', 'title', 'comparison', 'site_no', 'sample_no']].copy()
CO1_taxa1_diff_cleaned['combin'] = 'co1' #add marker so know this is the eDNA data when merging
CO1_taxa1_diff_cleaned = CO1_taxa1_diff_cleaned.dropna().copy()
CO1_taxa1_diff_cleaned = CO1_taxa1_diff_cleaned.rename(columns={'CO1_marker_ eDNA_taxa': 'taxa'}) #rename column

trad_taxa_CO2_diff = trad_taxa_CO2[['title', 'comparison', 'site_no', 'sample_no', 'tradno_taxa']].copy()
trad_taxa_CO2_diff['combin'] = 'trad' #add marker so know this is the eDNA data when merging
trad_taxa_CO2_diff_cleaned = trad_taxa_CO2_diff.dropna().copy()
trad_taxa_CO2_diff_cleaned = trad_taxa_CO2_diff_cleaned.rename(columns={'tradno_taxa': 'taxa'})  #rename column

merged_df1 = pd.merge(
    trad_taxa_CO2_diff_cleaned,
    CO1_taxa1_diff_cleaned,
    on=['taxa', 'title', 'comparison', 'site_no', 'sample_no'],
    how='outer',
    suffixes=('_df1', '_df2')
) #merge rows based on 'taxa', 'title', 'comparison', 'site_no', 'sample_no'. This will provide a clear presence/absence of taxa per sample/study whether detected by both methods or one (eDNA/Trad)

#See excel file: all_co1_taxa_withstudyID (produced from the output above after merged with 'all_taxa_present.xlsx' too link phyla)


###16S marker
S16_taxa1_diff = S16_taxa1[['16S_marker_eDNA_taxa', 'title', 'comparison', 'site_no', 'sample_no']]
S16_taxa1_diff_cleaned = S16_taxa1_diff.dropna().copy()
S16_taxa1_diff_cleaned.loc[:, 'combin'] = '16S'
S16_taxa1_diff_cleaned = S16_taxa1_diff_cleaned.rename(columns={'16S_marker_eDNA_taxa': 'taxa'})

trad_taxa_16S1_diff = trad_taxa_16S1[['trad_no_taxa', 'title', 'comparison', 'site_no', 'sample_no']]
trad_taxa_16S1_diff_cleaned = trad_taxa_16S1_diff.dropna().copy()
trad_taxa_16S1_diff_cleaned.loc[:, 'combin'] = 'trad'
trad_taxa_16S1_diff_cleaned = trad_taxa_16S1_diff_cleaned.rename(columns={'trad_no_taxa': 'taxa'})

merged_df2 = pd.merge(
    trad_taxa_16S1_diff_cleaned, S16_taxa1_diff_cleaned,
    on=['taxa', 'title', 'comparison', 'site_no', 'sample_no'],
    how='outer',
    suffixes=('_df1', '_df2')
)
#See excel file: all_16_taxa_withstudyID (produced from the output above after merged with 'all_taxa_present.xlsx' too link phyla)


###18S marker
S18_taxa1_diff = S18_taxa1[['18S_marker_Edna_taxa', 'title', 'comparison', 'site_no', 'sample_no']]
S18_taxa1_diff_cleaned = S18_taxa1_diff.dropna().copy()
S18_taxa1_diff_cleaned.loc[:, 'combin'] = '18s'
S18_taxa1_diff_cleaned = S18_taxa1_diff_cleaned.rename(columns={'18S_marker_Edna_taxa': 'taxa'})

trad_taxa_18S1_diff = trad_taxa_18S1[['trad_no_taxa', 'title', 'comparison', 'site_no', 'sample_no']]
trad_taxa_18S1_diff_cleaned = trad_taxa_18S1_diff.dropna().copy()
trad_taxa_18S1_diff_cleaned.loc[:, 'combin'] = 'trad'
trad_taxa_18S1_diff_cleaned = trad_taxa_18S1_diff_cleaned.rename(columns={'trad_no_taxa': 'taxa'})

merged_df3 = pd.merge(
    trad_taxa_18S1_diff_cleaned,
    S18_taxa1_diff_cleaned,
    on=['taxa', 'title', 'comparison', 'site_no', 'sample_no'],
    how='outer',
    suffixes=('_df1', '_df2')
)

#See excel file: all_18_taxa_withstudyID (produced from the output above after merged with 'all_taxa_present.xlsx' too link phyla)

#UNIQUE TAXA
#In EXCEL use unique funtion to identify the number of unique taxa detected by trad sampling overall and the number of these detected by eDNA

################################################################################


#SECTION 3 - log-ratio per study
##Each study is addressed seperately in the following section.

#set plot parameters
extract dataframes by title of study
sns.set(rc='axes', labelsize=16, titlesize=18)
sns.set(rc='xtick', labelsize=14)
sns.set(rc='ytick', labelsize=14)
sns.set(rc='legend', fontsize=14)
sns.set(rc={'axes.labelsize': 16, 'axes.titlesize': 18, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 14})



##df1 - SPAIN: extract
column_name = 'title'
target_value = 'Evaluating freshwater macroinvertebrates from eDNA metabarcoding: A river Nalón case study' #select study

# Create a new DataFrame with rows that have the target value in the specified column: for CO1 marker
df1_spain = df[df[column_name] == target_value].copy() #Log-ratio #select data (rows)
# Reset the index of the new DataFrame
df1_spain.reset_index(drop=True, inplace=True)


# Create a new DataFrame with rows that have the target value in the specified column: for 16s marker
column_name = 'Title'
target_value = 'Evaluating freshwater macroinvertebrates from eDNA metabarcoding: A river Nalón case study'
# Create a new DataFrame with rows that have the target value in the specified column
df1_spain_18 = df_18s[df_18s[column_name] == target_value].copy()
# Reset the index of the new DataFrame
df1_spain_18.reset_index(drop=True, inplace=True)


#Log ratio calaculation
#CO1
CO1_df1_spain= df1_spain[df1_spain['eDNA_no_taxa_marker_CO1'].notna()].copy()
CO1_df1_spain['CO1_df1_spain'] = CO1_df1_spain.apply(lambda row: ln(row['eDNA_no_taxa_marker_CO1'] / row['trad_no_taxa']) if not math.isnan(row['eDNA_no_taxa_marker_CO1']) and not math.isnan(row['trad_no_taxa']) else math.nan, axis=1)

#18S
df_18S_spain = df1_spain_18[df1_spain_18['eDNA no. taxa_ maker 18S'].notna()].copy()
df_18S_spain['df_18S_spain'] = df_18S_spain.apply(lambda row: ln(row['eDNA no. taxa_ maker 18S'] / row['Trad no. taxa']) if not math.isnan(row['eDNA no. taxa_ maker 18S']) and not math.isnan(row['Trad no. taxa']) else math.nan, axis=1)

##plot
#co1
sns.displot(CO1_df1_spain, x="CO1_df1_spain", kind="kde", fill = 'black', alpha=.5)
plt.title("Spain 1 CO1", loc = 'left')
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.tight_layout()
plt.show()
#plt.savefig('CO1_df_spain', dpi=330)

#18s
sns.displot(df_18S_spain, x="df_18S_spain", kind="kde", fill = 'black', alpha=.5)
plt.title("Spain 1 18S", loc = 'left')
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.tight_layout()
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.show()
#plt.savefig('18s_df_spain', dpi=330)


##df2 - Canada
column_name = 'title'
target_value = 'Assessment of stream macroinvertebrate communities with eDNA is not congruent with tissue-based metabarcoding'
# Create a new DataFrame with rows that have the target value in the specified column
df1_canada = df[df[column_name] == target_value].copy()
# Reset the index of the new DataFrame
df1_canada.reset_index(drop=True, inplace=True)

#Log ratio calaculation
CO1_df1_canada = df1_canada[df1_canada['eDNA_no_taxa_marker_CO1'].notna()].copy()
CO1_df1_canada['CO1_df1_canada'] = CO1_df1_canada.apply(lambda row: ln(row['eDNA_no_taxa_marker_CO1'] / row['trad_no_taxa']) if not math.isnan(row['eDNA_no_taxa_marker_CO1']) and not math.isnan(row['trad_no_taxa']) else math.nan, axis=1)

#plot
sns.displot(CO1_df1_canada, x="CO1_df1_canada", kind="kde", fill = 'black', alpha=.5)
plt.title("Canada 1 CO1", loc = 'left')
plt.xlabel("Log-ratio")
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.ylabel("Density")
plt.tight_layout()
plt.show()
#plt.savefig('CO1_canada_df_density', dpi=330)


##df 3 -South Africa
target_value = 'Metabarcoding unsorted kick-samples facilitates macroinvertebrate-based biomonitoring with increased taxonomic resolution, while outperforming environmental DNA'
# Create a new DataFrame with rows that have the target value in the specified column
df1_south_africa = df[df[column_name] == target_value].copy()
# Reset the index of the new DataFrame
df1_south_africa.reset_index(drop=True, inplace=True)

#Log ratio calaculation
CO1_df1_south_africa = df1_south_africa[df1_south_africa['eDNA_no_taxa_marker_CO1'].notna()].copy()
CO1_df1_south_africa['CO1_df1_south_africa'] = CO1_df1_south_africa.apply(lambda row: ln(row['eDNA_no_taxa_marker_CO1'] / row['trad_no_taxa']) if not math.isnan(row['eDNA_no_taxa_marker_CO1']) and not math.isnan(row['trad_no_taxa']) else math.nan, axis=1)

#plot
sns.displot(CO1_df1_south_africa, x="CO1_df1_south_africa", kind="kde", fill = 'black', alpha=.5)
plt.title("South Africa 1 co1", loc = 'left')
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.tight_layout()
plt.show()
#plt.savefig('CO1_southafrica_density', dpi=330)


##df 4 - greece

target_value = 'Assessment of hydrological barriers effect in river benthic fauna coupled with eDNA metabarcoding monitoring'
# Create a new DataFrame with rows that have the target value in the specified column
df1_greece = df[df[column_name] == target_value].copy()
# Reset the index of the new DataFrame
df1_greece.reset_index(drop=True, inplace=True)

#Log ratio calaculation
CO1_df1_greece = df1_greece[df1_greece['eDNA_no_taxa_marker_CO1'].notna()].copy()
CO1_df1_greece['CO1_df1_greece'] = CO1_df1_greece.apply(lambda row: ln(row['eDNA_no_taxa_marker_CO1'] / row['trad_no_taxa']) if not math.isnan(row['eDNA_no_taxa_marker_CO1']) and not math.isnan(row['trad_no_taxa']) else math.nan, axis=1)

#plot
sns.displot(CO1_df1_greece, x="CO1_df1_greece", kind="kde", fill = 'black', alpha=.5)
plt.title("Greece co1", loc = 'left')
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.tight_layout()
plt.show()
#plt.savefig('CO1_greece_density', dpi=330)


##df 6 - Belarus
column_name = 'Title'
target_value = 'Environmental DNA (eDNA) metabarcoding surveys show evidence of non-indigenous freshwater species invasion to new parts of Eastern Europe'
# Create a new DataFrame with rows that have the target value in the specified column
df1_belarus = df_16s[df_16s[column_name] == target_value].copy()
# Reset the index of the new DataFrame
df1_belarus.reset_index(drop=True, inplace=True)

#df1 - belarus: Log ratio calaculation
df1_belarus_16s= df1_belarus[df1_belarus['eDNA no. taxa_marker 16S'].notna()].copy()
df1_belarus_16s['df1_belarus_16s'] = df1_belarus_16s.apply(lambda row: ln(row['eDNA no. taxa_marker 16S'] / row['Trad no. taxa']) if not math.isnan(row['eDNA no. taxa_marker 16S']) and not math.isnan(row['Trad no. taxa']) else math.nan, axis=1)

#plot
sns.displot(df1_belarus_16s, x="df1_belarus_16s", kind="kde", fill = 'black', alpha=.5)
plt.title("Belarus 1 co1", loc = 'left')
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.tight_layout()
plt.show()
#plt.savefig('CO1_belarus_density', dpi=330)


##df 7 - spain
column-name= 'title'
target_value = 'How can eDNA contribute in riverine macroinvertebrate assessment? A metabarcoding approach in the Nalón River (Asturias, Northern Spain)'
# Create a new DataFrame with rows that have the target value in the specified column
df1_spain = df[df[column_name] == target_value].copy()
# Reset the index of the new DataFrame
df1_spain.reset_index(drop=True, inplace=True)

#Log ratio calaculation
CO1_df1_spain= df1_spain[df1_spain['eDNA_no_taxa_marker_CO1'].notna()].copy()
CO1_df1_spain['CO1_df1_spain'] = CO1_df1_spain.apply(lambda row: ln(row['eDNA_no_taxa_marker_CO1'] / row['trad_no_taxa']) if not math.isnan(row['eDNA_no_taxa_marker_CO1']) and not math.isnan(row['trad_no_taxa']) else math.nan, axis=1)

#plot
sns.displot(CO1_df1_spain, x="CO1_df1_spain", kind="kde", fill = 'black', alpha=.5)
plt.title("Spain 2 co1", loc = 'left')
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.tight_layout()
plt.show()
plt.savefig('CO1_spain_density', dpi=330)

##df 8 - SWITZERLAND 2
column_name = 'title'
target_value = 'Monitoring invasive alien macroinvertebrate species with environmental DNA'
# Create a new DataFrame with rows that have the target value in the specified column
df1_switz2 = df[df[column_name] == target_value].copy()
# Reset the index of the new DataFrame
df1_switz2.reset_index(drop=True, inplace=True)

#df1 - switz2: Log ratio calaculation
CO1_df1_switz2= df1_switz2[df1_switz2['eDNA_no_taxa_marker_CO1'].notna()].copy()
CO1_df1_switz2['CO1_df1_switz2'] = CO1_df1_switz2.apply(lambda row: ln(row['eDNA_no_taxa_marker_CO1'] / row['trad_no_taxa']) if not math.isnan(row['eDNA_no_taxa_marker_CO1']) and not math.isnan(row['trad_no_taxa']) else math.nan, axis=1)

#plot
sns.displot(CO1_df1_switz2, x="CO1_df1_switz2", kind="kde", fill = 'black', alpha=.5)
plt.title("Switzerland 2 co1", loc = 'left')
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.tight_layout()
plt.show()
#plt.savefig('CO1_switz2_density', dpi=330)


##df 10 - japan
target_value = 'Aquatic insect community structure revealed by eDNA metabarcoding derives indices for environmental assessment'
# Create a new DataFrame with rows that have the target value in the specified column
df1_japan = df[df[column_name] == target_value].copy()
# Reset the index of the new DataFrame
df1_japan.reset_index(drop=True, inplace=True)

#Log ratio calaculation
CO1_df1_japan= df1_japan[df1_japan['eDNA_no_taxa_marker_CO1'].notna()].copy()
CO1_df1_japan['CO1_df1_japan'] = CO1_df1_japan.apply(lambda row: ln(row['eDNA_no_taxa_marker_CO1'] / row['trad_no_taxa']) if not math.isnan(row['eDNA_no_taxa_marker_CO1']) and not math.isnan(row['trad_no_taxa']) else math.nan, axis=1)

#plot
sns.displot(CO1_df1_japan, x="CO1_df1_japan", kind="kde", fill = 'black', alpha=.5)
plt.title("Japan 1 co1", loc = 'left')
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.tight_layout()
plt.show()
#plt.savefig('CO1_japan_density', dpi=330)


##df 11- switz 3

target_value = 'Utility of environmental DNA for monitoring rare and indicator macroinvertebrate species'
# Create a new DataFrame with rows that have the target value in the specified column
df1_swiz3 = df[df[column_name] == target_value].copy()
# Reset the index of the new DataFrame
df1_swiz3.reset_index(drop=True, inplace=True)

#Log ratio calaculation
CO1_df1_swiz3= df1_swiz3[df1_swiz3['eDNA_no_taxa_marker_CO1'].notna()].copy()
CO1_df1_swiz3['CO1_df1_swiz3'] = CO1_df1_swiz3.apply(lambda row: ln(row['eDNA_no_taxa_marker_CO1'] / row['trad_no_taxa']) if not math.isnan(row['eDNA_no_taxa_marker_CO1']) and not math.isnan(row['trad_no_taxa']) else math.nan, axis=1)

#plot
sns.displot(CO1_df1_swiz3, x="CO1_df1_swiz3", kind="kde", fill = 'black', alpha=.5)
plt.title("Switzerland 3 co1", loc = 'left')
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.tight_layout()
plt.show()
#plt.savefig('CO1_switz3_density', dpi=330)


##df 12 - switz 4

target_value = 'Assessing different components of diversity across a river network using eDNA'
# Create a new DataFrame with rows that have the target value in the specified column
df1_switz4 = df[df[column_name] == target_value].copy()
# Reset the index of the new DataFrame
df1_switz4.reset_index(drop=True, inplace=True)

#Log ratio calaculation
CO1_df1_switz4= df1_switz4[df1_switz4['eDNA_no_taxa_marker_CO1'].notna()].copy()
CO1_df1_switz4['CO1_df1_switz4'] = CO1_df1_switz4.apply(lambda row: ln(row['eDNA_no_taxa_marker_CO1'] / row['trad_no_taxa']) if not math.isnan(row['eDNA_no_taxa_marker_CO1']) and not math.isnan(row['trad_no_taxa']) else math.nan, axis=1)

#plot
sns.displot(CO1_df1_switz4, x="CO1_df1_switz4", kind="kde", fill = 'black', alpha=.5)
plt.title("Switzerland 4 co1", loc = 'left')
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.tight_layout()
plt.show()
#plt.savefig('CO1_switz4_density', dpi=330)


##df 13 - wales
column_name = 'title'
target_value = 'Environmental DNA provides higher resolution assessment of riverine biodiversity and ecosystem function via spatio-temporal nestedness and turnover partitioning'
# Create a new DataFrame with rows that have the target value in the specified column
df1_wales = df[df[column_name] == target_value].copy()
# Reset the index of the new DataFrame
df1_wales.reset_index(drop=True, inplace=True)

#Log ratio calaculation
CO1_df1_wales= df1_wales[df1_wales['eDNA_no_taxa_marker_CO1'].notna()].copy()
CO1_df1_wales['CO1_df1_wales'] = CO1_df1_wales.apply(lambda row: ln(row['eDNA_no_taxa_marker_CO1'] / row['trad_no_taxa']) if not math.isnan(row['eDNA_no_taxa_marker_CO1']) and not math.isnan(row['trad_no_taxa']) else math.nan, axis=1)

#plot
sns.displot(CO1_df1_wales, x="CO1_df1_wales", kind="kde", fill = 'black', alpha=.5)
plt.title("Wales 1 co1", loc = 'left')
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.tight_layout()
plt.show()
#plt.savefig('CO1_wales_density', dpi=330)
################################################################################

#SECTION 4: log ratio for primers of CO1 marker
###For each marker of the eDNA primer CO1 we calculate the log-ratio. Here we find out which marker is the most effective at detecting riverine invertebrates.

# Calculate Ln(eDNA no. taxa_marker CO1 / Trad no. taxa) and create a new column

##Marker 1
BF2_BR2_co1_for_log['eDNA no. taxa_marker CO1'] = BF2_BR2_co1_for_log['eDNA no. taxa_marker CO1'].replace(0, 0.001)
BF2_CO1_df = BF2_BR2_co1_for_log[BF2_BR2_co1_for_log['eDNA no. taxa_marker CO1'].notna()].copy()
BF2_CO1_df['Ln_BF2_CO1'] = BF2_CO1_df.apply(lambda row: ln(row['eDNA no. taxa_marker CO1'] / row['Trad no. taxa']) if not math.isnan(row['eDNA no. taxa_marker CO1']) and not math.isnan(row['Trad no. taxa']) else math.nan, axis=1)

##Marker 2
fwhF2_EPTDr2n_co1_for_log['eDNA no. taxa_marker CO1'] = fwhF2_EPTDr2n_co1_for_log['eDNA no. taxa_marker CO1'].replace(0, 0.001)
fwhF2_CO1_df = fwhF2_EPTDr2n_co1_for_log[fwhF2_EPTDr2n_co1_for_log['eDNA no. taxa_marker CO1'].notna()].copy()
fwhF2_CO1_df['Ln_fwhF2_CO1'] = fwhF2_CO1_df.apply(lambda row: ln(row['eDNA no. taxa_marker CO1'] / row['Trad no. taxa']) if not math.isnan(row['eDNA no. taxa_marker CO1']) and not math.isnan(row['Trad no. taxa']) else math.nan, axis=1)

##Marker 3
m1COIintF_jgHCO2198_co1_for_log['eDNA no. taxa_marker CO1'] = m1COIintF_jgHCO2198_co1_for_log['eDNA no. taxa_marker CO1'].replace(0, 0.001)
m1COIintF_CO1_df = m1COIintF_jgHCO2198_co1_for_log[m1COIintF_jgHCO2198_co1_for_log['eDNA no. taxa_marker CO1'].notna()].copy()
m1COIintF_CO1_df['Ln_m1COIintF_CO1'] = m1COIintF_CO1_df.apply(lambda row: ln(row['eDNA no. taxa_marker CO1'] / row['Trad no. taxa']) if not math.isnan(row['eDNA no. taxa_marker CO1']) and not math.isnan(row['Trad no. taxa']) else math.nan, axis=1)

##Marker 4
LCO1490HCO2198_co1_for_log ['eDNA no. taxa_marker CO1'] = LCO1490HCO2198_co1_for_log ['eDNA no. taxa_marker CO1'].replace(0, 0.001)
LCO1490HCO2198_CO1_df = LCO1490HCO2198_co1_for_log [LCO1490HCO2198_co1_for_log ['eDNA no. taxa_marker CO1'].notna()].copy()
LCO1490HCO2198_CO1_df['LCO1490HCO2198_CO1'] = LCO1490HCO2198_CO1_df.apply(lambda row: ln(row['eDNA no. taxa_marker CO1'] / row['Trad no. taxa']) if not math.isnan(row['eDNA no. taxa_marker CO1']) and not math.isnan(row['Trad no. taxa']) else math.nan, axis=1)


#Plot log-ratios per marker

##LCO1490HCO2198 primer
stat, p_value = shapiro(LCO1490HCO2198_CO1_df['LCO1490HCO2198_CO1']) #check normality
print(f'Statistics={stat:.4f}, p-value={p_value:.4f}')

sns.set_style("whitegrid")
sns.set_theme()
plt.figure()
sns.displot(LCO1490HCO2198_CO1_df, x="LCO1490HCO2198_CO1", kind="kde", fill='black', alpha=.5)
plt.xticks(fontsize=14)  # Adjust font size as needed for x-axis
plt.yticks(fontsize=14)
plt.title("Density Plot CO1 LCO1490HCO2198", loc='left')
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.xlabel("Log-ratio", fontsize=16)  # Adjust font size as needed
plt.ylabel("Density", fontsize=16)  # Adjust font size as needed
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.tight_layout()
#plt.savefig('Density Plot CO1 LCO1490HCO2198.tiff',dpi = 330)
plt.show()

average_18S = LCO1490HCO2198_CO1_df['LCO1490HCO2198_CO1'].mean(skipna=True)  #average log-ratio
print(f"The average is: {average_18S}")


##BF2 primer
stat, p_value = shapiro(BF2_CO1_df['Ln_BF2_CO1']) #check normality
print(f'Statistics={stat:.4f}, p-value={p_value:.4f}')

sns.set_style("whitegrid")
sns.set_theme()
plt.figure()
sns.displot(BF2_CO1_df, x="Ln_BF2_CO1", kind="kde", fill='black', alpha=.5)
plt.xticks(fontsize=14)  # Adjust font size as needed for x-axis
plt.yticks(fontsize=14)
plt.title("Density Plot CO1 BF2", loc='left')
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.xlabel("Log-ratio", fontsize=16)  # Adjust font size as needed
plt.ylabel("Density", fontsize=16)  # Adjust font size as needed
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.tight_layout()
#plt.savefig('Density Plot CO1 BF2.tiff',dpi = 330)
plt.show()

average_18S = BF2_CO1_df['Ln_BF2_CO1'].mean(skipna=True)  #average log-ratio
print(f"The average is: {average_18S}")



##FwhF2 primer
stat, p_value = shapiro(fwhF2_CO1_df['Ln_fwhF2_CO1']) #check normality
print(f'Statistics={stat:.4f}, p-value={p_value:.4f}')

sns.set_style("whitegrid")
sns.set_theme()
plt.figure()
sns.displot(fwhF2_CO1_df, x="Ln_fwhF2_CO1", kind="kde", fill='black', alpha=.5)
plt.xticks(fontsize=14)  # Adjust font size as needed for x-axis
plt.yticks(fontsize=14)
plt.title("Density Plot CO1 fwhF2", loc='left')
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.xlabel("Log-ratio", fontsize=16)  # Adjust font size as needed
plt.ylabel("Density", fontsize=16) 
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.tight_layout()
plt.savefig('Density Plot CO1 fwhF2.tiff',dpi = 330)
plt.show()

average_18S = fwhF2_CO1_df['Ln_fwhF2_CO1'].mean(skipna=True) #average log-ratio
print(f"The average is: {average_18S}")



##m1CO1intF primer
stat, p_value = shapiro(m1COIintF_CO1_df['Ln_m1COIintF_CO1']) #check normality
print(f'Statistics={stat:.4f}, p-value={p_value:.4f}')
sns.histplot(merged_df['Ln_m1COIintF_CO1'], kde=True)
plt.show()

sns.set_style("whitegrid")
sns.set_theme()
plt.figure()
#sns.displot(m1COIintF_CO1_df, x="Ln_m1COIintF_CO1", kind="kde", fill='black', alpha=.5, hue="Taxa level")
sns.displot(m1COIintF_CO1_df, x="Ln_m1COIintF_CO1", kind="kde", fill='black', alpha=.5)
plt.xticks(fontsize=14)  # Adjust font size as needed for x-axis
plt.yticks(fontsize=14)
plt.title("Density Plot CO1 m1COIintF", loc='left')
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.xlabel("Log-ratio", fontsize=16)  # Adjust font size as needed
plt.ylabel("Density", fontsize=16) 
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.tight_layout()
#plt.savefig('Density Plot CO1 m1COIintF.tiff',dpi = 330)
plt.show()

average_18S = m1COIintF_CO1_df['Ln_m1COIintF_CO1'].mean(skipna=True)  #average log-ratio
print(f"The average is: {average_18S}")




#Determine if there is a significant difference between the mean log-ratios of the 4 eDNA markers of CO1 primer
###the four primer groups log-ratios are normally distrubuted (p values >0.05)

#ANOVA (ONE WAY) - test difference in means among groups
statistic, p_value = stats.f_oneway(BF2_CO1_df['Ln_BF2_CO1'], fwhF2_CO1_df['Ln_fwhF2_CO1'], m1COIintF_CO1_df['Ln_m1COIintF_CO1', LCO1490HCO2198_CO1_df['LCO1490HCO2198_CO1']]) 
# Print the results
df_between = 4 - 1  # k - 1, where k is the number of groups (4 in this case)
df_within = len(BF2_CO1_df) + len(fwhF2_CO1_df) + len(m1COIintF_CO1_df) + len(LCO1490HCO2198_CO1_df) - 4  # N (number of observations) - k (number of groups)
# Print the results
print(f'ANOVA F-statistic: {statistic:.4f}')
print(f'Degrees of freedom (between, within): {df_between}, {df_within}')
print(f'P-value: {p_value:.4f}'

# Check the significance level
alpha = 0.05
if p_value < alpha:
    print("Reject the null hypothesis. There are significant differences between groups.")
else:
        print("Fail to reject the null hypothesis. No significant differences between groups.") #not significant



################################################################################

#SECTION 5: Log-ratio for terrestrial markers

#1 - calculate log for each marker

df_1terres = df_terrestrial.drop(['Sample no.'], axis=1)#remove column

# Calculate Ln(eDNA no. taxa_marker CO1 / Trad no. taxa) and create a new column
df_1terres['eDNA no. taxa_marker CO1'] = df_1terres['eDNA no. taxa_marker CO1'].replace(0, 0.001)
CO1_df = df_1terres[df_1terres['eDNA no. taxa_marker CO1'].notna()].copy()
CO1_df['Ln_CO1'] = CO1_df.apply(lambda row: ln(row['eDNA no. taxa_marker CO1'] / row['Trad no. taxa']) if not math.isnan(row['eDNA no. taxa_marker CO1']) and not math.isnan(row['Trad no. taxa']) else math.nan, axis=1)

# Calculate Ln(eDNA no. taxa_marker 16S / Trad no. taxa) and create a new column
df_16S = df_1terres[df_1terres['eDNA no. taxa_marker 16S'].notna()].copy()
df_16S['Ln_16S'] = df_16S.apply(lambda row: math.log(row['eDNA no. taxa_marker 16S'] / row['Trad no. taxa']) if not math.isnan(row['eDNA no. taxa_marker 16S']) and not math.isnan(row['Trad no. taxa']) and row['Trad no. taxa'] != 0 else math.nan, axis=1)


#Density graphs
##CO1
sns.displot(CO1_df, x="Ln_CO1", kind="kde", fill = 'black', alpha=.5)
plt.title("density plot co1", loc = 'left')
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.tight_layout()
#plt.savefig('Density plot CO1 TERRESTRIAL.tiff',dpi = 330)
plt.show()
###AVERGAE value is:
CO1_average = CO1_df.loc[np.isfinite(CO1_df['Ln_CO1']), 'Ln_CO1'].mean()
print(f"The average is: {CO1_average}")

#Desity plot per group 
sns.displot(CO1_df, x="Ln_CO1", kind="kde", fill='black', alpha=.5, hue="Taxa level")
plt.title("Density Plot CO1", loc='left')
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.show()
#plt.savefig('Density plot CO1 TERRESTRIAL per group.tiff',dpi = 330)


##16S 
sns.displot(df_16S, x="Ln_16S", kind="kde", fill = 'black', alpha=.5)
plt.title("density plot 16s", loc = 'left')
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.tight_layout()
#plt.savefig('Density plot 16s terrestrial.tiff',dpi = 330) 
plt.show()
###AVERGAE value is:
average_16S = df_16S['Ln_16S'].mean(skipna=True)
print(f"The average is: {average_16S}")

#Desity plot per group
 sns.displot(df_16S, x="Ln_16S", kind="kde", fill='black', alpha=.5, hue="Taxa level")
plt.title("Density Plot CO1", loc='left')
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.show()
#plt.savefig('Density plot 16s terrestrial per group.tiff',dpi = 330)


################################################################################

#SECTION 6: Identifying missing terrestrial taxa


#eDNA
CO1_taxa = CO1_terrestrial[['Title','Comparison', 'Site no.', 'Sample no.', 'CO1 marker_ eDNA taxa (primer 1)']]
CO1_taxa = CO1_taxa.dropna(subset=['CO1 marker_ eDNA taxa (primer 1)'])
CO1_taxa1 = CO1_taxa.copy()
#CO1_taxa1.loc[:, 'method'] = 'CO1'

S16_taxa = S16_terrestrial[['Title','Comparison', 'Site no.', 'Sample no.', '16S marker_ eDNA taxa (primer 1)']]
S16_taxa = S16_taxa.dropna(subset=['16S marker_ eDNA taxa (primer 1)'])
S16_taxa1 = S16_taxa.copy()
#S16_taxa1.loc[:, 'method'] = '16S'

trad_taxa_CO1 = CO1_terrestrial[['Title','Comparison', 'Site no.', 'Sample no.', 'Trad no. taxa']]
trad_taxa_CO1 = trad_taxa_CO1.dropna(subset=['Trad no. taxa'])
trad_taxa_CO2 = trad_taxa_CO1.copy()
#trad_taxa_CO2.loc[:, 'method'] = 'Traditional'

trad_taxa_16S1 = S16_terrestrial[['Title','Comparison', 'Site no.', 'Sample no.', 'Trad no. taxa']]
trad_taxa_16S1 = trad_taxa_16S1.dropna(subset=['Trad no. taxa'])
trad_taxa_16S2 = trad_taxa_16S1.copy()
#trad_taxa_16S2.loc[:, 'method'] = 'Traditional'


##co1
CO1_taxa1_diff = CO1_taxa1[['CO1 marker_ eDNA taxa (primer 1)', 'Comparison', 'Site no.', 'Sample no.', 'Title']].copy()
CO1_taxa1_diff.loc[:, 'combin'] = 'co1'
CO1_taxa1_diff_cleaned = CO1_taxa1_diff.dropna().copy()
CO1_taxa1_diff_cleaned = CO1_taxa1_diff_cleaned.rename(columns={'CO1 marker_ eDNA taxa (primer 1)': 'taxa'})

trad_taxa_CO2_diff = trad_taxa_CO2[['Trad no. taxa', 'Comparison', 'Site no.', 'Sample no.', 'Title']]
trad_taxa_CO2_diff = trad_taxa_CO2_diff.copy()  # Ensure you are working with a copy
trad_taxa_CO2_diff.loc[:, 'combin'] = 'trad'
trad_taxa_CO2_diff_cleaned = trad_taxa_CO2_diff.dropna()
trad_taxa_CO2_diff_cleaned = trad_taxa_CO2_diff_cleaned.rename(columns={'Trad no. taxa': 'taxa'})

merged_df1 = pd.merge(
    CO1_taxa1_diff_cleaned,
    trad_taxa_CO2_diff_cleaned,
    on=['taxa', 'Comparison', 'Site no.', 'Sample no.', 'Title'],
    how='outer',
    suffixes=('_df1', '_df2')
)


##16s
S16_taxa1_diff = S16_taxa1[['16S marker_ eDNA taxa (primer 1)',  'Comparison', 'Site no.', 'Sample no.', 'Title']]
S16_taxa1_diff = S16_taxa1_diff.copy()  # Ensure you are working with a copy
S16_taxa1_diff.loc[:, 'combin'] = '16s'
S16_taxa1_diff_cleaned = S16_taxa1_diff.dropna()
S16_taxa1_diff_cleaned = S16_taxa1_diff_cleaned.rename(columns={'16S marker_ eDNA taxa (primer 1)': 'taxa'})

trad_taxa_16S1_diff = trad_taxa_16S1 [['Trad no. taxa',  'Comparison', 'Site no.', 'Sample no.', 'Title']]
trad_taxa_16S1_diff = trad_taxa_16S1_diff.copy()  # Ensure you are working with a copy
trad_taxa_16S1_diff.loc[:, 'combin'] = 'trad'
trad_taxa_16S1_diff_cleaned = trad_taxa_16S1_diff.dropna()
trad_taxa_16S1_diff_cleaned = trad_taxa_16S1_diff_cleaned.rename(columns={'Trad no. taxa': 'taxa'})

# Merge the DataFrames with indicator=True
merged_df2 = pd.merge(
    S16_taxa1_diff_cleaned,
    trad_taxa_16S1_diff_cleaned,
    on=['taxa', 'Comparison', 'Site no.', 'Sample no.', 'Title'],
    how='outer',
    suffixes=('_df1', '_df2')
)


#UNIQUE TAXA 
#In EXCEL use unique funtion to identify the number of unique taxa detected by trad sampling overall and the number of these detected by eDNA


################################################################################


##SECTION 7: log ratios per study


#Log for each study/area

##Set plot parameters
sns.set(rc='axes', labelsize=16, titlesize=18)
sns.set(rc='xtick', labelsize=14)
sns.set(rc='ytick', labelsize=14)
sns.set(rc='legend', fontsize=14)
sns.set(rc={'axes.labelsize': 16, 'axes.titlesize': 18, 'xtick.labelsize': 14, 'ytick.labelsize': 14, 'legend.fontsize': 14})


#extract dataframes by title of study
column_name = 'title'

##New zealand 1 (Soil)
target_value = 'Detecting invertebrate ecosystem service providers in orchards: traditional methods versus barcoding of environmental DNA in soil'
# Create a new DataFrame with rows that have the target value in the specified column
df1_newzeal1 = df_terrestrial[df_terrestrial[column_name] == target_value].copy()
# Reset the index of the new DataFrame
df1_newzeal1.reset_index(drop=True, inplace=True)

#Log ratio calaculation
df1_newzeal1['eDNA no. taxa_marker CO1'] = df1_newzeal1['eDNA no. taxa_marker CO1'].replace(0, 0.001)
CO1_df1_newzeal1 = df1_newzeal1[df1_newzeal1['eDNA no. taxa_marker CO1'].notna()].copy()
CO1_df1_newzeal1['CO1_df1_newzeal1'] = CO1_df1_newzeal1.apply(lambda row: ln(row['eDNA no. taxa_marker CO1'] / row['Trad no. taxa']) if not math.isnan(row['eDNA no. taxa_marker CO1']) and not math.isnan(row['Trad no. taxa']) else math.nan, axis=1)

#plot
sns.displot(CO1_df1_newzeal1, x="CO1_df1_newzeal1", kind="kde", fill = 'black', alpha=.5)
plt.title("New Zealand (1) co1 - soil", loc = 'left')
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.tight_layout()
plt.show()
#plt.savefig('CO1_df1_newzeal1_terr', dpi=330)


##New zealand 2 (Soil)
target_value = 'DNA metabarcoding as a tool for invertebrate community monitoring: a case study comparison with conventional techniques'
# Create a new DataFrame with rows that have the target value in the specified column
df1_newzeal1 = df_terrestrial[df_terrestrial[column_name] == target_value].copy()
# Reset the index of the new DataFrame
df1_newzeal1.reset_index(drop=True, inplace=True)

#Log ratio calaculation
CO1_df1_newzeal1 = df1_newzeal1[df1_newzeal1['eDNA no. taxa_marker CO1'].notna()].copy()
CO1_df1_newzeal1['CO1_df1_newzeal1'] = CO1_df1_newzeal1.apply(lambda row: ln(row['eDNA no. taxa_marker CO1'] / row['Trad no. taxa']) if not math.isnan(row['eDNA no. taxa_marker CO1']) and not math.isnan(row['Trad no. taxa']) else math.nan, axis=1)

#plot
sns.displot(CO1_df1_newzeal1, x="CO1_df1_newzeal1", kind="kde", fill = 'black', alpha=.5)
plt.title("New Zealand (2) co1 - soil", loc = 'left')
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.tight_layout()
plt.show()
#plt.savefig('CO1_df1_newzeal2_terr', dpi=330)


##Germany (rainwater)
target_value = 'Its raining species Rainwash eDNA metabarcoding as a minimally invasive method to assess tree canopy invertebrate diversity'
# Create a new DataFrame with rows that have the target value in the specified column
df1_newzeal1 = df_terrestrial[df_terrestrial[column_name] == target_value].copy()
# Reset the index of the new DataFrame
df1_newzeal1.reset_index(drop=True, inplace=True)

#Log ratio calaculation
CO1_df1_newzeal1 = df1_newzeal1[df1_newzeal1['eDNA no. taxa_marker CO1'].notna()].copy()
CO1_df1_newzeal1['CO1_df1_newzeal1'] = CO1_df1_newzeal1.apply(lambda row: ln(row['eDNA no. taxa_marker CO1'] / row['Trad no. taxa']) if not math.isnan(row['eDNA no. taxa_marker CO1']) and not math.isnan(row['Trad no. taxa']) else math.nan, axis=1)

#plot
plt.figure()
sns.displot(CO1_df1_newzeal1, x="CO1_df1_newzeal1", kind="kde", fill = 'black', alpha=.5)
plt.title("Germany co1 - rainwater", loc = 'left')
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.tight_layout()
plt.show()
#plt.savefig('CO1_df1_Germany1_terr', dpi=330)


##China (soil)
target_value = 'Diversity of soil faunal community as influenced by crop straw combined with different synthetic fertilizers in upland purple soil'
# Create a new DataFrame with rows that have the target value in the specified column
df1_newzeal1 = df_terrestrial[df_terrestrial[column_name] == target_value].copy()
# Reset the index of the new DataFrame
df1_newzeal1.reset_index(drop=True, inplace=True)

#Log ratio calaculation
CO1_df1_newzeal1 = df1_newzeal1[df1_newzeal1['eDNA no. taxa_marker CO1'].notna()].copy()
CO1_df1_newzeal1['CO1_df1_newzeal1'] = CO1_df1_newzeal1.apply(lambda row: ln(row['eDNA no. taxa_marker CO1'] / row['Trad no. taxa']) if not math.isnan(row['eDNA no. taxa_marker CO1']) and not math.isnan(row['Trad no. taxa']) else math.nan, axis=1)

#plot
sns.displot(CO1_df1_newzeal1, x="CO1_df1_newzeal1", kind="kde", fill = 'black', alpha=.5)
plt.title("China co1 - soil", loc = 'left')
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.tight_layout()
plt.show()
#plt.savefig('CO1_df1_China1_terr', dpi=330)


##England (soil)
target_value = 'Environmental DNA is more effective than soil-pit hand sorting in evaluating earthworm biodiversity responses to more regenerative agricultural management'
# Create a new DataFrame with rows that have the target value in the specified column
df1_newzeal1 = df_terrestrial[df_terrestrial[column_name] == target_value].copy()
# Reset the index of the new DataFrame
df1_newzeal1.reset_index(drop=True, inplace=True)

# Log ratio calaculation
df1_newzeal1['eDNA no. taxa_marker 16S'] = df1_newzeal1['eDNA no. taxa_marker 16S'].replace(0, 0.001)
df1_newzeal1['Trad no. taxa'] = df1_newzeal1['Trad no. taxa'].replace(0, 0.001)
CO1_df1_newzeal1 = df1_newzeal1[df1_newzeal1['eDNA no. taxa_marker 16S'].notna()].copy()
CO1_df1_newzeal1['CO1_df1_newzeal1'] = CO1_df1_newzeal1.apply(lambda row: ln(row['eDNA no. taxa_marker 16S'] / row['Trad no. taxa']) if not math.isnan(row['eDNA no. taxa_marker 16S']) and not math.isnan(row['Trad no. taxa']) else math.nan, axis=1)

#plot
sns.displot(CO1_df1_newzeal1, x="CO1_df1_newzeal1", kind="kde", fill = 'black', alpha=.5)
plt.title("England 16s - soil", loc = 'left')
plt.xlabel("Log-ratio")
plt.ylabel("Density")
plt.axvline(x=0, color='red', linestyle='--', linewidth=1.5)
plt.tight_layout()
plt.show()
#plt.savefig('CO1_df1_england1_terr', dpi=330)


################################################################################

#making binomial dataset

merged_df = pd.merge(merged_df_co1_corr, all_taxa_present_corr, on='taxa', how='outer')
merged_df.to_excel('all_co1_taxa_withstudyID.xlsx', index=False)

merged_df = pd.merge(merged_18S_df_corr, all_taxa_present_corr, on='taxa', how='outer')
merged_df.to_excel('all_18s_taxa_withstudyID.xlsx', index=False)

merged_df = pd.merge(merged_16S_df_corr, all_taxa_present_corr, on='taxa', how='outer')
merged_df.to_excel('all_16s_taxa_withstudyID.xlsx', index=False)


merged_df = pd.merge(merged_16S_df_terr, all_taxa_present_corr_terr, on='taxa', how='outer')
merged_df.to_excel('all_16s_taxa_withstudyID_terr.xlsx', index=False)

merged_df = pd.merge(merged_co1_df_terr, all_taxa_present_corr_terr, on='taxa', how='outer')
merged_df.to_excel('all_co1_taxa_withstudyID_terr.xlsx', index=False)

