# Semester Project II: Lipidomics Data Analysis

Overview of all R-scripts used during the project.  
To import all build in functions use the following code:  
 ################################################################
  '#  Libraries  
install.packages("devtools")  
library(devtools)  
 ################################################################################
  '#  source of all used functions  
source_url("https://raw.githubusercontent.com/AbrahamProject/Lipidomics_Project/main/Data_analysis_function.R" )  
source_url("https://raw.githubusercontent.com/AbrahamProject/Lipidomics_Project/main/Processing_function.R" )  
 ################################################################################


# I) Raw Data Analysis
### Raw_data_analysis_xcms.R
# II) MS-dial Data Analysis
## 1) Function Scripts
### a) Processing_function.R
### b) Data_analysis_function.R

## 2) Scripts using the functions
### a) MS_dial_processing_script.R
### b) MS_dial_Feature_script.R
### c) MS_dial_summary_plots_script.R
### d) MS_dial_differential_Analysis_script.R
### e) MS_dial_CID_vs_EAD_script.R

