##
# Janey Lienau
# Question 1: Qualifying Exams
# Analysis of Phillips et al. 2019 "Global distribution of earthworm diversity"
##

#----load data
# this data is downloaded from public repository from this link https://idata.idiv.de/ddm/Data/ShowData/1818?version=15
# the following csv files were downloaded on Monday Feb 24, 2025 at 12:45 pm. 
# 1818_2_sWormModelData_meta-data_correction_2020April_complete.csv
# 1818_2_sWormModelData_correction_2020April.csv

#set the working directory to this repository
setwd("/Users/janeylienau/Desktop/GitHubRepository/q1-earthworm-redo")
#load in data
metadata.df <- read.csv("raw-data/1818_15_1818_2_sWormModelData_meta-data_correction_2020April_complete.csv") 
modeldata.df <- read.csv("raw-data/1818_15_1818_2_sWormModelData_correction_2020April.csv")
#we are 