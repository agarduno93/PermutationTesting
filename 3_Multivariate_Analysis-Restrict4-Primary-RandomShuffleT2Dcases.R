#File Description: This file is the same as the primary analysis, but the sample ids are sampled.

#Multivariate Analysis
# Input
args <- commandArgs(trailingOnly = TRUE)

iteration=args[1]
temp_dir=args[2]


#Import libraries
library(tidyverse)
library(ggpubr)
library(dplyr)
library(lme4)
library(ggplot2)
library(rms)
library(mgcv)
library(tidymv)
library(tidyr)
library("ggthemes")
library(broom)
library(erer)
library(stringr)
library(R.utils)
library(data.table)
library("rqdatatable")

#Restricted analysis (only)

#Read the complete set (no imputation)
UKBB_AG2_m <- fread("/cellar/users/agarduno/jupyter/UKBB_AG2_12Jan21.txt", header = TRUE, na.strings=c("",".","NA")) %>% select(f.eid,T2D_status,GRS_WT_LIR,GRS_WT_IR,GRS_RAW_T2DIR,
                                 GRS_RAW_T2DIR2d2,GRS_WT_IR2d2,GRS_RAW_LIR2d2,GRS_WT_L5E8IR,
                                 GRS_WT_L5E8IRd,GRS_WT_L5E8IR2d2,GRS_WT_L1E5IR,GRS_WT_L1E5IRd,GRS_WT_L1E5IR2d2,ALBUMINERIA.0.0,
                                 GRS_RAW_T2DIRd,GRS_RAW_IR2d,GRS_WT_IR2d,GRS_RAW_LIR2d,GRS_WT_T2DIRd,GRS_WT_LIR2d2,
                                 ESKD.0.0,CKD.0.0,DN.0.0,ALL.0.0,NONESKD.0.0,DNCKD.0.0,
                                 CTRL_DNCKD.0.0,ACR.0.0,EGFR.0.0,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,SEX.0.0,
                                        IDEAL_DIET2.0.0,LIFESCORE,AGE.0.0,
                                        SES_TDI.0.0,BMI.0.0,EDUYEARS,SBP.0.0,HYP_POS1,STATIN,
                                        GRS_WT_LIR2d2s,GRS_WT_IRd2s,
                                        GRS_RAW_T2DIR2d2s,GRS_WT_L5E8IR2d2s,
                                        GRS_WT_L1E5IR2d2s,GRS_WT_IRf3,GRS_RAW_T2DIRf3,WHR.0.0,GRS_WT_IR53f3)
UKBB_AG2=as.data.frame(UKBB_AG2_m)
dim(UKBB_AG2)
rm(UKBB_AG2_m)

#Dichotomize Outcomes for Logistic Regression

#1_CKD
UKBB_AG2$CKD_only.0.0 <- factor(ifelse(UKBB_AG2$CKD.0.0=="CKD controls","CKD controls",
                              ifelse(UKBB_AG2$CKD.0.0=="CKD","CKD",NA)),
                levels = c("CKD controls", "CKD"))
#Set the refernece
UKBB_AG2$CKD_only.0.0 <- relevel(UKBB_AG2$CKD_only.0.0, ref = "CKD controls")

#2_CKD Extreme
UKBB_AG2$CKD_ex.0.0 <- factor(ifelse(UKBB_AG2$CKD.0.0=="CKD controls","CKD controls",
                              ifelse(UKBB_AG2$CKD.0.0=="CKD extreme","CKD extreme",NA)),
                levels = c("CKD controls", "CKD extreme"))
#Set the refernece
UKBB_AG2$CKD_ex.0.0 <- relevel(UKBB_AG2$CKD_ex.0.0, ref = "CKD controls")

#3_Micro
UKBB_AG2$micro.0.0 <- factor(ifelse(UKBB_AG2$ALBUMINERIA.0.0=="micro","micro",
                              ifelse(UKBB_AG2$ALBUMINERIA.0.0=="normo","normo",NA)),
                levels = c("normo", "micro"))
#Set the reference
UKBB_AG2$micro.0.0 <- relevel(UKBB_AG2$micro.0.0, ref = "normo")

#4_Macro
UKBB_AG2$macro.0.0 <- factor(ifelse(UKBB_AG2$ALBUMINERIA.0.0=="macro","macro",
                              ifelse(UKBB_AG2$ALBUMINERIA.0.0=="normo","normo",NA)),
                levels = c("normo", "macro"))
#Set the reference
UKBB_AG2$macro.0.0 <- relevel(UKBB_AG2$macro.0.0, ref = "normo")

#5_Macro
UKBB_AG2$macro.0.0 <- factor(ifelse(UKBB_AG2$ALBUMINERIA.0.0=="macro","macro",
                              ifelse(UKBB_AG2$ALBUMINERIA.0.0=="normo","normo",NA)),
                levels = c("normo", "macro"))

#6_ESKD vs. Macro
UKBB_AG2$ESKD_macro.0.0 <- factor(ifelse(UKBB_AG2$ESKD.0.0=="yes","ESKD",
                              ifelse(UKBB_AG2$ALBUMINERIA.0.0=="macro","macro",NA)),
                levels = c("macro","ESKD"))

#7_DNCKD vs. Control DNCKD
UKBB_AG2$DNCKD2.0.0 <- factor(ifelse(UKBB_AG2$DNCKD.0.0=="yes","DNCKD",
                              ifelse(UKBB_AG2$CTRL_DNCKD.0.0=="yes","DNCKD Control",NA)),
                levels = c("DNCKD Control","DNCKD"))

#8_ESKD vs. Normo, Macro, Micro
UKBB_AG2$ESKD_Albu.0.0 <- factor(ifelse(UKBB_AG2$ESKD.0.0=="yes","ESKD",
                              ifelse(UKBB_AG2$ALBUMINERIA.0.0 %in% c("normo","macro","micro"),"albu",NA)),
                levels = c("albu","ESKD"))

#Set the reference
UKBB_AG2$macro.0.0 <- relevel(UKBB_AG2$macro.0.0, ref = "normo")

#Summarize Counts of Disease Outcomes
#table(UKBB_AG2$CKD_only.0.0) #1
#table(UKBB_AG2$CKD_ex.0.0) #2
#table(UKBB_AG2$micro.0.0) #3
#table(UKBB_AG2$macro.0.0) #4
#table(UKBB_AG2$ESKD.0.0) #5
#table(UKBB_AG2$DN.0.0) #6
#table(UKBB_AG2$ALL.0.0) #7
#table(UKBB_AG2$ESKD.0.0) #8
#table(UKBB_AG2$ESKD_macro.0.0) #8
#table(UKBB_AG2$ESKD_Albu.0.0) #9
#table(UKBB_AG2$DNCKD2.0.0) #10

#3x2x1 - T2D x Disease x Exposure

#IR
#table(UKBB_AG2$T2D_status,UKBB_AG2$CKD_only.0.0,UKBB_AG2$GRS_WT_IR2d2) #1
#table(UKBB_AG2$T2D_status,UKBB_AG2$CKD_ex.0.0,UKBB_AG2$GRS_WT_IR2d2) #2
#table(UKBB_AG2$T2D_status,UKBB_AG2$micro.0.0,UKBB_AG2$GRS_WT_IR2d2) #3
#table(UKBB_AG2$T2D_status,UKBB_AG2$macro.0.0,UKBB_AG2$GRS_WT_IR2d2) #4
##table(UKBB_AG2$T2D_status,UKBB_AG2$ESKD.0.0,UKBB_AG2$GRS_WT_IR2d2) #5
#table(UKBB_AG2$T2D_status,UKBB_AG2$DN.0.0,UKBB_AG2$GRS_WT_IR2d2) #6
#table(UKBB_AG2$T2D_status,UKBB_AG2$ALL.0.0,UKBB_AG2$GRS_WT_IR2d2) #7
##table(UKBB_AG2$T2D_status,UKBB_AG2$ESKD.0.0,UKBB_AG2$GRS_WT_IR2d2) #8
##table(UKBB_AG2$T2D_status,UKBB_AG2$ESKD_macro.0.0,UKBB_AG2$GRS_WT_IR2d2) #8
##table(UKBB_AG2$T2D_status,UKBB_AG2$ESKD_Albu.0.0,UKBB_AG2$GRS_WT_IR2d2) #9
#table(UKBB_AG2$T2D_status,UKBB_AG2$DNCKD2.0.0,UKBB_AG2$GRS_WT_IR2d2) #10

#Cross-Tabulations of T2D x Disease x PGS (Overall)
#T2D
#table(UKBB_AG2$T2D_status,UKBB_AG2$CKD_only.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2) #1
#table(UKBB_AG2$T2D_status,UKBB_AG2$CKD_ex.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2) #2
#table(UKBB_AG2$T2D_status,UKBB_AG2$micro.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2) #3
#table(UKBB_AG2$T2D_status,UKBB_AG2$macro.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2) #4
##table(UKBB_AG2$T2D_status,UKBB_AG2$ESKD.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2) #5
#table(UKBB_AG2$T2D_status,UKBB_AG2$DN.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2) #6
#table(UKBB_AG2$T2D_status,UKBB_AG2$ALL.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2) #7
##table(UKBB_AG2$T2D_status,UKBB_AG2$ESKD.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2) #8
##table(UKBB_AG2$T2D_status,UKBB_AG2$ESKD_macro.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2) #8
##table(UKBB_AG2$T2D_status,UKBB_AG2$ESKD_Albu.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2) #9
#table(UKBB_AG2$T2D_status,UKBB_AG2$DNCKD2.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2) #10

#Cross-Tabulations of T2D x Disease x PGS (Stratum-Specific)
#T2D
#table(UKBB_AG2$T2D_status,UKBB_AG2$CKD_only.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2s) #1
#table(UKBB_AG2$T2D_status,UKBB_AG2$CKD_ex.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2s) #2
#table(UKBB_AG2$T2D_status,UKBB_AG2$micro.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2s) #3
#table(UKBB_AG2$T2D_status,UKBB_AG2$macro.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2s) #4
##table(UKBB_AG2$T2D_status,UKBB_AG2$ESKD.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2) #5
#table(UKBB_AG2$T2D_status,UKBB_AG2$DN.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2s) #6
#table(UKBB_AG2$T2D_status,UKBB_AG2$ALL.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2s) #7
##table(UKBB_AG2$T2D_status,UKBB_AG2$ESKD.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2) #8
##table(UKBB_AG2$T2D_status,UKBB_AG2$ESKD_macro.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2) #8
##table(UKBB_AG2$T2D_status,UKBB_AG2$ESKD_Albu.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2) #9
#table(UKBB_AG2$T2D_status,UKBB_AG2$DNCKD2.0.0,UKBB_AG2$GRS_RAW_T2DIR2d2s) #10

#Shuffling, first identify count of cases

#total T2D cases
T2D_cases <- table(UKBB_AG2$T2D_status)[2]
#print('T2D')
#print(T2D_cases)
#micro
T2Dmicro_cases <- table(UKBB_AG2$T2D_status,UKBB_AG2$micro.0.0)[2,2]
#print('micro')
#print(T2Dmicro_cases)
#micro_controls
T2Dmicro_controls <- table(UKBB_AG2$T2D_status,UKBB_AG2$micro.0.0)[2,1]
#print(T2Dmicro_controls)
#CKD
T2DCKD_cases <- table(UKBB_AG2$T2D_status,UKBB_AG2$CKD_only.0.0)[2,2]
#print('CKD')
#print(T2DCKD_cases)
#CKD controls
T2DCKD_controls <- table(UKBB_AG2$T2D_status,UKBB_AG2$CKD_only.0.0)[2,1]
#print(T2DCKD_controls)
#CKD extreme
#print('CKD ex')
T2DCKDex_cases <- table(UKBB_AG2$T2D_status,UKBB_AG2$CKD_ex.0.0)[2,2]
#print(T2DCKDex_cases)
#CKD extreme controls
T2DCKDex_controls <- table(UKBB_AG2$T2D_status,UKBB_AG2$CKD_ex.0.0)[2,1]
#print(T2DCKDex_controls)
#macro cases
T2Dmacro_cases <- table(UKBB_AG2$T2D_status,UKBB_AG2$macro.0.0)[2,2]
#print('macro')
#print(T2Dmacro_cases)
#macro_controls
T2Dmacro_controls <- table(UKBB_AG2$T2D_status,UKBB_AG2$macro.0.0)[2,1]
#print(T2Dmacro_controls)
#DN cases
T2D_DN_cases <- table(UKBB_AG2$T2D_status,UKBB_AG2$DN.0.0)[2,2]
#print('DN')
#print(T2D_DN_cases)
#DN controls
T2D_DN_controls <- table(UKBB_AG2$T2D_status,UKBB_AG2$DN.0.0)[2,1]
#print(T2D_DN_controls)
#All cases 
T2D_All_cases <- table(UKBB_AG2$T2D_status,UKBB_AG2$ALL.0.0)[2,2]
#print('All')
#print(T2D_All_cases)
#All controls
T2D_All_controls <- table(UKBB_AG2$T2D_status,UKBB_AG2$ALL.0.0)[2,1]
#print(T2D_All_controls)
#print("confirmed counts w. slides")
#DN-CKD cases
T2D_DNCKD_cases <- table(UKBB_AG2$T2D_status,UKBB_AG2$DNCKD2.0.0)[2,2]
#print('DN-CKD')
#print(T2D_DNCKD_cases)
#DN-CKD controls
T2D_DNCKD_controls <- table(UKBB_AG2$T2D_status,UKBB_AG2$DNCKD2.0.0)[2,1]
#print(T2D_DNCKD_controls)

#development section
# jj <- 'CKD_only.0.0'
# kidney_size <- T2DCKD_cases
#                    #identify sample ids for T2D cases
#                    T2D_cases <- UKBB_AG2 %>% filter(T2D_status==1 & (CKD_only.0.0=='CKD'|CKD_only.0.0=='CKD controls')) %>% 
#                                              select(f.eid) #20,120
# kidney_T2Dcases <- sample(T2D_cases$f.eid, size=kidney_size,replace = FALSE)
#                 #T2D kidney controls
#                 kidney_T2Dcontrols <- T2D_cases %>% filter(!f.eid %in% kidney_T2Dcases) 
# #dataset
# kidney_T2Dcases <- as.data.frame(kidney_T2Dcases)
# names(kidney_T2Dcases) <- 'f.eid'
# kidney_T2Dcases$shuffle_pheno <- 'cases'
# class(kidney_T2Dcases$shuffle_pheno)
# kidney_T2Dcontrols <- as.data.frame(kidney_T2Dcontrols)
# names(kidney_T2Dcontrols) <- 'f.eid'
# kidney_T2Dcontrols$shuffle_pheno <- 'controls'
# UKBB_shuffle <- natural_join(kidney_T2Dcases,kidney_T2Dcontrols,by="f.eid",jointype = "FULL") #full join
# head(UKBB_shuffle)
# table(UKBB_shuffle$shuffle_pheno)

# #shuffle
# UKBB_AG2_subset <- UKBB_AG2 %>% select(f.eid,T2D_status,AGE.0.0,SEX.0.0,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,all_of(jj)) 
# str(UKBB_AG2_subset)
# UKBB_shuffle <- merge(UKBB_shuffle,UKBB_AG2,by="f.eid",all=TRUE) %>% filter(shuffle_pheno=='cases'|
#                                                                                             shuffle_pheno=='controls'|
#                                                                                             !is.na(all_of(jj)))
# dim(UKBB_shuffle)
# UKBB_shuffle$jj <- 'CKD_only.0.0'
# UKBB_shuffle <- UKBB_shuffle %>% mutate(check22=case_when(T2D_status==1 & shuffle_pheno=='cases' & jj=='CKD_only.0.0' ~ 'CKD',
#                                           T2D_status==1 & shuffle_pheno=='controls' & jj=='CKD_only.0.0' ~ 'CKD controls',
#                                           T2D_status==0 & jj=='CKD_only.0.0' ~ as.character(CKD_only.0.0),
#                                           T2D_status==1 & shuffle_pheno=='cases' & jj=='micro.0.0' ~ 'micro',
#                                           T2D_status==1 & shuffle_pheno=='controls' & jj=='micro.0.0' ~ 'normo',
#                                           T2D_status==0 & jj=='micro.0.0' ~ as.character(micro.0.0)))


# #qa of this code
# table(UKBB_shuffle$T2D_status,UKBB_shuffle$check22)
# table(UKBB_shuffle$check22,UKBB_AG2$CKD_only.0.0,UKBB_shuffle$T2D_status)

##MODEL APPROACH #1 - Three Risk Groups

#Previously evaluated non-linearity in prior section#
#Current section evalautes the continuous form:

#counter
mm2 <- 1
#empty dataframes
coef_all4 <- data.frame()
coef_all4_t2d <- data.frame()
coef_all4_nd <- data.frame()

for (ii in c('relevel(as.factor(GRS_RAW_T2DIRf3),"low risk")','relevel(as.factor(GRS_WT_IR53f3),"low risk")')) {
    for(kk in c('Model1')){
          for (jj in c('CKD_only.0.0','CKD_ex.0.0','micro.0.0','macro.0.0','DN.0.0','ALL.0.0','DNCKD2.0.0')) {
              for (mm in c(iteration:iteration)) {
                  
                #generate new dataset
                RESULTS_CONT <- data.frame()
                RESULTS_OR <- data.frame()
            
                #set sampling size and generate the array with complete cases
                if (jj=='CKD_only.0.0') {
                   kidney_size <- T2DCKD_cases
                   #identify sample ids for T2D cases
                   T2D_cases <- UKBB_AG2 %>% filter(T2D_status==1 & (CKD_only.0.0=='CKD'|CKD_only.0.0=='CKD controls')) %>% 
                                             select(f.eid) #20,120
                } else if (jj=='micro.0.0') {
                    kidney_size <- T2Dmicro_cases
                    T2D_cases <- UKBB_AG2 %>% filter(T2D_status==1 & (micro.0.0=='normo'|micro.0.0=='micro')) %>% 
                                             select(f.eid) #20,291
                    
                } else if (jj=='macro.0.0') {
                    kidney_size <- T2Dmacro_cases
                    T2D_cases <- UKBB_AG2 %>% filter(T2D_status==1 & (macro.0.0=='normo'|macro.0.0=='macro')) %>% 
                                             select(f.eid)
                    
                } else if (jj=='CKD_ex.0.0') {
                    kidney_size <- T2DCKDex_cases
                    T2D_cases <- UKBB_AG2 %>% filter(T2D_status==1 & (CKD_ex.0.0=='CKD controls'|CKD_ex.0.0=='CKD extreme')) %>% 
                                             select(f.eid)
                    
                } else if (jj=='DN.0.0') {
                    kidney_size <- T2D_DN_cases
                    T2D_cases <- UKBB_AG2 %>% filter(T2D_status==1 & (DN.0.0=='yes'|DN.0.0=='no')) %>% 
                                             select(f.eid)
                    
                } else if (jj=='ALL.0.0') {
                    kidney_size <- T2D_All_cases
                    T2D_cases <- UKBB_AG2 %>% filter(T2D_status==1 & (ALL.0.0=='yes'|ALL.0.0=='no')) %>% 
                                             select(f.eid)
                    
                } else if (jj=='DNCKD2.0.0') {
                    kidney_size <- T2D_DNCKD_cases
                    T2D_cases <- UKBB_AG2 %>% filter(T2D_status==1 & (DNCKD2.0.0=='DNCKD'|DNCKD2.0.0=='DNCKD Control')) %>% 
                                             select(f.eid)
                    
                }
                  
                #Shift counter, needs to be based on the size
                mm2 <- as.numeric(iteration)*as.numeric(kidney_size) #shift the seed by the case count
                set.seed(mm2) #seed for shuffling 
                  
                #T2D kidney cases
                kidney_T2Dcases <- sample(T2D_cases$f.eid, size=kidney_size,replace = FALSE)
                #T2D kidney controls
                kidney_T2Dcontrols <- T2D_cases %>% filter(!f.eid %in% kidney_T2Dcases) 
                #combine into one dataset
                kidney_T2Dcases <- as.data.frame(kidney_T2Dcases)
                names(kidney_T2Dcases) <- 'f.eid'
                kidney_T2Dcases$shuffle_pheno <- 'cases'
                class(kidney_T2Dcases$shuffle_pheno)
                kidney_T2Dcontrols <- as.data.frame(kidney_T2Dcontrols)
                names(kidney_T2Dcontrols) <- 'f.eid'
                kidney_T2Dcontrols$shuffle_pheno <- 'controls'
                UKBB_shuffle <- natural_join(kidney_T2Dcases,kidney_T2Dcontrols,by="f.eid",jointype = "FULL") #full join, checked

                #select relevant variables from old dataset and merge
                UKBB_AG2_subset <- UKBB_AG2 %>% select(f.eid,AGE.0.0,SEX.0.0,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,all_of(jj)) 
                UKBB_shuffle <- merge(UKBB_shuffle,UKBB_AG2,by="f.eid",all=TRUE) %>% filter(shuffle_pheno=='cases'|
                                                                                            shuffle_pheno=='controls'|
                                                                                            !is.na(all_of(jj)))
                #recode cases based on shuffling for each phenotype
                UKBB_shuffle <- UKBB_shuffle %>% mutate(jj2:=case_when(T2D_status==1 & shuffle_pheno=='cases' & jj=='CKD_only.0.0' ~ 'CKD',
                                          T2D_status==1 & shuffle_pheno=='controls' & jj=='CKD_only.0.0' ~ 'CKD controls',
                                          T2D_status==0 & jj=='CKD_only.0.0' ~ as.character(CKD_only.0.0),
                                                                      
                                          T2D_status==1 & shuffle_pheno=='cases' & jj=='micro.0.0' ~ 'micro',
                                          T2D_status==1 & shuffle_pheno=='controls' & jj=='micro.0.0' ~ 'normo',
                                          T2D_status==0 & jj=='micro.0.0' ~ as.character(micro.0.0),
                                                                      
                                          T2D_status==1 & shuffle_pheno=='cases' & jj=='CKD_ex.0.0' ~ 'CKD extreme',
                                          T2D_status==1 & shuffle_pheno=='controls' & jj=='CKD_ex.0.0' ~ 'CKD controls',
                                          T2D_status==0 & jj=='CKD_ex.0.0' ~ as.character(CKD_ex.0.0),
                                                                      
                                          T2D_status==1 & shuffle_pheno=='cases' & jj=='macro.0.0' ~ 'macro',
                                          T2D_status==1 & shuffle_pheno=='controls' & jj=='macro.0.0' ~ 'normo',
                                          T2D_status==0 & jj=='macro.0.0' ~ as.character(macro.0.0),
                                                                      
                                          T2D_status==1 & shuffle_pheno=='cases' & jj=='DN.0.0' ~ 'yes',
                                          T2D_status==1 & shuffle_pheno=='controls' & jj=='DN.0.0' ~ 'no',
                                          T2D_status==0 & jj=='DN.0.0' ~ as.character(DN.0.0),
                                                       
                                          T2D_status==1 & shuffle_pheno=='cases' & jj=='ALL.0.0' ~ 'yes',
                                          T2D_status==1 & shuffle_pheno=='controls' & jj=='ALL.0.0' ~ 'no',
                                          T2D_status==0 & jj=='ALL.0.0' ~ as.character(ALL.0.0),
                  
                                          T2D_status==1 & shuffle_pheno=='cases' & jj=='DNCKD2.0.0' ~ 'DNCKD',
                                          T2D_status==1 & shuffle_pheno=='controls' & jj=='DNCKD2.0.0' ~ 'DNCKD Control',
                                          T2D_status==0 & jj=='DNCKD2.0.0'~ as.character(DNCKD2.0.0)))
                
                names(UKBB_shuffle)[which(names(UKBB_shuffle)=="jj2")] <- paste0(jj,"_v2") 
                  
                #Re-level, set ref
                if (jj=='CKD_only.0.0') {
                   UKBB_shuffle$CKD_only.0.0_v2 <- factor(UKBB_shuffle$CKD_only.0.0_v2)
                   UKBB_shuffle$CKD_only.0.0_v2 <- relevel(as.factor(UKBB_shuffle$CKD_only.0.0_v2),"CKD controls")
                } else if (jj=='micro.0.0') {
                   UKBB_shuffle$micro.0.0_v2 <- factor(UKBB_shuffle$micro.0.0_v2)
                   UKBB_shuffle$micro.0.0_v2 <- relevel(as.factor(UKBB_shuffle$micro.0.0_v2),"normo")
                } else if (jj=='macro.0.0') {
                   UKBB_shuffle$macro.0.0_v2 <- factor(UKBB_shuffle$macro.0.0_v2)
                   UKBB_shuffle$macro.0.0_v2 <- relevel(as.factor(UKBB_shuffle$macro.0.0_v2),"normo")
                } else if (jj=='CKD_ex.0.0') {
                   UKBB_shuffle$CKD_ex.0.0_v2 <- factor(UKBB_shuffle$CKD_ex.0.0_v2)
                   UKBB_shuffle$CKD_ex.0.0_v2 <- relevel(as.factor(UKBB_shuffle$CKD_ex.0.0_v2),"CKD controls")
                } else if (jj=='DN.0.0') {
                   UKBB_shuffle$DN.0.0_v2 <- factor(UKBB_shuffle$DN.0.0_v2)
                   UKBB_shuffle$DN.0.0_v2 <- relevel(as.factor(UKBB_shuffle$DN.0.0_v2),"no")
                } else if (jj=='ALL.0.0') {
                   UKBB_shuffle$ALL.0.0_v2 <- factor(UKBB_shuffle$ALL.0.0_v2)
                   UKBB_shuffle$ALL.0.0_v2 <- relevel(as.factor(UKBB_shuffle$ALL.0.0_v2),"no")
                } else if (jj=='DNCKD2.0.0') {
                   UKBB_shuffle$DNCKD2.0.0_v2 <- factor(UKBB_shuffle$DNCKD2.0.0_v2)
                   UKBB_shuffle$DNCKD2.0.0_v2 <- relevel(as.factor(UKBB_shuffle$DNCKD2.0.0_v2),"DNCKD Control")
                }   
                  
                #Used across formulas
                term <- ii   
                #Model 1 - Age, gender, PCI
                if(kk == "Model1"){

                fmla <- as.formula(paste0(paste0(jj,"_v2")," ~ ",term, "+ AGE.0.0 + SEX.0.0 + PC1 +
                                                PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"), env = environment()) }

                M1 <- glm(fmla, data=UKBB_shuffle, binomial(link="logit"))                                                                                
                M1_2 <-  M1 %>% summary()
                M1_3 <- anova(M1,test="LRT")

                #Model 1 - Sub-groups
                #Type 2 Diabetics
                M1_T2D <- UKBB_shuffle %>% filter(T2D_status == 1) %>% glm(formula=fmla,family=binomial(link="logit"))                                                                                
                M1_T2D2 <- M1_T2D  %>% summary()
                M1_T2D3 <- anova(M1_T2D ,test="LRT")

                #Non-Diabetes
                M1_ND <- UKBB_shuffle %>% filter(T2D_status == 0) %>% glm(formula=fmla, family=binomial(link="logit"))                                                                                
                M1_ND2 <- M1_ND  %>% summary()
                M1_ND3 <- anova(M1_ND ,test="LRT")

                TEMP<- list(shuffle_type=paste0("T2Dcases_it",mm2), var=jj,var2=ii, total=M1_2,t2d=M1_T2D2,nd=M1_ND2,lrt_tot=M1_3,
                            lrt_t2d=M1_T2D3,lrt_nd=M1_ND3)

                #Confidence Intervals
                #Entire Sample
                M1$dwn_conf <- exp(coefficients(TEMP$total)[,1]-1.96*coefficients(TEMP$total)[,2])
                M1$up_conf <- exp(coefficients(TEMP$total)[,1]+1.96*coefficients(TEMP$total)[,2])
                TABLE2 <- round(cbind(OR=exp(coefficients(TEMP$total)[,1]),CI=cbind(M1$dwn_conf,M1$up_conf)),3)   
                #Diabetics
                M1_T2D$dwn_conf <- exp(coefficients(TEMP$t2d)[,1]-1.96*coefficients(TEMP$t2d)[,2])
                M1_T2D$up_conf <- exp(coefficients(TEMP$t2d)[,1]+1.96*coefficients(TEMP$t2d)[,2])
                TABLE2_T2D <- round(cbind(OR=exp(coefficients(TEMP$t2d)[,1]),CI=cbind(M1_T2D$dwn_conf,M1_T2D$up_conf)),3)
                #Non-Diabetics
                M1_ND2$dwn_conf <- exp(coefficients(TEMP$nd)[,1]-1.96*coefficients(TEMP$nd)[,2])
                M1_ND2$up_conf <- exp(coefficients(TEMP$nd)[,1]+1.96*coefficients(TEMP$nd)[,2])
                TABLE2_ND <- round(cbind(OR=exp(coefficients(TEMP$nd)[,1]),CI=cbind(M1_ND2$dwn_conf,M1_ND2$up_conf)),3)

                #Summary Statistics, added shuffle # and shuffle_type, 6->7 entries long
                OR_CI <- list(shuffle_type="T2Dcases", var=jj,var2=ii,nd_ci=TABLE2_ND,t2d_ci=TABLE2_T2D,all_ci=TABLE2,shuffle_n=mm2)
                RESULTS_OR <- c(RESULTS_OR,OR_CI)    
                #Combined results
                RESULTS_CONT <- c(RESULTS_CONT,TEMP)
                  
                #Summary results - for export
                #Total Sample Score
                  coef_all3 <- data.frame()

                  #identifiers
                  iteration_id <- 1 #model id
                  model_id <- 2 #model id
                  i1 <- 4 #total models
                  i2 <- 5 #t2d models
                  i3 <- 6 #nd models
                  #lrt
                  i7 <- 7 #total models
                  i8 <- 8 #t2d models
                  i9 <- 9 #nd models 
                  term <- 1 - (1%/%10)*10 + 1 #remainder from subtracting from multiple of 30, depend on the model?
                  i4 <- 'Model 1'

                  ##################################ENTIRE SAMPLE 
                  #pull coefficients and convert to OR
                  coef_all <- data.frame(round(exp(coefficients(RESULTS_CONT[[i1]])),2)[c(2:3),1])
                  id_model <- names(round(exp(coefficients(RESULTS_CONT[[i1]])),2)[c(2:3),1])
                  model_id2 <- RESULTS_CONT[[model_id]]
                  rep_model <- rep(model_id2,dim(coef_all)[1])
                  rep_adj <- rep(i4,dim(coef_all)[1])

                  #pull coefficients and calculate 95% CI
                  up_coef <- data.frame(round(exp(coefficients(RESULTS_CONT[[i1]])[,1]+1.96*coefficients(RESULTS_CONT[[i1]])[,2]),2)[2:3])
                  down_coef <- data.frame(round(exp(coefficients(RESULTS_CONT[[i1]])[,1]-1.96*coefficients(RESULTS_CONT[[i1]])[,2]),2)[2:3])
                  #likelihood ratio
                  lrt_total <- RESULTS_CONT[[i7]][2,5]
                  lrt_total2 <- rep(as.character(lrt_total),dim(coef_all)[1])

                  coef_all2 <- cbind(rep_model,id_model) #model outcome to model var
                  coef_all2 <- cbind(coef_all2,coef_all) #coefficients
                  coef_all2 <- cbind(coef_all2,rep_adj)
                  coef_all2 <- cbind(coef_all2,down_coef) #95 CI
                  coef_all2 <- cbind(coef_all2,up_coef) #95 CI
                  coef_all2 <- cbind(coef_all2,lrt_total2)
                  coef_all3 <- rbind(coef_all3,coef_all2)

                    #reformat table
                    #rename
                    names(coef_all3) <- c("rep_model","id_model","OR","model_adj","lower","upper","lrt")
                    #combine HR and 95% CI
                    coef_all3$combo <- paste0(coef_all3$OR," (",coef_all3$lower,"-",coef_all3$upper,")")
                    coef_all3 <- subset(coef_all3, select = -c(3,5,6))
                    #high/low
                    coef_all3$category <- ifelse(str_detect(coef_all3$id_model,"high")==TRUE,"2_high","1_medium")
                    #substring
                    coef_all3$sub <- substr(coef_all3$id_model,19,32)
                    coef_all3$rep_model2 <- substr(coef_all3$rep_model,19,35)
                    #key
                    coef_all3$key <- paste0(coef_all3$sub,"-",coef_all3$rep_model2, model_id2,RESULTS_CONT[[iteration_id]])
        
                    #spread
                    coef_all3 <- subset(coef_all3, select = -c(1,2))
                    coeff_all3 <- spread(coef_all3,key=category,value=key)
                    coef_all4 <- rbind(coef_all4,coef_all3) 

                    ##############################T2D Sample Score
                    #Total Sample Score
                    coef_all3_t2d <- data.frame()
                    term <- 1 - (1%/%10)*10 + 1 #remainder from subtracting from multiple of 30, depend on the model?
                    i4 <- 'Model 1'

                      #ENTIRE SAMPLE 
                      #pull coefficients and convert to OR
                      coef_all <- data.frame(round(exp(coefficients(RESULTS_CONT[[i2]])),2)[c(2:3),1])
                      id_model <- names(round(exp(coefficients(RESULTS_CONT[[i2]])),2)[c(2:3),1])
                      model_id2 <- RESULTS_CONT[[model_id]]
                      rep_model <- rep(model_id2,dim(coef_all)[1])
                      rep_adj <- rep(i4,dim(coef_all)[1])

                      #pull coefficients and calculate 95% CI
                      up_coef <- data.frame(round(exp(coefficients(RESULTS_CONT[[i2]])[,1]+1.96*coefficients(RESULTS_CONT[[i2]])[,2]),2)[2:3])
                      down_coef <- data.frame(round(exp(coefficients(RESULTS_CONT[[i2]])[,1]-1.96*coefficients(RESULTS_CONT[[i2]])[,2]),2)[2:3])
                      #likelihood ratio
                      lrt_total <- RESULTS_CONT[[i8]][2,5]
                      lrt_total2 <- rep(as.character(lrt_total),dim(coef_all)[1])

                      coef_all2 <- cbind(rep_model,id_model) #model outcome to model var
                      coef_all2 <- cbind(coef_all2,coef_all) #coefficients
                      coef_all2 <- cbind(coef_all2,rep_adj)
                      coef_all2 <- cbind(coef_all2,down_coef) #95 CI
                      coef_all2 <- cbind(coef_all2,up_coef) #95 CI
                      coef_all2 <- cbind(coef_all2,lrt_total2)
                      coef_all3_t2d <- rbind(coef_all3_t2d,coef_all2)

                    #reformat table
                    #rename
                    names(coef_all3_t2d) <- c("rep_model","id_model","OR","model_adj","lower","upper","lrt")
                    #combine HR and 95% CI
                    coef_all3_t2d$combo <- paste0(coef_all3_t2d$OR," (",coef_all3_t2d$lower,"-",coef_all3_t2d$upper,")")
                    coef_all3_t2d <- subset(coef_all3_t2d, select = -c(3,5,6))
                    #high/low
                    coef_all3_t2d$category <- ifelse(str_detect(coef_all3_t2d$id_model,"high")==TRUE,"2_high","1_medium")
                    #substring
                    coef_all3_t2d$sub <- substr(coef_all3_t2d$id_model,19,32)
                    coef_all3_t2d$rep_model2 <- substr(coef_all3_t2d$rep_model,19,35)
                    #key
                    coef_all3_t2d$key <- paste0(coef_all3$sub,"-",coef_all3_t2d$rep_model2,model_id2,RESULTS_CONT[[iteration_id]])
                    #spread
                    coef_all3_t2d <- subset(coef_all3_t2d, select = -c(1,2))
                    coef_all3_t2d <- spread(coef_all3_t2d,key=category,value=key)
                    coef_all4_t2d <- rbind(coef_all4_t2d,coef_all3_t2d)

                    ####################################Non-Diabetic Sample Score
                    #Total Sample Score
                    coef_all3_nd <- data.frame()

                      #ENTIRE SAMPLE 
                      #pull coefficients and convert to OR
                      coef_all <- data.frame(round(exp(coefficients(RESULTS_CONT[[i3]])),2)[c(2:3),1])
                      id_model <- names(round(exp(coefficients(RESULTS_CONT[[i3]])),2)[c(2:3),1])
                      model_id2 <- RESULTS_CONT[[model_id]]
                      rep_model <- rep(model_id2,dim(coef_all)[1])
                      rep_adj <- rep(i4,dim(coef_all)[1])

                      #pull coefficients and calculate 95% CI
                      up_coef <- data.frame(round(exp(coefficients(RESULTS_CONT[[i3]])[,1]+1.96*coefficients(RESULTS_CONT[[i3]])[,2]),2)[2:3])
                      down_coef <- data.frame(round(exp(coefficients(RESULTS_CONT[[i3]])[,1]-1.96*coefficients(RESULTS_CONT[[i3]])[,2]),2)[2:3])
                      #likelihood ratio
                      lrt_total <- RESULTS_CONT[[i9]][2,5]
                      lrt_total2 <- rep(as.character(lrt_total),dim(coef_all)[1])

                      coef_all2 <- cbind(rep_model,id_model) #model outcome to model var
                      coef_all2 <- cbind(coef_all2,coef_all) #coefficients
                      coef_all2 <- cbind(coef_all2,rep_adj)
                      coef_all2 <- cbind(coef_all2,down_coef) #95 CI
                      coef_all2 <- cbind(coef_all2,up_coef) #95 CI
                      coef_all2 <- cbind(coef_all2,lrt_total2)
                      coef_all3_nd <- rbind(coef_all3_nd,coef_all2)

                    #reformat table
                    #rename
                    names(coef_all3_nd) <- c("rep_model","id_model","OR","model_adj","lower","upper","lrt")
                    #combine HR and 95% CI
                    coef_all3_nd$combo <- paste0(coef_all3_nd$OR," (",coef_all3_nd$lower,"-",coef_all3_nd$upper,")")
                    coef_all3_nd <- subset(coef_all3_nd, select = -c(3,5,6))
                    #high/low
                    coef_all3_nd$category <- ifelse(str_detect(coef_all3_nd$id_model,"high")==TRUE,"2_high","1_medium")
                    #substring
                    coef_all3_nd$sub <- substr(coef_all3_nd$id_model,19,32)
                    coef_all3_nd$rep_model2 <- substr(coef_all3_nd$rep_model,19,35)
                    #key
                    coef_all3_nd$key <- paste0(coef_all3_nd$sub,"-",coef_all3_nd$rep_model2,model_id2,RESULTS_CONT[[iteration_id]])
                    #spread
                    coef_all3_nd <- subset(coef_all3_nd, select = -c(1,2))
                    coef_all3_nd <- spread(coef_all3_nd,key=category,value=key)
                    coef_all4_nd <- rbind(coef_all4_nd,coef_all3_nd)
                  
                    #delete dataset to save space
                    rm(RESULTS_CONT,RESULTS_OR)
                  
                  }}}}


                    #Export the Sheet
                    write.csv(coef_all4,paste(temp_dir,"/All_10May21_T2DShuffle_it",iteration,".txt",sep=""))
                    write.csv(coef_all4_t2d,paste(temp_dir,"/T2D_10May21_T2DShuffle_it",iteration,".txt",sep=""))
                    write.csv(coef_all4_nd,paste(temp_dir,"/ND_10May21_T2DShuffle_it",iteration,".txt",sep=""))

