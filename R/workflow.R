library(tidyverse)
library(MSstats)
library(PeCorA)
library(RCurl)
library(here)

####getting the DDA.ipgr data#####
#https://github.com/tsunghengtsai/MSstats-feature-selection/blob/master/data/input.dda.iprg.pg.rda
load(here("data","input.dda.iprg.pg.rda"))

#PeCoRa analysis
ddaipgrfile<-input.dda.iprg.pg

import_preprocessed_for_PeCoRA <- function(x){
  dda.iprg <- dataProcess(x)
  dda_formated<-dda.iprg$ProcessedData  %>%
    select(PEPTIDE,PROTEIN, FEATURE, GROUP_ORIGINAL, RUN,ABUNDANCE) %>%
    dplyr::rename(Peptide=PEPTIDE,
                  Protein =PROTEIN,
                  Peptide.Modified.Sequence=FEATURE,
                  Condition=GROUP_ORIGINAL,
                  BioReplicate=RUN,
                  Normalized.Area= ABUNDANCE)
  return(dda_formated)
}
t<-import_preprocessed_for_PeCoRA(ddaipgrfile)
summary(t$Normalized.Area)

scaled_peptides <- PeCorA_preprocessing(t,
                                        area_column_name=6,
                                        threshold_to_filter=10,
                                        control_name="Condition1")

disagree_peptides <- PeCorA (scaled_peptides)

dim(disagree_peptides) #with threshold_to_filter=5
#9371 4 

dim(disagree_peptides) #with threshold_to_filter=1
#9372 4 

dim(disagree_peptides) #with threshold_to_filter=10
#9350 4 

write.table(disagree_peptides,here("files","DDAiPRG_Pecora.txt"), sep="\t", quote=F, row.names = F)

####Removing disagree peptides####

pecora_DDAiPRG_tofilter<-disagree_peptides %>%
  dplyr::filter(adj_pval<0.01) %>%
  dplyr::rename(ProteinName=protein) %>%
  separate(ProteinName,into=c("ProteinName","specie"),sep="_") %>%
  mutate(tofilter=paste(ProteinName,peptide,sep="_"))


input.dda_tofilter <- input.dda.iprg.pg %>%
  separate(ProteinName,into=c("ProteinName","specie"),sep="_") %>%
  unite(IsotopeLabelType,ProteinName,PeptideModifiedSequence,PrecursorCharge,FragmentIon,ProductCharge,sep="_") %>%
  mutate(tofilter=paste(IsotopeLabelType,"all",sep = "_"))
dim(input.dda_tofilter) 
#115752      7

new_input.dda.iprg.pg<-
  anti_join(input.dda_tofilter,pecora_DDAiPRG_tofilter,by="tofilter")

new_input.dda.iprg.pg<-new_input.dda.iprg.pg %>%
  separate(IsotopeLabelType,into=c("ProteinName","PeptideModifiedSequence", "PrecursorCharge","FragmentIon","ProductCharge"),sep="_") %>%
  mutate(IsotopeLabelType="L") %>%
  select(-tofilter) %>%
  select(-specie)

dim(new_input.dda.iprg.pg)  
#114348     10


####MSstats####

##All
DDA.proposed <- dataProcess(raw = input.dda.iprg.pg,
                            normalization = 'equalizeMedians',
                            summaryMethod = 'TMP',
                            censoredInt = "NA",
                            cutoffCensored = "minFeature",
                            MBimpute = TRUE,
                            maxQuantileforCensored=0.999)

#save(DDA.proposed, here("data",file="DDA.proposed.rda"))

##Pecora then MSstats
new_DDA.proposed <- dataProcess(raw = new_input.dda.iprg.pg,
                                normalization = 'equalizeMedians',
                                summaryMethod = 'TMP',
                                censoredInt = "NA",
                                cutoffCensored = "minFeature",
                                MBimpute = TRUE,
                                maxQuantileforCensored=0.999)

#save(new_DDA.proposed, here("data",file="new_DDA.proposed_F1.rda")) #threshol 1

##Feature selection
DDA.inf <- dataProcess(raw = input.dda.iprg.pg,
                       normalization = 'equalizeMedians',
                       summaryMethod = 'TMP',
                       featureSubset = "highQuality",
                       remove_uninformative_feature_outlier = TRUE)
#uninformative features does not seem to get removed in this version of MSstats!!
DDA.inf$ProcessedData %>% 
 # dplyr::filter(!(feature_quality== "Uninformative")) %>% 
  #       dplyr::filter(!(is_outlier == "TRUE")) %>% 
  count(feature_quality,is_outlier)

DDA.inf$ProcessedData[which(DDA.inf$ProcessedData$PROTEIN=="sp|P32898|CYM1_YEAST"),]

#save(DDA.inf, here("data",file="DDA.inf.rda"))


spike_in<-c("sp|P44015|VAC2_YEAST", "sp|P55752|ISCB_YEAST",
            "sp|P44374|SFG2_YEAST", "sp|P44983|UTR6_YEAST",
            "sp|P44683|PGA4_YEAST", "sp|P55249|ZRT4_YEAST")

mycomparison1<-matrix(c(-1,1,0,0),nrow=1)
mycomparison2<-matrix(c(-1,0,1,0),nrow=1)
mycomparison3<-matrix(c(-1,0,0,1),nrow=1)
mycomparison4<-matrix(c(0,-1,1,0),nrow=1)
mycomparison5<-matrix(c(0,-1,0,1),nrow=1)
mycomparison6<-matrix(c(0,0,-1,1),nrow=1)

mycomparison<-rbind(mycomparison1,mycomparison2, mycomparison3,mycomparison4,mycomparison5,mycomparison6)

row.names(mycomparison)<-c("C2-C1","C3-C1","C4-C1","C2-C3","C2-C4","C3-C4")

DDA.comparisons <- groupComparison(contrast.matrix = mycomparison, data = DDA.proposed)
DDA.comparisons$ComparisonResult %>% dim()
#18090    11

####count FP####
DDA.comparisons$ComparisonResult %>% 
  filter(adj.pvalue<0.01) %>%
  filter(!Protein %in% spike_in) %>% dim()
#106  11

fp_names<-DDA.comparisons$ComparisonResult %>% 
  filter(adj.pvalue<0.05) %>%
  filter(!Protein %in% spike_in) %>% pull(Protein)

#save(new_DDA.comparisons, here("data",file = "new_DDA.comparisons.RData"))

new_DDA.comparisons <- groupComparison(contrast.matrix = mycomparison, data = new_DDA.proposed)
new_DDA.comparisons$ComparisonResult %>% dim()
new_DDA.comparisons$ComparisonResult %>% 
  filter(adj.pvalue<0.01) %>%
  filter(!Protein %in% spike_in) %>%
  filter_all(all_vars(!is.infinite(.))) %>% 
  mutate(Protein=paste(Protein,"_YEAST",sep="")) %>%
 # filter(!Protein %in% fp_names) %>%
  dim()
# 86 11


DDA.inf.featureselected.comparisons <- groupComparison(contrast.matrix = mycomparison, data = DDA.inf)
DDA.inf.featureselected.comparisons$ComparisonResult %>% 
  filter(adj.pvalue<0.01) %>%
  filter(!Protein %in% spike_in) %>%
  #filter(!Protein %in% fp_names) %>%
  dim()

