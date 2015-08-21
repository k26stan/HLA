## Compile & Reformat HLA Reference Panels for SNP2HLA ##
## To be used for JnJ & HLI Projects ##
## February 23, 2015 ##
## Kristopher Standish ##

#############################################################
## GAME PLAN ################################################
#############################################################

## Data Sets
 # SNP2HLA - T1D: /projects/janssen/Tools/SNP2HLA_package_v1.0/
   # A, B, C, DRB1, DQA, DQB, DPA, DPB
 # SNP2HLA - CEU: /projects/janssen/Tools/SNP2HLA_package_v1.0/MakeReference/HAPMAP_CEU_HLA.ped
   # A, B, C, DRB1, DQA, DQB, DPA, DPB
 # 1KG - /home/kstandis/HLI/HLA_Typing/Data/1KG/*txt
   # A, B, C, DRB1, DQB1
 # Chinese Set (Mary Carrington) - /home/kstandis/HLI/HLA_Typing/Data/MARY/NPC_GWAS_samples_HLA.xlsx
   # A, B, C

## SNP2HLA: MakeReference
 # Input Files
   # bed/bim/fam files w/ SNPs
   # Plink formatted HLA type file (PED)
     # e.g., FID IID FATHER MOTHER SEX PHENOTYPE A B C DRB1 DQA DQB DPA DPB
     # Fill with 0 if no type is known/provided

## Overall Game Plan
 # Create single, complete reference panel w/ all available data sets
 # 10-fold cross validation, predicting HLA types
 # Check accuracy amongst cohorts and genes
 # Compare to typing w/ less diverse HLA reference cohort

## Immediate Game Plan
 # Reformat data into what would be appropriate to make a reference panel w/ SNP2HLA
 # Use all available reference data sets
 # Figure out how to deal with ambiguities

library(xlsx)

#############################################################
## LOAD DATA ################################################
#############################################################

## Set Paths to Data Sets & Save Locations (TSCC)
PathToT1D <- "/projects/janssen/Tools/SNP2HLA_package_v1.0/T1DGC/T1DGC_REF.HLA.bgl.phased"
PathToCEU <- "/projects/janssen/Tools/SNP2HLA_package_v1.0/MakeReference/HAPMAP_CEU_HLA.ped"
PathTo1KG <- "/home/kstandis/HLI/HLA_Typing/Data/1KG/20140702_hla_diversity.txt"
PathToMary <- "/home/kstandis/HLI/HLA_Typing/Data/MARY/NPC_GWAS_samples_HLA.txt"

## Set Paths to Data Sets & Save Locations (Mac)
# cd /home/kstandis/Data/HLI_Phase/20150223_HLA_Ref/
# head -4 T1DGC_REF.bgl.phased > T1DGC_REF.HLA.bgl.phased 
# cat T1DGC_REF.bgl.phased | grep HLA >> T1DGC_REF.HLA.bgl.phased
PathToT1D <- "/home/kstandis/Data/HLI_Phase/20150223_HLA_Ref/T1DGC_REF.HLA.bgl.phased"
PathToCEU <- "/home/kstandis/Data/HLI_Phase/20150223_HLA_Ref/HAPMAP_CEU_HLA.ped"
PathTo1KG <- "/home/kstandis/Data/HLI_Phase/20150223_HLA_Ref/20140702_hla_diversity.txt"
PathToMary <- "/home/kstandis/Data/HLI_Phase/20150223_HLA_Ref/NPC_GWAS_samples_HLA.txt"

## Load HLA Types
T1D.l <- read.table( PathToT1D, header=T, sep="" )
CEU.l <- read.table( PathToCEU, header=F, sep="" )
KG1.l <- read.table( PathTo1KG, header=T, sep="" )
MAR.l <- read.table( PathToMary, header=T, sep="\t" )
























#############################################################
## END OF DOC ###############################################
#############################################################
