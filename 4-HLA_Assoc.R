## Look for Association b/n HLA Amino Acids & Drug Response ##
## May 12, 2015 ##
## Updated Sept 25, 2015 ##
## Kristopher Standish ##

## GAME PLAN ##
 # Sept 25, 2015

## Brainstorming
 # Phenotypes
   # Drug Response (delta-DAS,lCRP,rSJC,rTJC)
   # Disease Severity
   # RF/ACPA Status
 # HLA Predictors
   # ANOVA over all Types for each Gene
     # 2- & 4-Digit Type  
   # Additive/Dominant/Recessive for each Type for each Gene
     # 2- & 4-Digit Type
   # Amino-Acid Level Tests for each Gene (4- or Best-Digit Precision)
     # ANOVA over all Amino Acids at each Position
     # Additive/Dominant/Recessive for each Type for each Gene
   # Collapsed Amino-Acid Haplotypes
     # DRB1
       # Pos 11,71,74
       # Pos 70-74 (Shared Epitope)
       # Pos 11,13,71,74
     # B
       # Pos 9
     # DPB1
       # Pos 9

library(gplots)

#############################################################
## LOAD DATA ################################################
#############################################################

## Set Date
DATE <- gsub("-","",Sys.Date())

## Set Paths to HLA Data Sets
PathToTypes <- "/Users/kstandis/Data/Burn/Data/HLA/20150512_SOAP_HLA_Types/20150511_HLA_Types.Rdata"
PathToTypes <- "/Users/kstandis/Data/Burn/Data/HLA/20151211_SOAP_HLA_Types/20151211_HLA_Types.Rdata"
PathToAA <- "/Users/kstandis/Data/Burn/Data/HLA/20150512_SOAP_HLA_Types/20150512_HLA_AA.Rdata"
PathToAA <- "/Users/kstandis/Data/Burn/Data/HLA/20151211_SOAP_HLA_Types/20151211_HLA_AA.Rdata"
PathTo1KG <- "/Users/kstandis/Data/HLI_Phase/20150223_HLA_Ref/20140702_hla_diversity.txt"
PathTo1KG.2 <- "/Users/kstandis/Data/TBL/Data/20141001/1KG/Panel_Key.txt"
PathToRefs <- "/Users/kstandis/Data/HLI_Phase/20150223_HLA_Ref/Alignments_Rel_3190/"
 # http://www.ebi.ac.uk/ipd/imgt/hla/nomenclature/alignments.html
PathToFT <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20141229_Full_Table.txt"
PathToOut <- "/Users/kstandis/Data/Burn/Data/HLA/20150512_SOAP_HLA_Types/"
PathToPlot <- paste("/Users/kstandis/Data/Burn/Plots/",DATE,"_HLA_Assoc/",sep="")
dir.create( PathToPlot )

## Load Janssen HLA Results
 # Types
load( PathToTypes )
TYPES.l <- COMPILE
 # Amino Acids
load( PathToAA )
AA.l <- COMPILE

## 1000 Genoems HLA Data
HLA.l <- read.table(PathTo1KG,header=T)
HLA.key <- read.table(PathTo1KG.2,header=T)
HLA <- merge( HLA.key, HLA.l, by.y="id",by.x="sample")

## Get Phenotype Info
FT.l <- read.table( PathToFT, sep="\t",header=T )

#############################################################
## ORGANIZE DATA ############################################
#############################################################

## Filter Phenotypic Data (based on length of participation in study)
RM.LT8 <- which( FT.l$IN < 8 )
RM.LT8.samps <- as.character( FT.l$ID[RM.LT8] )
FT <- FT.l[ -RM.LT8, ]

## Pull out Patient Data
PAT_DOS <- TYPES.l$GENES.2.list
PAT_TYP <- TYPES.l$GENES
PAT_AA <- AA.l$PAT_AA
 # Print Example
X <- 1:10
PAT_AA$DRB1[X,X]

## Patients & Genes
PATS <- colnames(PAT_DOS$A)
N.PATS <- length(PATS)
GENE_LIST <- names(AA.l$PAT_AA)
N.GENE <- length(GENE_LIST) # nrow(PAT_TYP)
# GENE_LIST <- rownames(TYPES.l$GENES)

## Cut Patient Types to 2 & 4 Digit Precision
 # PAT_TYP
PAT_TYP.4 <- PAT_TYP.2 <- PAT_TYP
PAT_TYP.4 <- apply( PAT_TYP, 1, function(x) apply( sapply( strsplit(x,":"), "[",1:2),2, function(y) paste(y,collapse=":") ) )
PAT_TYP.4 <- apply( PAT_TYP.4, 1, function(x) gsub("[A-z]","",x) )
PAT_TYP.2 <- apply( PAT_TYP, 1, function(x) sapply( strsplit(x,":"), "[",1) )
PAT_TYP.2 <- apply( PAT_TYP.2, 1, function(x) gsub("[A-z]","",x) )
 # PAT_DOS
PAT_DOS.4 <- PAT_DOS.2 <- list()
for ( g in 1:N.GENE ) { gene <- GENE_LIST[g]
	unique_haps <- unique( PAT_TYP.4[gene,] )
	PAT_DOS.4[[gene]] <- array( 0,c(length(unique_haps),N.PATS) )
	colnames(PAT_DOS.4[[gene]]) <- PATS ; rownames(PAT_DOS.4[[gene]]) <- unique_haps
	unique_haps <- na.omit( unique( PAT_TYP.2[gene,] ) )
	PAT_DOS.2[[gene]] <- array( 0,c(length(unique_haps),N.PATS) )
	colnames(PAT_DOS.2[[gene]]) <- PATS ; rownames(PAT_DOS.2[[gene]]) <- unique_haps
	for ( p in 1:N.PATS ) { pat <- PATS[p]
		which_haps <- PAT_TYP.4[ gene,paste(pat,1:2,sep="_") ]
		if ( length(which(duplicated(which_haps)))>0 ) {
			PAT_DOS.4[[gene]][which_haps,pat] <- 2
		}else{
			PAT_DOS.4[[gene]][which_haps,pat] <- 1
		}
		which_haps <- na.omit( PAT_TYP.2[ gene,paste(pat,1:2,sep="_") ] )
		if ( length(which(duplicated(which_haps)))>0 ) {
			PAT_DOS.2[[gene]][which_haps,pat] <- 2
		}else{
			PAT_DOS.2[[gene]][which_haps,pat] <- 1
		}
	}
}


#############################################################
## PLOT ALLELE FREQUENCIES ##################################
#############################################################

## Plot Allele Frequency of Each Haplotype
COLS <- c("seagreen3","mediumpurple3")
COLS.list <- c("firebrick2","chocolate2","gold1","chartreuse2","cadetblue2","dodgerblue2","slateblue2","magenta2")
for ( g in 1:N.GENE ) { gene <- GENE_LIST[g]
	FREQ.4 <- rowSums( PAT_DOS.4[[gene]] ) ; names(FREQ.4)[which(names(FREQ.4)==":")]<-"NA" ; FREQ.4 <- FREQ.4[order(names(FREQ.4))]
	# FREQ.4 <- table( PAT_TYP.4[gene,] ) ; FREQ.4 <- FREQ.4[order(names(FREQ.4))]
	FREQ.2 <- rowSums( PAT_DOS.2[[gene]] ) ; FREQ.2 <- FREQ.2[order(names(FREQ.2))]
	COLS.2 <- colorRampPalette(COLS.list)(length(FREQ.2)) ; names(COLS.2) <- names(FREQ.2)
	COLS.4 <- COLS.2[ sapply(strsplit(names(FREQ.4),":"),"[",1) ]
	png( paste(PathToPlot,"1-HapFreq.",gene,".png",sep=""), height=1200,width=1600,pointsize=30 )
	par(mfrow=c(2,1)) ; par(mar=c(4,4,3,2))
	# 4 digit
	barplot( FREQ.4, las=2, col=COLS.4, main=paste("HLA -",gene,": 4-Digit Haplotype Frequency"),xlab="Type (4-digit)",ylab="# Haplotypes" )
	barplot( FREQ.2, las=2, col=COLS.2, main=paste("HLA -",gene,": 2-Digit Haplotype Frequency"),xlab="Type (2-digit)",ylab="# Haplotypes" )
	dev.off()
}

## Plot AA Level Diversity within Cohort
COLS.list <- c("firebrick2","chocolate2","gold1","chartreuse2","cadetblue2","dodgerblue2","slateblue2","magenta2")
COLS.AA <- c(colorRampPalette(COLS.list)(26),"black","grey90","grey50") ; names(COLS.AA) <- c(LETTERS,"*",".","?")
for ( g in 1:length(GENE_LIST) ) {
	gene <- GENE_LIST[g]
	# Pull out Patient Data for Gene
	pat_typ <- PAT_TYP[gene,]
	pat_dos <- PAT_DOS[[gene]]
	pat_aa.r <- pat_aa <- PAT_AA[[gene]] # Raw Table
	ISNA <- which(is.na(pat_aa.r))
	if ( length(ISNA)>0 ) { pat_aa[ISNA] <- "?" } # After converting missing values to "?"

	## Get Amino Acid Frequencies
	 # What is the maximum number of AA's at one position?
	max.temp <- max( unlist(lapply( apply( pat_aa, 2, table ), length)) )
	 # Which Positions have >1 Known AA? (Use Raw Table before "?" were included)
	Find_Div <- function( aa_vector ) {
		TAB <- table( aa_vector )
		To_Ignore <- which(names(TAB) %in% c("?","*",".") )
		TAB <- TAB[-To_Ignore]
	}
	which.temp <- which( unlist(lapply( apply( pat_aa, 2, Find_Div ), length)) > 1 )
	 # Get AA Frequencies/Proportions & Names
	pat_tab <- sapply(apply( pat_aa, 2, table ), "[", 1:max.temp ) ; colnames(pat_tab) <- gsub("Pos_","",colnames(pat_tab))
	pat_tab.names <- sapply(apply( pat_aa, 2, function(x) names(table(x)) ), "[", 1:max.temp )
	pat_tab.prc <- pat_tab / nrow(pat_aa)
	 # Plot
	XLIM <- range(as.numeric(colnames(pat_tab)))
	png( paste(PathToPlot,"/2-AAfreq_",gene,".png",sep=""), height=1000,width=2000,pointsize=30 )
	layout( matrix(c(1,2), 1, 2, byrow = TRUE), widths=c(7,1) ) # layout( matrix(c(1,2,3,4), 2, 2, byrow = TRUE), widths=c(3,2), heights=c(1,1) )
	plot( 0,0,type="n",xlim=XLIM,ylim=c(0,1), main=paste("Amino Acid Frequency",gene),xlab="Amino Acid",ylab="Frequency",xaxt="n")
	axis( 1, at=seq(-1000,1000,20),las=2 )	
	# abline( v=seq(-1000,1000,20),lty=3,col="grey50" )
	abline( h=seq(0,1,.2),lty=3,col="grey50" ) ; abline( h=c(0,1) )
	allele_bar <- function(x) {
		xval <- as.numeric(colnames(pat_tab)[x])
		yvals <- cumsum(pat_tab.prc[,x])
		barplot( t(t(pat_tab.prc[,x])),add=T,xaxt="n",beside=F,border=NA,col=COLS.AA[pat_tab.names[,x]],space=c(xval,0),width=1)
	}
	DWAI <- lapply( which.temp, allele_bar )
	## Amino Acid Key
	Which_AA <- which(names(COLS.AA)%in% pat_aa )
	COLS.aa <- COLS.AA[Which_AA]
	barplot( matrix(rep(1,length(COLS.aa)),ncol=1), beside=F,col=COLS.aa,xaxt="n",yaxt="n",ylab="AA")
	axis( 2, at=1:length(COLS.aa)-.5, label=names(COLS.aa),las=2,cex.axis=.8 )
	dev.off()
}
	
#############################################################
## RUN ASSOCIATION TESTS ####################################
#############################################################

## Specify Phenos/Covs
PHENOS <- c("DEL_MNe_MN","DEL_lCRP_MNe_MN","DEL_rSJC_MNe_MN","DEL_rTJC_MNe_MN",
	"DAS_BL_MN","lCRP_BL_MN","rSJC_BL_MN","rTJC_BL_MN",
	"ACPA","RF" )
COVS <- c("DAS_BL_MN","lCRP_BL_MN","rSJC_BL_MN","rTJC_BL_MN",
	"RF_ACPA+log(DIS_DUR)","RF_ACPA+log(DIS_DUR)","RF_ACPA+log(DIS_DUR)","RF_ACPA+log(DIS_DUR)",
	"","")
PH.COV <- data.frame(PHENOS,COVS)

## Loop Through Genes & Test for Assoc
P.precise <- B.precise <- list()
for ( PRECISE in c(2,4) ) {
	P <- B <- list()
	P$TYP <- P$TYP_DOS <- P$AA_AOV <- P$AA_DOS <- list()
	B$TYP <- B$TYP_DOS <- B$AA_AOV <- B$AA_DOS <- list()
	# for ( g in 1:N.GENE ) {
	for ( g in 1:length(GENE_LIST) ) {
		gene <- GENE_LIST[g]
		print(paste("### Running Association on",gene))
		## Pull out Patient Haplotype Assignments for Gene
		pat_aa <- PAT_AA[[gene]]
		if ( PRECISE==4 ) {
			pat_typ <- PAT_TYP.4[gene,]
			pat_dos.a <- PAT_DOS.4[[gene]]
		}else{
			if ( PRECISE==2 ) {
				pat_typ <- PAT_TYP.2[gene,]
				pat_dos.a <- PAT_DOS.2[[gene]]	
			}else{
				PRECISE <- 0
				pat_typ <- PAT_TYP[gene,]
				pat_dos.a <- PAT_DOS[[gene]]
				pat_aa <- PAT_AA[[gene]]	
			}
		}

		## Filter Haplotypes to Common Haplotypes
		print("Filtering & Merging")
		HAP.freqs <- rowSums(pat_dos.a)
		which_common <- which(HAP.freqs>.05*N.PATS)
		which_common <- which(HAP.freqs>10)
		which_common.haps <- names(which_common)
		if ( length(which_common)==1 ) { next
		}else{
			pat_dos <- pat_dos.a[ which_common, -which(colnames(pat_dos.a)%in%RM.LT8.samps) ]
			rownames(pat_dos) <- paste("T",rownames(pat_dos),sep="")
			rownames(pat_dos) <- gsub(":","",rownames(pat_dos),fixed=T)
		}
		
		## Merge Files
		 # Merge Diplotype Tables w/ Pheno
		MG.DIP.1 <- merge( pat_aa, data.frame(PAT=names(pat_typ),TYP=pat_typ), by.x="row.names",by.y="PAT" )
		colnames(MG.DIP.1)[1] <- "DIP"
		MG.DIP.1 <- data.frame( ID=sapply(strsplit(MG.DIP.1[,1],"_"),"[",1), MG.DIP.1 )
		MG.DIP <- merge( FT, MG.DIP.1, by="ID")
		MG.DIP[,"ACPA"] <- MG.DIP[,"ACPA"]=="Positive"
		MG.DIP.com <- MG.DIP[ which(MG.DIP$TYP %in% which_common.haps), ]
		 # Merge w/ Phenotype Data
		MG.DOS <- merge( FT, t(pat_dos), by.x="ID",by.y="row.names" )

		## Run Association Tests
		print("ANOVA (Haplotypes)")
		 # ANOVA of Common Haplotypes vs Phenotypes
		B.TYP <- list()
		P.TYP <- numeric( length(PHENOS) )
		names(P.TYP) <- PHENOS
		for ( p in 1:length(PHENOS) ) {
			pheno <- PHENOS[p]
			cov <- COVS[p]
			formula <- as.formula(paste( pheno,"~",cov,"+TYP" ))
			if ( pheno %in% c("RF_ACPA","ACPA","RF") ) {
				MOD <- chisq.test( table( MG.DIP.com[,pheno],as.character(MG.DIP.com[,"TYP"]) ) )
				P.TYP[p] <- MOD$p.value
			}else{
				MOD <- lm( formula, data=MG.DIP.com )
				P.TYP[p] <- anova(MOD)["TYP","Pr(>F)"]
			}
			B.TYP[[pheno]] <- MOD
		}
		P$TYP[[gene]] <- P.TYP
		B$TYP[[gene]] <- B.TYP

		 # Additive Model for Haplotypes vs Phenotypes
		print("Additive (Haplotypes)")
		B.TYP_DOS <- list()
		P.TYP_DOS <- array( ,c(length(PHENOS),nrow(pat_dos)))
		colnames(P.TYP_DOS) <- rownames(pat_dos)
		rownames(P.TYP_DOS) <- PHENOS
		for ( p in 1:length(PHENOS) ) {
			pheno <- PHENOS[p]
			cov <- COVS[p]
			B.TYP_DOS[[pheno]] <- list()
			for ( d in 1:ncol(P.TYP_DOS) ) {
				dip <- colnames(P.TYP_DOS)[d]
				UNIQ <- unique(MG.DOS[,dip]) ; UNIQ <- UNIQ[which(!is.na(UNIQ))]
				if ( length(UNIQ)==1 ) {
					P.TYP_DOS[p,d] <- NA
					B.TYP_DOS[[pheno]][[dip]] <- "No_Mod"
				}else{
					# print(d)
					formula <- as.formula(paste( pheno,"~",cov,"+",dip ))
					if ( pheno %in% c("RF_ACPA","ACPA","RF") ) {
						MOD <- glm( formula, data=MG.DOS, family=binomial(logit) )
						P.TYP_DOS[p,d] <- summary(MOD)$coefficients[dip,"Pr(>|z|)"]
					}else{
						MOD <- lm( formula, data=MG.DOS )
						P.TYP_DOS[p,d] <- anova(MOD)[dip,"Pr(>F)"]
					}
				}
				B.TYP_DOS[[pheno]][[dip]] <- MOD
			}
		}
		P$TYP_DOS[[gene]] <- P.TYP_DOS
		B$TYP_DOS[[gene]] <- B.TYP_DOS

		 # ANOVA for all Amino Acids at each Position
		print("ANOVA (Amino Acids)")
		B.AA_AOV <- list()
		P.AA_AOV <- array( ,c(length(PHENOS),ncol(MG.DIP.1)-2))
		colnames(P.AA_AOV) <- colnames(MG.DIP.1)[3:ncol(MG.DIP.1)]
		rownames(P.AA_AOV) <- PHENOS
		for ( p in 1:length(PHENOS) ) {
			pheno <- PHENOS[p]
			cov <- COVS[p]
			B.AA_AOV[[pheno]] <- list()
			for ( d in 1:ncol(P.AA_AOV) ) {
				dip <- colnames(P.AA_AOV)[d]
				UNIQ <- unique(MG.DIP[,dip])
				 # Remove Unknown Amino Acids
				To_Ignore <- which( UNIQ%in%c("*","?",".") | is.na(UNIQ) )
				if ( length(To_Ignore)>0 ) { UNIQ <- UNIQ[-To_Ignore] } # UNIQ <- UNIQ[which(!is.na(UNIQ))]
				 # Subset Merged Table to only consider Known & Common (MAF>10) Amino Acids
				UNIQ.freq <- table(MG.DIP[,dip])
				UNIQ.rare <- names(UNIQ.freq)[which(UNIQ.freq<10)]
				UNIQ.com <- setdiff( UNIQ, UNIQ.rare )
				SUBSET <- which( MG.DIP[,dip]%in%c("*","?",".") | is.na(MG.DIP[,dip]) | MG.DIP[,dip]%in%UNIQ.rare )
				 # If more than 1 Known & Unique & Common Value exists
				if ( length(UNIQ.com)<=1 ) {
					P.AA_AOV[p,d] <- NA
					B.AA_AOV[[pheno]][[dip]] <- "No_Mod"
				}else{
					# print(d)
					formula <- as.formula(paste( pheno,"~",cov,"+",dip ))
					# MG.DIP <- MG.DIP[-SUBSET,]
					if ( pheno %in% c("RF_ACPA","ACPA","RF") ) {
						# MOD <- chisq.test( table( MG.DIP[,pheno],as.character(MG.DIP[,dip]) ) )
						if ( length(SUBSET)>0 ) {
							MOD <- chisq.test( table( MG.DIP[-SUBSET,pheno],as.character(MG.DIP[-SUBSET,dip]) ) )
						}else{
							MOD <- chisq.test( table( MG.DIP[,pheno],as.character(MG.DIP[,dip]) ) )
						}
						P.AA_AOV[p,d] <- MOD$p.value
					}else{
						if ( length(SUBSET)>0 ) {
							MOD <- lm( formula, data=MG.DIP, subset=-SUBSET )
						}else{
							MOD <- lm( formula, data=MG.DIP )
						}
						P.AA_AOV[p,d] <- anova(MOD)[dip,"Pr(>F)"]
					}
				}
				B.AA_AOV[[pheno]][[dip]] <- MOD
			} # Close Dip Loop
			# print(paste("Done with",pheno))
		} # Close Pheno Loo
		P$AA_AOV[[gene]] <- P.AA_AOV
		B$AA_AOV[[gene]] <- B.AA_AOV

		 # Additive Model for Amino Acids vs Phenotypes
		print("Additive (Amino Acids)")
		P.AA_DOS <- list()
		B.AA_DOS <- list()
		for ( p in 1:length(PHENOS) ) {
			pheno <- PHENOS[p]
			cov <- COVS[p]
			covs <- strsplit(cov,"+",fixed=T)[[1]] ; covs <- gsub("log(DIS_DUR)","DIS_DUR",covs,fixed=T)
			P.AA_DOS[[pheno]] <- list()
			B.AA_DOS[[pheno]] <- list()
			for ( d in 1:ncol(P.AA_AOV) ) {
				dip <- colnames(P.AA_AOV)[d]
				UNIQ <- unique(MG.DIP[,dip])
				To_Ignore <- which( UNIQ%in%c("*","?",".") | is.na(UNIQ) )
				if ( length(To_Ignore)>0 ) { UNIQ <- UNIQ[-To_Ignore] } # UNIQ <- UNIQ[which(!is.na(UNIQ))]
				 # Subset Merged Table to only consider Known Amino Acids
				UNIQ.freq <- table(MG.DIP[,dip])
				UNIQ.rare <- names(UNIQ.freq)[which(UNIQ.freq<10)]
				UNIQ.com <- setdiff( UNIQ, UNIQ.rare )
				# SUBSET <- which( MG.DIP[,dip]%in%c("*","?",".") | is.na(MG.DIP[,dip]) )# | MG.DIP[,dip]%in%UNIQ.rare )
				if ( length(UNIQ.com)<=1 ) {
					P.AA_DOS[[pheno]][[dip]] <- "NA"
					B.AA_DOS[[pheno]][[dip]] <- "NA"
				}else{				
					# print(d)
					P.AA_DOS[[pheno]][[dip]] <- numeric(length(UNIQ))
					names(P.AA_DOS[[pheno]][[dip]]) <- UNIQ
					B.AA_DOS[[pheno]][[dip]] <- list()
					for ( u in 1:length(UNIQ) ) {
						uniq <- as.character(UNIQ[u])
						tag <- paste(dip,uniq,sep="_") ; tag <- gsub( "*","star",tag,fixed=T ) ; tag <- gsub(":","",tag,fixed=T)
						dip_col <- aggregate( MG.DIP[,dip]==uniq, by=list(ID=MG.DIP$ID), FUN=function(x) length(which(x==T)) )
						TEMP <- merge( MG.DIP[,c("ID",pheno,covs)], dip_col, by="ID" ) ; colnames(TEMP)[ncol(TEMP)] <- tag
						TEMP <- TEMP[which(!duplicated(TEMP$ID)),]
						formula <- as.formula(paste( pheno,"~",cov,"+",tag ))
						if ( pheno %in% c("RF_ACPA","ACPA","RF") ) {
							MOD <- glm( formula, data=TEMP, family=binomial(logit) )
							P.AA_DOS[[pheno]][[dip]][uniq] <- summary(MOD)$coefficients[tag,"Pr(>|z|)"]
						}else{
							MOD <- lm( formula, data=TEMP )
							P.AA_DOS[[pheno]][[dip]][uniq] <- anova(MOD)[tag,"Pr(>F)"]
						}
						B.AA_DOS[[pheno]][[dip]][[uniq]] <- MOD
					}
				}
			} # Close Dip Loop # ; print(paste("Done with",pheno))
		} # Close Pheno Loop 
		P$AA_DOS[[gene]] <- P.AA_DOS
		B$AA_DOS[[gene]] <- B.AA_DOS

		# Moving On
		# print(paste("Done with",gene))
	}
	P.precise[[paste("Dig",PRECISE,sep="_")]] <- P
	B.precise[[paste("Dig",PRECISE,sep="_")]] <- B
}

## Save Output
save(P.precise,file=paste(PathToPlot,"0-P_Precise.Rdata",sep=""))
# load(paste(PathToPlot,"0-P_Precise.Rdata",sep=""))

#############################################################
## PLOT ASSOCIATION RESULTS #################################
#############################################################

## Write Function to Pull out all Results for a given Gene
 # Includes Results for all tested Phenotypes
 # Specify Precision (by inputting only List w/ specific Precision)
   # e.g., DAT = P.precise$Dig_2
PLOT_RESULTS <- function( DAT, gene ) {
	## Pull Results
	 # Pull TYP ANOVA Results
	TYP_AOV <- DAT$TYP[[gene]]
	 # Pull TYP Additive Results
	TYP_DOS <- DAT$TYP_DOS[[gene]]
	 # Pull AA ANOVA Results
	AA_AOV <- DAT$AA_AOV[[gene]]
	 # Pull AA Additive Results
	AA_DOS <- DAT$AA_DOS[[gene]]

	## Plot Results ##
	PHENOS <- rownames(TYP_DOS) ; N.PHENOS <- length(PHENOS)
	COLS.list <- c("firebrick1","chocolate1","gold2","springgreen2","cadetblue2","steelblue2","slateblue3")
	COLS.ph <- colorRampPalette(COLS.list)(N.PHENOS)
	COLS.aa <- rep(colorRampPalette(COLS.list)(4),3)[1:N.PHENOS]
	COLS.AA <- c(colorRampPalette(COLS.list)(26),"black","grey90","grey50") ; names(COLS.AA) <- c(LETTERS,"*",".","?")
	PCH.ph <- rep(1,N.PHENOS)
	PCH.ph[grep("BL_MN",PHENOS)] <- 2
	PCH.ph[grep("RF|ACPA",PHENOS)] <- 3

	## HLA-Type Results Plot
	TYP <- cbind( TYP_AOV, TYP_DOS ) ; colnames(TYP)[1] <- "ANOVA"
	TYP <- cbind( TYP_AOV, TYP_DOS[,order(colnames(TYP_DOS))] ) ; colnames(TYP)[1] <- "ANOVA"
	YLIM <- c(0, -log10(min(TYP/30)) )
	MAIN <- paste("HLA Haplotype Association:",gene)
	png( paste(PathToPlot,"3-HapAssoc_",gene,".png",sep=""),height=800,width=1000+50*ncol(TYP),pointsize=28 )
	barplot( -log10(TYP), beside=T, las=2,col=COLS.ph,ylim=YLIM,main=MAIN,ylab="-log10(p)")
	abline( h=seq(0,20,1),lty=3,col="grey50",lwd=1 )
	abline( h=-log10(.05/(10^(0:6))),lty=2,col="magenta2",lwd=2 )
	legend( "topright", legend=rownames(TYP),title="Phenotype",fill=COLS.ph, ncol=ceiling(nrow(TYP)/4),cex=.9 )
	barplot( -log10(TYP), beside=T, las=2,col=COLS.ph,add=T )
	dev.off()

	## Amino Acid Level Results Plot(s)
	ISNA <- which( apply(AA_AOV,2,function(x) all(is.na(x)) ) | colnames(AA_AOV)=="TYP" )
	 # ANOVA
	AA_AOV.2 <- AA_AOV[ ,-ISNA ]
	YLIM <- c(0, -log10(min(AA_AOV.2/30)) )
	XVALS <- gsub("Pos_","",colnames(AA_AOV.2))
	XVALS <- as.numeric( gsub(".","-",XVALS,fixed=T) )
	XLIM <- range(XVALS)
	MAIN <- paste("Amino Acid ANOVA: HLA",gene)
	 # Open File & Format Layout
	HEIGHT <- 800
	WIDTH <- 2000+5*ncol(AA_AOV.2)
	PLOT_RATIO <- c( (WIDTH-HEIGHT)/WIDTH, HEIGHT/WIDTH )
	png( paste(PathToPlot,"4-AA_ANOVA_",gene,".png",sep=""),height=HEIGHT,width=WIDTH,pointsize=28 )
	layout( matrix(c(1,2), 1, 2, byrow = TRUE), widths=PLOT_RATIO ) # layout( matrix(c(1,2,3,4), 2, 2, byrow = TRUE), widths=c(3,2), heights=c(1,1) )
	 # "Manhattan Style Plot"
	plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, main=MAIN,xlab="Amino Acid",ylab="-log10(p)",xaxt="n")
	axis( 1, at=seq(-1000,1000,20),las=2 )	
	abline( v=seq(-1000,1000,20),lty=3,col="grey50")
	abline( h=seq(0,20,1),lty=3,col="grey50",lwd=1 )
	abline( h=-log10(.05/(10^(0:6))),lty=2,col="magenta2",lwd=2 )
	legend( "topright", legend=rownames(AA_AOV.2),title="Phenotype",col=COLS.aa,pch=PCH.ph, ncol=2,cex=.8,pt.lwd=2 )
	for ( p in 1:N.PHENOS ) {
		pheno <- PHENOS[p]
		points( XVALS,-log10(AA_AOV.2[pheno,]), col=COLS.aa[p],pch=PCH.ph[p],lwd=2 )
	}
	 # QQ Plot
	plot( 0,0,type="n",xlim=YLIM,ylim=YLIM, main="QQ P-Vals",xlab="-log10(Exp)",ylab="-log10(Exp)")
	abline( h=seq(0,20,1),lty=3,col="grey50",lwd=1 )
	abline( v=seq(0,20,1),lty=3,col="grey50",lwd=1 )
	abline( h=-log10(.05/(10^(0:6))),lty=2,col="magenta2",lwd=2 )
	abline( 0,1, lwd=2, col="black" )
	legend( "bottomright", legend=rownames(AA_AOV.2),title="Phenotype",col=COLS.aa,pch=PCH.ph, ncol=2,cex=.5,pt.lwd=2 )
	for ( p in 1:N.PHENOS ) {
		pheno <- PHENOS[p]
		exp <- 1:ncol(AA_AOV.2) / ncol(AA_AOV.2)
		points( -log10(exp),-log10(sort(AA_AOV.2[pheno,])), col=COLS.aa[p],pch=PCH.ph[p],lwd=2 )
	}
	dev.off()
	 # Additive Model
	for ( p in 1:N.PHENOS ) {
		pheno <- PHENOS[p]
		ISNA <- which( unlist(lapply(AA_DOS[[pheno]],function(x) all(x=="NA") )) )
		AAP <- AA_DOS[[pheno]][-ISNA]
		AAP <- AAP[which(names(AAP)!="TYP")]
		YLIM <- c(0, -log10(as.numeric(Reduce( min, AAP ))/30) )
		XVALS <- gsub("Pos_","",names(AAP))
		XVALS <- as.numeric( gsub(".","-",XVALS,fixed=T) )
		XLIM <- range(XVALS)
		MAIN <- paste("Amino Acid Additive Model: HLA",gene,"-",pheno)
		HEIGHT <- 800
		WIDTH <- 2000+5*ncol(AA_AOV.2) + 200
		PLOT_RATIO <- c( (WIDTH-HEIGHT)/WIDTH, HEIGHT/WIDTH, 200/WIDTH )
		png( paste(PathToPlot,"5-AA_Add_",gene,"_",pheno,".png",sep=""),height=HEIGHT,width=WIDTH,pointsize=34 )
		layout( matrix(c(1,2,3), 1, 3, byrow = TRUE), widths=PLOT_RATIO ) # layout( matrix(c(1,2,3,4), 2, 2, byrow = TRUE), widths=c(3,2), heights=c(1,1) )
		 # Manhattan Style Plot
		plot( 0,0,type="n",xlim=XLIM,ylim=YLIM, main=MAIN,xlab="Amino Acid",ylab="Frequency",xaxt="n")
		axis( 1, at=seq(-1000,1000,20),las=2 )
		abline( v=seq(-1000,1000,20),lty=3,col="grey50")
		abline( h=seq(0,20,1),lty=3,col="grey50",lwd=1 )
		abline( h=-log10(.05/(10^(0:6))),lty=2,col="magenta2",lwd=2 )
		obs <- unlist(AAP)
		abline( h=-log10(.05/length(obs)),lty=3,col="firebrick2",lwd=3)
		lapply( 1:length(AAP), function(x) points(rep(XVALS[x],length(AAP[[x]])),-log10(AAP[[x]]), pch=names(AAP[[x]]),col=COLS.AA[names(AAP[[x]])] ) )
		 # QQ Plot
		plot( 0,0,type="n",xlim=YLIM,ylim=YLIM, main="QQ P-Vals",xlab="-log10(Exp)",ylab="-log10(Exp)")
		abline( h=seq(0,20,1),lty=3,col="grey50",lwd=1 )
		abline( v=seq(0,20,1),lty=3,col="grey50",lwd=1 )
		abline( h=-log10(.05/(10^(0:6))),lty=2,col="magenta2",lwd=2 )
		abline( 0,1, lwd=2, col="black" )
		exp <- 1:length(unlist(AAP)) / length(unlist(AAP))
		obs <- unlist(AAP)
		abline( h=-log10(.05/length(obs)),lty=3,col="firebrick2",lwd=3)
		# points( -log10(exp),-log10(sort(obs)), col=COLS.AA[unlist(lapply(AAP,names))][order(obs)],pch=PCH.ph[p],lwd=2 )
		points( -log10(exp),-log10(sort(obs)), col=COLS.AA[unlist(lapply(AAP,names))][order(obs)],pch=unlist(lapply(AAP,names))[order(obs)],lwd=2 )
		 # Amino Acid Key
		Which_AA <- which(names(COLS.AA)%in% pat_aa )
		COLS.aa <- COLS.AA[Which_AA]
		barplot( matrix(rep(1,length(COLS.aa)),ncol=1), beside=F,col=COLS.aa,xaxt="n",yaxt="n",ylab="AA")
		axis( 2, at=1:length(COLS.aa)-.5, label=names(COLS.aa),las=2,cex.axis=.8 )
		dev.off()
	}
}
# PLOT_RESULTS( P.precise$Dig_4, "C" )
# PLOT_RESULTS( P.precise$Dig_4, "DRB1" )

## Run Through Genes
WHICH_GENES <- c("A","B","C","DQB1","DRB1","DPA1","DPB1","DQA1")
WHICH_GENES <- c("A","B","C","DQB1","DRB1")
for ( gene in WHICH_GENES ) {
	PLOT_RESULTS( P.precise$Dig_4, gene )
}

## Run Through Genes
WHICH_GENES <- c("A","B","C","DQB1","DRB1")
for ( gene in WHICH_GENES ) {
	PLOT_RESULTS( P.precise$Dig_2, gene )
}

#############################################################
## POKE AROUND SPECIFIC ALLELES #############################
#############################################################

##########################################
## DRB1 ##################################
##########################################

## Function to do Haplotype Analysis of Specified Amino Acid Positions
HAP_AN <- function( Positions) {
	## Pull together Haplotypes from Specified Positions
	HAP <- apply( PAT_AA$DRB1[,paste("Pos",Positions,sep="_")], 1, function(x) paste(x,collapse="") )
	HAP <- gsub("NA","-",HAP)
	HAP.uniq <- sort(unique(HAP))
	N.HAP <- length(HAP.uniq)
	 # Get Haplotype Frequencies
	HAP.freq <- table(HAP)
	HAP.rare <- names(HAP.freq)[which(HAP.freq < 15 )]
	 # Convert to Array for Additive Analysis
	HAP.arr <- array( 0,c(N.PATS,length(HAP.uniq)) )
	colnames(HAP.arr) <- HAP.uniq
	rownames(HAP.arr) <- PATS
	for ( pat in PATS ) {
		HAP.pat <- HAP[ grep(pat,names(HAP)) ]
		if ( HAP.pat[1]==HAP.pat[2] ) { HAP.arr[pat,HAP.pat[1]] <- 2
		}else{ HAP.arr[pat,HAP.pat] <- 1 }
	}
	 # Plot Haplotype Frequency
	png( paste(PathToPlot,"/DRB1_",paste(Positions,collapse=""),"_1-HapFreq.png",sep=""), height=800,width=1600,pointsize=30 )
	barplot( table(HAP),las=2,col="dodgerblue2",main=paste("HLA-DRB1: Pos",paste(Positions,collapse=","),"Haplotype Frequency"),ylab="# Haplotypes")
	abline(h=seq(0,1000,20),lty=3,col="grey50")
	barplot( table(HAP),las=2,col="dodgerblue2",add=T)
	dev.off()

	## Combine Haplotype & Phenotype Data
	MG.HAP <- merge( HAP.arr, FT, by.x="row.names",by.y="ID" )
	MG.HAP2 <- merge( data.frame(Samp=sapply(strsplit(names(HAP),"_"),"[",1),HAP), FT, by.x="Samp",by.y="ID" )

	## Association w/ Additive by Haplotype
	LM.HAP.AA <- LM.HAP.ANV <- P.HAP.AA <- list()
	P.HAP.ANV <- numeric(length(PHENOS)) ; names(P.HAP.ANV) <- PHENOS
	for ( p in 1:length(PHENOS) ) {
		pheno <- PHENOS[p]
		cov <- COVS[p]
		# ANOVA
		formula <- as.formula(paste( pheno,"~",cov,"+HAP" ))
		if ( pheno %in% c("RF_ACPA","ACPA","RF") ) {
			LM.HAP.ANV[[pheno]] <- MOD <- chisq.test( table( MG.HAP2[,pheno],as.character(MG.HAP2[,"HAP"]) ) ) # glm( formula, data=MG.HAP2, family=binomial(logit), subset=which(!(HAP%in%HAP.rare)) )
			P.HAP.ANV[pheno] <- LM.HAP.ANV[[pheno]]$p.value
		}else{
			# ANOVA
			LM.HAP.ANV[[pheno]] <- MOD <- lm( formula, data=MG.HAP2, subset=which(!(HAP%in%HAP.rare)) )
			P.HAP.ANV[pheno] <- anova(MOD)["HAP","Pr(>F)"]
			# Plot Residuals
			formula <- as.formula(paste( pheno, "~", cov))
			RESID <- resid(lm( formula, data=MG.HAP2))
			png( paste(PathToPlot,"/DRB1_",paste(Positions,collapse=""),"_2-BoxPlot_",pheno,".png",sep=""), height=1400,width=1600,pointsize=30 )
			par(mfrow=c(2,1))
			barplot( table(HAP),las=2,col="dodgerblue2",main=paste("HLA-DRB1: Pos",paste(Positions,collapse=","),"Haplotype Frequency"),ylab="# Haplotypes")
			abline(h=seq(0,1000,20),lty=3,col="grey50")
			barplot( table(HAP),las=2,col="dodgerblue2",add=T)
			boxplot( RESID ~ MG.HAP2$HAP[as.numeric(names(RESID))], las=2,col="chartreuse2",main=paste("HLA-DRB1: Pos11,71,74 Haplotype vs",pheno),ylab=paste("Resid: vs",cov) )
			abline(h=seq(-10,10,1),lty=3,col="grey50")
			boxplot( RESID ~ MG.HAP2$HAP[as.numeric(names(RESID))],add=T,las=2,col="chartreuse2" )
			points( RESID ~ MG.HAP2$HAP[as.numeric(names(RESID))],pch="+" )
			dev.off()
		}
		# Additive
		LM.HAP.AA[[pheno]] <- list()
		P.HAP.AA[[pheno]] <- numeric(N.HAP) ; names(P.HAP.AA[[pheno]]) <- HAP.uniq
		for ( hap in HAP.uniq ) {
			if ( all( strsplit(hap,"")[[1]]=="-" ) ) { next }
			formula <- as.formula(paste( pheno,"~",cov,"+",hap ))
			if ( pheno %in% c("RF_ACPA","ACPA","RF") ) {
				formula <- as.formula(paste( pheno,"~",hap ))
				LM.HAP.AA[[pheno]][[hap]] <- MOD <- glm( formula, data=MG.HAP, family=binomial(logit) )
				P.HAP.AA[[pheno]][hap] <- summary(MOD)$coefficients[hap,"Pr(>|z|)"]
			}else{
				LM.HAP.AA[[pheno]][[hap]] <- MOD <- lm( formula, data=MG.HAP )
				P.HAP.AA[[pheno]][hap] <- anova(MOD)[hap,"Pr(>F)"]
			}
		} # Close Haplotype Loop

	} # Close Pheno Loop

	## Plot Results
	P.comp <- cbind( Reduce( rbind, P.HAP.AA ), P.HAP.ANV )
	colnames(P.comp)[ncol(P.comp)] <- "ANOVA"
	rownames(P.comp) <- names(P.HAP.AA)
	P.comp <- P.comp[,-which(colnames(P.comp)==paste(rep("-",length(Positions)),collapse=""))]

	COLS.list <- c("firebrick1","chocolate1","gold2","springgreen2","cadetblue2","steelblue2","slateblue3")
	COLS.ph <- colorRampPalette(COLS.list)(nrow(P.comp))
	YLIM <- c( 0,-log10(min(P.comp)/30) )
	png( paste(PathToPlot,"/DRB1_",paste(Positions,collapse=""),"_3-HapAssoc.png",sep=""), height=800,width=1600,pointsize=30 )
	barplot( -log10(P.comp), beside=T, col=COLS.ph,ylim=YLIM,las=2,main=paste("HLA-DRB1: Pos",paste(Positions,collapse=","),"Haplotype ANOVA"),ylab="-log10(p)")
	abline(h=seq(0,20,1),lty=3,col="grey50" )
	abline( h=-log10(.05/(10^(0:6))),lty=2,col="magenta2",lwd=2 )
	legend( "topright",legend=rownames(P.comp),fill=COLS.ph,ncol=3,cex=.8)
	barplot( -log10(P.comp), beside=T, col=COLS.ph,ylim=YLIM,las=2,add=T )
	dev.off()
	## Compile Outputs
	# OUT <- list( MOD.AA=LM.HAP.AA, MOD.ANV=LM.HAP.ANV, P.AA=P.HAP.AA, P.ANV=P.HAP.ANV )
	OUT <- list( MOD.AA=LM.HAP.AA, MOD.ANV=LM.HAP.ANV, P=P.comp )
}

OUT <- list()
##########################################
## POS 11, 71, 74 ##
 # Viatte, et al (2015)
Positions <- c(11,71,74)
OUT$p117174 <- HAP_AN(Positions)

##########################################
## POS 11, 13, 71, 74 ##
Positions <- c(11,13,71,74)
OUT$p11137174 <- HAP_AN(Positions)

##########################################
## POS 11, 13, 71, 74 ##
Positions <- 70:74
OUT$pSE <- HAP_AN(Positions)

BETAS <- SES <- PS <- list()
for ( hap in names(OUT) ) {
	PS[[hap]] <- matrix( unlist(lapply( PHENOS[1:8], function(x)unlist(lapply( names(OUT[[hap]]$MOD.AA[[x]]), function(y)summary(OUT[[hap]]$MOD.AA[[x]][[y]])$coefficients[y,4] )) )), byrow=F,ncol=8 )
	BETAS[[hap]] <- matrix( unlist(lapply( PHENOS[1:8], function(x)unlist(lapply( names(OUT[[hap]]$MOD.AA[[x]]), function(y)summary(OUT[[hap]]$MOD.AA[[x]][[y]])$coefficients[y,"Estimate"] )) )), byrow=F,ncol=8 )
	SES[[hap]] <- matrix( unlist(lapply( PHENOS[1:8], function(x)unlist(lapply( names(OUT[[hap]]$MOD.AA[[x]]), function(y)summary(OUT[[hap]]$MOD.AA[[x]][[y]])$coefficients[y,"Std. Error"] )) )), byrow=F,ncol=8 )
	colnames(BETAS[[hap]]) <- colnames(SES[[hap]]) <- colnames(PS[[hap]]) <- PHENOS[1:8]
	rownames(BETAS[[hap]]) <- rownames(SES[[hap]]) <- rownames(PS[[hap]]) <- names(OUT[[hap]]$MOD.AA[[1]])
	# pairs(BETAS[[hap]])
	for ( p in 1:4 ) {
		y_pheno <- PHENOS[p]
		x_pheno <- PHENOS[p+4]
		COLS.beta <- c("mediumpurple3","tomato2")
		png( paste(PathToPlot,"DRB1-BETA.",hap,".",x_pheno,".png",sep=""),height=1600,width=1600,pointsize=36 )
		XLIM <- extendrange(BETAS[[hap]][,x_pheno])
		YLIM <- extendrange(BETAS[[hap]][,y_pheno])
		plot( BETAS[[hap]][,x_pheno],BETAS[[hap]][,y_pheno], pch=21,col=COLS.beta[1],bg=COLS.beta[2], xlim=XLIM,ylim=YLIM,main=paste("Beta Estimates of",y_pheno,"vs",x_pheno,"-",hap),xlab=paste("Beta:",x_pheno),ylab=paste("Beta:",y_pheno) )
		abline( h=seq(-5,5,.5),v=seq(-5,5,.5),lty=3,col="grey50",lwd=1 )
		abline( h=0,v=0,lty=1,col="grey50",lwd=1 )
		arrows( BETAS[[hap]][,x_pheno]+SES[[hap]][,x_pheno],BETAS[[hap]][,y_pheno],BETAS[[hap]][,x_pheno]-SES[[hap]][,x_pheno],BETAS[[hap]][,y_pheno], code=3,angle=90,lwd=5,col=COLS.beta[1] )
		arrows( BETAS[[hap]][,x_pheno],BETAS[[hap]][,y_pheno]+SES[[hap]][,y_pheno],BETAS[[hap]][,x_pheno],BETAS[[hap]][,y_pheno]-SES[[hap]][,y_pheno], code=3,angle=90,lwd=5,col=COLS.beta[1] )
		text( BETAS[[hap]][,x_pheno],BETAS[[hap]][,y_pheno]+.02*diff(YLIM), label=rownames(BETAS[[hap]]), col=COLS.beta[1], pos=4,cex=1.2 )
		points( BETAS[[hap]][,x_pheno],BETAS[[hap]][,y_pheno], pch=21,col=COLS.beta[1],bg=COLS.beta[2], lwd=5,cex=1.2 )
		MOD <- lm( BETAS[[hap]][,y_pheno]~BETAS[[hap]][,x_pheno] )
		abline(MOD,lwd=6,lty=2,col=COLS.beta[2] )
		text( quantile(XLIM,.05),quantile(YLIM,.02), label=paste("p=",formatC(summary(MOD)$coefficients[length(coef(MOD)),4],digits=2,format="e"),sep=""), col=COLS.beta[2],cex=1.2 )
		dev.off()
	}	
}

BETAS.117174 <- matrix( unlist(lapply( PHENOS[1:8], function(x)unlist(lapply( names(OUT.117174$MOD.AA[[x]]), function(y)summary(OUT.117174$MOD.AA[[x]][[y]])$coefficients[y,1] )) )), byrow=F,ncol=8 )
colnames(BETAS.117174) <- PHENOS[1:8]
pairs(BETAS.117174)

BETAS.11137174 <- matrix( unlist(lapply( PHENOS[1:8], function(x)unlist(lapply( names(OUT.11137174$MOD.AA[[x]]), function(y)summary(OUT.11137174$MOD.AA[[x]][[y]])$coefficients[y,1] )) )), byrow=F,ncol=8 )
colnames(BETAS.11137174) <- PHENOS[1:8]
pairs(BETAS.11137174)

BETAS.SE <- matrix( unlist(lapply( PHENOS[1:8], function(x)unlist(lapply( names(OUT.SE$MOD.AA[[x]]), function(y)summary(OUT.SE$MOD.AA[[x]][[y]])$coefficients[y,"Estimate"] )) )), byrow=F,ncol=8 )
SES.SE <- matrix( unlist(lapply( PHENOS[1:8], function(x)unlist(lapply( names(OUT.SE$MOD.AA[[x]]), function(y)summary(OUT.SE$MOD.AA[[x]][[y]])$coefficients[y,"Std. Error"] )) )), byrow=F,ncol=8 )
colnames(BETAS.SE) <- colnames(SES.SE) <- PHENOS[1:8]
rownames(BETAS.SE) <- rownames(SES.SE) <- names(OUT.SE$MOD.AA[[1]])
pairs(BETAS.SE)















#############################################################
## END OF DOC ###############################################
#############################################################





# COLS.list <- c("firebrick2","chocolate2","gold1","springgreen2","steelblue2","slateblue3","black")
# COLS <- colorRampPalette(rev(COLS.list))(100)
# heatmap.2(-log10(P.DIP), scale="none",trace="none",Colv=NA,Rowv=NA,dendrogram="none",col=COLS)

