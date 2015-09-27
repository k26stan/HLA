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
PathToAA <- "/Users/kstandis/Data/Burn/Data/HLA/20150512_SOAP_HLA_Types/20150512_HLA_AA.Rdata"
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

## Get Phenotype Info
FT <- read.table( PathToFT, sep="\t",header=T )

#############################################################
## ORGANIZE DATA ############################################
#############################################################

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
for ( g in 1:N.GENE ) { gene <- GENE_LIST[g]
	FREQ.4 <- rowSums( PAT_DOS.4[[gene]] ) ; FREQ.4 <- FREQ.4[order(names(FREQ.4))]
	FREQ.2 <- rowSums( PAT_DOS.2[[gene]] ) ; FREQ.2 <- FREQ.2[order(names(FREQ.2))]
	png( paste(PathToPlot,"1-HapFreq.",gene,".png",sep=""), height=1200,width=1600,pointsize=30 )
	par(mfrow=c(2,1))
	# 4 digit
	barplot( FREQ.4, las=2, col=COLS[1], main=paste("HLA -",gene,": 4-Digit Haplotype Frequency") )
	barplot( FREQ.2, las=2, col=COLS[2], main=paste("HLA -",gene,": 2-Digit Haplotype Frequency") )
	dev.off()
}

## Plot AA Level Diversity within Cohort
COLS.list <- c("firebrick2","chocolate2","gold1","chartreuse2","cadetblue2","dodgerblue2","slateblue2","magenta2")
# COLS.2.list <- sample(COLS.list)
COLS.AA <- c(colorRampPalette(COLS.list)(26),"black","grey80") ; names(COLS.AA) <- c(LETTERS,"*",".")
for ( g in 1:length(GENE_LIST) ) {
	gene <- GENE_LIST[g]
	# Pull out Patient Data for Gene
	pat_typ <- PAT_TYP[gene,]
	pat_dos <- PAT_DOS[[gene]]
	pat_aa <- PAT_AA[[gene]]
	# Get Amino Acid Frequencies
	max.temp <- max( unlist(lapply( apply( pat_aa, 2, table ), length)) )
	which.temp <- which( unlist(lapply( apply( pat_aa, 2, table ), length)) > 1 )
	pat_tab <- sapply(apply( pat_aa, 2, table ), "[", 1:max.temp ) ; colnames(pat_tab) <- gsub("Pos_","",colnames(pat_tab))
	pat_tab.prc <- pat_tab / nrow(pat_aa)
	pat_tab.names <- sapply(apply( pat_aa, 2, function(x) names(table(x)) ), "[", 1:max.temp )
	 # Plot
	# which.temp.Pos <- as.numeric(as.character(gsub("Pos_","",names(which.temp))))
	# XLIM <- range(which.temp.Pos)
	XLIM <- range(as.numeric(colnames(pat_tab)))
	png( paste(PathToPlot,"/2-AAfreq_",gene,".png",sep=""), height=1000,width=2000,pointsize=30 )
	plot( 0,0,type="n",xlim=XLIM,ylim=c(0,1), main=paste("Amino Acid Frequency",gene),xlab="Amino Acid",ylab="Frequency",xaxt="n")
	axis( 1, at=seq(-1000,1000,20),las=2 )	
	abline( v=seq(-1000,1000,20),lty=3,col="grey50")
	allele_bar <- function(x) {
		xval <- as.numeric(colnames(pat_tab)[x])
		yvals <- cumsum(pat_tab.prc[,x])
		barplot( t(t(pat_tab.prc[,x])),add=T,xaxt="n",beside=F,border=NA,col=COLS.AA[pat_tab.names[,x]],space=c(xval,0),width=1)
	}
	DWAI <- lapply( which.temp, allele_bar )
	dev.off()
	# barplot( pat_tab[,which.temp]/nrow(pat_aa), col=COLS.2.list,border=NA,las=2, main=paste("Amino Acid Frequency",gene),xlab="Amino Acid",ylab="Frequency" )
	# barplot( pat_tab[,which.temp]/nrow(pat_aa), col=,border=NA,las=2, main=paste("Amino Acid Frequency",gene),xlab="Amino Acid",ylab="Frequency" )
}
	
#############################################################
## RUN ASSOCIATION TESTS ####################################
#############################################################

## Specify Phenos/Covs
PHENOS <- c("DEL_MNe_MN","DEL_lCRP_MNe_MN","DEL_rSJC_MNe_MN","DEL_rTJC_MNe_MN",
	"DAS_BL_MN","lCRP_BL_MN","rSJC_BL_MN","rTJC_BL_MN",
	"RF_ACPA","ACPA","RF" )
COVS <- c("DAS_BL_MN","lCRP_BL_MN","rSJC_BL_MN","rTJC_BL_MN",
	"RF_ACPA+log(DIS_DUR)","RF_ACPA+log(DIS_DUR)","RF_ACPA+log(DIS_DUR)","RF_ACPA+log(DIS_DUR)",
	"","","")
PH.COV <- data.frame(PHENOS,COVS)

## Loop Through Genes & Test for Assoc
P.precise <- list()
for ( PRECISE in c(2,4) ) {
	P <- list()
	P$DIP <- P$DOS <- P$AA <- list()
	for ( g in 1:N.GENE ) {
		gene <- GENE_LIST[g]
		## Pull out Patient Data for Gene
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
		HAP.freqs <- rowSums(pat_dos.a)
		which_common <- which(HAP.freqs>.05*N.PATS)
		which_common <- which(HAP.freqs>10)
		which_common.haps <- names(which_common)
		pat_dos <- pat_dos.a[ which_common, ]
		if ( nrow(pat_dos.a)==1 ) { next }

		## Merge Files
		 # Merge Diplotype Tables w/ Pheno
		MG.DIP.1 <- merge( pat_aa, data.frame(PAT=names(pat_typ),TYP=pat_typ), by.x="row.names",by.y="PAT" )
		colnames(MG.DIP.1)[1] <- "DIP"
		MG.DIP.1 <- data.frame( ID=sapply(strsplit(MG.DIP.1[,1],"_"),"[",1), MG.DIP.1 )
		MG.DIP <- merge( FT, MG.DIP.1, by="ID")
		MG.DIP <- MG.DIP[ which(MG.DIP$TYP %in% which_common.haps), ]
		MG.DIP[,"ACPA"] <- MG.DIP[,"ACPA"]=="Positive"
		 # Merge w/ Phenotype Data
		MG.DOS <- merge( FT, t(pat_dos), by.x="ID",by.y="row.names" )

		## Run Association Tests

		 # ANOVA of Common Haplotypes vs Phenotypes
		P.TYP <- numeric( length(PHENOS) )
		names(P.TYP) <- PHENOS
		for ( p in 1:length(PHENOS) ) {
			pheno <- PHENOS[p]
			cov <- COVS[p]
			formula <- as.formula(paste( pheno,"~",cov,"+TYP" ))
			if ( pheno %in% c("RF_ACPA","ACPA","RF") ) {
				MOD <- chisq.test( table( MG.DIP[,pheno],as.character(MG.DIP[,"TYP"]) ) )
				P.TYP[p] <- MOD$p.value
			}else{
				MOD <- lm( formula, data=MG.DIP )
				P.TYP[p] <- anova(MOD)["TYP","Pr(>F)"]
			}
		}

		 # ANOVA for all 
		P.DIP <- array( ,c(length(PHENOS),ncol(MG.DIP.1)-2))
		colnames(P.DIP) <- colnames(MG.DIP.1)[3:ncol(MG.DIP.1)]
		rownames(P.DIP) <- PHENOS
		for ( p in 1:length(PHENOS) ) {
			pheno <- PHENOS[p]
			cov <- COVS[p]
			for ( d in 1:ncol(P.DIP) ) {
				dip <- colnames(P.DIP)[d]
				UNIQ <- unique(MG.DIP[,dip]) ; UNIQ <- UNIQ[which(!is.na(UNIQ))]
				if ( length(UNIQ)==1 ) {
					P.DIP[p,d] <- NA
				}else{
					formula <- as.formula(paste( pheno,"~",cov,"+",dip ))
					print(d)
					MOD <- lm( formula, data=MG.DIP )
					P.DIP[p,d] <- anova(MOD)["dip","Pr(>F)"]	
				}
				# TEMP_ARR <- MG.DIP[, c(pheno,cov,dip) ]
				# colnames(TEMP_ARR) <- c("pheno","cov","dip")
				# TEMP_ARR <- TEMP_ARR[which(!is.na(TEMP_ARR$pheno)),]
				# TEMP_ARR <- TEMP_ARR[which(!is.na(TEMP_ARR$dip)),]
				# TEMP_ARR <- TEMP_ARR[which(TEMP_ARR$dip!="*"),]
				# UNIQ <- unique(TEMP_ARR$dip) ; UNIQ <- UNIQ[which(!is.na(UNIQ))]		
			}
		}
		P$DIP[[gene]] <- P.DIP
		 # Dosage Association
		P.DOS <- array( ,c(4,nrow(pat_dos)))
		colnames(P.DOS) <- rownames(pat_dos)
		rownames(P.DOS) <- PHENOS
		for ( p in 1:length(PHENOS) ) {
			pheno <- PHENOS[p]
			cov <- COVS[p]
			for ( d in 1:ncol(P.DOS) ) {
				dip <- colnames(P.DOS)[d]
				TEMP_ARR <- MG.DOS[, c(pheno,cov,dip) ]
				colnames(TEMP_ARR) <- c("pheno","cov","dip")
				TEMP_ARR <- TEMP_ARR[which(!is.na(TEMP_ARR$pheno)),]
				TEMP_ARR <- TEMP_ARR[which(!is.na(TEMP_ARR$dip)),]
				UNIQ <- unique(TEMP_ARR$dip) ; UNIQ <- UNIQ[which(!is.na(UNIQ))]
				if ( length(UNIQ)==1 ) {
					P.DOS[p,d] <- NA
				}else{
					# print(d)
					MOD <- lm( pheno ~ cov + dip, data=TEMP_ARR )
					P.DOS[p,d] <- anova(MOD)["dip","Pr(>F)"]	
				}
			}
		}
		P$DOS[[gene]] <- P.DOS
		# Moving On
		print(paste("Done with",gene))
	}
	P.precise[[paste("Dig",PRECISE,sep="_")]] <- P
}



for ( PRECISE in c(2,4) ) {
	P <- P.precise[[paste("Dig",PRECISE,sep="_")]]
	COLS <- c("firebrick1","gold3","chartreuse1","dodgerblue1")
	for ( g in 1:length(P$DOS) ) {
		gene <- names(P$DOS)[g]
		P_DIP <- P$DIP[[gene]]
		P_DOS <- P$DOS[[gene]]
		## Plot
		png( paste(PathToPlot,"/AS_",gene,".",PRECISE,".png",sep=""), height=1600,width=2000,pointsize=30 )
		par(mfrow=c(2,1))
		 # Dip
		XLIM <- c(1,ncol(P_DIP))
		YLIM <- c(0,-log10(min(P_DIP,na.rm=T)))
		XAXS <- gsub("Pos_","",colnames(P_DIP))[c(seq(1,XLIM[2],10),XLIM[2])]
		plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,ylab="-log10(p)",xlab="Amino Acid",main=paste("Diploid Association w/ Response Phenotypes -",gene),xaxt="n" )
		axis(1, at=c(seq(1,XLIM[2],10),XLIM[2]), label=XAXS, las=2 )
		for ( p in 1:length(PHENOS) ) {
			points( -log10(P_DIP[p,]), col=COLS[p],pch="+" )
		}
		abline( h=seq(0,10,1),lty=2,col="gray50",lwd=1)
		abline( h=-log10(.05),lty=2,col="magenta2",lwd=2)
		 # Dos
		XLIM <- c(1,ncol(P_DOS))
		YLIM <- c(0,-log10(min(P_DOS,na.rm=T)))
		plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,ylab="-log10(p)",xlab="HLA Type",main=paste("HLA Type Dosage Association w/ Response Phenotypes -",gene),xaxt="n" )
		axis(1, at=1:XLIM[2], label=gsub("Pos_","",colnames(P_DOS)), las=2 )
		for ( p in 1:length(PHENOS) ) {
			points( -log10(P_DOS[p,]), col=COLS[p],pch="+" )
		}
		abline( h=seq(0,10,1),lty=2,col="gray50",lwd=1)
		abline( h=-log10(.05),lty=2,col="magenta2",lwd=2)
		legend( "bottomleft",ncol=2,legend=c("DAS","lCRP","rSJC","rTJC"),col=COLS,pch="+" )
		dev.off()
	}
}

#############################################################
## POKE AROUND SPECIFIC ALLELES #############################
#############################################################

##########################################
## DRB1 ##################################
##########################################

## POS 11, 71, 74 ##
 # Viatte, et al (2015)
Positions <- c(11,71,74)
HAP <- apply( PAT_AA$DRB1[,paste("Pos",Positions,sep="_")], 1, function(x) paste(x,collapse="") )
HAP <- gsub("NA","-",HAP)
N.HAP <- length(unique(HAP))
HAP.arr <- array( 0,c(N.PATS,length(unique(HAP))) )
colnames(HAP.arr) <- unique(HAP)
rownames(HAP.arr) <- PATS
for ( pat in PATS ) {
	HAP.pat <- HAP[ grep(pat,names(HAP)) ]
	if ( HAP.pat[1]==HAP.pat[2] ) { HAP.arr[pat,HAP.pat[1]] <- 2
	}else{ HAP.arr[pat,HAP.pat] <- 1 }
}
 # Plot Haplotype Frequency
png( paste(PathToPlot,"/DRB1_1-HapFreq",paste(Positions,collapse=""),".png",sep=""), height=800,width=1600,pointsize=30 )
barplot( table(HAP),las=2,col="dodgerblue2",main="HLA-DRB1: Pos11,71,74 Haplotype Frequency",ylab="# Haplotypes")
abline(h=seq(0,1000,20),lty=3,col="grey50")
barplot( table(HAP),las=2,col="dodgerblue2",main="HLA-DRB1: Pos11,71,74 Haplotype Frequency",ylab="# Haplotypes",add=T)
dev.off()

MG.HAP <- merge( HAP.arr, FT, by.x="row.names",by.y="ID" )

LM.HAP <- LM.HAP.sev <- list()
for ( hap in unique(HAP) ) {
	formula <- as.formula( paste("DEL_MNe_MN ~",hap,"+ DAS_BL_MN") )
	# formula <- as.formula( paste("DEL_MNe_MN ~",hap,"+ DAS_BL_MN + RF_ACPA") )
	# formula <- as.formula( paste("DEL_MNe_MN ~",hap,"+ DAS_BL_MN + I(ACPA=='Positive')") )
	# formula <- as.formula( paste("DEL_MNe_MN ~",hap,"+ DAS_BL_MN + RF") )
	LM.HAP[[hap]] <- lm( formula, data=MG.HAP )
	formula <- as.formula( paste("DAS_BL_MN ~",hap,"+ RF_ACPA+log(DIS_DUR)") )
	LM.HAP.sev[[hap]] <- lm( formula, data=MG.HAP )
	# lm( DEL_MNe_MN ~ x + DAS_BL_MN, data=MG.HAP )
}
lapply( LM.HAP, anova )
Ps <- unlist(lapply( LM.HAP, function(x) anova(x)[1,"Pr(>F)"] ))
Ps <- unlist(lapply( LM.HAP.sev, function(x) anova(x)[1,"Pr(>F)"] ))
Ps.sev <- unlist(lapply( LM.HAP.sev, function(x) anova(x)[1,"Pr(>F)"] ))

Bs <- unlist(lapply( LM.HAP, function(x) summary(x)$coef[2,"Estimate"] ))
Bs.sev <- unlist(lapply( LM.HAP.sev, function(x) summary(x)$coef[2,"Estimate"] ))
# unlist(lapply( LM.HAP, function(x) anova(x)["RF_ACPA","Pr(>F)"] ))
# unlist(lapply( LM.HAP, function(x) anova(x)["RF","Pr(>F)"] ))
# unlist(lapply( LM.HAP, function(x) anova(x)['I(ACPA == "Positive")',"Pr(>F)"] ))
plot( -log10(1:17/17), -log10(sort(Ps)) )
abline(0,1)

## POS 11, 13, 71, 74 ##
Positions <- c(11,13,71,74)
HAP <- apply( PAT_AA$DRB1[,paste("Pos",Positions,sep="_")], 1, function(x) paste(x,collapse="") )
HAP <- gsub("NA","-",HAP)
N.HAP <- length(unique(HAP))
HAP.arr <- array( 0,c(N.PATS,length(unique(HAP))) )
colnames(HAP.arr) <- unique(HAP)
rownames(HAP.arr) <- PATS
for ( pat in PATS ) {
	HAP.pat <- HAP[ grep(pat,names(HAP)) ]
	if ( HAP.pat[1]==HAP.pat[2] ) {
		HAP.arr[pat,HAP.pat[1]] <- 2
	}else{
		HAP.arr[pat,HAP.pat] <- 1
	}
}
 # Plot Haplotype Frequency
png( paste(PathToPlot,"/DRB1_1-HapFreq",paste(Positions,collapse=""),".png",sep=""), height=800,width=1600,pointsize=30 )
barplot( table(HAP),las=2,col="dodgerblue2",main="HLA-DRB1: Pos11,71,74 Haplotype Frequency",ylab="# Haplotypes")
abline(h=seq(0,1000,20),lty=3,col="grey50")
barplot( table(HAP),las=2,col="dodgerblue2",main="HLA-DRB1: Pos11,71,74 Haplotype Frequency",ylab="# Haplotypes",add=T)
dev.off()









#############################################################
## END OF DOC ###############################################
#############################################################





# COLS.list <- c("firebrick2","chocolate2","gold1","springgreen2","steelblue2","slateblue3","black")
# COLS <- colorRampPalette(rev(COLS.list))(100)
# heatmap.2(-log10(P.DIP), scale="none",trace="none",Colv=NA,Rowv=NA,dendrogram="none",col=COLS)

