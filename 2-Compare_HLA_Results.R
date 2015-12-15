## Compare HLA Results Across Platforms ##
## May 22, 2015 ##
## Kristopher Standish ##

## Compare Results of HLA typing from multiple platforms
 # SOP = Reads -> SOAP-HLA
 # CHP = SNP Chip -> SNP2HLA
 # SEQ = HaplotypeCaller SNPs -> SNP2HLA
 # LAB = Lab Typing (Gold Standard)

#############################################################
## LOAD DATA ################################################
#############################################################
library(xlsx)
library(stringr)

## Set Date
DATE <- gsub("-","",Sys.Date())

## Set Paths to Data Sets & Save Locations (Mac)
PathToSOAP <- "/Users/kstandis/Data/Burn/Data/HLA/20151211_SOAP_HLA_Types/20151211_HLA_Types.Rdata"
PathToS2H <- "/Users/kstandis/Data/Burn/Data/HLA/20140924_SNP2HLA_Types/"
PathToLAB <- "/Users/kstandis/Data/Burn/Data/HLA/20150522_Lab_HLA/LabCorp HLA typing data_ART3001.xls"
PathToPheno <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20150520_Full_Table.txt"
PathToPlot <- paste("/Users/kstandis/Data/Burn/Plots/",DATE,"_HLA/",sep="")
PathToRefs <- "/Users/kstandis/Data/HLI_Phase/20150223_HLA_Ref/Alignments_Rel_3190/"
if ( !file.exists(PathToPlot) ) { dir.create( PathToPlot ) }

## Load HLA Types
 # SOAP-HLA
load( PathToSOAP )
SOP.l <- COMPILE
 # SNP2HLA
CHP.2.l <- read.table( paste(PathToS2H,"CHP_2.txt",sep=""),sep="\t",header=T,colClasses="character" )
CHP.4.l <- read.table( paste(PathToS2H,"CHP_4.txt",sep=""),sep="\t",header=T,colClasses="character" )
SEQ.2.l <- read.table( paste(PathToS2H,"HC2_2.txt",sep=""),sep="\t",header=T,colClasses="character" )
SEQ.4.l <- read.table( paste(PathToS2H,"HC2_4.txt",sep=""),sep="\t",header=T,colClasses="character" )
 # Lab Typing
LAB.l <- read.xlsx( PathToLAB, sheetIndex=1, rowIndex=1:101, header=T, colIndex=1:11)

## Load Phenotypes & Sample Lists
FT <- read.table( PathToPheno, sep="\t",header=T )

## Load Amino Acid Data
FILES <- list.files(PathToRefs)
gsub("_prot.txt","",FILES[grep("_prot",FILES)])

#############################################################
## ORGANIZE PATIENT DATA ####################################
#############################################################

## Set Platform Names
PLAT.names <- c("SOP","CHP","SEQ","LAB")

#### SAMPLE LISTS ####

## Get Sample List for each Platform
SAMP.all <- union( as.character(FT$ID), colnames(SOP.l$GENES.2.list$DRB1) )
SAMP.sop <- colnames( SOP.l$GENES.2.list$DRB1 )
SAMP.chp <- rownames( CHP.2.l )
SAMP.seq <- rownames( SEQ.2.l )
SAMP.lab <- unique(as.character( LAB.l$Accession[which( LAB.l$Accession %in% SAMP.all )] ))

## LAB: Specify Rows by Sample Name
Which_Samp <- which( LAB.l$Accession %in% SAMP.lab )
Subject_Key <- data.frame( ID=LAB.l$Accession[Which_Samp], NUM=LAB.l$Subject[Which_Samp] )
Subject_Key <- Subject_Key[ which(!duplicated(Subject_Key$ID)), ]
LAB.2 <- merge( Subject_Key, LAB.l, by.x="NUM",by.y="Subject" )
 # Problem w/ Subject Names W367072 & W367073
 # Listed as SAME sample in LAB typing, but different for phenotypes
LAB.2 <- LAB.2[ -which(LAB.2$ID=="W367073"), ]
SAMP.lab <- as.character(unique( LAB.2$ID ))

## Get Samples that Intersect ALL Platforms
SAMP.int <- Reduce( intersect, list(SAMP.sop,SAMP.chp,SAMP.seq,SAMP.lab) )

#### GENE LISTS ####

## Get Gene List for each Platform
GENE.sop <- sort( names( SOP.l$GENES.2.list ) )
GENE.chp <- sort( unique( sapply( colnames(CHP.2.l), function(x) substr(x,1,nchar(x)-1) ) ) ) # sort( gsub("HLA","", unique( sapply( colnames(CHP.2.l), function(x) substr(x,1,nchar(x)-1) ) )) )
GENE.seq <- sort( GENE.chp )
GENE.lab <- sort( unique(as.character(LAB.2$TestName)) ) # sort( gsub("HLA-","", unique(as.character(LAB.2$TestName)) ) )
GENE.all <- gsub("[0-9]", "", GENE.sop[2:length(GENE.sop)] ) ; GENE.all[grep("TAP",GENE.all)] <- c("TAP1","TAP2")
GENE.all <- gsub("[0-9]", "", GENE.sop )
GENE.tab <- array(, c(length(GENE.all),4) )
rownames(GENE.tab) <- GENE.all ; colnames(GENE.tab) <- PLAT.names
 # Fill Table
for ( row in 1:nrow(GENE.tab) ) {
	name <- rownames(GENE.tab)[row]
	# if ( row < 4 ) {
	which_sop <- union( which(GENE.sop==name | GENE.sop==paste("HLA",name,sep="") | GENE.sop==paste("HLA",name,sep="-") ), grep( paste(name,"[0-9]",sep=""), GENE.sop ) )
	which_chp <- union( which(GENE.chp==name | GENE.chp==paste("HLA",name,sep="") | GENE.chp==paste("HLA",name,sep="-") ), grep( paste(name,"[0-9]",sep=""), GENE.chp ) )
	which_seq <- union( which(GENE.seq==name | GENE.seq==paste("HLA",name,sep="") | GENE.seq==paste("HLA",name,sep="-") ), grep( paste(name,"[0-9]",sep=""), GENE.seq ) )
	which_lab <- union( which(GENE.lab==name | GENE.lab==paste("HLA",name,sep="") | GENE.lab==paste("HLA",name,sep="-") ), grep( paste(name,"[0-9]",sep=""), GENE.lab ) )
	GENE.tab[row,"SOP"] <- ifelse( length(which_sop)>0, GENE.sop[ which_sop ], NA )
	GENE.tab[row,"CHP"] <- ifelse( length(which_chp)>0, GENE.chp[ which_chp ], NA )
	GENE.tab[row,"SEQ"] <- ifelse( length(which_seq)>0, GENE.seq[ which_seq ], NA )
	GENE.tab[row,"LAB"] <- ifelse( length(which_lab)>0, GENE.lab[ which_lab ], NA )
}
GENE.tab["C","LAB"] <- "HLA-CW"

#############################################################
## CREATE FLAT AA TABLE FOR EACH GENE #######################
#############################################################

## Load Amino Acid Reference Data
AA <- list()
for ( g in 1:length(GENE.all) ) {
	gene <- GENE.all[g]
	PathToAA <- paste(PathToRefs,gene,"_prot.txt",sep="")
	if ( !file.exists(PathToAA) ) {
		PathToAA <- paste(PathToRefs,gsub("1","",gene),"_prot.txt",sep="")
		if ( !file.exists(PathToAA) ) {
			next
		}
	}
	# Get File Info
	META <- read.table( pipe(paste("cat",PathToAA,"| grep -n Prot") ),fill=T )
	closeAllConnections()
	SKIP_LINES <- as.numeric(gsub(":","",META[2:nrow(META),1]))
	NUM_LINES <- SKIP_LINES[2:length(SKIP_LINES)] - SKIP_LINES[2:length(SKIP_LINES)-1]
	PROT_POS <- as.numeric(as.character( META[2:nrow(META),3] ))
	## Load Table & Pull Name/Sequences
	AA.1 <- AA.split <- AA.type.1 <- AA.prot.1 <- AA.prot.2 <- AA.frame.1 <- list()
	for ( n in 1:length(SKIP_LINES) ) {
		if ( n < length(SKIP_LINES) ) { AA.1[[n]] <- read.table( PathToAA, sep="\t",header=F, skip=SKIP_LINES[n]+1, nrow=NUM_LINES[n]-3, fill=T, colClasses="character" )[,1] }
		if ( n == length(SKIP_LINES) ) {
			AA.1[[n]] <- read.table( PathToAA, sep="\t",header=F, skip=SKIP_LINES[n]+1, fill=T, colClasses="character" )[,1]
			AA.1[[n]] <- AA.1[[n]][1:(length(AA.1[[n]])-1)]
		}
		AA.split[[n]] <- strsplit( AA.1[[n]], "  " )
		AA.type.1[[n]] <- gsub(paste(" ",gene,"*",sep=""),"",unlist(lapply( AA.split[[n]], function(x) head(x,1) )),fixed=T)
		AA.prot.1[[n]] <- gsub(" ","",unlist(lapply( AA.split[[n]], function(x) x[which(nchar(x[2:length(x)])>0)[1]+1] )) )
		if ( n == 1 ) {
			AA.paste <- AA.prot.1[[n]]
		}else{
			AA.paste <- paste( AA.paste, AA.prot.1[[n]], sep="" )
		}
	}
	## Create Table w/ Type & Each Amino Acid
	MAX.nCHAR <- max( sapply(AA.paste,nchar) )
	START_POS <- PROT_POS[2]-max(sapply(AA.prot.1[[1]],nchar))
	STOP_POS <- START_POS + MAX.nCHAR - 1
	AA.2 <- t(sapply(strsplit(AA.paste,""),"[",1:MAX.nCHAR))
	colnames(AA.2) <- paste("Pos",START_POS:STOP_POS,sep="_")
	rownames(AA.2) <- AA.type.1[[1]]
	AA.3 <- AA.2 ; for ( col in 1:ncol(AA.3) ) { AA.3[which(AA.3[,col]=="-"),col] <- AA.3[1,col] }
	AA[[gene]] <- AA.3
	## Move to next Gene
	print(paste("Done with",gene))
}

## Compile Possibilities from Ambiguous Lab Data
OPTIONS.1 <- sapply(strsplit( gsub("[[:space:]]","",as.character(LAB.2$Code_Translation_1,fixed=T)), "="),"[",2)
OPTIONS.2 <- sapply(strsplit( gsub("[[:space:]]","",as.character(LAB.2$Code_Translation_2,fixed=T)), "="),"[",2)
Sentence <- "The following rare alleles could not be ruled out by the SBT procedure and their presumed frequency is very low: "
LOW_FREQ.1 <- gsub(",","/",gsub(Sentence,"",as.character(LAB.2$Comment_1)))
LOW_FREQ.2 <- gsub(",","/",gsub(Sentence,"",as.character(LAB.2$Comment_2)))
LOW_FREQ.2[grep("are identical",LOW_FREQ.2)] <- NA
LAB.3  <- data.frame( LAB.2[,1:7],OPTIONS.1,LOW_FREQ.1,Allele_2=LAB.2[,"Allele_2"],OPTIONS.2,LOW_FREQ.2 )

#############################################################
## COMPILE TYPES ############################################
#############################################################

N.plat <- length(PLAT.names)
N.samp <- length(SAMP.all)
N.gene <- nrow(GENE.tab)

#### RAW/BEST VALUES ####
 # Make 3-D table for Gene, Type, &
TYPE.r.colnames <- c(PLAT.names,"SOP.alt","LAB.opt","LAB.lik") 
TYPE.r <- array( , c(N.samp,length(TYPE.r.colnames),N.gene) )
rownames(TYPE.r) <- SAMP.all
colnames(TYPE.r) <- TYPE.r.colnames
dimnames(TYPE.r)[[3]] <- GENE.all
for ( s in 1:N.samp ) {
	samp <- SAMP.all[s]
	for ( g in 1:N.gene ) {
		gene <- GENE.all[g]
		## Pull SOAP-HLA Type
		gene_sop <- GENE.tab[gene,"SOP"]
		TEMP_TAB <- SOP.l$GENES.list[[gene_sop]]
		TEMP_TAB <- TEMP_TAB[,grep(samp,colnames(TEMP_TAB))]
		TYPE.sop <- rownames(TEMP_TAB)[which(rowSums(TEMP_TAB)>0)]
		TYPE.sop <- ifelse( length(TYPE.sop)==2, paste(TYPE.sop,collapse="///"), ifelse( length(TYPE.sop)==1,paste(rep(TYPE.sop,2),collapse="///"),"NA///NA" ) )
		TYPE.r[samp,"SOP",gene] <- TYPE.sop
		 # Alternate SOAP-HLA Type
		TEMP_ALT <- SOP.l$CONF.list[[gene_sop]]
		TEMP_ALT <- TEMP_ALT[,grep(samp,colnames(TEMP_ALT))]
		TYPE.sop.alt <- rownames(TEMP_ALT)[apply( TEMP_ALT, 2, function(x) order(x[which(x>0)],decreasing=T)[2] )]
		TYPE.sop.alt <- ifelse( length(TYPE.sop.alt)==2, paste(TYPE.sop.alt,collapse="///"), ifelse( length(TYPE.sop.alt)==1,paste(rep(TYPE.sop.alt,2),collapse="///"),"NA///NA" ) )
		# TYPE.sop.alt <- ifelse( length(TYPE.sop.alt)>1, paste(TYPE.sop.alt,collapse="///"), paste(rep(TYPE.sop.alt,2),collapse="///") )
		TYPE.r[samp,"SOP.alt",gene] <- TYPE.sop.alt
		## Pull CHIP Type
		gene_chp <- GENE.tab[gene,"CHP"]
		if ( !is.na(gene_chp) ) {
			TYPE.chp.1 <- CHP.4.l[samp,paste(gene_chp,"1",sep="")]
			TYPE.chp.2 <- CHP.4.l[samp,paste(gene_chp,"2",sep="")]
			TYPE.chp <- paste(TYPE.chp.1,TYPE.chp.2,sep="///")
			TYPE.r[samp,"CHP",gene] <- TYPE.chp
		}else{
			TYPE.r[samp,"CHP",gene] <- paste(NA,NA,sep="///")
		}
		## Pull SEQ Type
		gene_seq <- GENE.tab[gene,"SEQ"]
		if ( !is.na(gene_seq) ) {
			TYPE.seq.1 <- SEQ.4.l[samp,paste(gene_seq,"1",sep="")]
			TYPE.seq.2 <- SEQ.4.l[samp,paste(gene_seq,"2",sep="")]
			TYPE.seq <- paste(TYPE.seq.1,TYPE.seq.2,sep="///")
			TYPE.r[samp,"SEQ",gene] <- TYPE.seq
		}else{
			TYPE.r[samp,"SEQ",gene] <- paste(NA,NA,sep="///")
		}
		## Pull LAB Type
		gene_lab <- GENE.tab[gene,"LAB"]
		which_row <- which( LAB.2$TestName==gene_lab & LAB.2$ID==samp )
		if ( length(which_row)>0 ) {
			# Types as Stated
			TYPE.lab.1 <- as.character( LAB.2$Allele_1[which_row] )
			TYPE.lab.2 <- as.character( LAB.2$Allele_2[which_row] )
			# Types w/ Ambiguity
			ALT.lab.1 <- as.character( LAB.3$OPTIONS.1[which_row] )
			ALT.lab.2 <- as.character( LAB.3$OPTIONS.2[which_row] )
			# Most Likely
			LOWf.lab.1 <- as.character( LAB.3$LOW_FREQ.1[which_row] )
			LOWf.lab.2 <- as.character( LAB.3$LOW_FREQ.2[which_row] )
			LIK.lab.1 <- setdiff( strsplit(ALT.lab.1,"/")[[1]], strsplit(LOWf.lab.1,"/")[[1]] ) ; if ( length(LIK.lab.1)==0 ) { LIK.lab.1 <- NA }
			LIK.lab.2 <- setdiff( strsplit(ALT.lab.2,"/")[[1]], strsplit(LOWf.lab.2,"/")[[1]] ) ; if ( length(LIK.lab.2)==0 ) { LIK.lab.2 <- NA }
			if ( TYPE.lab.1=="-" ) { TYPE.lab.1 <- TYPE.lab.2 ; ALT.lab.1 <- ALT.lab.2 ; LIK.lab.1 <- LIK.lab.2 }
			if ( TYPE.lab.2=="-" ) { TYPE.lab.2 <- TYPE.lab.1 ; ALT.lab.2 <- ALT.lab.1 ; LIK.lab.2 <- LIK.lab.1 }
			TYPE.lab <- paste(TYPE.lab.1,TYPE.lab.2,sep="///")
			TYPE.r[samp,"LAB",gene] <- TYPE.lab
			ALT.lab <- paste(ALT.lab.1,ALT.lab.2,sep="///")
			TYPE.r[samp,"LAB.opt",gene] <- ALT.lab
			LIK.lab <- paste( paste(LIK.lab.1,collapse="/"), paste(LIK.lab.2,collapse="/"), sep="///")
			TYPE.r[samp,"LAB.lik",gene] <- LIK.lab
			# if ( LIK.lab!="///" ) { TYPE.r[samp,"LAB.lik",gene] <- LIK.lab
			# }else{ TYPE.r[samp,"LAB.lik",gene] <- paste(NA,NA,sep="///") }
		}else{
			TYPE.r[samp,"LAB",gene] <- TYPE.r[samp,"LAB.opt",gene] <- TYPE.r[samp,"LAB.lik",gene] <- paste(NA,NA,sep="///")
		}
	}
	## Spit out Update
	if ( s%%20 == 0 ) {
		print(paste( "Done with",s,"of",N.samp,"Samples" ))
	}
}
for ( g in 1:N.gene ) { TYPE.r[which(TYPE.r[,1,g]==""),1,g] <- paste(NA,NA,sep="///") }

#### FOUR-DIGIT ####
TYPE.4 <- TYPE.r
for ( g in 1:N.gene ) {
	gene <- GENE.all[g]
	for ( col in c("CHP","SEQ") ) {
		# TYPE.4[,col,gene] <- gsub(":","",TYPE.4[,col,gene] )
		# TYPE.4[,col,gene] <- gsub("[A-Z]","",TYPE.4[,col,gene] )
		split.1 <- strsplit( TYPE.4[,col,gene], "///" )
		split.2 <- lapply( split.1, function(x) substr( x,1,4 ) )
		split.3 <- unlist(lapply( split.2, function(x) paste( x,collapse="///" ) ))
		TYPE.4[,col,gene] <- split.3
	}
	for ( col in c("SOP","LAB","SOP.alt") ) {
		split.1 <- strsplit( TYPE.4[,col,gene], "///" )
		split.2 <- lapply( split.1, function(x) unlist(lapply( strsplit(x,":"),function(y)paste(y[1:min(2,length(y))],collapse="") )) )
		split.3 <- unlist(lapply( split.2, function(x) paste( x,collapse="///" ) ))
		TYPE.4[,col,gene] <- split.3
	}
	for ( col in c("LAB.opt","LAB.lik") ) {
		split.1 <- strsplit( TYPE.4[,col,gene], "///" )
		split.1b <- lapply( split.1, function(x) strsplit(x,"/") )
		split.2 <- lapply( split.1b, function(x) unlist(lapply( x,function(y) paste(unlist(lapply(strsplit(y,":"),function(z)paste(z[1:min(2,length(z))],collapse="") )),collapse="/") )) )
		split.3 <- unlist(lapply( split.2, function(x) paste( x,collapse="///" ) ))
		TYPE.4[,col,gene] <- split.3	
	}
}
TYPE.4[SAMP.lab,,"DRB"]

#### TWO-DIGIT ####
TYPE.2 <- TYPE.r
for ( g in 1:N.gene ) {
	gene <- GENE.all[g]
	for ( col in c("CHP","SEQ") ) {
		# TYPE.2[,col,gene] <- gsub(":","",TYPE.2[,col,gene] )
		# TYPE.2[,col,gene] <- gsub("[A-Z]","",TYPE.2[,col,gene] )
		split.1 <- strsplit( TYPE.4[,col,gene], "///" )
		split.2 <- lapply( split.1, function(x) substr( x,1,2 ) )
		split.3 <- unlist(lapply( split.2, function(x) paste( x,collapse="///" ) ))
		TYPE.2[,col,gene] <- split.3
	}
	for ( col in c("SOP","LAB","SOP.alt") ) {
		split.1 <- strsplit( TYPE.2[,col,gene], "///" )
		split.2 <- lapply( split.1, function(x) unlist(lapply( strsplit(x,":"),function(y)y[1] )) )
		split.3 <- unlist(lapply( split.2, function(x) paste( x,collapse="///" ) ))
		TYPE.2[,col,gene] <- split.3
	}
	for ( col in c("LAB.opt","LAB.lik") ) {
		split.1 <- strsplit( TYPE.2[,col,gene], "///" )
		split.1b <- lapply( split.1, function(x) strsplit(x,"/") )
		split.2 <- lapply( split.1b, function(x) unlist(lapply( x,function(y) paste(unlist(lapply(strsplit(y,":"),function(z)z[1] )),collapse="/") )) )
		split.3 <- unlist(lapply( split.2, function(x) paste( x,collapse="///" ) ))
		TYPE.2[,col,gene] <- split.3	
	}
}
TYPE.2[SAMP.lab,,"DRB"]

#############################################################
## SUMMARIZE TYPING BY PATIENT & PLATFORM ###################
#############################################################

## Split Patients to Different Cohorts
PATS.3k1 <- grep("-",rownames(TYPE.r),invert=T)
PATS.3k2 <- grep("-",rownames(TYPE.r))

#############################################################
## Number of Patients Typed by Each Platform (by Gene) ######
N.3k1.na2 <- apply( TYPE.r[PATS.3k1,,], c(2,3), function(x)length(grep("NA///NA",x,invert=T)) )
N.3k2.na2 <- apply( TYPE.r[PATS.3k2,,], c(2,3), function(x)length(grep("NA///NA",x,invert=T)) )

## Plot it
COLS.plat <- c("chartreuse1","orange1","deepskyblue1","firebrick1","mediumpurple1","yellow1","seagreen1")
names(COLS.plat) <- colnames(TYPE.r)
c <- 1:3
 # ART3001
png( paste(PathToPlot,"1-Patients_Typed.3k1.na2.",paste(c,collapse=""),".png",sep=""),height=1400,width=2000,pointsize=36 )
YLIM <- range( N.3k1.na2 ) * c(1,1.2)
barplot( N.3k1.na2[c,], beside=T, col=COLS.plat[c], ylim=YLIM, main="HLA Typing of ART3001 by Gene/Platform",xlab="Gene",ylab="# Patients")
abline( h=seq(0,800,25),lty=3,col="grey50",lwd=1 )
barplot( N.3k1.na2[c,], beside=T, col=COLS.plat[c], add=T )
legend("topright",fill=COLS.plat[c],legend=rownames(N.3k1.na2)[c],ncol=2,bg="white")
dev.off()
 # ART3002
png( paste(PathToPlot,"1-Patients_Typed.3k2.na2.",paste(c,collapse=""),".png",sep=""),height=1400,width=2000,pointsize=36 )
YLIM <- range( N.3k2.na2 ) * c(1,1.2)
barplot( N.3k2.na2[c,], beside=T, col=COLS.plat[c], ylim=YLIM, main="HLA Typing of ART3002 by Gene/Platform",xlab="Gene",ylab="# Patients")
abline( h=seq(0,800,25),lty=3,col="grey50",lwd=1 )
barplot( N.3k2.na2[c,], beside=T, col=COLS.plat[c], add=T )
legend("topright",fill=COLS.plat[c],legend=rownames(N.3k2.na2)[c],ncol=2,bg="white")
dev.off()

## How many have BOTH alleles determined by Comp Method?
N.3k1.na1 <- apply( TYPE.r[PATS.3k1,,], c(2,3), function(x)length(grep("NA",x,invert=T)) )
## Plot it
COLS.plat <- c("chartreuse1","orange","deepskyblue2","firebrick1","mediumpurple2","yellow2","lightseagreen")
c <- 1:3
png( paste(PathToPlot,"1-Patients_Typed.3k1.na1.",paste(c,collapse=""),".png",sep=""),height=1400,width=2000,pointsize=36 )
 # ART3001
YLIM <- range( N.3k1.na2 ) * c(1,1.2)
barplot( N.3k1.na1[c,], beside=T, col=COLS.plat[c], ylim=YLIM, main="HLA Typing of ART3001 by Gene/Platform",xlab="Gene",ylab="# Patients")
abline( h=seq(0,800,25),lty=3,col="grey50",lwd=1 )
barplot( N.3k1.na1[c,], beside=T, col=COLS.plat[c], add=T )
legend("topright",fill=COLS.plat[c],legend=rownames(N.3k1.na1)[c],ncol=2,bg="white")
dev.off()

## Show with LAB
c <- 1:4
 # ART3001
png( paste(PathToPlot,"1-Patients_Typed.3k1.na2.",paste(c,collapse=""),".png",sep=""),height=1400,width=2000,pointsize=36 )
YLIM <- range( N.3k1.na2 ) * c(1,1.2)
barplot( N.3k1.na2[c,], beside=T, col=COLS.plat[c], ylim=YLIM, main="HLA Typing of ART3001 by Gene/Platform",xlab="Gene",ylab="# Patients")
abline( h=seq(0,800,25),lty=3,col="grey50",lwd=1 )
barplot( N.3k1.na2[c,], beside=T, col=COLS.plat[c], add=T )
legend("topright",fill=COLS.plat[c],legend=rownames(N.3k1.na2)[c],ncol=2,bg="white")
dev.off()

#############################################################
## At What Precision Are Methods Determined? ################

Genes <- c("A","B","C","DQB","DRB")
## Better than 4-digit Precision
N.3k1.PrGT4.2 <- apply( TYPE.r[PATS.3k1,,Genes], c(2,3), function(x) length(which( unlist(lapply(strsplit(x,"///"),function(y) length(which(str_count(y,":")>1)) )) ==2)) )
N.3k1.PrGT4.1 <- apply( TYPE.r[PATS.3k1,,Genes], c(2,3), function(x) length(which( unlist(lapply(strsplit(x,"///"),function(y) length(which(str_count(y,":")>1)) )) ==1)) )
N.3k1.PrGT4.0 <- N.3k1.na2[,Genes] - (N.3k1.PrGT4.2+N.3k1.PrGT4.1)
 # Plot it
png( paste(PathToPlot,"2a-Patients_Typed.3k1.PrGT4.png",sep=""),height=1400,width=2000,pointsize=36 )
par(mfrow=c(2,2))
par(mar=c(4,4,3,2))
for ( which_plat in 1:4 ) {
	To_Plot <- rbind( N.3k1.PrGT4.2[which_plat,],N.3k1.PrGT4.1[which_plat,],N.3k1.PrGT4.0[which_plat,] ) ; rownames(To_Plot) <- paste(2:0,"Alleles",sep="_")
	COLS.temp <- colorRampPalette(c(COLS.plat[which_plat],"black"))(10)[c(1,4,7)]
	 # ART3001
	YLIM <- range( N.3k1.na2 ) * c(1,1.2)
	barplot( To_Plot, beside=F, col=COLS.temp, ylim=YLIM, main="HLA Typing of ART3001 by Gene: >4-Digit Precision",xlab="Gene",ylab="# Patients")
	abline( h=seq(0,800,25),lty=3,col="grey50",lwd=1 )
	barplot( To_Plot, beside=F, col=COLS.temp, add=T)
	if ( which_plat==4 ) { legend("topright",fill=COLS.plat[1:4],legend=rownames(N.3k1.na1)[1:4],ncol=2,bg="white",title="Platform") }
	if ( which_plat==3 ) { legend("topleft",fill=COLS.temp,legend=rownames(To_Plot),ncol=3,bg="white",title=">4-digit Precision") }
}
dev.off()
## 4-digit Precision
N.3k1.Pr4.2 <- apply( TYPE.4[PATS.3k1,,Genes], c(2,3), function(x) length(which( unlist(lapply(strsplit(x,"///"),function(y) length(which(nchar(gsub("[A-Z]","",y))>=4)) )) ==2)) )
N.3k1.Pr4.1 <- apply( TYPE.4[PATS.3k1,,Genes], c(2,3), function(x) length(which( unlist(lapply(strsplit(x,"///"),function(y) length(which(nchar(gsub("[A-Z]","",y))>=4)) )) ==1)) )
N.3k1.Pr4.0 <- N.3k1.na2[,Genes] - (N.3k1.Pr4.2+N.3k1.Pr4.1)
 # Plot it
png( paste(PathToPlot,"2b-Patients_Typed.3k1.Pr4.png",sep=""),height=1400,width=2000,pointsize=36 )
par(mfrow=c(2,2))
par(mar=c(4,4,3,2))
for ( which_plat in 1:4 ) {
	To_Plot <- rbind( N.3k1.Pr4.2[which_plat,],N.3k1.Pr4.1[which_plat,],N.3k1.Pr4.0[which_plat,] ) ; rownames(To_Plot) <- paste(2:0,"Alleles",sep="_")
	COLS.temp <- colorRampPalette(c(COLS.plat[which_plat],"black"))(10)[c(1,4,7)]
	YLIM <- range( N.3k1.na2 ) * c(1,1.2)
	barplot( To_Plot, beside=F, col=COLS.temp, ylim=YLIM, main="HLA Typing of ART3001 by Gene: 4-Digit Precision",xlab="Gene",ylab="# Patients")
	abline( h=seq(0,800,25),lty=3,col="grey50",lwd=1 )
	barplot( To_Plot, beside=F, col=COLS.temp, add=T)
	if ( which_plat==4 ) { legend("topright",fill=COLS.plat[1:4],legend=rownames(N.3k1.na1)[1:4],ncol=2,bg="white",title="Platform") }
	if ( which_plat==3 ) { legend("topleft",fill=COLS.temp,legend=rownames(To_Plot),ncol=3,bg="white",title=">4-digit Precision") }
}
dev.off()
## 2-digit Precision (aka, "any")
 # Plot it
png( paste(PathToPlot,"2c-Patients_Typed.3k1.Pr2.png",sep=""),height=1400,width=2000,pointsize=36 )
par(mfrow=c(2,2))
par(mar=c(4,4,3,2))
for ( which_plat in 1:4 ) {
	To_Plot <- rbind( N.3k1.na1[which_plat,Genes],N.3k1.na2[which_plat,Genes]-N.3k1.na1[which_plat,Genes] ) ; rownames(To_Plot) <- paste(2:1,"Alleles",sep="_")
	COLS.temp <- colorRampPalette(c(COLS.plat[which_plat],"black"))(10)[c(1,4)]
	YLIM <- range( N.3k1.na2 ) * c(1,1.2)
	barplot( To_Plot, beside=F, col=COLS.temp, ylim=YLIM, main="HLA Typing of ART3001 by Gene: 2-Digit Precision",xlab="Gene",ylab="# Patients")
	abline( h=seq(0,800,25),lty=3,col="grey50",lwd=1 )
	barplot( To_Plot, beside=F, col=COLS.temp, add=T)
	if ( which_plat==4 ) { legend("topright",fill=COLS.plat[1:4],legend=rownames(N.3k1.na1)[1:4],ncol=2,bg="white",title="Platform") }
	if ( which_plat==3 ) { legend("topleft",fill=COLS.temp,legend=rownames(To_Plot),ncol=2,bg="white",title="2-digit Precision") }
}
dev.off()

#############################################################
## CALCULATE CONCORDANCE ####################################
#############################################################

## Comparisons
 # LAB vs:
   # SOP
   # CHP
   # SEQ
 # SOP vs:
   # CHP
   # SEQ
 # CHP vs:
   # SEQ

## FCT: Compare One Platform to Several Others
Calc_Conc <- function( TYPE.arr, Platform, Prim_Comp, Samples, Genes ) {
	## Determine Number of Shared Alleles per Person per Gene
	TYPE.plat <- TYPE.arr[ Samples,,Genes ]
	TYPE.plat.conc <- TYPE.plat
	for ( g in 1:length(Genes) ) {
		gene <- Genes[g]
		for ( c in 1:length(Prim_Comp) ) {
			comp <- Prim_Comp[c]
			# Split Alleles
			split.plat <- strsplit( TYPE.plat[,Platform,gene], "///" )
			split.comp <- strsplit( TYPE.plat[,comp,gene], "///" )
			for ( s in 1:nrow(TYPE.plat) ) {
				samp <- rownames(TYPE.plat)[s]
				temp.plat <- split.plat[[samp]] ; temp.plat <- temp.plat[which(temp.plat!="NA")]
				temp.comp <- split.comp[[samp]] ; temp.comp <- temp.comp[which(temp.comp!="NA")]
				if ( length(temp.plat)==0 ) { 
					TYPE.plat.conc[samp,comp,gene] <- NA
					next
				}
				if ( length(temp.comp)>0 ) {
					temp_conc.fwd <- length(which( temp.plat == temp.comp ))
					temp_conc.rev <- length(which( temp.plat == rev(temp.comp) ))
					TYPE.plat.conc[samp,comp,gene] <- max( temp_conc.fwd, temp_conc.rev )
				}else{
					TYPE.plat.conc[samp,comp,gene] <- NA
				}
			}
		}
	}

	## Calculate % Concordance for Each Gene & Platform
	TYPE.plat.c.a <- array( ,c( length(Prim_Comp), dim(TYPE.plat)[3] ))
	colnames(TYPE.plat.c.a) <- dimnames(TYPE.plat)[[3]]
	rownames(TYPE.plat.c.a) <- Prim_Comp
	 # ...and determine # of Patients in Intersection of Platforms
	TYPE.plat.c.b <- TYPE.plat.c.a
	 # Number of Patients w/ 0,1,2 correct alleles
	TYPE.plat.c.c <- array( ,c( 3, dim(TYPE.plat)[3], length(Prim_Comp) ))
	rownames(TYPE.plat.c.c) <- 0:2
	colnames(TYPE.plat.c.c) <- dimnames(TYPE.plat)[[3]]
	dimnames(TYPE.plat.c.c)[[3]] <- Prim_Comp
	 # Loop through Genes & Patients
	for ( g in 1:length(Genes) ) {
		gene <- Genes[g]
		TYPE.plat.c.a[,gene] <- apply( TYPE.plat.conc[,Prim_Comp,gene], 2, function(x) mean(as.numeric(x),na.rm=T)/2 )
		TYPE.plat.c.b[,gene] <- apply( TYPE.plat.conc[,Prim_Comp,gene], 2, function(x) length(which(!is.na(x))) )
		TYPE.plat.c.c[,gene,] <- rbind( apply(TYPE.plat.conc[,Prim_Comp,gene],2,function(x) length(which(x==0))), apply(TYPE.plat.conc[,Prim_Comp,gene],2,function(x) length(which(x==1))), apply(TYPE.plat.conc[,Prim_Comp,gene],2,function(x) length(which(x==2))) )
	}
	## Compile Outputs
	COMPILE <- list( Type=TYPE.plat, Conc=TYPE.plat.conc, A=TYPE.plat.c.a, B=TYPE.plat.c.b, C=TYPE.plat.c.c )
	return(COMPILE)
}

## FCT: Compare One Platform to Several Others (Lab-Best Precision)
Calc_Conc.B <- function( TYPE.arr, Platform, Prim_Comp, Samples, Genes ) {
	TYPE.plat <- TYPE.4[ Samples,,Genes ]
	TYPE.plat.conc <- TYPE.plat
	for ( g in 1:length(Genes) ) {
		gene <- Genes[g]
		for ( c in 1:length(Prim_Comp) ) {
			comp <- Prim_Comp[c]
			# Split Alleles
			split.plat <- strsplit( TYPE.plat[,"LAB",gene], "///" )
			split.comp <- strsplit( TYPE.plat[,comp,gene], "///" )
			split.lik <- strsplit( TYPE.plat[,"LAB.lik",gene], "///" )
			for ( s in 1:nrow(TYPE.plat) ) {
				samp <- rownames(TYPE.plat)[s]
				temp.plat <- split.plat[[samp]] ; temp.plat <- temp.plat[which(temp.plat!="NA")]
				temp.comp <- split.comp[[samp]] ; temp.comp <- temp.comp[which(temp.comp!="NA")]
				if ( plat=="LAB" & any(grepl("[A-Z]",temp.plat)) & gene=="DRB" ) {
					temp.plat[grepl("[A-Z]",temp.plat)] <- split.lik[[samp]][grepl("[A-Z]",temp.plat)] ; temp.plat <- temp.plat[which(temp.plat!="NA")]
				}else{ temp.plat <- gsub("[A-Z]","",temp.plat) }
				if ( length(temp.plat)==0 ) { 
					TYPE.plat.conc[samp,comp,gene] <- NA
					next
				}
				if ( length(temp.comp)>0 ) {
					prec <- nchar( gsub("[A-Z]","",temp.plat) )
					temp_conc.fwd.1 <- length(which( temp.plat == c( substr(temp.comp[1],1,prec[1]),substr(temp.comp[2],1,prec[2]) ) ))
					temp_conc.fwd.2 <- length(which( temp.plat == c( substr(temp.comp[1],1,prec[2]),substr(temp.comp[2],1,prec[1]) ) ))
					temp_conc.rev.1 <- length(which( temp.plat == rev(c( substr(temp.comp[1],1,prec[1]),substr(temp.comp[2],1,prec[2]) )) ))
					temp_conc.rev.2 <- length(which( temp.plat == rev(c( substr(temp.comp[1],1,prec[2]),substr(temp.comp[2],1,prec[1]) )) ))
					# temp_conc.rev <- length(which( temp.plat == rev(temp.comp) ))
					TYPE.plat.conc[samp,comp,gene] <- max( temp_conc.fwd.1, temp_conc.fwd.2, temp_conc.rev.1, temp_conc.rev.2  )
				}else{ TYPE.plat.conc[samp,comp,gene] <- NA }
			}
		}
	}

	## Calculate % Concordance for Each Gene & Platform
	TYPE.plat.c.a <- array( ,c( length(Prim_Comp), dim(TYPE.plat)[3] ))
	colnames(TYPE.plat.c.a) <- dimnames(TYPE.plat)[[3]]
	rownames(TYPE.plat.c.a) <- Prim_Comp
	 # ...and determine # of Patients in Intersection of Platforms
	TYPE.plat.c.b <- TYPE.plat.c.a
	 # Number of Patients w/ 0,1,2 correct alleles
	TYPE.plat.c.c <- array( ,c( 3, dim(TYPE.plat)[3], length(Prim_Comp) ))
	rownames(TYPE.plat.c.c) <- 0:2
	colnames(TYPE.plat.c.c) <- dimnames(TYPE.plat)[[3]]
	dimnames(TYPE.plat.c.c)[[3]] <- Prim_Comp
	 # Loop through Genes & Patients
	for ( g in 1:length(Genes) ) {
		gene <- Genes[g]
		TYPE.plat.c.a[,gene] <- apply( TYPE.plat.conc[,Prim_Comp,gene], 2, function(x) mean(as.numeric(x),na.rm=T)/2 )
		TYPE.plat.c.b[,gene] <- apply( TYPE.plat.conc[,Prim_Comp,gene], 2, function(x) length(which(!is.na(x))) )
		TYPE.plat.c.c[,gene,] <- rbind( apply(TYPE.plat.conc[,Prim_Comp,gene],2,function(x) length(which(x==0))), apply(TYPE.plat.conc[,Prim_Comp,gene],2,function(x) length(which(x==1))), apply(TYPE.plat.conc[,Prim_Comp,gene],2,function(x) length(which(x==2))) )
	}
	## Compile Outputs
	COMPILE <- list( Type=TYPE.plat, Conc=TYPE.plat.conc, A=TYPE.plat.c.a, B=TYPE.plat.c.b, C=TYPE.plat.c.c )
	return(COMPILE)
}

###############################################
## TWO-DIGIT ##################################
COMP <- list()
COMP$P2 <- list()

####  LAB vs ... ####
Platform <- "LAB"
Prim_Comp <- c("SOP","CHP","SEQ","SOP.alt")
Samples <- SAMP.lab
COMP$P2$lab <- Calc_Conc( TYPE.2, Platform, Prim_Comp, Samples, Genes )
####  SOP vs ... ####
Platform <- "SOP"
Prim_Comp <- c("CHP","SEQ","LAB","LAB.lik")
Samples <- SAMP.sop
COMP$P2$sop <- Calc_Conc( TYPE.2, Platform, Prim_Comp, Samples, Genes )
####  SEQ vs ... ####
Platform <- "SEQ"
Prim_Comp <- c("CHP","SOP.alt","LAB","LAB.lik") # PLAT.names[1:3]
Samples <- SAMP.seq
COMP$P2$seq <- Calc_Conc( TYPE.2, Platform, Prim_Comp, Samples, Genes )

###############################################
## FOUR-DIGIT #################################
COMP$P4 <- list()

####  LAB vs ... ####
Platform <- "LAB"
Prim_Comp <- c("SOP","CHP","SEQ","SOP.alt")
Samples <- SAMP.lab
COMP$P4$lab <- Calc_Conc( TYPE.4, Platform, Prim_Comp, Samples, Genes )
####  SOP vs ... ####
Platform <- "SOP"
Prim_Comp <- c("CHP","SEQ","LAB","LAB.lik")
Samples <- SAMP.sop
COMP$P4$sop <- Calc_Conc( TYPE.4, Platform, Prim_Comp, Samples, Genes )
####  SEQ vs ... ####
Platform <- "SEQ"
Prim_Comp <- c("CHP","SOP.alt","LAB","LAB.lik") # PLAT.names[1:3]
Samples <- SAMP.seq
COMP$P4$seq <- Calc_Conc( TYPE.4, Platform, Prim_Comp, Samples, Genes )

# TEMP <- Reduce( rbind, lapply(COMP$P4,function(x)x$A) ) ; data.frame( COMP=rownames(TEMP), PLAT=rep(names(COMP$P4),unlist(lapply(COMP$P4,function(x)nrow(x$A))) ), TEMP )
# lapply( Genes, function(x)COMP$P2$lab$Type[which(COMP$P2$lab$Conc[,"SOP",x]<2),,x] )
# lapply( Genes, function(x)COMP$P4$lab$Type[which(COMP$P4$lab$Conc[,"SOP",x]<2),,x] )
# COMP$P4$lab$Conc[,,"DRB"]
# COMP$P4$sop$Conc[,,"DRB"]
# COMP$P4$sop$Type[which(COMP$P4$sop$Conc[,"LAB","DRB"]<2),,"DRB"]

###############################################
## LAB-BEST ###################################

COMP$PB <- list()

####  LAB vs ... ####
Platform <- "LAB"
Prim_Comp <- c("SOP","CHP","SEQ","SOP.alt")
Samples <- SAMP.lab
COMP$PB$lab <- Calc_Conc.B( TYPE.4, Platform, Prim_Comp, Samples, Genes )
####  SOP vs ... ####
Platform <- "SOP"
Prim_Comp <- c("CHP","SEQ","LAB","LAB.lik")
Samples <- SAMP.sop
COMP$PB$sop <- Calc_Conc.B( TYPE.4, Platform, Prim_Comp, Samples, Genes )
####  SEQ vs ... ####
Platform <- "SEQ"
Prim_Comp <- c("CHP","SOP.alt","LAB","LAB.lik") # PLAT.names[1:3]
Samples <- SAMP.seq
COMP$PB$seq <- Calc_Conc.B( TYPE.4, Platform, Prim_Comp, Samples, Genes )

#############################################################
## PLOT CONCORDANCE #########################################
#############################################################

###############################################
## Non-LAB vs ... #############################
 # 2-digit
COMP.nlab.2 <- rbind( COMP$P2$sop$A[1:2,], COMP$P2$seq$A[1,] )
COMP.nlab.2.b <- rbind( COMP$P2$sop$B[1:2,], COMP$P2$seq$B[1,] )
rownames(COMP.nlab.2) <- rownames(COMP.nlab.2.b) <- c("SOPvCHP","SOPvSEQ","SEQvCHP")
 # 4-digit
COMP.nlab.4 <- rbind( COMP$P4$sop$A[1:2,], COMP$P4$seq$A[1,] )
COMP.nlab.4.b <- rbind( COMP$P4$sop$B[1:2,], COMP$P4$seq$B[1,] )
rownames(COMP.nlab.4) <- rownames(COMP.nlab.4.b) <- c("SOPvCHP","SOPvSEQ","SEQvCHP")

## Plot It
YLIM <- c(0,max(COMP.nlab.2)) * c(1,1.2)
png( paste(PathToPlot,"3-Compare_NonLab.png",sep=""),height=1200,width=2000,pointsize=36 )
par(mfrow=c(1,2))
 # 2-digit
TEMP <- barplot( COMP.nlab.2, beside=T, col=COLS.plat[sapply(strsplit(rownames(COMP.nlab.2),"v"),"[",1)],ylim=YLIM,yaxt="n",main="Concordance at 2-Digit Precision",xlab="Gene",ylab="Fraction of Patients" )
axis(2,at=seq(0,1,.1) ) ; abline(h=seq(0,1,.1),lty=3,col="grey50",lwd=1)
barplot( COMP.nlab.2, beside=T, col=COLS.plat[sapply(strsplit(rownames(COMP.nlab.2),"v"),"[",2)],density=30,yaxt="n",add=T )
legend( "topright",legend=rownames(COMP.nlab.2),fill=COLS.plat[sapply(strsplit(rownames(COMP.nlab.2),"v"),"[",1)], title="Comparison",ncol=3,cex=.7)
legend( "topright",legend=rownames(COMP.nlab.2),fill=COLS.plat[sapply(strsplit(rownames(COMP.nlab.2),"v"),"[",2)],density=30, title="Comparison",bg=NA,ncol=3,cex=.7)
text( TEMP, .2, srt=90,label=paste("n=",round(c(COMP.nlab.2)*2*c(COMP.nlab.2.b)),"/",2*c(COMP.nlab.2.b),sep="") )
 # 4-digit
TEMP <- barplot( COMP.nlab.4, beside=T, col=COLS.plat[sapply(strsplit(rownames(COMP.nlab.4),"v"),"[",1)],ylim=YLIM,yaxt="n",main="Concordance at 4-Digit Precision",xlab="Gene",ylab="Fraction of Patients" )
axis(2,at=seq(0,1,.1) ) ; abline(h=seq(0,1,.1),lty=3,col="grey50",lwd=1)
barplot( COMP.nlab.4, beside=T, col=COLS.plat[sapply(strsplit(rownames(COMP.nlab.4),"v"),"[",2)],density=30,yaxt="n",add=T )
legend( "topright",legend=rownames(COMP.nlab.4),fill=COLS.plat[sapply(strsplit(rownames(COMP.nlab.4),"v"),"[",1)], title="Comparison",ncol=3,cex=.7)
legend( "topright",legend=rownames(COMP.nlab.4),fill=COLS.plat[sapply(strsplit(rownames(COMP.nlab.4),"v"),"[",2)],density=30, title="Comparison",bg=NA,ncol=3,cex=.7)
text( TEMP, .2, srt=90,label=paste("n=",round(c(COMP.nlab.4)*2*c(COMP.nlab.4.b)),"/",2*c(COMP.nlab.4.b),sep="") )
dev.off()

###############################################
## LAB vs ... #################################
 # 2-digit
COMP.lab.2 <- COMP$P2$lab$A[1:3,]
COMP.lab.2.b <- COMP$P2$lab$B[1:3,]
rownames(COMP.lab.2) <- rownames(COMP.lab.2.b) <- paste("LABv",rownames(COMP.lab.2),sep="")
 # 4-digit
COMP.lab.4 <- COMP$P4$lab$A[1:3,]
COMP.lab.4.b <- COMP$P4$lab$B[1:3,]
rownames(COMP.lab.4) <- rownames(COMP.lab.4.b) <- paste("LABv",rownames(COMP.lab.4),sep="")

## Plot It
YLIM <- c(0,max(COMP.lab.2)) * c(1,1.2)
png( paste(PathToPlot,"3-Compare_Lab.png",sep=""),height=1200,width=2000,pointsize=36 )
par(mfrow=c(1,2))
 # 2-digit
TEMP <- barplot( COMP.lab.2, beside=T, col=COLS.plat[sapply(strsplit(rownames(COMP.lab.2),"v"),"[",1)],ylim=YLIM,yaxt="n",main="Concordance at 2-Digit Precision",xlab="Gene",ylab="Fraction of Patients" )
axis(2,at=seq(0,1,.1) ) ; abline(h=seq(0,1,.1),lty=3,col="grey50",lwd=1)
barplot( COMP.lab.2, beside=T, col=COLS.plat[sapply(strsplit(rownames(COMP.lab.2),"v"),"[",2)],density=30,yaxt="n",add=T )
legend( "topright",legend=rownames(COMP.lab.2),fill=COLS.plat[sapply(strsplit(rownames(COMP.lab.2),"v"),"[",1)], title="Comparison",ncol=3,cex=.7)
legend( "topright",legend=rownames(COMP.lab.2),fill=COLS.plat[sapply(strsplit(rownames(COMP.lab.2),"v"),"[",2)],density=30, title="Comparison",bg=NA,ncol=3,cex=.7)
text( TEMP, .2, srt=90,label=paste("n=",round(c(COMP.lab.2)*2*c(COMP.lab.2.b)),"/",2*c(COMP.lab.2.b),sep="") )
 # 4-digit
TEMP <- barplot( COMP.lab.4, beside=T, col=COLS.plat[sapply(strsplit(rownames(COMP.lab.4),"v"),"[",1)],ylim=YLIM,yaxt="n",main="Concordance at 4-Digit Precision",xlab="Gene",ylab="Fraction of Patients" )
axis(2,at=seq(0,1,.1) ) ; abline(h=seq(0,1,.1),lty=3,col="grey50",lwd=1)
barplot( COMP.lab.4, beside=T, col=COLS.plat[sapply(strsplit(rownames(COMP.lab.4),"v"),"[",2)],density=30,yaxt="n",add=T )
legend( "topright",legend=rownames(COMP.lab.4),fill=COLS.plat[sapply(strsplit(rownames(COMP.lab.4),"v"),"[",1)], title="Comparison",ncol=3,cex=.7)
legend( "topright",legend=rownames(COMP.lab.4),fill=COLS.plat[sapply(strsplit(rownames(COMP.lab.4),"v"),"[",2)],density=30, title="Comparison",bg=NA,ncol=3,cex=.7)
text( TEMP, .2, srt=90,label=paste("n=",round(c(COMP.lab.4)*2*c(COMP.lab.4.b)),"/",2*c(COMP.lab.4.b),sep="") )
dev.off()

 # 4-digit (w/ "Likely LAB Haplotypes")
COMP.lab.4 <- COMP$P4$lab$A[1:3,]
COMP.lab.4["SOP","DRB"] <- COMP.lab.4["SOP","DRB"] + COMP$P4$sop$A["LAB.lik","DRB"]
COMP.lab.4.b <- COMP$P4$lab$B[1:3,]
rownames(COMP.lab.4) <- rownames(COMP.lab.4.b) <- paste("LABv",rownames(COMP.lab.4),sep="")
## Plot It
YLIM <- c(0,max(COMP.lab.2)) * c(1,1.2)
png( paste(PathToPlot,"3-Compare_Lab.lik.png",sep=""),height=1200,width=1000,pointsize=36 )
 # 4-digit
TEMP <- barplot( COMP.lab.4, beside=T, col=COLS.plat[sapply(strsplit(rownames(COMP.lab.4),"v"),"[",1)],ylim=YLIM,yaxt="n",main="Concordance at 4-Digit Precision",xlab="Gene",ylab="Fraction of Patients" )
axis(2,at=seq(0,1,.1) ) ; abline(h=seq(0,1,.1),lty=3,col="grey50",lwd=1)
barplot( COMP.lab.4, beside=T, col=COLS.plat[sapply(strsplit(rownames(COMP.lab.4),"v"),"[",2)],density=30,yaxt="n",add=T )
legend( "topright",legend=rownames(COMP.lab.4),fill=COLS.plat[sapply(strsplit(rownames(COMP.lab.4),"v"),"[",1)], title="Comparison",ncol=3,cex=.7)
legend( "topright",legend=rownames(COMP.lab.4),fill=COLS.plat[sapply(strsplit(rownames(COMP.lab.4),"v"),"[",2)],density=30, title="Comparison",bg=NA,ncol=3,cex=.7)
text( TEMP, .2, srt=90,label=paste("n=",round(c(COMP.lab.4)*2*c(COMP.lab.4.b)),"/",2*c(COMP.lab.4.b),sep="") )
dev.off()

 # Best-Lab (w/ "Likely LAB Haplotypes")
COMP.lab.B <- COMP$PB$lab$A[1:3,]
# COMP.lab.B["SOP","DRB"] <- COMP.lab.B["SOP","DRB"] + COMP$PB$sop$A["LAB.lik","DRB"]
COMP.lab.B.b <- COMP$PB$lab$B[1:3,]
rownames(COMP.lab.B) <- rownames(COMP.lab.B.b) <- paste("LABv",rownames(COMP.lab.B),sep="")
## Plot It
YLIM <- c(0,max(COMP.lab.B)) * c(1,1.2)
png( paste(PathToPlot,"3-Compare_Lab.Best.png",sep=""),height=1200,width=1000,pointsize=36 )
 # B-digit
TEMP <- barplot( COMP.lab.B, beside=T, col=COLS.plat[sapply(strsplit(rownames(COMP.lab.B),"v"),"[",1)],ylim=YLIM,yaxt="n",main="Concordance at Best Lab Precision",xlab="Gene",ylab="Fraction of Patients" )
axis(2,at=seq(0,1,.1) ) ; abline(h=seq(0,1,.1),lty=3,col="grey50",lwd=1)
barplot( COMP.lab.B, beside=T, col=COLS.plat[sapply(strsplit(rownames(COMP.lab.B),"v"),"[",2)],density=30,yaxt="n",add=T )
legend( "topright",legend=rownames(COMP.lab.B),fill=COLS.plat[sapply(strsplit(rownames(COMP.lab.B),"v"),"[",1)], title="Comparison",ncol=3,cex=.7)
legend( "topright",legend=rownames(COMP.lab.B),fill=COLS.plat[sapply(strsplit(rownames(COMP.lab.B),"v"),"[",2)],density=30, title="Comparison",bg=NA,ncol=3,cex=.7)
text( TEMP, .2, srt=90,label=paste("n=",round(c(COMP.lab.B)*2*c(COMP.lab.B.b)),"/",2*c(COMP.lab.B.b),sep="") )
dev.off()



# ###############################################
# ## LAB vs ... #################################

# ## CONCORDANCE By Gene ##
# N.Gene.c <- length(Genes) # ncol(Which_Dat$A)
# # COLS <- c("tomato2","slateblue3","chartreuse2")
# # COLS <- c("chocolate2","slateblue3","dodgerblue1")
# # COLS <- c("tomato2","slateblue3","dodgerblue1")
# XLIM <- c(1,N.Gene.c)
# YLIM <- c(0,1)
# png( paste(PathToPlot,"1-Lab_Concordance.png",sep=""),height=800,width=2400,pointsize=36 )
# par(mfrow=c(1,3))
#  # 2/4 Digit Precision
# plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,xlab="Gene",ylab="% Concordant",main="Concordance w/ Lab HLA Types by Platform",xaxt="n")
# axis( 1,at=1:N.Gene.c,label=colnames(TYPE.B.lab.c.a),las=1 )
# abline( h=seq(0,1,.1),lty=3,col="grey50",lwd=1 )
# for ( i in 1:nrow(TYPE.B.lab.c.a) ) {
# 	points( 1:N.Gene.c, TYPE.2.lab.c.a[i,], col=COLS[i],type="o",lty=2,lwd=3,pch=1 )
# 	points( 1:N.Gene.c, TYPE.4.lab.c.a[i,], col=COLS[i],type="o",lty=4,lwd=3,pch=2 )
# }
# legend( "bottomright", lty=c(1,2,3),pch=c(19,1,2),legend=c("Lab Best","2-Digit","4-Digit"),title="Precision",bg="white",lwd=3 )
#  # Lab Best Precision
# # YLIM <- c(0.5,1)
# plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,xlab="Gene",ylab="% Concordant",main="Concordance w/ Lab HLA Types by Platform",xaxt="n")
# axis( 1,at=1:N.Gene.c,label=colnames(TYPE.B.lab.c.a),las=1 )
# abline( h=seq(0,1,.1),lty=3,col="grey50",lwd=1 )
# for ( i in 1:nrow(TYPE.B.lab.c.a) ) {
# 	points( 1:N.Gene.c, TYPE.B.lab.c.a[i,], col=COLS[i],type="o",lty=1,lwd=3,pch=19 )
# }
# legend( "bottomleft", fill=COLS,legend=rownames(TYPE.2.lab.c.a),title="Platform",bg="white" )
# ## Number of Patients By Gene
# YLIM <- c(0,20)
#  # 2/4 Digit Precision
# barplot( TYPE.B.lab.c.b, beside=T, col=COLS )
# abline( h=seq(0,20,5),lty=3,col="grey50",lwd=1 )
# barplot( TYPE.B.lab.c.b, beside=T, col=COLS, add=T,main="Number of Patients Typed by Platform",xlab="Gene",ylab="# Patients") # ,legend=T,args.legend=c(x="topright",bg="white") )
# dev.off()

# ## TABLE OF CONCORDANCE By Gene ##
# COLS.tab <- sapply( COLS, function(x) colorRampPalette(c("white",x,"black"))(6)[2:4] )
# XLIM <- c(0,20)
# YLIM <- c(0,1)
# png( paste(PathToPlot,"1-Lab_Concordance_Table.png",sep=""),height=800,width=2400,pointsize=36 )
# par(mfrow=c(1,3))
#  # Best Lab Digit Precision
# barplot( prop.table(TYPE.B.lab.c.c,c(2,3))[,,"SOP"], space=c(0,rep(3,4)), beside=F, col=COLS.tab[,1], ylim=YLIM,xaxt="n",xlim=XLIM,xlab="Gene",main="Concordant Alleles per Person: Best Lab Precision",ylab="Fraction Patients w/ 0,1,2 Concordant Alleles" )
# barplot( prop.table(TYPE.B.lab.c.c,c(2,3))[,,"CHP"], space=c(1,rep(3,4)), beside=F, col=COLS.tab[,2], add=T,xaxt="n" )
# barplot( prop.table(TYPE.B.lab.c.c,c(2,3))[,,"SEQ"], space=c(2,rep(3,4)), beside=F, col=COLS.tab[,3], add=T,xaxt="n" )
# axis( 1, at=seq(1.5,20,4), labels=colnames(TYPE.B.lab.c.c) )
#  # 2-Digit Precision
# barplot( prop.table(TYPE.2.lab.c.c,c(2,3))[,,"SOP"], space=c(0,rep(3,4)), beside=F, col=COLS.tab[,1], ylim=YLIM,xaxt="n",xlim=XLIM,xlab="Gene",main="Concordant Alleles per Person: 2-Digit Precision",ylab="Fraction Patients w/ 0,1,2 Concordant Alleles" )
# barplot( prop.table(TYPE.2.lab.c.c,c(2,3))[,,"CHP"], space=c(1,rep(3,4)), beside=F, col=COLS.tab[,2], add=T,xaxt="n" )
# barplot( prop.table(TYPE.2.lab.c.c,c(2,3))[,,"SEQ"], space=c(2,rep(3,4)), beside=F, col=COLS.tab[,3], add=T,xaxt="n" )
# axis( 1, at=seq(1.5,20,4), labels=colnames(TYPE.2.lab.c.c) )
# legend( "topright",title="Conc Alleles",legend=2:0,fill=paste("grey",c(20,50,80),sep="") )
#  # 4-Digit Precision
# barplot( prop.table(TYPE.4.lab.c.c,c(2,3))[,,"SOP"], space=c(0,rep(3,4)), beside=F, col=COLS.tab[,1], ylim=YLIM,xaxt="n",xlim=XLIM,xlab="Gene",main="Concordant Alleles per Person: 4-Digit Precision",ylab="Fraction Patients w/ 0,1,2 Concordant Alleles" )
# barplot( prop.table(TYPE.4.lab.c.c,c(2,3))[,,"CHP"], space=c(1,rep(3,4)), beside=F, col=COLS.tab[,2], add=T,xaxt="n" )
# barplot( prop.table(TYPE.4.lab.c.c,c(2,3))[,,"SEQ"], space=c(2,rep(3,4)), beside=F, col=COLS.tab[,3], add=T,xaxt="n" )
# axis( 1, at=seq(1.5,20,4), labels=colnames(TYPE.4.lab.c.c) )
# dev.off()

#############################################################
## WRITE RESULT TABLES ######################################
#############################################################

## Table of Alleles (LabBest)
TempPath <- paste(PathToPlot,"Types_LabBest.xlsx",sep="")
write.xlsx( TYPE.B.lab[,,1], TempPath, append=F, sheetName=dimnames(TYPE.B.lab)[[3]][1], col.names=T,row.names=T, showNA=T )
for ( g in 2:5 ) {
	gene <- dimnames(TYPE.B.lab)[[3]][g]
	write.xlsx( TYPE.B.lab[,,gene], TempPath, append=T, sheetName=gene, col.names=T,row.names=T, showNA=T )
}
## Table of # Concordant Alleles
TempPath <- paste(PathToPlot,"Concordance_LabBest.xlsx",sep="")
write.xlsx( TYPE.B.lab.conc[,,1], TempPath, append=F, sheetName=dimnames(TYPE.B.lab.conc)[[3]][1], col.names=T,row.names=T, showNA=T )
for ( g in 2:5 ) {
	gene <- dimnames(TYPE.B.lab.conc)[[3]][g]
	write.xlsx( TYPE.B.lab.conc[,,gene], TempPath, append=T, sheetName=gene, col.names=T,row.names=T, showNA=T )
}
## Summary Statistics for Concordant Alleles
TempPath <- paste(PathToPlot,"SummaryTables_LabBest.xlsx",sep="")
write.xlsx( TYPE.B.lab.c.a, TempPath, append=F, sheetName="Percent_Conc", col.names=T,row.names=T, showNA=T )
write.xlsx( TYPE.B.lab.c.c[,,1], TempPath, append=T, sheetName="Correct_Alleles_Num_SOP", col.names=T,row.names=T, showNA=T )
write.xlsx( TYPE.B.lab.c.c[,,2], TempPath, append=T, sheetName="Correct_Alleles_Num_CHP", col.names=T,row.names=T, showNA=T )
write.xlsx( TYPE.B.lab.c.c[,,3], TempPath, append=T, sheetName="Correct_Alleles_Num_SEQ", col.names=T,row.names=T, showNA=T )






















#############################################################
## END OF DOC ###############################################
#############################################################
