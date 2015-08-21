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

## Set Date
DATE <- gsub("-","",Sys.Date())

## Set Paths to Data Sets & Save Locations (TSCC)
PathToSOAP <- "/Users/kstandis/Data/Burn/Data/HLA/20150512_SOAP_HLA_Types/20150511_HLA_Types.Rdata"
PathToS2H <- "/Users/kstandis/Data/Burn/Data/HLA/20140924_SNP2HLA_Types/"
PathToLAB <- "/Users/kstandis/Data/Burn/Data/HLA/20150522_Lab_HLA/LabCorp HLA typing data_ART3001.xls"
PathToPheno <- "/Users/kstandis/Data/Burn/Data/Phenos/Full_Tables/20150520_Full_Table.txt"
PathToPlot <- paste("/Users/kstandis/Data/Burn/Plots/",DATE,"_HLA/",sep="")
dir.create( PathToPlot )

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

#############################################################
## ORGANIZE DATA ############################################
#############################################################

## Set Platform Names
PLAT.names <- c("SOP","CHP","SEQ","LAB")

#### SAMPLE LISTS ####

## Get Sample List for each Platform
SAMP.all <- as.character( FT$ID )
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
## COMPILE TYPES ############################################
#############################################################

N.plat <- length(PLAT.names)
N.samp <- length(SAMP.all)
N.gene <- nrow(GENE.tab)

#### RAW/BEST VALUES ####
 # Make 3-D table for Gene, Type, & 
TYPE.r <- array( , c(N.samp,N.plat,N.gene) )
rownames(TYPE.r) <- SAMP.all
colnames(TYPE.r) <- PLAT.names
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
		TYPE.sop <- ifelse( length(TYPE.sop)>1, paste(TYPE.sop,collapse="///"), paste(rep(TYPE.sop,2),collapse="///") )
		TYPE.r[samp,"SOP",gene] <- TYPE.sop
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
			TYPE.lab.1 <- LAB.2$Allele_1[which_row]
			TYPE.lab.2 <- LAB.2$Allele_2[which_row]
			if ( TYPE.lab.1=="-" ) { TYPE.lab.1 <- TYPE.lab.2 }
			if ( TYPE.lab.2=="-" ) { TYPE.lab.2 <- TYPE.lab.1 }
			TYPE.lab <- paste(TYPE.lab.1,TYPE.lab.2,sep="///")
			TYPE.r[samp,"LAB",gene] <- TYPE.lab
		}else{
			TYPE.r[samp,"LAB",gene] <- paste(NA,NA,sep="///")
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
	TYPE.4[,,gene] <- gsub(":","",TYPE.4[,,gene] )
	TYPE.4[,,gene] <- gsub("[A-Z]","",TYPE.4[,,gene] )
	for ( p in 1:N.plat ) {
		plat <- PLAT.names[p]
		split.1 <- strsplit( TYPE.4[,plat,gene], "///" )
		split.2 <- lapply( split.1, function(x) substr( x,1,4 ) )
		split.3 <- unlist(lapply( split.2, function(x) paste( x,collapse="///" ) ))
		TYPE.4[,plat,gene] <- split.3
	}
}

#### TWO-DIGIT ####
TYPE.2 <- TYPE.r
for ( g in 1:N.gene ) {
	gene <- GENE.all[g]
	TYPE.2[,,gene] <- gsub(":","",TYPE.2[,,gene] )
	TYPE.2[,,gene] <- gsub("[A-Z]","",TYPE.2[,,gene] )
	for ( p in 1:N.plat ) {
		plat <- PLAT.names[p]
		split.1 <- strsplit( TYPE.2[,plat,gene], "///" )
		split.2 <- lapply( split.1, function(x) substr( x,1,2 ) )
		split.3 <- unlist(lapply( split.2, function(x) paste( x,collapse="///" ) ))
		TYPE.2[,plat,gene] <- split.3
	}
}

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

###############################################
## LAB vs ... #################################

#### TWO-DIGIT ####

## Determine Number of Shared Alleles per Person per Gene
 # For each Comparison Platform
TYPE.2.lab.genes <- GENE.all[which(!is.na(GENE.tab[,"LAB"]))]
TYPE.2.lab <- TYPE.2[ SAMP.lab,,TYPE.2.lab.genes ]
To_Compare <- PLAT.names[1:3]
TYPE.2.lab.conc <- TYPE.2.lab
for ( g in 1:length(TYPE.2.lab.genes) ) {
	gene <- TYPE.2.lab.genes[g]
	for ( c in 1:length(To_Compare) ) {
		comp <- To_Compare[c]
		# Split Alleles
		split.lab <- strsplit( TYPE.2.lab[,"LAB",gene], "///" )
		split.comp <- strsplit( TYPE.2.lab[,comp,gene], "///" )
		for ( s in 1:nrow(TYPE.2.lab) ) {
			samp <- rownames(TYPE.2.lab)[s]
			if ( length(split.comp[[samp]])>0 ) {
				temp_conc.fwd <- length(which( split.lab[[samp]] == split.comp[[samp]] ))
				temp_conc.rev <- length(which( split.lab[[samp]] == rev(split.comp[[samp]]) ))
				TYPE.2.lab.conc[samp,comp,gene] <- max( temp_conc.fwd, temp_conc.rev )
			}else{
				TYPE.2.lab.conc[samp,comp,gene] <- NA
			}
		}
	}
}

## Calculate % Concordance for Each Gene & Platform
TYPE.2.lab.c.a <- array( ,c( length(To_Compare), dim(TYPE.2.lab)[3] ))
colnames(TYPE.2.lab.c.a) <- dimnames(TYPE.2.lab)[[3]]
rownames(TYPE.2.lab.c.a) <- To_Compare
 # ...and determine # of Patients in Intersection of Platforms
TYPE.2.lab.c.b <- array( ,c( length(To_Compare), dim(TYPE.2.lab)[3] ))
colnames(TYPE.2.lab.c.b) <- dimnames(TYPE.2.lab)[[3]]
rownames(TYPE.2.lab.c.b) <- To_Compare
 # Number of Patients w/ 0,1,2 correct alleles
TYPE.2.lab.c.c <- array( ,c( 3, dim(TYPE.2.lab)[3], length(To_Compare) ))
rownames(TYPE.2.lab.c.c) <- 0:2
colnames(TYPE.2.lab.c.c) <- dimnames(TYPE.2.lab)[[3]]
dimnames(TYPE.2.lab.c.c)[[3]] <- To_Compare
 # Loop through Genes & Patients
for ( g in 1:length(TYPE.2.lab.genes) ) {
	gene <- TYPE.2.lab.genes[g]
	TYPE.2.lab.c.a[,gene] <- apply( TYPE.2.lab.conc[,To_Compare,gene], 2, function(x) mean(as.numeric(x),na.rm=T)/2 )
	TYPE.2.lab.c.b[,gene] <- apply( TYPE.2.lab.conc[,To_Compare,gene], 2, function(x) length(which(!is.na(x))) )
	TYPE.2.lab.c.c[,gene,] <- rbind( apply(TYPE.2.lab.conc[,To_Compare,gene],2,function(x) length(which(x==0))), apply(TYPE.2.lab.conc[,To_Compare,gene],2,function(x) length(which(x==1))), apply(TYPE.2.lab.conc[,To_Compare,gene],2,function(x) length(which(x==2))) )
}

#### FOUR-DIGIT ####

## Determine Number of Shared Alleles per Person per Gene
 # For each Comparison Platform
TYPE.4.lab.genes <- GENE.all[which(!is.na(GENE.tab[,"LAB"]))]
TYPE.4.lab <- TYPE.4[ SAMP.lab,,TYPE.4.lab.genes ]
To_Compare <- PLAT.names[1:3]
TYPE.4.lab.conc <- TYPE.4.lab
for ( g in 1:length(TYPE.4.lab.genes) ) {
	gene <- TYPE.4.lab.genes[g]
	for ( c in 1:length(To_Compare) ) {
		comp <- To_Compare[c]
		# Split Alleles
		split.lab <- strsplit( TYPE.4.lab[,"LAB",gene], "///" )
		split.comp <- strsplit( TYPE.4.lab[,comp,gene], "///" )
		for ( s in 1:nrow(TYPE.4.lab) ) {
			samp <- rownames(TYPE.4.lab)[s]
			if ( length(split.comp[[samp]])>0 ) {
				temp_conc.fwd <- length(which( split.lab[[samp]] == split.comp[[samp]] ))
				temp_conc.rev <- length(which( split.lab[[samp]] == rev(split.comp[[samp]]) ))
				TYPE.4.lab.conc[samp,comp,gene] <- max( temp_conc.fwd, temp_conc.rev )
			}else{
				TYPE.4.lab.conc[samp,comp,gene] <- NA
			}
		}
	}
}

## Calculate % Concordance for Each Gene & Platform
TYPE.4.lab.c.a <- array( ,c( length(To_Compare), dim(TYPE.4.lab)[3] ))
colnames(TYPE.4.lab.c.a) <- dimnames(TYPE.4.lab)[[3]]
rownames(TYPE.4.lab.c.a) <- To_Compare
 # ...and determine # of Patients in Intersection of Platforms
TYPE.4.lab.c.b <- array( ,c( length(To_Compare), dim(TYPE.4.lab)[3] ))
colnames(TYPE.4.lab.c.b) <- dimnames(TYPE.4.lab)[[3]]
rownames(TYPE.4.lab.c.b) <- To_Compare
 # Number of Patients w/ 0,1,2 correct alleles
TYPE.4.lab.c.c <- array( ,c( 3, dim(TYPE.4.lab)[3], length(To_Compare) ))
rownames(TYPE.4.lab.c.c) <- 0:2
colnames(TYPE.4.lab.c.c) <- dimnames(TYPE.4.lab)[[3]]
dimnames(TYPE.4.lab.c.c)[[3]] <- To_Compare
 # Loop through Genes & Patients
for ( g in 1:length(TYPE.4.lab.genes) ) {
	gene <- TYPE.4.lab.genes[g]
	TYPE.4.lab.c.a[,gene] <- apply( TYPE.4.lab.conc[,To_Compare,gene], 2, function(x) mean(as.numeric(x),na.rm=T)/2 )
	TYPE.4.lab.c.b[,gene] <- apply( TYPE.4.lab.conc[,To_Compare,gene], 2, function(x) length(which(!is.na(x))) )
	TYPE.4.lab.c.c[,gene,] <- rbind( apply(TYPE.4.lab.conc[,To_Compare,gene],2,function(x) length(which(x==0))), apply(TYPE.4.lab.conc[,To_Compare,gene],2,function(x) length(which(x==1))), apply(TYPE.4.lab.conc[,To_Compare,gene],2,function(x) length(which(x==2))) )
}

#### LAB-BEST PRECISION ####

## Determine Number of Shared Alleles per Person per Gene
 # For each Comparison Platform
TYPE.B.lab.genes <- GENE.all[which(!is.na(GENE.tab[,"LAB"]))]
TYPE.B.lab <- TYPE.4[ SAMP.lab,,TYPE.B.lab.genes ]
To_Compare <- PLAT.names[1:3]
TYPE.B.lab.conc <- TYPE.B.lab
for ( g in 1:length(TYPE.B.lab.genes) ) {
	gene <- TYPE.B.lab.genes[g]
	for ( c in 1:length(To_Compare) ) {
		comp <- To_Compare[c]
		# Split Alleles
		split.lab <- strsplit( TYPE.B.lab[,"LAB",gene], "///" )
		split.comp <- strsplit( TYPE.B.lab[,comp,gene], "///" )
		for ( s in 1:nrow(TYPE.B.lab) ) {
			samp <- rownames(TYPE.B.lab)[s]
			if ( length(split.comp[[samp]])>0 ) {
				prec <- nchar( split.lab[[samp]] )
				temp_conc.fwd.1 <- length(which( split.lab[[samp]] == c( substr(split.comp[[samp]][1],1,prec[1]),substr(split.comp[[samp]][2],1,prec[2]) ) ))
				temp_conc.fwd.2 <- length(which( split.lab[[samp]] == c( substr(split.comp[[samp]][1],1,prec[2]),substr(split.comp[[samp]][2],1,prec[1]) ) ))
				temp_conc.rev.1 <- length(which( split.lab[[samp]] == rev(c( substr(split.comp[[samp]][1],1,prec[1]),substr(split.comp[[samp]][2],1,prec[2]) )) ))
				temp_conc.rev.2 <- length(which( split.lab[[samp]] == rev(c( substr(split.comp[[samp]][1],1,prec[2]),substr(split.comp[[samp]][2],1,prec[1]) )) ))
				# temp_conc.rev <- length(which( split.lab[[samp]] == rev(split.comp[[samp]]) ))
				TYPE.B.lab.conc[samp,comp,gene] <- max( temp_conc.fwd.1, temp_conc.fwd.2, temp_conc.rev.1, temp_conc.rev.2  )
			}else{
				TYPE.B.lab.conc[samp,comp,gene] <- NA
			}
		}
	}
}

## Calculate % Concordance for Each Gene & Platform
TYPE.B.lab.c.a <- array( ,c( length(To_Compare), dim(TYPE.B.lab)[3] ))
colnames(TYPE.B.lab.c.a) <- dimnames(TYPE.B.lab)[[3]]
rownames(TYPE.B.lab.c.a) <- To_Compare
 # ...and determine # of Patients in Intersection of Platforms
TYPE.B.lab.c.b <- array( ,c( length(To_Compare), dim(TYPE.B.lab)[3] ))
colnames(TYPE.B.lab.c.b) <- dimnames(TYPE.B.lab)[[3]]
rownames(TYPE.B.lab.c.b) <- To_Compare
 # Number of Patients w/ 0,1,2 correct alleles
TYPE.B.lab.c.c <- array( ,c( 3, dim(TYPE.B.lab)[3], length(To_Compare) ))
rownames(TYPE.B.lab.c.c) <- 0:2
colnames(TYPE.B.lab.c.c) <- dimnames(TYPE.B.lab)[[3]]
dimnames(TYPE.B.lab.c.c)[[3]] <- To_Compare
 # Loop through Genes & Patients
for ( g in 1:length(TYPE.B.lab.genes) ) {
	gene <- TYPE.B.lab.genes[g]
	TYPE.B.lab.c.a[,gene] <- apply( TYPE.B.lab.conc[,To_Compare,gene], 2, function(x) mean(as.numeric(x),na.rm=T)/2 )
	TYPE.B.lab.c.b[,gene] <- apply( TYPE.B.lab.conc[,To_Compare,gene], 2, function(x) length(which(!is.na(x))) )
	TYPE.B.lab.c.c[,gene,] <- rbind( apply(TYPE.B.lab.conc[,To_Compare,gene],2,function(x) length(which(x==0))), apply(TYPE.B.lab.conc[,To_Compare,gene],2,function(x) length(which(x==1))), apply(TYPE.B.lab.conc[,To_Compare,gene],2,function(x) length(which(x==2))) )
}

#############################################################
## PLOT CONCORDANCE #########################################
#############################################################

###############################################
## LAB vs ... #################################

## CONCORDANCE By Gene ##
N.Gene.c <- ncol(TYPE.B.lab.c.a)
COLS <- c("tomato2","slateblue3","chartreuse2")
COLS <- c("chocolate2","slateblue3","dodgerblue1")
COLS <- c("tomato2","slateblue3","dodgerblue1")
XLIM <- c(1,N.Gene.c)
YLIM <- c(0,1)
png( paste(PathToPlot,"1-Lab_Concordance.png",sep=""),height=800,width=2400,pointsize=36 )
par(mfrow=c(1,3))
 # 2/4 Digit Precision
plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,xlab="Gene",ylab="% Concordant",main="Concordance w/ Lab HLA Types by Platform",xaxt="n")
axis( 1,at=1:N.Gene.c,label=colnames(TYPE.B.lab.c.a),las=1 )
abline( h=seq(0,1,.1),lty=3,col="grey50",lwd=1 )
for ( i in 1:nrow(TYPE.B.lab.c.a) ) {
	points( 1:N.Gene.c, TYPE.2.lab.c.a[i,], col=COLS[i],type="o",lty=2,lwd=3,pch=1 )
	points( 1:N.Gene.c, TYPE.4.lab.c.a[i,], col=COLS[i],type="o",lty=4,lwd=3,pch=2 )
}
legend( "bottomright", lty=c(1,2,3),pch=c(19,1,2),legend=c("Lab Best","2-Digit","4-Digit"),title="Precision",bg="white",lwd=3 )
 # Lab Best Precision
# YLIM <- c(0.5,1)
plot( 0,0,type="n",xlim=XLIM,ylim=YLIM,xlab="Gene",ylab="% Concordant",main="Concordance w/ Lab HLA Types by Platform",xaxt="n")
axis( 1,at=1:N.Gene.c,label=colnames(TYPE.B.lab.c.a),las=1 )
abline( h=seq(0,1,.1),lty=3,col="grey50",lwd=1 )
for ( i in 1:nrow(TYPE.B.lab.c.a) ) {
	points( 1:N.Gene.c, TYPE.B.lab.c.a[i,], col=COLS[i],type="o",lty=1,lwd=3,pch=19 )
}
legend( "bottomleft", fill=COLS,legend=rownames(TYPE.2.lab.c.a),title="Platform",bg="white" )
## Number of Patients By Gene
YLIM <- c(0,20)
 # 2/4 Digit Precision
barplot( TYPE.B.lab.c.b, beside=T, col=COLS )
abline( h=seq(0,20,5),lty=3,col="grey50",lwd=1 )
barplot( TYPE.B.lab.c.b, beside=T, col=COLS, add=T,main="Number of Patients Typed by Platform",xlab="Gene",ylab="# Patients") # ,legend=T,args.legend=c(x="topright",bg="white") )
dev.off()

## TABLE OF CONCORDANCE By Gene ##
COLS.tab <- sapply( COLS, function(x) colorRampPalette(c("white",x,"black"))(6)[2:4] )
XLIM <- c(0,20)
YLIM <- c(0,1)
png( paste(PathToPlot,"1-Lab_Concordance_Table.png",sep=""),height=800,width=2400,pointsize=36 )
par(mfrow=c(1,3))
 # Best Lab Digit Precision
barplot( prop.table(TYPE.B.lab.c.c,c(2,3))[,,"SOP"], space=c(0,rep(3,4)), beside=F, col=COLS.tab[,1], ylim=YLIM,xaxt="n",xlim=XLIM,xlab="Gene",main="Concordant Alleles per Person: Best Lab Precision",ylab="Fraction Patients w/ 0,1,2 Concordant Alleles" )
barplot( prop.table(TYPE.B.lab.c.c,c(2,3))[,,"CHP"], space=c(1,rep(3,4)), beside=F, col=COLS.tab[,2], add=T,xaxt="n" )
barplot( prop.table(TYPE.B.lab.c.c,c(2,3))[,,"SEQ"], space=c(2,rep(3,4)), beside=F, col=COLS.tab[,3], add=T,xaxt="n" )
axis( 1, at=seq(1.5,20,4), labels=colnames(TYPE.B.lab.c.c) )
 # 2-Digit Precision
barplot( prop.table(TYPE.2.lab.c.c,c(2,3))[,,"SOP"], space=c(0,rep(3,4)), beside=F, col=COLS.tab[,1], ylim=YLIM,xaxt="n",xlim=XLIM,xlab="Gene",main="Concordant Alleles per Person: 2-Digit Precision",ylab="Fraction Patients w/ 0,1,2 Concordant Alleles" )
barplot( prop.table(TYPE.2.lab.c.c,c(2,3))[,,"CHP"], space=c(1,rep(3,4)), beside=F, col=COLS.tab[,2], add=T,xaxt="n" )
barplot( prop.table(TYPE.2.lab.c.c,c(2,3))[,,"SEQ"], space=c(2,rep(3,4)), beside=F, col=COLS.tab[,3], add=T,xaxt="n" )
axis( 1, at=seq(1.5,20,4), labels=colnames(TYPE.2.lab.c.c) )
legend( "topright",title="Conc Alleles",legend=2:0,fill=paste("grey",c(20,50,80),sep="") )
 # 4-Digit Precision
barplot( prop.table(TYPE.4.lab.c.c,c(2,3))[,,"SOP"], space=c(0,rep(3,4)), beside=F, col=COLS.tab[,1], ylim=YLIM,xaxt="n",xlim=XLIM,xlab="Gene",main="Concordant Alleles per Person: 4-Digit Precision",ylab="Fraction Patients w/ 0,1,2 Concordant Alleles" )
barplot( prop.table(TYPE.4.lab.c.c,c(2,3))[,,"CHP"], space=c(1,rep(3,4)), beside=F, col=COLS.tab[,2], add=T,xaxt="n" )
barplot( prop.table(TYPE.4.lab.c.c,c(2,3))[,,"SEQ"], space=c(2,rep(3,4)), beside=F, col=COLS.tab[,3], add=T,xaxt="n" )
axis( 1, at=seq(1.5,20,4), labels=colnames(TYPE.4.lab.c.c) )
dev.off()

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
