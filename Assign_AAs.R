## Get Amino Acids for each patient ##
## May 12, 2015 ##
## Kristopher Standish ##

#############################################################
## LOAD DATA ################################################
#############################################################

## Set Date
DATE <- gsub("-","",Sys.Date())

## Set Paths to HLA Data Sets
PathToTypes <- "/Users/kstandis/Data/Burn/Data/HLA/20150512_SOAP_HLA_Types/20150511_HLA_Types.Rdata"
PathToRefs <- "/Users/kstandis/Data/HLI_Phase/20150223_HLA_Ref/Alignments_Rel_3190/"
 # http://www.ebi.ac.uk/ipd/imgt/hla/nomenclature/alignments.html
PathToOut <- "/Users/kstandis/Data/Burn/Data/HLA/20150512_SOAP_HLA_Types/"
dir.create( PathToOut )

## Load Janssen HLA Results
load( PathToTypes )

## Get Gene List
GENE_LIST <- names(COMPILE$GENES.list)
gene <- GENE_LIST[1]

FILES <- list.files(PathToRefs)
gsub("_prot.txt","",FILES[grep("_prot",FILES)])

#############################################################
## CREATE FLAT AA TABLE FOR EACH GENE #######################
#############################################################

## Load Reference Data
AA <- list()
for ( g in 1:length(GENE_LIST) ) {
	gene <- GENE_LIST[g]
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

#############################################################
## ASSIGN AMINO ACIDS to PATIENTS ###########################
#############################################################

## Pull out Individual Types per Patient
DAT <- COMPILE$GENES

## Set up
PAT_AA <- AMB <- list()
for ( g in 1:length(AA) ) {
	gene <- names(AA)[g]
	GENE <- AA[[gene]]
	TEMP_TAB <- array( ,c(ncol(DAT),ncol(GENE)) )
	colnames(TEMP_TAB) <- colnames(GENE)
	rownames(TEMP_TAB) <- colnames(DAT)
	AMB[[gene]] <- array(,c(0,2))
	for ( row in 1:nrow(TEMP_TAB) ) {
		which_type <- DAT[gene,row]
		PREC <- nchar(which_type)
		if (is.na(which_type)) { next }
		which_row <- which( rownames(GENE)==which_type )
		if ( length(which_row)==0 ) {
			which_row <- which( substr(rownames(GENE),1,PREC)==which_type )
			if (length(which_row)>1) {
				which_row <- which_row[1]
				if ( PREC<5 ) {
					print(paste(gene,"-",row,"-",which_type))
					# AMB[[gene]] <- c( AMB[[gene]], paste(rownames(TEMP_TAB)[row],which_type,sep=":") )
					AMB[[gene]] <- rbind( AMB[[gene]], c(rownames(TEMP_TAB)[row],which_type) )
				}
			}
			if ( length(which_row)==0 ) {
				print(paste(gene,"-",row,"-",which_type,"- SKIP"))
				# AMB[[gene]] <- c( AMB[[gene]], paste(rownames(TEMP_TAB)[row],which_type,sep=":") )
				AMB[[gene]] <- rbind( AMB[[gene]], c(rownames(TEMP_TAB)[row],which_type) )
				next
			}
		}
		TEMP_TAB[row,] <- GENE[which_row,]
	}
	PAT_AA[[gene]] <- TEMP_TAB
}

# MAX <- max( apply(PAT_AA$DRB1,2,function(x) length(unique(x)) ) )
# TEMP <- t(sapply( apply(PAT_AA$DRB1,2,table), "[",1:MAX ))

#############################################################
## SAVE DATA ################################################
#############################################################

## Save Patients' Amino Acid Assignments
COMPILE <- list( AA, PAT_AA, AMB )
names(COMPILE) <- c("AA","PAT_AA","AMB")
save( COMPILE, file=paste(PathToOut,"20150512_HLA_AA.Rdata",sep="") )





















#############################################################
## END OF DOC ###############################################
#############################################################