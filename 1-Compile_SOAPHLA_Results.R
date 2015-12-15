## Compile & Reformat SOAP-HLA Results for Janssen Cohort ##
## February 23, 2015 ##
## Updated Dec 11, 2015 ##
## Kristopher Standish ##

library(xlsx)

#############################################################
## LOAD DATA ################################################
#############################################################

## Set Date
DATE <- gsub("-","",Sys.Date())

## Set Paths to Data Sets & Save Locations (TSCC)
PathToEurList <- "/projects/janssen/HLA/HLA_SOAP_Samp.txt"
PathToNonEurList <- "/projects/janssen/HLA/Non_EUR_Samp.txt"
PathToEurResults <- "/projects/janssen/HLA/20140917_SOAP_HLA/Test_Output/"
PathToNonEurResults <- "/projects/janssen/HLA/20141110_SOAP_HLA_Non_EUR/Output/"
PathToART3002Results <- "/projects/janssen/HLA/20150827_ART3002_SOAP_HLA/HLA_Types/"
PathToHLADB <- "/projects/janssen/HLA/Alignments_Rel_3190/"
PathToSave <- paste("/projects/janssen/HLA/SOAP_HLA_Types/",DATE,sep="")
dir.create( PathToSave )

## Load Janssen Sample Lists
Eur_List <- read.table( PathToEurList,colClasses="character" )[,1]
NonEur_List <- read.table( PathToNonEurList,colClasses="character" )[,1]
ART3002_List <- list.files( PathToART3002Results )

## Loop through Samples & Load Relevant Data
HLA.list <- list()
 # Europeans
for ( samp in Eur_List ) {
	Path <- paste( PathToEurResults, samp,"/",samp,".type", sep="" )
	HLA.list[[samp]] <- read.table( Path, sep="\t",header=F, fill=T, col.names=c("Type","Alt","Conf",paste("Ex",1:5,sep="_")), colClasses=c("character","character","numeric","character","character","character","character") )
	colnames(HLA.list[[samp]])[1:3] <- c("Type","Alt","Conf")
}
 # Non-Europeans
for ( samp in NonEur_List ) {
	Path <- paste( PathToNonEurResults, samp,"/",samp,".type", sep="" )
	HLA.list[[samp]] <- read.table( Path, sep="\t",header=F, fill=T, col.names=c("Type","Alt","Conf",paste("Ex",1:5,sep="_")), colClasses=c("character","character","numeric","character","character","character","character") )
	colnames(HLA.list[[samp]])[1:3] <- c("Type","Alt","Conf")
}
 # ART3002
for ( samp in ART3002_List ) {
	Path <- paste( PathToART3002Results, samp,"/",samp,".type", sep="" )
	HLA.list[[samp]] <- read.table( Path, sep="\t",header=F, fill=T, col.names=c("Type","Alt","Conf",paste("Ex",1:5,sep="_")), colClasses=c("character","character","numeric","character","character","character","character") )
	colnames(HLA.list[[samp]])[1:3] <- c("Type","Alt","Conf")
}

#############################################################
## ORGANIZE DATA ############################################
#############################################################

## Compile all gene names typed in any list
Compile_Gene_Names <- unique( unlist(lapply( HLA.list, function(x) sapply(strsplit( as.character(x[,"Type"]),"*",fixed=T),"[",1) )) )
Ind_Gene_Names <- lapply( HLA.list, function(x) sapply(strsplit( as.character(x[,"Type"]),"*",fixed=T),"[",1) )
Ind_Gene_Types <- lapply( HLA.list, function(x) sapply(strsplit( as.character(x[,"Type"]),"*",fixed=T),"[",2) )
Ind_Gene_Alt <- lapply( HLA.list, function(x) sapply(strsplit( as.character(x[,"Alt"]),"*",fixed=T),"[",2) )
Ind_Gene_Conf <- lapply( HLA.list, function(x) as.numeric(x[,"Conf"]) )
Nov_Gene_Vars <- lapply( HLA.list, function(x) apply(x[,4:8],1,function(y) paste(y,collapse="_")) )

## Specify Number of Samples
N_Samp <- length(HLA.list)
Names_Samp <- names(HLA.list)

#############################################################
## Create Table for each Gene (Num_Type rows X 2*Num_Samp columns)
Gene_Names <- Compile_Gene_Names[c(1:8,10:14)]

GENES.list <- CONF.list <- NOVEL.list <- list()
for ( g in 1:length(Gene_Names) ) {
	gene <- Gene_Names[g]
	WHICH <- lapply( Ind_Gene_Names, function(x) which(x==gene) )
	Gene_Types <- list()
	Novel_Vars <- list()
	Gene_Types <- lapply( Names_Samp, function(x)Ind_Gene_Types[[x]][ WHICH[[x]] ] )
	Gene_Alts <- lapply( Names_Samp, function(x)Ind_Gene_Alt[[x]][ WHICH[[x]] ] )
	Gene_Types.conf <- lapply( Names_Samp, function(x)Ind_Gene_Conf[[x]][ WHICH[[x]] ] )
	Gene_Alts.conf <- lapply( Gene_Types.conf, function(x) 100-x )
	Novel_Vars <- lapply( Names_Samp, function(x)Nov_Gene_Vars[[x]][ WHICH[[x]] ] )
	names(Gene_Types) <- names(Gene_Alts) <- names(Gene_Types.conf) <- names(Gene_Alts.conf) <- names(Novel_Vars) <- Names_Samp
	# Novel_Vars[ which(unlist(lapply(Novel_Vars,function(x)any(x!="____")))) ]
	# Gene_Types[ which(unlist(lapply(Novel_Vars,function(x)any(x!="____")))) ]
	## Compile Union of all "Types" in Cohort
	Union_Gene_Types <- sort( Reduce( union, Gene_Types ) )
	Union_Gene_Alts <- sort( Reduce( union, Gene_Alts ) )
	Compile_Gene_Types <- sort(union( Union_Gene_Types,Union_Gene_Alts ))
	## Compile HLA Type Table for Gene
	GENES.list[[gene]] <- array( 0, c(length(Compile_Gene_Types),2*N_Samp) )
	rownames(GENES.list[[gene]]) <- Compile_Gene_Types
	colnames(GENES.list[[gene]]) <- paste( rep(Names_Samp,rep(2,N_Samp)), 1:2,sep="_" )
	NOVEL.list[[gene]] <- CONF.list[[gene]] <- GENES.list[[gene]]
	# NOVEL.list[[gene]] <- array( , c(length(Compile_Gene_Types),2*N_Samp) )
	# rownames(NOVEL.list[[gene]]) <- Compile_Gene_Types
	# colnames(NOVEL.list[[gene]]) <- paste( rep(Names_Samp,rep(2,N_Samp)), 1:2,sep="_" )
	for ( s in 1:N_Samp ) {
		samp <- Names_Samp[s]
		Col_Ind <- grep(samp,colnames(GENES.list[[gene]]))
		Temp.Types <- Gene_Types[[samp]][ order(Gene_Types[[samp]]) ]
		Temp.Types.conf <- Gene_Types.conf[[samp]][ order(Gene_Types[[samp]]) ]
		Temp.Alts <- Gene_Alts[[samp]][ order(Gene_Types[[samp]]) ]
		Temp.Alts.conf <- Gene_Alts.conf[[samp]][ order(Gene_Types[[samp]]) ]
		Temp.Novel <- Novel_Vars[[samp]][ order(Gene_Types[[samp]]) ]
		if ( length(which(!is.na(Temp.Types)))>0 ) {
			GENES.list[[gene]][Temp.Types[1],Col_Ind[1]] <- 1
			CONF.list[[gene]][Temp.Types[1],Col_Ind[1]] <- Temp.Types.conf[1]
			NOVEL.list[[gene]][Temp.Types[1],Col_Ind[1]] <- Temp.Novel[1]
		}
		if ( !is.na(Temp.Alts[1]) ) {
			CONF.list[[gene]][Temp.Alts[1],Col_Ind[1]] <- Temp.Alts.conf[1]
		}
		if ( length(which(!is.na(Temp.Types)))>1 ) {
			GENES.list[[gene]][Temp.Types[2],Col_Ind[2]] <- 1
			NOVEL.list[[gene]][Temp.Types[2],Col_Ind[2]] <- Temp.Novel[2]
			CONF.list[[gene]][Temp.Types[2],Col_Ind[2]] <- Temp.Types.conf[2]
		}
		if ( !is.na(Temp.Alts[2]) ) {
			CONF.list[[gene]][Temp.Alts[2],Col_Ind[2]] <- Temp.Alts.conf[2]
		}
	}
	NOVEL.list[[gene]][which(NOVEL.list[[gene]]=="____",arr.ind=T)] <- 0
	print(paste("Done with",gene))
}

#############################################################
## Convert to Var Dosage Table (0/1/2 for each Type)
GENES.2.list <- CONF.2.list <- list()
for ( g in 1:length(Gene_Names) ) {
	gene <- Gene_Names[g]
	GENES.2.list[[gene]] <- GENES.list[[gene]][,seq(1,2*N_Samp,2)] + GENES.list[[gene]][,seq(2,2*N_Samp,2)]
	CONF.2.list[[gene]] <- CONF.list[[gene]][,seq(1,2*N_Samp,2)] + CONF.list[[gene]][,seq(2,2*N_Samp,2)]
	# GENES.2.list[[gene]] <- array( ,c(length(Gene_Names),N_Samp) )
	colnames(GENES.2.list[[gene]]) <- Names_Samp
}

#############################################################
## Compile Types by Name into Single Table for All Genes
GENES <- array( ,c(length(Gene_Names),2*N_Samp) )
colnames(GENES) <- paste( rep(Names_Samp,rep(2,N_Samp)), 1:2,sep="_" )
rownames(GENES) <- Gene_Names
 # Create Table to Compile Novel Variants
NOVEL <- CONF <- GENES
# NOVEL <- array( ,c(length(Gene_Names),2*N_Samp) )
# colnames(NOVEL) <- paste( rep(Names_Samp,rep(2,N_Samp)), 1:2,sep="_" )
# rownames(NOVEL) <- Gene_Names
for ( g in 1:length(Gene_Names) ) {
	gene <- Gene_Names[g]
	## Compile Types for each Gene/Sample
	temp <- unlist(apply( GENES.list[[gene]], 2, function(x) rownames(GENES.list[[gene]])[which(x==1)]  ))
	GENES[gene,names(temp)] <- temp
	## Compile Novel Variants for each Gene/Sample
	temp <- unlist(apply( NOVEL.list[[gene]], 2, function(x) unname(x[which(!is.na(x))]) )) # rownames(NOVEL.list[[gene]])[which(!is.na(x))]  ))
	NOVEL[gene,names(temp)] <- temp
	## Compile
}


#############################################################
## WRITE DATA ###############################################
#############################################################

## Compile Data into Single List
COMPILE <- list( GENES.list, GENES.2.list, GENES, NOVEL.list, NOVEL, CONF.list, CONF.2.list )
names(COMPILE) <- c( "GENES.list", "GENES.2.list", "GENES", "NOVEL.list", "NOVEL", "CONF.list", "CONF.2.list" )

## Write Rdata w/ Gene Arrays
save( COMPILE, file=paste(PathToSave,"/",DATE,"_HLA_Types.Rdata",sep="") )








# PROT <- list()
# gene <- "A"
# PROT[[gene]] <- read.table( paste("/projects/janssen/HLA/Alignments_Rel_3190/")




#############################################################
## END OF DOC ###############################################
#############################################################
