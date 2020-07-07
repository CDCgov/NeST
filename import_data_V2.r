#### provide genotype datasheet
#data = read.csv("Example_haplotype_data_sheet.txt",sep = "\t")
#print(data[4])
setwd('/home/momo/Desktop/CDC/Bayesian2/pyamd')
data = read.csv("Example_haplotype_data_sheet.txt",skip=0,stringsAsFactors = FALSE, sep = "\t")
#print(data)
#print(data[4])

## provide names of loci as per the genotype datasheet provided
locinames = c("CDC1","CDC2","CDC3","CDC4","X378_PART_A","X378_PART_B","X378_PART_C",
		    "X360i2_PART_A","X360i2_PART_B","X360i2_PART_C","X360i2_PART_D",
		    "Junction","MSR_Left","MSR_Right")
locinames_base = c("CDC1","CDC2","CDC3","CDC4","X378","X360","Junction","MSR")

## provide ploidy of each locus - ordered the same as the locinames variable
ploidy = c(2,2,2,2,2,2,2,2,2,2,2,1,1,1)

### is.na means not missing data
data = data[!is.na(data$Seq_ID) & data$Seq_ID != "",]
#print(data)
ids = data$Seq_ID
nids = length(ids)
nloci = length(locinames)

### transform data 

newdata = c()
datacompleteness = c()
# MaxMOI for each locus 
#print(nloci)
#print(locinames)
for (j in 1:nloci) {
	#print(paste(locinames[j],"",sep=""))
	#print(colnames(data))
	locicolumns = grepl(paste(locinames[j],"",sep=""),colnames(data))
	###why true and fasle????
	#print(locicolumns)
	#print(locicolumns)
	raw_alleles = data[,locicolumns]
	#print(raw_alleles)
	#print("working")
	#break
	#print(raw_alleles)
	maxMOI = max(c(2,rowSums(raw_alleles == "X",na.rm=TRUE)))
	#print(rowSums(raw_alleles == "X",na.rm=TRUE))
	#print(maxMOI)
	#break
	MOI = rowSums(raw_alleles == "X",na.rm=TRUE)
	#print(MOI)
	nalleles = sum(locicolumns,na.rm=TRUE)
	#print(nalleles)
	newdatatemp = rbind(matrix(NA,nids,maxMOI))
	#print(newdatatemp)
	sapply(1:nids,function (x) if (MOI[x] > 0) { newdatatemp[x,1:MOI[x]] <<- paste("Hap_",which(raw_alleles[x,] == "X"),sep="")})
	#print(sapply(1:nids,function (x) if (MOI[x] > 0) { newdatatemp[x,1:MOI[x]] <<- paste("Hap_",which(raw_alleles[x,] == "X"),sep="")}))
	if (ploidy[j] > 1) {
		sapply(1:nids,function (x) if (MOI[x] == 1) { newdatatemp[x,1:2] <<- paste("Hap_",which(raw_alleles[x,] == "X"),sep="")})
		#print(sapply(1:nids,function (x) if (MOI[x] == 1) { newdatatemp[x,1:2] <<- paste("Hap_",which(raw_alleles[x,] == "X"),sep="")}))
	
	}
	#break
	#print(sapply(1:nids,function (x) if (MOI[x] == 1) { newdatatemp[x,1:2] <<- paste("Hap_",which(raw_alleles[x,] == "X"),sep="")}))
	colnames(newdatatemp ) = paste(1:maxMOI,"_Hap_",locinames[j],sep="")
	#print(colnames(newdatatemp ))
	newdata = cbind(newdata,newdatatemp )
	#print(newdatatemp)
	#print(newdata)
	datacompleteness  = cbind(datacompleteness,MOI)
	#print(datacompleteness)
}

data = cbind(ids,newdata)
#print(data)
data = data.frame(data)

# calculate locus-specific amplification success
		       
datacompleteness_bylocus = sapply(locinames_base,function (x) rowSums(cbind(datacompleteness[,grepl(x,locinames)]))>0)

#print(datacompleteness)
#print(datacompleteness_bylocus)

#print(rowSums(datacompleteness_bylocus))
print("Number of loci sequenced")
table(rowSums(datacompleteness_bylocus))



# define inclusion criteria below - filters genotype datasheet to exclude specimens that fail to meet these criteria

cleandata = data[(rowSums(datacompleteness_bylocus[,c(5,7,8)]) == 3 & rowSums(datacompleteness_bylocus) >= 4) | rowSums(datacompleteness_bylocus) >= 5 | (rowSums(datacompleteness_bylocus[,c(6,7,8)]) == 3 & rowSums(datacompleteness_bylocus) >= 4) | (rowSums(datacompleteness_bylocus[,c(5,6,7)]) == 3 & rowSums(datacompleteness_bylocus) >= 4) | (rowSums(datacompleteness_bylocus[,c(5,6,8)]) == 3 & rowSums(datacompleteness_bylocus) >= 4),] 

data = cleandata
ids = data$ids
nids = length(ids)

####write.csv(data,"transformed_data.csv")
