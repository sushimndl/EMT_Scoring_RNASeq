##==================================================================================================
## Raw read count to TPM conversion code
## Author: Priyanka Chakraborty
## Date: 06-08-2020
## References : 
##  https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
##  http://robpatro.com/blog/?p=235
##  https://www.biostars.org/p/390038/
##====================================================================================================



countToTpm <- function(rawCount, GSEID) {    

	geneLengthMap = read.table(file="../Annotation/mm10_gene_length_kb.txt", sep = "\t",header = TRUE,stringsAsFactors = F)
	geneCol = 2
	commonGenes = intersect(rawCount[,1], geneLengthMap[,1])
	
	if(length(commonGenes) == 0) {
		geneLengthMap = read.delim("../Annotation/hg38_gene_length_kb.txt")
		commonGenes = intersect(rawCount[,1], geneLengthMap[,1])
		geneCol = 2
	}


	idx1 = match(commonGenes,rawCount[,1])
	rawCountSubset = apply(rawCount[idx1,-c(1,2)], 2, as.numeric)
	idx2 = match(commonGenes,geneLengthMap[,1])
	geneLengthSubset = geneLengthMap[idx2, ]
	geneSymbol = geneLengthSubset[ ,geneCol]
	featureLength = geneLengthSubset[ ,3]

	rownames(rawCountSubset) = geneSymbol

	    # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(rawCountSubset))
  
  # Compute effective lengths of features in each library.
  effLen <- featureLength
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(rawCountSubset), function(i) {
    rate = (rawCountSubset[,i])/effLen
    rate/sum(rate) * 1e6
  }))
  
  log2TPM = log2(tpm+1)
  # Copy the row and column names from the original matrix.
  colnames(log2TPM) <- colnames(rawCountSubset)
  rownames(log2TPM) <- toupper(rownames(rawCountSubset))

	  outFile = paste("../Data_generated/", paste(GSEID, "_TPM.tsv"), sep = "")
	  if(!file.exists(outFile)) file.create(outFile)
	  write.table(log2TPM,file = outFile)

  return(log2TPM)
  
}

rnaToMA = function(TPMVal){
	MA = (0.57 +(TPMVal* 0.37 ))
	return(MA)
}


