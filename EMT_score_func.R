EMT76GS = function(log2Exp){


cat("Calculating EMT Score by using 76 gene signatures ....\n")
   

  remIdx = which(apply(log2Exp,1,function(x) any(x == NaN | x == -Inf)) == TRUE)
  if(length(remIdx) > 0 ) log2Exp = log2Exp[-remIdx, ]

  log2Exp[is.na(log2Exp)] = 0


  sampleNum = ncol(log2Exp)
  genes = rownames(log2Exp)
  exp = apply(log2Exp,2,as.numeric)

  EMTSignature = data.frame(read_excel("../Gene_signatures/76GS/EMT_signature_76GS.xlsx",col_names = TRUE))
  EMTIdx = match(unique(na.omit(EMTSignature[,2])),genes)
  geneFound = length(na.omit(EMTIdx))
  cat(paste(geneFound ,"gene's expression values found \n",sep = " "))
  
  ## Add sd
  sdVal = rnorm(ncol(exp), mean=0, sd= 0.01)
 
  EMTMat = t(apply(exp[na.omit(EMTIdx),], 1, function(x)  x + sdVal))
  row.names(EMTMat) = genes[na.omit(EMTIdx)]


 
  ## get the weights for each gene
  ecadhExp = grep("^CDH1$",row.names(EMTMat))
  if(length(ecadhExp) == 0 ){
    cat("CDH1 gene not found- 76 GS EMT score cannot be calculated\n") 
    EMTScoreStd = rep(0, ncol(exp))
  } else{
    ecadhExpVal = EMTMat[ecadhExp, ]
    weightVal = apply(EMTMat,1,function(x) cor(x,ecadhExpVal))
    EMTMatWt = weightVal*EMTMat
    EMTScore = apply(EMTMatWt,2,sum)
    EMTScoreMean = mean(EMTScore)
    EMTScoreStd = EMTScore-EMTScoreMean

  }  
  emtWrite = cbind(colnames(exp),EMTScoreStd)
  colnames(emtWrite) = c("Sample_name","GS76")
  return(emtWrite)

}


KSScore = function(expMat){

  cat("Calculating EMT Score by KS score method ....")

  genes = rownames(expMat)
  exp = apply(expMat,2,as.numeric)
  
  EMTSignature = data.frame(read_excel("../Gene_signatures/KS/EM_gene_signature_cellLine_KS.xlsx",col_names = FALSE))
  commonSig = intersect(EMTSignature[,1],genes)
  EMTExpIdx = match(commonSig,genes)
  EMTExp = exp[EMTExpIdx, ]
  EMTGrpIdx = match(commonSig,EMTSignature[,1])
  geneCat = EMTSignature[EMTGrpIdx,2]
  epiIdx = which(geneCat == "Epi")
  mesIdx = which(geneCat == "Mes")

    ## Perform KS test
    sampleScore2 = matrix(0,nrow=ncol(EMTExp),ncol=6)
    rownames(sampleScore2) = colnames(EMTExp)
    for(i in 1:ncol(EMTExp)){
        ## Two sided test
        ksTwoSided =  ks.test(EMTExp[mesIdx,i],EMTExp[epiIdx,i])
        ## One sided test: ecdf(Mes) > ecdf(Epi)
        ksResGrt = ks.test(EMTExp[mesIdx,i],EMTExp[epiIdx,i],alternative = "greater")
        ## One sided test: ecdf(Epi) > ecdf(Mes)
        ksResLess = ks.test(EMTExp[epiIdx,i],EMTExp[mesIdx,i],alternative = "greater")
        sampleScore2[i, ] = c(ksTwoSided$statistic,ksTwoSided$p.value,
                  ksResGrt$statistic,ksResGrt$p.value,
                  ksResLess$statistic,ksResLess$p.value)
    }

    ## Assign signs to EMT score of sample based on test statistic and pvalue
    finalScore = matrix(0,nrow = nrow(sampleScore2),ncol = 1)
    for (i in 1:nrow(sampleScore2)) {
        if(sampleScore2[i,4] < 0.05){
          finalScore[i, ] = c(-1 * sampleScore2[i,3])
        } else if (sampleScore2[i,6] < 0.05){
          finalScore[i, ] = sampleScore2[i,5]
        } else {

          if(sampleScore2[i,5] == max(c(sampleScore2[i,3],sampleScore2[i,5]))){
            finalScore[i, ] = max(c(sampleScore2[i,3],sampleScore2[i,5]))
          } else {
            finalScore[i, ] = (-1 * max(c(sampleScore2[i,3],sampleScore2[i,5])))
          }
          
      }
  }


  ksOut = cbind(colnames(EMTExp),finalScore)
  colnames(ksOut) = c("sampleName", "KS_score")
  return(finalScore)
}


emtPredExp = function(geneExpMat,gseID){

  #Predictor + Normalizer gene expression to run MLR model for EMT prediction 
  geneList =  scan(file = "../Gene_signatures/MLR/genes_for_EMT_score.txt",sep = '\n',what = "vector")
  idx = match(geneList,rownames(geneExpMat))
  geneExpSubMat = geneExpMat[na.omit(idx), ]

  notFound = setdiff(geneList,rownames(geneExpMat))
  print(paste("Gene not found in the dataset", paste(notFound, collapse = ",")))
  sampleNum = ncol(geneExpSubMat)
  
  fileName = trimws(paste(gseID,"EMT_gene_explevels.txt",sep="_"))
  outFileName = paste("../Data_generated",fileName,sep = "/")
  write.table(geneExpSubMat,file = outFileName ,sep = '\t',quote = F,row.names = T)
  return(list(geneExpSubMat,fileName))
}

getNCIindices = function(EMTexp,seriesID){
  nci60_data = read.table(file = "../Gene_signatures/MLR/GPL570-55999.txt",sep = '\t',header = T,stringsAsFactors = F, fill = T, quote = "")
  nci60_indices = na.omit(match(rownames(EMTexp),nci60_data[,11]))
  nci60_out = paste(seriesID,"_nci60_use_probe.txt",sep = "")
  nci60_indices_file = paste("../Data_generated",nci60_out,sep = "/")
  write(nci60_indices,file = nci60_indices_file ,sep = '\n')
  return(list(nci60_indices,nci60_out))
}


mlrEMTPred = function(expSet,gse_id){

  cat("EMT predictors expression file generated\n")
  emtPredExpOut = emtPredExp(expSet,gse_id)
  emtExpMat = emtPredExpOut[[1]]
  emtWrite = emtPredExpOut[[2]]
  

  NCIindicesOut = getNCIindices(emtExpMat,gse_id)
  nci60Indices = NCIindicesOut[[1]]
  nci60Write = NCIindicesOut[[2]]
  
  mlrOutfile  = paste(gse_id,"_emt_score_MLR.txt",sep  = "")
  writeOut = mlrOutfile
  write(rbind(emtWrite,nci60Write,writeOut),file = "../MLR3_Code/matlab_input_temp.txt", sep = '\n')

  cat("running MLR EMT predictions (matlab code)...\n")
  run_matlab_script("../MLR3_Code/MLR3_automated.m")  

  mlrOut = read.table(file  = paste("../Data_generated/",writeOut, sep = ""),header = T,stringsAsFactors = F)
  file.remove("../MLR3_Code/matlab_input_temp.txt")
  return(mlrOut)
}

all_scoreCor = function(scoreList){
  
  count = 0
  corVal = NULL
  for(i in 1:2){
    for(j in (i+1):3){
      count = count + 1
      corEst = cor.test(as.numeric(scoreList[[i]]),as.numeric(scoreList[[j]]))
      corVal = c(corVal,c(corEst$estimate,corEst$p.value))
    }
  }

  return(corVal)
}


xlsx.addEMTscores<-function(sheetName,GS76, KS, MLR, correlationVal){
  EMTDF = data.frame("Sample_name" = GS76[,1], 
                "GS76" = GS76[,2],
                "MLR" = MLR[,2],
                "KS" = KS[,1])
  write.xlsx(EMTDF, "../Output/all_EMT_Scores.xlsx",sheetName, row.names = F, append = TRUE)
  
}


