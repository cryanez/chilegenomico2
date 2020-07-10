args <- commandArgs(TRUE)
name_set <- args[1]
name_popinfo <- args[2]
modo <- args[3]
if(modo==2){
  Q5 <- args[4]
  startValue <- args[5]
  workDirectory <- args[6]
  path_output <- file.path( paste(workDirectory, "/output/"  ,sep = ""))
  path_data  <- file.path(paste(workDirectory, "/output/"  ,sep = ""))
  path_input  <- file.path(paste(workDirectory, "/output/"  ,sep = ""))
  path_fam <- file.path(paste(workDirectory, "/output/"  ,sep = ""))
  CLG2_Ancestry <- args[7]
}else{
  workDirectory <- args[4]
  path_output <- file.path( paste(workDirectory, "/temporal/"  ,sep = ""))
  path_data  <- file.path(paste(workDirectory, "/temporal/"  ,sep = ""))
  path_input  <- file.path(paste(workDirectory, "/temporal/"  ,sep = ""))
  path_fam <- file.path(paste(workDirectory, "/temporal/"  ,sep = ""))
}
# Set working directory
if(basename(getwd())=="code"){
  if(!file.exists(path_data))
    stop("Input directory doesn't exist: ", path_data)
  if(!file.exists(path_output))
    dir.create(path_output, mode = "0755", recursive=T)
}
sortAncestryComps <- function(x, pops, poporder, rename=T) {
    if(! all(pops %in% poporder)) {
      stop("poporder is incomplete. It must contain the name of all populations in vector pops, in the desired order.")
    }
    poporder <- poporder[poporder %in% pops]
    anc.mean <- apply(x, 2, tapply, pops, mean)
    anc.mean <- anc.mean[match(poporder, row.names(anc.mean)),]
    coli <- 1:ncol(anc.mean)
    usedi <- NULL
    anc.out <- NULL
    outcol  <- NULL
    outpop  <- NULL
    for(j in 1:nrow(anc.mean)) {
        if(length(coli)<1) break
        jmax <- which.max(anc.mean[j,coli])
        anc.out <- cbind(anc.out, x[,coli[jmax]])
        outcol <- c(outcol, colnames(x)[coli[jmax]])
        outpop <- c(outpop, poporder[j])
        usedi <- c(usedi, coli[jmax])
        coli <- setdiff(coli, usedi)
    }
    if(rename) {
      colnames(anc.out) <- outpop
    } else {
      colnames(anc.out) <- outcol
    }
    if(length(coli)>0)
      anc.out <- cbind(anc.out, x[,coli])
    rownames(anc.out) <- rownames(x)
    anc.out
}
adx_files <- dir(path_data, patt=".Q$", full.names=T)
errors  <- read.table(file.path(paste(path_data,"CV_", name_set,".txt", sep = "")))
popinfo <- read.table(file.path(path_data, name_popinfo),header=T,sep="\t")
fam     <- read.table(file.path(path_fam, paste(name_set,".fam", sep = "")),sep=" ")
samps <- fam[,2]
popinfo <- popinfo[popinfo[,"IndID"] %in% samps,]
cverr <- errors[,2]
K     <- errors[,1]
adx_files <- adx_files[order(K)]
cverr <- cverr[order(K)]
K <- K[order(K)]
adx_files <- paste( paste(path_data, name_set, sep = ""), K, "Q", sep=".")
popinfo$slab <- sprintf("%-10s %8s", popinfo$IndID, popinfo$Ancestry)
x_pop <- tapply((1:nrow(popinfo))*1.2, popinfo$Ancestry, mean)
m_pop <- tapply((1:nrow(popinfo))*1.2, popinfo$Ancestry, max)
x_anc <- tapply((1:nrow(popinfo))*1.2, list(popinfo$Ancestry, popinfo$Ancestry), mean)
x_anc_ls <- list()
for(anc in rownames(x_anc)) {
  x_anc_ls[[anc]] <- x_anc[anc, !is.na(x_anc[anc,])]
}
x_anc_vec <- unlist(x_anc_ls)
names(x_anc_vec) <- sub("\\..*$", "", names(x_anc_vec))
pops  <- unique(as.character(popinfo$Ancestry))
N_pop <- length(pops)
comps_palette <- c("brown3", "forestgreen", "darkslategray2", "darkorange1", "mediumorchid4", "lightpink", "darkblue","yellow1", "darkolivegreen2", "lightpink1", "goldenrod1", "cadetblue3", "burlywood", "lightskyblue3", "rosybrown4")
palette(comps_palette)
pdf(file.path(path_output, paste("CV_",name_set,".pdf", sep = "")), width=6.5, height=4)
plot(cverr~K, ylab="error VC", type="b", col="darkblue", las=1, main="Error de validacion cruzada ADMIXTURE")
dev.off()
pdf(file.path(path_output,  paste("ADMIXTURE_K_4_",name_set,".pdf", sep = "")),width=20,height=5)
    for (i in 1) {
    qtbl <- read.table(adx_files[i])
    qtbl <- qtbl[match(popinfo[,"IndID"], samps),]
    qtbl <- sortAncestryComps(qtbl, popinfo$Ancestry, unique(popinfo$Ancestry))
    qtbl <- t(as.matrix(qtbl))#[,rev(1:ncol(qtbl))])
    barplot(qtbl, main=sprintf("K=%i", K[i]), xlab="sample", ylab="Ancestry", col=1:nrow(qtbl), border=NA, names.arg=popinfo$IndID, las=2,xaxt='n')
    text(x_pop, 1.05, names(x_pop), srt=45, font=2, cex=.7, adj=0, xpd=NA)
    text((1:ncol(qtbl))+cumsum(c(0,rep(.2,ncol(qtbl)-1)))-.5, par("usr")[3]-0.05, srt = 90, adj = 1 , labels = popinfo$slab, xpd = TRUE, cex=0.2)
    abline(v=m_pop, col="blue", lwd=.2)
  }
dev.off()
if(modo==2){
  Q5_table <- read.table(file = paste(path_data,Q5,sep = ""), header = F)
  startTail = as.integer(startValue)-1
  Q5_table <- tail(Q5_table, n=-(startTail))
  avg_AFR <- mean(Q5_table$V1)
  avg_EUR <- mean(Q5_table$V2)
  avg_AYM <- mean(Q5_table$V3)
  avg_Chile_Sur <- mean(Q5_table$V4)
  avg_AMR <- mean(Q5_table$V5)
  sd_AFR <- sd(Q5_table$V1)
  sd_EUR <- sd(Q5_table$V2)
  sd_AYM <- sd(Q5_table$V3)
  sd_Chile_Sur <- sd(Q5_table$V4)
  sd_AMR <- sd(Q5_table$V5)
  get_SE <- function(x,y){
    value <- x/sqrt(y)
    return(value)
  } 
  sE_AFR <- get_SE(sd_AFR, length(Q5_table$V1))
  sE_EUR <- get_SE(sd_EUR, length(Q5_table$V2))
  sE_AYM <- get_SE(sd_AYM, length(Q5_table$V3))
  sE_Chile_Sur <- get_SE(sd_Chile_Sur, length(Q5_table$V4))
  sE_AMR <- get_SE(sd_AMR, length(Q5_table$V5))
  CLG2_An <- read.table(CLG2_Ancestry, sep = "\t", header = T)   
  CLG_avg_AFR <- mean(CLG2_An$AFR)
  CLG_avg_EUR <- mean(CLG2_An$EUR)
  CLG_avg_AYM <- mean(CLG2_An$AYMARA)
  CLG_avg_Chile_Sur <- mean(CLG2_An$CHILE_SUR)
  CLG_avg_AMR <- mean(CLG2_An$AMERINDIO)
  CLG_sd_AFR <- sd(CLG2_An$AFR)
  CLG_sd_EUR <- sd(CLG2_An$EUR)
  CLG_sd_AYM <- sd(CLG2_An$AYMARA)
  CLG_sd_Chile_Sur <- sd(CLG2_An$CHILE_SUR)
  CLG_sd_AMR <- sd(CLG2_An$AMERINDIO)
  get_SE <- function(x,y){
    value <- x/sqrt(y)
    return(value)
  } 
  CLG_sE_AFR <- get_SE(CLG_sd_AFR, length(CLG2_An$AFR))
  CLG_sE_EUR <- get_SE(CLG_sd_EUR, length(CLG2_An$EUR))
  CLG_sE_AYM <- get_SE(CLG_sd_AYM, length(CLG2_An$AYMARA))
  CLG_sE_Chile_Sur <- get_SE(CLG_sd_Chile_Sur, length(CLG2_An$CHILE_SUR))
  CLG_sE_AMR <- get_SE(CLG_sd_AMR, length(CLG2_An$AMERINDIO))
  colnames(Q5_table)[colnames(Q5_table)=="V1"] <- "AFR"
  colnames(Q5_table)[colnames(Q5_table)=="V2"] <- "EUR"
  colnames(Q5_table)[colnames(Q5_table)=="V3"] <- "AYMARA"
  colnames(Q5_table)[colnames(Q5_table)=="V4"] <- "CHILE_SUR"
  colnames(Q5_table)[colnames(Q5_table)=="V5"] <- "AMERINDIO"
  ancestrias <- colnames(Q5_table)
  list_pval_dif <- NA
  for (an in ancestrias) {
    print(an)
    pvalue_vardif <- var.test(Q5_table[,an] ,   CLG2_An[,an], alternative = "two.sided" )
    pvalue_vardif <- as.numeric(pvalue_vardif$p.value)
    if(pvalue_vardif >= 0.05){
      pvalue_tdif <- t.test(Q5_table[,an] ,   CLG2_An[,an], paired = F, var.equal = T)
    }else{
      pvalue_tdif <- t.test(Q5_table[,an] ,   CLG2_An[,an], paired = F, var.equal = F)
    }
    pvalue_tdif <- as.numeric(pvalue_tdif$p.value)
    print(pvalue_tdif)
    list_pval_dif <- append(list_pval_dif,pvalue_tdif)
  }
  list_pval_dif <- list_pval_dif[-1]
  print("Ancestry;Average;SD;SE;Average_CLG;SD_CLG;SE_CLG;Dif_Valor_P")
  print(paste("Africano",round(avg_AFR, 4),round(sd_AFR,4),round(sE_AFR,4), round(CLG_avg_AFR,4), round(CLG_sd_AFR,4), round(CLG_sE_AFR,4), list_pval_dif[1] ,sep=";"))
  print(paste("Europeo",round(avg_EUR,4),round(sd_EUR,4),round(sE_EUR,4), round(CLG_avg_EUR,4), round(CLG_sd_EUR,4), round(CLG_sE_EUR,4), list_pval_dif[2] ,sep=";"))
  print(paste("Chile_Norte",round(avg_AYM,4),round(sd_AYM,4),round(sE_AYM,4), round(CLG_avg_AYM,4), round(CLG_sd_AYM,4), round(CLG_sE_AYM,4), list_pval_dif[3] ,sep=";"))
  print(paste("Chile_Sur",round(avg_Chile_Sur,4),round(sd_Chile_Sur,4),round(sE_Chile_Sur,4), round(CLG_avg_Chile_Sur,4), round(CLG_sd_Chile_Sur,4), round(CLG_sE_Chile_Sur,4), list_pval_dif[4] ,sep=";"))
  print(paste("Amerindio",round(avg_AMR,4),round(sd_AMR,4),round(sE_AMR,4), round(CLG_avg_AMR,4), round(CLG_sd_AMR,4), round(CLG_sE_AMR,4), list_pval_dif[5] ,sep=";"))
  listPrint <- "Ancestry;Average;SD;SE;Average_CLG;SD_CLG;SE_CLG;Dif_Valor_P"
  listPrint <- append(listPrint,   paste("Africano",round(avg_AFR, 4),round(sd_AFR,4),round(sE_AFR,4), round(CLG_avg_AFR,4), round(CLG_sd_AFR,4), round(CLG_sE_AFR,4), list_pval_dif[1] ,sep=";")    )
  listPrint <- append(listPrint,   paste("Europeo",round(avg_EUR,4),round(sd_EUR,4),round(sE_EUR,4), round(CLG_avg_EUR,4), round(CLG_sd_EUR,4), round(CLG_sE_EUR,4), list_pval_dif[2] ,sep=";")    )
  listPrint <- append(listPrint,   paste("Chile_Norte",round(avg_AYM,4),round(sd_AYM,4),round(sE_AYM,4), round(CLG_avg_AYM,4), round(CLG_sd_AYM,4), round(CLG_sE_AYM,4), list_pval_dif[3] ,sep=";")    )
  listPrint <- append(listPrint,   paste("Chile_Sur",round(avg_Chile_Sur,4),round(sd_Chile_Sur,4),round(sE_Chile_Sur,4), round(CLG_avg_Chile_Sur,4), round(CLG_sd_Chile_Sur,4), round(CLG_sE_Chile_Sur,4), list_pval_dif[4] ,sep=";")    )
  listPrint <- append(listPrint,   paste("Amerindio",round(avg_AMR,4),round(sd_AMR,4),round(sE_AMR,4), round(CLG_avg_AMR,4), round(CLG_sd_AMR,4), round(CLG_sE_AMR,4), list_pval_dif[5] ,sep=";")    )
  write.table(listPrint, file = paste(path_output, "Summary_ancestry_",name_set,".csv" ,sep = ""), quote = F, row.names = F , col.names = F)
}
