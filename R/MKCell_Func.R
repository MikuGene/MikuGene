## MKcell Functions. ##
# Sun Jun 21 09:54:52 2020 ------------------------------

MKrcpp = F
MKtime = T

MK_time <- function(){
  mtime = gsub(" ", "", gsub(":", "", date()))
  mtime = sub("^...", " ", mtime)
  if(!MKtime){mtime = NULL}
  return(mtime)
}

message("  MKCell time system: From ", date(), " to", MK_time())

# Need nnls;Matrix;Seurat;ggplot2;harmony; if MKrcpp, neead Rtools.

options(stringsAsFactors = F)
if(!any(installed.packages() %in% "Matrix")){
  install.packages("Matrix")
}
suppressMessages(library(Matrix))

## MK_Tree 8a03a29901b31176e32928321b1349e6 ##
#
MK_Tree <- function(Bulk,Sigl,Cluster,SuClust,SuType = c("Cancer cell","Immune cell")){

  # check bulk #
  Check0 <- apply(Bulk,2,function(i) sum(i != 0)) > 0
  if(sum(Check0) == 0){stop(message("!!! Error, Bulk reads all 0 !!!", MK_time()))}
  Bulk <- Bulk[,Check0,drop = F]
  rm(Check0)
  if(ncol(data.frame(Bulk)) < 2){Bulk <- cbind(Bulk,Bulk)}

  # tree 1 #
  Re <- MK_Cell(Bulk,Sigl,Cluster)

  # tree 2 #
  Re_Tree <- list()
  Re_Type <- c()

  # save index #
  Te <- 0
  if(is.null(SuType)){
    message("Only print Tree 1.", MK_time())
    return(data.frame(Re[[1]]))
  }

  for (Type in SuType) {

    Te <- Te + 1

    # extract one type #
    Re_Bulk <- sapply(Re[[2]], function(i) i[,colnames(i) == Type])
    Re_Sigl <- Sigl[,Cluster == Type]
    colnames(Re_Bulk) <- colnames(Bulk)

    # check bulk #
    Check0 <- apply(Re_Bulk,2,function(i) sum(i != 0)) > 0
    if(sum(Check0) < 1){
      message("No tree 2, only print tree 1", MK_time())
      return(data.frame(Re[[1]]))
    }

    Re_Bulk <- Re_Bulk[,Check0]
    if(ncol(data.frame(Re_Bulk)) < 2){Re_Bulk <- cbind(Re_Bulk,Re_Bulk)}

    # rerun #
    Re_Sub <- MK_Cell(Re_Bulk,Re_Sigl,SuClust[Cluster == Type])

    # tree 2 proportion #
    Fix <- NULL
    if(sum(!Check0) > 0){
      Fix <- matrix(0,nrow = nrow(Re_Sub[[1]]),ncol = sum(!Check0))
      colnames(Fix) <- colnames(Re[[1]])[!Check0]}
    Re_Sub[[1]] <- cbind(Re_Sub[[1]],Fix)
    Re_Tree[[Te]] <- sweep(Re_Sub[[1]][,match(colnames(Re_Sub[[1]]),colnames(Re[[1]]))],2,Re[[1]][Type,],FUN = "*")
    Re_Type[Te] <- Type

    rm(Re_Sub)}

  # tree all proportion #
  ReAll <- data.frame(rbind(Re[[1]],do.call(rbind,Re_Tree)))

  # tree tag #
  ReAll$Tree <- c(rep("Tree 1",nrow(Re[[1]])),
                  rep(paste0("Tree 2 in ",Re_Type[1]),nrow(Re_Tree[[1]])),
                  rep(paste0("Tree 2 in ",Re_Type[2]),nrow(Re_Tree[[2]])))

  # cell proportion #
  rm(Re,Re_Tree,Re_Type,Te)
  gc()
  return(ReAll)}
#
## 8a03a29901b31176e32928321b1349e6 ##

## MK_Cell 8a03a29901b31176e32928321b1349e6 ##
#
MK_Cell <- function(Bulk,Sigl,Cluster) {

  # process scRNA #
  Sigl <- Sigl[Matrix::rowSums(Sigl) != 0,]

  Meta1 <- apply(Sigl, 2, sum)
  Meta2 <- apply(Sigl, 2, function(i) i/sum(i))

  # abudent variance library and design #
  Cluster <- as.character(Cluster)
  Abud <- sapply(unique(Cluster), function(i) rowMeans(Meta2[,Cluster %in% i]))
  Vari <- sapply(unique(Cluster), function(i) apply(Meta2[,Cluster %in% i],1,var))
  Libr <- sapply(unique(Cluster), function(i) mean(Meta1[Cluster %in% i]))
  Dsig <- sweep(Abud,2,Libr,FUN = "*")

  rm(Meta1,Meta2,Abud)
  gc()

  # common gene and match #
  Gene <- intersect(rownames(Dsig), rownames(Bulk))
  if (length(Gene) < 0.5*min(nrow(Bulk), nrow(Sigl))) {write.csv(c("Lower than 50% common genes!",MK_time()),"Error.csv")}
  Dsig <- Dsig[match(Gene, rownames(Dsig)),]
  Vari <- Vari[match(Gene, rownames(Dsig)),]
  Bulk <- apply(Bulk[match(Gene, rownames(Bulk)),], 2, function(i) i/sum(i))
  rm(Gene)
  gc()

  # weight-nnls #
  ReLm <- list()
  for (i in 1:ncol(Bulk)) {

    # process bulk #
    Tag <- (Bulk[, i] != 0)

    CoDsig <- Dsig[Tag,]
    CoVari <- Vari[Tag,]
    CoBulk <- Bulk[Tag,i]

    # center and normalize #
    CoDsig <- (CoDsig - mean(CoDsig)) / sd(CoDsig)
    CoBulk <- (CoBulk - mean(CoBulk)) / sd(CoBulk)

    # weight-nnls #
    ReLm[[i]] <- MK_Lm(CoDsig,CoBulk,CoVari,Libr)
    rm(Tag,CoDsig,CoVari,CoBulk)+gc()}

  # cell proportion #
  Cell <- sapply(ReLm, function(i) i$x/sum(i$x))
  rownames(Cell) <- colnames(Dsig)
  colnames(Cell) <- colnames(Bulk)

  # cell matrix #
  MaCell <- lapply(ReLm, function(i) sweep(Dsig,2,i$x,FUN = "*"))

  # results #
  rm(Dsig,ReLm,Vari,Libr)+gc()
  return(list(Cell,MaCell))}
#
## 8a03a29901b31176e32928321b1349e6 ##

## MK_Lm 8a03a29901b31176e32928321b1349e6 ##
#
MK_Lm <- function(CoDsig,CoBulk,CoVari,Libr){

  # nnls #
  if(!any(installed.packages() %in% "nnls")){
    install.packages("nnls")
  }
  library(nnls)
  Lm <- nnls(CoDsig,CoBulk)
  R <- resid(Lm)

  # weight gene #
  WgGene <- as.numeric(1/(1e-04 + R^2 + colSums((Lm$x * Libr)^2 * t(CoVari))))

  rm(Lm,R)
  gc()

  # weight bulk and design #
  WgCoBulk <- CoBulk * sqrt(WgGene)
  WgCoDsig <- CoDsig * sqrt(WgGene)

  # weight-nnls #
  WgLm <- nnls(WgCoDsig, WgCoBulk)
  WgP <- WgLm$x/sum(WgLm$x)
  WgR <- resid(WgLm)

  # iterate #
  for (j in 1:1000){

    WgGene <- as.numeric(1/(1e-04 + WgR^2 + colSums((WgLm$x * Libr)^2 * t(CoVari))))

    WgCoBulk <- CoBulk * sqrt(WgGene)
    WgCoDsig <- CoDsig * sqrt(WgGene)

    WgLm <- nnls(WgCoDsig, WgCoBulk)
    WgPn <- WgLm$x/sum(WgLm$x)
    WgRn <- resid(WgLm)

    # distance lower than 0.01 #
    if (sum(abs(WgPn - WgP)) < 0.01){
      rm(WgGene,WgCoBulk,WgCoDsig,WgP,WgR,WgPn,WgRn)
      gc()
      return(WgLm)}

    WgP <- WgPn
    WgR <- WgRn}

  # max iteration return #
  rm(WgGene,WgCoBulk,WgCoDsig,WgP,WgR,WgPn,WgRn)
  gc()
  return(WgLm)}
#
## 8a03a29901b31176e32928321b1349e6 ##

## MK_read 8a03a29901b31176e32928321b1349e6 ##
MK_read <- function(path){
  MTX = Matrix::readMM(paste0(path, "_mt.mtx"))
  MTX_cell = read.csv(paste0(path, "_cell.csv"), header = T, row.names = 1)$x
  MTX_gene = read.csv(paste0(path, "_gene.csv"), header = T, row.names = 1)$x
  MTX@Dimnames = list(MTX_gene, MTX_cell)
  rm(MTX_cell, MTX_gene)
  gc()
  return(MTX)}
#
MK_reads = function(path, IDin = NULL, verbose = T){
  name = gsub(".*/", "", path)
  if(is.null(IDin)){
    IDin = unique(gsub("_.*", "", gsub(".* ", "", list.files(path))))
  }
  MKfiles <- list()
  for (i in 1:max(as.numeric(IDin))) {
    if(verbose){message(" Read MM ", i, MK_time())}
    MKfile = MK_read(paste0(path, "/", name, " ", i))
    MKfiles[[i]] = MKfile
    rm(MKfile)
  }
  rm(IDin)

  ## return ##
  MK_comb = MKfiles[[1]]
  if(i == 1){
    rm(MKfiles)
    return(MK_comb)
  }else{
    for (i in 2:length(MKfiles)) {
      MK_comb = MK_cbind_s(MK_comb, MKfiles[[i]])
      message(" Cbinding ", i, MK_time())
    }
    rm(MKfiles)
    return(MK_comb)
  }
}
#
## 8a03a29901b31176e32928321b1349e6 ##

## MK_pie 8a03a29901b31176e32928321b1349e6 ##
#
MK_pie <- function(x,Title = "Cell Composition",Size = 3.5){
  library(ggplot2)
  Iden <- data.frame(table(x))
  Iden$Var2 <- factor(c(1:nrow(Iden)))
  Iden_label <- paste0(as.vector(Iden[,1]), " (", round(Iden$Freq / sum(Iden$Freq) * 100, 2), "%)")
  pie <- ggplot(Iden, aes(x = "", y = Iden$Freq, fill = Iden$Var2))+geom_bar(stat = "identity")+coord_polar(theta = "y")+labs(x = "", y = "", title = Title)+theme_light()+
    theme(axis.ticks = element_blank(),axis.text.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    scale_fill_discrete(breaks = Iden$Var2,labels = Iden_label)+theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5))+
    geom_text(aes(x = 1.5,y = rev(c(0,cumsum(rev(Iden$Freq))[-length(Iden$Freq)])+rev(Iden$Freq/2)),label = as.vector(Iden$Var2)),size = Size)
  return(pie)}
#
## 8a03a29901b31176e32928321b1349e6 ##

## MK_list 8a03a29901b31176e32928321b1349e6 ##
#
MK_box <- function(x, Size = 3) {
  library(ggplot2)
  Data <- list()
  for (i in 1:ncol(x)) {
    Datam <- data.frame(Valu = as.numeric(x[,i]))
    Datam$Name <- colnames(x)[i]
    Datam$Cells <- factor(rownames(x),levels = rownames(x))
    Data[[i]] <- Datam
    rm(Datam)}
  Data <- na.omit(do.call(rbind,Data))
  Box <- ggplot(Data,aes(Cells,Valu))+geom_boxplot(aes(fill = Cells))+geom_point(aes(color = Name),size = Size)+
    xlab("Cell Type")+ylab("Cell Composition")+scale_x_discrete(labels = substr(unique(Data$Cells),1,6))+theme_light()
  return(Box)}
#
## 8a03a29901b31176e32928321b1349e6 ##

## MK_cbind 8a03a29901b31176e32928321b1349e6
#
MK_cbind_i <- function(F1, F2){
  comfr = intersect(rownames(F1), rownames(F2))
  F1 = F1[comfr,]
  F2 = F2[comfr,]
  rm(comfr)
  return(cbind(F1, F2))}
#
MK_cbind_s <- function(F1, F2){

  ## Miss empty ##
  if(any(dim(F1) == 0)){return(F2)}
  if(any(dim(F2) == 0)){return(F1)}

  ## Each set ##
  rowall = c(rownames(F1), rownames(F2))
  if(length(setdiff(rowall, rownames(F1))) > 0){
    SF1r = matrix(0, nrow = length(setdiff(rowall, rownames(F1))), ncol = ncol(F1))
    rownames(SF1r) = setdiff(rowall, rownames(F1))
    colnames(SF1r) = colnames(F1)
    F1 = rbind(F1, SF1r)
    rm(SF1r)}
  if(length(setdiff(rowall, rownames(F2))) > 0){
    SF2r = matrix(0, nrow = length(setdiff(rowall, rownames(F2))), ncol = ncol(F2))
    rownames(SF2r) = setdiff(rowall, rownames(F2))
    colnames(SF2r) = colnames(F2)
    F2 = rbind(F2, SF2r)
    rm(SF2r)}
  F2 = F2[rownames(F1),]
  return(cbind(F1, F2))}
#
## 8a03a29901b31176e32928321b1349e6

## MK_fix 8a03a29901b31176e32928321b1349e6
#
MK_fix <- function(x,Fix){
  Fname = setdiff(rownames(Fix),rownames(x))
  fr = matrix(0,ncol = ncol(x),nrow = length(Fname),dimnames = list(Fname,colnames(x)))
  fr = rbind(x,fr)
  return(fr)}
#
## 8a03a29901b31176e32928321b1349e6

## MK_toMM 8a03a29901b31176e32928321b1349e6
#
MK_toMM <- function(x, HK_bm = F, Mito_rm = T, AC_rm = T, RP_rm = T, RPLS_rm = T, MIR_rm = T, ATP_rm = T, IGKV_rm = T, verbose = F, name = "temp"){

  # Change sign #
  if(verbose){print(grep("\\.", rownames(x), value = T)[1:6])}
  rownames(x) = gsub("\\.", "-", rownames(x))
  rownames(x) = gsub("_", "-", rownames(x))
  
  # Same gene #
  if(sum(rownames(x) %in% c("IGJ", "JCHAIN")) > 1){
    a = which(rownames(x) %in% c("IGJ", "JCHAIN"))
    JCHAIN = apply(x[a,], 2, sum)
    x = rbind(x[-a,], JCHAIN)
    rm(a, JCHAIN)
  }

  # Rm MT #
  if(Mito_rm){
  
  # MT-, MTATP, MTRNR, MTND, MTCO, MTCYB, MTT #
    if(verbose){
      message("Removing MT-; MTATP; MTRNR; MTND; MTCO; MTCYB; MTT ...", MK_time())
      print(grep("^MT-", rownames(x), value = T)[1:6])
      print(grep("^MTATP", rownames(x), value = T)[1:6])
      print(grep("^MTRNR", rownames(x), value = T)[1:6])
      print(grep("^MTND", rownames(x), value = T)[1:6])
      print(grep("^MTCO", rownames(x), value = T)[1:6])
      print(grep("^MTCYB", rownames(x), value = T)[1:6])
      print(grep("^MTT", rownames(x), value = T)[1:6])
    }
  x = x[!grepl("^MT-", rownames(x)),]
  x = x[!grepl("^MTATP", rownames(x)),]
  x = x[!grepl("^MTRNR", rownames(x)),]
  x = x[!grepl("^MTND", rownames(x)),]
  x = x[!grepl("^MTCO", rownames(x)),]
  x = x[!grepl("^MTCYB", rownames(x)),]
  x = x[!grepl("^MTT", rownames(x)),]
  }

  # Rm AC #
  if(AC_rm){
    if(verbose){
      message("Removing AX; LINC; YRNA; LOC; CTX-; XX-; -AS -IT -OT ...", MK_time())
      print(grep("^A[A-Z][0-9][0-9]", rownames(x), value = T)[1:6])
      print(grep("^LINC[0-9]", rownames(x), value = T)[1:6])
      print(rownames(x)[grepl("^Y-RNA", rownames(x)) | grepl("^Y_RNA", rownames(x))][1:6])
      print(grep("^LOC[0-9]", rownames(x), value = T)[1:6])
      print(grep("^CT[A-Z]-", rownames(x), value = T)[1:6])
      print(rownames(x)[grepl("^XX-", rownames(x)) | grepl("^XX[a-z]", rownames(x))][1:6])
      print(grep("-AS", rownames(x), value = T)[1:6])
      print(grep("-IT", rownames(x), value = T)[1:6])
      print(grep("-OT", rownames(x), value = T)[1:6])
    }
    x = x[!grepl("^A[A-Z][0-9][0-9]", rownames(x)),]
    x = x[!grepl("^LINC[0-9]", rownames(x)),]
    x = x[!grepl("^LOC[0-9]", rownames(x)),]
    x = x[!grepl("^CT[A-Z]-", rownames(x)),]
    x = x[!(grepl("^Y-RNA", rownames(x)) | grepl("^Y_RNA", rownames(x))),]
    x = x[!(grepl("^XX-", rownames(x)) | grepl("^XX[a-z]", rownames(x))),]
    x = x[!grepl("-AS", rownames(x)),]
    x = x[!grepl("-IT", rownames(x)),]
    x = x[!grepl("-OT", rownames(x)),]
  }

  # Rm RP #
  if(RP_rm){
    if(verbose){
      message("Removing RP-; orf; APXXX; ENSGX ...", MK_time())
      print(rownames(x)[grepl("^RP", rownames(x)) & grepl("-", rownames(x))][1:6])
      print(grep("orf[0-9]", rownames(x), value = T)[1:6])
      print(grep("^AP[0-9][0-9][0-9]", rownames(x), value = T)[1:6])
      print(grep("^ENSG[0-9]", rownames(x), value = T)[1:6])
    }
    x = x[!(grepl("^RP", rownames(x)) & grepl("-", rownames(x))),]
    x = x[!grepl("orf[0-9]", rownames(x)),]
    x = x[!grepl("^AP[0-9][0-9][0-9]", rownames(x)),]
    x = x[!grepl("^ENSG[0-9]", rownames(x)),]
  }

  # Rm RPL and RPS #
  if(RPLS_rm){
    message("Removing RPL; RPS ...", MK_time())
    if(verbose){print(rownames(x)[grepl("^RPL", rownames(x)) | grepl("^RPS", rownames(x))][1:6])}
    x = x[!(grepl("^RPL", rownames(x)) | grepl("^RPS", rownames(x))),]
  }
  
  # Rm MIR #
  if(MIR_rm){
    message("Removing MIR ...", MK_time())
    if(verbose){print(rownames(x)[grepl("^MIR[0-9]", rownames(x))][1:6])}
    x = x[!grepl("^MIR[0-9]", rownames(x)),]
  }
  
  # Rm ATP #
  if(ATP_rm){
    message("Removing ATP ...", MK_time())
    if(verbose){print(rownames(x)[grepl("^ATP", rownames(x))][1:6])}
    x = x[!grepl("^ATP", rownames(x)),]
  }
  
  # Rm IGXV #
  if(IGKV_rm){
    message("Removing IGXVxx-x ...", MK_time())
    if(verbose){
      print(rownames(x)[grepl("^IG[A-Z]V", rownames(x)) & grepl("-", rownames(x))][1:6])
    }
    x = x[!(grepl("^IG[A-Z]V", rownames(x)) & grepl("-", rownames(x))),]
  }
  
  # Rm 0 #
  x = x[Matrix::rowSums(x) != 0,]

  # HK batch remove #
  if(HK_bm){
    if(verbose){print(x[1:7, 1:7])}

    # Cell library #
    lib <- apply(x, 2, sum)

    # HK genes #
    HKs <- c("CFL1","RPS15A","RPL38","RPL28","FAU","RPL35A","RPL39","RPL23","RPLP2","RPS23")
    HKbm <- x[HKs,]

    # Predict with library #
    HKbm <- sweep(HKbm,2,lib,FUN = "/")
    for (i in 1:nrow(HKbm)) {
      if(verbose){print(paste(rownames(HKbm)[i],i,sum(as.numeric(HKbm[i,]) == 0),median(as.numeric(HKbm[i,HKbm[i,] != 0]))))}

      # By median #
      HKbm[i,as.numeric(HKbm[i,]) == 0] <- median(as.numeric(HKbm[i,HKbm[i,] != 0]))

      if(verbose){print(paste(i,sum(as.numeric(HKbm[i,]) == 0),median(as.numeric(HKbm[i,HKbm[i,] != 0]))))}}

    # Focus by mean #
    HKbm <- apply(HKbm, 2, mean)

    # Normal matrix #
    x <- sweep(x,2,HKbm/mean(HKbm),FUN = "/")

    if(verbose){print(x[1:7,1:7])}
    rm(HKbm,HKs,lib)
    gc()}

  # Save by dgTMatrix #
  col = colnames(x)
  row = rownames(x)
  if(MKrcpp){
    x = as(MK_asMatr(x), "dgTMatrix")
  }else{
    x = as(as.matrix(x), "dgTMatrix")
  }

  write.csv(col, paste0(name, "_cell.csv"))
  write.csv(row, paste0(name, "_gene.csv"))
  Matrix::writeMM(x, paste0(name, "_mt.mtx"))

  rm(x)
  gc()}
#
## 8a03a29901b31176e32928321b1349e6

## MK_read10x 8a03a29901b31176e32928321b1349e6 ##
#
MK_read10X <- function(MKdir = getwd(), IDin = NULL, Barfile = "barcode", Genefile = "gene", Exprfile = "matrix", View = T, Save = T){

  ## Check dir ##
  if(!dir.exists(MKdir)){
    stop("There is no target dir exist!", MK_time())
  }

  ## IDin ##
  if(is.null(IDin)){
    IDin = gsub("_.*", "", list.files(MKdir))
    IDin = names(table(IDin))[table(IDin) > 2]
  }

  ## IDname ##
  IDname = c(Barfile, Genefile, Exprfile)

  ## Change wd ##
  Odir = getwd()
  setwd(MKdir)
  if(View){message("MK read10X! Now we working on: ", MKdir, MK_time())}

  ## Check File names ##
  File_name = list()
  for (i in 1:3) {
    File_name[[i]] = grep(IDname[i], list.files(), value = T)
    if(!length(File_name[[i]]) > 0){stop("No files matched by your ", IDname[i], " name !", MK_time())}
  }

  ## Check IDin ##
  File_data = list()
  for (i in 1:length(IDin)) {

    ## In each ID ##
    if(View){message(i, " Now we process: ", IDin[i], MK_time())}

    ## Check files ##
    File_nameI = list()
    for (j in 1:3) {
      File_nameI[[j]] = grep(IDin[i], File_name[[j]], value = T)
      if(length(File_nameI[[j]]) != 1){stop("Your < ", IDin[i], " > in ID is not unique or cannot match for files !", MK_time())}
    }

    ## Read Barcode ##
    Barcode = readLines(File_nameI[[1]])

    ## Read Feature ##
    Gene = read.delim(File_nameI[[2]], header = F)
    Ig = NULL
    if(ncol(Gene) > 2){
      Ig = Gene$V3
      if(Save){write.csv(Ig, paste0(IDin[i], MK_time(), "_Ig.csv"))}
    }

    ## Read Expr ##
    Expr = Matrix::readMM(File_nameI[[3]])

    ## Process ##
    colnames(Expr) = paste0(IDin[i], "_", Barcode)
    rownames(Expr) = make.unique(Gene$V2)

    ## Save ##
    File_data[[i]] = Expr

    ## Clean ##
    rm(File_nameI, Barcode, Gene, Ig, Expr)
    gc()
  }

  ## Back wd ##
  setwd(Odir)
  rm(Odir)

  return(File_data)
}
#
## 8a03a29901b31176e32928321b1349e6 ##

## MK_scRNA 8a03a29901b31176e32928321b1349e6 ##
#
MK_scRNA <- function(x, name = NULL, Reso = 0.6, nGene = c(200, Inf), nCount = c(200, Inf), nVar = 3, Dim = 2, SCT = F, BatchRemove = F, Umap = F, Plot = T, Norm = T, save = T){
  if(!any(installed.packages() %in% "Seurat")){
    install.packages("Seurat")
  }
  suppressMessages(library(Seurat))

  if(is.null(name)){name = "temp"}
  name = as.character(name)
  
  ## Creat Seurat v3.2 ##
  x = CreateSeuratObject(x, name, min.features = nGene[1])
  if(Plot){
    message("Matrix dim: ", paste(dim(x), collapse = " "), MK_time())
    print(VlnPlot(x, c("nFeature_RNA", "nCount_RNA"), ncol = 2, pt.size = 0.2))
    message("Input nGene Min-cutoff: ", MK_time())
    nGene[1] = scan()
    message("Input nGene Max-cutoff: ", MK_time())
    nGene[2] = scan()
    message("Nice Your nGene cut off: ", paste(nGene, collapse = " "), MK_time())
    message("Input nCount Min-cutoff: ", MK_time())
    nCount[1] = scan()
    message("Input nCount Max-cutoff: ", MK_time())
    nCount[2] = scan()
    message("Nice Your nCount cut off: ", paste(nCount, collapse = " "), MK_time())
  }
  x = x[, x$nFeature_RNA >= nGene[1] & x$nFeature_RNA <= nGene[2] & x$nCount_RNA >= nCount[1] & x$nCount_RNA <= nCount[2]]
  if(Plot){message("Matrix dim: ", paste(dim(x), collapse = " "), MK_time())}

  ## If SCTransform ##
  if(SCT){
    x = SCTransform(x, verbose = Plot, variable.features.n = 1000*nVar)

  }else{
    ## Standard process ##
    ## Normalize ##
    if(Norm){x = NormalizeData(x, verbose = Plot)}

    ## Find varible ##
    x = FindVariableFeatures(x, nfeatures = 1000*nVar)
    if(Plot){print(LabelPoints(VariableFeaturePlot(x), points = head(VariableFeatures(x), 10), repel = T))}

    ## Scale data ##
    x = ScaleData(x, features = rownames(x), verbose = Plot)
  }

  ## Cluster ##
  x = RunPCA(x, features = VariableFeatures(x), verbose = F)
  x = x[, !duplicated(x@reductions$pca@cell.embeddings)]

  if(BatchRemove){

    ## Harmony v1.0 ##
    if(!any(installed.packages() %in% "harmony")){
      if(!any(installed.packages() %in% "devtools")){
        install.packages("devtools")
      }
      devtools::install_github("immunogenomics/harmony")
    }
    suppressMessages(library(harmony))
    AData = ifelse(SCT, "SCT", "RNA")
    x = RunHarmony(x, "orig.ident", max.iter.harmony = 20, max.iter.cluster = 50, assay.use = AData, verbose = Plot, plot_convergence = Plot)
  }

  ## TSNE and UMAP ##
  ARedu = ifelse(BatchRemove, "harmony", "pca")
  x = RunTSNE(x, dims = 1:50, dim.embed = Dim, reduction = ARedu)
  if(Umap){x = RunUMAP(x, dims = 1:50, n.components = Dim, reduction = ARedu, verbose = Plot)}

  ## SNN ##
  x = FindNeighbors(x, dims = 1:50, reduction = ARedu, verbose = Plot)
  x = FindClusters(x, resolution = Reso, verbose = Plot)

  if(Plot){print(DimPlot(x, label = T, pt.size = 0.9))}

  ## Save ##
  if(save){
    dir.create("backup")
    saveRDS(x, paste0("backup/", name, MK_time(), "_backup.rds"))}
  gc()
  message("MK_scRNA process well done !!!", MK_time())

  return(x)
}
#
## 8a03a29901b31176e32928321b1349e6 ##

## MK_clear 8a03a29901b31176e32928321b1349e6 ##
#
MK_clear <- function(x,name = "temp",Save = T){
  library(Seurat)

  ## Main data ##
  DataB <- Matrix::rowSums(x@assays$RNA@counts)

  ## Clear ##
  x@assays$RNA@counts <- matrix()
  x@assays$RNA@scale.data <- matrix()
  x@graphs <- list()
  x@reductions$pca <- list()
  x@reductions$harmony <- list()
  x@meta.data <- data.frame()
  x@commands <- list()
  x@assays$RNA@var.features <- list()
  x@assays$RNA@meta.features <- data.frame()
  if(Save){
    saveRDS(x,paste0(name,"_clear.rds"))
    write.csv(DataB,paste0(name,"_bulk.csv"))
  }
  rm(DataB)
  gc()
  return(x)}
#
## 8a03a29901b31176e32928321b1349e6 ##

## MK_Enrich 8a03a29901b31176e32928321b1349e6 ##
#
MK_Enrich <- function(x, EnID = "temp", CutP = 0.01, Save = T, Wid = 8, Hig = 8.3){
  
  if(!any(installed.packages() %in% "clusterProfiler")){
    if(!any(installed.packages() %in% "BiocManager")){
      install.packages("BiocManager")
    }
    BiocManager::install("clusterProfiler")
  }
  if(!any(installed.packages() %in% "ReactomePA")){
    if(!any(installed.packages() %in% "BiocManager")){
      install.packages("BiocManager")
    }
    BiocManager::install("ReactomePA")
  }
  if(!any(installed.packages() %in% "ggplot2")){
    install.packages("ggplot2")
  }
  
  suppressMessages(library(clusterProfiler))
  suppressMessages(library(ReactomePA))
  suppressMessages(library(ggplot2))

  ## set enrich ##
  path = getwd()
  dir.create("Enrich")
  setwd("Enrich")

  ## wdir gene-id ##
  EnID = as.character(EnID)
  dir.create(EnID)
  setwd(EnID)
  message("Now working in ", getwd(), MK_time())
  GeneID = bitr(gsub("^MT\\.","MT-",x), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  if(Save){write.csv(GeneID, paste0(EnID, MK_time(), "_GeneID.csv"))}

  ## GO ##
  GoRE = data.frame(tryCatch(
    enrichGO(gene = as.character(GeneID[,2]), "org.Hs.eg.db", ont = "ALL", pAdjustMethod = "BH", minGSSize = 3, pvalueCutoff = CutP, readable = TRUE),
    error = function(e) data.frame()))
  if(Save){write.csv(GoRE, paste0(EnID, MK_time(), "_GoRE.csv"))}
  if(nrow(GoRE) > 0){
    GoRE = na.omit(GoRE[c(which(GoRE$ONTOLOGY == "BP")[1:10], which(GoRE$ONTOLOGY == "CC")[1:10], which(GoRE$ONTOLOGY == "MF")[1:10]),])
    GoRE = GoRE[order(GoRE$Count),]

    ## plot ##
    if(!any(installed.packages() %in% "stringr")){
      install.packages("stringr")
    }
    suppressMessages(library(stringr))
    GoPlot = ggplot(GoRE,
                     aes(GoRE$Count/length(as.character(GeneID[,2])), factor(GoRE$Description, levels = GoRE$Description)))+
      geom_point(aes(size = GoRE$Count, color = -1*log10(GoRE$qvalue), shape = GoRE$ONTOLOGY))+
      scale_colour_gradient(low = "blue", high = "red")+
      labs(color = expression(-log[10](Qvalue)), size = "Gene number", shape = "Ontology", x = "GeneRatio", y = NULL, title = "GO enrichment")+
      theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text.y = element_text(size = 12))+
      scale_y_discrete(labels = function(x) str_wrap(x, width = 50))
    if(Save){ggsave(paste0(EnID, MK_time(), ".tiff"), plot = GoPlot, device = "tiff", width = Wid, height = Hig)}
    rm(GoPlot)
  }
  gc()

  ## KEGG ##
  KeRE = data.frame(tryCatch(
    enrichKEGG(gene = as.character(GeneID[,2]), organism = "hsa", pvalueCutoff = CutP),
    error = function(e) data.frame()))
  if(nrow(KeRE) > 0){
    Ncol = ncol(KeRE)
    for (i in 1:nrow(KeRE)) {
      KeRE[i,Ncol+1] = paste(as.character(GeneID[,1])[GeneID[,2] %in% strsplit(KeRE$geneID,"/")[[i]]],collapse = "/")
    }
    rm(Ncol,i)

    ## plot ##
    if(Save){
      if(!any(installed.packages() %in% "pathview")){
        if(!any(installed.packages() %in% "BiocManager")){
          install.packages("BiocManager")
        }
        BiocManager::install("pathview")
      }
      suppressMessages(library(pathview))
      tryCatch(pathview(gene.data = as.character(GeneID[,2]), pathway.id = KeRE$ID, species = "hsa"))
      rm(bods, gene.idtype.bods, gene.idtype.list, cpd.simtypes, korg)
    }
  }
  if(Save){write.csv(KeRE, paste0(EnID, MK_time(), "_KeRE.csv"))}
  gc()

  ## ReactPA ##
  PaRE = data.frame(tryCatch(
    enrichPathway(gene = as.character(GeneID[,2]), pvalueCutoff = CutP, organism = "human", readable = T)@result,
    error = function(e) data.frame()))
  if(Save){write.csv(PaRE, paste0(EnID, MK_time(), "_PaRE.csv"))}
  gc()

  ## Return ##
  setwd(path)
  return(list(GoRE, KeRE, PaRE))}
#
## 8a03a29901b31176e32928321b1349e6 ##

## MK_toMMs 8a03a29901b31176e32928321b1349e6 ##
#
MK_toMMs <- function(x, name = "temp", Cells = 10, verbose = T, HK_bm = F, Mito_rm = T, AC_rm = T, RP_rm = T, RPLS_rm = T){
  name = as.character(name)
  dir.create(name)

  i = 1
  ## more 10k*Cells ##
  if(ncol(x) >= (Cells*1000)){
    for (i in 1:floor(ncol(x)/(Cells*1000))) {
      message("Save MM ", i, MK_time())
      MK_toMM(x[, ((i-1)*(Cells*1000)+1):(i*(Cells*1000))], name = paste0(name, "/", name, " ", i), verbose = verbose,
              HK_bm = HK_bm, Mito_rm = Mito_rm, AC_rm = AC_rm, RP_rm = RP_rm, RPLS_rm = RPLS_rm)
    }
    i = i + 1
  }

  ## remain ##
  if(ncol(x) %% (Cells*1000) != 0){
    message("Save MM ex ", i, MK_time())
    MK_toMM(x[, ((i-1)*(Cells*1000)+1):ncol(x)], name = paste0(name, "/", name, " ", i), verbose = verbose,
            HK_bm = HK_bm, Mito_rm = Mito_rm, AC_rm = AC_rm, RP_rm = RP_rm, RPLS_rm = RPLS_rm)
  }
  message("MK_toMMs done !!!", MK_time())
}
#
## 8a03a29901b31176e32928321b1349e6 ##

## MK_TCGA 8a03a29901b31176e32928321b1349e6 ##
MK_Tcga_CanerFilt <- function(x,cancer = T,save = F,name = "temp"){
  if(sum(duplicated(substr(colnames(x),1,12)))>0){
    x <- x[,colnames(x) %in% sort(colnames(x))[!duplicated(substr(sort(colnames(x)),1,12))]]
  }
  Can<-x[,which(substr(colnames(x),14,14) == "0")]
  Nor<-x[,which(substr(colnames(x),14,14) == "1")]
  colnames(Can) <- gsub("\\.","_",substr(toupper(colnames(Can)),1,12))
  if(save){
    write.csv(Can,paste0(name,"_Tcga_Cancer.csv"))
    write.csv(Nor,paste0(name,"_Tcga_Normal.csv"))
  }
  if(cancer){
    rm(Nor)
    return(Can)
  }else{
    rm(Can)
    return(Nor)
  }
}
#
MK_DEGs <- function(x, y, filt = T, log2FC = 2, padj = 0.01, pval = 0.01, save = T, Order = T, name = "temp"){
  options(stringsAsFactors = F)
  if(!any(installed.packages() %in% "limma")){
    if(!any(installed.packages() %in% "BiocManager")){
      install.packages("BiocManager")
    }
    BiocManager::install("limma")
  }
  suppressMessages(library(limma))
  Data = cbind(x,y)
  Group = data.frame(row.names = colnames(Data), Group1 = c(rep(1, ncol(x)), rep(0, ncol(y))), Group2 = c(rep(0, ncol(x)), rep(1, ncol(y))))
  Sig = makeContrasts("Group1-Group2", levels = Group)
  fit = eBayes(contrasts.fit(lmFit(Data, Group), Sig))
  OP = na.omit(topTable(fit, number = nrow(Data), coef = 1, adjust = "BH"))
  rm(Data, Group, Sig, fit)
  colnames(OP) = c("LogFC","AveExpr","T","P_val","P_adj","B")
  OP$MeanP = apply(x,1,mean)
  OP$MedianP = apply(x,1,median)
  OP$MeanN = apply(y,1,mean)
  OP$MedianN = apply(y,1,median)
  if(filt){OP = subset(OP, P_adj < padj & abs(LogFC) > log2FC)}
  OP$Sig = "NO"
  OP$Sig[OP$P_val < 0.05] = "Yes"
  OP$State = "None"
  OP$State[OP$LogFC > 0] = "Up"
  OP$State[OP$LogFC < 0] = "Down"
  if(Order){OP = OP[order(abs(OP$LogFC), decreasing = T),]}
  if(save){write.csv(OP, paste0(name, MK_time(), "_DEGs.csv"))}
  return(OP)
}
#
MK_DEGplot <- function(x,Title = "temp",pvalue = 0.01,log2FC = 2,plimit = 30,log2limit = 5,color = 3,Lima = F,adj = T){
  library(ggplot2)
  if(Lima){
    if(adj){colnames(x) <- c("log2FoldChange","AveExpr","t","p","padj","B")}
    else{colnames(x) <- c("log2FoldChange","AveExpr","t","padj","p","B")}}
  x$Legend <- as.factor(ifelse(x$padj < pvalue & abs(x$log2FoldChange) >= log2FC, ifelse(x$log2FoldChange > log2FC,'Up','Down'),'Not'))
  if(color == 3){colornum <- c("blue", "black", "red")}
  if(color == 2){colornum <- c("black", "red")}
  DEp <- ggplot(data = x,aes(x = log2FoldChange, y = -log10(padj),colour = Legend))+ggtitle(Title)+
    xlab("log2 Foldchange")+ylab("-log10 Padj")+geom_vline(xintercept = c(-log2FC,log2FC),lty = 6,col = "grey",lwd = 0.5)+
    geom_hline(yintercept = -log10(pvalue),lty = 4,col = "grey",lwd = 0.5)+scale_color_manual(values = colornum)+theme(legend.position = "right")+theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.title = element_blank())+xlim(-log2limit,log2limit)+ylim(0,plimit)+
    theme(plot.title = element_text(hjust = 0.5))+geom_point(alpha=0.4, size=1.2)
  return(DEp)
}
##

## MK_rem0 8a03a29901b31176e32928321b1349e6 ##
#
MK_rem0 <- function(x, Rem0 = 0.1, raito = T, MKrcpp = F){
  if(MKrcpp){
    Mat = MK_asMatr(x)
  }else{
    Mat = as.matrix(x)
  }
  Mat = apply(Mat, 2, as.numeric)
  r = apply(Mat, 1, function(i) sum(i == 0))
  if(raito){ok = which(r <= dim(Mat)[2]*Rem0)
  }else{
    ok = which(r <= Rem0)
  }
  rm(r)
  x = x[ok,]
  rm(ok)
  return(x)
}
#
## 8a03a29901b31176e32928321b1349e6 ##

## MK_WGCNA 8a03a29901b31176e32928321b1349e6 ##
#
MK_WG_Tom <- function(x, name = "temp", nGene = 10000, Save = T){
  if(!any(installed.packages() %in% "WGCNA")){
    install.packages("WGCNA")
  }
  suppressMessages(library(WGCNA))
  if(nGene > 11000){stop("Sorry Memory Size Crush when nGene over 11000.")}

  ## Filt virable genes ##
  if(nGene < nrow(x)){
    x = t(x[order(apply(x, 1, mad), decreasing = T)[1:nGene],])
  }else{
    x = t(x)
  }

  ## Check genes and samples ##
  gsg = goodSamplesGenes(x)
  x = x[gsg$goodSamples, gsg$goodGenes]
  rm(gsg)

  ## Estimate power ##
  sft = pickSoftThreshold(x, blockSize = 12000)
  okpower = sft$powerEstimate
  rm(sft)
  if(is.na(okpower)){stop("No okpower !")}

  ## Make net ##
  net = blockwiseModules(x, power = okpower, numericLabels = T, saveTOMs = T, saveTOMFileBase = paste(name, "WGTOM"), maxBlockSize = 12000)
  moduleColors = labels2colors(net$colors)
  geneTree = net$dendrograms[[1]]
  MEs = orderMEs(moduleEigengenes(x, moduleColors)$eigengenes)
  rm(net)
  load(grep(paste(name, "WGTOM"), list.files(), value = T))
  gc()

  ## Return ##
  RE = list(x, MEs, moduleColors, okpower, TOM)

  if(Save){
    dir.create("backup")
    saveRDS(RE,paste0("backup/", name, MK_time(), "_WGtom_backup.rds"))

    ## Plot TOM ##
    plotTOM = -(1 - as.matrix(TOM)) ^ 7
    diag(plotTOM) = NA
    sizeGrWindow(9,9)
    tiff(paste(name,"TOM.tiff"), width = 800, height = 850, pointsize = 30, compression = "lzw")
    TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap of all modules")
    dev.off()
    rm(plotTOM)
    gc()
  }

  rm(MEs, moduleColors, okpower, geneTree)
  return(RE)
}
#
MK_WG_CliIn <- function(Clin,WG_Tom,name = "temp",color = "ALL",classname = "ALL",Cys = T){
  if(!any(installed.packages() %in% "WGCNA")){
    install.packages("WGCNA")
  }
  suppressMessages(library(WGCNA))
  if(any(color == "ALL")){color = substring(names(WG_Tom[[2]]), 3)}
  if(any(classname == "ALL")){classname = colnames(Clin)}

  ## MM and GS ##
  weight <- data.frame(Clin[,classname])
  names(weight) <- classname
  modNames <- substring(names(WG_Tom[[2]]), 3)
  MM <- cor(WG_Tom[[1]], WG_Tom[[2]], use = "p")
  MM_P <- corPvalueStudent(MM, nrow(Clin))
  colnames(MM)<- paste("MM", modNames)
  colnames(MM_P) <- paste("P-MM", modNames)
  GS <- cor(WG_Tom[[1]], weight, use = "p")
  GS_P <- corPvalueStudent(as.matrix(GS), nrow(Clin))
  colnames(GS) <- paste("GS", colnames(weight))
  colnames(GS_P) <- paste("P-GS", colnames(weight))

  ## Save ##
  write.csv(MM,paste(name,"WG_MM.csv"))
  write.csv(MM_P,paste(name,"WG_MM_P.csv"))
  write.csv(GS,paste(name,"WG_GS.csv"))
  write.csv(GS_P,paste(name,"WG_GS_P.csv"))
  rm(weight,MM_P,GS_P)

  ## MMGS Plot ##
  ModGenes <- list()
  for (i in 1:length(color)) {
    column <- match(color[i], modNames)
    for (j in 1:length(classname)) {
      tiff(paste(name,color[i],classname[j],"WG_MMGS.tiff"),width = 700,height = 550,pointsize = 15,compression = "lzw")
      verboseScatterplot(abs(MM[WG_Tom[[3]] == color[i], column]), abs(GS[WG_Tom[[3]] == color[i], j]),
                         xlab = paste("Module Membership in", color[i], "module"), ylab = paste("Gene significance for", classname[j]),
                         abline = 1,abline.lty = 1,abline.color = "red",main = paste("Module membership vs. Gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")
      dev.off()
    }

    ## Cys ##
    inModule <- is.finite(match(WG_Tom[[3]], color[i]))
    modGenes <- colnames(WG_Tom[[1]])[inModule]
    if(Cys){
      modTOM <- as.matrix(WG_Tom[[5]])[inModule, inModule]
      dimnames(modTOM) <- list(modGenes,modGenes)
      cys <- exportNetworkToCytoscape(modTOM,threshold = 0.05,edgeFile = paste(name,"edges",color[i],"WG_Cys.txt"),nodeFile = paste(name,"nodes",color[i],"WG_Cys.txt"),
                                      nodeNames = modGenes, nodeAttr = WG_Tom[[3]][inModule])
      modGenes <- cys$nodeData$nodeName
      rm(modTOM,cys)
    }
    ModGenes <- c(ModGenes,list(modGenes))
    rm(inModule,column,modGenes)
  }

  ## Return ##
  rm(modNames)
  return(ModGenes)
}
#
## 8a03a29901b31176e32928321b1349e6 ##

## MK_asNum 8a03a29901b31176e32928321b1349e6 ##
#
MK_asNum = function(x, nsep = 5, verbose = T){

  ## save names ##
  coln = colnames(x)
  rown = rownames(x)

  ## check ncol ##
  if(ncol(x) < nsep*3){

    x = apply(x, 2, as.numeric)
    colnames(x) = coln

    return(x)
  }

  ## step ##
  sp = floor(ncol(x)/nsep)

  ## seprate ##
  Re = list()

  for (i in 1:nsep) {
    a = (i-1) * sp + 1
    z = i * sp
    mt = apply(x[,a:z], 2, as.numeric)
    Re[[i]] = as(as.matrix(mt), "dgCMatrix")
    if(verbose){message("Num: ", a, " to ", z, MK_time())}
    rm(a,z,mt)
  }
  i = i + 1

  ## remain ##
  if(ncol(x) %% sp != 0){
    a = (i-1) * sp + 1
    mt = apply(x[,a:ncol(x)], 2, as.numeric)
    Re[[i]] = as(as.matrix(mt), "dgCMatrix")
    if(verbose){message("Num: ex ", a, " to ", ncol(x), MK_time())}
    rm(a,mt)
  }

  ## cbind ##
  Re = do.call(cbind, Re)
  colnames(Re) = coln
  rownames(Re) = rown
  rm(sp,coln,rown)

  return(Re)
}
#
## 8a03a29901b31176e32928321b1349e6 ##

## MK_Exct 8a03a29901b31176e32928321b1349e6 ##
#
MK_Exct <- function(x, filed = 1, exct = "\\|", verbose = F){
  if(verbose){message("Before: ", x[1], MK_time())}
  x <- sapply(strsplit(x, exct), function(i) i[[filed]])
  if(verbose){message("After: ", x[1], MK_time())}
  return(x)
}
#
## 8a03a29901b31176e32928321b1349e6 ##

## MK_Download 8a03a29901b31176e32928321b1349e6 ##
#
MK_Download <- function(x, name = NULL, sleep = 2, outdir = getwd()){

  # check name #
  if(is.null(name)){
    name = 1:length(x)
  }

  for (i in 1:length(x)) {

    # set index #
    write.csv("ok", "Erro.csv", row.names = F)

    # check file #
    if(any(list.files(outdir) %in% name[i])){
      message(i, " Already here for ", as.character(name)[i], MK_time())
    }else{
      tryCatch(download.file(as.character(x[i]),
                             destfile = paste0(outdir, "/", as.character(name)[i]),
                             quiet = T, cacheOK = T, method = "libcurl"),
               error = function(e){
                 write.csv("MikuGene download erro", "Erro.csv", row.names = F)
                 message(i," Download Failed for ", as.character(name)[i], MK_time())
               }
      )

      # check index #
      Erro <- as.character(read.csv("Erro.csv", header = T)[1,])
      if(Erro == "ok"){
        message(i," Download OK for ", as.character(name)[i], MK_time())
      }
      while (Erro == "MikuGene download erro"){
        write.csv("ok", "Erro.csv", row.names = F)
        tryCatch(download.file(as.character(x[i]),
                               destfile = paste0(outdir, "/", as.character(name)[i]),
                               quiet = T,cacheOK = T,method = "libcurl"),
                 error = function(e){
                   write.csv("MikuGene download erro", "Erro.csv", row.names = F)
                   message(i," Download Failed for ", as.character(name)[i], MK_time())
                 }
        )
        Erro <- as.character(read.csv("Erro.csv", header = T)[1,])
        if(Erro == "ok"){
          message(i," Download OK for ", as.character(name)[i], MK_time())
        }
        Sys.sleep(sleep)}
    }
  }
  Sys.sleep(0.5)
}
#
## 8a03a29901b31176e32928321b1349e6 ##

## MK_Microbiology 8a03a29901b31176e32928321b1349e6 ##
#
MK_BuildVirusRef <- function(version = "2020.3", OutVs = "default", verbose = T){
  options(scipen = 200)

  # OutV #
  if(!is.null(OutVs)){
    if(any(OutVs %in% "default")){
      if(verbose){
        message("Remove some out-virus ...", MK_time())
      }
      OutVs = c("NC_041925", "NC_032111", "NC_022518", "NC_018464",
                "NC_001506", "NC_038858", "NC_043329", "NC_035797",
                "NC_015253", "NC_043329", "NC_043314", "NC_041831")
    }
  }

  # Download from virusite.org #
  if(!any(list.files() %in% "MikuGene_virus.fasta.zip")){
    if(verbose){
      message("Downloading genes.fasta ...", MK_time())
    }
    MK_Download(paste0("http://www.virusite.org/archive/", version, "/genes.fasta.zip"), "MikuGene_virus.fasta.zip")
  }
  if(!any(list.files() %in% "MikuGenome_virus.fasta.zip")){
    if(verbose){
      message("Downloading genomes.fasta ...", MK_time())
    }
    MK_Download(paste0("http://www.virusite.org/archive/", version, "/genomes.fasta.zip"), "MikuGenome_virus.fasta.zip")
  }
  unzip("MikuGene_virus.fasta.zip")
  unzip("MikuGenome_virus.fasta.zip")

  if(!any(installed.packages() %in% "Biostrings")){
    if(!any(installed.packages() %in% "BiocManager")){
      install.packages("BiocManager")
    }
    BiocManager::install("Biostrings")
  }
  suppressMessages(library(Biostrings))

  # Ref from Genome #
  if(verbose){
    message("Processing genomes.fasta ...", MK_time())
  }
  Geno = readBStringSet("genomes.fasta")
  Genome = Geno@ranges@NAMES
  # ID  Length  Name #
  Geno_ID = MK_Exct(Genome, 2)
  # Rem OutVs #
  if(!is.null(OutVs)){
    Geno = Geno[!Geno_ID %in% OutVs]
    Genome = Geno@ranges@NAMES
    Geno_ID = MK_Exct(Genome, 2)
  }
  Geno_siz = MK_Exct(Genome, 3)
  Geno_name = MK_Exct(Genome, 4)
  # Save genome meta ##
  Virus_genome = data.frame(ID = Geno_ID, Length = Geno_siz, Name = Geno_name)
  write.csv(Virus_genome, "Virus_genome.csv")
  # Make GTF ##
  Gtf = data.frame(SeqID = Geno_ID, source = "MikuGenome", feature = "exon",
                   start = 1, end = as.numeric(gsub("nt","",Geno_siz)),
                   sore = ".", strand = "+", frame = ".",
                   attributes = paste0("gene_id ", Geno_ID, ";", "transcript_id ", Geno_ID, ";"))
  write.table(Gtf, "MikuGenome_virus.gtf", row.names = F, col.names = F, sep = "\t", quote = F)
  # Process fasta name #
  Geno@ranges@NAMES = Geno_ID
  writeXStringSet(Geno, "MikuGenome_virus.fasta", append = F)
  rm(Geno) + gc()

  # Ref from Gene #
  if(verbose){
    message("Processing genes.fasta ...", MK_time())
  }
  Gene = readBStringSet("genes.fasta")
  Gtf = Gene@ranges@NAMES
  # ID Location Name #
  Gtf_Attr = MK_Exct(Gtf, 2)
  if(!is.null(OutVs)){
    Gene = Gene[!Gtf_Attr %in% OutVs]
    Gtf = Gene@ranges@NAMES
    Gtf_Attr = MK_Exct(Gtf, 2)
  }
  Gtf_local = MK_Exct(Gtf, 3)
  Gtf_name = MK_Exct(Gtf, 4)
  # Save gene meta #
  Virus_gene = data.frame(ID = Gtf_Attr, Location = Gtf_local, Name = Gtf_name)
  write.csv(Virus_gene, "Virus_gene.csv")
  # Filt char #
  Gtf_local = gsub("complement\\(", "", gsub("\\)", "", Gtf_local))
  # Remove join and >/< and order #
  GetID = which(!(
    grepl("join", Gtf_local) | grepl("<", Gtf_local) |
      grepl(">", Gtf_local) | grepl("order", Gtf_local)
  ))
  # Start and End #
  Gtf_local = Gtf_local[GetID]
  Gtf_start = as.numeric(gsub("\\..*", "", Gtf_local))
  Gtf_end = as.numeric(gsub(".*\\.", "", Gtf_local))
  # Make Gtf attribute #
  Gtf_Attr = Gtf_Attr[GetID]
  Gtf_gene = paste("gene_id", paste0("MK", GetID))
  Gtf_tran = paste("transcript_id", paste0("MK", GetID))
  # Make GTF #
  Gtf = data.frame(SeqID = Gtf_Attr, source = "MikuGene", feature = "exon",
                   start = Gtf_start, end = Gtf_end,
                   sore = ".", strand = "+", frame = ".",
                   attributes = paste0(Gtf_gene, ";", Gtf_tran, ";"))
  write.table(Gtf, "MikuGene_virus.gtf", row.names = F, col.names = F, sep = "\t", quote = F)
  rm(list = ls()) + gc()

  message("MK_virusref build done .", MK_time())
}
##
MK_VirMap <- function(path_r1, path_r2, name = NULL, maxMiss = 3, GTF = T){
  options(scipen = 200)
  if(is.null(name)){
  name = "temp"
  }
  if(!any(installed.packages() %in% "Rsubread")){
    if(!any(installed.packages() %in% "BiocManager")){
      install.packages("BiocManager")
    }
    BiocManager::install("Rsubread")
  }
  suppressMessages(library(Rsubread))

  # Ref check #
  if(!any(grepl("MikuVirusref", list.files()) & !grepl("MikuVirusref.log", list.files()))){

    # build ref-index #
    buildindex("MikuVirusref", "MikuGenome_virus.fasta")
  }

  # Read Geno/Gene meta #
  GenoID = read.csv("Virus_genome.csv", row.names = 1)
  GeneID = read.csv("Virus_gene.csv",row.names = 1)

  # Align #
  align(index = "MikuVirusref",
        readfile1 = path_r1,
        readfile2 = path_r2,
        output_file = paste0(name, ".BAM"),
        nBestLocations = nrow(GenoID),
        maxMismatches = maxMiss,
        nthreads = 6)

  # fea-count #
  if(GTF){
    Fea1 = featureCounts(paste0(name, ".BAM"),
                         annot.ext = "MikuGene_virus.gtf",
                         isGTFAnnotationFile = T,
                         verbose = F)
    Fea2 = featureCounts(paste0(name, ".BAM"),
                         annot.ext = "MikuGenome_virus.gtf",
                         isGTFAnnotationFile = T,
                         verbose = F)

    # process gene #
    Re_Gene = data.frame(Fea1$counts)
    Re_Gene$Name = GeneID$Name[as.numeric(gsub("MK", "", rownames(Re_Gene)))]
    Re_Gene$Local = GeneID$Location[as.numeric(gsub("MK", "", rownames(Re_Gene)))]
    Re_Gene = Re_Gene[as.numeric(Re_Gene[,1]) != 0,]
    if(nrow(Re_Gene) > 0){
      Re_Gene = Re_Gene[order(Re_Gene[,1], decreasing = T),]
    }
    write.csv(Re_Gene, paste(name, "Re_Gene.csv"))

    # process genome #
    Re_Genome = data.frame(Fea2$counts)
    Re_Genome$Name = GenoID$Name[match(rownames(Re_Genome), GenoID$ID)]
    Re_Genome = Re_Genome[as.numeric(Re_Genome[,1]) != 0,]
    if(nrow(Re_Genome) > 0){
      Re_Genome = Re_Genome[order(Re_Genome[,1], decreasing = T),]
    }
    write.csv(Re_Genome, paste(name, "Re_Genome.csv"))
  }
  rm(list = ls())
  message("MK_VirMap done ...", MK_time())
}
#
MK_BuildBacterRef <- function(verbose = T){
  
  if(!any(list.files() %in% "MikuGenome_bacter.fasta")){
    
    # Meta from NCBI #
    MK_Download("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt",
                "MikuGene_BacRef.txt")
    
    # Process #
    if(!any(installed.packages() %in% "tidyverse")){
      install.packages("tidyverse")
    }
    library(tidyverse)
    Down = read.table("MikuGene_BacRef.txt", sep = "\t", fill = T)
    Down = apply(Down, 1, function(i){
      str_extract_all(i, "ftp:.*",simplify = T)
    })
    Down = do.call(c, Down)
    Down = str_remove(Down[Down != ""], "\t.*")
    
    # Download #
    BacID = gsub(".*/", "", Down)
    Bac_Genome_fa = paste0(BacID, "_genomic.fna.gz")
    Bac_Genome_gtf = paste0(BacID, "_genomic.gtf.gz")
    Bac_dowm_fa = paste(Down, Bac_Genome_fa, sep = "/")
    Bac_dowm_gtf = paste(Down, Bac_Genome_gtf, sep = "/")
    dir.create("MikuBacRef")
    if(verbose){
      message("Downloading genomes.fasta ...", MK_time())
    }
    MK_Download(c(Bac_dowm_fa, Bac_dowm_gtf),
                name = c(Bac_Genome_fa, Bac_Genome_gtf),
                outdir = "MikuBacRef")
    
    # Ref from Genome #
    if(verbose){
      message("Processing genomes.fasta ...", MK_time())
    }
    if(!any(installed.packages() %in% "Biostrings")){
      if(!any(installed.packages() %in% "BiocManager")){
        install.packages("BiocManager")
      }
      BiocManager::install("Biostrings")
    }
    suppressMessages(library(Biostrings))
    Fa = grep("fna.gz", list.files("MikuBacRef"), value = T)
    Fa_info = list()
    for (i in 1:floor(length(Fa)/1000)) {
      a = (i-1)*1000 +1
      z = i*1000
      if(verbose){
        message("Reading Fasta From ", a, " to ", z)
      }
      Fna = readBStringSet(paste0("MikuBacRef/", Fa[a:z]))
      Fna = Fna[Fna@ranges@width < 1e+7]
      FaID = MK_Exct(Fna@ranges@NAMES, exct = " ")
      Fa_geno = data.frame(ID = FaID, Length = Fna@ranges@width, Name = Fna@ranges@NAMES)
      Fa_info[[i]] = Fa_geno
      Fna@ranges@NAMES = FaID
      writeXStringSet(Fna, "MikuGenome_bacter.fasta", append = T)
      rm(Fna, a, z, FaID) + gc()
    }
    i = i +1
    if(length(Fa) %% 1000 != 0){
      a = (i-1)*1000 +1
      if(verbose){
        message("Reading Fasta From ", a, " to ", length(Fa))
      }
      Fna = readBStringSet(paste0("MikuBacRef/", Fa[a:length(Fa)]))
      Fna = Fna[Fna@ranges@width < 1e+7]
      FaID = MK_Exct(Fna@ranges@NAMES, exct = " ")
      Fa_geno = data.frame(ID = FaID, Length = Fna@ranges@width, Name = Fna@ranges@NAMES)
      Fa_info[[i]] = Fa_geno
      Fna@ranges@NAMES = FaID
      writeXStringSet(Fna, "MikuGenome_bacter.fasta", append = T)
      rm(Fna, a, FaID) + gc()
    }
    Fa_info = do.call(rbind, Fa_info)
    write.csv(Fa_info, "Bacter_genome.csv")
    
    # Make GTF ##
    Gtf = data.frame(SeqID = Fa_info$ID, source = "MikuGenome", feature = "exon",
                     start = 1, end = Fa_info$Length,
                     sore = ".", strand = "+", frame = ".",
                     attributes = paste0("gene_id ", Fa_info$ID, ";", "transcript_id ", Fa_info$ID, ";"))
    write.table(Gtf, "MikuGenome_bacter.gtf", row.names = F, col.names = F, sep = "\t", quote = F)
    
    # Ref from Gene #
    if(verbose){
      message("Processing gtfs (waiting this function)...", MK_time())
    }
    rm(list = ls())
  }
}
#
## 8a03a29901b31176e32928321b1349e6 ##

## MK_asMatr 8a03a29901b31176e32928321b1349e6 ##
#
if(MKrcpp){
  Rcpp::sourceCpp(code = '
                  #include <Rcpp.h>
                  using namespace Rcpp;
                  //[[Rcpp::export]]
                  IntegerMatrix MKrc_asmat(NumericVector rp, NumericVector cp, NumericVector z,
                  int nrows, int ncols){
                  int k = z.size() ;
                  IntegerMatrix  mat(nrows, ncols);

                  for (int i = 0; i < k; i++){
                  mat(rp[i],cp[i]) = z[i];
                  }
                  return mat;
                  }
                  ')

  MK_asMatr <- function(mat){

    ## extract ##
    rowP <- mat@i
    colP <- findInterval(seq(mat@x)-1, mat@p[-1])

    Mat <- MKrc_asmat(rp = rowP, cp = colP, z = mat@x,
                      nrows =  mat@Dim[1], ncols = mat@Dim[2])

    ## return ##
    row.names(Mat) <- mat@Dimnames[[1]]
    colnames(Mat) <- mat@Dimnames[[2]]
    return(Mat)
  }
}
##
message("  Welcome to MikuGene Bioinformatics Ecological Community !!! --- Lianhao Song (CodeNight) 2020-11-08 10:55.")
