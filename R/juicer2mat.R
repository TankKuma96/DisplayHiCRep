#' calculate the reproducibility of HiChIP data in .hic format and display all inter-chromosome and intra-chromosome in data frame and heat map
#'
#' @param dic path of the contant maps exist
#' @param chromfile_path a file that contain chromosome length information that is used for calculating bin size
#' @param pat1 pattern of files for chromosome 1. Default is "_r1".
#' @param pat2 pattern of files for chromosome 2. Default is "_r2".
#' @param len amount of chromosome that will be compared about, which count from chr1
#' @param maxrange An integer indicating the maximum distance of interaction that is considered. Default is 5000000.
#' @param resolution An integer indicating the resolution level of the HiChIP matrix. Default is 500000.
#' @param mapHeat If need to display in heat map or not. Default is TRUE
#' @return matrix giving number of all inter-chromosome and intra-chromosome comparasion
#' @export

juicer2mat = function(dic, chromfile_path, pat1 = "_r1", pat2 = "_r2", len, maxrange = 5000000, resolution = 500000, mapHeat = TRUE){

  setwd(dic)

  chrom_file = read.table(chromfile_path)

  rep1 <- list.files(pattern = pat1)
  rep2 <- list.files(pattern = pat2)
  nameorder <- gsub("(_chr)([1-9])(_)","_chr0\\2_",rep1)

  rep1files <- rep1[str_order(nameorder)]
  rep2files <- rep2[str_order(nameorder)]
  rep1d <- lapply(rep1files, function(X) read_table2(X, col_names = FALSE))
  rep2d <- lapply(rep2files, function(X) read_table2(X, col_names = FALSE))
  filepat <- gsub("_r1.txt","",rep1files) %>%
    gsub(pattern = "[[:alpha:]]*_[[:digit:]]*kb_", replacement = "")


  chr_f = sapply(filepat, function(X) gsub("_VS_chr[[:alnum:]]*","",X))
  chr_b = sapply(filepat, function(X) gsub("chr[[:alnum:]]*_VS_","",X))
  pktb = tibble(chr_f,chr_b)

  B <- chrom_file %>%
    set_names(c("chr", "len")) %>%
    filter(chr %in% pktb$chr_f) %>%
    mutate(len = ceiling(len/resolution))
  points <- cumsum(B$len)
  parafile <- data.frame("chr" = B$chr,
                         "startpt" = c(0,points[1:(length(points)-1)])+1,
                         "endpt" = points)
  fulld <- list()
  for(i in seq_along(rep1d)){
    chr_f <- pktb$chr_f[i]
    chr_b <- pktb$chr_b[i]
    if(chr_f == chr_b){
      sd1 <- sparseMatrix(i=(rep1d[[i]]$X1)/resolution,
                          j=(rep1d[[i]]$X2)/resolution,
                          x = rep1d[[i]]$X3,
                          dims = c(filter(B, chr == chr_f)[2],
                                   filter(B, chr == chr_b)[2]),
                          index1 = FALSE)
      full1 <- as.matrix(sd1) + t(as.matrix(sd1))
      diag(full1) <- diag(full1)/2

      sd2 <- sparseMatrix(i=(rep2d[[i]]$X1)/resolution,
                          j=(rep2d[[i]]$X2)/resolution,
                          x = rep2d[[i]]$X3,
                          dims = c(filter(B, chr == chr_f)[2],
                                   filter(B, chr == chr_b)[2]),
                          index1 = FALSE)
      full2 <- as.matrix(sd2) + t(as.matrix(sd2))
      diag(full2) <- diag(full2)/2

      fulld[[i]] <- list("rep1" = full1, "rep2" = full2)
    }
    else{
      sd1 <- sparseMatrix(i=(rep1d[[i]]$X1)/resolution,
                          j=(rep1d[[i]]$X2)/resolution,
                          x = rep1d[[i]]$X3,
                          dims = c(filter(B, chr == chr_f)[2],
                                   filter(B, chr == chr_b)[2]),
                          index1 = FALSE)
      full1 <- as.matrix(sd1)

      sd2 <- sparseMatrix(i=(rep2d[[i]]$X1)/resolution,
                          j=(rep2d[[i]]$X2)/resolution,
                          x = rep2d[[i]]$X3,
                          dims = c(filter(B, chr == chr_f)[2],
                                   filter(B, chr == chr_b)[2]),
                          index1 = FALSE)
      full2 <- as.matrix(sd2)

      fulld[[i]] <- list("rep1" = full1, "rep2" = full2)
    }
  }

  superMatrix1 <- superMatrix2 <- matrix(0, ncol = max(points), nrow = max(points))
  for(i in seq_along(fulld)){
    chr1 <- pktb$chr_f[i]
    chr2 <- pktb$chr_b[i]

    if(chr1 == chr2){
      idx <- which(parafile$chr == as.character(chr1))
      x <- parafile[idx,] %>% mutate(chr = NULL)
      superMatrix1[x$startpt:x$endpt, x$startpt:x$endpt] <- fulld[[i]]$rep1
      superMatrix2[x$startpt:x$endpt, x$startpt:x$endpt] <- fulld[[i]]$rep2
    } else
    {
      idx1 <- which(parafile$chr == as.character(chr1))
      idx2 <- which(parafile$chr == as.character(chr2))
      x <- parafile[idx1,] %>% mutate(chr = NULL)
      y <- parafile[idx2,] %>% mutate(chr = NULL)

      superMatrix1[x$startpt:x$endpt, y$startpt:y$endpt] <- fulld[[i]]$rep1
      superMatrix1[y$startpt:y$endpt, x$startpt:x$endpt] <- t(fulld[[i]]$rep1)

      superMatrix2[x$startpt:x$endpt, y$startpt:y$endpt] <- fulld[[i]]$rep2
      superMatrix2[y$startpt:y$endpt, x$startpt:x$endpt] <- t(fulld[[i]]$rep2)
    }
  }

  SMmat1 = fast.mean.filter(superMatrix1, 1)
  SMmat2 = fast.mean.filter(superMatrix2, 1)

  ###################################################
  # Break super matrix back into pieces for hicrep to take

  matPiece1 <- matPiece2 <- list()
  for(i in seq_along(fulld)){
    chr1 <- pktb$chr_f[i]
    chr2 <- pktb$chr_b[i]

    if(chr1 == chr2){
      idx <- which(parafile$chr == as.character(chr1))
      x <- parafile[idx,] %>% mutate(chr = NULL)
      matPiece1[[i]] <- SMmat1[x$startpt:x$endpt, x$startpt:x$endpt]
      matPiece2[[i]] <- SMmat2[x$startpt:x$endpt, x$startpt:x$endpt]
    } else
    {
      idx1 <- which(parafile$chr == as.character(chr1))
      idx2 <- which(parafile$chr == as.character(chr2))
      x <- parafile[idx1,] %>% mutate(chr = NULL)
      y <- parafile[idx2,] %>% mutate(chr = NULL)

      matPiece1[[i]] <- SMmat1[x$startpt:x$endpt, y$startpt:y$endpt]
      matPiece2[[i]] <- SMmat2[x$startpt:x$endpt, y$startpt:y$endpt]
    }
  }
  names(matPiece1) <- names(matPiece2) <- paste0(pktb$chr_f, "VS", pktb$chr_b)

  cor_out <- rep(NA, nrow(pktb))
  for(i in seq_along(matPiece1)){
    chr1 <- pktb$chr_f[i]
    chr2 <- pktb$chr_b[i]

    if(chr1 == chr2){
      M1 <- as.matrix(matPiece1[[i]])
      M2 <- as.matrix(matPiece2[[i]])
      hicscc = get.scc(M1, M2, resol = resolution, 0, lbr = 0, ubr = maxrange)
      cor_out[i] = hicscc$scc
    }
    else{
      M1 <- as.matrix(matPiece1[[i]])
      M2 <- as.matrix(matPiece2[[i]])
      cor_out[i] = cor(c(M1), c(M2), method = "spearman")
    }
  }
  mattb <- data.frame(pktb, "cor_out" = cor_out)

  namelist <- unique(pktb$chr_f)
  DisplayMat <- matrix(0, nrow = len, ncol= len,
                       dimnames = list(namelist, namelist))

  for (i in seq(nrow(mattb))){
    chr1 <- pktb$chr_f[i]
    chr2 <- pktb$chr_b[i]

    if(chr1 == chr2){
      ridx <- which(namelist == chr1)
      DisplayMat[ridx, ridx] = mattb[i,"cor_out"]
    } else{
      ridx <- which(namelist == chr1)
      cidx <- which(namelist == chr2)
      DisplayMat[ridx, cidx] = mattb[i,"cor_out"]
      DisplayMat[cidx, ridx] = mattb[i,"cor_out"]
    }
  }
  if (mapHeat == TRUE){
    result = DisplayMat
    print(result)
    heatmap(result, Rowv=NA, Colv=NA, scale='none',main = "Heatmap for result",
            xlab = "chrom #2", ylab = "chrom #1",RowSideColors = heat.colors(nrow(result)))
  }else{
    return(DisplayMat)
  }
}


