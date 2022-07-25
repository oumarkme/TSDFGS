

#' PEV score
#' 
#' This function calculate prediction error variance (PEV) score <doi:10.1186/s12711-015-0116-6> by given training set and test set. This function is imported from STPGA package developed by Deniz Akdemir <https://CRAN.R-project.org/package=STPGA>.
#' 
#' @import STPGA
#' 
#' @param X A numeric matrix. The training set genotypic information matrix can be given as genotype matrix (coded as -1, 0, 1) or principle component matrix (row: sample; column: marker). 
#' @param X0 A numeric matrix. The test set genotypic information matrix can be given as genotype matrix (coded as -1, 0, 1) or principle component matrix (row: sample; column: marker).
#' @param lambda A float number. Default lambda = 1.
#' @return A floating-point number, PEV score.
#' 
#' @export
#' @examples 
#' data(geno)
#' \dontrun{pev_score(geno[1:50,], geno[51:100,])}
#' 
pev_score = function(X, X0, lambda = 1){
    
    if(! is.matrix(X)) stop("X is not a matrix")
    if(! is.matrix(X0)) stop("X0 is not a matrix")
    if(ncol(X) != ncol(X0))stop("X and X0 matrix have different column size. Please check your input matrix.")
    
    PMAT = rbind(X, X0)
    rownames(PMAT) = 1:nrow(PMAT)
    
    pev_score = STPGA::PEVMAX(1:nrow(X), (nrow(X) + 1): nrow(PMAT), PMAT, lambda, NULL)
    return(pev_score)
}


#' CD score
#' 
#' This function calculate CD-score doi:10.1186/1297-9686-28-4-359 by given training set and test set. This function is imported from STPGA package developed by Deniz Akdemir <https://CRAN.R-project.org/package=STPGA>.
#' 
#' @import STPGA
#' 
#' @param X A numeric matrix. The training set genotypic information matrix can be given as genotype matrix (coded as -1, 0, 1) or principle component matrix (row: sample; column: marker). 
#' @param X0 A numeric matrix. The test set genotypic information matrix can be given as genotype matrix (coded as -1, 0, 1) or principle component matrix (row: sample; column: marker).
#' @param lambda A float number. Default lambda = 1.
#' @param C A numeric matrix. Default as NULL.
#' @return A floating-point number, CD score.
#' 
#' @export
#' @examples 
#' data(geno)
#' \dontrun{cd_score(geno[1:50,], geno[51:100,])}
#' \dontrun{cd_score(geno[1:50,], matrix(geno[51,], nrow = 1))}
#' 
cd_score = function(X, X0, lambda = 1, C = NULL){
    
    if(! is.matrix(X)) stop("X is not a matrix")
    if(! is.matrix(X0)) stop("X0 is not a matrix")
    if(ncol(X) != ncol(X0))stop("X and X0 matrix have different column size. Please check your input matrix.")
    
    PMAT = rbind(X, X0)
    rownames(PMAT) = 1:nrow(PMAT)
    
    cd_score = STPGA::CDMAX(1:nrow(X), (nrow(X) + 1): nrow(PMAT), PMAT, lambda, C)
    return(cd_score)
}