#' returns tissue_positions_list and count matrix given outs folder.
#'
#' @param outs_dir spaceranger outs directory
#'
#' @param ensembl_or_gene format for gene names in count matrix
#'
#' @return list with objects 'tpl' and 'raw_matrix'
#'
#' @keywords spatial
#'
#' @export
#'
#' @examples
#' R code here showing how your function works
.datatable.aware=TRUE
rawmatr_and_tpl <- function(outs_dir, ensembl_or_gene = 'ensembl') {
    raw_matrix <- as.matrix(Matrix::readMM(paste0(outs_dir, "/raw_feature_bc_matrix/matrix.mtx.gz")))
    bcs <- data.table::fread(paste0(outs_dir, "/raw_feature_bc_matrix/barcodes.tsv.gz"),
                                col.names = 'barcode', header = F)
    genes <- data.table::fread(paste0(outs_dir, "/raw_feature_bc_matrix/features.tsv.gz"),
                                        select = c(1:2), col.names = c('ensembl', 'gene'), header = F)
    colnames(raw_matrix) <- bcs$barcode
    if (ensembl_or_gene == 'ensembl') {
       rownames(raw_matrix) <- genes$ensembl
    } else {
       rownames(raw_matrix) <- genes$gene
    }
    umi_cnts <- apply(raw_matrix, 2, sum)

    scales <- rjson::fromJSON(file = paste0(outs_dir, "/spatial/scalefactors_json.json"))
    tpl <- data.table::fread(paste0(outs_dir, "/spatial/tissue_positions_list.csv"),
                    col.names=c('barcode','tissue','row','col','imagerow','imagecol'))
    # scale tissue coordinates for lowres image
    tpl[, imagerow := imagerow * scales$tissue_lowres_scalef]
    tpl[, imagecol := imagecol * scales$tissue_lowres_scalef]
    tpl[, umis := umi_cnts[barcode]]

    return(list('tpl' = tpl, 'raw_matrix' = raw_matrix))
}
