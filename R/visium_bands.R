#' Creates bands inside and outside of tissue.
#'
#' Adds a "band" column to the tissue position list table.
#'
#' @param tpl The tissue position list table, a data.table
#'
#' @return tpl with "band" column, a data.table
#'
#' @keywords spatial
#'
#' @export
#'
#' @examples
#' R code here showing how your function works
visium_bands <- function(tpl) {
    in_tissue_spots <- tpl[tissue == 1]
    out_tissue_spots <- tpl[tissue == 0]
    in_tissue_mat <- as.matrix(in_tissue_spots[, list(imagerow, imagecol)])
    rownames(in_tissue_mat) <- in_tissue_spots$barcode
    out_tissue_mat <- as.matrix(out_tissue_spots[, list(imagerow, imagecol)])
    rownames(out_tissue_mat) <- out_tissue_spots$barcode
    dists <- pdist(in_tissue_mat, out_tissue_mat)
    dist_matr <- as.matrix(dists)
    rownames(dist_matr) <- rownames(in_tissue_mat)
    colnames(dist_matr) <- rownames(out_tissue_mat)
    in_tissue_spots[, min_dist := apply(dist_matr, 1, min)]
    out_tissue_spots[, min_dist := apply(dist_matr, 2, min)]
    band_width <- (1/4)*(max(in_tissue_spots$min_dist)-min(in_tissue_spots$min_dist))
    if (band_width < 30) { # units are pixels.
        band_width <- 30
    }

    # positive bands, outside of tissue.
    bp_ctr <- 1
    while (0 < nrow(out_tissue_spots[(bp_ctr-1)*band_width <= min_dist & min_dist < bp_ctr*band_width])) {
        command <- paste0("out_tissue_spots[, bp", bp_ctr, " := 0]")
        eval(parse(text = command))
        command <- paste0("out_tissue_spots[(bp_ctr-1)*band_width <= min_dist & min_dist < bp_ctr*band_width, bp", bp_ctr, " := 1]")
        eval(parse(text = command))
        command <- paste0("in_tissue_spots[, bp", bp_ctr, " := 0]")
        eval(parse(text = command))
        bp_ctr <- bp_ctr+1
    }
    last_bp <- bp_ctr-1
    if (last_bp == 0) {
        # no spots fall in any of the bands.
        out_tissue_spots[, bp1 := 1]
        in_tissue_spots[, bp1 := 0]
        last_bp <- 1
    }

    # negative bands, inside of tissue.
    bm_ctr <- 1
    while (0 < nrow(in_tissue_spots[(bm_ctr-1)*band_width <= min_dist & min_dist < bm_ctr*band_width])) {
        command <- paste0("in_tissue_spots[, bm", bm_ctr, " := 0]")
        eval(parse(text = command))
        command <- paste0("in_tissue_spots[(bm_ctr-1)*band_width <= min_dist & min_dist < bm_ctr*band_width, bm", bm_ctr, " := 1]")
        eval(parse(text = command))
        command <- paste0("out_tissue_spots[, bm", bm_ctr, " := 0]")
        eval(parse(text = command))
        bm_ctr <- bm_ctr+1
    }
    last_bm <- bm_ctr-1
    if (last_bm == 0) {
        # no spots fall in any of the bands.
        in_tissue_spots[, bm1 := 1]
        out_tissue_spots[, bm1 := 0]
        last_bm <- 1
    }

    bp_vec <- paste0('bp', 1:last_bp)
    bm_vec <- paste0('bm', last_bm:1)

    args_to_paste <- as.list(bm_vec)
    args_to_paste[['sep']] <- ','
    command <- paste0("band_index <- apply(as.matrix(in_tissue_spots[, list(", do.call(paste, args_to_paste), ")]), 1, which.max)")
    eval(parse(text = command))
    in_tissue_spots[, band := bm_vec[band_index]]

    args_to_paste <- as.list(bp_vec)
    args_to_paste[['sep']] <- ','
    command <- paste0("band_index <- apply(as.matrix(out_tissue_spots[, list(", do.call(paste, args_to_paste), ")]), 1, which.max)")
    eval(parse(text = command))
    out_tissue_spots[, band := bp_vec[band_index]]

    return(rbind(in_tissue_spots,
                 out_tissue_spots)[, list(barcode,
                                          tissue,
                                          row,
                                          col,
                                          imagerow,
                                          imagecol,
                                          umis,
                                          band = factor(band,
                                                         levels = c(bm_vec, bp_vec),
                                                         ordered = T))])
}
