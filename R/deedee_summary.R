deedee_summary <- function(deedee_list,
                           pthresh = 0.05,
                           scatter_select1 = 1,
                           scatter_select2 = 2,
                           scatter_color_by = "pval1",
                           heatmap_show_first = 25,
                           heatmap_show_gene_names = FALSE,
                           heatmap_dist = "euclidean",
                           heatmap_clust = "average",
                           heatmap_show_na = FALSE,
                           venn_mode = "both",
                           upset_mode = "both_colored",
                           upset_min_setsize = 10,
                           qqmult_ref = 1,
                           cat_ref = 1,
                           cat_maxrank = 1000,
                           cat_mode = "up") {

  # ----------------------------- argument check ------------------------------
  checkmate::assert_list(deedee_list, type = "data.frame", min.len = 2)
  for (i in 1:length(deedee_list)) {
    checkmate::assert_data_frame(deedee_list[[i]], type = "numeric")
  }
  checkmate::assert_number(pthresh, lower = 0, upper = 1)

  checkmate::assert_number(scatter_select1, lower = 1, upper = length(deedee_list))
  checkmate::assert_number(scatter_select2, lower = 1, upper = length(deedee_list))
  choices <- c("pval1", "pval2")
  checkmate::assert_choice(scatter_color_by, choices)

  checkmate::assert_number(heatmap_show_first, lower = 1)
  checkmate::assert_logical(heatmap_show_gene_names)
  choices1 <- c("euclidean", "manhattan", "pearson", "spearman")
  checkmate::assert_choice(heatmap_dist, choices1)
  choices2 <- c("single", "complete", "average", "centroid")
  checkmate::assert_choice(heatmap_clust, choices2)
  checkmate::assert_logical(heatmap_show_na)

  choices <- c("up", "down", "both")
  checkmate::assert_choice(venn_mode, choices)

  choices <- c("up", "down", "both", "both_colored")
  checkmate::assert_choice(upset_mode, choices)
  checkmate::assert_number(upset_min_setsize, lower = 0)

  checkmate::assert_number(qqmult_ref, lower = 1, upper = length(deedee_list))

  checkmate::assert_number(cat_ref, lower = 1, upper = length(deedee_list))
  checkmate::assert_number(cat_maxrank, lower = 1)
  choices <- c("up", "down", "both")
  checkmate::assert_choice(cat_mode, choices)

  # -------------------------- calling the functions ---------------------------
  sc <- DeeDee::deedee_scatter(data = deedee_list,
                         pthresh = pthresh,
                         select1 = scatter_select1,
                         select2 = scatter_select2,
                         color_by = scatter_color_by)

  hm <- DeeDee::deedee_heatmap(data = deedee_list,
                         pthresh = pthresh,
                         show_first = heatmap_show_first,
                         show_gene_names = heatmap_show_gene_names,
                         dist = heatmap_dist,
                         clust = heatmap_clust,
                         show_na = heatmap_show_na)

  vn <- DeeDee::deedee_venn(data = deedee_list,
                      pthresh = pthresh,
                      mode = venn_mode)

  us <- DeeDee::deedee_upset(data = deedee_list,
                       pthresh = pthresh,
                       mode = upset_mode,
                       min_setsize = upset_min_setsize)

  qq <- DeeDee::deedee_qqmult(data = deedee_list,
                        pthresh = pthresh,
                        ref = qqmult_ref)

  ct <- DeeDee::deedee_cat(data = deedee_list,
                     pthresh = pthresh,
                     ref = cat_ref,
                     maxrank = cat_maxrank,
                     mode = cat_mode)

  pdf("~/Desktop/deedee_summary.pdf")
    print(sc)
    print(hm)
    print(vn)
    print(us)
    print(qq)
    print(ct)
  dev.off()

}
