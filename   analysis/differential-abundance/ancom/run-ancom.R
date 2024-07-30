
# ANCOM example -------
library(phyloseq)
library(ANCOMBC)

# ANCOM with just diabetes yes/no, all samples ---------

if (ntaxa(ps) == 0 || nsamples(ps) == 0) {
  stop("The phyloseq object has no taxa or samples after subsetting.")
}

ancom_output <- ancombc2(
  ps,
  fix_formula = "Danska.Mice..Diabetic.",
  #fix_formula = "Danska.Mice..Fate.Primary",
  p_adj_method = "holm",
  prv_cut = 0.10,
  lib_cut = 1000,
  group = "Danska.Mice..Fate.Primary",
  struc_zero = TRUE,
  neg_lb = TRUE,
  #tol = 1e-5,
  #max_iter = 100,
  #conserve = TRUE,
  alpha = 0.05,
  global = TRUE
)

ancom_output <- tryCatch({
  ancombc2(
    ps,
    fix_formula = "Danska.Mice..Diabetic.",
    p_adj_method = "holm",
    prv_cut = 0.10,
    lib_cut = 1000,
    group = "Danska.Mice..Diabetic.",
    struc_zero = TRUE,
    neg_lb = TRUE,
    #tol = 1e-5,
    #max_iter = 100,
    #conserve = TRUE,
    alpha = 0.05,
    global = TRUE
  )
}, error = function(e) {
  print(paste("ANCOM-BC2 Error:", e$message))
  NULL
})

res_ancom <- ancom_output$res

# ANCOM with "early" or "late diabetes" -------

ancom_output_J17 <- tryCatch({
  ancombc2(
    ps2,
    fix_formula = "DiabetesType",
    p_adj_method = "holm",
    prv_cut = 0.10,
    lib_cut = 1000,
    group = "DiabetesType",
    struc_zero = TRUE,
    neg_lb = TRUE,
    #tol = 1e-5,
    #max_iter = 100,
    #conserve = TRUE,
    alpha = 0.05,
    global = TRUE
  )
}, error = function(e) {
  print(paste("ANCOM-BC2 Error:", e$message))
  NULL
})

## Parameters:
# fix_formula:
# group:
# normalized? Y, 

July4 <- tryCatch({
  ancombc2(
    ps3,
    #tax_level = "BLAST_species.(as.done.by.EAV)",
    #tax_level = "Sample_ID",
    #fix_formula = "DiabetesType",
    fix_formula = "Diabetes",
    p_adj_method = "holm",
    prv_cut = 0.05, # No.Samples X prv_cut value (66 * 0.1 is 6.6)
    lib_cut = 0,
    #group = "DiabetesType",
    group = "Diabetes",
    struc_zero = TRUE,
    #struc_zero = FALSE,
    neg_lb = TRUE,
    #tol = 1e-5,
    #max_iter = 100,
    #conserve = TRUE,
    alpha = 0.05
  )
}, error = function(e) {
  print(paste("ANCOM-BC2 Error:", e$message))
  NULL
})

## Manually group (aka. collapse) in the ps object

meta2[rownames(meta2) %in% AF.colnames, ]
meta2[rownames(meta2) %in% AF.colnames, colnames(meta2) %in% c("CollectionTime", "MouseID", "Danska Mice::Age At Retirement", "DiabetesType")]

a.res <- a[["res"]]
a.res$significance <- ifelse(a.res$q_DiabetesTypelate < 0.05, "Significant", "Not Significant")

p <- ggplot(a.res, aes(x = lfc_DiabetesTypelate, y = -log10(q_DiabetesTypelate), color = significance)) +
  geom_point() +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot of Differentially Abundant Taxa",
    x = "Log Fold Change",
    y = "-log10(Adjusted P-value)",
    color = "Significance"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

p2 <- ggplot(a.res, aes(x = lfc_DiabetesTypelate, y = W_DiabetesTypelate, color = significance)) +
  geom_point() +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot of Differentially Abundant Taxa",
    x = "Log Fold Change",
    y = "-log10(Adjusted P-value)",
    color = "Significance"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )
