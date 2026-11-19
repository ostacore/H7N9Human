install.packages("usethis")
library(usethis)


# Configure your Git identity (use your GitHub email)
use_git_config(user.name = "Your Name", user.email = "your_email@example.com")



# WHICH QUESTION IS ANSWERED WHERE
# --------------------------------
#   Question 1 (H7N9 similar/different to other viruses?)  -> PART 8
#   Question 2 (which time point is most dysregulated?)    -> PART 7
#   Question 3 (which pathways, and how do they differ?)   -> PART 9



library(limma)            # statistics for comparing gene expression between groups
library(ggplot2)          # making nice plots
library(clusterProfiler)  # gene enrichment tools (also converts gene ID types)
library(org.Hs.eg.db)     # database of human genes (used to convert gene names)
library(ReactomePA)       # tests for enrichment of Reactome pathways
library(enrichplot)       # plots for enrichment results

setwd("/Users/emmakomers/Desktop/projectB/")
mydata <- read.delim("Project-H7N9Human.txt", check.names = FALSE)


# The viruses and time points in the experiment.
# They must match the names used in your sample names exactly.
viruses <- c("H7N9", "H7N7", "H3N2", "VN1203")
times   <- c("03h", "07h", "12h", "24h")
padj_cutoff  <- 0.05
logfc_cutoff <- 1

# SAFETY CHECKS: stopifnot() halts the script with an error if something is wrong.
# This is better than getting confusing results later.
stopifnot(all(sapply(mydata, is.numeric)))      # every column holds numbers
stopifnot(!anyNA(mydata))                       # no missing values
stopifnot(all(as.matrix(mydata) > 0))           # all values positive (needed for log2)
stopifnot(!anyDuplicated(rownames(mydata)))     # each gene name appears only once


# =============================================================================
# PART 3: CREATE A SAMPLE INFORMATION TABLE ("metadata")
# =============================================================================
# R needs to know which virus and time point each sample belongs to. We get
# this from the sample names, which look like:  H7N9_24h_1
#                                               virus_time_replicate

# Start a table with one row per sample.
metadata <- data.frame(sample = colnames(mydata))

# sub() finds a pattern in text and replaces it. Here we use it to cut out
# the pieces of each sample name.
metadata$virus     <- sub("_.*", "", metadata$sample)                  # text before first "_"
metadata$timepoint <- sub("^[^_]+_([^_]+)_.*$", "\\1", metadata$sample) # text between the "_"s
metadata$replicate <- sub("^.*_", "", metadata$sample)                 # text after last "_"

# Turn timepoint into an ordered factor so plots show 3h, 7h, 12h, 24h in order
# (otherwise R sorts alphabetically).
metadata$timepoint <- factor(metadata$timepoint, levels = times)

# Show how many samples (replicates) we have for each virus/time combination.
print(table(metadata$virus, metadata$timepoint))

# More safety checks: all time points were recognised, and all virus names are expected.
stopifnot(!anyNA(metadata$timepoint))
stopifnot(all(metadata$virus %in% c("mock", viruses)))


# =============================================================================
# PART 4: LOG-TRANSFORM THE DATA AND DO QUALITY CHECKS
# =============================================================================

# Take log2 of every value. Gene expression values are very skewed; the log
# makes them more symmetric, which the statistics below assume.
# NOTE: only do this if your data are NOT already log-transformed!
expression_log <- log2(as.matrix(mydata))

# --- Boxplot of every sample ---------------------------------------------------
# All boxes should look roughly similar. A very different box = a possible
# problem sample.
pdf("01_sample_boxplots.pdf", width = 14, height = 6)   # open a PDF file
boxplot(expression_log, las = 2, cex.axis = 0.5,
        main = "Log2 expression by sample", ylab = "Log2 expression")
dev.off()                                               # close (save) the PDF

# --- PCA plot --------------------------------------------------------------------
# PCA squeezes thousands of genes into 2 axes (PC1, PC2) so we can SEE which
# samples are similar. Samples close together behave similarly.

# Genes that never change (zero variation) break the scaling step, so drop them.
variable_genes <- apply(expression_log, 1, sd) > 0

# prcomp() needs samples as rows, so t() flips the table (transposes it).
pca <- prcomp(t(expression_log[variable_genes, , drop = FALSE]), scale. = TRUE)

# Put the first two PCA axes next to the sample information for plotting.
pca_plot <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], metadata)

# Colour = virus, shape = time point.
pca_figure <- ggplot(pca_plot, aes(PC1, PC2, color = virus, shape = timepoint)) +
  geom_point(size = 2.5) +
  theme_minimal() +
  labs(title = "PCA of samples")
print(pca_figure)                                        # show in RStudio
ggsave("02_PCA.png", pca_figure, width = 8, height = 6)  # save as an image
# Look for odd samples that sit far from their replicates before going further.


# =============================================================================
# PART 5: BUILD THE STATISTICAL MODEL (limma)
# =============================================================================
# Idea: treat every virus + time combination (e.g. "H7N9_24h") as one group,
# then compare groups to each other.

# Make a "group" label for each sample, e.g. "H7N9_24h" or "mock_03h".
metadata$group <- factor(paste(metadata$virus, metadata$timepoint, sep = "_"))

# The design matrix tells limma which group each sample belongs to.
# "~ 0 + group" means: one column per group, no intercept.
design <- model.matrix(~ 0 + group, data = metadata)
colnames(design) <- sub("^group", "", colnames(design))  # tidy the column names

# Safety check: the samples in metadata are in the same order as in the data.
stopifnot(identical(metadata$sample, colnames(expression_log)))

# Fit a linear model to every gene at once.
fit <- lmFit(expression_log, design)


# =============================================================================
# PART 6: DEFINE THE COMPARISONS ("contrasts") AND FIND DEGs
# =============================================================================
# A contrast is a comparison written as "groupA - groupB". We build two sets:
#   (a) each virus vs. mock at the SAME time point   (what the virus changed)
#   (b) H7N9 vs. each other virus at each time point (how H7N9 differs)
# Two empty "shopping lists" to fill up in the loops below:
comparison_formulas <- character(0)  # the formulas limma needs, e.g. "H7N9_24h - mock_24h"
comparison_names    <- character(0)  # friendly names, e.g. "H7N9_vs_mock_24h"

# (a) Each virus vs. mock. A "for" loop repeats the same steps for every item.
for (v in viruses) {                 # for each virus...
  for (tp in times) {                # ...and each time point:
    comparison_formulas <- c(comparison_formulas, paste0(v, "_", tp, " - mock_", tp))
    comparison_names    <- c(comparison_names,    paste0(v, "_vs_mock_", tp))
  }
}

# (b) H7N9 vs. each of the other viruses.
for (v in c("H7N7", "H3N2", "VN1203")) {
  for (tp in times) {
    comparison_formulas <- c(comparison_formulas, paste0("H7N9_", tp, " - ", v, "_", tp))
    comparison_names    <- c(comparison_names,    paste0("H7N9_vs_", v, "_", tp))
  }
}

# Give limma the comparisons, then calculate statistics.
# eBayes() "borrows strength" across genes, which makes results more reliable
# when there are few replicates.
contrast_matrix <- makeContrasts(contrasts = comparison_formulas, levels = design)
colnames(contrast_matrix) <- comparison_names
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Collect the full results table for each comparison into a list.
# Each table has, for every gene:
#   logFC     = log2 fold change (positive = higher in the FIRST group)
#   P.Value   = raw p-value
#   adj.P.Val = p-value corrected for testing thousands of genes (use this one!)
results <- list()
for (name in comparison_names) {
  results[[name]] <- topTable(fit2, coef = name, number = Inf, adjust.method = "BH")
  results[[name]]$gene <- rownames(results[[name]])   # keep gene names as a column
}

# Count the DEGs in every comparison and store them in a summary table.
summary_table <- data.frame(comparison = comparison_names,
                            total_DEGs = 0, upregulated = 0, downregulated = 0)

for (i in seq_along(comparison_names)) {
  x <- results[[comparison_names[i]]]
  # TRUE/FALSE for each gene: significant AND big enough change?
  up   <- x$adj.P.Val < padj_cutoff & x$logFC >=  logfc_cutoff
  down <- x$adj.P.Val < padj_cutoff & x$logFC <= -logfc_cutoff
  # sum() of TRUE/FALSE counts the TRUEs. na.rm = TRUE ignores missing values.
  summary_table$upregulated[i]   <- sum(up, na.rm = TRUE)
  summary_table$downregulated[i] <- sum(down, na.rm = TRUE)
  summary_table$total_DEGs[i]    <- sum(up | down, na.rm = TRUE)
}
print(summary_table)


# =============================================================================
# PART 7: QUESTION 2 - AT WHICH TIME POINT IS THE RESPONSE MOST DYSREGULATED?
# =============================================================================
# Approach: count DEGs (virus vs. mock) at each time point. The time point
# with the most DEGs is the most dysregulated by this measure.

# Keep only the "virus vs mock" rows of the summary table.
mock_summary <- summary_table[grepl("_vs_mock_", summary_table$comparison), ]

# Split the comparison name back into virus and time point for plotting.
mock_summary$virus     <- sub("_vs_mock_.*", "", mock_summary$comparison)
mock_summary$timepoint <- sub(".*_vs_mock_", "", mock_summary$comparison)
mock_summary$timepoint <- factor(mock_summary$timepoint, levels = times)

# Line plot: one line per virus, number of DEGs over time.
plot_counts <- ggplot(mock_summary,
                      aes(timepoint, total_DEGs, color = virus, group = virus)) +
  geom_line() + geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Differentially expressed genes vs mock",
       x = "Timepoint", y = "Number of DEGs")
print(plot_counts)
ggsave("04_DEG_counts_by_time.png", plot_counts, width = 8, height = 5)

# CAUTION: more DEGs means more genes passed the cutoffs, not necessarily a
# bigger biological effect.


# =============================================================================
# PART 8: QUESTION 1 - IS H7N9 SIMILAR TO OR DIFFERENT FROM OTHER VIRUSES?
# =============================================================================
# We look at the 24h time point. Four comparisons matter here:
comparisons_24h <- c("H7N9_vs_mock_24h",    # what H7N9 changes compared to no virus
                     "H7N9_vs_H7N7_24h",    # H7N9 directly vs each other virus
                     "H7N9_vs_H3N2_24h",
                     "H7N9_vs_VN1203_24h")

# Show the DEG counts for these four comparisons.
# FEW DEGs in "H7N9 vs other virus" = the two viruses behave SIMILARLY.
# MANY DEGs = they behave DIFFERENTLY.
print(summary_table[summary_table$comparison %in% comparisons_24h, ])

# Reminder: positive logFC = higher in H7N9; negative logFC = lower in H7N9.

# --- Volcano plot: H7N9 vs mock at 24h ------------------------------------------
# Each dot is a gene. Far left/right = big change; high up = strong evidence.
volcano <- results[["H7N9_vs_mock_24h"]]

# Label each gene as higher, lower, or not a DEG (used for colouring).
volcano$status <- "Not DEG"
volcano$status[volcano$adj.P.Val < padj_cutoff & volcano$logFC >=  logfc_cutoff] <- "Higher in H7N9"
volcano$status[volcano$adj.P.Val < padj_cutoff & volcano$logFC <= -logfc_cutoff] <- "Lower in H7N9"

# -log10 of the p-value makes small p-values into big numbers (higher = stronger).
# pmax() avoids taking the log of exactly zero.
volcano$minus_log10_p <- -log10(pmax(volcano$adj.P.Val, .Machine$double.xmin))

volcano_plot <- ggplot(volcano, aes(logFC, minus_log10_p, color = status)) +
  geom_point(alpha = 0.5, size = 1) +
  theme_minimal() +
  labs(title = "H7N9 vs mock at 24h",
       x = "Log2 fold change", y = "-log10 adjusted p-value")
print(volcano_plot)
ggsave("05_H7N9_24h_volcano.png", volcano_plot, width = 8, height = 6)


# =============================================================================
# PART 9: QUESTION 3 - WHICH PATHWAYS ARE DYSREGULATED, AND HOW DO THEY DIFFER?
# =============================================================================
# "Pathway enrichment" asks: do our DEGs contain more genes from a given
# biological pathway (e.g. "Interferon signalling") than expected by chance?
# We use the Reactome pathway database, again at 24h.

# Step 1: Reactome uses Entrez ID numbers, but our data uses gene SYMBOLS
# (e.g. "TP53"). bitr() converts between them. Some names (e.g. LOC/XLOC
# entries) won't convert and are dropped.
all_converted <- clusterProfiler::bitr(rownames(mydata), fromType = "SYMBOL",
                                       toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Step 2: The "background" (universe) = all genes we could convert. Enrichment
# compares our DEGs against this list, so it must be the same conversion we use below.
background_entrez <- unique(all_converted$ENTREZID)
cat("Mapped symbols:", length(unique(all_converted$SYMBOL)), "of", nrow(mydata), "\n")

# Step 3: Run the enrichment for each of the four 24h comparisons, separately
# for genes that are higher ("up") and lower ("down") in H7N9.
reactome_results <- list()   # results are kept inside R (no CSV files)

for (comparison in comparisons_24h) {
  reactome_results[[comparison]] <- list()
  cat("Reactome analysis:", comparison, "\n")   # progress message

  x <- results[[comparison]]
  up_symbols   <- x$gene[x$adj.P.Val < padj_cutoff & x$logFC >=  logfc_cutoff]
  down_symbols <- x$gene[x$adj.P.Val < padj_cutoff & x$logFC <= -logfc_cutoff]

  for (direction in c("up", "down")) {
    # Choose the up or down gene list.
    selected_symbols <- if (direction == "up") up_symbols else down_symbols

    # Convert those gene symbols to Entrez IDs using our conversion table.
    selected_entrez <- unique(all_converted$ENTREZID[all_converted$SYMBOL %in% selected_symbols])

    # If no genes are left, skip ahead ("next") to the next direction.
    if (length(selected_entrez) == 0) {
      cat("No mapped genes for", comparison, direction, "\n")
      next
    }

    # The actual enrichment test.
    pathway_result <- ReactomePA::enrichPathway(
      gene = selected_entrez,       # our DEGs
      universe = background_entrez, # all genes tested
      organism = "human",
      pvalueCutoff = 0.05,          # keep significant pathways only
      pAdjustMethod = "BH",         # correct for testing many pathways
      readable = TRUE               # show gene symbols, not ID numbers
    )

    # Store the result as a normal table inside our list.
    pathway_table <- as.data.frame(pathway_result)
    reactome_results[[comparison]][[direction]] <- pathway_table

    # If any pathways were found, draw a dot plot of the top 15 and save it.
    # (Dot size = number of genes; colour = significance.)
    if (nrow(pathway_table) > 0) {

      # Make a clean figure title WITHOUT underscores:
      # gsub("_", " ", ...) replaces every "_" with a space,
      # so "H7N9_vs_mock_24h" becomes "H7N9 vs mock 24h".
      clean_comparison <- gsub("_", " ", comparison)
      clean_direction  <- if (direction == "up") "upregulated" else "downregulated"

      # label_format = 35 wraps long pathway names onto a new line after
      # about 35 characters, so they do not run off the figure.
      # theme(axis.text.y = ...) makes the pathway names smaller (size 9) and
      # puts a little extra space between lines (lineheight), so that wrapped
      # names do not touch the names above and below them.
      p <- enrichplot::dotplot(pathway_result, showCategory = 15,
                               label_format = 35) +
        ggtitle(paste(clean_comparison, "-", clean_direction, "genes")) +
        theme(axis.text.y = element_text(size = 9, lineheight = 0.9))
      print(p)

      # The saved figure is made taller (height = 10) to give every pathway
      # name enough room. If names still touch, increase height (e.g. 12).
      ggsave(paste0("Reactome_", comparison, "_", direction, ".png"),
             p, width = 11, height = 10)
    }
  }
}
# =============================================================================
# LIMITATIONS TO KEEP IN MIND
# =============================================================================
# - Not every gene name converts to an Entrez ID, so some genes are lost in Part 9.
# - We assumed the data were normalized; this wasn't verified.
# - The comparisons test differences in gene EXPRESSION, not pathway ACTIVITY.
