# YinYangPlot

**YinYangPlot** is an R visualization package designed to generate a new class of differential expression plot: the **Yin–Yang plot**. This plot is intended as a powerful *complement* to the traditional volcano plot, helping researchers focus directly on biologically meaningful signals rather than being overwhelmed by thousands of non-informative points.

The package also provides publication-ready volcano plots and enhanced pirate plots, all optimized for RNA-seq / differential expression workflows (e.g. DESeq2 results).


<img width="700" height="700" alt="image" src="https://github.com/user-attachments/assets/5f34a8cf-4625-4f32-a57a-50e79ffa6c6e" />


---

## ✨ Why Yin–Yang Plot?

Traditional volcano plots display *all* genes simultaneously. While comprehensive, they often suffer from two limitations:

1. **Visual clutter** – tens of thousands of points obscure key signals
2. **Reduced focus** – significant up- and down-regulated genes are not clearly separated

The **Yin–Yang plot** solves this by:

- Displaying **only statistically significant genes**
- Explicitly separating **up-regulated** and **down-regulated** genes into two hemispheres
- Encoding effect size (|log2FC|) as radius
- Preserving symmetry and biological intuition

This design helps scientists:
- Rapidly identify research targets
- Improve exploratory efficiency
- Communicate differential expression patterns more intuitively

---

## 📦 Installation

```r
# install.packages("devtools")
devtools::install_github("LingzhangMeng/YinYangPlot")
```

Load the package:

```r
library(YinYangPlot)
```

---

## 📊 Input Data Requirements

All plotting functions in **YinYangPlot** expect a data frame similar to a **DESeq2 results table**, containing at least:

| Column name        | Description |
|--------------------|-------------|
| `log2FoldChange`   | log2 fold change |
| `padj`             | adjusted p-value |
| `gene_symbol`      | gene symbol (recommended) |
| `GENEBIOTYPE`      | gene biotype (optional but recommended) |

Row names should correspond to gene IDs (e.g. Ensembl IDs).

---

## 1️⃣ Yin–Yang Plot (`yinyang()`)

### Purpose

Generate a Yin–Yang circular plot that highlights **significant up- and down-regulated genes only**.

### Key Features

- Up-regulated genes → **right hemisphere**
- Down-regulated genes → **left hemisphere**
- Radius = |log2FoldChange|
- Optional gene labeling (top up & top down)
- Clean, symmetric, publication-ready design

### Basic Usage

```r
yinyang(
  data = res_df,
  padj_threshold = 0.05,
  lfc_threshold  = 0.5,
  top_n_labels   = 10,
  title = "Yin–Yang Plot of Differential Expression"
)
```

### Output


<img width="800" height="800" alt="Yiinyang plot-2" src="https://github.com/user-attachments/assets/11a84254-ee89-425e-9cb7-2190af3b8336" />




---

## 2️⃣ Enhanced Volcano Plot (`yy_volcano()`)

### Purpose

Provide a cleaner, more informative volcano plot with **biotype-aware coloring** and intelligent gene labeling.

### Key Features

- Automatic DEG classification
- Customizable thresholds
- Top gene labeling
- Compatible with the same input as `yinyang()`

### Example

```r
yy_volcano(
  data = res_df,
  padj_threshold = 0.05,
  lfc_threshold  = 0.5,
  top_n_labels   = 15,
  title = "Differential Expression Volcano Plot"
)
```

<img width="720" height="650" alt="Volcanoplot" src="https://github.com/user-attachments/assets/c4bf4ac8-7c57-4f5c-ab6c-7c3ecddc5774" />


---

## 3️⃣ Pirate Plots for Gene-Level Visualization

### `yy_pirateplot()`

Designed for **multi-gene comparison**, combining:

- Violin distribution
- Beeswarm points
- Mean and confidence intervals

Example:

```r
yy_pirateplot(
  dds = dds,
  genes = c("Aloxe3", "Cxcl13", "Igkc", "Sni1", "Cd163", "Wnt3"),
  group_var = "condition"
)
```

<img width="750" height="650" alt="PiratePlot" src="https://github.com/user-attachments/assets/92b16dd2-85cc-4008-98ab-953169bc6d57" />


---



## 🔬 Typical Workflow

1. Perform RNA-seq differential expression analysis (e.g. DESeq2)
2. Annotate genes (symbol, biotype)
3. Generate:
   - `yy_volcano()` for global overview
   - `yinyang()` for focused biological interpretation
   - `yy_pirateplot()` for target gene validation

---

## 📘 Full Tutorial: End-to-End Test Script

This section converts the **complete testing script** used during package development into a **step-by-step tutorial**. You can copy–paste this directly into your own workflow to reproduce the analysis and figures.



* * *

RNA-Seq Differential Expression Analysis Workflow
=================================================

This repository contains R scripts for processing gene count data, performing differential expression analysis using **DESeq2**, and visualizing results. The pipeline covers everything from data cleaning and metadata preparation to statistical testing.

1\. Environment Setup and Dependencies
--------------------------------------

First, we increase the RStudio console limit to ensure all outputs are captured and load the necessary bioinformatics and visualization libraries.

```
# Extend Rstudio Console to output 10000 lines
rstudioapi::writeRStudioPreference("console_max_lines", 10000L)

# Load Required Packages
library(YinYangPlot)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(apeglm)
library(org.Rn.eg.db)
library(AnnotationDbi)
library(EnsDb.Rnorvegicus.v79)
library(biomaRt)
library(tibble)
library(tidyr)
library(scales)
library(stringr)
library(ggbeeswarm)
library(patchwork)
library(plotly)
```

* * *

2\. Data Loading and Cleaning
-----------------------------

We load the raw counts generated by `featureCounts`. The column names are initially messy (containing full file paths), so we sanitize them and map them to meaningful group names (Control vs. Treated).

```
# Read featureCounts output
counts <- read.table(
   "//path...//gene_counts.txt",   # Locate the counts file on your computer
   header = TRUE,
   row.names = 1,
   comment.char = "#",
   check.names = FALSE,
   stringsAsFactors = FALSE
)

# Remove annotation columns (keeping only count columns)
annotation_cols <- 1:5
counts <- counts[, -annotation_cols]

# Clean column names (removing file paths and extensions)
colnames(counts) <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(counts))
colnames(counts) <- sub("^aligned/", "", colnames(counts))

# Rename samples names with meaningful group names
name_mapping <- c(
   "L1MKL1609037-M_1" = "Control_1",
   "L1MKL1609038-M_2" = "Control_2", 
   "L1MKL1609039-M_3" = "Control_3",
   "L1MKL1609040-M_4" = "Control_4",
   "L1MKL1609041-M_5" = "Control_5",
   "L1MKL1609042-P_1" = "Treated_1",
   "L1MKL1609043-P_2" = "Treated_2",
   "L1MKL1609044-P_3" = "Treated_3",
   "L1MKL1609045-P_4" = "Treated_4",
   "L1MKL1609046-P_5" = "Treated_5"
)

colnames(counts) <- sapply(colnames(counts), function(x) {
   if (x %in% names(name_mapping)) name_mapping[x] else x
})
```

### Data Preview

The resulting count matrix contains Ensembl Gene IDs as rows and sample names as columns:

```
> print(colnames(counts))
 [1] "Control_1" "Control_2" "Control_3" "Control_4" "Control_5" "Treated_1" "Treated_2" "Treated_3" "Treated_4" "Treated_5"

> head(counts, 5)
                    Control_1 Control_2 Control_3 Control_4 Control_5 Treated_1 Treated_2 Treated_3 Treated_4 Treated_5
ENSRNOG00000066169          1         5         0         4         2         3         5         4         1         4
ENSRNOG00000070168          0         0         0         0         0         0         0         0         0         0
ENSRNOG00000070901       5580      3708      6034      5033      4218      3770      4732      5964      5108      4484
ENSRNOG00000018029        182       247       332        57        51       161       181       188        54       296
ENSRNOG00000031391         14        27        13         8         4         0        24        17         4        13
```

* * *

3\. Metadata Preparation
------------------------

Before running DESeq2, define the experimental design by creating a metadata table (`colData`) that assigns each sample to a condition.

```
sample_names <- colnames(counts)
condition <- factor(c(rep("Control", 5), rep("Treated", 5)),
                    levels = c("Control", "Treated"))

colData <- data.frame(
   row.names = sample_names,
   condition = condition,
   stringsAsFactors = TRUE
)
```

### Metadata Summary

```
          condition
Control_1   Control
Control_2   Control
Control_3   Control
Control_4   Control
Control_5   Control
Treated_1   Treated
Treated_2   Treated
Treated_3   Treated
Treated_4   Treated
Treated_5   Treated
```

* * *

4\. DESeq2 Analysis and Filtering
---------------------------------

Initialize the `DESeqDataSet` object and perform a pre-filtering step to remove genes with low counts across all samples, which improves statistical power.

```
# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
   countData = counts,
   colData = colData,
   design = ~ condition
)

# Pre-filter: keep genes with at least 10 total counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
```

**Output:**

```
> class(dds)
[1] "DESeqDataSet"
attr(,"package")
[1] "DESeq2"

Number of genes after filtering: [The count of genes remaining]
```

* * *



* * *

5\. Differential Expression Analysis
------------------------------------

Execute the DESeq2 pipeline, which includes size factor estimation, dispersion estimation, and model fitting.  Apply the `apeglm` shrinkage method to the  $log_{2}$  Fold Changes, to improve visualization and ranking of genes.

```
# Run DESeq2
dds <- DESeq(dds)

# Generate results for Treated vs Control
res <- results(dds, contrast = c("condition", "Treated", "Control"))

# Shrink log2 fold changes for better visualization (removes noise from low-count genes)
resLFC <- lfcShrink(dds, 
                     coef = "condition_Treated_vs_Control", 
                     type = "apeglm")

# Convert to dataframe for further processing
res_df <- as.data.frame(resLFC)
res_df$gene_id <- rownames(res_df)
```

### Analysis Summary

The summary below indicates the distribution of differentially expressed genes at a default alpha of 0.1:

```
out of 19623 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 101, 0.51%
LFC < 0 (down)     : 166, 0.85%
outliers [1]       : 110, 0.56%
low counts [2]     : 2262, 12%
(mean count < 5)
```

* * *

6\. Comprehensive Gene Annotation
---------------------------------

To make the results biologically interpretable, we map Ensembl IDs to symbols and biotypes using a tiered fallback approach: **EnsDb** (local)  $\to$  **biomaRt** (online)  $\to$  **Basic ID retrieval**.

```
# Remove version numbers from Ensembl IDs if present
ensembl_ids <- gsub("\\..*", "", rownames(res_df))

# METHOD 1: Using EnsDb package
annotations <- tryCatch({
   edb <- EnsDb.Rnorvegicus.v79
   ensembldb::select(edb,
                     keys = ensembl_ids,
                     columns = c("GENEID", "SYMBOL", "GENEBIOTYPE", "GENENAME", 
                                 "ENTREZID", "SEQNAME", "GENESEQSTART", "GENESEQEND"),
                     keytype = "GENEID")
}, error = function(e) { NULL })

# METHOD 2: Using biomaRt (Fallback if EnsDb is unavailable)
if (is.null(annotations) || nrow(annotations) == 0) {
   mart <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
   annotations <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", 
                                       "gene_biotype", "description", "entrezgene_id",
                                       "chromosome_name", "start_position", "end_position"),
                         filters = "ensembl_gene_id",
                         values = ensembl_ids,
                         mart = mart)
   colnames(annotations) <- c("GENEID", "SYMBOL", "GENEBIOTYPE", "GENENAME", 
                               "ENTREZID", "SEQNAME", "GENESEQSTART", "GENESEQEND")
}

# Merge annotations back into the results dataframe
res_df <- merge(res_df, annotations, by.x = "gene_id", by.y = "GENEID", all.x = TRUE, sort = FALSE)

# Use Ensembl ID as fallback for missing symbols
res_df$gene_symbol <- ifelse(is.na(res_df$SYMBOL) | res_df$SYMBOL == "", 
                              res_df$gene_id, res_df$SYMBOL)
```

* * *

7\. Categorization and Labeling for Plots
-----------------------------------------

Categorize genes into functional groups (Protein-coding, lncRNA, miRNA, etc.) and define significance thresholds ( $p_{adj}<0.05$  and  $∣log_{2}FC∣>0.3$ ) for downstream visualization like Volcano or Yin-Yang plots.

```
# Define non-coding RNA categories
ncRNA_biotypes <- c("miRNA", "snoRNA", "snRNA", "scaRNA", "vault_RNA",
                     "rRNA", "tRNA", "lncRNA", "misc_RNA", "ribozyme",
                     "sRNA", "piRNA", "siRNA", "guide_RNA")

# Create intelligent labels for plot display
res_df$volcano_label <- sapply(1:nrow(res_df), function(i) {
   row <- res_df[i, ]
   if (!is.na(row$SYMBOL) && row$SYMBOL != row$gene_id && !grepl("^ENSRNOG", row$SYMBOL)) return(row$SYMBOL)
   return(paste0(row$GENEBIOTYPE, "_", substr(row$gene_id, 10, 14)))
})

# Categorize genes into broad groups
res_df$gene_category <- sapply(1:nrow(res_df), function(i) {
   biotype <- res_df[i, ]$GENEBIOTYPE
   if (is.na(biotype)) return("Unknown")
   if (biotype == "protein_coding") return("Protein-coding")
   if (biotype %in% ncRNA_biotypes) return("ncRNA")
   if (grepl("pseudogene", biotype)) return("Pseudogene")
   return("Other")
})

# Add significance status
res_df$significance <- ifelse(
   res_df$padj < 0.05 & !is.na(res_df$padj) & abs(res_df$log2FoldChange) > 0.3,
   ifelse(res_df$log2FoldChange > 0, "Upregulated", "Downregulated"),
   "Not significant"
)
```

* * *



## Response:
In this final section, refine the statistical significance flags and implement an automated file management system to export the results. The workflow concludes by segregating the large dataset into specific files based on gene biotype, facilitating targeted downstream analysis for protein-coding genes, miRNAs, and lncRNAs.

* * *

8\. Significance Flagging and Results Inspection
------------------------------------------------

To identify the most biological relevant changes, apply a stricter significance threshold ( $P<0.05$  and  $∣log_{2}FC∣>1$ ) and sort the data by significance.

```
# Add significance flag for high-confidence genes
res_df$is_significant <- with(res_df, 
                               pvalue < 0.05 & abs(log2FoldChange) > 1 & !is.na(pvalue))

# Sort results by p-value (most significant first)
res_df <- res_df[order(res_df$pvalue, na.last = TRUE), ]
```

### Top 10 Differentially Expressed Genes

The following table shows the top-ranked genes, including `Cxcl13` (highly upregulated) and `Cd163` (significantly downregulated).

```
       gene_id            baseMean  log2FC      pvalue       padj      SYMBOL   GENEBIOTYPE
12013  ENSRNOG00000024899  260.27   3.34    5.14e-15    8.86e-11   Cxcl13   protein_coding
3349   ENSRNOG00000010253  409.53  -1.58    2.91e-09    2.51e-05   Cd163    protein_coding
12204  ENSRNOG00000050996  700.61   1.01    2.35e-08    1.35e-04   Kctd4    protein_coding
7535   ENSRNOG00000005758 1434.42  -0.61    1.77e-07    7.65e-04   Btbd11   protein_coding
2149   ENSRNOG00000018186 2305.23   0.33    6.06e-07    1.85e-03   Lamtor5  protein_coding
```

* * *

9\. Automated Data Export by Biotype
------------------------------------

Because RNA-Seq datasets often contain a mix of different RNA species, the script automatically creates an output directory and saves separate CSV files for each gene biotype (e.g., Protein-coding, miRNA, Pseudogenes).

```
# Define and create output directory
output_dir <- "//path..//Yinyang_CSV files"
if (!dir.exists(output_dir)) {
   dir.create(output_dir, recursive = TRUE)
}

# Save the full master results list
write.csv(res_df, file.path(output_dir, "All_Genes_DESeq2_Results.csv"), row.names = FALSE)

# Clean biotype names and separate into individual files
unique_biotypes <- unique(na.omit(res_df$GENEBIOTYPE))

for (biotype in unique_biotypes) {
   biotype_df <- res_df[!is.na(res_df$GENEBIOTYPE) & res_df$GENEBIOTYPE == biotype, ]
   if (nrow(biotype_df) == 0) next
   
   clean_name <- gsub("[^[:alnum:]]", "_", biotype)
   filepath <- file.path(output_dir, paste0("DESeq2_Results_", clean_name, ".csv"))
   write.csv(biotype_df, filepath, row.names = FALSE)
}
```

### Export Summary

The export process successfully categorized the genes as follows:

| Biotype | Gene Count |
| --- | --- |
| **Protein Coding** | 15,448 |
| **Processed Pseudogene** | 87  |
| **Pseudogene** | 79  |
| **miRNA** | 56  |
| **processed\_transcript** | 32  |
| **lincRNA** | 11  |
| **Other (Mt\_RNA, etc.)** | 13  |

* * *

10\. Summary of Findings
------------------------

The final distribution shows that **75.55%** of the detected transcriptome consists of protein-coding genes, while approximately **23.09%** remains unannotated (`NA`) in the current reference database, highlighting potential areas for novel transcript discovery.



## note:
In this final downstreaming phase, the workflow focuses on high-stringency filtering and organizational export. The script isolates significantly differentially expressed genes and categorizes them into specialized subdirectories for focused biological interpretation.

* * *

11\. High-Confidence Filtering and Sub-sorting
----------------------------------------------

To prioritize genes with the strongest evidence of differential expression, apply a more localized filtering strategy. This step separates genes based on both statistical significance ( $p_{adj}<0.05$ ) and biological magnitude ( $∣log_{2}FC∣>0.3$ ).

```
# Define customized thresholds
padj_cutoff <- 0.05   
lfc_cutoff  <- 0.3    

# Create specialized directory for significant results
base_dir <- file.path(output_dir, "Significant_By_GeneBiotype")
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
```

* * *

12\. Automated Export Pipeline
------------------------------

The script iterates through every detected gene biotype (including those with `NA` annotations) to generate a comprehensive folder structure. Each biotype folder contains:

1.  **All Significant Genes:** Combined list meeting the  $p_{adj}$  and  $LFC$  criteria.
2.  **Upregulated Genes:** Subset with positive fold changes.
3.  **Downregulated Genes:** Subset with negative fold changes.
    
```
for (bt in biotypes) {
   bt_name <- ifelse(is.na(bt), "NA", bt)
   
   # Subset and Filter
   sig_bt <- res_df %>% 
      dplyr::filter(
        if(is.na(bt)) is.na(GENEBIOTYPE) else GENEBIOTYPE == bt,
        padj < padj_cutoff,
        abs(log2FoldChange) > lfc_cutoff
      ) %>%
      dplyr::arrange(padj)
   
   if (nrow(sig_bt) == 0) next
   
   # Create Biotype-specific subfolder
   clean_bt <- gsub("[^[:alnum:]]+", "_", bt_name)
   bt_dir <- file.path(base_dir, clean_bt)
   dir.create(bt_dir, showWarnings = FALSE)
   
   # Save files
   write.csv(sig_bt, file.path(bt_dir, paste0("All_Significant_", clean_bt, ".csv")), row.names = FALSE)
   write.csv(sig_bt %>% filter(log2FoldChange > 0), file.path(bt_dir, "Upregulated.csv"), row.names = FALSE)
   write.csv(sig_bt %>% filter(log2FoldChange < 0), file.path(bt_dir, "Downregulated.csv"), row.names = FALSE)
}
```

### Final Export Log

The console output confirms the successful isolation of high-interest genes:

```
======================================================================
Processing gene biotype: protein_coding 
======================================================================
Total genes: 15448 
Significant genes: 243 

Saved significant protein_coding:      :   40 genes (LFC > 1)
Saved significant pseudogene:          :    1 genes (LFC > 1)
Saved top 20 significant protein_coding:   20 genes
Saved genes with NA biotype: 4722 genes
```

* * *

13\. Project Structure Summary
------------------------------

The final output directory `F:/Chunxia/r_analysis/4. Processing/Yinyang_CSV files` is organized as follows:

*   **`All_Genes_DESeq2_Results.csv`**: The master table of all genes.
*   **`Gene_Biotype_Summary.csv`**: Statistical breakdown of RNA species.
*   **`Significant_By_GeneBiotype/`**:
    *   **`protein_coding/`**: (All\_Significant, Upregulated, Downregulated)
    *   **`miRNA/`**: (Biotype-specific filtered results)
    *   **`NA/`**: (Unannotated significant features)

* * *



## note:
This final section of the workflow transitions from data processing to advanced data visualization. We utilize the `YinYangPlot` package alongside `ggplot2` to create both traditional Volcano plots and innovative Yin-Yang plots, specifically focusing on the protein-coding gene subset.

* * *

14\. Advanced Visualization and Plotting
----------------------------------------

After segregating the data, focus on the **protein-coding** subset to generate publication-quality figures. We utilize two primary visualization methods: a customized Volcano Plot and the unique Yin-Yang Plot.

### Data Preparation for Plotting

We isolate the protein-coding genes and ensure the significance labels are correctly mapped for the plotting functions.

```
# Abstract protein_coding genes for downstream plotting
protein_coding_df <- res_df %>%
   filter(GENEBIOTYPE == "protein_coding")

# Update significance levels for visualization
protein_coding_df <- protein_coding_df %>%
   mutate(
     significance = case_when(
       padj < 0.05 & log2FoldChange > 0.3 ~ "Upregulated",
       padj < 0.05 & log2FoldChange < -0.3 ~ "Downregulated",
       TRUE ~ "NotSignificant"
     )
   )
```

* * *

15\. The Yin-Yang Volcano Plot
------------------------------

The `yy_volcano()` function provides a refined version of the traditional volcano plot, allowing for custom color schemes, label boxing, and automated labeling of top-tier genes.

```
# Generate an enhanced Volcano Plot
yy_volcano(
   df = protein_coding_df,
   fc_cutoff = 0.3,           # Fold change threshold
   p_cutoff = 0.05,           # Adjusted p-value threshold
   top_n_labels = 20,         # Label top 20 DEGs
   boxed_labels = TRUE,       # Box labels for readability
   col_up = "darkgreen",      # Upregulated color
   col_down = "purple",       # Downregulated color
   label_size = 4             
)
```

* * *

16\. The Yin-Yang Plot (Circular DEGs Representation)
-----------------------------------------------------

The `yinyang()` plot offers a novel circular layout where the **radius** represents the magnitude of  $∣log_{2}FC∣$  and the **point size** represents the significance ( $-log_{10}\left(p_{adj}\right)$ ). This provides a visually striking summary of experimental results.

```
# Full customization of the Yin-Yang Plot
yinyang(protein_coding_df,
         title = "Yin-Yang Plot for DEGs (Protein-Coding)",
         subtitle = "Top upregulated & downregulated genes\nRadius = |log2FC|, Point size = -log10(padj)",
         
         # Selection & Labeling
         show_all = TRUE,         
         top_n_up = 10,           
         top_n_down = 10,         
         gene_label_size = 5,     
         
         # Aesthetics
         point_size_range = c(2, 8),
         yin_color = "#90EE90",    # Light Green
         yang_color = "#FFB6C1",   # Light Pink
         point_yin_color = "#006400",
         point_yang_color = "#8B0000",
         border_color = "gold",
         
         # Scaling
         min_radius = 0.4,
         max_radius = 0.9,
         seed = 42
)
```

### Interpretation of Results

In the protein-coding analysis:

*   **Total significant genes identified:** 104
*   **Upregulated:** 37 (e.g., `Cxcl13`, `Wnt3`, `Ighg`)
*   **Downregulated:** 67 (e.g., `Cd163`, `Epha2`)

| Gene Symbol | $log_{2}FC$ | $p_{adj}$ | Description |
| --- | --- | --- | --- |
| **Cxcl13** | 3.34 | $8.86\times 10^{-11}$ | C-X-C motif chemokine ligand 13 |
| **Cd163** | \-1.58 | $2.51\times 10^{-5}$ | Scavenger receptor cysteine-rich type 1 |
| **Wnt3** | 6.71 | $3.03\times 10^{-3}$ | Wnt family member 3 |

* * *

17\. Final Summary of Exported Files
------------------------------------

Upon completion, all statistical tables and summaries are saved to: `//path...//Yinyang_CSV files/`

The exported files include:

*   **`Top50_by_padj.csv`**: Quick reference for the most significant genes per biotype.
*   **`Upregulated.csv` / `Downregulated.csv`**: Direction-specific gene lists.
*   **`All_Significant_..._padj0.05_lfc0.3.csv`**: Comprehensive significant DEG tables.

* * *



## note:
This final section of the workflow utilizes the `yy_pirateplot()` function to visualize the normalized expression levels of specific genes of interest. Pirate plots provide a more transparent alternative to bar charts by combining raw data points, density distributions, and central tendency measures.

* * *

18\. Individual Gene Expression Visualization (Pirate Plots)
------------------------------------------------------------

To validate the findings from the differential expression analysis, we can visualize the raw count distribution for specific candidate genes (e.g., `Cxcl13`, `Cd163`, `Wnt3`). This helps confirm that the statistical significance is driven by consistent trends across all biological replicates.

```
# Define a vector of top-priority genes
my_genes <- c("Aloxe3", "Cxcl13", "Igkc", "Snai1", "Cd163", "Wnt3")

# Generate customized pirate plots
pirate_plot <- yy_pirateplot(
   dds_object = dds,
   results_df = res_df,
   gene_list = my_genes,
   condition_col = "condition",
   group_names = c("Control" = "Control", "Treated" = "Treated"),
   group_colors = c("Control" = "#F781BF", "Treated" = "#92C5DE"),
   plot_style = "pirate",
   layout_strategy = "auto"
)

# Extract and display the combined multi-panel plot
combined_plot <- pirate_plot[["combined_plot"]]
print(combined_plot)
```

### Key Visualization Features

The `yy_pirateplot` generates a comprehensive layout for each gene, including:

*   **Raw Data Points:** Individual sample expression levels to identify outliers.
*   **Density Bean:** Shows the distribution shape of the expression data.
*   **Mean/Median Line:** Clear indicator of the group's central tendency.
*   **Automatic Layout:** The `layout_strategy = "auto"` automatically arranges the six genes into an optimized grid (e.g., a  $2\times 3$  or  $3\times 2$  matrix).

* * *

Summary of Completed Pipeline
-----------------------------

This repository provides a robust end-to-end R workflow for RNA-Seq analysis:

1.  **Preprocessing:** Data cleaning and sample metadata mapping.
2.  **Statistical Analysis:** DESeq2 modeling and `apeglm` LFC shrinkage.
3.  **Annotation:** Multi-tier mapping of Ensembl IDs to symbols and biotypes.
4.  **Data Management:** Automated export of results categorized by biotype.
5.  **Visualization:** High-quality Volcano, Yin-Yang, and Pirate plots for biological interpretation.

* * *



---


This completes a **full reproducible workflow** from raw counts to Yin–Yang visualization.


## 🧠 Design Philosophy

YinYangPlot is built with the following principles:

- **Biology first** – highlight meaningful signals, not noise
- **Robustness** – safe handling of missing annotations
- **Reproducibility** – no fragile dependencies
- **Publication quality** – clean defaults, flexible customization

---

## 📄 Citation

If you use **YinYangPlot** in your research, please cite:

> Meng L. *YinYangPlot: A Yin–Yang visualization framework for differential gene expression analysis*.

(Preprint / manuscript in preparation)

---

## 🤝 Contributions

Issues, feature requests, and pull requests are welcome.

---

## 📬 Contact

**Prof. Dr. rer. nat. Lingzhang Meng** 
E-mail: lzmeng@gxams.org.cn
Center for Medical Big Data and Artificial Intelligence Research
Guangxi Academy of Medical Sciences & People’s Hospital of Guangxi Zhuang Autonomous Region

---

Happy plotting ☯️

