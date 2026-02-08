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

> **Note**: File paths are examples and should be adapted to your local system.

---

## Step 1. Load Required Packages

```r
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

---

## Step 2. Load and Prepare Count Data

```r
counts <- read.table(
  "/path/gene_counts.txt",
  header = TRUE,
  row.names = 1,
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Remove annotation columns
counts <- counts[, -c(1:5)]
counts <- counts[, -annotation_cols]

# Clean sample names
colnames(counts) <- sub("\.Aligned\.sortedByCoord\.out\.bam$", "", colnames(counts))
```


[output]
```
> colnames(counts)
[1] "Treated_1" "Treated_2" "Treated_3" "Treated_4" "Treated_5"
```


```r
# Clean column names
colnames(counts) <- sub("\\.Aligned\\.sortedByCoord\\.out\\.bam$", "", colnames(counts))
colnames(counts) <- sub("^aligned/", "", colnames(counts))

# Check column names again
colnames(counts)
```

[output]
```
> colnames(counts)
[1] "Treated_1" "Treated_2" "Treated_3" "Treated_4" "Treated_5"
```

```r
head(counts, 5)
```

[output]
```
> head(counts,5)
                   Treated_1 Treated_2 Treated_3 Treated_4 Treated_5
ENSRNOG00000066169         3         5         4         1         4
ENSRNOG00000070168         0         0         0         0         0
ENSRNOG00000070901      3770      4732      5964      5108      4484
ENSRNOG00000018029       161       181       188        54       296
ENSRNOG00000031391         0        24        17         4        13
```

```r
name_mapping <- c(
  "L1MKL1609037-M_1" = "Control_1",
  "L1MKL1609042-P_1" = "Treated_1"
)

colnames(counts) <- sapply(colnames(counts), function(x) {
  if (x %in% names(name_mapping)) name_mapping[x] else x
})
```

---

## Step 3. Sample Metadata

```r
condition <- factor(c(rep("Control", 5), rep("Treated", 5)),
                    levels = c("Control", "Treated"))

colData <- data.frame(
  row.names = colnames(counts),
  condition = condition
)
```

---

## Step 4. Differential Expression Analysis (DESeq2)

```r
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = colData,
  design = ~ condition
)

# Filter low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "Treated", "Control"))

# Shrink fold changes
resLFC <- lfcShrink(
  dds,
  coef = "condition_Treated_vs_Control",
  type = "apeglm"
)

res_df <- as.data.frame(resLFC)
res_df$gene_id <- rownames(res_df)
```

---

## Step 5. Gene Annotation (EnsDb + biomaRt Fallback)

```r
ensembl_ids <- gsub("\..*", "", res_df$gene_id)

annotations <- ensembldb::select(
  EnsDb.Rnorvegicus.v79,
  keys = ensembl_ids,
  columns = c("GENEID", "SYMBOL", "GENEBIOTYPE"),
  keytype = "GENEID"
)

res_df <- merge(
  res_df,
  annotations,
  by.x = "gene_id",
  by.y = "GENEID",
  all.x = TRUE
)

res_df$gene_symbol <- ifelse(
  is.na(res_df$SYMBOL),
  res_df$gene_id,
  res_df$SYMBOL
)
```

---

## Step 6. Define Significance for Visualization

```r
res_df$significance <- ifelse(
  res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.3,
  ifelse(res_df$log2FoldChange > 0, "Upregulated", "Downregulated"),
  "Not significant"
)
```

---

## Step 7. Volcano Plot (yy_volcano)

```r
protein_coding_df <- res_df %>%
  filter(GENEBIOTYPE == "protein_coding")

yy_volcano(
  df = protein_coding_df,
  fc_cutoff = 0.3,
  p_cutoff = 0.05,
  top_n_labels = 20,
  col_up = "darkgreen",
  col_down = "purple"
)
```

*(Insert yy_volcano output figure here)*

---

## Step 8. Yin–Yang Plot

### Basic Usage

```r
yinyang(protein_coding_df)
```

### Fully Customized Example

```r
yinyang(
  protein_coding_df,
  title = "Yin–Yang Plot of Protein-Coding DEGs",
  subtitle = "Top up- and down-regulated genes",
  top_n_up = 15,
  top_n_down = 15,
  yin_color = "#90EE90",
  yang_color = "#FFB6C1",
  seed = 42
)
```

*(Insert Yin–Yang plot here)*

---

## Step 9. Pirate Plots for Target Genes

```r
my_genes <- c("Aloxe3", "Cxcl13", "Igkc", "Snai1", "Cd163", "Wnt3")

pirate_plot <- yy_pirateplot(
  dds_object = dds,
  results_df = res_df,
  gene_list = my_genes,
  condition_col = "condition"
)

pirate_plot
```

*(Insert pirate plot figure here)*

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

**Lingzhang Meng, PhD**  
Center for Medical Big Data and Artificial Intelligence Research  
Guangxi Academy of Medical Sciences & People’s Hospital of Guangxi Zhuang Autonomous Region

---

Happy plotting ☯️

