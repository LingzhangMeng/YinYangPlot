# YinYangPlot

**YinYangPlot** is an R visualization package designed to generate a new class of differential expression plot: the **Yin–Yang plot**. This plot is intended as a powerful *complement* to the traditional volcano plot, helping researchers focus directly on biologically meaningful signals rather than being overwhelmed by thousands of non-informative points.

The package also provides publication-ready volcano plots and enhanced pirate plots, all optimized for RNA-seq / differential expression workflows (e.g. DESeq2 results).

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
devtools::install_github("<your_github_username>/YinYangPlot")
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

- A `ggplot2` object
- Ready for further customization or direct export

*(Insert Yin–Yang plot image here)*

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

*(Insert volcano plot image here)*

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
  genes = c("IL6", "TNF", "CXCL8"),
  group_var = "condition"
)
```

*(Insert combined pirate plot image here)*

---

### `yy_pirateplot_simple()`

A lightweight wrapper for **single-gene** or rapid exploratory plotting.

```r
yy_pirateplot_simple(
  dds = dds,
  gene = "IL6",
  group_var = "condition"
)
```

---

## 🔬 Typical Workflow

1. Perform RNA-seq differential expression analysis (e.g. DESeq2)
2. Annotate genes (symbol, biotype)
3. Generate:
   - `yy_volcano()` for global overview
   - `yinyang()` for focused biological interpretation
   - `yy_pirateplot()` for target gene validation

---

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

