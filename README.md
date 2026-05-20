# ImmuneSigR 🧬

[![CRAN status](https://www.r-pkg.org/badges/version/ImmuneSigR)](https://CRAN.R-project.org/package=ImmuneSigR)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/grand-total/ImmuneSigR)](https://CRAN.R-project.org/package=ImmuneSigR)
**ImmuneSigR** is a comprehensive and lightweight R package designed for immune cell signature retrieval, GMT file management, and dependency-free single-cell expression scoring. It provides a robust, literature-derived database of immune cell markers to streamline your single-cell RNA-seq downstream analysis.

## 🌟 Key Features

* **🔍 Signature Retrieval:** Easily search and retrieve highly curated immune cell signatures based on specific cell types or literature sources (PMIDs).
* **📂 GMT File Management:** Natively export signatures to standard Gene Matrix Transposed (GMT) formats, or create custom GMT files for downstream tools like GSEA.
* **📊 Single-Cell Scoring:** Dependency-free scoring of gene-by-cell expression matrices using rank-based or mean-expression methods.
* **📖 Literature Backed:** Core signatures are rigorously curated from landmark papers, including the lung cell atlas (Travaglini et al., 2020) and pan-cancer B cell atlas (Fitzsimons et al., 2024).

## 🚀 Installation

You can install the released version of **ImmuneSigR** directly from [CRAN](https://CRAN.R-project.org/package=ImmuneSigR):

```R
install.packages("ImmuneSigR")
```

Alternatively, you can install the latest development version from GitHub:

```R
# install.packages("devtools")
devtools::install_github("YingYanLaing/ImmuneSigR")
```

## 💡 Quick Start

```R
library(ImmuneSigR)

# 1. Search for specific immune cell signatures
b_cell_sigs <- Search_ImmuneSigR(keyword = "B cell")

# 2. Export signatures to a GMT file
Export_ImmuneSigR_GMT(cell_type = "B cell", output_dir = tempdir())

# 3. Score a single-cell expression matrix
# Assuming 'expr_mat' is your normalized expression matrix
scores <- Score_ImmuneSigR(expr_mat, target_cells = "Plasma cell", method = "mean")
```

For detailed usage and real-world validation (including UMAP visualizations), please refer to the package [Vignette](https://www.google.com/search?q=https://CRAN.R-project.org/package%3DImmuneSigR/vignettes/ImmuneSigR-tutorial.html).

## 📝 License

This project is licensed under the GPL-3 License.
