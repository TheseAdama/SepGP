# 📦 SepGP
The **SepGP** package provides tools for the modeling and prediction of spatio-temporal functions defined as:

$$
f : (x, t) \in \mathbb{R}^{d} \times \mathbb{R} \ \longmapsto \ f(x, t) \in \mathbb{R}
$$

In this framework, the function $f(x, t)$ is observed at discrete time points $t_1, \dots, t_{N_t}$ for all spatial locations $x \in \mathbb{R}^d$.

This package is designed to efficiently handle and analyze spatio-temporal data where observations are collected over space and time.

## 📥 Installation

You can install the latest version of the package manually or directly from GitHub.

### Option 1: Install from GitHub

Make sure you have the `devtools` package installed, then use:

```r
install.packages("devtools")
devtools::install_github("TheseAdama/SepGP")
```
### Option 2: Manual Installation (Download & Install ZIP)

1. **Download the ZIP or TAR.GZ file**
   
   Download the latest version of the package in ZIP or TAR.GZ format.

   - For **Windows**: `SepGP_x.y.z.zip`
   - For **Linux/macOS**: `SepGP_x.y.z.tar.gz`

3. **Install the package manually in R**

   Open your R session and run one of the following commands, replacing the file path with where you downloaded the archive:

   - **On Windows**:
     ```r
     install.packages("path/to/SepGP_x.y.z.zip", repos = NULL, type = "win.binary")
     ```

   - **On Linux/macOS**:
     ```r
     install.packages("path/to/SepGP_x.y.z.tar.gz", repos = NULL, type = "source")
     ```

---

After installation, load the package:

```r
library(SepGP)
```
