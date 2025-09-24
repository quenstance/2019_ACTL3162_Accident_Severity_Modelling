# Accident Severity Modelling & Ruin Probability

### Overview

This project, an individual assignment for the **ACTL3162 General Insurance course**, tackles two fundamental areas of actuarial science: **accident severity modeling** and the analysis of **ruin probability**. The primary goal was to apply statistical and actuarial theories to analyze a dataset of motor insurance claims. All analysis, modeling, and computations were performed using the **R programming language**.

For more information, visit the [project web page](https://quenstance.pages.dev/projects/accident-severity-model/).

---

### How to Use

To replicate this project and run the analysis yourself, you will need to have **R** and **RStudio** installed.

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/quenstance/2019_ACTL3162_Accident_Severity_Modelling.git](https://github.com/quenstance/2019_ACTL3162_Accident_Severity_Modelling.git)
    ```
2.  **Navigate to the project directory and set your working directory in RStudio:**
    ```R
    #Set working directory to the location of the cloned repository
    setwd("C:/path/to/your/cloned/repo") 
    ```
3.  **Install the required R packages:** Run the following commands in the R console to install the necessary libraries:
    ```R
    install.packages("fitdistrplus")
    install.packages("actuar")
    install.packages("moments")
    install.packages("rSymPy")
    install.packages("rJava")
    install.packages("rootSolve")
    ```
4.  **Load packages:** After installation, load the libraries by running:
    ```R
    library(fitdistrplus)
    library(actuar)
    library(moments)
    library(rSymPy)
    library(rJava)
    library(rootSolve)
    ```
5.  **Run the script:** Execute the R script files in the correct order as described in the report to reproduce the analysis and findings. The `data.csv` file should be included in the repository for the script to run.

---

### Methodology

The methodology for this project was structured into two main tasks:

* **Task 1: Accident Severity Modeling**: This involved using **Maximum Likelihood Estimation (MLE)** to fit various distributions (e.g., Gamma, Lognormal, and Burr) to claims data. The model selection process was comprehensive, combining a graphical approach (using Cullen and Frey graphs, Q-Q plots, and P-P plots) with statistical tests like the Anderson-Darling, Kolmogorov-Smirnov, and Cramér-von Mises tests, along with the AIC and BIC information criteria.
* **Task 2: Ruin Theory**: This task involved applying ruin theory to a simple surplus model, examining scenarios with and without insurance, as well as with both proportional and non-proportional reinsurance. A key part of this task also involved implementing the Panjer recursion method to perform convolution and approximate the distribution of aggregated claims.

---

### Key Finding

The project yielded several important findings:

* The **Gamma distribution** was identified as the most suitable model for the claims data, outperforming other distributions in statistical tests and effectively capturing the crucial tail risk of the distribution.
* A deeper analysis of the more complex Burr distribution revealed a high correlation in bootstrapped parameters, suggesting that its additional parameters might not offer significant value.
* The analysis demonstrated that non-proportional reinsurance was more effective than proportional reinsurance at increasing the adjustment coefficient, thus improving solvency.
* The implementation of Panjer recursion successfully produced numerical results that were very close to the true values, validating the chosen approximation method.

---

### Limitation

The project was subject to several constraints:

* The dataset consisted of just over 1,000 claims from a single motor insurance product, which may limit the generalizability of the findings.
* The process of finding appropriate initial conditions for the MLE optimization was noted as being challenging.
* The continuous model for ruin probability proved to be computationally intensive and susceptible to numerical instability when calculating values for larger time horizons.

---

### Author

* **Quenstance Lau**
* **Web Portfolio**: [https://quenstance.pages.dev/projects/accident-severity-model/](https://quenstance.pages.dev/projects/accident-severity-model/)
* **LinkedIn**: [https://www.