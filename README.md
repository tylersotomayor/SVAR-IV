# IV-SVAR Tutorial: Identifying Monetary Policy Shocks with External Instruments

## 📌 Introduction
Welcome to this step-by-step tutorial on implementing **Instrumental Variable Structural Vector Autoregression (IV-SVAR)** in MATLAB. I started this tutorial in Spring of 2024 to teach my classmates how to use SVAR in MATLAB. Disclaimer: I am in the process of updating this tutorial to be even better in 2025. This tutorial is designed for beginners with **minimal background in MATLAB, coding, or SVAR/macroeconometrics**. By the end of this tutorial, you will understand how to:

✅ Load and preprocess macroeconomic data
✅ Estimate a **Vector Autoregression (VAR)** model
✅ Identify **monetary policy shocks** using external instruments
✅ Bootstrap confidence intervals to assess robustness
✅ Compare alternative identification strategies
✅ Replicate empirical findings from **Gertler & Karadi (2015)** and **Miranda-Agrippino & Ricco (2021)**

# IV-SVAR Tutorial: Identifying Monetary Policy Shocks with External Instruments

## 📂 Repository Structure
This repository contains:

### **1. Main MATLAB Scripts**
These files contain the core implementation of the IV-SVAR model and related econometric procedures:
- `IV_VAR.m` → Main script implementing IV-SVAR
- `VAR.m` → Function for estimating a standard VAR
- `woldirf.m` → Computes Wold IRFs
- `choleskyIRF.m` → Implements Cholesky IRFs
- `bootstrapChol.m`, `bootstrapIV_corrected.m` → Bootstrapping confidence bands
- `plotirf_partial.m`, `plot_vardec.m` → Functions for plotting results
- `myols.m` → Helper function for OLS regression
- `remove_bias.m` → Bias correction function

### **2. Data Files**
These datasets are necessary to replicate the results:
- `Proxydata.xlsx` → Macro dataset (FRED data)
- `instruments_GK2015.csv` → External instruments for identification
- `OilDataM.mat`, `OilSurprisesMLog.mat`, `OilSurprisesMLogControl.mat` → Used for **Känzig (2021) oil supply news replication**

### **3. Documentation**
- `README.md` → This tutorial walkthrough
- `LICENSE` (Optional) → Defines open-source terms (if applicable)

## 🛠 Prerequisites
This tutorial assumes no prior coding experience. However, you will need:
- MATLAB installed on your computer
- **Basic knowledge** of economics (inflation, interest rates, monetary policy)
- An interest in **macroeconometric modeling**

If you are completely new to MATLAB, consider reviewing **[MATLAB Onramp](https://matlabacademy.mathworks.com/details/matlab-onramp/getting-started)** before starting.

## 📥 Step 1: Setting Up Your MATLAB Environment
1. Open MATLAB
2. Set your **working directory** to the folder containing this repository:
   ```matlab
   cd 'path_to_your_repository'
   ```
3. Ensure that `Proxydata.xlsx` and `instruments_GK2015.csv` are in your working directory.

## 📁 Step 2: Navigating the Files
Once you have cloned the repository, here’s how to navigate the key files:
- **Start with `IV_VAR.m`**: This is the main script that runs the IV-SVAR model.
- **Data files (`.xlsx`, `.mat`, `.csv`)**: These are preloaded datasets used for estimation.
- **Helper functions (`.m` files)**: These perform operations like estimating VARs, computing impulse responses, and bootstrapping.
- **Figures will be generated**: Running the script will create IRF and variance decomposition plots automatically.

To execute the full analysis, open MATLAB and run:
```matlab
run('IV_VAR.m')
```
This will generate the impulse response functions and other outputs.

## 📊 Step 3: Estimating a Cholesky VAR
Before implementing IV-SVAR, we estimate a **standard Cholesky VAR** as a benchmark.
```matlab
p = 12;  % Lag length
c = 1;   % Include constant term
[beta, residuals] = VAR(finaldata, p, c);
```

### Computing Impulse Responses
```matlab
hor = 48; % Forecast horizon
wold = woldirf(beta, c, p, hor);
sigma = (residuals' * residuals) ./ (T - 1 - p - N*p);
S = chol(sigma, 'lower');
cholirf = choleskyIRF(wold, S);
```

## 🎯 Step 4: IV-SVAR Identification (Using External Instruments)
We now use an **instrumental variable** to isolate monetary policy shocks.

### Load the Instrument
```matlab
ivGK = csvread('instruments_GK2015.csv',1);
ff4 = ivGK(127:end,4);  % January 1990 - June 2012
```

### IV Estimation
```matlab
u_p = residuals(:,3);  % Policy residuals
ff4_final = ff4(1:end-30);
u_p_final = u_p(193:end);
u_q_final = residuals(193:end,[1,2,4]);
[~, uhat, ~] = myols(u_p_final, ff4_final,0);
sq_sp = (uhat'*uhat)\uhat'*u_q_final;
s = [sq_sp(1) sq_sp(2) 1 sq_sp(3)]';
```

## 🔄 Step 5: Bootstrapping Confidence Intervals
```matlab
nboot = 1000;
prc = 68;
[upper, lower, meanirf, medianirf] = bootstrapIV_corrected(finaldata,p,c,beta,residuals,ff4, nboot, 2000, [1,240], [193,432], 3, hor, prc);
```

## 📖 Step 6: Visualizing and Interpreting Results
```matlab
plotirf_partial(meanirf,upper,lower,{'IP', 'CPI', 'SR', 'EBP'}, 'Monetary Policy Shock', prc);
```

## 📚 Further Reading & References
- **Gertler & Karadi (2015)**: "Monetary Policy Surprises, Credit Spreads, and Economic Activity"
- **Miranda-Agrippino & Ricco (2021)**: "The Transmission of Monetary Policy Shocks"
- **Stock & Watson (2012)**: "Disentangling the Channels of Monetary Policy"

## 🛠 Contributing
If you find any issues or have suggestions, feel free to **open a GitHub issue** or **submit a pull request**.

## 📧 Contact
For questions, email me at: **[your email]**

---

This tutorial is designed to be an **accessible introduction** to IV-SVAR. If you're stuck at any step, feel free to reach out!