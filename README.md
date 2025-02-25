## High Dimensional Model Representation (HDMR) High with Uncorrelated and/or Correlated Input Variables: MATLAB and Python Toolboxes

## Description

High-Dimensional Model Representation High (HDMR-H) for $d > 20$ using B-spline functions for variance-based global sensitivity analysis (GSA) with correlated and uncorrelated inputs. This function uses as input a $N \times d$ matrix **X** of $N$ different $d$-vectors of model inputs (factors/parameters) and a $N \times 1$ vector of corresponding model outputs **y** and returns to the user each factor's first, second, and third order sensitivity coefficient (separated in total, structural and correlative contributions), an estimate of their 95% confidence intervals (from bootstrap method) and the coefficients of the significant B-spline basis functions that govern output, $y$ (determined by an F-test of the error residuals of the HDMR model (emulator) with/without a given first, second and/or third order B-spline). These coefficients define an emulator that can be used to predict the output, $y$, of the original (CPU-intensive?) model for any $d$-vector of model inputs. For uncorrelated model inputs (columns of **X** are independent), the HDMR-H sensitivity indices reduce to a single index (= structural contribution), consistent with their values derived from commonly used variance-based GSA methods.

## Getting Started

### Installing: MATLAB

* Download and unzip the zip file 'MATLAB_code_HDMR_H_V2.0.zip' in a directory 'HDMR-H'
* Add the toolbox to your MATLAB search path by running the script 'install_HDMR_H.m' available in the root directory
* You are ready to run the examples.

### Executing program

* After intalling, you can simply direct to each example folder and execute the local 'example_X.m' file.
* Please make sure you read carefully the instructions (i.e., green comments) in 'install_HDMR_H.m'
  
### Installing: Python

* Download and unzip the zip file 'Python_code_HDMR_H_V2.0.zip' to a directory called 'HDMR-H'.

### Executing program

* Go to Command Prompt and directory of example_X in the root of 'HDMR-H'
* Now you can execute this example by typing 'python example_X.py'
* Instructions can be found in the file 'HDMR_H.py' 
  
## Authors

* Vrugt, Jasper A. (jasper@uci.edu)

## Literature
1. Gao, Y., C. Guilloteau, E. Foufoula-Georgiou, C. Xu, X. Sun, and J.A. Vrugt (2024), Soil moisture-cloud-precipitation feedback in the lower atmosphere from functional decomposition of satellite observations. Geophysical Research Letters, 51, e2024GL110347, https://doi.org/10.1029/2024GL110347
2. Gao, Y., A. Sahin, and J.A. Vrugt (2023), Probabilistic Sensitivity Analysis With Dependent Variables: Covariance-Based Decomposition of Hydrologic Models, Water Resources Research, 59 (4), e2022WR0328346, https://doi.org/10.1029/2022WR032834
3. Chastaing, G., F. Gamboa, and C. Prieur (2012), Generalized Hoeffding-Sobol decomposition for dependent variables - application to sensitivity analysis, _Electronic Journal of Statistics_, 6, pp. 2420–2448, https://doi.org/10.1214/12-EJS749
4. Li, G. H. Rabitz, P.E. Yelvington, O.O. Oluwole, F. Bacon, C.E. Kolb, and J. Schoendorf (2010), Global sensitivity analysis for systems with independent and/or correlated inputs, _Journal of Physical Chemistry A_, 114 (19), pp. 6022-6032

## Version History

* 1.0
    * Initial Release
* 2.0
    * Python implementation
    * New built-in case studies

## Built-in Case Studies

1. Example 1: Multivariate normal benchmark study
2. Example 2: Multivariate normal benchmark study with correlated variables, $\Sigma \neq I$dentity matrix
3. Example 3: Multivariate normal benchmark study with correlated variables, $\Sigma \neq I$dentity matrix
4. Example 4: Ishigami function with uncorrelated variables
5. Example 5: Ishigami function with correlated variables $\Sigma \neq I$dentity matrix
6. Example 6: Function with $\Sigma \neq I$dentity matrix
7. Example 7: Sobol $g$ function
8. Example 8: Multivariate normal with correlated variables
9. Example 9: Example 1 of Chastaing et al. (2012)
10. Example 10: Example 2 of Chastaing et al. (2012)
11. Example 21: Soil temperature modeling
12. Example 22: Rainfall runoff modeling: hmodel
13. Example 23: Rainfall runoff modeling: SAC-SMA model
14. Example 24: One-predator-one-prey model
15. Example 25: Two-predators-two-preys model
16. Example 26: Two-predators-two-preys model with measured data
17. Example 27: Simple multivariate function

## Acknowledgments

The Python toolbox is a translation of the MATLAB code. 

