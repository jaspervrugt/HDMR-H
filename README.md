# High Dimensional Model Representation (HDMR) High with Uncorrelated and/or Correlated Input Variables: MATLAB and Python Toolboxes

## Description

High-Dimensional Model Representation High (HDMR-H) for $d > 20$ using B-spline functions for variance-based global sensitivity analysis (GSA) with correlated and uncorrelated inputs. This function uses as input a N x d matrix of N different d-vectors of model inputs (factors/parameters) and a N x 1 vector of corresponding model outputs and returns to the user each factor's first, second, and third order sensitivity coefficient (separated in total, structural and correlative contributions), an estimate of their 95% confidence intervals (from bootstrap method) and the coefficients of the significant B-spline basis functions that govern output, y (determined by an F-test of the error residuals of the HDMR model (emulator) with/without a given first, second and/or third order B-spline). These coefficients define an emulator that can be used to predict the output, y, of the original (CPU-intensive?) model for any d-vector of model inputs. For uncorrelated model inputs (columns of X are independent), the HDMR-H sensitivity indices reduce to a single index (= structural contribution), consistent with their values derived from commonly used variance-based GSA methods.

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

## Version History

* 1.0
    * Initial Release
* 2.0
    * Python implementation
    * New built-in case studies

## Acknowledgments
The Python toolbox is based on the MATLAB code. 

