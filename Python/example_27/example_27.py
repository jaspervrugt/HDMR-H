# ####################################################################### #
#                                                                         #
#  EEEEEE  XX  XX   AAAA   MM   MM  PPPPPP  LL      EEEEEE   222222 77777 #
#  EE       XXXX   AA  AA  MMM MMM  PP  PP  LL      EE       22 22     77 #
#  EEEEE     XX    AA  AA  MMMMMMM  PPPPPP  LL      EEEEE      22     77  #
#  EE       XXXX   AAAAAA  MM   MM  PP      LL      EE        22     77   #
#  EEEEEE  XX  XX  AA  AA  MM   MM  PP      LLLLLL  EEEEEE   222222 77    #
#                                                                         #
# ####################################################################### #

import sys
import os
import numpy as np

# Get the current working directory
current_directory = os.getcwd()
# Go up one directory
parent_directory = os.path.abspath(os.path.join(current_directory, '..'))
# add this to path
sys.path.append(parent_directory)
# Add another directory
misc_directory = os.path.abspath(os.path.join(parent_directory, 'miscellaneous'))
# add this to path
sys.path.append(misc_directory)

from HDMR_H import HDMR_H

# Define the number of samples
N = int(1e5)

# Sample from U(0,1) N times (uniform distribution)
X = np.random.rand(N, 3)

# Define the function f(x)
def f(x):
    return x[:, 0] + 2 * x[:, 1] + 3 * x[:, 2] + x[:, 0] * x[:, 1]

# Compute the function output y
y = f(X)

# HDMR-H options dictionary
options = {
    'graphics': 1,
    'maxorder': 3,
    'maxiter': 100,
    'bf1': 1,
    'm': 2,
    'K': 10,
    'R': 300,
    'method': 1,
    'alfa': 0.01,
    'lambda': 0.10,
    'vartol': 1e-3,
    'refit': 1
}

if __name__ == '__main__':
	# Run the HDMR-H toolbox
	S, Ss, Fx, Em, Xy = HDMR_H(X, y, options)

	# Output the results (optional)
	print(S, Ss, Fx, Em, Xy)
