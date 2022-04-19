#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Read out the 31P NMR chemical shift tensor from an ORCA
calculation and then compute isotropy, eta, and zeta.

Run script in the folder containing the op.orca output 
file from your ORCA run.  
'''

import re
from sympy import Matrix, Trace, eye

p = re.compile('\d+.\d+') # regex for finding numbers w/ decimals

def get_cs_tensor(filename='op.orca'):
    
    # read the output file
    f = open(filename, "r")
    lines = f.readlines()
    f.close()
    
    tensor_line = 0
    # find the lines with the chemical shift tensor
    for i, line in enumerate(lines):
        line=line.rstrip()
        if 'Total shielding tensor' in line: # WARNING: takes first line this is true for, if calc mult nmr shifts check which tensor you want to read
            tensor_line = i
            break 
    # check to make sure it found cs tensor
    if tensor_line==0:
        raise ValueError("Could not find chemical shift tensor in file")
    
    # Convert the values in file into list of floats
    cs_tensor = []
    for j in range(1,4):
        line_as_string = lines[tensor_line + j]
        line_as_list = p.findall(line_as_string)
        line_as_floats = [float(k) for k in line_as_list]
        cs_tensor.append(line_as_floats)

    return(cs_tensor)

def convert_cs_tensor(X):
    # Compute traceless symmetric part
    trace = Trace(X).simplify()
    S = 0.5 * (X + X.T) - (trace * eye(3)) / 3

    # Compute lambda values.
    # Sort diagonal values by absolute value,
    # then assign true value to lambda variables.

    # Diagonalize the matrix
    M, D = S.diagonalize()
    # get the values on the diagonal
    lambda_list = [D[i, i] for i in range(0, 3)]
    # sort by absolute value from smallest to largest and assign
    lambda_yy, lambda_xx, lambda_zz = sorted(lambda_list, key=abs)
    # Compute asymmetry parameter
    asym = (lambda_yy - lambda_xx) / lambda_zz

    # Print out values
    print("The isotropy is:", trace/3.)
    print("The anisotropy (zeta) is:", lambda_zz)
    print("The asymmetry parameter (eta) is:", asym)
    print(trace/3.)
    print(lambda_zz)
    print(asym)

if __name__ == '__main__':
    C = get_cs_tensor() # list of lists
    M = Matrix(C)       # sympy Matrix
    convert_cs_tensor(M)
