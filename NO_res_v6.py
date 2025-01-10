#!/usr/bin/env python3
import sys
import numpy as np
import math
import os
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt

######Version 6.0########
# - Added --help function


class bcolors :
    Yes = '\033[96m'
    No = '\033[0m'
###Help 
if sys.argv[1].casefold() == 'help' or sys.argv[1].casefold() == '--help':
    print(bcolors.Yes + '-- NO_res Script Version 6.0 --')
    print('To run please provide common prefix of file name')
    print('Ensure all .out files are present in current folder')
    print('Script Currently only works for molpro files.')
    print('Tested to molpro v2023.3')
    print(bcolors.No)
    quit()
####å


def convert_sci(string):
    string_split = string.split('D')
    coefficent = float(string_split[0])
    power = float(string_split[1])
    result = coefficent * (10 ** power)
    return result

def readfile(file,start):
    cont = open(file,'r').readlines()
    for i, n in enumerate(cont):
        check = n.strip().split()
        check = " ".join(check)
        if "DISTANCE ANG" in check:
            end=i   

    R_val, E_val = [], []
    for i in cont[end:]:
        if len(i.split())==0:
            break
        elif "DISTANCE" in i:
            continue
        else:
            res = (i.strip().split())
            R_val.append(res[0])
            E_val.append(res[2])

    if start == True: return E_val, R_val
    else: return E_val

def CSV_to_text_output(matrix):
    if len(matrix[0,:]) == 2:
        quit()
    else:
        #build matrix
        x_vals = matrix[1:,0]
        len_x = len(x_vals)
        y_vals = matrix[0,1:]
        len_y = len(y_vals)
        new_len = len_x * len_y
        new_mat = np.zeros((new_len,3))
        # input values
        x = 0
        y = 0
        while x < new_len:
            p = 0
            while p < int(len_x):
                new_mat[x,0] = x_vals[p]
                new_mat[x,1] = y_vals[y]
                new_mat[x,2] = matrix[p+1,y+1]
                p += 1
                x += 1
            y += 1
    return new_mat
    
#Select file
input = sys.argv[1]
cwd = os.getcwd()
files, files_sort, angles = [], [], []
angles = []
for i in os.listdir(cwd):
    if '.out_' in i:
        continue
    elif '.out' in i:
        if input in i:
            files.append(i)
            angles.append(int((i.split('.')[0]).split('_')[-1]))
angles.sort()

#Sort files in ascending order
for i in angles:
    for n in files:
        if ('_' + str(i) + '.') in n:
            files_sort.append(n)

#Pre-build results matrix
#(requires getting results from the first file first)
res_0, R_vals = readfile(files_sort[0],True)
res_mat = np.zeros((len(R_vals)+1,len(angles)+1))
nR_vals = len(R_vals)
for i in range(len(angles)):
    res_mat[0,i+1] = angles[i]
for i in range(nR_vals):
    res_mat[i+1,0] = R_vals[i]
    if 'D' in res_0[i]:
        temp = convert_sci(res_0[i])
        res_mat[i+1,1] = temp
    else:
        res_mat[i+1,1] = res_0[i]

#populate matrix
for n, i in enumerate(files_sort[1:]):
    res = readfile(i,False)
    for p in range(nR_vals):
        if 'D' in res[p]:
            temp = convert_sci(res[p])
            res_mat[p+1,n+2] = temp
        else:
            res_mat[p+1,n+2] = res[p]

#Print results
new_name = input.rsplit('_',1)[0] + '.csv'
np.savetxt(new_name, res_mat, delimiter=', ')
text_format = CSV_to_text_output(res_mat)
np.savetxt(new_name.replace('.csv','.txt'), text_format, delimiter=' ')





        
        
       
