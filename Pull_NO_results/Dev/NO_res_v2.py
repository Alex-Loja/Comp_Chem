#!/usr/bin/env python3
import sys
import numpy as np
import math
import os
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt

#Select file
file = sys.argv[1]
cont = open(file,'r').readlines()

#Pull file contents
for i, n in enumerate(cont):
    check = n.strip().split()
    check = " ".join(check)
    if "DISTANCE ANG" in check:
        end=i
    
results = []
for i in cont[end:]:
    results.append(i.strip().split())
    if len(i.split())==0:
        break

final = ",".join(results[1])+','
for i in results[2:]:
    if len(i) == 0:
        break
    else:
        final += (",".join(i))+','

final = final.split(',')

#####Pre-build matrix
#total number of values
units = int(len(final))-1
#quantity of lengths and angles
lengths, angles = [], []
a = float(final[0])
lengths.append(a)
n_lenghts = 1
b = int(float(final[1]))
angles.append(b)
n_angles = 1

for i in range(3,units,3):
    x = float(final[i])
    while x < a:
        n_lenghts += 1
        lengths.append(x)
        a = x
    y = int(float(final[i+1]))
    if y != angles[-1]:
        angles.append(y)
        n_angles += 1
    

#Build matrix
res_mat = np.zeros((n_lenghts+1,n_angles+1))
for i,n in enumerate(lengths):
    res_mat[i+1,0] = n
for i,n in enumerate(angles):
    res_mat[0,i+1] = n

#A = column, B = row
A = 1
B = 1
for i in range(2,units,3):
    if 'D' in final[i]:
        val = final[i].split('D')
        value = float(val[0])*(10**(float(val[1])))
    else:
        value = final[i]
    res_mat[B,A] = value
    B += 1
    if B == n_lenghts +1:
        B = 1
        A += 1

#Print at specific decimals
# print(np.around(res_mat,5))

# dataframe = pd.DataFrame(res_mat)
# dataframe = pd.DataFrame.to_csv(dataframe)
new_name = file.replace('.out','.csv')
# open((new_name),'w').writelines(dataframe)

np.savetxt(new_name, res_mat, delimiter=', ')




        
        
       
