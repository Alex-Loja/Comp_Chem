#!/usr/bin/env python3
import sys
import numpy as np
import math
import os
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt

######Version 8.0########
# - Additionally produces text file with all coordinates 


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
####

def Parse_Angle(num):
    if "0.00" in num:
        return 0.0
    elif "179" in num:
        return 180.0
    else: return num

def convert_sci(string):
    string_split = string.split('D')
    coefficent = float(string_split[0])
    power = float(string_split[1])
    result = coefficent * (10 ** power)
    return result

def readfile(file):
    Ang_NO = 0; Ang_NN = 0; Ang_H2 = 0; Dih = 0
    Partner_mol = " "
    cont = open(file,'r').readlines()
    for n, i in enumerate(cont):
        if "angno=" in i.casefold():
            Ang_NO = Parse_Angle(i.split("=")[1].split()[0])
        if "angnn=" in i.casefold():
            Ang_NN = Parse_Angle(i.split("=")[1].split()[0])
            Partner_mol = "Ang_NN"
        if "angh2=" in i.casefold():
            Ang_H2 = Parse_Angle(i.split("=")[1].split()[0])
            Partner_mol = "Ang_H2"
        if "dih2=" in i.casefold():
            Dih = Parse_Angle(i.split("=")[1].split()[0])
        check = i.strip().split()
        check = " ".join(check)
        if "DISTANCE ANG" in check:
            end=n   

    R_val, E_val = [], []
    for i in cont[end:]:
        if len(i.split())==0:
            break
        elif "DISTANCE" in i:
            for value, header in enumerate(i.strip().split()):
                if '_cp' in header.casefold() and 'bsse' in file.casefold():
                    E_column = value
                elif '_e' in header.casefold() or 'efinal' in header.casefold() and 'scan' in file.casefold():
                    E_column = value
            continue
        else:
            res = (i.strip().split())
            R_val.append(res[0])
            E_val.append(res[E_column])
    
    val_dict = {
        "energy"  :   E_val,
        "R"         :   R_val,
        "Ang_NO"    :   [],
        "Ang_NN"    :   [],
        "Ang_H2"    :   [],
        "Dih"       :   [],
        "Partner"   :   Partner_mol
    }
    for i in range(len(R_val)):
        val_dict["Ang_NN"].append(Ang_NN)
        val_dict["Ang_H2"].append(Ang_H2)
        val_dict["Ang_NO"].append(Ang_NO)
        val_dict["Dih"].append(Dih)
    return val_dict

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
        
#populate matrix
all_Ang_NO, all_Ang_other, all_Dih, all_R_vals, all_E_vals = [], [], [], [], []
for n, i in enumerate(files_sort):
    res_dict = readfile(i)
    #Build text file
    for x in res_dict["energy"]: all_E_vals.append(x)
    for x in res_dict["R"]: all_R_vals.append(x)
    for x in res_dict["Ang_NO"]: all_Ang_NO.append(x)
    for x in res_dict["Dih"]:  all_Dih.append(x)
    if res_dict["Partner"] == "Ang_NN":
        for x in res_dict["Ang_NN"]:   all_Ang_other.append(x)
    else:
        for x in res_dict["Ang_H2"]:   all_Ang_other.append(x)
    #Pre-build matrix on first file
    if n == 0:
        res_mat = np.zeros((len(res_dict["R"])+1,len(angles)+1))
        nR_vals = len(res_dict["R"])
        for x in range(len(angles)):
            res_mat[0,x+1] = angles[x]
        for x in range(nR_vals):
            res_mat[x+1,0] = res_dict["R"][x]
            if 'D' in res_dict["energy"][x]:
                temp = convert_sci(res_dict["energy"][x])
                res_mat[x+1,1] = temp
            else:
                res_mat[x+1,1] = res_dict["energy"][x]     
    else:
        for p in range(nR_vals):
            if 'D' in res_dict["energy"][p]:
                temp = convert_sci(res_dict["energy"][p])
                res_mat[p+1,n+2] = temp
            else:
                res_mat[p+1,n+1] = res_dict["energy"][p]
            
#Print results
new_name = input.rsplit('_',1)[0] + '_.csv'
np.savetxt(new_name, res_mat, delimiter=', ')

#text-column format
if res_dict["Partner"] == "Ang_NN":
    output = ["R_val, Ang_NO, Ang_NN, Dih, Energy/Ha\n"]
else:
    output = ["R_val, Ang_NO, Ang_H2, Dih, Energy/Ha\n"]
for i in range(len(all_R_vals)):
    output.append(str(all_R_vals[i]) + ", " + str(all_Ang_NO[i]) + ", " +
                  str(all_Ang_other[i]) + ", " + str(all_Dih[i]) + ", " +
                  str(all_E_vals[i]) + '\n')
new_name = new_name.replace(".csv", ".txt")
open(new_name, 'w').writelines(output)





        
        
       
