#!/usr/bin/env python3
import sys
import numpy as np
import math
import os
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt

######Version 10.0########
# - Takes values more directly from where F12 values are printed if F12 is used
# - Takes into account default F12 requests both F12a and F12b
# - Will only take both F12a and F12b values
# - Attempts to smooth code by introducing classes


class bcolors :
    Yes = '\033[96m'
    No = '\033[0m'
###Help 
if sys.argv[1].casefold() == 'help' or sys.argv[1].casefold() == '--help':
    print(bcolors.Yes + '-- NO_res Script Version 10 --')
    print('To run please provide common prefix of file name')
    print('Ensure all .out files are present in current folder')
    print('Script Currently only works for molpro files.')
    print('Tested to molpro v2024.1')
    print(bcolors.No)
    quit()
####

def Parse_Angle(num):
    if "0.00" in num:
        return 0.0
    if "1e-05" in num:
        return 0.0
    if "1e-06" in num:
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
 
class CCSD_Scan:
    """Handling an output file of NO scans"""
    
    def __init__(self, name):
        """ initialise type of file """
        self.name = name
        self.angle = name.strip().split('_')[-1].split('.')[0]
        self.read()
        
    def read(self):
        """reads the file"""
        contents = open(self.name,'r').readlines() #Read in the entire file
        r_vals, e_vals = [], []                 #Prep arrays for lengths and energies
        self.partner = ' '                         #Prep variable for name of the partner molecule
        collect_energy_flag = False                #Prep flag for picking up energies
        for line in contents:
            #Section for Collecting F12a energies
            if collect_energy_flag == True:        #Energies are taken as the line following the flag
                collect_energy_flag = False
                e_vals.append(line.strip().split(' ')[-1])
            if "RHF-RCCSD-T energy" in line:   #Flag to show energy is on the following line
                collect_energy_flag = True
            #Section to collect lengths
            if "setting r " in line.casefold():
                value = line.strip().split("=")[1].strip()
                if 'angstrom' in value.casefold(): continue #Ensure only numerical values are taken
                if float(value) > 50: continue              #Ensure the BSSE 100Ang isn't picked up
                if len(r_vals) == 0: r_vals.append(value)   #Negates the double initial setting R for BSSE
                if value == r_vals[0]: continue             #Makes sure lengths appear in numerical ordering
                if float(value) >= float(r_vals[-1]): 
                    print('error, inconsistent R values, check output'); quit()
                r_vals.append(value)
            #Section for angle collections
            if "angno=" in line.casefold():
                self.ang_no = Parse_Angle(line.split("=")[1].split()[0])
            if "angn2=" in line.casefold() or "angnn=" in line.casefold():
                self.ang_other = Parse_Angle(line.split("=")[1].split()[0])
                if self.partner == '': self.partner = 'n2'
            if "angh2=" in line.casefold():
                self.ang_other = Parse_Angle(line.split("=")[1].split()[0])
                if self.partner == '': self.partner == 'h2'
            if "dih2=" in line.casefold():
                self.ang_dih = Parse_Angle(line.split("=")[1].split()[0])
        self.e_vals = e_vals
        self.r_vals = r_vals
        self.n_r_vals = len(r_vals)
        
class CCSD_BSSE(CCSD_Scan):
    """For handeling BSSE files specifically"""
    def __init__(self, name):   
        super().__init__(name)  #Take overlapping class values from the explicit correlation scan class
        self.cp = float(self.e_vals[0])
        energies = self.e_Fa_vals[1:]
        new_e_vals = []
        for energy in energies: new_e_vals.append(float(energy)-self.cp)
        self.e_vals = new_e_vals
    
class Expcorr_Scan:
    """Handling an output file of NO scans"""
    
    def __init__(self, name):
        """ initialise type of file """
        self.name = name
        self.angle = name.strip().split('_')[-1].split('.')[0]
        self.read()
        
    def read(self):
        """reads the file"""
        contents = open(self.name,'r').readlines() #Read in the entire file
        r_vals, e_Fa_vals, e_Fb_vals = [], [], []                    #Prep arrays for lengths and energies
        self.partner = ''                         #Prep variable for name of the partner molecule
        collect_energy_Fa_flag = False                #Prep flag for picking up energies
        collect_energy_Fb_flag = False
        for line in contents:
            #Section for Collecting F12a energies
            if collect_energy_Fa_flag == True:        #Energies are taken as the line following the flag
                collect_energy_Fa_flag = False
                e_Fa_vals.append(line.strip().split(' ')[-1])
            if "RHF-RCCSD-T-F12a energy" in line:   #Flag to show energy is on the following line
                collect_energy_Fa_flag = True
            #Section for Collecting F12b energies
            if collect_energy_Fb_flag == True:        #Energies are taken as the line following the flag
                collect_energy_Fb_flag = False
                e_Fb_vals.append(line.strip().split(' ')[-1])
            if "RHF-RCCSD-T-F12b energy" in line:   #Flag to show energy is on the following line
                collect_energy_Fb_flag = True
            #Section to collect lengths
            if "setting r " in line.casefold():
                value = line.strip().split("=")[1].strip()
                if 'angstrom' in value.casefold(): continue #Ensure only numerical values are taken
                if float(value) > 50: continue              #Ensure the BSSE 100Ang isn't picked up
                if len(r_vals) == 0: r_vals.append(value)   #Negates the double initial setting R for BSSE
                if value == r_vals[0]: continue             #Makes sure lengths appear in numerical ordering
                if float(value) >= float(r_vals[-1]): 
                    print('error, inconsistent R values, check output'); quit()
                r_vals.append(value)
            #Section for angle collections
            if "angno=" in line.casefold():
                self.ang_no = Parse_Angle(line.split("=")[1].split()[0])
            if "angn2=" in line.casefold() or "angnn=" in line.casefold():
                self.ang_other = Parse_Angle(line.split("=")[1].split()[0])
                if self.partner == '': self.partner = 'n2'
            if "angh2=" in line.casefold():
                self.ang_other = Parse_Angle(line.split("=")[1].split()[0])
                if self.partner == '': self.partner = 'h2'
            if "dih2=" in line.casefold():
                self.ang_dih = Parse_Angle(line.split("=")[1].split()[0])
        self.e_Fa_vals = e_Fa_vals
        self.e_Fb_vals = e_Fb_vals
        self.r_vals = r_vals
        self.n_r_vals = len(r_vals)
    
class Expcorr_BSSE(Expcorr_Scan):
    """For handeling BSSE files specifically"""
    def __init__(self, name):   
        super().__init__(name)  #Take overlapping class values from the explicit correlation scan class
        self.cp_Fa, self.e_Fa_vals = self.Calc_CP('Fa')
        self.cp_Fb, self.e_Fb_vals = self.Calc_CP('Fb')
        
    def Calc_CP(self, type):
        if type == 'Fa': 
            cp = float(self.e_Fa_vals[0])
            energies = self.e_Fa_vals[1:]
        else: 
            cp = float(self.e_Fb_vals[0])
            energies = self.e_Fb_vals[1:]
        new_e_vals = []
        for energy in energies: new_e_vals.append(float(energy)-cp)
        return cp, new_e_vals
        
#Select file
input = sys.argv[1]
cwd = os.getcwd()
files, files_sort, angles = [], [], []
for i in os.listdir(cwd):   #Pulls file names and angles from said names
    if '.out_' in i:
        continue
    elif '.out' in i:
        if input in i:
            files.append(i)
            angles.append(int((i.split('.')[0]).split('_')[-1]))
angles.sort()       #Sort angles into numerical order
for i in angles:    #Sort files in ascending order by assigning angles
    for n in files:
        if ('_' + str(i) + '.') in n:
            files_sort.append(n)

method = ''
all_Ang_NO, all_Ang_other, all_Dih, all_R_vals = [], [], [], []
all_E_vals, all_Fa_E_vals, all_Fb_E_vals = [], [], []
for n, i in enumerate(files_sort):
    #Need to read the initial file to figure out if its an F12 calculation or not
    if n == 0:
        with open(i) as file_object:
            for x, line in enumerate(file_object):
                if x == 200: break
                if 'f12' in line.casefold(): method = 'F12'
                
    #Pull file info depending on the type of file it is
    if method == 'F12':
        if 'bsse' in i.casefold(): file_info = Expcorr_BSSE(i)
        elif 'scan' in i.casefold(): file_info = Expcorr_Scan(i)
        else: quit()
    else: 
        if 'bsse' in i.casefold(): file_info = CCSD_BSSE(i)
        elif 'scan' in i.casefold(): file_info = CCSD_Scan(i)
        else: quit()
        
    
        
    
    #Construct empty matrix
    if n == 0:
        n_ang_vals = len(angles)
        res_mat = np.zeros((file_info.n_r_vals+1, n_ang_vals+1))
        for x in range(n_ang_vals): res_mat[0 , x+1] = angles[x]
        if method == 'F12':
            Fa_mat = res_mat.copy()
            Fb_mat = res_mat.copy()
        for x in range(file_info.n_r_vals): 
            if method == 'F12':
                Fa_mat[x+1 , 0] = file_info.r_vals[x]
                Fb_mat[x+1 , 0] = file_info.r_vals[x]
            else:
                res_mat[x+1 , 0] = file_info.r_vals[x]
        
            
    if method == 'F12':   #Populate matrix - F12
        for x in range(file_info.n_r_vals):
            Fa_mat[x+1 , n+1] = file_info.e_Fa_vals[x]
            Fb_mat[x+1 , n+1] = file_info.e_Fb_vals[x]
            #Add file-consistent info to text-output arrays
            all_Ang_NO.append(file_info.ang_no)
            all_Ang_other.append(file_info.ang_other)
            all_Dih.append(file_info.ang_dih)
            all_R_vals.append(file_info.r_vals[x])
            all_Fa_E_vals.append(file_info.e_Fa_vals[x])
            all_Fb_E_vals.append(file_info.e_Fb_vals[x])

    else:   #Populate matrix - non-F12
        for x in range(file_info.n_r_vals):
            res_mat[x+1 , n+1] = file_info.e_vals[x]
            #Add file-consistent info to text-output arrays
            all_Ang_NO.append(file_info.ang_no)
            all_Ang_other.append(file_info.ang_other)
            all_Dih.append(file_info.ang_dih)
            all_R_vals.append(file_info.r_vals[x])
            all_E_vals.append(file_info.e_vals[x])

if method == 'F12':
    #Print matrix results
    new_name_Fa = input.rsplit('_',1)[0] + '_Fa.csv'
    np.savetxt(new_name_Fa, Fa_mat, delimiter=', ')
    new_name_Fb = input.rsplit('_',1)[0] + '_Fb.csv'
    np.savetxt(new_name_Fb, Fb_mat, delimiter=', ')
    #Print text results
    if file_info.partner == "n2":
        output = ["R_val, Ang_NO, Ang_NN, Dih, Energy_Fa/Ha, Energy_Fb/Ha\n"]
    elif file_info.partner == "h2":
        output = ["R_val, Ang_NO, Ang_H2, Dih, Energy_Fa/Ha, Energy_Fb/Ha\n"]
    else:
        print('partner molecule type not supported'); quit()
    for i in range(len(all_R_vals)):
        output.append(str(all_R_vals[i]) + ", " + str(all_Ang_NO[i]) + ", " +
                    str(all_Ang_other[i]) + ", " + str(all_Dih[i]) + ", " +
                    str(all_Fa_E_vals[i]) + ", " + str(all_Fb_E_vals[i]) + '\n')
    new_name = new_name_Fa.replace("_Fa.csv", ".txt")
    open(new_name, 'w').writelines(output)

    
else:
    #Print matrix results
    new_name = input.rsplit('_',1)[0] + '_.csv'
    np.savetxt(new_name, res_mat, delimiter=', ')
    #Print text results
    if file_info.partner == "n2":
        output = ["R_val, Ang_NO, Ang_NN, Dih, Energy/Ha\n"]
    elif file_info.partner == "h2":
        output = ["R_val, Ang_NO, Ang_H2, Dih, Energy/Ha\n"]
    else:
        print('partner molecule type not supported'); quit()
    for i in range(len(all_R_vals)):
        output.append(str(all_R_vals[i]) + ", " + str(all_Ang_NO[i]) + ", " +
                    str(all_Ang_other[i]) + ", " + str(all_Dih[i]) + ", " +
                    str(all_E_vals[i]) + '\n')
    new_name = new_name.replace(".csv", ".txt")
    open(new_name, 'w').writelines(output)
quit()