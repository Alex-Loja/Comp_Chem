#!/usr/bin/env python3
import sys
import numpy as np
import math
import os
import pandas as pd
import re
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import colors as c
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from math import pi
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import matplotlib.ticker as ticker
from matplotlib import gridspec
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import FormatStrFormatter
from datetime import datetime
from matplotlib.patches import Arc

##### NO_plot > Polar_Heatmap_v007_unix.py #####
### By Alexandre De Matos Loja
### 
### Script for plotting NO + N2 PES data as a polar heatmap
### Script will also convert from Ha units to selected options given (Ha, cm^-1, eV)
### Logging function is present to keep track of what has been plotted and how
###
### Current Version : 007
### Changes:
###     change starting point of plot - Make it so R starts at 0
###     Tick lines more consistent
###     Font of R values is reduced
###     Ticks for R values are more reasonable (whole numbers)

### To run ###
## Have a parent folder containing subfolders for scans and BSSE
## run from within the Scans folder and make sure the BSSE folders are in the parent directory

class bcolors : 
    DIRC = '\033[96m'
    ENDC = '\033[0m'
    PROMPT = '\033[95m'
########################

####Help section
if sys.argv[1].casefold() == 'Help' or sys.argv[1].casefold() == '--help' or sys.argv[1].casefold() == '-help':
    print(bcolors.PROMPT + '--- NO scan polar heatmap plotter ---')
    print(' -- Version 5.0 --')
    print('To run provide a .csv file with Y-axis = R and X-axis = Angle')
    print(bcolors.ENDC)
    quit()
#####

###Options####
options = {
    'save'              :   True,
    'units'             :   'cm',
    'rot_angle'         :   'NO',   #NO/NN/dih
    'show'              :   False,
    'up_lim'            :   1000,
    'up_lim_ticks'      :   20,
    'down_lim'          :   -240.0,
    'down_lim_ticks'    :   16,
    'R_ticks'           :   [3, 8],
    'R_limits'          :   [0, 8],
    'R_ticks_step'      : 1,
    'R_ticks_integers'  :   True,
    'inc. title'        :   False,
    'font_size'         :   15,
    'BSSE'              :   True,
    'partner'           :   'N2',
    'save_plotmat'      :   True,
    }
##############

class csv_file:
    """Handling commeon file reading and conversion to np.mat"""
    
    def __init__(self, file_name):
        self.name = file_name
        self.DF = pd.read_csv(file_name)
        self.mat = pd.DataFrame.to_numpy(self.DF)
        self.rvals = self.mat[:,0]
        self.angles = self.DF.columns.values[1:]
        self.EV = self.mat[:,1:]
        
def match_list(list_points,keys):
    matches = list_points.copy()
    for i in list_points:
        for key in keys:
            if key not in i.casefold() and i in matches: matches.remove(i)
    if len(matches) > 1: print('!!! error too many file matches with BSSE !!!'); quit()
    return matches[0]
        

def Counterpoise(file_info, options):
    #Find folders for relevant BSSE files
    cwd = os.getcwd().split('/')[:-1]
    parent_dir = "/".join(cwd)
    folders = os.listdir(parent_dir)
    bsse_no_folder = parent_dir + '/BSSE_NO'
    bsse_partner_folder = parent_dir + '/' + 'BSSE_' + options['partner']
    
    #Breakdown current file name to match with BSSE files:
    keys = file_info.name.casefold().split('.')[0].split('scan')[-1].split('_')
    if keys[0] == '': keys = keys[1:]
    
    #Match files to BSSE files and extract info
    #BSSE NO
    bsse_no_match = match_list(os.listdir(bsse_no_folder), keys)
    bsse_no_info = csv_file(bsse_no_folder + '/' + bsse_no_match)
    
    #BSSE partner
    bsse_partner_match = match_list(os.listdir(bsse_partner_folder), keys)
    bsse_partner_info = csv_file(bsse_partner_folder + '/' + bsse_partner_match)
    
    #Check to see everything matches
    if np.shape(file_info.EV) != np.shape(bsse_no_info.EV) and np.shape(file_info.EV) != np.shape(bsse_partner_info.EV):
        print("!!! error matrix sizes don't match!!"); quit()
    
    #subtract CP corrections from matrix
    correction_mat = bsse_no_info.EV + bsse_partner_info.EV
    file_info.EV = file_info.EV - (correction_mat)
    file_info.maxCP = np.amax(correction_mat)
    return file_info
        
####Run  
file_info = csv_file(sys.argv[1])
print('Working with file: ' + file_info.name)

if options['BSSE'] == True:
        print("Running with BSSE corrections")
        Counterpoise(file_info, options)

### Unit conversion
if options['units'] == 'cm':    #option for wavenumbers
    unit_mat = file_info.EV * 219474.6
elif options['units'] == 'ev':
    unit_mat = file_info.EV * 27.2114
    
#limits
mat_min = np.amin(unit_mat)
mat_max = np.amax(unit_mat)
zero_point = np.amin(unit_mat[0,:])
min_val = mat_min - zero_point
max_val = mat_max - zero_point
unit_mat =  unit_mat - zero_point    #Change to value of mat to diff

print('Positive limit is equal to: ' + str(max_val))
print('Negative limit is equal to: ' + str(min_val))    

###Title
NO_angle_name = '$\\Theta_{NO}$='
if options['partner'] == 'N2': partner_angle_name = '$\\Theta_{N_2}$='
else: partner_angle_name =  '$\\Theta_{D_2}$='
Dih_angle_name = '$\\psi$='
title_prefix = 'NO(A $\\Sigma^+$) + N$_2$(X $\\Sigma^+$) collision complex PES cut at fixed angles of \n'
title_suffix_1 = 'rotating over '
file_angle_1 = re.sub('\D', '', file_info.name.split('_')[0])
file_angle_2 = re.sub('\D', '', file_info.name.split('_')[1].split('.csv')[0])
if options['rot_angle'] == 'NO':
    title = title_prefix + partner_angle_name + str(file_angle_1) + ' ' + Dih_angle_name + str(file_angle_2) + ' ' \
        + title_suffix_1 + NO_angle_name + '0 to ' + NO_angle_name + '180'
elif options['rot_angle'] == 'NN':
    title = title_prefix + NO_angle_name + str(file_angle_1) + ' ' + Dih_angle_name + str(file_angle_2) + ' ' \
        + title_suffix_1 + partner_angle_name + '0 to ' + partner_angle_name + '180'
elif options['rot_angle'] == 'dih':
    title = title_prefix + NO_angle_name + str(file_angle_1) + ' ' + partner_angle_name + str(file_angle_2) + ' ' \
        + title_suffix_1 + Dih_angle_name + '0 to ' + Dih_angle_name + '180'
if options['inc. title'] == False:
    title = ""

###Plotting
    if options['units'] == 'cm': units = "$cm^{-1}$"
    angles = np.radians(np.linspace(0,180, len(file_info.angles)))
    len_r = len(file_info.angles)
    
    angle_mesh = np.zeros(np.shape(unit_mat))
    for i in range(np.shape(angle_mesh)[0]):
        angle_mesh[i,:] = np.array(angles).T   


    R_mesh = np.zeros(np.shape(unit_mat))
    for i in range(np.shape(R_mesh)[1]):
        R_mesh[:,i] = np.array(file_info.rvals)
 
    new_unit_mat = unit_mat.copy()
    #Extend unit_mat to within plot limits
    if R_mesh[-1, -1] != options['R_ticks'][0]:
        ext_R = options['R_ticks'][0]
        new_angle_row = angle_mesh[-1, :]
        new_R_row = np.full(len(angle_mesh[-1, :]), ext_R)
        new_unit_row = unit_mat[-1, :]
        
        angle_mesh = np.vstack([angle_mesh, new_angle_row])
        R_mesh = np.vstack([R_mesh, new_R_row])
        new_unit_mat = np.vstack([new_unit_mat, new_unit_row])
    
        
    my_gradient = LinearSegmentedColormap.from_list('my_gradient', (
    (0.000, (0.000, 0.114, 0.576)),
    (0.250, (0.000, 0.204, 0.996)),
    (0.475, (0.992, 0.992, 1.000)),
    (0.500, (1.000, 1.000, 1.000)),
    (0.525, (1.000, 1.000, 1.000)),
    (1.000, (1.000, 0.000, 0.047))))
    
    low_levels = np.linspace(options['down_lim'], 0, options['down_lim_ticks'])[0:-2]
    high_levels = np.linspace(0, options['up_lim'], options['up_lim_ticks'])
    levels = np.append(low_levels,high_levels)
    norm = TwoSlopeNorm(vmin = options['down_lim'], vcenter = 0, vmax = options['up_lim'])
    
    fig, axis = plt.subplots(subplot_kw=dict(projection='polar'))  
    matplotlib.rcParams.update({'font.size': options['font_size']})
    contourplot = axis.contourf(angle_mesh, R_mesh, new_unit_mat, cmap = my_gradient, norm=norm, levels=levels, extend='max')
    
    #labels
    axis.set_title(title)
    axis.set_xticks(np.arange(0, 2.0 * np.pi, np.pi / 6))
    axis.set_thetamin(0)
    axis.set_thetamax(180)
    axis.set_rmin(options['R_limits'][0])
    axis.set_rmax(options['R_limits'][1])
    r_ticks = np.arange(options['R_ticks'][0], (options['R_ticks'][1]+1), options['R_ticks_step'])
    r_ticks = np.round(r_ticks,1)
    r_ticks = np.insert(r_ticks, 0, 0)
    if options['R_ticks_integers'] == True: r_ticks = r_ticks.astype(int)
    axis.set_rticks(r_ticks)
    axis.tick_params(axis='y', which='major', labelsize=12)

    #side bar
    cbar = plt.colorbar(contourplot, shrink=.6, pad=0.08)
    cbar.set_label('$\\Delta$E (' + units + ')')

    label_position = axis.get_rlabel_position()
    axis.set_xlabel('R ($\\AA$)')
    axis.xaxis.set_label_coords(0.85, 0.17)
    axis.set_title(title, y=0.85)
        
    min_level = min(levels)
    max_level = max(levels)
    majors = [min_level, min_level/1.5, min_level/3, 0, max_level/3, max_level/1.5, max_level]
    majors = np.round(majors,0)
    labels = [str(i) for i in majors]
    cbar.ax.yaxis.set_major_locator(ticker.FixedLocator(majors))
    cbar.ax.set_yticklabels(labels)
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    #file handeling
    new_file = file_info.name.replace('.csv','.png')
    new_csv = file_info.name.replace('.csv', '_plotted.csv')

    if options['save'] == True: plt.savefig(new_file, dpi=1600)
    if options['show'] == True: plt.show()
    
    #Save matrix that is plotted
    plot_mat = np.zeros((np.shape(unit_mat)[0]+1, np.shape(unit_mat)[1]+1))
    angles = file_info.angles.copy().astype(float).astype(int)
    angles[-1] = 180
    plot_mat[0,1:] = angles
    plot_mat[1:,0] = file_info.rvals.copy().astype(float)
    plot_mat[1:,1:] = unit_mat.copy().astype(float)
    if options['save_plotmat'] == True: np.savetxt(new_csv, plot_mat, delimiter=', ')

now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
if 'plotting_log.csv' not in os.listdir(os.getcwd()):
    log = 'name, time/date, +ve limit, -ve limit, Saved?, save_name, BSSE?, units, save_plotmat? \n'
    log = log + file_info.name + ', ' + dt_string + ', ' + '{0:.5f}'.format(max_val) + \
    ', ' + '{0:.5f}'.format(min_val) + ', ' + "{0}".format(options['save']) + ', ' +  \
    new_file + ', ' + "{0}".format(options['BSSE']) + ', ' + options['units'] + ', ' +  "{0}".format(options['save_plotmat']) + '\n'
    
else: log = file_info.name + ', ' + dt_string + ', ' + '{0:.5f}'.format(max_val) + \
    ', ' + '{0:.5f}'.format(min_val) + ', ' + "{0}".format(options['save']) + ', ' + \
        new_file + ', ' + "{0}".format(options['BSSE']) + ', ' + options['units'] + ', ' +  "{0}".format(options['save_plotmat']) + '\n'
file = open('plotting_log.csv', 'a+').writelines(log)
