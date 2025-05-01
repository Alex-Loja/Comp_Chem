#!/usr/bin/env python3
import sys
import math
import os
from math import floor

######Molpro_Time#####
#Ver 1.0
#Simple script to calculate the average calculation time of a subset of files
#run in Molpro


class bcolors :
    Yes = '\033[96m'
    No = '\033[0m'
###Help 
if sys.argv[1].casefold() == 'help' or sys.argv[1].casefold() == '--help':
    print(bcolors.Yes + '-- Molpro time Script Version 0.0 --')
    print('Collect and average run times for all .out files')
    print('To run please provide common prefix of out files')
    print('Ensure all .out files are present in current folder')
    print('Tested to molpro v2024.1')
    print(bcolors.No)
    quit()
####

#Select files
input = sys.argv[1]
cwd = os.getcwd()
files = []
for i in os.listdir(cwd):
    if '.out_' in i:
        continue
    elif '.out' in i:
        if input in i:
            files.append(i)
            
def readfile(file):
    cont = open(file,'r').readlines()
    times = []
    for i in cont:
        if 'REAL TIME ' in i:
            times.append(i.strip().split()[-2])
    time = float(times[-1])
    return time
         
all_times = []
final_time = 0
for i in files:
    file_time = readfile(i)
    all_times.append(file_time)
    final_time += file_time
final_time = final_time/len(all_times)


second = 0
minute = 0
hour = 0
if final_time >= 60:
    minute = floor(final_time/60)
    second = floor(((final_time/60) - minute) * 60)
    if minute >= 60:
        hour = floor(minute/60)
        rem_min = floor(((minute/60) - hour) * 60)
        minute = rem_min
    
    
result = str(hour) + ':' + str(minute) + ':' + str(second)     

print(result)
