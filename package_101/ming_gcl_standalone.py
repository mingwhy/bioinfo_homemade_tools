"""
https://github.com/HamdiTools/GCL_Tool_BIU
GCL - Global Coordination Level Python Tool.
November 2021.
Code written by Omer Hamdi - omerhamdilf2@gmail.com.
August 2022. Editted by Ming
"""
"""
Requirements:
- Python 3 and above (3.9.1 is recommended).
- Packages for python: (if you don't have one of the following, right click inside python on the missing import -> Show Context Actions -> Install Package):
numpy, math, random, scipy, threading. for the example you will also need: matplotlib, sys, pathlib, os.

How to include the tool in your own python project:
Put the gcl_library.py file inside your python project.
Add the import: import gcl_library as gcl_lib
"""

# Imports:
import gcl_library as gcl_lib

# check if input file/s exists:
import sys
from pathlib import Path
from os import path
# for threading:
import threading
# math tools:
import numpy as np
from scipy.spatial.distance import cdist
import math
import random
# gcl library import:
from matplotlib import pyplot as plt

files_arr = 'old.csv,young.csv'.split(',')
files_arr = ['old.csv', 'young.csv']
jack_knifes=70
num_divisions=10
jack_knife_percentage=0.8 
task='jackknife'

"""
The main gcl start function to generate and present the histograms.
param files_arr: Array of .CSV files.
param jack_knifes: Number of jackknives for calculation.
param num_divisions: Number of random gene division for calculation.
param jack_knife_percentage: Percentage of cells to choose for jackknife realization.
param task: String, either 'jackknife' or 'regular_calc' for the requested task.
return: none.
"""

result_arr, file_names = [], []
for file in files_arr:
    csv_mat = np.genfromtxt(file, delimiter=',')
    file_names.append(Path(file).stem)
    if task == 'jackknife':
        result_arr.append(gcl_lib.jackknife(csv_mat, jack_knifes, jack_knife_percentage, num_divisions))
    elif task == 'regular_calc':
        result_arr.append(gcl_lib.gcl(csv_mat, num_divisions))
        print('GCL value of ' + file_names[-1] + ' is: ' + str(result_arr[-1]))
    else:
        raise NameError('Invalid Task!')

# plotting the histogram/s:
plt.figure(figsize=(10, 10))
if task == 'jackknife':
    for result in range(len(result_arr)):
        plt.hist(result_arr[result], 8, density=False, edgecolor='black', label=file_names[result], alpha=.8)
    plot_title = 'GCL ' + (
        'jackknife ' if task == 'jackknife' else 'Regular Calculation ') + 'Hist. with: ' + str(
        jack_knifes) + ' realizations' + ', ' + str(jack_knife_percentage * 100) + '%'
    plt.title(plot_title, fontsize=26)
# plot settings:
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('GCL', fontsize=32)
plt.ylabel('Iterations', fontsize=32)
plt.rc('legend', fontsize=26)
plt.legend(loc='upper right')
plt.show()


