#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 15:48:18 2020

@author: millerlab

Description: This script allows a user to simulate a DIC experiment with 
continuous updating from a previously recorded DIC experiment. This script is
mainly used for debugging purposes.

Usage: Run this script first and leave running in the background. Run 
RealTimeDIC.m selecting "new_dic.par" in the par file selection section. Once
RealTimeDIC.m reaches function call ContinuousUpdate.m, the user can then 
"simulate" new DIC measurements by hitting enter in this Python script.
RealTimeDIC.m will automatically update the stress-strain curve with the new
measurement.
"""

# *****************************************************************************
#%% IMPORTS
# *****************************************************************************
import os


# *****************************************************************************
#%% DECLARE PATH AND FILENAMES
# *****************************************************************************
# path to dic image folder
full_image_folder_path = 'C:/path/to/image/folder'

# filename of the complete par file stored in the full_image_folder_path
full_dic_par_name = 'dic.par'

# filename of the new/updating par file stored in the full_image_folder_path
new_dic_par_name = 'new_dic.par'


# *****************************************************************************
#%% MAIN FUNCTION
# *****************************************************************************
# open complete par file and read lines
full_par = open(os.path.join(full_image_folder_path, full_dic_par_name), 'r')
full_par_lines = full_par.readlines()

# check if the new par file exists, if so remove it
if os.path.isfile(os.path.join(full_image_folder_path, new_dic_par_name)):
    os.remove(os.path.join(full_image_folder_path, new_dic_par_name))

# set the number of images to be processed in first for loop of MATLAB DIC code
num_images_to_write_first = 5

# for each image/line in the complete par file, write it to the new par file
for i, line in enumerate(full_par_lines):
    # check the number of images to be proessed first before moving to sim exp
    if i == num_images_to_write_first:
        print('Start simulating experiment')
    if i > num_images_to_write_first:
        my_input = input('Press Enter for next image (q then enter to quit): ')
        if 'q' in my_input:
            break;
    
    # split line into info
    line_split = line.split(' ')
    image_num = line_split[5]
    load = line_split[6]
    screw_pos = line_split[7]
    
    # print all information
    print('Image #: %s \t Load: %s \t Screw Pos: %s' %(image_num, load, screw_pos))
    
    # open new par, write line, then close
    new_par = open(os.path.join(full_image_folder_path, new_dic_par_name), 'a')
    new_par.write(line)
    new_par.close()

print('Done!')