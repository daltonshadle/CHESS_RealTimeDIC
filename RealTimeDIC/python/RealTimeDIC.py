#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 14:57:43 2022

@author: djs522
"""

# CONSTANT_NAMES = all uppercase with underscores
# FunctionNames() = camel case with parentheses and without underscores
# variable_names = all lowercase with underscores

#*****************************************************************************
#%% IMPORTS
import numpy as np

from CPCorrFunctions import process_correlations
from RTDIC_Classes import dic_parameters, dic_paths, dic_matrices
from RTDIC_Widgets import dic_parameters_selector_widget, dic_continuous_update_widget, dic_field_value_widget

#*****************************************************************************
#%% USER INPUT

raw_dir = "/Users/jacksonearls/documents/GitHub/CHESS_RealTimeDIC/example/"
aux_dir = "/Users/jacksonearls/documents/GitHub/CHESS_RealTimeDIC/example/"
output_fname = 'dictest_output.txt'
img_template = 'dic_%06i.tiff'
sample_width = 1 # (mm) cross-sectional width of the sample for macro stress calculations
sample_thickness = 1 # (mm) cross-sectional thickness of the sample for macro stress calculations
sample_orientation_in_img = 'v' # describes sample orientaiton in image, 
                                # can be 'v' = vertical, 'h'=horizontal

#*****************************************************************************
#%% INITIALIZE OBJECTS
exp_dic_params = dic_parameters(sample_width=sample_width, 
                                sample_thickness=sample_thickness, 
                                sample_orientation=sample_orientation_in_img)
exp_dic_paths = dic_paths(base_dir=raw_dir, dic_json_dir=raw_dir, dic_par_dir=raw_dir, img_dir=raw_dir, 
                          img_fname_template=img_template, output_dir=aux_dir, 
                          output_fname=output_fname)
exp_dic_mats = dic_matrices()


#*****************************************************************************
#%% OPEN DIC PAR FILE AND CORRESPONDING JSON FILE
exp_dic_paths.open_dic_par_file(exp_dic_mats)

#*****************************************************************************
#%% SELECT FIRST IMAGE / REFERENCE IMAGE
exp_dic_paths.open_first_image()

#*****************************************************************************
#%% SET CONTROL POINT GRID AND DIC PARAMTERS
dpsw = dic_parameters_selector_widget(exp_dic_paths, exp_dic_params, exp_dic_mats, adjust_grid=True)
[exp_dic_paths, exp_dic_params, exp_dic_mats] = dpsw.get_all_dic_objects()

#*****************************************************************************
#%% PROCESS CURRENT IMAGES

first_img_dir = exp_dic_paths.get_first_img_dir()
img_nums = exp_dic_mats.dic_par_mat[:, 0]
num_imgs = img_nums.size
exp_dic_mats.reset_output_mat()

# for each file to process...
for i, cur_img_num in enumerate(img_nums):
    # process image index for full image and name
    cur = exp_dic_mats.get_dic_par_mat_at_img_num(cur_img_num)
    cur_img_num = cur[0]
    cur_force = cur[1]
    cur_screw = cur[2]
    cur_img_dir = exp_dic_paths.get_img_num_dir(cur_img_num)
    
    if i == 0:
        # process second image
        exp_dic_mats.cur_points = process_correlations(first_img_dir, cur_img_dir, 
                                                         exp_dic_mats.ref_points,
                                                         exp_dic_mats.ref_points,
                                                         exp_dic_params)
    else:
        # process the rest of the images
        exp_dic_mats.cur_points = process_correlations(first_img_dir, cur_img_dir, 
                                                         exp_dic_mats.ref_points,
                                                         exp_dic_mats.cur_points,
                                                         exp_dic_params)
    
    # do stress and strain calculations and write to output
    [cur_stress, cur_strain_xx, cur_strain_yy, cur_strain_xy] = exp_dic_mats.get_cur_stress_strain(cur_force, exp_dic_params)
    exp_dic_mats.add_to_output(np.array([cur_img_num, 
                                         cur_stress,
                                         cur_strain_xx,
                                         cur_strain_yy,
                                         cur_strain_xy,
                                         cur_screw]))
    
    # display to screen
    if exp_dic_params.is_sample_horizontal():
        print('Image: %i / %i \t DIC Image Number: %i \t Stress: %0.2f MPa \t Strain_XX: %0.3f %%' 
              %(i+1, num_imgs, cur_img_num, cur_stress, cur_strain_xx * 100))
    else:
        print('Image: %i / %i \t DIC Image Number: %i \t Stress: %0.2f MPa \t Strain_YY: %0.3f %%' 
              %(i+1, num_imgs, cur_img_num, cur_stress, cur_strain_yy * 100))

exp_dic_mats.save_output_mat_to_file(exp_dic_paths.get_output_full_dir())

#*****************************************************************************
#%% START REAL TIME PROCESSING
dcuw = dic_continuous_update_widget(exp_dic_paths, exp_dic_params, exp_dic_mats, prev_img_num=cur_img_num)
