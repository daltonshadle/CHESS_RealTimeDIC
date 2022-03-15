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
from RTDIC_Widgets import dic_parameters_selector_widget, dic_continuous_update_widget

#*****************************************************************************
#%% USER INPUT

base_dir = "/home/djs522/additional_sw/RealTimeDIC/CHESS_RealTimeDIC/example/dic_images/"
img_template = 'dic_%06i.tiff'
sample_width = 1 # (mm) cross-sectional width of the sample for macro stress calculations
sample_thickness = 1 # (mm) cross-sectional thickness of the sample for macro stress calculations

#*****************************************************************************
#%% INITIALIZE OBJECTS
exp_dic_params = dic_parameters()
exp_dic_paths = dic_paths(base_dir=base_dir)
exp_dic_mats = dic_matrices()

exp_dic_paths.set_img_fname_template(img_template)
exp_dic_params.set_sample_width(sample_width)
exp_dic_params.set_sample_thickness(sample_thickness)

#*****************************************************************************
#%% OPEN DIC PAR FILE
exp_dic_paths.open_dic_par_file(exp_dic_mats)

#*****************************************************************************
#%% SELECT FIRST IMAGE / REFERENCE IMAGE
exp_dic_paths.open_first_image(exp_dic_mats)

#*****************************************************************************
#%% SET CONTROL POINT GRID AND DIC PARAMTERS
dpsw = dic_parameters_selector_widget(exp_dic_paths, exp_dic_params, exp_dic_mats, adjust_grid=True)
[exp_dic_paths, exp_dic_params, exp_dic_mats] = dpsw.get_all_dic_objects()

#*****************************************************************************
#%% PROCESS CURRENT IMAGES
first_img_idx = exp_dic_paths.get_first_img_idx()
last_img_idx = exp_dic_paths.get_last_img_idx()
img_indices = np.arange(first_img_idx, last_img_idx, step=1)
num_imgs = img_indices.size
first_img_dir = exp_dic_paths.get_first_img_dir()

exp_dic_mats.reset_output_mat()


# for each file to process...
for i in range(num_imgs):
    # process image index for full image and name
    [cur_img_num, cur_force, cur_screw] = exp_dic_mats.get_dic_par_mat()[img_indices[i], :]
    cur_img_dir = exp_dic_paths.get_img_num_dir(cur_img_num)
    
    if i == 0:
        # process second image
        exp_dic_mats.set_cur_points(process_correlations(first_img_dir, cur_img_dir, 
                                                         exp_dic_mats.get_ref_points(),
                                                         exp_dic_mats.get_ref_points(),
                                                         exp_dic_params))
    else:
        # process the rest of the images
        exp_dic_mats.set_cur_points(process_correlations(first_img_dir, cur_img_dir, 
                                                         exp_dic_mats.get_ref_points(),
                                                         exp_dic_mats.get_cur_points(),
                                                         exp_dic_params))
    
    # do stress and strain calculations and write to output
    [cur_stress, cur_strain_xx, cur_strain_yy, cur_strain_xy] = exp_dic_mats.get_cur_stress_strain(cur_force, exp_dic_params)
    exp_dic_mats.add_to_output(np.array([cur_img_num, 
                                         cur_stress,
                                         cur_strain_xx,
                                         cur_strain_yy,
                                         cur_strain_xy,
                                         cur_screw]))
    
    # display to screen
    print('Image: %i / %i \t DIC Image Number: %i \t Stress: %0.2f MPa \t Strain_XX: %0.3f %%' 
          %(i+1, num_imgs, cur_img_num, cur_stress, cur_strain_xx * 100))

exp_dic_mats.save_output_mat_to_file(exp_dic_paths.get_output_full_dir())

#*****************************************************************************
#%% START REAL TIME PROCESSING

dcuw = dic_continuous_update_widget(exp_dic_paths, exp_dic_params, exp_dic_mats, prev_img_idx=last_img_idx)

