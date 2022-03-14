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
from RTDIC_Widgets import dic_paramteres_selector_widget

#*****************************************************************************
#%% USER INPUT

base_dir = "/home/djs522/additional_sw/RealTimeDIC/CHESS_RealTimeDIC/example/dic_images/"


#*****************************************************************************
#%% INITIALIZE OBJECTS
exp_dic_params = dic_parameters()
exp_dic_paths = dic_paths(base_dir=base_dir)
exp_dic_mats = dic_matrices()

#*****************************************************************************
#%% OPEN DIC PAR FILE
exp_dic_paths.open_dic_par_file(exp_dic_mats)

#*****************************************************************************
#%% SELECT FIRST IMAGE / REFERENCE IMAGE
exp_dic_paths.open_first_image()

#*****************************************************************************
#%% SET CONTROL POINT GRID AND DIC PARAMTERS
dpsw = dic_paramteres_selector_widget(exp_dic_paths, exp_dic_params, exp_dic_mats, adjust_grid=True)
[exp_dic_paths, exp_dic_params, exp_dic_mats] = dpsw.get_all_dic_objects()

#*****************************************************************************
#%% PROCESS CURRENT IMAGES
first_img_idx = exp_dic_mats.get_img_idx_from_num(exp_dic_paths.get_first_img_num())
last_img_idx = exp_dic_mats.get_dic_par_mat().shape[0]
img_indices = np.arange(first_img_idx, last_img_idx, step=1)
num_imgs = img_indices.size
first_img_dir = exp_dic_paths.get_img_num_dir(exp_dic_paths.get_first_img_num())

# img_num, stress, strain_xx, strain_yy, strain_xy, screw_pos
output = np.zeros([num_imgs, 6])

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
    
    
    '''
    output[i, :] = np.array([cur_img_num, 
                             cur_stress,
                             cur_strain_xx,
                             cur_strain_yy,
                             cur_strain_xy,
                             cur_screw])
    np.savetxt(os.path.join(base_dir, output_filename), output, fmt='%0.6f', delimiter='\t')
    '''
    
    # display to screen
    print('Image: %i / %i \t DIC Image Number: %i \t Stress: %0.2f MPa \t Strain_XX: %0.3f %%' 
          %(i+1, num_imgs, cur_img_num, cur_stress, cur_strain_xx * 100))



'''

#%%

# paths and filenames
base_dir = "/home/djs522/additional_sw/RealTimeDIC/CHESS_RealTimeDIC/example/"

# dic parameters
grid_spacing = 60 # spacing of correlation points
fixed_corr_dimen = [150, 80] # (pixels) fixed correlation area dimensions for control point [width, height]
moving_corr_dimen = [40, 40] # (pixels) moving subregion area dimensions for control point [width, height]

# image processing options
img_processing_gap = 1

# sample geometry
sample_width = 1 # (mm) cross-sectional width of the sample for macro stress calculations
sample_thickness = 1 # (mm) cross-sectional thickness of the sample for macro stress calculations
dic_matssample_area = sample_width * sample_thickness # (mm^2) calculate sample cross-sectional area for stress calculations

# bool parameter options
save_stress_strain = True
reference_process = False
double_process = False
debug = False


#%% preprocessing
dic_params = dic_parameters(grid_spacing=grid_spacing, 
                            fixed_corr_dimen=fixed_corr_dimen, 
                            moving_corr_dimen=moving_corr_dimen)
output_dir = os.path.join(base_dir, output_filename)

#%% select the dic.par file for current sample, contains image index, force, and screw position
[dic_par_mat, dic_par_dir] = process_dic_par_file(base_dir=None)

#%% select the reference image
[img_dir, first_img_fname, first_img_num, first_img_idx] = select_first_image(dic_par_mat, base_dir=None)
img_dir_template = os.path.join(img_dir, img_template)
first_img_dir = img_dir_template %(first_img_num)

#%% generate grid from the first reference image
t = grid_selector_widget(first_img_dir, save_dir=base_dir, dic_params=dic_params)    
dic_params = t.get_dic_params()
ref_coords = t.get_grid_coords()

#%%
if debug:
    ref_coords = np.load(os.path.join(base_dir, 'grid_coords.npy')).astype(float)
    print(ref_coords)

#%% calculate stress and strain for images that are already in dic.par first
last_img_idx = dic_par_mat.shape[0]
img_indices = np.arange(first_img_idx, last_img_idx, img_processing_gap)
num_imgs = img_indices.size

# img_num, stress, strain_xx, strain_yy, strain_xy, screw_pos
output = np.zeros([num_imgs, 6])

# for each file to process...
for i in range(num_imgs):#range(num_imgs):
    # process image index for full image and name
    curr_img_num, curr_force, curr_screw = dic_par_mat[img_indices[i], :]
    curr_img_dir = img_dir_template %(curr_img_num)
    
    # special case for second image
    if i == 0:
        # process second image
        curr_coords = process_correlations(first_img_dir, curr_img_dir, ref_coords, 
                         ref_coords, dic_params=dic_params)
    else:
        # process the rest of the images, double processing, if necessary
        
        if reference_process:
            curr_coords = process_correlations(first_img_dir, curr_img_dir, ref_coords, 
                         ref_coords, dic_params=dic_params)
        if double_process:
            curr_coords = process_correlations(first_img_dir, curr_img_dir, ref_coords, 
                         curr_coords, dic_params=dic_params)

        curr_coords = process_correlations(first_img_dir, curr_img_dir, ref_coords, 
                         curr_coords, dic_params=dic_params)
    
    # do stress and strain calculations and write to output
    curr_stress = curr_force / sample_area
    curr_strain_xx = fit_strain(ref_coords[:, 0], (curr_coords[:, 0]-ref_coords[:, 0]))
    curr_strain_yy = fit_strain(ref_coords[:, 1], (curr_coords[:, 1]-ref_coords[:, 1]))
    curr_strain_xy = (fit_strain(ref_coords[:, 0], (curr_coords[:, 1]-ref_coords[:, 1])) + fit_strain(ref_coords[:, 1], (curr_coords[:, 0]-ref_coords[:, 0]))) / 2
    
    output[i, :] = np.array([curr_img_num, 
                             curr_stress,
                             curr_strain_xx,
                             curr_strain_yy,
                             curr_strain_xy,
                             curr_screw])
    np.savetxt(os.path.join(base_dir, output_filename), output, fmt='%0.6f', delimiter='\t')
    
    # display to screen
    print('Image: %i / %i \t DIC Image Number: %i \t Stress: %0.2f MPa \t Strain_XX: %0.3f %%' 
          %(i+1, num_imgs, curr_img_num, output[i, 1], output[i, 2] * 100))

if debug:
    fig = plt.figure()
    plt.scatter(ref_coords[:, 0], ref_coords[:, 1])
    plt.scatter(curr_coords[:, 0], curr_coords[:, 1])

    ## intial plotting
    stress_strain_figure = plt.figure();
    plt.scatter(output[:, 2] * 100, output[:, 1], c='b');
    plt.xlabel('Macroscopic Strain (%)')
    plt.ylabel('Macroscopic Stress (MPa)')
    plt.show()

#%% continuous updating as new dic images are collected


window = tk.Tk()
start = continuous_update_widget(window, dic_par_dir, dic_params, sample_area, 
                         img_dir_template, first_img_idx, img_indices[i], 
                         ref_coords, curr_coords, output, output_dir)
def on_closing():
    window.destroy()
window.protocol("WM_DELETE_WINDOW", on_closing)
window.mainloop()




'''










