# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

#%% ***************************************************************************
# IMPORTS
import os
import numpy as np

import pickle

import pandas as pd

import tkinter as tk

# ***************************************************************************
# CLASS DECLARATION
class dic_matrices():
    def __init__(self, 
                 dic_par_mat=np.array([]), 
                 ref_points=np.array([]), 
                 cur_points=np.array([]),
                 output_mat=np.array([])):
        
        # initialize class variables
        self.dic_par_mat = dic_par_mat
        self.ref_points = ref_points
        self.cur_points = cur_points
        self.output_mat = output_mat
    
    # setters and getters
    def get_dic_par_mat(self):
        return self.dic_par_mat
    def set_dic_par_mat(self, dic_par_mat):
        self.dic_par_mat = dic_par_mat
    
    def get_ref_points(self):
        return self.ref_points
    def set_ref_points(self, ref_points):
        self.ref_points = ref_points
    
    def get_cur_points(self):
        return self.cur_points
    def set_cur_points(self, cur_points):
        self.cur_points = cur_points
    
    def get_output_mat(self):
        return self.output_mat
    def set_output_mat(self, output_mat):
        self.output_mat = output_mat
        
    # extra functions
    def process_dic_par_file(self, dic_par_dir):
        # day, month, day number, time, year, image number, load, screw position, extra
        df = pd.read_csv(dic_par_dir, sep=" ", header=None)
        self.dic_par_mat = np.array(df)[:, [5, 6, 7]]
    
    def load_dic_matrices_from_file(self, dic_matrices_dir):
        with open(dic_matrices_dir, "rb") as input_file:
             e = pickle.load(input_file)
             self.set_dic_par_mat(e.get_dic_par_mat())
             self.set_ref_points(e.get_ref_points())
             self.set_cur_points(e.get_cur_points())
             self.set_output_mat(e.get_output_mat())
             
    def save_dic_matrices_to_file(self, dic_matrices_dir):
        with open(dic_matrices_dir, "wb") as output_file:
            pickle.dump(self, output_file)
            
    def save_output_mat_to_file(self, output_dir):
        np.savetxt(output_dir, self.output_mat, fmt='%0.6f', delimiter='\t')
    
    def calc_grid_spacing(self):
        if len(self.ref_points.shape) == 2:
            xs = np.unique(self.ref_points[:, 0])
            ys = np.unique(self.ref_points[:, 1])
            
            return [np.abs(xs[0] - xs[1]), np.abs(ys[0] - ys[1])]
    
    def calc_bounding_box(self):
        if len(self.ref_points.shape) == 2:
            return np.array([[np.min(self.ref_points[:, 0]), np.min(self.ref_points[:, 1])],
                             [np.max(self.ref_points[:, 0]), np.max(self.ref_points[:, 1])]])
        else:
            return np.zeros([2,2])
    
    def fit_strain(self, ref_coords, curr_disp):
        # fit_strain - This function does a linear fit to calculate strain from displacements
        # 
        #   INPUT:
        # 
        #   ref_coords is a n x 2 matrix
        #      reference coodinates given as (x, y) pairs
        # 
        #   curr_disp is a n x 2 matrix
        #      current displacements to calculate strain from
        # 
        #   OUTPUT:
        # 
        #   p[0] is a float
        #      The linear coefficient of fitting displacements for strains
        #      
        #   Notes:
        #      Written and tested Spring 2017
        #
        #   REVISION:
        #      Carmen Fang, Chris Budrow (2017) - Originally written and tested
        #
        #      Tyrone Bracker, Dalton Shadle (2020) - Modified to include changing
        #      correlation and subregion size
        #
        #      Dalton Shadle (2022) - Converted to Python
        
        p = np.polyfit(ref_coords, curr_disp, 1);
        return p[0]
    
    def get_cur_stress_strain(self, cur_force, dic_params):
        cur_stress = cur_force / dic_params.get_sample_area()
        cur_strain_xx = self.fit_strain(self.ref_points[:, 0], (self.cur_points[:, 0]-self.ref_points[:, 0]))
        cur_strain_yy = self.fit_strain(self.ref_points[:, 1], (self.cur_points[:, 1]-self.ref_points[:, 1]))
        cur_strain_xy = (self.fit_strain(self.ref_points[:, 0], (self.cur_points[:, 1]-self.ref_points[:, 1])) 
                         + self.fit_strain(self.ref_points[:, 1], (self.cur_points[:, 0]-self.ref_points[:, 0]))) / 2
        return [cur_stress, cur_strain_xx, cur_strain_yy, cur_strain_xy]
    
    def add_to_output(self, cur_array):
        # img_num, stress, strain_xx, strain_yy, strain_xy, screw_pos
        cur_array = np.atleast_2d(cur_array)
        if cur_array.shape[1] != 6:
            cur_array = cur_array.T
        if self.output_mat.size == 0:
            self.output_mat = cur_array
        else:
            self.output_mat = np.vstack([self.output_mat, cur_array])
    
    def replace_in_output(self, cur_array, img_num):
        # img_num, stress, strain_xx, strain_yy, strain_xy, screw_pos
        cur_array = np.atleast_2d(cur_array)
        if cur_array.shape[1] != 6:
            cur_array = cur_array.T
        try:
            [o, img_idx] = self.get_output_mat_at_img_num(img_num, return_idx=True)
            self.output_mat[img_idx, :] = cur_array
        except Exception as e:
            print(e)
    
    def reset_output_mat(self):
        self.output_mat = np.array([])
    
    def reset_cur_points(self):
        self.set_cur_points(self.get_ref_points())
        
    def get_dic_par_mat_at_img_num(self, img_num, return_idx=False):
        img_idx = np.where(self.dic_par_mat[:, 0] == img_num)[0]
        if return_idx:
            return self.dic_par_mat[img_idx, :].flatten(), img_idx
        else:
            return self.dic_par_mat[img_idx, :].flatten()
        
    def get_output_mat_at_img_num(self, img_num, return_idx=False):
        img_idx = np.where(self.output_mat[:, 0] == img_num)[0]
        if return_idx:
            return self.output_mat[img_idx, :].flatten(), img_idx
        else:
            return self.output_mat[img_idx, :].flatten()
    
    # str and rep
    def __repr__(self):
        return "dic_matrices()\n" + self.__str__()
    def __str__(self):
        class_dict = {'dic_par_mat' : self.dic_par_mat,
        'ref_points' : self.ref_points,
        'cur_points' : self.cur_points,
        'output_mat' : self.output_mat}
        
        return str(class_dict)
        

class dic_parameters():
    def __init__(self, 
                 grid_spacing=[40, 40], 
                 fixed_corr_dimen=[100, 100], 
                 moving_corr_dimen=[40, 40],
                 sample_width=1,
                 sample_thickness=1):
        
        # initialize class variables
        self.grid_spacing = grid_spacing
        self.fixed_corr_dimen = fixed_corr_dimen
        self.moving_corr_dimen = moving_corr_dimen
        self.sample_width = sample_width
        self.sample_thickness = sample_thickness
    
    # getters and setters
    def get_grid_spacing(self):
        return self.grid_spacing
    def set_grid_spacing(self, grid_spacing):
        self.grid_spacing = grid_spacing
    
    def get_fixed_corr_dimen(self):
        return self.fixed_corr_dimen
    def set_fixed_corr_dimen(self, fixed_corr_dimen):
        self.fixed_corr_dimen = fixed_corr_dimen
    
    def get_moving_corr_dimen(self):
        return self.moving_corr_dimen
    def set_moving_corr_dimen(self, moving_corr_dimen):
        self.moving_corr_dimen = moving_corr_dimen
        
    def get_sample_width(self):
        return self.sample_width
    def set_sample_width(self, sample_width):
        self.sample_width = sample_width
    
    def get_sample_thickness(self):
        return self.sample_thickness
    def set_sample_thickness(self, sample_thickness):
        self.sample_thickness = sample_thickness
    
    def get_sample_area(self):
        return self.sample_width * self.sample_thickness
    
    # extra fuctions
    def load_dic_parameters_from_file(self, dic_params_dir):
        with open(dic_params_dir, "rb") as input_file:
             e = pickle.load(input_file)
             self.set_grid_spacing(e.get_grid_spacing())
             self.set_fixed_corr_dimen(e.get_fixed_corr_dimen())
             self.set_moving_corr_dimen(e.get_moving_corr_dimen())
             self.set_sample_width(e.get_sample_width())
             self.set_sample_thickness(e.get_sample_thickness())
             
    def save_dic_parameters_to_file(self, dic_params_dir):
        with open(dic_params_dir, "wb") as output_file:
            pickle.dump(self, output_file)
    
    # str and rep
    def __repr__(self):
        return "dic_parameters()\n" + self.__str__()
    def __str__(self):
        class_dict = {'grid_spacing' : self.grid_spacing,
        'fixed_corr_dimen' : self.fixed_corr_dimen,
        'moving_corr_dimen' : self.moving_corr_dimen,
        'sample_width' : self.sample_width,
        'sample_thickness' : self.sample_thickness}
        
        return str(class_dict)

class dic_paths():
    def __init__(self, 
                 base_dir=os.getcwd()):
        
        # intialize class variables
        self.base_dir = base_dir
        self.dic_par_dir = base_dir
        self.dic_par_fname = 'dic.par'
        self.img_dir = base_dir
        self.img_fname_template = 'dic_%06i.tiff'
        self.output_dir = base_dir
        self.output_fname = 'output.txt'
        self.first_img_num = 0
        
    # setters and getters
    def get_base_dir(self):
        return self.base_dir
    def set_base_dir(self, base_dir):
        self.base_dir = base_dir
    
    def get_dic_par_dir(self):
        return self.dic_par_dir
    def set_dic_par_dir(self, dic_par_dir):
        self.dic_par_dir = dic_par_dir
    
    def get_dic_par_fname(self):
        return self.dic_par_fname
    def set_dic_par_fname(self, dic_par_fname):
        self.dic_par_fname = dic_par_fname
    
    def get_img_dir(self):
        return self.img_dir
    def set_img_dir(self, img_dir):
        self.img_dir = img_dir
    
    def get_img_fname_template(self):
        return self.img_fname_template
    def set_img_fname_template(self, img_fname_template):
        self.img_fname_template = img_fname_template
    
    def get_output_dir(self):
        return self.output_dir
    def set_output_dir(self, output_dir):
        self.output_dir = output_dir
    
    def get_output_fname(self):
        return self.output_fname
    def set_output_fname(self, output_fname):
        self.output_fname = output_fname
    
    def get_first_img_num(self):
        return self.first_img_num
    def set_first_img_num(self, first_img_num):
        self.first_img_num = first_img_num 
    
    # extra functions
    def open_dic_par_file(self, dic_mats):
        root = tk.Tk()
        root.withdraw()
        
        dic_par_dir = tk.filedialog.askopenfilename(initialdir=self.base_dir, 
                                                    defaultextension='.par',
                                                    filetypes=[("dic par files", "*.par")],
                                                    title='Select dic.par File')
        
        if not dic_par_dir:
            quit()
        else:
            try:
                self.set_dic_par_dir(os.path.dirname(dic_par_dir))
                self.set_dic_par_fname(os.path.basename(dic_par_dir))
                dic_mats.process_dic_par_file(dic_par_dir)
            except Exception as e:
                print(e)
                self.open_dic_par_file(dic_mats)
    
    def open_first_image(self):
        root = tk.Tk()
        root.withdraw()
        
        first_img_dir = tk.filedialog.askopenfilename(initialdir=self.base_dir, 
                                                    defaultextension='.tiff',
                                                    filetypes=[("TIFF Files", "*.tiff")],
                                                    title='Select First DIC Image File')
        
        if not first_img_dir:
            quit()
        else:
            try:
                self.set_img_dir(os.path.dirname(first_img_dir))
                first_img_fname = os.path.basename(first_img_dir)
                self.set_first_img_num(int((first_img_fname.split('_')[-1]).split('.')[0]))
            except:
                self.open_first_image()
    
    def get_img_num_dir(self, img_num):
        return os.path.join(self.img_dir, self.img_fname_template %(img_num))
    def get_first_img_dir(self):
        return self.get_img_num_dir(self.get_first_img_num())
    
    def get_dic_par_full_dir(self):
        return os.path.join(self.dic_par_dir, self.dic_par_fname)
    def get_output_full_dir(self):
        return os.path.join(self.output_dir, self.output_fname)
    
    # str and rep
    def __repr__(self):
        return "dic_paths()\n" + self.__str__()
    def __str__(self):
        class_dict = {'base_dir' : self.base_dir,
        'dic_par_dir' : self.dic_par_dir,
        'dic_par_fname' : self.dic_par_fname,
        'img_dir' : self.img_dir,
        'img_fname_template' : self.img_fname_template,
        'output_dir' : self.output_dir,
        'output_fname' : self.output_fname,
        'first_img_num' : self.first_img_num}
        
        return str(class_dict)