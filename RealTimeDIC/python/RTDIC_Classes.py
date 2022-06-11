# -*- coding: utf-8 -*-

# -*- coding: utf-8 -*-

#%% ***************************************************************************
# IMPORTS
import os

import numpy as np

import pickle

import json

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
        self._dic_par_mat = dic_par_mat
        self._ref_points = ref_points
        self._cur_points = cur_points
        self._output_mat = output_mat
    
    # setters and getters
    @property
    def dic_par_mat(self):
        return self._dic_par_mat
    @dic_par_mat.setter
    def dic_par_mat(self, dic_par_mat):
        self._dic_par_mat = dic_par_mat
    
    @property
    def ref_points(self):
        return self._ref_points
    @ref_points.setter
    def ref_points(self, ref_points):
        self._ref_points = ref_points
    
    @property
    def cur_points(self):
        return self._cur_points
    @cur_points.setter
    def cur_points(self, cur_points):
        self._cur_points = cur_points
    
    @property
    def output_mat(self):
        return self._output_mat
    @output_mat.setter
    def output_mat(self, output_mat):
        self._output_mat = output_mat
        
    # extra functions
    def process_dic_par_file(self, dic_par_dir, column_nums):
        # day, month, day number, time, year, image number, load, screw position, extra
        # {"0": "date", "1": "time", "2": "epoch", "3": "filenumber", "4": "TEN", "5": "Load1", "6": "SAM_BOT_ROT"}
        # can eventually point to a .json file that has the data labeled in dic.par
        df = pd.read_csv(dic_par_dir, sep=" ", header=None)
        # image number, load, screw position
        self._dic_par_mat = np.array(df)[:, column_nums]
        
        # old dic.par format for image number, load, screw position
        #self._dic_par_mat = np.array(df)[:, [5, 6, 7]]
    
    def load_dic_matrices_from_file(self, dic_matrices_dir):
        with open(dic_matrices_dir, "rb") as input_file:
             e = pickle.load(input_file)
             self._dic_par_mat = e.dic_par_mat
             self._ref_points = e.ref_points
             self._cur_points = e.cur_points
             self._output_mat = e.output_mat
             
    def save_dic_matrices_to_file(self, dic_matrices_dir):
        with open(dic_matrices_dir, "wb") as output_file:
            pickle.dump(self, output_file)
            
    def save_output_mat_to_file(self, output_dir):
        np.savetxt(output_dir, self._output_mat, fmt='%0.6f', delimiter='\t')
    
    def calc_grid_spacing(self):
        if len(self._ref_points.shape) == 2:
            xs = np.unique(self._ref_points[:, 0])
            ys = np.unique(self._ref_points[:, 1])
            
            return [np.abs(xs[0] - xs[1]), np.abs(ys[0] - ys[1])]
    
    def calc_bounding_box(self):
        if len(self._ref_points.shape) == 2:
            return np.array([[np.min(self._ref_points[:, 0]), np.min(self._ref_points[:, 1])],
                             [np.max(self._ref_points[:, 0]), np.max(self._ref_points[:, 1])]])
        else:
            return np.zeros([2, 2])
    
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
        cur_strain_xx = self.fit_strain(self._ref_points[:, 0], (self._cur_points[:, 0]-self._ref_points[:, 0]))
        cur_strain_yy = self.fit_strain(self._ref_points[:, 1], (self._cur_points[:, 1]-self._ref_points[:, 1]))
        cur_strain_xy = (self.fit_strain(self._ref_points[:, 0], (self._cur_points[:, 1]-self._ref_points[:, 1])) 
                         + self.fit_strain(self._ref_points[:, 1], (self._cur_points[:, 0]-self._ref_points[:, 0]))) / 2
        return [cur_stress, cur_strain_xx, cur_strain_yy, cur_strain_xy]
    
    def add_to_output(self, cur_array):
        # img_num, stress, strain_xx, strain_yy, strain_xy, screw_pos
        cur_array = np.atleast_2d(cur_array)
        if cur_array.shape[1] != 6:
            cur_array = cur_array.T
        if self._output_mat.size == 0:
            self._output_mat = cur_array
        else:
            self._output_mat = np.vstack([self._output_mat, cur_array])
    
    def replace_in_output(self, cur_array, img_num):
        # img_num, stress, strain_xx, strain_yy, strain_xy, screw_pos
        cur_array = np.atleast_2d(cur_array)
        if cur_array.shape[1] != 6:
            cur_array = cur_array.T
        try:
            [o, img_idx] = self.get_output_mat_at_img_num(img_num, return_idx=True)
            self._output_mat[img_idx, :] = cur_array
        except Exception as e:
            print(e)
    
    def reset_output_mat(self):
        self._output_mat = np.array([])
    
    def reset_cur_points(self):
        self._cur_points = self._ref_points
        
    def get_dic_par_mat_at_img_num(self, img_num, return_idx=False):
        img_idx = np.where(self._dic_par_mat[:, 0] == img_num)[0]
        if return_idx:
            return self._dic_par_mat[img_idx, :].flatten(), img_idx
        else:
            return self._dic_par_mat[img_idx, :].flatten()
        
    def get_output_mat_at_img_num(self, img_num, return_idx=False):
        img_idx = np.where(self._output_mat[:, 0] == img_num)[0]
        if return_idx:
            return self._output_mat[img_idx, :].flatten(), img_idx
        else:
            return self._output_mat[img_idx, :].flatten()
    
    # str and rep
    def __repr__(self):
        return "dic_matrices()\n" + self.__str__()
    def __str__(self):
        class_dict = {'dic_par_mat' : self._dic_par_mat,
        'ref_points' : self._ref_points,
        'cur_points' : self._cur_points,
        'output_mat' : self._output_mat}
        
        return str(class_dict)
        

class dic_parameters():
    def __init__(self, 
                 grid_spacing=[40, 40], 
                 fixed_corr_dimen=[100, 100], 
                 moving_corr_dimen=[40, 40],
                 sample_width=1,
                 sample_thickness=1,
                 sample_orientation='horizontal'):
        
        # declare constants
        self._vert_const = ['v', 'vertical']
        self._horz_const = ['h', 'horizontal']
        
        # initialize class variables
        self._grid_spacing = grid_spacing
        self._fixed_corr_dimen = fixed_corr_dimen
        self._moving_corr_dimen = moving_corr_dimen
        self._sample_width = sample_width
        self._sample_thickness = sample_thickness
        self._sample_orientation = sample_orientation
    
    # getters and setters
    @property
    def grid_spacing(self):
        return self._grid_spacing
    @grid_spacing.setter
    def grid_spacing(self, grid_spacing):
        self._grid_spacing = grid_spacing
    
    @property
    def fixed_corr_dimen(self):
        return self._fixed_corr_dimen
    @fixed_corr_dimen.setter
    def fixed_corr_dimen(self, fixed_corr_dimen):
        self._fixed_corr_dimen = fixed_corr_dimen
    
    @property
    def moving_corr_dimen(self):
        return self._moving_corr_dimen
    @moving_corr_dimen.setter
    def moving_corr_dimen(self, moving_corr_dimen):
        self._moving_corr_dimen = moving_corr_dimen
    
    @property
    def sample_width(self):
        return self._sample_width
    @sample_width.setter
    def sample_width(self, sample_width):
        self._sample_width = sample_width
    
    @property
    def sample_thickness(self):
        return self._sample_thickness
    @sample_thickness.setter
    def sample_thickness(self, sample_thickness):
        self._sample_thickness = sample_thickness
    
    @property
    def sample_orientation(self):
        return self._sample_orientation
    @sample_orientation.setter
    def sample_orientation(self, sample_orientation):
        sample_orientation = str(sample_orientation).lower()
        if sample_orientation in self._vert_const or sample_orientation in self._horz_const:
            self._sample_orientation = sample_orientation
        else:
            raise ValueError("Sample orientation '%s' needs to be set as:\n\
                  vertical: %s\n horizontal: %s" %(sample_orientation, self._vert_const, self._horz_const))
    
    def is_sample_horizontal(self):
        return self._sample_orientation in self._horz_const
    
    def get_sample_area(self):
        return self._sample_width * self._sample_thickness
    
    # extra fuctions
    def load_dic_parameters_from_file(self, dic_params_dir):
        with open(dic_params_dir, "rb") as input_file:
             e = pickle.load(input_file)
             self._grid_spacing = e.grid_spacing
             self._fixed_corr_dimen = e.fixed_corr_dimen
             self._moving_corr_dimen = e.moving_corr_dimen
             self._sample_width = e.sample_width
             self._sample_thickness = e.sample_thickness
             
    def save_dic_parameters_to_file(self, dic_params_dir):
        with open(dic_params_dir, "wb") as output_file:
            pickle.dump(self, output_file)
    
    # str and rep
    def __repr__(self):
        return "dic_parameters()\n" + self.__str__()
    def __str__(self):
        class_dict = {'grid_spacing' : self._grid_spacing,
        'fixed_corr_dimen' : self._fixed_corr_dimen,
        'moving_corr_dimen' : self._moving_corr_dimen,
        'sample_width' : self._sample_width,
        'sample_thickness' : self._sample_thickness}
        
        return str(class_dict)

class dic_paths():
    def __init__(self, 
                 base_dir=os.getcwd(),
                 dic_json_dir=os.getcwd(),
                 dic_json_fname='dic.json',
                 json_col_nums=[0]*3,
                 dic_par_dir=os.getcwd(),
                 dic_par_fname='dic.par',
                 img_dir=os.getcwd(),
                 img_fname_template='dic_%06i.tif',
                 output_dir=os.getcwd(),
                 output_fname='output.txt'
                 ):
        
        # intialize class variables
        self._base_dir = base_dir
        self._dic_json_dir = dic_json_dir
        self._dic_json_fname = dic_json_fname
        self._json_col_nums = json_col_nums
        self._dic_par_dir = dic_par_dir
        self._dic_par_fname = dic_par_fname
        self._img_dir = img_dir
        self._img_fname_template = img_fname_template
        self._output_dir = output_dir
        self._output_fname = output_fname
        self._first_img_num = 0
        
    # setters and getters  
    @property
    def base_dir(self):
        return self._base_dir
    @base_dir.setter
    def base_dir(self, base_dir):
        if os.path.exists(base_dir):
            self._base_dir = base_dir
        else:
            raise ValueError("Base directory '%s' does not exists" %(base_dir))
    
    @property
    def dic_json_dir(self):
        return self._dic_json_dir
    @dic_json_dir.setter
    def dic_json_dir(self, dic_json_dir):
        if os.path.exists(dic_json_dir):
            self._dic_json_dir = dic_json_dir
        else:
            raise ValueError("Dic json directory '%s' does not exists" %(dic_json_dir))
    
    @property
    def json_col_nums(self):
        return self._json_col_nums
    @json_col_nums.setter
    def json_col_nums(self, json_col_nums):
        if len(json_col_nums) == 3:
            self._json_col_nums = json_col_nums
        else:
            raise ValueError("Only three columns can be extracted from the json file")

    @property
    def dic_json_fname(self):
        return self._dic_json_fname
    @dic_json_fname.setter
    def dic_json_fname(self, dic_json_fname):
        if dic_json_fname.endswith('.json'):
            self._dic_json_fname = dic_json_fname
        else:
            self._dic_json_fname = dic_json_fname + '.json'

    @property
    def dic_par_dir(self):
        return self._dic_par_dir
    @dic_par_dir.setter
    def dic_par_dir(self, dic_par_dir):
        if os.path.exists(dic_par_dir):
            self._dic_par_dir = dic_par_dir
        else:
            raise ValueError("Dic par directory '%s' does not exists" %(dic_par_dir))
    
    @property
    def dic_par_fname(self):
        return self._dic_par_fname
    @dic_par_fname.setter
    def dic_par_fname(self, dic_par_fname):
        if dic_par_fname.endswith('.par'):
            self._dic_par_fname = dic_par_fname
        else:
            self._dic_par_fname = dic_par_fname + '.par'
    
    @property
    def img_dir(self):
        return self._img_dir
    @img_dir.setter
    def img_dir(self, img_dir):
        if os.path.exists(img_dir):
            self._img_dir = img_dir
        else:
            raise ValueError("Img directory '%s' does not exists" %(img_dir))
    
    @property
    def img_fname_template(self):
        return self._img_fname_template
    @img_fname_template.setter
    def img_fname_template(self, img_fname_template):
        self._img_fname_template = img_fname_template
    
    @property
    def output_dir(self):
        return self._output_dir
    @output_dir.setter
    def output_dir(self, output_dir):
        if os.path.exists(output_dir):
            self._output_dir = output_dir
        else:
            raise ValueError("Output directory '%s' does not exists" %(output_dir))
    
    @property
    def output_fname(self):
        return self._output_fname
    @output_fname.setter
    def output_fname(self, output_fname):
        if output_fname.endswith('.txt'):
            self._output_fname = output_fname
        else:
            self._output_fname = output_fname + '.txt'
    
    @property
    def first_img_num(self):
        return self._first_img_num
    @first_img_num.setter
    def first_img_num(self, first_img_num):
        self._first_img_num = first_img_num 
    
    # extra functions
    def open_dic_par_file(self, dic_mats):
        root = tk.Tk()
        root.withdraw()
        # user selects json file to use
        dic_json_dir = tk.filedialog.askopenfilename(initialdir=self._dic_json_dir,
                                                    defaultextension=".json",
                                                    filetypes=[("dic json files", "*.json"),
                                                               ("All Files", "*.*")],
                                                    title="Select dic.json File")
        # user selects par file to use
        dic_par_dir = tk.filedialog.askopenfilename(initialdir=self._dic_par_dir, 
                                                    defaultextension='.par',
                                                    filetypes=[("dic par files", "*.par"),
                                                               ("All Files", "*.*")],
                                                    title='Select dic.par File')
        if not dic_par_dir:
            quit()
        else:
            try:
                self.dic_json_dir = os.path.dirname(dic_json_dir)
                self.dic_json_fname = os.path.basename(dic_json_dir)
                json_file = open(dic_json_dir)
                json_dict = json.load(json_file)
                i = 0
                print("\nColumn and Key Pairs for dic.par")
                print(*json_dict.items())
                
                prompt_item_list = ['Image Number', 'Load', 'Screw Position']
                while 3 > i:
                    prompt = input("\nChose %s column to import. Press q to quit. \n" %(prompt_item_list[i]))
                    if input == "q":
                        break
                    self.json_col_nums[i] = int(prompt)
                    i+=1
                self.dic_par_dir = os.path.dirname(dic_par_dir)
                self.dic_par_fname = os.path.basename(dic_par_dir)
                dic_mats.process_dic_par_file(dic_par_dir, self.json_col_nums)
            except Exception as e:
                print(e)
                self.open_dic_par_file(dic_mats)
    
    def open_first_image(self):
        root = tk.Tk()
        root.withdraw()
        
        first_img_dir = tk.filedialog.askopenfilename(initialdir=self._img_dir, 
                                                    defaultextension='.tiff',
                                                    filetypes=[("TIFF Files", "*.tiff"),
                                                               ("TIFF Files", "*.tif"),
                                                               ("All Files", "*.*")],
                                                    title='Select First DIC Image File')
        
        if not first_img_dir:
            quit()
        else:
            try:
                self.img_dir = os.path.dirname(first_img_dir)
                first_img_fname = os.path.basename(first_img_dir)
                self.first_img_num = int((first_img_fname.split('_')[-1]).split('.')[0])
            except Exception as e:
                print(e)
                self.open_first_image()
    
    def get_img_num_dir(self, img_num):
        return os.path.join(self._img_dir, self._img_fname_template %(img_num))
    def get_first_img_dir(self):
        return self.get_img_num_dir(self._first_img_num)
    def get_dic_json_full_dir(self):
        return os.path.join(self._dic_json_dir, self._dic_json_fname)
    def get_json_col_nums(self):
        return self.json_col_nums
    def get_dic_par_full_dir(self):
        return os.path.join(self._dic_par_dir, self._dic_par_fname)
    def get_output_full_dir(self):
        return os.path.join(self._output_dir, self._output_fname)
    
    # str and rep
    def __repr__(self):
        return "dic_paths()\n" + self.__str__()
    def __str__(self):
        class_dict = {'base_dir' : self._base_dir,
        'dic_par_dir' : self._dic_par_dir,
        'dic_par_fname' : self._dic_par_fname,
        'img_dir' : self._img_dir,
        'img_fname_template' : self._img_fname_template,
        'output_dir' : self._output_dir,
        'output_fname' : self._output_fname,
        'first_img_num' : self._first_img_num}
        
        return str(class_dict)