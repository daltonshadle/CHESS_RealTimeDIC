#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 14:57:43 2022

@author: djs522
"""

# Main function for DIC processing. Keeps running until "control+c" 
# Input: 1)user specified reference image
#        2)user specified area of interest and gridlines
#        3)user specified sample width and thickness for stress calculation
# Output:1)"gridx.dat, gridy.dat", XY coordinates specified by gridlines.
#        2)"StressStrain.dat", image number, stress value, strain value
#        3) plots of current stress and strain 

# Restart DIC processing when process has been interupted.  Will process
# all existing images then continue to look for changes in the par file to
# process any new images

# CONSTANT_NAMES = all uppercase with underscores
# FunctionNames() = camel case with parentheses and without underscores
# variable_names = all lowercase with underscores

## user input paramters


#%% ***************************************************************************
# IMPORTS
import os
import numpy as np

import pickle

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
from matplotlib.widgets import Slider, Button
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

import pandas as pd

import cv2

import tkinter as tk

from CPCorrFunctions import cpcorr


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
    
    def calc_grid_spacing(self):
        if len(self.ref_points.shape) == 2:
            xs = np.unique(self.ref_points[:, 0])
            ys = np.unique(self.ref_points[:, 1])
            
            return [np.abs(xs[0] - xs[1]), np.abs(ys[0] - ys[1])]
    
    def calc_bounding_box(self):
        if len(self.ref_points.shape) == 2:
            return np.array([[np.min(self.ref_points[:, 0]), np.min(self.ref_points[:, 1])],
                             [np.max(self.ref_points[:, 0]), np.max(self.ref_points[:, 1])]])
        
    # str and rep
    def __repr__(self):
        print("dic_matrices()")
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
        print("dic_parameters()")
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
            except:
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
    
    # str and rep
    def __repr__(self):
        print("dic_paths()")
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

  




class dic_paramteres_selector_widget():
    def __init__(self, dic_paths, dic_params, dic_mats, adjust_grid=True):
        
        self.dic_paths = dic_paths
        self.dic_params = dic_params
        self.dic_mats = dic_mats
        
        self.window = tk.Tk()
        self.window.geometry("1000x800")
        
        self.first_img = cv2.imread(dic_paths.get_img_num_dir(dic_paths.get_first_img_num()), 0)
        
        self.fig = Figure(figsize=(6,6))
        self.fig.suptitle("Use the left mouse button to select a region of interest")
        self.first_img_ax = self.fig.add_subplot(111)
        self.first_img_ax.imshow(self.first_img, cmap='Greys_r')
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.window)
        self.canvas.get_tk_widget().place(x=0, y=0, height=800, width=800)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.window)
        
        # handle mouse cursor button for figure grid selection
        self.num_clicks = 0
        self.bound_box_array = np.zeros([2,2])
        def on_fig_cursor_press(event):
            if event.inaxes is not None and str(event.button) == 'MouseButton.LEFT' and str(self.toolbar.mode) != 'zoom rect':
                if self.num_clicks == 0:
                    self.bound_box_array[0, :] = [event.xdata, event.ydata]
                    self.first_img_ax.scatter(event.xdata, event.ydata, c='r', s=200, marker='+')
                    self.canvas.draw()
                    # update plot here with on point scatter
                else:
                    self.bound_box_array[1, :] = [event.xdata, event.ydata]
                    self.update_plot()
                    # update plot here with on point scatter and with grid
                
                self.num_clicks = (self.num_clicks + 1) % 2
        
        if adjust_grid:
            self.fig.canvas.mpl_connect('button_press_event', on_fig_cursor_press)
        
        
        
        
        # Add a button for loading dic paramters
        def load_params_on_click():
            print("Load DIC Parameters")
            root = tk.Tk()
            root.withdraw()
            
            dic_params_dir = tk.filedialog.askopenfilename(initialdir=self.dic_paths.get_output_dir(),
                                                           title='Select DIC Parameters File')
            
            if not dic_params_dir:
                pass
            else:
                try:
                    self.dic_paths.set_output_dir(os.path.dirname(dic_params_dir))
                    self.dic_params.load_dic_parameters_from_file(dic_params_dir)
                    
                    update_sliders()
                    if self.dic_mats.get_ref_points().size != 0:
                        self.update_plot()
                except:
                    print("Something failed when loading the params file")
                    pass
            
            # update sliders
        self.load_params_button = tk.Button(self.window, text="Load DIC Parameters", command=load_params_on_click)
        self.load_params_button.place(x=820, y=150, height=40, width=160)
        
        # Add a button for saving dic parameters
        def save_params_on_click():
            print("Save DIC Parameters")
            root = tk.Tk()
            root.withdraw()
            
            dic_params_dir_fname = tk.filedialog.asksaveasfilename(confirmoverwrite=False,
                                                                 initialdir=self.dic_paths.get_output_dir(),
                                                                 title='Save DIC Parameters')
            
            if not dic_params_dir_fname:
                pass
            else:
                try:
                    self.dic_paths.set_output_dir(os.path.dirname(dic_params_dir_fname))
                    self.dic_params.save_dic_parameters_to_file(dic_params_dir_fname)
                except:
                    print("Something failed when saving the params file")
                    pass
            
        self.save_params_button = tk.Button(self.window, text="Save DIC Parameters", command=save_params_on_click)
        self.save_params_button.place(x=820, y=200, height=40, width=160)
        
        # Add a button for loading grid points
        def load_grid_on_click():
            print("Load Grid")
            root = tk.Tk()
            root.withdraw()
            
            dic_grid_dir = tk.filedialog.askopenfilename(initialdir=self.dic_paths.get_output_dir(),
                                                           title='Select DIC Grid File')
            
            if not dic_grid_dir:
                pass
            else:
                try:
                    self.dic_mats.load_dic_matrices_from_file(dic_grid_dir)
                    
                    self.dic_params.set_grid_spacing(self.dic_mats.calc_grid_spacing())
                    self.bound_box_array = self.dic_mats.calc_bounding_box()
                    
                    update_sliders()
                    if self.dic_mats.get_ref_points().size != 0:
                        self.update_plot()
                except:
                    print("Something failed when loading grid file")
                    pass
            # update grid spacing and sliders
        self.load_grid_button = tk.Button(self.window, text="Load Grid", command=load_grid_on_click)
        self.load_grid_button.place(x=820, y=250, height=40, width=160)
        
        # Add a button for saving grid points
        def save_grid_on_click():
            print("Save Grid")
            root = tk.Tk()
            root.withdraw()
            
            dic_grid_dir_fname = tk.filedialog.asksaveasfilename(confirmoverwrite=False,
                                                                 initialdir=self.dic_paths.get_output_dir(),
                                                                 title='Save DIC Grid')
            
            if not dic_grid_dir_fname:
                pass
            else:
                try:
                    self.dic_paths.set_output_dir(os.path.dirname(dic_grid_dir_fname))
                    self.dic_mats.save_dic_matrices_to_file(dic_grid_dir_fname)
                except:
                    print("Something failed when saving the grid file")
                    pass
        self.save_grid_button = tk.Button(self.window, text="Save Grid", command=save_grid_on_click)
        self.save_grid_button.place(x=820, y=300, height=40, width=160)
        
        
        # add slider for grid spacing
        def grid_spacing_slider_change(event):
            self.dic_params.set_grid_spacing([self.grid_width_slider.get(), 
                                              self.grid_height_slider.get()])
            if self.dic_mats.get_ref_points().size != 0:
                self.update_plot()
        
        self.grid_width_slider = tk.Scale(self.window, label='Grid Spacing Width', 
                                          command=grid_spacing_slider_change,
                                          orient='horizontal',
                                          from_=20, to=200)
        self.grid_width_slider.place(x=820, y=350, height=50, width=160)
        
        self.grid_height_slider = tk.Scale(self.window, label='Grid Spacing Height', 
                                           command=grid_spacing_slider_change,
                                           orient='horizontal',
                                           from_=20, to=200)
        self.grid_height_slider.place(x=820, y=400, height=50, width=160)
        
        # add slider for fixed corr dimen
        def fixed_slider_change(event):
            self.dic_params.set_fixed_corr_dimen([self.fixed_width_slider.get(), 
                                                  self.fixed_height_slider.get()])
            if self.dic_mats.get_ref_points().size != 0:
                self.update_plot()
        
        self.fixed_width_slider = tk.Scale(self.window, label='Fixed Box Width', 
                                          command=fixed_slider_change,
                                          orient='horizontal',
                                          from_=20, to=200)
        self.fixed_width_slider.place(x=820, y=450, height=50, width=160)
        
        self.fixed_height_slider = tk.Scale(self.window, label='Fixed Box Height', 
                                           command=fixed_slider_change,
                                           orient='horizontal',
                                           from_=20, to=200)
        self.fixed_height_slider.place(x=820, y=500, height=50, width=160)
        
        # add slider for moving corr dimen
        def moving_slider_change(event):
            self.dic_params.set_moving_corr_dimen([self.moving_width_slider.get(), 
                                                   self.moving_height_slider.get()])
            if self.dic_mats.get_ref_points().size != 0:
                self.update_plot()
        
        self.moving_width_slider = tk.Scale(self.window, label='Moving Box Width', 
                                          command=moving_slider_change,
                                          orient='horizontal',
                                          from_=20, to=200)
        self.moving_width_slider.place(x=820, y=550, height=50, width=160)
        
        self.moving_height_slider = tk.Scale(self.window, label='Moving Box Height', 
                                           command=moving_slider_change,
                                           orient='horizontal',
                                           from_=20, to=200)
        self.moving_height_slider.place(x=820, y=600, height=50, width=160)
        
        
        def update_sliders():
            self.grid_width_slider.set(dic_params.get_grid_spacing()[0])
            self.grid_height_slider.set(dic_params.get_grid_spacing()[1])
            self.fixed_width_slider.set(dic_params.get_fixed_corr_dimen()[0])
            self.fixed_height_slider.set(dic_params.get_fixed_corr_dimen()[1])
            self.moving_width_slider.set(dic_params.get_moving_corr_dimen()[0])
            self.moving_height_slider.set(dic_params.get_moving_corr_dimen()[1])
        
        update_sliders()
        
        
        
        # Add a button for saving grid points
        def on_closing():
            self.window.destroy()
        def quit_on_click():
            print("Quit")
            on_closing()
            
        self.quit_button = tk.Button(self.window, text="Quit", command=quit_on_click)
        self.quit_button.place(x=820, y=700, height=40, width=160)
        
        self.window.protocol("WM_DELETE_WINDOW", on_closing)
        self.window.mainloop()
    
        
    def update_plot(self):
        # initialize local variables
        bbc = self.bound_box_array
        fcd = self.dic_params.get_fixed_corr_dimen()
        mcd = self.dic_params.get_moving_corr_dimen()
        gs = self.dic_params.get_grid_spacing()
        
        # generate grid points
        x = np.linspace(np.min(bbc[:, 0]), np.max(bbc[:, 0]), 
                        int(np.abs(bbc[0, 0] - bbc[1, 0]) / gs[0]) + 1, 
                        endpoint=True)
        y = np.linspace(np.min(bbc[:, 1]), np.max(bbc[:, 1]), 
                        int(np.abs(bbc[0, 1] - bbc[1, 1]) / gs[1]) + 1, 
                        endpoint=True)
        
        xv, yv = np.meshgrid(x, y)
        grid_pts = np.vstack([xv.flatten(), yv.flatten()]).T.astype(float)
        self.dic_mats.set_ref_points(grid_pts)
            
        # generate rectangle patches        
        fixed_rect = patches.Rectangle((grid_pts[0, 0] - fcd[0] / 2, grid_pts[0, 1] - fcd[1] / 2), 
                                 fcd[0], fcd[1], 
                                 linewidth=1, edgecolor='g', facecolor='none')
        
        moving_rect = patches.Rectangle((grid_pts[0, 0] - mcd[0] / 2, grid_pts[0, 1] - mcd[1] / 2), 
                                 mcd[0], mcd[1], 
                                 linewidth=1, edgecolor='c', facecolor='none')
        
        
        # Add the patch to the Axes
        self.first_img_ax.cla()
        self.fig.suptitle("Use the left mouse button to select a region of interest \n - Control Point Fixed Box in Green \n - Control Point Moving Search Box in Cyan")
        self.first_img_ax.imshow(self.first_img, cmap='Greys_r')
        self.first_img_ax.add_patch(fixed_rect)#, edgecolor='g')
        self.first_img_ax.add_patch(moving_rect)#, edgecolor='c')
        self.first_img_ax.scatter(grid_pts[:, 0], grid_pts[:, 1], c='b', s=15)
        self.canvas.draw()
        




base_dir = "/home/djs522/additional_sw/RealTimeDIC/CHESS_RealTimeDIC/example/dic_images/"
exp_dic_params = dic_parameters()
exp_dic_paths = dic_paths(base_dir=base_dir)
exp_dic_mats = dic_matrices()

exp_dic_paths.open_dic_par_file(exp_dic_mats)
exp_dic_paths.open_first_image()

dic_paramteres_selector_widget(exp_dic_paths, exp_dic_params, exp_dic_mats, adjust_grid=True)







#%%

class grid_selector_widget():
    # GenerateGrid - Function used to create or load existing grid for points
    #                of interest (seed points)
    #
    #   INPUT:
    #
    #   first_image_name is a string
    #      Complete name (including folder extension) of the first image.
    #
    #   grid_spacing is an integer
    #      Spacing between points of interest in the area of interst.
    #
    #   fixed_corr_dimen is a 1 x 2 vector
    #      (pixels) fixed correlation area dimensions for control point [height, width]
    #
    #   moving_sub_dimen is a 1 x 2 vector
    #      (pixels) moving subregion area dimensions for control point [height, width]
    #
    #   OUTPUT:
    #
    #   grid_x is a 1 x n vector
    #      List of x-coordinates of a seed box.  Will be written to the file:
    #      'gridx.dat'
    #
    #   grid_y is a 1 x n vector
    #      List of y-coordinates of a seed box.  Will be written to the file:
    #      'gridy.dat'
    #
    #   fixed_corr_dimen is a 1 x 2 vector
    #      (pixels) fixed correlation area dimensions for control point [height, width]
    #
    #   moving_sub_dimen is a 1 x 2 vector
    #      (pixels) moving subregion area dimensions for control point [height, width]
    #
    #   NOTES:
    #      Written and tested Spring 2017
    #
    #      'GridX' & 'GridY' can be loaded in from existing files rather than
    #      regenerated.
    #
    #      Always generates a rectangular grid.
    #
    #   REVISION:
    #      Carmen Fang (2017) - Originally written and tested
    #
    #      Tyrone Bracker and Dalton Shadle (2020) - Modified for displaying
    #      and adjusting grid spacing, correlation area, and subregion area
    #      sizes on-the-fly
    #
    #      Dalton Shadle (2022) - Converted to Python
    
    def __init__(self, first_img_dir, save_dir=None, dic_params=dic_parameters(), grid_coords=[], adjust_grid=True):
        
        self.dic_params = dic_params
        self.grid_coords = grid_coords
        self.img = cv2.imread(first_img_dir, 0)
        self.fig = plt.figure(figsize=plt.figaspect(0.5))
        self.gs = GridSpec(12,3)
        self.img_ax = self.fig.add_subplot(self.gs[:,:2])
        self.save_dir = save_dir
        
        self.img_ax.imshow(self.img, cmap='Greys_r')
        
        if adjust_grid:
            self.fig.canvas.mpl_connect('button_press_event', self.onpress)
            self.fig.canvas.mpl_connect('button_release_event', self.onrelease)
        self.fig.canvas.mpl_connect('close_event', self.on_close)
        
        # Add a button for loading grid
        load_grid_button = Button(self.fig.add_subplot(self.gs[3,2]), 'Load Grid', hovercolor='0.975')
        def load_grid_button_on_click(mouse_event):
            print('To Be Implemented')
        load_grid_button.on_clicked(load_grid_button_on_click)
        
        # Add a button for generating grid
        save_grid_button = Button(self.fig.add_subplot(self.gs[4,2]), 'Save Grid', hovercolor='0.975')
        def save_grid_button_on_click(mouse_event):
            if self.save_dir is not None:
                np.save(os.path.join(self.save_dir, 'grid_coords.npy'), self.grid_coords)
            plt.close(self.fig)
        save_grid_button.on_clicked(save_grid_button_on_click)
        
        
        # Define an axes area and draw a slider in it
        grid_spacing_slider = Slider(self.fig.add_subplot(self.gs[5,2]), 'Spacing', 5, 200, 
                                     valinit=self.dic_params.get_grid_spacing(), valfmt='%i')
        fixed_x_slider = Slider(self.fig.add_subplot(self.gs[6,2]), 'Fixed X', 40, 200, 
                                valinit=self.dic_params.get_fixed_corr_dimen()[0], valfmt='%i')
        fixed_y_slider = Slider(self.fig.add_subplot(self.gs[7,2]), 'Fixed Y', 40, 200, 
                                valinit=self.dic_params.get_fixed_corr_dimen()[1], valfmt='%i')
        move_x_slider = Slider(self.fig.add_subplot(self.gs[8,2]), 'Move X', 20, 150, 
                               valinit=self.dic_params.get_moving_corr_dimen()[0], valfmt='%i')
        move_y_slider = Slider(self.fig.add_subplot(self.gs[9,2]), 'Move Y', 20, 150, 
                               valinit=self.dic_params.get_moving_corr_dimen()[1], valfmt='%i')
        
        # Define an action for modifying the line when any slider's value changes
        def grid_spacing_slider_change(val):
            self.dic_params.set_grid_spacing(int(grid_spacing_slider.val))
            self.update_plot()
        def fixed_slider_change(val):
            self.dic_params.set_fixed_corr_dimen([int(fixed_x_slider.val), int(fixed_y_slider.val)])
            self.update_plot()
        def moving_slider_change(val):
            self.dic_params.set_moving_corr_dimen([int(move_x_slider.val), int(move_y_slider.val)])
            self.update_plot()
        
        grid_spacing_slider.on_changed(grid_spacing_slider_change)
        fixed_x_slider.on_changed(fixed_slider_change)
        fixed_y_slider.on_changed(fixed_slider_change)
        move_x_slider.on_changed(moving_slider_change)
        move_y_slider.on_changed(moving_slider_change)
        
        if len(self.grid_coords) == 0:
            self.bound_box_coords = np.zeros([2, 2]) 
        else:
            self.bound_box_coords = np.vstack([self.grid_coords[0, :], self.grid_coords[-1, :]])
            self.update_plot()
        
        plt.show()
    
    def onpress(self, event):
        if str(event.button) == 'MouseButton.RIGHT':
            ix, iy = event.xdata, event.ydata
            self.bound_box_coords[0, :] = [ix, iy]
            
    def onrelease(self, event):
        if str(event.button) == 'MouseButton.RIGHT':
            ix, iy = event.xdata, event.ydata
            self.bound_box_coords[1, :] = [ix, iy]
            self.bound_box_coords = self.bound_box_coords.astype(int)
            
            self.update_plot()
    
    def on_close(self, event):
        np.save(os.path.join(self.save_dir, 'grid_coords.npy'), self.grid_coords)
    
    def update_plot(self):
        bbc = self.bound_box_coords
        fcd = self.dic_params.get_fixed_corr_dimen()
        mcd = self.dic_params.get_moving_corr_dimen()
        
        # generate grid points
        x = np.linspace(np.min(bbc[:, 0]), np.max(bbc[:, 0]), 
                        int(np.abs(bbc[0, 0] - bbc[1, 0]) / self.dic_params.get_grid_spacing()) + 1, 
                        endpoint=True)
        y = np.linspace(np.min(bbc[:, 1]), np.max(bbc[:, 1]), 
                        int(np.abs(bbc[0, 1] - bbc[1, 1]) / self.dic_params.get_grid_spacing()) + 1, 
                        endpoint=True)
        
        xv, yv = np.meshgrid(x, y)
        xv = xv.astype(int)
        yv = yv.astype(int)
        
        self.grid_coords = np.vstack([xv.flatten(), yv.flatten()]).T.astype(float)
            
        # generate rectangle patches
        bounding_rect = patches.Rectangle((bbc[:, 0].min(), bbc[:, 1].min()), 
                                 np.abs(bbc[0, 0] - bbc[1, 0]), np.abs(bbc[0, 1] - bbc[1, 1]), 
                                 linewidth=1, edgecolor='r', facecolor='none')
        
        fixed_rect = patches.Rectangle((self.grid_coords[0, 0] - fcd[0] / 2, self.grid_coords[0, 1] - fcd[1] / 2), 
                                 fcd[0], fcd[1], 
                                 linewidth=1, edgecolor='g', facecolor='none')
        
        moving_rect = patches.Rectangle((self.grid_coords[0, 0] - mcd[0] / 2, self.grid_coords[0, 1] - mcd[1] / 2), 
                                 mcd[0], mcd[1], 
                                 linewidth=1, edgecolor='c', facecolor='none')
        
        
        # Add the patch to the Axes
        self.img_ax.cla()
        self.img_ax.imshow(self.img, cmap='Greys_r')
        self.img_ax.add_patch(bounding_rect)
        self.img_ax.add_patch(fixed_rect)
        self.img_ax.add_patch(moving_rect)
        self.img_ax.scatter(self.grid_coords[:, 0], self.grid_coords[:, 1], c='b')
        self.fig.canvas.draw_idle()
    
    def get_dic_params(self):
        return self.dic_params
    
    def get_grid_coords(self):
        return self.grid_coords

class continuous_update_widget():
    def __init__(self, window, dic_par_dir, dic_params, sample_area,
                 img_dir_template, first_img_idx, prev_img_idx, 
                 ref_coords, curr_coords, output, output_dir):
        
        self.dic_par_dir = dic_par_dir
        self.dic_params = dic_params
        self.sample_area = sample_area
        self.img_dir_template = img_dir_template
        self.first_img_idx = first_img_idx
        self.prev_img_idx = prev_img_idx
        self.ref_coords = ref_coords
        self.curr_coords = curr_coords
        self.output = output
        self.output_dir = output_dir
        
        self.window = window
        self.window.geometry("1000x800")
        
        # Add a button for reprocessing
        self.reprocess = False
        def repreocess_image_on_click():
            self.reprocess = True
        self.reprocess_button = tk.Button(window, text="Reprocess Image", command=repreocess_image_on_click)
        self.reprocess_button.place(x=820, y=100, height=40, width=160)
        
        # Add a button for resetting control point coordinates
        def reset_coords_on_click():
            self.curr_coords = self.ref_coords
            self.reprocess = True
        self.reset_coords_button = tk.Button(window, text="Reset CP Coordinates", command=reset_coords_on_click)
        self.reset_coords_button.place(x=820, y=150, height=40, width=160)
        
        # Add a button for adjusting dic params
        def adjust_dic_params_on_click():
            t = grid_selector_widget(self.img_dir_template %(self.dic_par_mat[self.first_img_idx, 0]), 
                                     dic_params=self.dic_params, grid_coords=self.ref_coords, adjust_grid=False)    
            self.dic_params = t.get_dic_params()
            self.reprocess = True
        self.adjust_dic_params_button = tk.Button(window, text="Adjust Correlation Dimensions", command=adjust_dic_params_on_click)
        self.adjust_dic_params_button.place(x=820, y=200, height=40, width=160)
        
        self.fig = plt.figure(figsize=(6,6))
        self.stress_strain_ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.window)
        self.canvas.get_tk_widget().place(x=0, y=0, height=800, width=800)
        
        self.update_plot()
        
        self.window.after(500, self.check_for_new_image)
    
    def check_for_new_image(self):
        self.dic_par_mat = open_dic_par_file(self.dic_par_dir)
        self.curr_img_idx = self.dic_par_mat.shape[0] - 1
        
        if (self.curr_img_idx != self.prev_img_idx) or (self.reprocess):
            self.process_img()
            self.reprocess = False
        
        self.window.after(500, self.check_for_new_image)
    
    def process_img(self):
        curr_img_num, curr_force, curr_screw = self.dic_par_mat[self.curr_img_idx, :]
                
        curr_coords = process_correlations(self.img_dir_template %(self.dic_par_mat[self.first_img_idx, 0]), 
                                           self.img_dir_template %(curr_img_num), 
                                           self.ref_coords, 
                                           self.curr_coords, 
                                           self.dic_params)
        
        # do stress and strain calculations and write to output
        curr_stress = curr_force / self.sample_area
        curr_strain_xx = fit_strain(ref_coords[:, 0], (curr_coords[:, 0]-ref_coords[:, 0]))
        curr_strain_yy = fit_strain(ref_coords[:, 1], (curr_coords[:, 1]-ref_coords[:, 1]))
        curr_strain_xy = (fit_strain(ref_coords[:, 0], (curr_coords[:, 1]-ref_coords[:, 1])) + fit_strain(ref_coords[:, 1], (curr_coords[:, 0]-ref_coords[:, 0]))) / 2
        
        # add to output
        curr_info = np.array([curr_img_num, curr_stress, curr_strain_xx,
                              curr_strain_yy, curr_strain_xy, curr_screw])
        if self.reprocess:
            self.output[self.curr_img_idx, :] = curr_info
        else:
            self.output = np.vstack([self.output, curr_info])
        np.savetxt(self.output_dir, self.output, fmt='%0.6f', delimiter='\t')
        
        # display to screen
        print('DIC Image Number: %i \t Stress: %0.2f MPa \t Strain_XX: %0.3f %%' 
              %(curr_img_num, self.output[self.curr_img_idx, 1], self.output[self.curr_img_idx, 2] * 100))
        
        self.update_plot()
        
        # update prev_img_idx
        self.prev_img_idx = self.curr_img_idx
        
    def update_plot(self):
        # Add the patch to the Axes
        self.stress_strain_ax.cla()
        self.stress_strain_ax.scatter(self.output[:, 2] * 100, self.output[:, 1], c='b');
        self.stress_strain_ax.set_xlabel('Macroscopic Strain (%)')
        self.stress_strain_ax.set_ylabel('Macroscopic Stress (MPa)')
        
        self.canvas.draw()


#*****************************************************************************
#%% FUNCTION DECLARATIONS

def process_dic_par_file(base_dir=None):
    # process_dic_par_file - Used to load and extract data from .par file.
    #   INPUT:
    # 
    #   base_dir is a string
    #      path of the base directory (default = None)
    #    
    #   OUTPUT:
    # 
    #   par_file is a structure
    #      .name - name of the file (ex: 'dic.par')
    #      .dir - Directory the file is stored in.
    #      .full_name - Full extension of the file (Directory and Name)
    # 
    #   file_mat is a n x 3 array
    #      First column is a list of image numbers
    #      Second Column is a list of forces corresponding to the first column
    #      Third Column is a list of screw positions corresponding to the first column
    # 
    #   NOTES:
    #      Written and tested Spring 2017
    # 
    #      This function process "dic.par" file to obtain a list of image names
    #      and the force during loading
    #
    #   REVISION:
    #      Carmen Fang, Chris Budrow (2017) - Originally written and tested
    #
    #      Tyrone Bracker, Dalton Shadle (2020) - Modified to include screw
    #      positions in return file
    #
    #      Dalton Shadle (2022) - Converted to Python
    
    if base_dir is not None:
        base_dir = os.getcwd()
    
    root = tk.Tk()
    root.withdraw()
    
    dic_par_dir = tk.filedialog.askopenfilename(initialdir=base_dir, defaultextension='.par', title='Select dic.par File')
    
    return [open_dic_par_file(dic_par_dir), dic_par_dir]

def open_dic_par_file(dic_par_dir):    
    # day, month, day number, time, year, image number, load, screw position, extra
    df = pd.read_csv(dic_par_dir, sep=" ", header=None)
    dic_par_mat = np.array(df)
    
    return dic_par_mat[:, [5, 6, 7]]


def select_first_image(dic_par_mat, base_dir=None):
    # SelectFirstImage - Function used to extract the name of user specified 
    #                    reference image
    # 
    #   INPUT:
    # 
    #   file_mat is a n x 3 array
    #      First column is a list of image numbers
    #      Second Column is a list of forces corresponding to the first column
    #      Third Column is a list of screw positions corresponding to the first column
    #
    #   first_image_path is a string variable
    #      Contains the string directory path to DIC images
    #
    #   OUTPUT:
    # 
    #   image_struct is a structure
    # 
    #      .leader - Text leader before image number.  Found by searching the 
    #          filename of the image selected, 'FirstImageName', for an '_'.
    #      .extension - Extension of the image file.  Allows for selection of 
    #          .tif or .tiff, typically equals '.tiff'.
    # 
    #   first_image_index is a string
    #      Location (row number) of the find file that contains the selected
    #      image name
    # 
    #   first_image_struct is a structure
    #      .name - sting, name of selected image.
    #      .folder - sting, folder location of the selected image.
    #      .full_name - Full extension of the file (Folder and Name)
    # 
    #   NOTES:
    #      Written and tested Spring 2017
    # 
    #      Problems occur when the image has multiple underscores ('_' in the 
    #      file name.
    # 
    #
    #   REVISION:
    #      Carmen Fang, Chris Budrow (2017) - Originally written and tested
    #
    #      Tyrone Bracker, Dalton Shadle (2020) - Modified to include
    #      first_image_path
    #
    #      Dalton Shadle (2022) - Converted to Python

    
    # First image information, only allow .tif and .tiff files to be selected
    if base_dir is not None:
        base_dir = os.getcwd()
    
    root = tk.Tk()
    root.withdraw()
        
    first_img_dir = tk.filedialog.askopenfilename(initialdir=base_dir, defaultextension='.tiff', title='Select First DIC Image')
    
    # process first image
    first_img_fname = first_img_dir.split('/')[-1]
    img_dir = first_img_dir.replace(first_img_fname, '')
    first_img_num = int((first_img_fname.split('_')[-1]).split('.')[0])
    first_img_idx = np.where(dic_par_mat[:, 0] == first_img_num)[0]
    
    return [img_dir, first_img_fname, first_img_num, first_img_idx]  


def process_correlations(ref_img_dir, curr_img_dir, grid_coords, 
                         prev_valid_coords, dic_params):
    # ProcessCorrelations - This function selects and defines processing modes
    #                       for the correlation.
    # 
    #   INPUT:
    # 
    #   ref_image_dir is a string
    #      The reference image, the image that was selected.  This image will
    #      be used as the zero strain comparative image
    # 
    #   curr_image_dir is a string
    #      The image being correlated. The image where we want to know the
    #      strain.
    # 
    #   grid_x is a 1 x n vector
    #      List of x-coordinates of a seed box.  
    # 
    #   grid_y is a 1 x n vector
    #      List of y-coordinates of a seed box.  
    # 
    #   prev_valid_x is a 1 x n vector
    #      List of x-coordinates of a reference box. An initial guess to aid in
    #      corelation.  Will change from image to image. (Initial guess of
    #      NewValidX)
    # 
    #   prev_valid_y is a 1 x n vector
    #      List of y-coordinates of a reference box. An initial guess to aid in
    #      corelation.  Will change from image to image. (Initial guess of
    #      NewValidY)
    #
    #   fixed_corr_dimen is a 1 x 2 vector
    #      (pixels) fixed correlation area dimensions for control point [height, width]
    #
    #   moving_corr_dimen is a 1 x 2 vector
    #      (pixels) moving subregion area dimensions for control point [height, width]
    # 
    #   OUTPUT:
    # 
    #   curr_valid_x is a 1 x n vector
    #      The x-coordinates that the seed point has translated to.
    # 
    #   curr_valid_y is a 1 x n vector
    #      The y-cooridnates that the seed point has translated to.
    #      
    #   Notes:
    #      Written and tested Spring 2017
    # 
    #      No filter applied to images, base code that this was derived from
    #      had the ability to apply a filter to the images.
    #      Refer to first image as the baseline instead of previous image.
    #        - Not sure what this means, it is a note from Carmen. - CJB
    #        - I dont think this holds true any more
    # 
    #      GridX & GridY are equal to PreValidX & PreValidY for the first
    #      correlation.  After that they can begin to deviate when there are
    #      better guesses for where PreValidX & PreValidY should be
    # 
    #      Need to give 'CPcorr' input files
    #
    #   REVISION:
    #      Carmen Fang, Chris Budrow (2017) - Originally written and tested
    #
    #      Tyrone Bracker, Dalton Shadle (2020) - Modified to include changing
    #      correlation and subregion size
    #
    #      Dalton Shadle (2022) - Converted to Python
    
    # read in images
    ref_img = cv2.imread(ref_img_dir, 0)
    curr_img = cv2.imread(curr_img_dir, 0)
    
    # Use  function 'CPcorr', this is the main correlation function.
    # Derived from Matlabs 'cpcorr' function but allowing for larger scan area.
    # It will be slower but results in better corelation.
    [curr_valid_coords, StdX, StdY, CorrCoef, errorInfos] = cpcorr(np.copy(prev_valid_coords), 
                                                           np.copy(grid_coords), 
                                                           curr_img, ref_img, 
                                                           dic_params.get_fixed_corr_dimen(), 
                                                           dic_params.get_moving_corr_dimen())
    
    return curr_valid_coords

def fit_strain(ref_coords, curr_disp):
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

#%%

# paths and filenames
base_dir = "/home/djs522/additional_sw/RealTimeDIC/CHESS_RealTimeDIC/example/"
output_filename = "dp718-1_nov2020.txt"
img_template = 'dic_%06i.tiff'

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















