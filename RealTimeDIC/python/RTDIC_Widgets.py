# -*- coding: utf-8 -*-
#%% ***************************************************************************
# IMPORTS
import os
import numpy as np

import tkinter as tk

import cv2

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

# ***************************************************************************
# CLASS DECLARATION
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
        
        
        self.bound_box_array = self.dic_mats.calc_bounding_box()
        
        # handle mouse cursor button for figure grid selection
        if adjust_grid:
            self.num_clicks = 0
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
        
        
        if adjust_grid:
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
            if adjust_grid:
                self.grid_width_slider.set(dic_params.get_grid_spacing()[0])
                self.grid_height_slider.set(dic_params.get_grid_spacing()[1])
            self.fixed_width_slider.set(dic_params.get_fixed_corr_dimen()[0])
            self.fixed_height_slider.set(dic_params.get_fixed_corr_dimen()[1])
            self.moving_width_slider.set(dic_params.get_moving_corr_dimen()[0])
            self.moving_height_slider.set(dic_params.get_moving_corr_dimen()[1])
        
        update_sliders()
        
        if not adjust_grid and self.dic_mats.get_ref_points().size != 0:
            self.update_plot()
        
        # Add a button for quitting
        def on_closing(root):
            root.destroy()
            root.quit()            
        self.quit_button = tk.Button(self.window, text="Quit", command=lambda root=self.window:on_closing(root))
        self.quit_button.place(x=820, y=700, height=40, width=160)
        self.window.protocol("WM_DELETE_WINDOW", lambda root=self.window:on_closing(root))
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
    
    def get_all_dic_objects(self):
        return [self.dic_paths, self.dic_params, self.dic_mats]


#%%
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


