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

from CPCorrFunctions import process_correlations

# ***************************************************************************
# CLASS DECLARATION
class dic_parameters_selector_widget():
    def __init__(self, dic_paths, dic_params, dic_mats, adjust_grid=True):
        
        self.dic_paths = dic_paths
        self.dic_params = dic_params
        self.dic_mats = dic_mats
        
        self.window = tk.Tk()
        self.window.geometry("1000x800")
        self.window.title('DIC Parameter Selector')
        self.first_img = cv2.imread(dic_paths.get_img_num_dir(dic_paths.first_img_num), 0)
        
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
                    if self.dic_mats.ref_points.size != 0:
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
                        
                        self.dic_params.grid_spacing = self.dic_mats.calc_grid_spacing()
                        self.bound_box_array = self.dic_mats.calc_bounding_box()
                        
                        update_sliders()
                        if self.dic_mats.ref_points.size != 0:
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
                self.dic_params.grid_spacing = [self.grid_width_slider.get(), 
                                                self.grid_height_slider.get()]
                if self.dic_mats.ref_points.size != 0:
                    self.update_plot()
            
            self.grid_width_slider = tk.Scale(self.window, label='Grid Spacing Width', 
                                              command=grid_spacing_slider_change,
                                              orient='horizontal',
                                              from_=5, to=200)
            self.grid_width_slider.place(x=820, y=350, height=50, width=160)
            
            self.grid_height_slider = tk.Scale(self.window, label='Grid Spacing Height', 
                                               command=grid_spacing_slider_change,
                                               orient='horizontal',
                                               from_=5, to=200)
            self.grid_height_slider.place(x=820, y=400, height=50, width=160)
        
        # add slider for fixed corr dimen
        def fixed_slider_change(event):
            self.dic_params.fixed_corr_dimen = [self.fixed_width_slider.get(), 
                                                self.fixed_height_slider.get()]
            if self.dic_mats.ref_points.size != 0:
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
            self.dic_params.moving_corr_dimen = [self.moving_width_slider.get(), 
                                                 self.moving_height_slider.get()]
            if self.dic_mats.ref_points.size != 0:
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
                self.grid_width_slider.set(self.dic_params.grid_spacing[0])
                self.grid_height_slider.set(self.dic_params.grid_spacing[1])
            self.fixed_width_slider.set(self.dic_params.fixed_corr_dimen[0])
            self.fixed_height_slider.set(self.dic_params.fixed_corr_dimen[1])
            self.moving_width_slider.set(self.dic_params.moving_corr_dimen[0])
            self.moving_height_slider.set(self.dic_params.moving_corr_dimen[1])
        
        update_sliders()
        
        if not adjust_grid and self.dic_mats.ref_points.size != 0:
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
        fcd = self.dic_params.fixed_corr_dimen
        mcd = self.dic_params.moving_corr_dimen
        gs = self.dic_params.grid_spacing
        
        # generate grid points
        x = np.linspace(np.min(bbc[:, 0]), np.max(bbc[:, 0]), 
                        int(np.abs(bbc[0, 0] - bbc[1, 0]) / gs[0]) + 1, 
                        endpoint=True)
        y = np.linspace(np.min(bbc[:, 1]), np.max(bbc[:, 1]), 
                        int(np.abs(bbc[0, 1] - bbc[1, 1]) / gs[1]) + 1, 
                        endpoint=True)
        
        xv, yv = np.meshgrid(x, y)
        grid_pts = np.vstack([xv.flatten(), yv.flatten()]).T.astype(float)
        self.dic_mats.ref_points = grid_pts
            
        # generate rectangle patches        
        fixed_rect = patches.Rectangle((grid_pts[0, 0] - fcd[0] / 2, grid_pts[0, 1] - fcd[1] / 2), 
                                 fcd[0], fcd[1], 
                                 linewidth=1, edgecolor='g', facecolor='none')
        
        moving_rect = patches.Rectangle((grid_pts[0, 0] - mcd[0] / 2, grid_pts[0, 1] - mcd[1] / 2), 
                                 mcd[0], mcd[1], 
                                 linewidth=1, edgecolor='c', facecolor='none')
        
        
        # Add the patch to the Axes
        self.first_img_ax.cla()
        self.fig.suptitle("Use the left mouse button to select two corners of a rectangular region of interest \n - Control Point Fixed Box in Green \n - Control Point Moving Search Box in Cyan")
        self.first_img_ax.imshow(self.first_img, cmap='Greys_r')
        self.first_img_ax.add_patch(fixed_rect)#, edgecolor='g')
        self.first_img_ax.add_patch(moving_rect)#, edgecolor='c')
        self.first_img_ax.scatter(grid_pts[:, 0], grid_pts[:, 1], c='b', s=15)
        self.canvas.draw()
    
    def get_all_dic_objects(self):
        return [self.dic_paths, self.dic_params, self.dic_mats]

#%%
class dic_continuous_update_widget():
    def __init__(self, dic_paths, dic_params, dic_mats, prev_img_num):
        
        self.dic_paths = dic_paths
        self.dic_params = dic_params
        self.dic_mats = dic_mats
        self.prev_img_num = prev_img_num
        
        self.window = tk.Tk()
        self.window.geometry("1000x800")
        self.window.title('RealTimeDIC Stress Strain Curve')
        
        self.fig = plt.figure(figsize=(6,6))
        self.stress_strain_ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.window)
        self.canvas.get_tk_widget().place(x=0, y=0, height=800, width=800)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self.window)
        
        
        # Add a button for reprocessing
        self.reprocess = False
        def repreocess_image_on_click():
            self.reprocess = True
        self.reprocess_button = tk.Button(self.window, text="Reprocess Image", command=repreocess_image_on_click)
        self.reprocess_button.place(x=820, y=150, height=40, width=160)
        
        # Add a button for resetting control point coordinates
        def reset_coords_on_click():
            self.dic_mats.reset_cur_points()
            self.reprocess = True
        self.reset_coords_button = tk.Button(self.window, text="Reset Grid Point Pos.", command=reset_coords_on_click)
        self.reset_coords_button.place(x=820, y=200, height=40, width=160)
        
        # Add a button for adjusting dic params
        def adjust_dic_params_on_click():
            dpsw = dic_parameters_selector_widget(self.dic_paths, self.dic_params, self.dic_mats, adjust_grid=False)
            [self.dic_paths, self.dic_params, self.dic_mats] = dpsw.get_all_dic_objects()
            self.reprocess = True
        self.adjust_dic_params_button = tk.Button(self.window, text="Adjust DIC Parameters", command=adjust_dic_params_on_click)
        self.adjust_dic_params_button.place(x=820, y=250, height=40, width=160)
        
        # Add a button for adjusting dic params
        def calc_backoff_stress_on_click():
            # process image index for full image and name
            cur = self.dic_mats.get_dic_par_mat_at_img_num(self.cur_img_num)
            cur_force = cur[1]
            
            # do stress and strain calculations and write to output
            [cur_stress, cur_strain_xx, cur_strain_yy, cur_strain_xy] = self.dic_mats.get_cur_stress_strain(cur_force, self.dic_params)
            
            # calc 90% backoff stress
            cur_backoff_stress = cur_stress * 0.9
            self.backoff_stress_label['text'] = '%0.1f MPa' %(cur_backoff_stress)
            
            
        self.backoff_stress_button = tk.Button(self.window, text="Calc Backoff Stress", command=calc_backoff_stress_on_click)
        self.backoff_stress_button.place(x=820, y=400, height=40, width=160)
        self.backoff_stress_label = tk.Label(self.window, text='- MPa', font=("Arial", 16))
        self.backoff_stress_label.place(x=820, y=450, height=40, width=160)
        
        self.update_plot()
        
        # Add a button for quitting
        def on_closing(root):
            MsgBox = tk.messagebox.askquestion ('Exit RealTimeDIC','Are you sure you want to exit RealTimeDIC? Output is saved, but all real-time processing will be lost.',icon = 'warning')
            if MsgBox == 'yes':
               root.destroy()
               root.quit()  
                
                      
        self.quit_button = tk.Button(self.window, text="Quit", command=lambda root=self.window:on_closing(root))
        self.quit_button.place(x=820, y=700, height=40, width=160)
        self.window.protocol("WM_DELETE_WINDOW", lambda root=self.window:on_closing(root))
        
        self.window.after(500, self.check_for_new_image)
        
        self.window.mainloop()
    
    def check_for_new_image(self):
        self.dic_mats.process_dic_par_file(self.dic_paths.get_dic_par_full_dir(), self.dic_paths.get_json_col_nums())
        self.cur_img_num = self.dic_mats.dic_par_mat[-1, 0]
        
        if (self.cur_img_num != self.prev_img_num) or (self.reprocess):
            self.process_img()
            self.reprocess = False
        
        self.window.after(500, self.check_for_new_image)
    
    def process_img(self):
        # process image index for full image and name
        cur = self.dic_mats.get_dic_par_mat_at_img_num(self.cur_img_num)
        cur_img_num = cur[0]
        cur_force = cur[1]
        cur_screw = cur[2]
        cur_img_dir = self.dic_paths.get_img_num_dir(self.cur_img_num)
        
        # process the rest of the images
        self.dic_mats.cur_points = process_correlations(self.dic_paths.get_first_img_dir(),
                                                         cur_img_dir, 
                                                         self.dic_mats.ref_points,
                                                         self.dic_mats.cur_points,
                                                         self.dic_params)
        
        # do stress and strain calculations and write to output
        [cur_stress, cur_strain_xx, cur_strain_yy, cur_strain_xy] = self.dic_mats.get_cur_stress_strain(cur_force, self.dic_params)
        cur_info = np.array([cur_img_num, 
                            cur_stress,
                            cur_strain_xx,
                            cur_strain_yy,
                            cur_strain_xy,
                            cur_screw])
        
        
        if self.reprocess:
            self.dic_mats.replace_in_output(cur_info, cur_img_num)
        else:
            self.dic_mats.add_to_output(cur_info)
        self.dic_mats.save_output_mat_to_file(self.dic_paths.get_output_full_dir())
        
        # display to screen
        if self.dic_params.is_sample_horizontal():
            print('DIC Image Number: %i \t Stress: %0.2f MPa \t Strain_XX: %0.3f %%' 
                  %(cur_img_num, cur_stress, cur_strain_xx * 100))
        else:
            print('DIC Image Number: %i \t Stress: %0.2f MPa \t Strain_YY: %0.3f %%' 
                  %(cur_img_num, cur_stress, cur_strain_yy * 100))
        
        self.update_plot()
        
        # update prev_img_num
        self.prev_img_num = self.cur_img_num
        
    def update_plot(self):
        # Add the patch to the Axes
        fontsize = 20
        self.stress_strain_ax.cla()
        # output = [img num, stress, strain_xx, strain_yy, strain_xy, screw]
        output = self.dic_mats.output_mat
        
        if self.dic_params.is_sample_horizontal():
            self.stress_strain_ax.scatter(output[:, 2] * 100, output[:, 1], c='b')
            self.fig.suptitle('Current Stress: %0.2f MPa     Current Strain_XX: %0.3f %%' 
                  %(output[-1, 1], output[-1, 2] * 100))
        else:
            self.stress_strain_ax.scatter(output[:, 3] * 100, output[:, 1], c='b')
            self.fig.suptitle('Current Stress: %0.2f MPa     Current Strain_YY: %0.3f %%' 
                  %(output[-1, 1], output[-1, 3] * 100), fontsize=fontsize)
                
        self.stress_strain_ax.set_xlabel('Macroscopic Strain (%)', fontsize=fontsize)
        self.stress_strain_ax.set_ylabel('Macroscopic Stress (MPa)', fontsize=fontsize)
        self.stress_strain_ax.tick_params(axis='both', which='major', labelsize=fontsize)
        
        self.canvas.draw()

