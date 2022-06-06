
function ContinuousUpdate(par_file_full_name, image_struct, ...
                      grid_x, grid_y, first_image_struct, valid_ref_x, valid_ref_y, ...
                      valid_x, valid_y, sample_area, ...
                      output_filename, plot_figure,...
                      fixed_corr_dimen, moving_sub_dimen)
% MassPlot - Function used to continuously read and plot DIC data
% 
%   AUTOR: Carman Fang, Chris Budrow
% 
%   USAGE:
% 
%   ContinuousUpdate(par_file_full_name, image_struct, ...
%                      grid_x, grid_y, first_image_struct, valid_ref_x, valid_ref_y, ...
%                      valid_x, valid_y, sample_area, ...
%                      output_filename, plot_figure,...
%                      fixed_corr_dimen, moving_sub_dimen)
% 
%   INPUT:
% 
%   par_file_full_name is a string
%      Full Name of .par file.
% 
%   image_struct is a structure
%      .leader - Text leader before image number.  Found by searching the 
%          filename of the image selected, 'FirstImageName', for an '_'.
%      .extension - Extension of the image file.  Allows for selection of 
%          .tif or .tiff, typically equals '.tiff'.
% 
%   grid_x is a 1 x n vector
%      List of x-coordinates of a seed box.  
% 
%   grid_y is a 1 x n vector
%      List of y-coordinates of a seed box.  
% 
%   first_image_struct is a structure
%      .name - sting, name of selected image.
%      .folder - sting, folder location of the selected image.
%      .full_name - Full extension of the file (Folder and Name)
% 
%   valid_ref_x is a 1 x n vector
%      List of x-coordinates of a reference box. Initial guess to aid in
%      corelation for first itteration of loop (Initial guess of valid_x).
% 
%   valid_ref_y is a 1 x n vector
%      List of y-coordinates of the reference box. Initial guess to aid in
%      corelation for first itteration of loop (Initial guess of valid_y).
% 
%   valid_x is a 1 x n vector
%      List of x-coordinates of a reference box. An initial guess to aid in
%      corelation.  Will change from image to image. 
% 
%   valid_y is a 1 x n vector
%      List of y-coordinates of a reference box. An initial guess to aid in
%      corelation.  Will change from image to image. 
%
%   sample_area is a double
%      Cross-sectional area of the sample in mm^2
%
%   corrsize is a double
%      edge size of correlation area
%
%   output_filename is a string
%      output filename for stress strain data
%
%   plot_figure is a figure object
%      figure associated with stress-strain curve from RealTimeDIC
%
%   fixed_corr_dimen is a 1 x 2 vector
%      (pixels) fixed correlation area dimensions for control point [height, width]
%
%   moving_sub_dimen is a 1 x 2 vector
%      (pixels) moving subregion area dimensions for control point [height, width]
%
%
% 
%   OUTPUT:
% 
%   None
% 
%   NOTES:
%      Written and tested Spring 2017
% 
%      Ouputs are saved in file: 'stressStrain.dat'
%        Image number, stress, Strain
% 
%      PostProcessing called here.
% 
%      Function to plot stress and strain in real time. No termination
%      action implemented, program should be stopped by "control+c".
% 
%      Compare newest image with previous newest image, if different, 
%      function will send the current newest to correlation, if not, 
%      function pauses for 1s then grabs the current newest image again.
%      This is done though seaching the last line of the .par file.
% 
%      loads in the file 'ForLoopData.mat' which has information from
%      previous load step.
% 

% START INIT CONTINUOUS ***************************************************
% init image string variables
pre_number_name = image_struct.leader;
image_extension = image_struct.extension;
first_image_full_name = first_image_struct.full_name;
image_folder = first_image_struct.folder;
first_img_index = first_image_struct.file_mat_index;

% init dic.par file information
[file_mat] = OpenDicParFile(par_file_full_name);
prev_index = size(file_mat,1);
run=1;

% init boundary storage structures
ref_struct = struct();
ref_struct.valid_x = valid_ref_x;
ref_struct.valid_y = valid_ref_y;
ref_struct.screw_pos = file_mat(first_img_index,3);
ref_struct.stress = 0;
ref_struct.strain_xx = 0;

ten_struct = struct();
ten_struct.valid_x = valid_x;
ten_struct.valid_y = valid_y;
ten_struct.screw_pos = file_mat(first_img_index,3);
ten_struct.stress = 0;
ten_struct.strain_xx = 0;
ten_struct.is_set = false;

comp_struct = struct();
comp_struct.valid_x = valid_x;
comp_struct.valid_y = valid_y;
comp_struct.screw_pos = file_mat(first_img_index,3);
comp_struct.stress = 0;
comp_struct.strain_xx = 0;
comp_struct.is_set = false;

upper_screw_bound = file_mat(first_img_index,3);
lower_screw_bound = file_mat(first_img_index,3);

hold on
fprintf('Ref Screw Pos: %0.6f\n', ref_struct.screw_pos);
% END INIT CONTINUOUS *****************************************************


% START UI CONTROLS *******************************************************
% initialize useful constants
BOUNDARY_MODE = 1;
INCREMENT_MODE = 0;

% initialize useful variables
mode = INCREMENT_MODE; % 0 for incremental mode, 1 for boundary mode

% create figure for plotting and UI control
set(gca,'units','normalized','position',[0.1, 0.1, 0.65, 0.8]);
set(gcf,'units','normalized','position',[0.1, 0.1, 0.8, 0.8]);

% ui location vectors
loc_btn_reprocess = [.8 .85 .15 .04];
loc_btn_reset_valid_xy = [.8 .8 .15 .04];
loc_btn_show_cycle = [.8 .75 .15 .04];
loc_btn_ten_tip = [.8 .7 .15 .04];
loc_txt_ten_tip = [.8 .6 .15 .09];
loc_btn_comp_tip = [.8 .55 .15 .04];
loc_txt_comp_tip = [.8 .45 .15 .09];
loc_btn_group_mode = [.8 .35 .15 .09];
loc_btn_calc_backoff = [0.775 0.20 0.20 0.10];
loc_slider_group = [0.775 0.04 0.20 0.15];


% Add ui widgets for reprocessing image
reprocess_btn_name = 'Reprocess Image';
reprocess_image = false;
btn_reporcess_image = uicontrol('Parent',plot_figure,...
                         'Style','pushbutton',...
                         'Units','normalized',...
                         'Position',loc_btn_reprocess,...
                         'String',reprocess_btn_name,...
                         'Callback', @BtnReprocessImageFunctionCB);
      function BtnReprocessImageFunctionCB(src, event)
          % set reprocess image boolean to true
          reprocess_image = true;
          disp('Reprocessing image...');
      end
  
% Add ui widgets for reseting valid_x and valid_y to reference
reset_btn_name = 'Reset Estimate';
btn_reset_valid_xy = uicontrol('Parent',plot_figure,...
                         'Style','pushbutton',...
                         'Units','normalized',...
                         'Position',loc_btn_reset_valid_xy,...
                         'String',reset_btn_name,...
                         'Callback', @BtnResetValidXYCB);
      function BtnResetValidXYCB(src, event)
          % set reset estimate for control point positions
          valid_x = valid_ref_x;
          valid_y = valid_ref_y;
          
          % reprocess image
          BtnReprocessImageFunctionCB(0, 0);
      end

% Add ui widget for cycling mode
show_cycle = true;
btn_hide_show_cycle = uicontrol('Parent',plot_figure,...
                        'Style','pushbutton',...
                        'Units','normalized',...
                        'Position',loc_btn_show_cycle,...
                        'String','Show Cycling Options',...
                        'Callback', @BtnShowHideCycleCB);
      function BtnShowHideCycleCB(src, event)
          % set variable and check what needs to be done
          if show_cycle
              btn_ten_tip.Visible = 'off';
              btn_comp_tip.Visible = 'off';
              txt_ten_tip.Visible = 'off';
              txt_comp_tip.Visible = 'off';
              
              btn_group_mode.Visible = 'off';
              
              btn_hide_show_cycle.String = 'Show Cycling Options';
          else
              btn_ten_tip.Visible = 'on';
              btn_comp_tip.Visible = 'on';
              txt_ten_tip.Visible = 'on';
              txt_comp_tip.Visible = 'on';
              
              % check if mode buttons can be displayed
              if comp_struct.is_set && ten_struct.is_set
                  btn_group_mode.Visible = 'on';
              end
              
              btn_hide_show_cycle.String = 'Hide Cycling Options';
          end
          
          show_cycle = ~show_cycle;
      end

% Add ui widgets for tension tip control
ten_tip_btn_name = 'Lock Tension Tip';
is_ten_tip_locked = false;
btn_ten_tip = uicontrol('Parent',plot_figure,...
                        'Style','pushbutton',...
                        'Units','normalized',...
                        'Position',loc_btn_ten_tip,...
                        'String',ten_tip_btn_name,...
                        'Callback', @BtnLockTenTipCB);
      function BtnLockTenTipCB(src, event)
          % check if tip is locked and if it needs to be reset
          if is_ten_tip_locked
              [ten_struct, upper_screw_bound, ten_str] = ...
                SaveTipData(valid_x, valid_y, curr_screw_pos, curr_stress, curr_strain_xx, ref_struct.screw_pos, txt_ten_tip);
              
              % set text for confirmation of recorded variables
              disp(['Ten Tip Locked @ ', ten_str]);
          end
          
          % lock tension tip
          is_ten_tip_locked = true;
          
          % set string to reset
          btn_ten_tip.String = 'Reset Tension Tip';
          
          % check if mode buttons can be displayed
          if comp_struct.is_set && ten_struct.is_set
              btn_group_mode.Visible = 'on';
          end
          
      end
txt_ten_tip = uicontrol('Parent',plot_figure,...
                        'Style','text',...
                        'Units','normalized',...
                        'Position',loc_txt_ten_tip,...
                        'String','');

% Add ui widgets for compression tip control
comp_tip_btn_name = 'Lock Compression Tip';
is_comp_tip_locked = false;
btn_comp_tip = uicontrol('Parent',plot_figure,...
                         'Style','pushbutton',...
                         'Units','normalized',...
                         'Position',loc_btn_comp_tip,...
                         'String',comp_tip_btn_name,...
                         'Callback', @BtnLockCompTipBtn);
      function BtnLockCompTipBtn(src, event)
          % check if tip is locked and if it needs to be reset
          if is_comp_tip_locked
              [comp_struct, lower_screw_bound, comp_str] = ...
                SaveTipData(valid_x, valid_y, curr_screw_pos, curr_stress, curr_strain_xx, ref_struct.screw_pos, txt_comp_tip);
              
              % set text for confirmation of recorded variables
              disp(['Comp Tip Locked @ ', comp_str]);
          end
          
          % lock compression tip
          is_comp_tip_locked = true;
          
          % set string to reset
          btn_comp_tip.String = 'Reset Compression Tip';
          
          % check if mode buttons can be displayed
          if comp_struct.is_set && ten_struct.is_set
              btn_group_mode.Visible = 'on';
          end
          
      end
txt_comp_tip = uicontrol('Parent',plot_figure,...
                        'Style','text',...
                        'Units','normalized',...
                        'Position',loc_txt_comp_tip,...
                        'String','');
  
% Add toggle buttons for mode control
boundary_btn_name = 'Boundary Mode'; 
incrememtal_btn_name = 'Incremental Mode';

% create button group 
btn_group_mode = uibuttongroup('Parent',plot_figure,...
                  'Visible','off',...
                  'Units','normalized',...
                  'Position',loc_btn_group_mode,...
                  'SelectionChangedFcn',@BtnSwitchModeCB);
              
% create two toggle buttons in the button group for switching modes
tbtn_incremental_mode = uicontrol(btn_group_mode,...
                  'Style','togglebutton',...
                  'String',incrememtal_btn_name,...
                  'Units','normalized',...
                  'Position',[.1 .6 .8 .3],...
                  'HandleVisibility','off');
tbtn_boundary_mode = uicontrol(btn_group_mode,...
                      'Style','togglebutton',...
                      'String',boundary_btn_name,...
                      'Units','normalized',...
                      'Position',[.1 .2 .8 .3],...
                      'HandleVisibility','off');

              
% Make the uibuttongroup visible after creating child objects
btn_group_mode.Visible = 'off';
set(btn_group_mode,'SelectedObject',tbtn_incremental_mode);
function BtnSwitchModeCB(src, ~)
   mode = tbtn_boundary_mode.Value;
   if mode == BOUNDARY_MODE
       mode_str = boundary_btn_name;
   else
       mode_str = incrememtal_btn_name;
   end
   disp(['Mode: ' mode_str]);
end

% create sliders group for handling fixed and moving rectangle changes
adjust_corr_btn_name = 'Adjust Correlation Dimensions';
group_sliders = uibuttongroup('Parent',plot_figure,...
                              'Units','normalized',...
                              'Position',loc_slider_group);

btn_adjust_corr_dimen = uicontrol(group_sliders,...
                                'style', 'pushbutton',...
                                'String', adjust_corr_btn_name,...
                                'units', 'normalized',...
                                'position', [0.1 0.70 0.80 0.2],...
                                'callback',@BtnAdjustCorrDimen);
txt_fixed_height = uicontrol(group_sliders,...
                                'style', 'text',...
                                'units', 'normalized',...
                                'position', [0.0 0.40 0.50 0.20],...
                                'String',sprintf('Fix Height: %i', fixed_corr_dimen(1)));
txt_fixed_width = uicontrol(group_sliders,...
                                'style', 'text',...
                                'units', 'normalized',...
                                'position', [0.5 0.40 0.50 0.20],...
                                'String',sprintf('Fix Width: %i', fixed_corr_dimen(2)));
txt_moving_height = uicontrol(group_sliders,...
                                'style', 'text',...
                                'units', 'normalized',...
                                'position', [0.0 0.10 0.50 0.20],...
                                'String',sprintf('Move Height: %i', moving_sub_dimen(1)));
txt_moving_width = uicontrol(group_sliders,...
                                'style', 'text',...
                                'units', 'normalized',...
                                'position', [0.5 0.10 0.50 0.20],...
                                'String',sprintf('Move Width: %i', moving_sub_dimen(2)));

function BtnAdjustCorrDimen(src, event)
    [fixed_corr_dimen, moving_sub_dimen] = ...
    AdjustCorrDimen(first_image_struct.full_name, grid_x, grid_y, fixed_corr_dimen, moving_sub_dimen);
    
    txt_fixed_height.String = sprintf('Fix Height: %i', fixed_corr_dimen(1));
    txt_fixed_width.String = sprintf('Fix Width: %i', fixed_corr_dimen(2));
    txt_moving_height.String = sprintf('Move Height: %i', moving_sub_dimen(1));
    txt_moving_width.String = sprintf('Move Width: %i', moving_sub_dimen(2));
end

% create button group for handling backoff stress calculations
backoff_calc_btn_name = 'Calculate 90% Backoff Stress';
group_backoff_btn = uibuttongroup('Parent',plot_figure,...
                              'Units','normalized',...
                              'Position',loc_btn_calc_backoff);

btn_backoff_calc = uicontrol(group_backoff_btn,...
                                'style', 'pushbutton',...
                                'String', backoff_calc_btn_name,...
                                'units', 'normalized',...
                                'position', [0.1 0.55 0.80 0.4],...
                                'callback',@BtnCalcBackoffStress);
txt_backoff_stress = uicontrol(group_backoff_btn,...
                                'style', 'text',...
                                'units', 'normalized',...
                                'position', [0.1 0.05 0.80 0.4],...
                                'String',sprintf('Backoff Stress: --- MPa'));

function BtnCalcBackoffStress(src, event)
    txt_backoff_stress.String = sprintf('Backoff Stress: %0.3f MPa', curr_stress * 0.9);
    disp(fprintf('Backoff Stress: %0.3f MPa', curr_stress * 0.9));
end



% run this once to set show cycle info
BtnShowHideCycleCB(0, 0);
% END UI CONTROLS *********************************************************

while(run==1)
    [file_mat] = OpenDicParFile(par_file_full_name);
    new_index = size(file_mat,1); 
    
    while(new_index == prev_index && ~reprocess_image)
        pause(0.2)
        [file_mat] = OpenDicParFile(par_file_full_name);
        prev_index = new_index;
        new_index = size(file_mat,1);
        curr_screw_pos = file_mat(new_index,3);
    end
    reprocess_image = false;
    
    new_image = [image_folder pre_number_name num2str(sprintf('%06d', file_mat(end,1))) image_extension];
    
    if mode == INCREMENT_MODE
        [valid_x,valid_y]=ProcessCorrelations(first_image_full_name,new_image,...
                           grid_x,grid_y,valid_x,valid_y,...
                           fixed_corr_dimen, moving_sub_dimen);
    else
        if curr_screw_pos > upper_screw_bound
            [valid_x,valid_y]=ProcessCorrelations(first_image_full_name,new_image,...
                                grid_x,grid_y,ten_struct.valid_x,ten_struct.valid_y,...
                                fixed_corr_dimen, moving_sub_dimen);
        elseif curr_screw_pos < lower_screw_bound
            [valid_x,valid_y]=ProcessCorrelations(first_image_full_name,new_image,...
                                grid_x,grid_y,comp_struct.valid_x,comp_struct.valid_y,...
                                fixed_corr_dimen, moving_sub_dimen);
        else
            [valid_x,valid_y]=ProcessCorrelations(first_image_full_name,new_image,...
                                grid_x,grid_y,ref_struct.valid_x,ref_struct.valid_y,...
                                fixed_corr_dimen, moving_sub_dimen);
        end
    end

    curr_strain_xx = FitStrain(valid_ref_x, (valid_x-valid_ref_x));
    curr_strain_yy = FitStrain(valid_ref_y, (valid_y-valid_ref_y));
    curr_strain_xy = (FitStrain(valid_ref_x, (valid_y-valid_ref_y))+FitStrain(valid_ref_y, (valid_x-valid_ref_x)))/2;
    curr_stress = file_mat(new_index,2)/(sample_area);
    
    % write to output file
    output=[file_mat(new_index,1), curr_stress, curr_strain_xx, curr_strain_yy, curr_strain_xy, curr_screw_pos];
    dlmwrite(output_filename,output,'delimiter','\t','precision','%.6f','-append');
    
    prev_index = new_index;
    drawnow
    plot(curr_strain_xx, curr_stress,'d','MarkerFaceColor','b','MarkerEdgeColor','k');
    title(sprintf('Current Stress = %0.3fMPa    |    Current Strain = %0.3f%%    |    Current Screw = %0.5f', curr_stress, curr_strain_xx*100, curr_screw_pos), 'FontSize', 15);
    disp(fprintf('Current Stress = %0.3fMPa | Current Strain = %0.3f%% | Current Screw = %0.5f', curr_stress, curr_strain_xx*100, curr_screw_pos));
    drawnow
    clear title  
    
    % check if tension tip is locked and update automatically
    if (curr_screw_pos > ten_struct.screw_pos) && (~is_ten_tip_locked)
        [ten_struct, upper_screw_bound, ~] = ...
            SaveTipData(valid_x, valid_y, curr_screw_pos, curr_stress, curr_strain_xx, ref_struct.screw_pos, txt_ten_tip);
    end
    
    % check if tension tip is locked and update automatically
    if (curr_screw_pos < comp_struct.screw_pos) && (~is_comp_tip_locked)
        [comp_struct, lower_screw_bound, ~] = ...
            SaveTipData(valid_x, valid_y, curr_screw_pos, curr_stress, curr_strain_xx, ref_struct.screw_pos, txt_comp_tip);
    end
end

end

function [tip_struct, tip_bound, tip_str] = ...
    SaveTipData(valid_x, valid_y, curr_screw_pos, curr_stress, curr_Strain_xx, ref_screw_pos, txt_box)
    % set variables for structure
    tip_struct.valid_x = valid_x;
    tip_struct.valid_y = valid_y;
    tip_struct.screw_pos = curr_screw_pos;
    tip_struct.stress = curr_stress;
    tip_struct.strain_xx = curr_Strain_xx;
    tip_struct.is_set = true;
    
    % calculate bound for screw position
    tip_bound = (tip_struct.screw_pos + ref_screw_pos) / 2;
    
    % update text each time
    tip_str = sprintf('Stress: %0.2f \n Strain: %0.3d \n Screw Pos: %0.5f \n Boundary Pos: %0.5f', ...
        tip_struct.stress, tip_struct.strain_xx, tip_struct.screw_pos, tip_bound);
    txt_box.String = tip_str;
end