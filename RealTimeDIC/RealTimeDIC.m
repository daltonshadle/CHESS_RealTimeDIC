% Main function for DIC processing. Keeps running until "control+c" 
% Input: 1)user specified reference image
%        2)user specified area of interest and gridlines
%        3)user specified sample width and thickness for stress calculation
% Output:1)"gridx.dat, gridy.dat", XY coordinates specified by gridlines.
%        2)"StressStrain.dat", image number, stress value, strain value
%        3) plots of current stress and strain 

% Restart DIC processing when process has been interupted.  Will process
% all existing images then continue to look for changes in the par file to
% process any new images

% CONSTANT_NAMES = all uppercase with underscores
% FunctionNames() = camel case with parentheses and without underscores
% variable_names = all lowercase with underscores

%% start new envorinment
clc; clear;

%% user input paramters

% DIC Parameters
grid_spacing = 40; % spacing of correlation points
fixed_corr_dimen = [60, 120]; % (pixels) fixed correlation area dimensions for control point [height, width]
moving_sub_dimen = [40, 40]; % (pixels) moving subregion area dimensions for control point [height, width]
output_filename='dp_718_2_fatigue.txt';

% Sample Geometry
sample_width = 1; % (mm) cross-sectional width of the smaple for macro stress calculations
sample_thickness = 1; % (mm) cross-sectional thickness of the smaple for macro stress calculations

% Saving
save_stress_strain = true;

%% select the dic.par file for current sample, contains image index, force, and screw position
[par_file, file_mat] = ProcessDicParFile();

%% select the reference image
[image_struct, first_image_index, first_image_struct] = SelectFirstImage(file_mat, par_file.dir);

%% generate grid from the first reference image
[grid_x, grid_y, fixed_corr_dimen, moving_sub_dimen] = ...
    GenerateGrid(first_image_struct.full_name, grid_spacing, fixed_corr_dimen, moving_sub_dimen);

%% calculate stress and strain
% Need to process up until the second to last image. Then save all of that
% data off so it can be loaded in by the real time program.

% calculate smaple cross-sectional area for stress calculations
sample_area = sample_width*sample_thickness;

% index images that are already in dic.par to process first
first_img_idx = first_image_index;
last_img_idx = size(file_mat,1);
num_image_gap = 1;
image_indices = first_img_idx:num_image_gap:last_img_idx;
num_images = length(image_indices);

% initialize stress and strain matrices
strain_xx = zeros(num_images,1);
strain_yy = zeros(num_images,1);
strain_xy = zeros(num_images,1);
stress = zeros(num_images,1);

% for each file to process...
tic;
for i=1:num_images
    % process image index for full image and name
    curr_image = [image_struct.leader num2str(sprintf('%06d',file_mat(image_indices(i),1))) image_struct.extension];
    curr_image_full_name = [first_image_struct.folder curr_image];
    
    % special case for second image
    if i==1
        % process second image
        [valid_ref_x, valid_ref_y] = ProcessCorrelations(...
            first_image_struct.full_name, curr_image_full_name,...
            grid_x, grid_y, grid_x, grid_y, fixed_corr_dimen, moving_sub_dimen);  
        valid_x = valid_ref_x;
        valid_y = valid_ref_y;
    else
        % process the rest of the images
        [valid_x, valid_y] = ProcessCorrelations(...
            first_image_struct.full_name, curr_image_full_name, ...
            grid_x, grid_y, valid_x, valid_y, fixed_corr_dimen, moving_sub_dimen);
    end
    
    % do stress and strain calculations
    strain_xx(i) = FitStrain(valid_ref_x, (valid_x-valid_ref_x));
    strain_yy(i) = FitStrain(valid_ref_y, (valid_y-valid_ref_y));
    strain_xy(i) = (FitStrain(valid_ref_x, (valid_y-valid_ref_y))+FitStrain(valid_ref_y, (valid_x-valid_ref_x)))/2;
    stress(i) = file_mat(image_indices(i),2)/(sample_area);
    
    % write to output file
    output=[file_mat(image_indices(i),1),stress(i),strain_xx(i),strain_yy(i),strain_xy(i)];
    dlmwrite(output_filename,output,'delimiter','\t','precision','%.6f','-append')
    
    % display to screen
    disp(['File: ' num2str(i) ', Image: ' num2str(file_mat(image_indices(i),1))])
    disp(['Strain_XX: ' num2str(strain_xx(i))])
    disp(['Strain_YY: ' num2str(strain_yy(i))])
    disp(['Strain_XY: ' num2str(strain_xy(i))])
end
toc;

%% intial plotting
stress_strain_figure = figure();
plot(strain_xx,stress,'bd','MarkerFaceColor','b','MarkerEdgeColor','k');
xlabel('Macroscopic Strain')
ylabel('Macroscopic Stress (MPa)')
grid on

%% process the rest of the images in a while loop --> built-in while loop in "Continuous Update"
% ASSUMES THAT YOU HAVE ALREADY PROCESSED ALL THE DATA THAT YOU HAVE
hold on
drawnow
ContinuousUpdate(par_file.full_name, image_struct, grid_x, grid_y, first_image_struct, ...
                 valid_ref_x, valid_ref_y, valid_x, valid_y, sample_area, ...
                 output_filename, stress_strain_figure, ...
                 fixed_corr_dimen, moving_sub_dimen);
grid on

%% Save stress and strain values
if save_stress_strain
    save('strain_xx.mat', 'strain_xx')
    save('strain_xy.mat', 'strain_xy')
    save('strain_yy.mat', 'strain_yy')
    save('stress.mat', 'stress')
end
