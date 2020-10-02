function [image_struct, first_image_index, first_image_struct] = SelectFirstImage(file_mat, first_image_path)
% SelectFirstImage - Function used to extract the name of user specified 
%                    reference image
% 
%   AUTOR: Carman Fang, Chris Budrow, 
% 
%   USAGE:
% 
%   [image_struct, first_image_index, first_image_struct] = 
%       SelectFirstImage(file_mat, first_image_path)
% 
%   INPUT:
% 
%   file_mat is a n x 3 array
%      First column is a list of image numbers
%      Second Column is a list of forces corresponding to the first column
%      Third Column is a list of screw positions corresponding to the first column
%
%   first_image_path is a string variable
%      Contains the string directory path to DIC images
%
%   OUTPUT:
% 
%   image_struct is a structure
% 
%      .leader - Text leader before image number.  Found by searching the 
%          filename of the image selected, 'FirstImageName', for an '_'.
%      .extension - Extension of the image file.  Allows for selection of 
%          .tif or .tiff, typically equals '.tiff'.
% 
%   first_image_index is a string
%      Location (row number) of the find file that contains the selected
%      image name
% 
%   first_image_struct is a structure
%      .name - sting, name of selected image.
%      .folder - sting, folder location of the selected image.
%      .full_name - Full extension of the file (Folder and Name)
% 
%   NOTES:
%      Written and tested Spring 2017
% 
%      Problems occur when the image has multiple underscores ('_' in the 
%      file name.
% 
%
%   REVISION:
%      Carmen Fang, Chris Budrow (2017) - Originally written and tested
%
%      Tyrone Bracker, Dalton Shadle (2020) - Modified to include
%      first_image_path
%

drawnow

% First image information, only allow .tif and .tiff files to be selected
[first_image_struct.name, first_image_struct.folder] = uigetfile([first_image_path, '*.tiff'], 'Open First Image'); 
first_image_struct.full_name = [first_image_struct.folder first_image_struct.name];

% Ouput file name information as separate variables 
dash_pos = strfind(first_image_struct.name, '_');
dot_pos = strfind(first_image_struct.name, '.');

% Outputs
image_struct.leader = first_image_struct.name(1 : dash_pos);
image_struct.extension = first_image_struct.name(dot_pos : end);

% Find the location in 'File' of the selected image
image_num = first_image_struct.name(dash_pos+1 : dot_pos-1);
first_image_index = find(file_mat(:,1) == str2double(image_num));
end
