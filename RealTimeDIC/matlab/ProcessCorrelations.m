function [new_valid_x, new_valid_y] = ProcessCorrelations(base_image, new_image,...
                              grid_x, grid_y, prev_valid_x, prev_valid_y,...
                              fixed_corr_dimen, moving_sub_dimen)
% ProcessCorrelations - This function selects and defines processing modes
%                       for the correlation.
% 
%   AUTOR: Carman Fang, Chris Budrow, modified from code package
% 
%   USAGE:
% 
%   [new_valid_x, new_valid_y] = ProcessCorrelations(base_image,...
%           new_image, grid_x, grid_y, prev_valid_x, prev_valid_y,...
%           fixed_corr_dimen, moving_sub_dimen)
% 
%   INPUT:
% 
%   base_image is a string
%      The reference image, the image that was selected.  This image will
%      be used as the zero strain comparative image
% 
%   new_image is a string
%      The image being correlated. The image where we want to know the
%      strain.
% 
%   grid_x is a 1 x n vector
%      List of x-coordinates of a seed box.  
% 
%   grid_y is a 1 x n vector
%      List of y-coordinates of a seed box.  
% 
%   prev_valid_x is a 1 x n vector
%      List of x-coordinates of a reference box. An initial guess to aid in
%      corelation.  Will change from image to image. (Initial guess of
%      NewValidX)
% 
%   prev_valid_y is a 1 x n vector
%      List of y-coordinates of a reference box. An initial guess to aid in
%      corelation.  Will change from image to image. (Initial guess of
%      NewValidY)
%
%   fixed_corr_dimen is a 1 x 2 vector
%      (pixels) fixed correlation area dimensions for control point [height, width]
%
%   moving_sub_dimen is a 1 x 2 vector
%      (pixels) moving subregion area dimensions for control point [height, width]
% 
%   OUTPUT:
% 
%   new_valid_x is a 1 x n vector
%      The x-coordinates that the seed point has translated to.
% 
%   new_valid_y is a 1 x n vector
%      The y-cooridnates that the seed point has translated to.
%      
%   Notes:
%      Written and tested Spring 2017
% 
%      No filter applied to images, base code that this was derived from
%      had the ability to apply a filter to the images.
%      Refer to first image as the baseline instead of previous image.
%        - Not sure what this means, it is a note from Carmen. - CJB
%        - I dont think this holds true any more
% 
%      GridX & GridY are equal to PreValidX & PreValidY for the first
%      correlation.  After that they can begin to deviate when there are
%      better guesses for where PreValidX & PreValidY should be
% 
%      Need to give 'CPcorr' input files
%
%   REVISION:
%      Carmen Fang, Chris Budrow (2017) - Originally written and tested
%
%      Tyrone Bracker, Dalton Shadle (2020) - Modified to include changing
%      correlation and subregion size
%


% assemble x,y points into one matrix for base and input
prev_points = [prev_valid_x, prev_valid_y];
base_points = [grid_x, grid_y];

% process image data type for base (reference) and input (current) images
data_field_type = GetImageDataType(base_image);
if (isempty(data_field_type))
    warning('Illegal bit depth');
    return
end
base_image = data_field_type(mean(double(imread(base_image)),3));
input_image = data_field_type(mean(double(imread(new_image)),3));

% Use  function 'CPcorr', this is the main correlation function.
% Derived from Matlabs 'cpcorr' function but allowing for larger scan area.
% It iwill be slower but results in better corelation.

InputCorr=CPcorr(round(prev_points),round(base_points),input_image,base_image,...
                 fixed_corr_dimen,moving_sub_dimen);
new_valid_x=InputCorr(:,1); 
new_valid_y=InputCorr(:,2);
end