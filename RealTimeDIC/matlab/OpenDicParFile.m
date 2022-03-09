function [file_mat] = OpenDicParFile(full_name)
% OpenDicParFile - Used to load .par file and store data in 'File'
% 
%   AUTHOR: Carmen Fang, Chris Budrow
% 
%   USAGE:
% 
%   [file_mat] = OpenDicParFile(full_name)
% 
%   INPUT:
% 
%   full_name is a string
%      Name of .par file to load.  Should include the path before the file
%      name.
%    
%   OUTPUT:
% 
%   file_mat is a n x 3 array
%      First column is a list of image numbers
%      Second Column is a list of forces corresponding to the first column
%      Third Column is a list of screw positions corresponding to the first column
% 
%   NOTES:
%      Written and tested Spring 2017
% 
%      This is the file that should be changed if the format of the par
%      file changes
%
%   REVISION:
%      Carmen Fang, Chris Budrow (2017) - Originally written and tested
%
%      Tyrone Bracker, Dalton Shadle (2020) - Modified to include screw
%      positions in return file
%

% open dic.par file and read info
fid = fopen(full_name);
%info = textscan(fid,'%s %s %s %s %s %f %f %f %s %s');
info = textscan(fid,'%s %s %s %f %f %f %f');
% save off image index (6), force (7), screw position (8)
file_mat = [info{4} info{6} info{5}];
fclose(fid);
end