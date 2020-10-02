function [par_file, file_mat]=ProcessDicParFile()
% ProcessDicParFile - Used to load and extract data from .par file.
% 
%   AUTHOR:  Carmen Fang, Chris Budrow
% 
%   USAGE:
% 
%   [par_file, file_mat]=ProcessDicParFile()
% 
%   INPUT:
% 
%      None
%    
%   OUTPUT:
% 
%   par_file is a structure
%      .name - name of the file (ex: 'dic.par')
%      .dir - Directory the file is stored in.
%      .full_name - Full extension of the file (Directory and Name)
% 
%   file_mat is a n x 3 array
%      First column is a list of image numbers
%      Second Column is a list of forces corresponding to the first column
%      Third Column is a list of screw positions corresponding to the first column
% 
%   NOTES:
%      Written and tested Spring 2017
% 
%      This function process "dic.par" file to obtain a list of image names
%      and the force during loading
%
%   REVISION:
%      Carmen Fang, Chris Budrow (2017) - Originally written and tested
%
%      Tyrone Bracker, Dalton Shadle (2020) - Modified to include screw
%      positions in return file
%

% open UI file direcotry for selecting dic.par file
[par_file.name, par_file.dir] = uigetfile('*.par','Open dic.par File');
par_file.full_name = [par_file.dir par_file.name];

% Load in .par file and save information to File
[file_mat] = OpenDicParFile(par_file.full_name);
end



