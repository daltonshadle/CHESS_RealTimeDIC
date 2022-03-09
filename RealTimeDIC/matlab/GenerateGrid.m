function [grid_x, grid_y, fixed_corr_dimen, moving_sub_dimen] = ...
    GenerateGrid(first_image_name, grid_spacing, fixed_corr_dimen, moving_sub_dimen)
% GenerateGrid - Function used to create or load existing grid for points
%                of interest (seed points)
%
%   AUTHOR: Carman Fang, and code package
%
%   USAGE:
%
%   [grid_x, grid_y, fixed_corr_dimen, moving_sub_dimen] = ...
%       GenerateGrid(first_image_name, grid_spacing, fixed_corr_dimen,moving_sub_dimen)
%
%   INPUT:
%
%   first_image_name is a string
%      Complete name (including folder extension) of the first image.
%
%   grid_spacing is an integer
%      Spacing between points of interest in the area of interst.
%
%   fixed_corr_dimen is a 1 x 2 vector
%      (pixels) fixed correlation area dimensions for control point [height, width]
%
%   moving_sub_dimen is a 1 x 2 vector
%      (pixels) moving subregion area dimensions for control point [height, width]
%
%   OUTPUT:
%
%   grid_x is a 1 x n vector
%      List of x-coordinates of a seed box.  Will be written to the file:
%      'gridx.dat'
%
%   grid_y is a 1 x n vector
%      List of y-coordinates of a seed box.  Will be written to the file:
%      'gridy.dat'
%
%   fixed_corr_dimen is a 1 x 2 vector
%      (pixels) fixed correlation area dimensions for control point [height, width]
%
%   moving_sub_dimen is a 1 x 2 vector
%      (pixels) moving subregion area dimensions for control point [height, width]
%
%   NOTES:
%      Written and tested Spring 2017
%
%      'GridX' & 'GridY' can be loaded in from existing files rather than
%      regenerated.
%
%      Always generates a rectangular grid.
%
%   REVISION:
%      Carmen Fang (2017) - Originally written and tested
%
%      Tyrone Bracker and Dalton Shadle (2020) - Modified for displaying
%      and adjusting grid spacing, correlation area, and subregion area
%      sizes on-the-fly

% initialize GridX and GridY
grid_x=[];
grid_y=[];

% read reference image
base_image = imread(first_image_name);

% call select grid menu
[grid_x, grid_y, fixed_corr_dimen, moving_sub_dimen] = ...
    SelectGridType(first_image_name, base_image, grid_x, grid_y, grid_spacing,...
    fixed_corr_dimen,moving_sub_dimen);
close all;
end
% End GenerateGrid Function


% Select which type of grid you want to create
function [grid_x, grid_y, fixed_corr_dimen, moving_sub_dimen] = ...
    SelectGridType(base_image_name, base_image, grid_x, grid_y, grid_spacing,...
    fixed_corr_dimen, moving_sub_dimen)

% display grid generator menu and process choices
close_grid_menu = false;
all_rects = false;

% reset plotting area by setting hold off and plotting grid
hold off
set(gca,'units','normalized','position',[0.05, 0.1, 0.65, 0.8]);
set(gcf,'units','normalized','position',[0.1, 0.1, 0.8, 0.8]);
imshow(base_image,'InitialMagnification',100);
axis on
hold on
    
while(~close_grid_menu)
    % reset all plots
    delete(findall(gca, 'type', 'Scatter'));
    delete(findall(gca, 'type', 'Rectangle'));
    delete(findall(gcf,'type','annotation'));
    
    % check if GridX and GridY have been intialized, if so draw correlation region
    if ~isempty(grid_x)
        scatter(grid_x,grid_y,'+r');

        % reshape grid points into nx2 matrix and calculate rectangles
        grid_xy = [reshape(grid_x,[],1),reshape(grid_y,[],1)];
        corr_rects = CalcRects(grid_xy,fixed_corr_dimen);
        sub_rects = CalcRects(grid_xy,moving_sub_dimen);

        % all_rects = true: displays correlation region for ALL control points
        % all_rects = false: displays correlation region for ONE control point
        if all_rects
            for i = 1:size(corr_rects,1)
                rectangle('Position', corr_rects(i,:), 'EdgeColor','g');
                rectangle('Position', sub_rects(i,:), 'EdgeColor','b');
            end
        else
            rectangle('Position', corr_rects(1,:), 'EdgeColor','g');
            rectangle('Position', sub_rects(1,:), 'EdgeColor','b');
        end
        str = sprintf('Red = Control Points\nGreen = Correlation Area (Fixed)\nBlue = Subregion Area (Move)');
        annotation('textbox',[0.4, 0.75, 0.1, 0.1],'String',str,'FitBoxToText','on','Units','normalized','Color','white');
    end  
    
    % create ui controls
    btn_load_grid = uicontrol('style', 'pushbutton',...
                              'String', 'Load Grid',...
                              'units', 'normalized',...
                              'position', [0.75 0.85 0.15 0.04],...
                              'callback', @BtnLoadGridCB);
                          
    btn_rect_grid = uicontrol('style', 'pushbutton',...
                              'String', 'Create Grid',...
                              'units', 'normalized',...
                              'position', [0.75 0.80 0.15 0.04],...
                              'callback', @BtnRectGridCB);

    btn_save_grid = uicontrol('style', 'pushbutton',...
                              'String', 'Save Grid',...
                              'units', 'normalized',...
                              'position', [0.75 0.75 0.15 0.04],...
                              'callback', @BtnSaveGridCB);

    btn_close = uicontrol('style', 'pushbutton',...
                          'String', 'Close',...
                          'units', 'normalized',...
                          'position', [0.75 0.65 0.15 0.04],...
                          'callback', @BtnCloseCB);  
    
    % toggle button for showing all correlation areas
    tbtn_disp_rects = uicontrol('style', 'togglebutton',...
                              'String', 'Show All Rectangles',...
                              'units', 'normalized',...
                              'position', [0.75 0.70 0.15 0.04],...
                              'callback', @TBtnDispRects);
    
    % create sliders group for handling fixed and moving rectangle changes
    group_sliders = uibuttongroup('Visible','off',...
                                  'Units','normalized',...
                                  'Position',[0.675 0.30 0.3 0.25]);
    slider_offset = 10;             
    slider_fixed_height = uicontrol(group_sliders,...
                                    'style', 'slider',...
                                    'units', 'normalized',...
                                    'position', [0.05 0.70 0.6 0.19],...
                                    'min', moving_sub_dimen(1) + slider_offset,...
                                    'max', 250,...
                                    'value', fixed_corr_dimen(1),...
                                    'callback',@SliderFixedHeightCB);
    txt_fixed_height = uicontrol(group_sliders,...
                                    'style', 'text',...
                                    'units', 'normalized',...
                                    'position', [0.7 0.70 0.3 0.19],...
                                    'String',sprintf('Fixed H (Green): %i', fixed_corr_dimen(1)));
    slider_fixed_width = uicontrol(group_sliders,...
                                    'style', 'slider',...
                                    'units', 'normalized',...
                                    'position', [0.05 0.50 0.6 0.19],...
                                    'min', moving_sub_dimen(2) + slider_offset,...
                                    'max', 250,...
                                    'value', fixed_corr_dimen(2),...
                                    'callback',@SliderFixedWidthCB);
    txt_fixed_width = uicontrol(group_sliders,...
                                    'style', 'text',...
                                    'units', 'normalized',...
                                    'position', [0.7 0.50 0.3 0.19],...
                                    'String',sprintf('Fixed W (Green): %i', fixed_corr_dimen(2)));
    slider_moving_height = uicontrol(group_sliders,...
                                    'style', 'slider',...
                                    'units', 'normalized',...
                                    'position', [0.05 0.30 0.6 0.19],...
                                    'min', 10,...
                                    'max', fixed_corr_dimen(1)-slider_offset,...
                                    'value', moving_sub_dimen(1),...
                                    'callback',@SliderMovingHeightCB);
    txt_moving_height = uicontrol(group_sliders,...
                                    'style', 'text',...
                                    'units', 'normalized',...
                                    'position', [0.7 0.30 0.3 0.19],...
                                    'String',sprintf('Move H (Blue): %i', moving_sub_dimen(1)));
    slider_moving_width = uicontrol(group_sliders,...
                                    'style', 'slider',...
                                    'units', 'normalized',...
                                    'position', [0.05 0.10 0.6 0.19],...
                                    'min', 10,...
                                    'max', fixed_corr_dimen(2)-slider_offset,...
                                    'value', moving_sub_dimen(2),...
                                    'callback',@SliderMovingWidthCB);
    txt_moving_width = uicontrol(group_sliders,...
                                    'style', 'text',...
                                    'units', 'normalized',...
                                    'position', [0.7 0.10 0.3 0.19],...
                                    'String',sprintf('Move W (Blue): %i', moving_sub_dimen(2)));
                                
    
                                
    % handle visibility of ui controls
    if isempty(grid_x)
        btn_save_grid.Visible = 'off';
        btn_close.Visible = 'off';
        group_sliders.Visible = 'off';
        tbtn_disp_rects.Visible = 'off';
    else
        btn_save_grid.Visible = 'on';
        btn_close.Visible = 'on';
        group_sliders.Visible = 'on';
        tbtn_disp_rects.Visible = 'on';
    end
    
    % GUI while loop to keep processing until ui control is selected and processed
    is_btn_pressed = false;
    while(~is_btn_pressed)
        pause(0.2);
    end
end

function BtnLoadGridCB(src, event)
    drawnow
    [grid_x_name, grid_x_path] = uigetfile('*.dat','Open gridx.dat');
    if grid_x_name == 0
        disp('You did not select a file!')
    end
    grid_x = importdata([grid_x_path grid_x_name],'\t');
    drawnow
    [grid_y_name, grid_y_path] = uigetfile('*.dat','Open gridy.dat');
    if grid_y_name == 0
        disp('You did not select a file!')
    end
    grid_y = importdata([grid_y_path grid_y_name],'\t');
    
    is_btn_pressed = true;
end

function BtnRectGridCB(src, event)
    
    [grid_x,grid_y] = CreateRectangularGrid(grid_spacing);
    
    is_btn_pressed = true;
end

function BtnSaveGridCB(src, event)
    SaveGrid(grid_x,grid_y);
    disp('Saved Grid Data');
    
    is_btn_pressed = true;
end

function BtnCloseCB(src, event)
    close_grid_menu = true;
    
    is_btn_pressed = true;
end

function SliderFixedHeightCB(src, event)
    fixed_corr_dimen(1) = round(get(slider_fixed_height,'value'));
    
    is_btn_pressed = true;
end

function SliderFixedWidthCB(src, event)
    fixed_corr_dimen(2) = round(get(slider_fixed_width,'value'));
    
    is_btn_pressed = true;
end

function SliderMovingHeightCB(src, event)
    moving_sub_dimen(1) = round(get(slider_moving_height,'value'));
    
    is_btn_pressed = true;
end

function SliderMovingWidthCB(src, event)
    moving_sub_dimen(2) = round(get(slider_moving_width,'value'));
    
    is_btn_pressed = true;
end

function TBtnDispRects(src, event)
    all_rects = ~all_rects;
    
    is_btn_pressed = true;
end

end
% End SelectGridType Function

% Select a rectangular area
function [grid_x, grid_y] = CreateRectangularGrid(grid_spacing)

% Select rect (lower left and upper right)
title(sprintf('Define the region of interest. Pick (single click) a point in the LOWER LEFT region of the gage section.\n  Do the same for a point in the UPPER RIGHT portion of the gage section.'));

% draw lines for selecting region of interest
[x(1,1), y(1,1)] = ginput(1);
hold on
scatter(x(1,1), y(1,1),'+b');
[x(2,1), y(2,1)] = ginput(1);
hold on
scatter(x(2,1), y(2,1),'+b');
drawnow
x_min = min(x);
x_max = max(x);
y_min = min(y);
y_max = max(y);

x_spacing = grid_spacing;
y_spacing = grid_spacing;

% Round XMin,XMax and YMin,YMax "up" based on selected spacing
num_x_elements = ceil((x_max-x_min)/x_spacing)-1;
num_y_elements = ceil((y_max-y_min)/y_spacing)-1;
x_min_new = (x_max+x_min)/2-((num_x_elements/2)*x_spacing);
x_max_new = (x_max+x_min)/2+((num_x_elements/2)*x_spacing);
y_min_new = (y_max+y_min)/2-((num_y_elements/2)*y_spacing);
y_max_new = (y_max+y_min)/2+((num_y_elements/2)*y_spacing);

% Create the analysis grid and show user
[x,y] = meshgrid(x_min_new:x_spacing:x_max_new,y_min_new:y_spacing:y_max_new);

% Overwrite GridX and GridY with X and Y
grid_x = reshape(x,[],1);
grid_y = reshape(y,[],1);

% remove blue plot points
title(sprintf('Region defined with %i control points.', size(grid_x,1)));
delete(findall(gca, 'type', 'scatter'));
end
% End CreateRectangularGrid Function

% Save grid (GridX, GridY)
function SaveGrid(grid_x, grid_y)
save gridx.dat grid_x -ascii -tabs
save gridy.dat grid_y -ascii -tabs
end
% End SaveGrid Function

function rect = CalcRects(xy, rect_dimen)
% Calculate rectangles so imcrop will return image with xy coordinate
% inside center pixel

% xy specifies center of rectangle, need upper left
default_height = rect_dimen(1);
default_width = rect_dimen(2);
upperleft = round(xy) - [default_width/2, default_height/2];

% break apart upper and left coords and add in heights and widths
upper = upperleft(:,2);
left = upperleft(:,1);
width = default_width * ones(size(upper));
height = default_height * ones(size(upper));

% return left and upper coords and width and height of each rectangle
% centered at xy
rect = [left upper width height];
end
% End CalcRects Function