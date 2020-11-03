function [fixed_corr_dimen, moving_sub_dimen] = ...
    AdjustCorrDimen(curr_image_name, grid_x, grid_y, old_fixed_corr_dimen, old_moving_sub_dimen)
% AdjustCorrDimen - Function used to adjust existing grid dimensions for points
%                of interest (seed points)
%
%   AUTHOR: Carman Fang, and code package
%
%   USAGE:
%
%   [fixed_corr_dimen, moving_sub_dimen] = ...
%        AdjustCorrDimen(curr_image_name, grid_x, grid_y, fixed_corr_dimen, moving_sub_dimen)
%
%   INPUT:
%
%   curr_image_name is a string
%      Complete name (including folder extension) of the current image.
%
%   grid_x is a 1 x n vector
%      List of x-coordinates of a seed box.
%
%   grid_y is a 1 x n vector
%      List of y-coordinates of a seed box.
%
%   old_fixed_corr_dimen is a 1 x 2 vector
%      (pixels) fixed correlation area dimensions for control point [height, width]
%
%   old_moving_sub_dimen is a 1 x 2 vector
%      (pixels) moving subregion area dimensions for control point [height, width]
%
%   OUTPUT:
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

% read reference image
funct_fig = figure('name','adjust_corr');
base_image = imread(curr_image_name);

% call select grid menu
[fixed_corr_dimen, moving_sub_dimen] = ...
    SelectGridType(base_image, grid_x, grid_y,...
    old_fixed_corr_dimen, old_moving_sub_dimen);
close(funct_fig);
end
% End GenerateGrid Function


% Select which type of grid you want to create
function [fixed_corr_dimen, moving_sub_dimen] = ...
    SelectGridType(base_image, grid_x, grid_y,...
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
        pause(0.5);
        
        if isempty(findobj('type','figure','name','adjust_corr'))
            is_btn_pressed = true;
            close_grid_menu = true;
        end
    end
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