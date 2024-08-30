function [output_spline_lengths, output_spline_wire_lengths, output_coordinates] = wire_function_v8(starting_points, ranges, points_changing, s, varargin)
% Authors:
% Fatima Ben Haj, Davis Senior High School
% Felipe Godinez, PhD, University of California Davis
% Date Created: September 2022
%{
% Code Last Modified: Jan 24, 2024

%%%%%<<<<<<<<<<<<<<<<<<<<<<REQUIRED INPUTS>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%starting_points: the original set of points, which should correspond with
%                 the electric field axis being passed in. If you want the
%                 spline to be displayed in the spatial field plot, the 
%                 points in the x direction should be between -2 and 2, the
%                 points in the y direction should be between -1 and 1, and 
%                 the points in the z direction should be between -1 and 1. 
%
%         
%         ranges: a matrix with the number of rows equal to the number of 
%                 fixed points and odd columns are minimum range and even 
%                 columns are the maximum range for each direction. 
%                 Ex: ranges = [[p1x_min, p1x_max;...
%                                p2x_min, p2x_max],...
%                               [p1y_min, p1y_max;...
%                                p2y_min, p2y_max],...
%                               [p1z_min, p1z_max;...
%                                p2z_min, p2z_max]]; 
%                 If the point does not move, the mins and maxs should be 
%                 the same.
% 
% 
%points_changing: a matrix that contains the indices of points that will be
%                 moving. For example: [3 4 5] means that the ranges of points 
%                 3, 4, and 5 will be used to move the spline. 
%
%              s: electric field data structure (should correspond with the
                  inputted points)

%%%%%<<<<<<<<<<<<<<<<<OPTIONAL INPUTS WITH VARARGIN>>>>>>>>>>>>>>>>>>>>>>>>
%      num_steps: the number of steps between ranges. For Example if there is a
%                 range [-1 2] and num_steps is 3, the function will use 3 
%                 values in between -1 and 2. Default value: 5
%
%     nom_length: the desired length of the wire in meters. The input set of 
%                 points will be modified to fit this length (only shortened). 
%                 The wire is modified by taking off points from the proximal 
%                 end (the side entering the body). The end modified lengths 
%                 should be within 0.2 of the input length, as shown by the 
%                 graph if show_graph is inputted. nom_length should always be 
%                 less than the unmodified length of the line because it will 
%                 be hard to keep added points on the trajectory. Default value: 2 
%
%     show_graph: show the lengths of the modified lines compared to the lengths
%                 of the original line on a histogram. 
%
%      show_plot: show the plotted splines (modified and unmodified) along with
%                 their lengths in the 3d. 
%
% intersecting_index: The index of the point that enters the body (which is
%                     shown as the pink line at y = 0.3 and z = 0.3 on the graph). 
%                     If this is inputed, the plot will show the length of the 
%                     modified wire outside of the body. If the values are negative, 
%                     it shows how far from the outside of the body the modified 
%                     line is. For example: if 2 is inputted, then the 2nd
%                     point of starting_points is where the wire enters the body. 
%                     Default value: 0 (the distance outside the body will also 
%                     be portrayed as 0 if there is no input or if the input is 0
%
%   give_warnings: shows the user a warning when their inputted spline points
%                  fall out of range of their inputted electric field axis values.
%
%       direction: The direction from which points are removed from the spline.
%                  This should be the proximal side of the wire. For
%                  example, if the user inputs "Right", it will be assumed that 
%                  the right side of the spline is the proximal end and
%                  points will be removed from it. Options are "Right" and
%                  "Left". Default value is "Left".
%
%   show_electric: This allows the user to turn off the electric field
%                  plot. The electric field points are still generated, just not 
%                  displayed.


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<OUTPUTS>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%          output_spline_lengths: all of the lengths of the unmodified wire
%         output_modified_spline_lengths: all of the lengths of the modified wire
%         output_coordinates: A matrix with 3 columns that contains the coordinates of all 
%                  the points on the possible modified lines that are converted to 
%                  electric field coordinates that correspond to the eField axes 
%                  that are passed in.
%}
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<CODE BODY>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Checks for correct number of inputs
if nargin < 4
error('wire_function_modified:NotEnoughInputs', "Requires at least 8 inputs")
end
if nargin > 22
error('wire_function_modified:TooManyInputs', "Requires at most 22 inputs")
end

% sets default values for varargin inputs
show_electric = true
num_points = length(starting_points);
num_steps = 5;
nom_length = 1;
show_graph = false;
show_plot = true;
intersecting_index = 0;
give_warnings = true;
eField0 = s.axis0;
eField1 = s.axis1;
eField2 = s.axis2;
eField = s.E;
direction = "Left";
show_electric = true;

% checks for optional inputs using varargin to set values to variables
for n=1:nargin-4
    if (strcmp(varargin{n}, 'num_steps'))
        num_steps = varargin{n+1};
    end

    if (strcmp(varargin{n}, 'nom_length'))
        nom_length = varargin{n+1};
    end

    if (strcmp(varargin{n}, 'show_graph'))
        show_graph = varargin{n+1};
    end

    if (strcmp(varargin{n}, 'show_plot'))
        show_plot = varargin{n+1};
    end

    if (strcmp(varargin{n}, 'intersecting_index'))
        intersecting_index = varargin{n+1}; 
    end

    if (strcmp(varargin{n}, 'give_warnings'))
        give_warnings = varargin{n+1}; 
    end

    if (strcmp(varargin{n}, 'direction'))
        direction = varargin{n+1}; 
    end

    if (strcmp(varargin{n}, 'show_electric'))
        show_electric = varargin{n+1}; 
    end

end
 
% set range of random points for each coordinate 
x_ranges = ranges(1:num_points, 1:2); 
y_ranges = ranges(1:num_points, 3:4);
z_ranges = ranges(1:num_points, 5:6);

% determine number of elements in the field for each dimension
nx = length(eField0(1:end-1));
ny = length(eField1(1:end-1));
nz = length(eField2(1:end-1));

% initiate variables
points = starting_points;
wire_lengths = [];
modified_lengths = [];
c = 1;

list_vals = [];

f1 = figure; %f1 is for the electric field plot
f2 = figure; %f2 is for the spatial field plot

%iterate over points, changing the point according to the ranges variable if it is in the input points_changing.
for x = 1:num_points
if ismember(x, points_changing)

    % compute 3 random matrices of length num_steps of points in between the 
    % mins and maxs of the inputted ranges. 
    current_xrange = x_ranges(x,1) + rand(1,num_steps)*(x_ranges(x,2)-x_ranges(x,1));
    current_yrange = y_ranges(x,1) + rand(1,num_steps)*(y_ranges(x,2)-y_ranges(x,1));
    current_zrange = z_ranges(x,1) + rand(1,num_steps)*(z_ranges(x,2)-z_ranges(x,1));

    % iterate over number of points in x direction
    for i=1:num_steps
        points(1, x) = current_xrange(i);
            
        % iterate over number of points in y direction
        for j = 1:num_steps
            points(2, x) = current_yrange(j);
            
            % iterate over number of points in z direction
            for k = 1:num_steps
                points(3, x) = current_zrange(k);
               
                % Compute trajectory/create a spline
                all_points = fnplt(cscvn(points));

                %interpolation of points on the wire. The default number of
                %points on the wire is 500

                %JUST COMMENTED
                pts_x = interp1(all_points(1,:), linspace(1, length(all_points), 500), "spline");
                pts_y = interp1(all_points(2,:), linspace(1, length(all_points), 500), "spline");
                pts_z = interp1(all_points(3,:), linspace(1, length(all_points), 500), "spline"); 

                all_points = [pts_x; pts_y; pts_z];

                %JUST COMMENTED BETWEEN THESE

                % compute length of inserted wire section
                distance_inside = getInsertedLength(all_points, points, intersecting_index);

                % calculate the length of the wire and add it to the
                % wire_lengths list
                wire_length = arclength(all_points(1,:), all_points(2,:), all_points(3,:));
                wire_lengths(c) = wire_length;
                c = c + 1;
                
                % Length modification -- shortens the spline to input_length
                [modified_points, modified_lengths(c)] = modifyWireLength(all_points, wire_length, nom_length);
                
                outside_length = modified_lengths(c) - distance_inside;

                % plotting the spline in the spatial field
                figure(f2);
                plotfunc(show_plot,all_points,points,modified_points, wire_length, modified_lengths(c),distance_inside, outside_length,num_steps)

                % extract E-field indices that correspond to wire trajectory 
                nn = length(modified_points);
                mm = 1:.025:nn;
                l = interp1(1:nn,modified_points',mm); 

                ddx = interp1(eField0(1:end-1),1:nx,l(:,1),'nearest');
                ddy = interp1(eField1(1:end-1),1:ny,l(:,2),'nearest');
                ddz = interp1(eField2(1:end-1),1:nz,l(:,3),'nearest');

                %warn user when points fall out side of electric field range
                if not(isnan(ddx) & isnan(ddy) & isnan(ddx))
                    list_vals{c} = [ddx, ddy, ddz]';
                    indice_all = [ddx ddy ddz];
                    eField(abs(eField) > 5) = 5; %changed from 5 to 7
                    vol = zeros(length(eField0)-1,length(eField1)-1,length(eField2)-1);
                    for nn=1:length(indice_all)
                        vol(indice_all(nn,1),indice_all(nn,2),indice_all(nn,3)) = 1;
                    end

                    % Display the electric field plot
                    if show_electric
                        t = tiledlayout(f1,2,1);
                        nexttile(t)
                        
                        %display the wire in the coronal view

                        coronal = imfuse (squeeze (sum (vol, 1) ), squeeze(abs(eField(:,75,:,1)) ), 'Scaling', 'independent') ;



                        imshow(coronal);
    
                        axis on
                        title('coronal view')
                        ylabel('Y')
                        xlabel('Z')
    
                        nexttile(t)
                       
                        %display the wire in the saggital view (on the same
                        %figure)
                        

                        saggital = imfuse (squeeze (sum (vol, 2) ), squeeze(abs(eField(75,:,:,1)) ), 'Scaling', 'independent');
                        imshow(saggital);
    
                        axis on
                        title('saggital')
                        ylabel('X')
                        xlabel('Z')
                    end

                    
                else 
                    % The user can choose to be warned if the inputted spatial points lie outside of the
                    % range of the electric field 
                    if (give_warnings)
                        warning("inputted spline points are outside the inputted electric field range.")
                    end
                end
               drawnow
            end
        end
    end 
end
end

%histogram
if show_graph
    edges=linspace(0.2,2,150); %sets the number of boxes 

    figure;
    histogram(wire_lengths, edges)
    hold on,
    histogram(modified_lengths, edges)
    legend('unmodified', 'modified');
    xlim([0.2 2]) %sets the bounds

    %axis labels
    xlabel('line lengths')
    ylabel('number of lines')
end

%set the outputs
output_spline_lengths = wire_lengths;
output_spline_wire_lengths = modified_lengths; 
output_coordinates = list_vals;


%% functions

%calculates the length of the wire in meters that is inside the body (past
%the point that enters the body) 
function [distance_inside] = getInsertedLength(all_modified_points, points, intersecting_index)
distance_inside = 0;  
if intersecting_index ~= 0
    points_for_coor = points'; 
    coordinate = points_for_coor(intersecting_index, :);
    coordinate_list = points_for_coor(intersecting_index, :);
    points_inside = all_modified_points';
    list_x = all_modified_points(1,:) - coordinate_list(1);
    list_y = all_modified_points(2,:) - coordinate_list(2);
    list_z = all_modified_points(3,:) - coordinate_list(3);
    [~, idx] = min(abs(list_x));
    [~, idy] = min(abs(list_y));
    [~, idz] = min(abs(list_z));

    coordinate = [all_modified_points(1,idx), all_modified_points(2,idy), all_modified_points(3,idz)];


    while points_inside(end, :) ~= coordinate 
        points_inside(end, :) = [];
    end
    distance_inside = arclength(points_inside(1,:), points_inside(2,:), points_inside(3,:));
    inside_points = points_inside;
    coordinate = points_for_coor(intersecting_index, :);
end
end

function [modified_points,new_length] = modifyWireLength(all_points, wire_length, nom_length)
% Subtract from the back end / outside the body (the dismal) by deleting points from
% the beginning of modified points. Check with a while loop to continue shortening until
% the wire is the right length. Plot the modified points. All wires should be within
% 0.02m of nom_length. 

rounded_length = round(wire_length, 1);
rounded_input = round(nom_length, 1);
modified_points = all_points;
p = 0;

%removes points depending on the direction that the user inputted
if direction == "Left"
    while rounded_length > rounded_input %while the line length is too long, remove points
        modified_points(:,1) = []; 
        rounded_length = round(arclength(modified_points(1,:), modified_points(2,:), modified_points(3,:)), 2);
        p = p + 1;
        if (p > 100)
            break;
        end
    end
    new_length = arclength(modified_points(1,:), modified_points(2,:), modified_points(3,:));
    end
if direction == "Right"
    while rounded_length > rounded_input 
        modified_points(:,end) = []; 
        rounded_length = round(arclength(modified_points(1,:), modified_points(2,:), modified_points(3,:)), 2);
        p = p + 1;
        if (p > 100)
            fprintf(sprintf('length of rounded wire %1.2f',rounded_length))
            break;
        end
    end
    new_length = arclength(modified_points(1,:), modified_points(2,:), modified_points(3,:));
else
    fprintf('Invalid input for direction\n');
end
end

%plot the spline in the sptial field
function plotfunc(show_plot,all_points,points,modified_points, wire_length,new_length,distance_inside, outside_length,num_steps)
if show_plot
    
    %Set the bounds for the axes 
    x_axis = [-2 2];
    y_axis = [-1 1];
    z_axis = [-1 1];
    
    %interpolate the spline points so that they are close to evenly spaced
    spline_x = interp1(all_points(1,:), linspace(1, length(all_points)), "spline");
    spline_y = interp1(all_points(2,:), linspace(1, length(all_points)), "spline");
    spline_z = interp1(all_points(3,:), linspace(1, length(all_points)), "spline");

    %plot the unmodified and modified splines with the control points
    spline = plot3(spline_z, spline_y, spline_x);
    spline.LineWidth = 1;
    hold on;
    plot3(points(3,:), points(2,:), points(1,:), 'ro');
    hold off;
    grid on;
    view(-50, 10)
    title('local wire MRI (units are in meters)')
    hold on;
    new_line = plot3(modified_points(3,:), modified_points(2,:), modified_points(1,:), ".",'markersize',12);
    set(new_line, "Color", '#EDB120')
    hold off;
    %display the lengths of the modified and unmodified splines 
    txt1 = {'Original length of wire:', num2str(wire_length)};
    txt = {'Modified length of wire:', num2str(new_length)};
    t = text(-0.8,0.9, 0.7,txt);
    t1 = text(-0.8, 0.9, 0.5, txt1);
    t.FontSize = 10;
    t1.FontSize = 10;
    set(t, 'Color', '#EDB120');
    set(t1, 'Color', 'b');
    txt2 = {'Distance of wire inside body:', num2str(distance_inside)};
    t2 = text(0.5, 0.3, 0.7, txt2);
    t2.FontSize = 10;
    txt3 = {'Distance of modified wire outside body:', num2str(outside_length)};
    t3 = text(0.5, 0.3, 0.5, txt3);
    t3.FontSize = 10;
    
    %plot the line that represents the body. 
    x_line = [-2,2];
    y_line = [0, 0];
    z_line = [0, 0];
    ll = line(x_line, y_line, z_line, 'linewidth', 3);
    set(ll, 'Color', 'm')
    xlim(x_axis)
    ylim (y_axis)
    zlim(z_axis)

    %adding legend and labels to axis. In the mri the z axis is the x axis and vice versa. The z
    %axis is the length of the cylinder/scanner bore
    xlabel('z')
    ylabel('y')
    zlabel('x')
    axis square  
    legend('original spline', 'control points', 'modified spline', 'body');
    
end
end

end