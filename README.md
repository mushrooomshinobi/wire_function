# wire_function
This is a function that generates wire trajectories in an electric field based on user inputs. For more information on the function see [citation]. 

Variables descriptions with examples: 

**Starting_points** contains a set of coordinates for a 3d spline, organized in 3 x the number of points array. It is the original set of points, which should correspond with the electric field spatial axis that are passed in. The function can display a 3d plot of the spline in the spatial field. The size of this plot is fixed so points may fall out of the range of view. This isn’t a problem for plotting on the electric field as long as starting_points correspond to the domain of the electric field volume. For the 3d spatial plot, the body (which is displayed as a pink line) is 2 meters and it is at a y and x value of 0. For reference, to make a line that stays within the spatial field axes, the x axis ranges from -1 to 1 m, the y axis ranges from -1 to 1 m, and the z axis ranges from -2 to 2 m.


**Ranges** contains the minimum and maximum values for each of the coordinates in the x, y, and z directions. The number of rows should equal the number of points. The odd columns are minimum values and the even columns are maximum values.  If the node points are to remain static, an integer placeholder value should be in the place of the range values. The ranges can be created by defining 3 subarrays (one for each of the 3 dimensions) which contain the ranges for each point and concatenating them into a single array. For example: x_ranges = [p1x_min, p1x_max; p2x_min, p2x_max],  y_ranges = [p1y_min, p1y_max; p2y_min, p2y_max], z_ranges = [p1z_min, p1z_max; p2z_min, p2 z_max], and then inputting the variable all_ranges = [x_ranges, y_ranges, z_ranges]. This is an example of how it would be set up if starting_points contained only 2 points:
              	Ex: ranges = [[p1x_min, p1x_max;...
                             	p2x_min, p2x_max],...
                            	[p1y_min, p1y_max;...
                             	p2y_min, p2y_max],...
                            	[p1z_min, p1z_max;...
                             	p2z_min, p2z_max]];
 

**Points_changing**, is a one-dimensional array of the indices (row numbers) of all the points (from starting_points) that will be moving. For example, if [3 4 5] is selected,  the ranges for the coordinates of starting_points at indices (rows) 3, 4, and 5 will be used to move the spline. 

We will continue with the optional inputs, which are inputted using varargin. The user inputs the name of the variable as a string with single quotation marks and the value of the variable in the next argument. 

Num_steps is the number of steps between ranges. For example if there is a range [-1 2] and num_steps is 3, the function will use 3 random (using the function rand()) values between -1 and 2. The default value for this variable is 5. 

The input nom_length holds the user’s desired length of the wire in meters. The inputted set of points will be modified to fit this length (only if it is too long it will be shortened). The default value for this variable is two meters. 

There are three optional boolean (true or false) variables. **show_graph** should be true (default is false) if the user wants a histogram of the lengths of the modified lines compared to the original lines displayed after the program finishes running. 

**show_plot** should be true (default is false) if the user wants to show the plotted splines (modified and unmodified) along with their lengths in the 3d. 

When the variable **show_warnings** is set to true (default is false), the program displays a warning to the user whenever their inputted spline points fall out of range of their inputted electric field axis. These cases are ignored when displaying the graphs.

**Intersecting_index** contains the index (of starting_points) where the wire enters the body. This is used to calculate the distance of the wire inside and outside the body. To ensure accuracy, the point at this index should be fixed. The default value is 0, in which the lengths (inside and outside) will be displayed as 0. 

**Direction**  is the direction from which points will be removed from the spline. This should be the proximal side of the wire. 

There are 3 outputs of the function besides the graphs. 

**output_wire_lengths** is a list of all the lengths of the unmodified spline,
**output_modified_wire_lengths** contains all the lengths of the modified wire, and 
**output_coordinates** is a matrix with 3 columns that contains all the points on the modified lines that have been converted to points in the electric field. 
