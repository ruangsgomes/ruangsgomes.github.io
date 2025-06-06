function[h] = Fill_Ranges(x,y_low, y_up, color)
% Fill the ranges with a given color

% Input: x = vector with independent variables
%        y_low = vector with lower range o dependent variables
%        y_up  = vector with upper range o dependent variables
%        color = filling color

% Create a vector x1
X = [x(:); flipud(x(:)); x(1)];

% Corresponding matrix y1
Y = [y_up(:); flipud(y_low(:)); y_up(1)];

% Now fill area with "fill" function
h = fill(X,Y,color);

%set color
set (h, 'edgecolor',[0 0 0]);