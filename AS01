% J.Williams
% University of Wisconsin-Madison
% EPD629: Powertrain Systems and Controls
% Assignement 01
% 2022-01-26



% Q1:
% Find intersection of two lines (2 points)
% 
% Define slope and y-intercept for line 1 as: m1 = 2 and b1 = 10
% Define x-breakpoints for line 1 as: xVector1 = -10:10;
% Calculate corresponding y points as: yVector1 = m1 * xVector1 + b1;
% Plot line 1.
% 
% Define slope and y-intercept for line 2 as: m2 = -1 and b2 = 20
% Define x-breakpoints for line 2 as: xVector2 = -10:10;
% Calculate corresponding y points as: yVector2 = m2 * xVector2 + b2;
% Plot line 2 (overlay line 2 on top of plot for line 1 using the “hold on” command prior to the second plot command).
% 
% These two lines can be represented in matrix form as:
% A * z = B
% or
% -m11-m21 * xy = b1b2
% 
% Create matrix A and B in your script
% 
% Find intersection of these two lines using Cramer’s rule:
% xIntersect = det([b1 1;b2 1]) / det(A)
% yIntersect = det([-m1 b1;-m2 b2]) / det(A)
% 
% Add the intersection point to your plot using the following command
% 	plot(xIntersect,yIntersect,'*')
% 
% Label the plot axes using the xlabel land ylabel commands
% xlabel('x')
% ylabel('y')
% 
% The intersection can also be found by multiplying both sided of the matrix equation by the inverse of A;
% 	A-1 * A * z = A-1 * B
% or
% 	z = A-1 * B
% 
% Use the MATLAB inv() function to solve for vector z  
% 	z = inv(A) * B
% 
% Use the MATLAB left matrix divide “\” function to solve for vector z
% 	z = A \ B
% 
% The first term of vector z is the x-intersect and the second term of vector z is the y-intersect. Make sure the result using the inverse is the same as the result using Cramer’s rule.

clear;clc; % Clean workspace and command window
fprintf("Q1:\n")

m1 = 2; % slope of line 1
b1 = 10; % y-intercept of line 1
xVector1 = [-10:10]; % x-breakpoints of line 1 (I assume this is the valid domain)
yVector1 = m1 * xVector1 + b1; % Calculate y-points of line 1 (y=mx+b) 
plot(xVector1, yVector1)

m2 = -1; % slope of line 2
b2 = 20; % y-intercept of line 2
xVector2 = [-10:10]; % domain of line 2
yVector2 = m2 * xVector2 + b2; % range of line 2
hold on; plot(xVector2, yVector2) % plot line 2 on same plot as line 1

A = [-m1, 1; -m2, 1]; 
B = [b1, b2]'; % y-intercept vector (column)

xIntersect = det([b1 1;b2 1]) / det(A); % Find x position of intersection point
yIntersect = det([-m1 b1;-m2 b2]) / det(A); % Find y position of intersection point
disp("Point is at: (" + xIntersect + ", " + yIntersect + ")")
plot(xIntersect, yIntersect,"*") % Add the intersection point to the plot
xlabel('x')
ylabel('y')

% Use different ways of finding the intersection point ...
z1 = [xIntersect, yIntersect];
z2 = inv(A) * B;
z3 = A\B;

%%% Below did not work, would have to adjust because of row vs
%%% column vector...
% if z1 == z2
%     z_12_same = 1
% else
%     z_12_same = 0
% end
% 
% if z2 == z3
%     z_23_same = 1
% else
%     z_12_same = 0
% end

disp("z1 = (" + z1(1) + "," + z1(2) + ")")
disp("z2 = (" + z2(1) + "," + z2(2) + ")")
disp("z3 = (" + z3(1) + "," + z3(2) + ")")

hold off



% Q2:
% Find equation of line through two points (2 points)
% 
% Define two points as follows:
% x1 = -5
% y1 = -5
% x2 = 5
% y2 = 10
% 
% We want an equation in the form of a line.
% 	C1 + C2 * x = y
% We are given two points for x and y so we can solve for C1 and C2 by first writing the equations in matrix form.
% A * C = B
% or
% 1x11x2 * C1C2 = y1y2
% 
% Solve for the coefficients as follows:
% 	C = A-1 * B
% 
% Use the MATLAB inv() function to solve for vector C  
% 	C = inv(A) * B
% 
% Use the MATLAB left matrix divide “\” function to solve for vector C
% 	C = A \ B
% 
% Define x-breakpoints for line as: xVector = -10:10;
% Calculate corresponding y points as: yVector = C(1) + C(2) * xVector;
% Plot the line and two original points on the same plot.
% Label the plot axes.

% Clean the workspace of variables, store given points
clear;
fprintf("\nQ2:\n")
x1 = -5;
y1 = -5;
x2 = 5;
y2 = 10;

% Formula of line: (C_1 + C_2 * x = y)
% This is asking for us to solve system of equations using matrices
% Recall that this is setup where the matrix form is:
% [coefficients_array]*[Variables_col_vector]=[constants_col_vector]
A = [1, x1; 1, x2]; % Coefficient matrix array
% C = [C1; C2]; % Variables array to solve for
B = [y1; y2]; % 'Constants' array

% Find variables matrix solution, compare different methods
C1 = A^(-1)*B;
C2 = inv(A)*B;
C3 = A\B; % The 'B' must be on 'top'
disp("C1 = (" + C1(1) + "," + C1(2) + ")")
disp("C2 = (" + C2(1) + "," + C2(2) + ")")
disp("C3 = (" + C3(1) + "," + C3(2) + ")")
C = C2;

xVector = [-10, 10]; % define domain of line 2
yVector = C(1) + C(2) * xVector; % define range of line 2
% plot line and two original given points on same plot
figure % Create a new figure window, so the old one from Q1 does not get overwritten
plot(xVector, yVector)
hold on
plot(x1,y1,"*"); plot(x2,y2,"*")
xlabel('x')
ylabel('y')
hold off



% Q3:
% Find best fit equation of line through three points (2 points)
% 
% Define three points as follows:
% x1 = -5
% y1 = -5
% x2 = 5
% y2 = 10
% x3 = 1
% y3 = 0
% 
% We want an equation in the form of a line.
% 	C1 + C2 * x = y
% We are given three points for x and y so we can solve for C1 and C2 by first writing the equations in matrix form.
% A * C = B
% or
% 1x11x21x3 * C1C2 = y1y2y3
% 
% In this case the inverse function does not work because matrix A is not square. Instead we use the pseudoinverse pinv() function.  
% Use the MATLAB pinv() function to solve for vector C  
% 	C = pinv(A) * B
% 
% Use the MATLAB left matrix divide “\” function to solve for vector C
% 	C = A \ B
% 
% In cases where there is a unique solution “pinv” and “\” provide a solution that minimizes the sum of the errors squared (method of least squares).  
% If there are not enough linear independent equations to find all the coefficients, then infinite solutions exist.  For example, if you try fitting a line to data points that all have the same x values, the slope could be any value.  In this case the equations are said to be rank deficient (the “rank” of A is less than the number of coefficients in C).  For rank deficient equations, the left matrix divide “\” function will assign zero to some of the coefficients while the pseudoinverse will find a solution that may assigns numbers to all the coefficients while minimizing the “norm” function. 
% Define x-breakpoints for line as: xVector = -10:10;
% Calculate corresponding y points as: yVector = C(1) + C(2) * xVector;
% Plot the line and three original points on the same plot.
% Label the plot axes.


clear; % Clean the workspace of variables in previous problem
fprintf("\nQ3:\n")

% Define given points
x1 = -5;
y1 = -5;
x2 = 5;
y2 = 10;
x3 = 1;
y3 = 0;

A = [1, x1; 1, x2; 1, x3]; % Define array A as 'coefficient' matrix
B = [y1; y2; y3]; % Define column vector B as 'constants' matrix
C = A\B; % Find the variables of C1 and C2
C_alt = pinv(A)*B; % Use alternative method to check above; use 'psuedoinverse' instead of inverse because matrix not square....
C1 = C(1);
C2 = C(2);
C_check = (round(C(1),8) == round(C_alt(1),8));
disp("Check if C and C_alt calculate to the same: " + C_check)

% Plot points that we need to find best fit line for
figure
hold on
plot(x1, y1, "*")
plot(x2, y2, "*")
plot(x3, y3, "*")
xlim([-10,10]) %%%
ylim([-15,15]) %%%

% Plot best fit line through the points
xVector = [-10:10]; % domain to plot
yVector = C1 + C2*xVector;
plot(xVector, yVector)
hold off





% Q4:
% Find equation of parabola through three points (2 points)
% 
% Define three points as follows:
% x1 = -5
% y1 = -5
% x2 = 5
% y2 = 10
% x3 = 1
% y3 = 0
% 
% We want an equation in the form of a parabola:
% 	C1 + C2 * x + C3 * x2 = y
% We are given three points for x and y so we can solve for C1, C2 and C3 by first writing the equations in matrix form.
% A * C = B
% or
% 1x1x121x2x221x3x32 * C1C2C3 = y1y2y3
% 
% Use the MATLAB inv() function to solve for vector 
% 	C = inv(A) * B
% 
% Use the MATLAB left matrix divide “\” function to solve for vector C
% 	C = A \ B
% 
% Define x-breakpoints for line as: xVector = -10:10;
% Calculate corresponding y points as: yVector = C(1) + C(2) * xVector + C(3) * xVector.^2;
% 	where the dot notation “.^2” means to square each term of the array
% Plot the parabola and three original points on the same plot.
% Label the plot axes.


clear; % Clean the workspace of variables in previous problem
fprintf("\nQ4:\n")

% Define given points
x1 = -5;
y1 = -5;
x2 = 5;
y2 = 10;
x3 = 1;
y3 = 0;

% Plot the given points
figure
hold on
plot(x1, y1, "*")
plot(x2, y2, "*")
plot(x3, y3, "*")
xlim([-10, 10])
ylim([-10, 15])

% Find variables (C1, C2, C3)
A = [1, x1, x1^2; % Coefficient matrix (coefficients in front of variable C1, C2 etc)
     1, x2, x2^2;
     1, x3, x3^3];
B = [y1; y2; y3]; % Constants matrix
C = inv(A)*B; % Find variable values (C1, C2, C3)
C_alt = A\B; % Find variable values using different method to check

% Check if alternative method agrees with original method
if round(C(1),8) == round(C_alt(1),8) & round(C(2),8) == round(C_alt(2),8) & round(C(3),8) == round(C_alt(3),8)
    C_check = true;
else
    C_check = false;
end
disp("Check if C and C_alt calculate to the same: " + C_check)

xVector = [-10, 10];
yVector = C(1) + C(2) * xVector + C(3) * xVector.^2;
plot(xVector, yVector)
hold on
xlabel('x')
ylabel('y')
hold off











% Q5:
% Find equation of plane through three points (2 points)
% 
% Define three points as follows:
% x1 = -10
% y1 = 10
% z1 = 0
% x2 = 10
% y2 = -10
% z2 = 0
% x3 = -10
% y3 = -10
% z3 = -30
% 
% We want an equation of the form of a plane:
% 	C1 + C2 * x + C3 * y = z
% We are given three points for x, y and z so we can solve for C1, C2 and C3 by first writing the equations in matrix form. 
% A * C = B
% or
% 1x1y11x2y21x3y3 * C1C2C3 = z1z2z3
% 
% Use the MATLAB inv() function to solve for vector C
% 	C = inv(A) * B
% 
% Use the MATLAB left matrix divide “\” function to solve for vector C
% 	C = A \ B
% 
% Define x-breakpoints for line as: xVector = -10:2:10;
% Define y-breakpoints for line as: yVector = -10:10;
% Define x and y arrays using the meshgrid() function.
% 	[xArray,yArray] = meshgrid(xVector,yVector) 
% Calculate corresponding z points as: 
% zArray = C(1) + C(2) * xArray + C(3) * yArray; 
% Plot the plane surface and three original points on the same plot.
% surf(xArray,yArray,zArray)
% hold on
% plot3(x1,y1,z1,'.r')
% plot3(x2,y2,z2,'.r')
% plot3(x3,y3,z3,'.r') 
% Label the plot axes.

fprintf("\nQ5:\n")
clear;

% Define the given points
x1 = -10;
y1 = 10;
z1 = 0;
x2 = 10;
y2 = -10;
z2 = 0;
x3 = -10;
y3 = -10;
z3 = -30;

% Find the variable matrix
A = [1, x1, y1; % Coefficient matrix
     1, x2, y2;
     1, x3, y3];
B = [z1, z2, z3]'; % Constants matrix
C = inv(A)*B; % Variable matrix
C_alt = A\B; % Alternative method to find variable matrix

% Check if alternative method agrees with original method
if round(C(1),8) == round(C_alt(1),8) & round(C(2),8) == round(C_alt(2),8) & round(C(3),8) == round(C_alt(3),8)
    C_check = true;
else
    C_check = false;
end
disp("Check if C and C_alt calculate to the same: " + C_check)

% Solve for Z
xVector = [-10:2:10]; % domain
yVector = [-10:10]; % range
[xArray, yArray] = meshgrid(xVector, yVector); % believe this just copies thex vector how ever many times the y vector has points?
zArray = C(1) + C(2)*xArray + C(3)*yArray;

% Plot plane and points
figure
surf(xArray, yArray, zArray) % Plot plane (surface)
hold on
plot3(x1, y1, z1, ".r") % Plot first point
plot3(x2, y2, z2, ".r") % Plot second point
plot3(x3, y3, z3, ".r") % Plot third point
xlabel('x')
ylabel('y')
hold off
