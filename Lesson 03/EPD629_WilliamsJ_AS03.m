% J.Williams
% University of Wisconsin-Madison
% EPD629: Powertrain Systems and Controls
% Assignement 03
% 2022-02-12 (updated 2022-03-05 after noticing missing information)
  
% Part 4 is 'final product' of the assingment, may want to skip to part 4
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parts 1-2 of the assignment were turned in via .pdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 3: Flywheel without damping

% Recall that
% Torque = (Moment of Inertia)*(angular acceleration)
% T(t) = J*d/dt[omega]
% T(t) = J*alpha(t)

format compact
clear;clc;
close all

J = 10; % Moment of inertia
Trq = 100; % Torque

%%% Initial speed is zero

G = tf(1,[J, 0]) % Flywheel transfer function. [Control System Toolbox function]...
%% G(s) is the laplace transfer function
% G(s) = input/output = omega(s)/T(s) = 1/(J*s)
% In MATLAB, tf() is a tranfer function. tf(numerator/denominator)
% Here, the [J, 0] represents the denominator. The '0' is just a
% placeholder, the denominator is actually J s ... if damping was considered, 0
% should be replaced by the damping
% placeholder, it will
figure(1)
step(Trq*G) % [Control System Toolbox function] The step function is for a unit step, so need to multiply G by Trq


syms s % Defines 's' as a symbol for use with Symbolic Math. 

T = Trq/s % % Input function. Laplace transform for torque step (Function for input torque. Step function with magnitude of 'Trq') Recall 1/s represents a step function (u(t))
G = 1/(J*s) % % Transfer function. Flywheel symbolic transfer function

omega_s = T*G;
omega_t = ilaplace(omega_s) % Inverse laplace transform [Symbolic Math Toolbox function]

% Move on to another case...
%%% Initial speed is 50 rad/s

omega_0 = 50; % initial angular speed (rad/s)
omega_s = Trq/(J*s^2) + omega_0/s % Same as T*G + omega_0/s
omega_t = ilaplace(omega_s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 4: Flywheel with damping and initial conditions
close all; % Will ignore 
clear; clc;
fprintf('Previous part of assignment cleared... \n\n')

% Defined Inputs
J = 10; % Moment of inertia
Trq = 100; % Torque (N*m)
b = 1; % Damping
omega_0 = 50; % rad/s
fprintf("Inputs were:\n" + ...
    " J = %f\n" + ...
    " Trq = %f\n" + ...
    " b = %f\n" + ...
    " omega_0 = %f\n\n", ...
    J,Trq,b,omega_0)

% G = tf(1,[J, 0]) % Flywheel transfer function
% figure(1) % open a new figure window, figure 1
% step(Trq*G) % plots the step response, omega(t) 

syms s t % define symbolic variable

% Governing equation
% T - b*omega(t) = J * d/dt[omega(t)]
% Torque input - frictional losses (same as damping * shaft velocity) = moment of inertia * shaft acceleration

T_s = Trq/s % define T(s) function
fprintf('\n')
G_s = 1/(J*s+b)+(J*omega_0)/((J*s+b)*(T_s)) % define G(s) function - takes into account damping, and initial conditions
fprintf('\n')

omega_s = T_s*G_s % throttle shaft speed in terms of frequency
fprintf('\n')
omega_t = ilaplace(omega_s) % throttle shaft speed in terms of time
fplot(omega_t)

