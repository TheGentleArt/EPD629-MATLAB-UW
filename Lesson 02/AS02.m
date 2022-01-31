% J.Williams
% University of Wisconsin-Madison
% EPD629: Powertrain Systems and Controls
% Assignement 02
% 2022-01-31

% Link to textbook end of chapter problems:
% https://bcs.wiley.com/he-bcs/Books?action=mininav&bcsId=11568&itemId=1119474221&assetId=473646&resourceId=45905&newwindow=true

% Homework assignment (copy/pasted) below:
% % This homework assignment is to create a vehicle longitudinal dynamic model for a truck.  Operation of the truck will be simulated in Simulink.  The free body diagram for the truck is shown below.   
% %  The governing equations for this model are as follows:
% % 	Faero = ½  Cd  A    v2
% % 	Frolling = Cr  m  g
% %  	Fincline = m  g  sin()
% % 	∑ F = m  a
% % 	a = (Ftrac - Faero - Frolling - Fincline) / m
% % 	v =  a dt
% % 	x =  v dt
% % In this model tractive force is the input.  Vehicle speed and position are the outputs.
% % 
% % 1)	Create an initialization script with the following commands (1 points)
% % clear
% % close all
% % clc
% %  
% % % Ambient Conditions
% %  
% % AmbientPress_kPa = 100;
% % AmbientTemp_C = 20;
% %  
% % % Constants
% %  
% % AirGasConstant_kJpkgK = 0.278;
% % CelsiusToKelvinConstant = 273.15;
% % GravitationalAccel_mps2 = 9.8067;
% %  
% % % Road Grade
% %  
% % PitchAngleData_deg = [2 2];
% % PitchAnglePositionBP_m = [0 1e6];
% %  
% % % Vehicle
% %  
% % VehicleMass_kg = 30000;
% % AeroArea_m2 = 10;
% % AeroDragCoef = 0.6;
% % RollingDragCoef = 0.005;
% % 
% % 2)	In Simulink create the following subsystem model for Road Load (2 points)
% % 
% %  
% % 
% % Configure the pitch angle table lookup as shown below.
% %  
% %  
% % 
% % 3)	In Simulink create the following subsystem model for Vehicle Dynamics (2 points)
% % 
% %  
% % 
% % Configure the two Integrator blocks as shown below.
% % 
% %   Vehicle Speed				     Vehicle Position
% %  
% % 
% % 4)	Attach the two subsystems as shown below. (1 point)
% %  
% % Configure the Step bock as shown below.
% % 
% %  
% % 
% % 
% % 5)	Set the Stop Time as 2000 seconds. Run the simulation. Verify that your results are the same as the following using the Scope. (1 point)
% %  
% % 
% % 6)	Complete chapter 1 textbook problem 13a, 13b, and 13c. (3 points)
% % 
% % 7)	Complete chapter 1 textbook problem 14a, 14b, and 14c. (3 points)

% Q1:
close all;
clear; clc;
fprintf("\nQ1\n")

% Ambient Conditions
AmbientPress_kPa = 100; % Ambient air pressure
AmbientTemp_C = 20; % Ambient air temperature
 
% Constants
AirGasConstant_kJpkgK = 0.278; % Gas constant of air
CelsiusToKelvinConstant = 273.15; 
GravitationalAccel_mps2 = 9.8067; % Acceleration due to gravity
 
% Road Grade
PitchAngleData_deg = [2 2]; 
PitchAnglePositionBP_m = [0 1e6];
 
% Vehicle
VehicleMass_kg = 30000; % Vehicle mass
AeroArea_m2 = 10; % Frontal Area of vehicle
AeroDragCoef = 0.6; % Vehicle aerodynamic drag coefficient
RollingDragCoef = 0.005; % Vehicle rolling resistance coefficient

% % Formulas
% F_sum = m*a 
% a = (F_trac - F_aero - F_rolling - F_incline)/m 
% v = integral(a dt)
% F_aero = (1/2)*C_d*A*rho*v.^2 % Aerodynamic drag resistance force
% F_rolling = C_r*m*g % Rolling resistance force
% F_incline = m*g_*sin(theta) % Force to overcome component of gravity on hill









%Q2 (6):
%From textbook:
%1-13:
% % Given the electric network shown in Figure P1.5.
% % [Review]
% % a. Write the differential equation for the network if
% % v(t) = u(t), a unit step.
% % b. Solve the differential equation for the current, i(t), if
% % there is no initial energy in the network.
% % c. Make a plot of your solution if R/L=1.
%
%Attempt:
% If V_source-V_resistor-V_inductor=0
% V_source = V_resistor + V_inductor
% V_source = i*R + L*(di/dt)
% % % dV/dt = R*(di/dt) + L*(di/dt)^2 %%% Part a
% % % dV/dt = (R+L*(di/dt))*(di/dt)
% L*(di/dt) = V_source - i*R
% di/dt = (V_source - i*R)/L
% di/(V_source - i*R) = dt/L
% Integrate
% -ln(V - i*R)/R = 1/L*t + Constant
% solving for constant:
% Constant = -ln(V - i*R)/R - 1/L*t
% if i(0)=0 IE current is initially zero
% Constant = -ln(V)/R
% then
% -ln(V - i*R)/R = (1/L)*t -ln(V)/R
% -ln(V - i*R) = (R/L)*t - ln(V)
% ln((V - i*R)/(V) = -(R/L)*t
% (V - i*R)/V = e^(-(R/L)/t)
% 1 - i*R/V = e^(-(R/L)/t)
% -(R/V)*i = e^(-(R/L)/t) - 1
% i =  (-1 + e^(-(R/L)/t))/(-R/V)
% i = (1 + e^(-(R/L)/t)))*(V/R)




%Q3 (7)
% From textbook:
%1-14:
% % Repeat Problem 13 using the network shown in
% % Figure P1.6. Assume R=1Ω, L=0.5H, and 1/LC=16.



