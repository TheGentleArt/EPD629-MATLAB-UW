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

%Attempt:
% If V_source-V_resistor-V_inductor=0
% V_source = V_resistor + V_inductor
% V_source = i*R + L*(di/dt)
% % % dV/dt = R*(di/dt) + L*(di/dt)^2 %%% Part a
% % % dV/dt = (R+L*(di/dt))*(di/dt)
% L*(di/dt) = V_source - i*R
% di/dt = (V_source - i*R)/L
% di/dt + (R/L)*i = V_source/L
% Linear 1st order diff eq
% Form is dy/dx + P(x)*y = Q(x)
% Rho(x) = e^(int[P(x),x])
% Here, P(x) = (R/L), so then Rho(x) = exp(int[(R/L),t])
% Rho(x) = exp((R/L)*t)
% This is our missing piece
% exp((R/L)*t)*(di/dt + (R/L)*i = V/L)
% exp((R/L)*t)*(di/dt) + (R/L)*exp((R/L)*t)*i = (V/L)*exp((R/L)*t)
% The left side can be reduced to a derivative then
% d/dt[(exp((R/L)*t)*i] = (V/L)*exp((R/L)*t)
% Then integrate with respect to t,
% (exp((R/L)*t)*i = (V/L)*(L/R)*exp((R/L)*t) + C
% i*exp((R/L)*t) = (V/R)*exp((R/L)*t) + C
% i(t) = (V/R) + C*exp(-(R/L)*t)
% Initial Conditions are i(0) = 0
% i(0) = (V/R) + C*exp(0)
% C = -V/R
% Then the particular solution:
% i = (V/R) + (-V/R)*exp(-(R/L)*t)

% Plot a solution
R = 110; % Ohms
V = 14.5; % Volts
Ratio_RL = 1; % R/L
L = R/Ratio_RL; % Henries (Inductance)
t = linspace(0,20,1000) % Time vector
i = (V/R) + (-V/R)*exp(-(R/L).*t)
figure
plot(t, i) % Plot current vs time





%Q3 (7)
% From textbook:
%1-14:
% % Repeat Problem 13 using the network shown in
% % Figure P1.6. Assume R=1Ω, L=0.5H, and 1/LC=16.

% Attempt:
% If V_source - V_resistor - V_inductor - V_capacitor = 0
% q is charge ; i = dq/dt
% V_source = V_resistor + V_inductor + V_capacitor
% V(t) = i*R + L(di/dt) + q/C
% q/C can be expressed as (1/C)*int[i dtau,0,t]
% V(t) = i*R + L(di/dt) + (1/C)*int[i dtau,0,t]
% V(t) = i*(R/L) + (di/dt) + (1/(L*C))*int[i dtau,0,t]
% at time t=0, say V(0) = 0 (Homogenous portion)
% 0 = i*(R/L) + (di/dt) + (1/(L*C))*int[i dtau,0,t]
% 0 = (di/dt) + i*(R/L) + (1/(L*C))*int[i dtau,0,t]
% Differentiate with resepect to t,
% 0 = (d^2i/dt^2) + (R/L)*(di/dt) + (1/(L*C))*i
% The 'auxillary equation' is then: m^2 + (R/L)*m - (!/(L*C))=0
% aux_coefficients = [1, R/L, 1/(L*C)]
aux_coefficients = [1, R/L, 16];
roots(aux_coefficients);
% This is 'Scenario 3', where roots are complex
% General solution for Scenario 3 is:
% y = e^(alpha*x)*(A*cos(beta*x)+B*sin(beta*x))




% Previous work discarded below:
% V(t) = i*R + L(di//dt) + q/C
%%
% V(t) = L(di//dt) + i*R + q/C
% dV(t)/dt = L(d^2i/dt^2) + R*(di/dt) + (1/C)*dq/dt
% dV(t)/dt = L(d^2i/dt^2) + R*(di/dt) + (1/C)*i
% (1/L)*(dV(t)/dt) = d^2i/dt^2 + (R/L)*(di/dt) + (1/(L*C))*i
% Must find natural response and forced response to get total response
% For natural response, say V(t) = 0
% Initial conditions i(o) = 0; di/dt(0) = 0
% Assume i(t) = A*e^(s*t)     IE i(t) = A*exp(s*t)
%%
% Characteristic equation
% s^2 + (R/L)*s + 1/(L*C) = 0
% Roots are:
% s = (-1/2)*(R/L) +/- ((R/(2*L))^2-1/(L*C))^(1/2)
% Recall i(t) = A*exp(s*t)
% 
% Since i = dq/dt, then
% V(t) = (dq/dt)*R + L*(d^2q/dt^2) + q/C
% V(t) = L*(d^2q/dt^2) + (dq/dt)*R + q/C
% V(t)/L = (d^2q/dt^2) + (R/L)*(dq/dt) + (1/(C*L))*q