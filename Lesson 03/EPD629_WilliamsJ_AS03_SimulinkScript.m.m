% J.Williams
% University of Wisconsin-Madison
% EPD629: Powertrain Systems and Controls
% Assignement 03
% 2022-02-12

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lesson 02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear; clc;
format compact;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lesson 03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
ConversionFactor_mile_per_km = 0.621371;
ConversionFactor_Nm_per_ftlb = 1.3558179;
 
% Engine
% Cummins X15 Performance Series 605 hp
% https://mart.cummins.com/imagelibrary/data/assetfiles/0033023.pdf
EngineSpeedBreakpoints_rpm = [900, 1000, 1100, 1160, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000]; % (1/min [rpm])
EngineTorqueMaxData_ftlb = [1050, 1520, 1880, 2050, 2050, 2050, 2050, 2050, 1980, 1870, 1760, 1660, 1560]; % (lbf*ft)
EngineTorqueMaxData_Nm = EngineTorqueMaxData_ftlb * ConversionFactor_Nm_per_ftlb; % (N*m)

% Plot engine torque vs engine speed (english)
figure(1)
plot(EngineSpeedBreakpoints_rpm,EngineTorqueMaxData_ftlb)
xlabel('Engine Speed (rpm)')
ylabel('Engine Brake Torque (ft*lb)')
grid on

% Plot engine torque vs engine speed (metric)
figure(2)
plot(EngineSpeedBreakpoints_rpm,EngineTorqueMaxData_Nm)
xlabel('Engine Speed (rpm)')
ylabel('Engine Brake Torque (N*m)')
grid on

% Estimate engine friction
X = [1, 1000, 1000^2; 1, 1500, 1500^2; 1, 2000, 2000^2];  % speed matrix for parabolic friction curve fit
Y = [130, 175, 240]'; % nominal friction estimate
C = X \ Y;
EngineFrictionData_Nm = C(1) + C(2) * EngineSpeedBreakpoints_rpm + C(3) * EngineSpeedBreakpoints_rpm .^2; % (N*m)

% Plot engine friction vs engine speed
figure(3)
plot(EngineSpeedBreakpoints_rpm,EngineFrictionData_Nm)
xlabel('Engine Speed (rpm)')
ylabel('Engine Friction (N*m)')
grid on

EngineTorqueTimeConstant_s = 0.1; % (N*m)

% Transmission
% Eaton Endurant HD
% https://www.eatoncummins.com/us/en-us/catalog/transmissions/endurant.specifications.html
TransGearRatioData = [14.43, 11.05, 8.44, 6.46, 4.95, 3.79, 2.91, 2.23, 1.70, 1.30, 1.00, 0.77];  

% Axle 
AxleRatio = 3.25;

% Tires
TireRevsPerMile = 500; % (rev/mi)
TireRevsPerKilometer = TireRevsPerMile * ConversionFactor_mile_per_km; % (rev/km)
