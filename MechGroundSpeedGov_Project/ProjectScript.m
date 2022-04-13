% J.Williams
% University of Wisconsin-Madison
% EPD629: Powertrain Systems and Controls
% Project Draft
% 2022-04-13

% This script and the simulink file associated with it are used to model the governor system on a small car.
% This will help in showing others trends when wanting to change things in the governor system.
% It will also be helpful to approximate values needed to achieve desired results,
% hopefully reducing prototyping and testing time.

close all;
clear; clc;

% General Environment Info
g = 9.80665; % m/s^2 gravitational acceleration
rho_air = 1.224991; % kg/m^3 air density
GradeIncline_pcnt = 0 % Grade of hill in percentage

% General Vehicle Info
TireDia_in = 18; % in
CTireRollRadius = 0.944; % Rolling Radius Factor for Tire .965 later
CurbWeight_lbm = 800; % lbm
VehiclePayload_lbm = 500; % lbm
MassVehicle_kg = (CurbWeight_lbm + VehiclePayload_lbm)*0.4535924 % Total Vehicle Weight
BrakeDrag_N = 2.5; % N of brake drag resistance
Crr = 0.0138796; % Rolling resistance coefficient estimate of RXV
Cd = 0.7691142;
VehicleFrontalArea_m2 = 1.8; % m^2

% Accelerator Pedal Assembly Info
SpringRate_lbfpin = 12.6 % lbf/in was 12.6 7.1091, 11.05266 
PedalPosMax_deg = 30 % deg
CableStrokeVsPedalPos_mmpdeg = 0.995519 % mm/deg
PedalTravelIntTension_deg = 10 % Pedal position at which tension is started
                               % within the accelerator cable, if less than
                               % 10 risk opening throttle before pedal 
                               % switch closure

% Governor System Info
LengthGovArmUpper_mm = 67.5; % mm
LengthGovArmLower_mm = 90.4;
GovArmRange_deg = 14.7; % deg of travel of gov arm allowed, was 14.7, throttle limiting estimated around 13.4
% Coefficients for governor torque vs speed
%Flyweights closed
FlyCoeffAClosed = 5.31512605042017e-8;
FlyCoeffBClosed = 0.0000162815126050419;
FlyCoeffCClosed = -0.0102941176470588;
% Flyweights open
FlyCoeffAOpen = 2.27450980392157e-7;
FlyCoeffBOpen = -0.000070868347338936;
FlyCoeffCOpen = 0.00588235294117689;

% Powertrain Info
% Gearbox/Axle Info
AxleRatio = 11.47;
GearboxRatio = 1.95;
% CVT Info
CVTEngagementSpd_rpm = 1400; % rpm
CVTEff = 0.35 % CVT Efficiency
% CVT ratio as function of vehicle speed (valid for limited scenarios)
CVTRatioCoeffVsMphA = 0.000494471902965804
CVTRatioCoeffVsMphB = -0.0203171948865201
CVTRatioCoeffVsMphC = 0.301637634575451
CVTRatioCoeffVsMphD = -2.09463337350964
CVTRatioCoeffVsMphE = 7.47628950179633
% Engine Info
% Directly below info is original test data
% EngPowerTableData_hp = [0,3.64, 4.85, 5.08, 5.39;
%                      0,3.94, 5.65, 6.1, 6.44;
%                      0, 4.22, 6.13, 6.75,7.16;
%                      0, 4.38, 6.48, 7.48, 7.75;
%                      0, 4.64, 6.8, 7.76,8.33;
%                      0,4.92, 7.12, 8.47, 8.9;
%                      0, 5.04, 7.27, 8.77, 9.37;
%                      0,5.25, 7.23, 8.81, 9.52;
%                      0, 4.81, 6.6, 8.31, 9.04];
% EngPowerTableSpdBreakpoints_rpm = (3500:500:7500);
% EngPowerTableThrottleBreakpoints_prcnt = [0:.25:1];
EngPowerTableData_hp = [0, 0, 0, 0, 0;
                     0, 0.1, 0.19, 0.25, 0.28;
                     0, 0.25, 0.6, 0.72, 0.8;
                     0, 0.7, 1.1, 1.3, 1.4;
                     0, 1.15, 1.7, 2.05, 2.2;
                     0, 1.9, 2.6, 2.9, 3.1;
                     0, 2.75, 3.5, 3.95, 4.1; % last line of added
                     0, 3.64, 4.85, 5.08, 5.39;
                     0, 3.94, 5.65, 6.1, 6.44;
                     0, 4.22, 6.13, 6.75,7.16;
                     0, 4.38, 6.48, 7.48, 7.75;
                     0, 4.64, 6.8, 7.76,8.33;
                     0,4.92, 7.12, 8.47, 8.9;
                     0, 5.04, 7.27, 8.77, 9.37;
                     0,5.25, 7.23, 8.81, 9.52;
                     0, 4.81, 6.6, 8.31, 9.04];
EngPowerTableSpdBreakpoints_rpm = (0:500:7500);
EngPowerTableThrottleBreakpoints_prcnt = [0:.25:1];

% Throttle Body Info
ThrottleReturnSpringRate_Nmmpdeg = 0.786881; % N*mm/deg
ThrottleReturnSpringPreload_Nmm = 163.7711; % N*mm
ThrottleCableDia_mm = 1.190625; % mm
ThrottleLeverArmRadius_mm = 14.25; % mm
ThrottlePosMax_deg = 82; % deg
DistCableSinkThrottle_mm = 0.25; % Distance the throttle cable sinks below
                                % the throttle lever arm due to the pulley
                                % being made of an s-shape plate bolted to
                                % a flat plate



