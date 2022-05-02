% J.Williams
% University of Wisconsin-Madison
% EPD629: Powertrain Systems and Controls
% Project Draft
% 2022-04-29

close all;
clear; clc;

% General Environment Info
g = 9.80665; % m/s^2 gravitational acceleration
rho_air = 1.224991; % kg/m^3 air density
GradeIncline_pcnt = 0; % Grade of hill in percentage

% General Vehicle Info
TireDia_in = 18; % in
CTireRollRadius = 0.944; % Rolling Radius Factor for Tire .965 later
CurbWeight_lbm = 800; % lbm
VehiclePayload_lbm = 360; % lbm
MassVehicle_kg = (CurbWeight_lbm + VehiclePayload_lbm)*0.4535924; % Total Vehicle Weight
BrakeDrag_N = 2.5; % N of brake drag resistance
Crr = 0.0138796; % Rolling resistance coefficient estimate of RXV
Cd = 0.7691142;
VehicleFrontalArea_m2 = 1.65; % m^2
DesCruisingSpeed_mph = 12.5; % mph % Desired flat-ground steady-state veh spd
MnfgSpdAdjGain = .925; % Amount of adjustment to make on the end-of-line 
                      % governor dyno stand for each iteration of
                      % adjustment.
MnfgSpdAdjItr = 3; % Amount of adjustments to make on the end-of-line
                  % governor dyno stand to adjust vehicle speed to desired
                  % flat-ground steady-state cruising speed.

% Accelerator Pedal Assembly Info
SpringRate_lbfpin = 11.5; % lbf/in was 12.6 7.1091, 11.05266 
PedalPosMax_deg = 30; % deg
CableStrokeVsPedalPos_mmpdeg = 0.995519; % mm/deg
Pedal_Change_After_Dyno = 0; % Initial Value set, simulink should overwrite
PedalTravelIntTension_deg = 10; % Pedal position at which tension is started
                               % within the accelerator cable, if less than
                               % 10 risk opening throttle before pedal 
                               % switch closure.

% Governor System Info
LengthGovArmUpper_mm = 67.5; % mm
LengthGovArmLower_mm = 90.4; % mm
GovArmRange_deg = 14.7; % deg of travel of gov arm allowed, was 14.7, throttle limiting estimated around 13.4
% Coefficients for governor torque vs speed (from test data)
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
CVTEngagementSpd_rpm = 1500; % rpm
CVTEff = 0.375; % CVT Efficiency
% CVT ratio as function of vehicle speed (valid for limited scenarios)
CVTRatioCoeffVsMphA = 0.000494471902965804;
CVTRatioCoeffVsMphB = -0.0203171948865201;
CVTRatioCoeffVsMphC = 0.301637634575451;
CVTRatioCoeffVsMphD = -2.09463337350964;
CVTRatioCoeffVsMphE = 7.47628950179633;
% CVT Efficiency Lookup Table Info

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

% Below replicates the manufacturing line adjusting the position of the 
% accelerator cable (thus increasing/decreasing tension/slack in the
% cable), to get the desired cruising speed with the given accelerator
% cable spring setup. This will consist of estimating the amount of slack
% to introduce or remove, multiplied by a gain, and repeating.
% This should result in the model being set near the desired cruising
% speed.
GovernorSet = 0; % Indicates that the governor has not been set yet.
StopTime = 35; % Sets the amount of time to run the simulation
CaptureTime = StopTime - 5;
DoSim = sim('EPD629_WilliamsJ_Project_.slx');
MnfgLineEOLTestSpeedList = DoSim.SimResults.VehicleSpeed_mph.data;
for n = 1:MnfgSpdAdjItr
    Pedal_Change_After_Dyno = Pedal_Change_After_Dyno + DoSim.SimResults.Pedal_Change_Itr.data*MnfgSpdAdjGain;
    DoSim.SimResults.Pedal_Change_Itr.data*MnfgSpdAdjGain;
    DoSim = sim('EPD629_WilliamsJ_Project_.slx');
    MnfgLineEOLTestSpeedList = [MnfgLineEOLTestSpeedList, DoSim.SimResults.VehicleSpeed_mph.data];
end
GovernorSet = 1; % Supposed to stop the governor from being adjusted after this but...

% Show results after getting off of governor stand.
disp(' ')
disp('Once getting off of the manufacturing line governor stand...')
VarToDisp = ['   Cable tension starts at: ~',num2str(DoSim.SimResults.PedalPosWhenTensionStarts_deg.data),' deg of pedal travel.'];
disp(VarToDisp)
VarToDisp = ['   Cruising vehicle speed set to: ~', num2str(DoSim.SimResults.VehicleSpeed_mph.data),' mph.'];
disp(VarToDisp)
VarToDisp = ['   Engine Speed while cruising is: ~', num2str(DoSim.SimResults.EngineSpeed_rpm.data), ' rpm.'];
disp(VarToDisp)
VarToDisp = ['   Throttle Position while cruising is: ~', num2str(DoSim.SimResults.ThrottlePos_deg.data),' deg.'];
disp(VarToDisp)
disp(' ')
VarToDisp = ['Vehicle Speeds for manufacturing line dyno were: ',num2str(MnfgLineEOLTestSpeedList)];
disp(VarToDisp)
disp(' ')
disp('Now moving onto simulating testing scenarios...')
disp(' ')

% Now the model SHOULD be close to how the vehicle should be off the line
% Can now do whatever other simulations needed...
disp('Simulations are done, homie.')

%
%StopTime = 90; % Set the amount of simulation time
DoSim = sim('EPD629_WilliamsJ_Project_.slx');
disp(' Sim ran again')
