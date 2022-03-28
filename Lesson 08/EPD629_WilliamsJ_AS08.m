% J.Williams
% University of Wisconsin-Madison
% EPD629: Powertrain Systems and Controls
% Assignement 08
% 2022-03-27
  
% No Changes since AS07

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment:
% This homework assignment is to create an idle speed control model of a 
% spark-ignition engine with electronic throttle control.
% Complete the following assignment and submit your MATLAB script 
% file (.m or .mlx file type) and the Simulink model (.slx file type).
% Shown below is the top-level view of the model. Equations for this model 
% are shown in Guzzella Appendix B.1. Calibrations are shown 
% in Guzzella Appendix B.2.3.
%  
% Part 1: Model Initialization File (1 points)
% Copy the following commands to your model initialization (calibration) file:
%
% The remainder of the assignment is done via simulink.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that this was largely a copy-paste from the homework assignment


% This script file generates calibrations for throttle area, throttle flow, 
% volumetric efficiency, engine airflow, and engine torque. 
 
close all
clear
clc
 
% Define Common Calibration Table Breakpoints
 
EngineSpeedBreakpoints_rpm = 400:200:5000;
ManifoldPressureBreakpoints_kPa = 5:5:100;
 
% Throttle Position to Area
 
ThrottlePositionBreakpoints_pct = 0:100;
 
ThrottleAngleMin_deg = 7.9;
ThrottleAngleMax_deg = 90;
ThrottleDiameter_mm = 58.7;
ThrottleLeakageArea_mm2 = 5.3;
 
ThrottleAlpha_deg = ThrottleAngleMin_deg + ...
    (ThrottleAngleMax_deg - ThrottleAngleMin_deg) * ...
    ThrottlePositionBreakpoints_pct / 100;
 
ThrottleArea_mm2 = pi * ThrottleDiameter_mm^2 / 4 * ...
    (1 - cosd(ThrottleAlpha_deg) / cosd(ThrottleAngleMin_deg)) + ...
    ThrottleLeakageArea_mm2;
 
WideOpenThrottleArea_mm2 = ThrottleArea_mm2(end);
ThrottleAreaData_pct = 100 * ThrottleArea_mm2 / WideOpenThrottleArea_mm2;
 
figure(1)
plot(ThrottlePositionBreakpoints_pct,ThrottleAreaData_pct)
xlabel('Throttle Position (%)')
ylabel('Throttle Area (%)')
axis([0 100 0 100])
grid on

% Throttle Area to Position

ThrottleAreaBreakpoints_pct  = ThrottleAreaData_pct;
ThrottlePositionData_pct = ThrottlePositionBreakpoints_pct;

% Throttle Flow Function

ThrottlePressureRatioBreakpoints_frac = 0:0.01:1;

DischargeCoefficient_frac = 0.7;
AirGasConstant_kJpkgK = 0.287;  % kJ/kg*K
AmbientPressure_kPa = 98;
InletAirTemperature_K = 298;
RatioOfSpecificHeats_frac = 1.35;  % Cp/Cv

PR_min = (2/(RatioOfSpecificHeats_frac+1))^ ...
    (RatioOfSpecificHeats_frac/(RatioOfSpecificHeats_frac-1));

PR = max(PR_min,ThrottlePressureRatioBreakpoints_frac);
ThrottleFlowFunctionData_frac = PR.^(1/RatioOfSpecificHeats_frac) .* ...
    sqrt( 2 * RatioOfSpecificHeats_frac / (RatioOfSpecificHeats_frac - 1) .* ...
    (1 - PR.^((RatioOfSpecificHeats_frac - 1) / RatioOfSpecificHeats_frac)) );

figure(2)
plot(ThrottlePressureRatioBreakpoints_frac,ThrottleFlowFunctionData_frac)
xlabel('Pressure Ratio (%/100)')
ylabel('Throttle Flow Function')
axis([0 1 0 1])
grid on

% Intake Manifold
 
IntakeManifoldVolume_L = 5.8;
ManifoldAirTemperature_K = 340;
 
% Engine Volumetric Efficiency
 
Gamma_0 = 0.45;     % -
Gamma_1 = 3.42e-3;  % s
Gamma_2 = -7.7e-6;  % s^2
EngineDisplacement_L = 2.77;
ClearanceVolumeTotal_L = 0.277; %%% Assummed a Comp Ratio
ExhaustManifoldPressure_kPa = 108;
AirToFuelEquivalenceRatio_frac = 1;  % A/F Equivalence Ratio  =  AFR/StoichAFR StoichiometricAFR_frac = 14.7;  % Stoich A/F ratio, Sigma_o in Guzzella
StoichiometricAFR_frac = 14.7;  % Stoich A/F ratio, Sigma_o in Guzzella

[EngineSpeedArray,ManifoldPressureArray] = meshgrid(EngineSpeedBreakpoints_rpm, ...
    ManifoldPressureBreakpoints_kPa);
 
for r = 1:length(EngineSpeedArray(:,1))
    for c = 1:length(EngineSpeedArray(1,:))
        Omega_e = EngineSpeedArray(r,c) * pi / 30;
        VolumetricEfficiencyData_frac(r,c) =  ...
            (Gamma_0 + Gamma_1*Omega_e + Gamma_2*Omega_e^2) * ...
            ((ClearanceVolumeTotal_L + EngineDisplacement_L) / EngineDisplacement_L - ...
            ClearanceVolumeTotal_L / EngineDisplacement_L * ...
            (ExhaustManifoldPressure_kPa / ManifoldPressureArray(r,c))^ ...
            (1/RatioOfSpecificHeats_frac));
        VolumetricEfficiencyData_frac(r,c) = max(VolumetricEfficiencyData_frac(r,c),0);
        EngineMassFlow_gps(r,c) = ManifoldPressureArray(r,c) / ...
          (AirGasConstant_kJpkgK * ManifoldAirTemperature_K) * ...
          VolumetricEfficiencyData_frac(r,c) * EngineDisplacement_L * Omega_e / ...
          (4*pi)/(1+1/StoichiometricAFR_frac);
   end
end

figure(3)
surf(EngineSpeedArray,ManifoldPressureArray,VolumetricEfficiencyData_frac)
colorbar
xlabel('Engine Speed (rpm)')
ylabel('Intake Manifold Pressure (kPa)')
zlabel('Volumetric Efficiency (%/100)')
axis([0 5000 0 100 0 1])

figure(4)
contourf(EngineSpeedArray,ManifoldPressureArray,VolumetricEfficiencyData_frac)
colorbar
xlabel('Engine Speed (rpm)')
ylabel('Intake Manifold Pressure (kPa)')
title('Volumetric Efficiency (%/100)')
axis([0 5000 0 100])
grid on

figure(5)
surf(EngineSpeedArray,ManifoldPressureArray,EngineMassFlow_gps)
colorbar
xlabel('Engine Speed (rpm)')
ylabel('Intake Manifold Pressure (kPa)')
zlabel('Engine Mass Flow (g/s)')
axis([0 5000 0 100 0 60])

figure(6)
contourf(EngineSpeedArray,ManifoldPressureArray,EngineMassFlow_gps)
colorbar
xlabel('Engine Speed (rpm)')
ylabel('Intake Manifold Pressure (kPa)')
title('Engine Mass Flow (g/s)')
axis([0 5000 0 100])
grid on

% Engine Torque Model

FuelLowerHeatingValue_kJpkg = 44000;  % Gasoline Lower Heating Value
CrankAngleIVCToTDC_deg = 170;  % Delay A/F entering cylinder to top-center

% Thermocynamic Efficiency  =  Eta_0 + Eta_1*Omega_e + Eta_2*Omega_e^2
% where Omega_e is the engine speed in rad/sec 
% see Fig 2.32 in Guzzella

Eta_0 = 3.7313e-001;  % -                     
Eta_1 = 3.2500e-004;  % s         
Eta_2 = -5.6250e-007; % s^2       
Omega_e = EngineSpeedBreakpoints_rpm * pi / 30; 
ThermodynamicEfficiencyData_frac = Eta_0 + Eta_1*Omega_e + Eta_2*Omega_e.^2;

figure(7), plot(EngineSpeedBreakpoints_rpm,ThermodynamicEfficiencyData_frac)
axis([0 5000 0.35 0.45])
grid on
xlabel('Engine Speed (rpm)')
ylabel('Thermodynamic Efficiency (%/100)')

% Engine Friction and Accessory Torque  =  Beta_0 + Beta_2*Omega_e^2
 
Beta_0 = 15.6;        % Nm
Beta_2 = 0.175e-3;    % Nm*s^2
EngineMotoringFrictionData_Nm  =  Beta_0 + Beta_2*Omega_e.^2;
 
% Engine Inertia
 
EngineInertia_kgm2 = 0.2;  % kg*m^2




% Model Correction
 
ThrottleMassFlowCorrIntGain = 10;
VECorrectionFactorIntGain = 0.02;
 
% Spark Advance
 
NumberOfCylinders = 6;
 
CylinderAirMass_mg = 120000 * EngineMassFlow_gps / NumberOfCylinders ./ EngineSpeedArray;
 
figure(8)
surf(EngineSpeedArray,ManifoldPressureArray,CylinderAirMass_mg)
colorbar
xlabel('Engine Speed (rpm)')
ylabel('Intake Manifold Pressure (kPa)')
zlabel('Cylinder Air Mass (mg)')
axis([0 5000 0 100 0 400])
 
figure(9)
contourf(EngineSpeedArray,ManifoldPressureArray,CylinderAirMass_mg)
colorbar
xlabel('Engine Speed (rpm)')
ylabel('Intake Manifold Pressure (kPa)')
title('Cylinder Air Mass (mg)')
axis([0 5000 0 100])
grid on
 
CylinderAirMassBreakpoints_mg = 0:25:400;
 
Ne = [600 2000 2000 2000 3400]';
CylAir = [200 50 200 350 200]';
SpkAdv = [25 31 27 21 29]';
R = [ones(length(Ne),1) Ne Ne.^2 CylAir CylAir.^2];
C = R\SpkAdv;
[X,Y] = meshgrid(EngineSpeedBreakpoints_rpm,CylinderAirMassBreakpoints_mg);
for r = 1:length(X(:,1))
    for c = 1:length(Y(1,:))
        MBTSparkAdvanceData_deg(r,c) = [1 X(r,c) X(r,c)^2 Y(r,c) Y(r,c)^2] * C;
    end
end
 
figure(10)
contourf(X,Y,MBTSparkAdvanceData_deg)
colorbar
xlabel('Engine Speed (rpm)')
ylabel('Cylinder Air Mass (mg)')
title('MBT Spark Advance (deg)')
axis([0 5000 0 400])
grid on
 
Ne = [600 2000 2000 2000 3400]';
CylAir = [200 50 200 350 200]';
SpkAdv = [24 40 27 12 29]';
R = [ones(length(Ne),1) Ne Ne.^2 CylAir CylAir.^2];
C = R\SpkAdv;
[X,Y] = meshgrid(EngineSpeedBreakpoints_rpm,CylinderAirMassBreakpoints_mg);
for r = 1:length(X(:,1))
    for c = 1:length(Y(1,:))
        BorderlineSparkAdvData_deg(r,c) = [1 X(r,c) X(r,c)^2 Y(r,c) Y(r,c)^2] * C;
    end
end
 
figure(11)
contourf(X,Y,BorderlineSparkAdvData_deg)
colorbar
xlabel('Engine Speed (rpm)')
ylabel('Cylinder Air Mass (mg)')
title('Borderline Knock Spark Advance (deg)')
axis([0 5000 0 400])
grid on
 
BaseSparkAdvanceData_deg = min(MBTSparkAdvanceData_deg,BorderlineSparkAdvData_deg);
 
figure(12)
surf(X,Y,BaseSparkAdvanceData_deg)
colorbar
xlabel('Engine Speed (rpm)')
ylabel('Cylinder Air Mass (mg)')
zlabel('Base Spark Advance (deg)')
axis([0 5000 0 400])
 
figure(13)
contourf(X,Y,BaseSparkAdvanceData_deg)
colorbar
xlabel('Engine Speed (rpm)')
ylabel('Cylinder Air Mass (mg)')
title('Base Spark Advance (deg)')
axis([0 5000 0 400])
grid on
 
SparkTorqFactorRetFromMBTBP_deg = 0:5:45;
SparkTorqFactorData_frac = 1 - (SparkTorqFactorRetFromMBTBP_deg.^2) / ...
    (SparkTorqFactorRetFromMBTBP_deg(end)^2); 
 
n = length(SparkTorqFactorRetFromMBTBP_deg);
for i = 1:n
    SparkTorqFactorRetFromMBTData_deg(i) = SparkTorqFactorRetFromMBTBP_deg(n+1-i);
    SparkTorqFactorBreakpoints_frac(i) = SparkTorqFactorData_frac(n+1-i); 
end
 
figure(14)
plot(SparkTorqFactorRetFromMBTBP_deg,SparkTorqFactorData_frac)
xlabel('Spark Retard From MBT (deg)')
ylabel('Spark Retard From MBT Torque Factor (%/100)')
axis([0 45 0 1])
grid on
 
figure(15)
plot(SparkTorqFactorBreakpoints_frac,SparkTorqFactorRetFromMBTData_deg)
xlabel('Spark Retard From MBT Torque Factor (%/100)')
ylabel('Spark Retard From MBT (deg)')
axis([0 1 0 45])
grid on
 
% Fuel Injection and Oxygen Sensing Model
 
FuelInjImpactFactor_frac = 0.3;
FuelInjEvapTimeConstant_s = 0.8;
 
DelayFuelInjectionToIVC_rev = 0.75;
DelayIVCToOxygenSensor_rev = 7;
 
OxygenSensorLambdaBP_frac = [0 0.9 0.92 0.94 0.96 0.985 0.999 1.001 1.015 1.04 1.06 1.08 1.1 2];
OxygenSensorVoltageData_mV = [800 800 794 788 775 750 700 200 150 125 112 106 100 100];
 
figure(16)
plot(OxygenSensorLambdaBP_frac,OxygenSensorVoltageData_mV)
xlabel('Air to Fuel Equivalence Ratio, Lambda')
ylabel('Oxygen Sensor Voltage (mV)')
axis([0 2 0 1000])
grid on
 
% Fuel Injection Control
 
DelayFuelInjectionToIVCEst_rev = 0.75;
 
FuelInjImpactFactorEst_frac = 0.3;
FuelInjEvapTimeConstantEst_s = 0.8;
 
OxygenSensorSwitchSetpoint_mV = 450;
 
FuelControlProportionalGain = 5e-6;
FuelControlIntegralGain = 2e-5;
