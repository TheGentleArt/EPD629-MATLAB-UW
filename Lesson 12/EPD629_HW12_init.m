
clear
close all
clc

% Turbocharger

load compressor_large_cmp.mat
MXRVOL = 0.8*MXRVOL;
suffix = 'cmp';
eval(['APPR_',suffix,' = APPR;'])
eval(['ARPM_',suffix,' = ARPM;'])
eval(['DFAC_',suffix,' = DFAC;'])
% eval(['DIAW_',suffix,' = DIAW;'])
eval(['EFFI_',suffix,' = EFFI;'])  
eval(['MASS_',suffix,' = MASS;'])  
eval(['MXRSP_',suffix,' = MXRSP;'])  
eval(['MXRVOL_',suffix,' = MXRVOL;']) 
% eval(['NPR_',suffix,' = NPR;']) 
% eval(['NRPM_',suffix,' = NRPM;']) 
eval(['PMAX_',suffix,' = PMAX;']) 
eval(['PMIN_',suffix,' = PMIN;'])     
clear APPR ARPM DFAC DIAW EFFI MASS MXRSP MXRVOL NPR NRPM PMAX PMIN suffix

load turbine_med_trb.mat
MXRVOL = 1.2*MXRVOL;
suffix = 'trb';
eval(['APPR_',suffix,' = APPR;'])
eval(['ARPM_',suffix,' = ARPM;'])
% eval(['DIAW_',suffix,' = DIAW;'])
eval(['EFFI_',suffix,' = EFFI;'])  
eval(['MASS_',suffix,' = MASS;'])  
eval(['MXRSP_',suffix,' = MXRSP;'])  
eval(['MXRVOL_',suffix,' = MXRVOL;']) 
% eval(['NPR_',suffix,' = NPR;']) 
% eval(['NRPM_',suffix,' = NRPM;']) 
% eval(['NUMRACK_',suffix,' = NUMRACK;']) 
eval(['PMAX_',suffix,' = PMAX;']) 
eval(['PMIN_',suffix,' = PMIN;'])     
clear APPR ARPM DIAW EFFI MASS MXRSP MXRVOL NPR NRPM NUMRACK PMAX PMIN suffix

TurbochargerInertia_kgm2 = 1e-4;

% Gas Properties

Cp_cmp = 1.004;   % kJ/(kg*K)  @ 300 K
Cv_cmp = 0.717;   % kJ/(kg*K)  @ 300 K
R_cmp = Cp_cmp - Cv_cmp;   % kJ/(kg*K)
k_cmp = Cp_cmp/Cv_cmp;

Cp_trb = 1.073;   % kJ/(kg*K)  @ 700 K
Cv_trb = 0.786;   % kJ/(kg*K)  @ 700 K
R_trb = Cp_trb - Cv_trb;   % kJ/(kg*K)
k_trb = Cp_trb/Cv_trb;

Cp_im = 1.013;   % kJ/(kg*K)  @ 400 K
Cv_im = 0.725;   % kJ/(kg*K)  @ 400 K
IntakeManifoldGasConstant_kJpkgK = Cp_im - Cv_im;   % kJ/(kg*K)
IntakeManifoldRatioOfSpecificHeats_frac = Cp_im/Cv_im;

ExhaustManifoldGasConstant_kJpkgK = R_trb;   % kJ/(kg*K)

AirMassFracO2 = 0.233;

% EGR Valve Position to Area

EGRValvePositionBreakpoints_pct = 0:100;

EGRValveAngleMin_deg = 7.9;
EGRValveAngleMax_deg = 90;
EGRValveDiameter_mm = 58.7;
EGRValveLeakageArea_mm2 = 5.3;

EGRValveAlpha_deg = EGRValveAngleMin_deg + ...
    (EGRValveAngleMax_deg - EGRValveAngleMin_deg) * ...
    EGRValvePositionBreakpoints_pct / 100;

EGRValveArea_mm2 = pi * EGRValveDiameter_mm^2 / 4 * ...
    (1 - cosd(EGRValveAlpha_deg) / cosd(EGRValveAngleMin_deg)) + ...
    EGRValveLeakageArea_mm2;

EGRValveWideOpenArea_mm2 = EGRValveArea_mm2(end);
EGRValveAreaData_pct = 100 * EGRValveArea_mm2 / EGRValveWideOpenArea_mm2;

% figure(1)
% plot(EGRValvePositionBreakpoints_pct,EGRValveAreaData_pct)
% xlabel('EGR Valve Position (%)')
% ylabel('EGR Valve Area (%)')
% axis([0 100 0 100])
% grid on

% EGR Valve Area to Position

EGRValveAreaBreakpoints_pct = EGRValveAreaData_pct;
EGRValvePositionData_pct = EGRValvePositionBreakpoints_pct;

% Back Pressure Valve Position to Area

BPVPositionBreakpoints_pct = EGRValvePositionBreakpoints_pct;
BPVAreaData_pct = EGRValveAreaData_pct;
BPVDiameter_mm = 100;
BPVWideOpenArea_mm2 = pi * BPVDiameter_mm^2 / 4;

% EGR Valve Area to Position

BPVAreaBreakpoints_pct = BPVAreaData_pct;
BPVPositionData_pct = BPVPositionBreakpoints_pct;

% Define Common Calibration Table Breakpoints

EngineSpeedBreakpoints_rpm = 600:200:3000;
ManifoldPressureBreakpoints_kPa = 80:20:500;
DesiredFuelBreakpoints_mg = 0:20:500;

% Manifold Volume

CACVolume_L = 40;
EGRCoolerToValveVolume_L = 2;
ExhaustManifoldVolume_L = 8;
IntakeManifoldVolume_L = 8;
TurbineToBPVVolume_L = 10;

% Engine Volumetric Efficiency

EngineDisplacement_L = 8;
ClearanceVolumeTotal_L = 8/15;
StoichiometricAFR_frac = 14.7;  % Stoich A/F ratio, Sigma_o in Guzzella

[EngineSpeedArray,ManifoldPressureArray] = meshgrid(EngineSpeedBreakpoints_rpm, ...
    ManifoldPressureBreakpoints_kPa);

for r = 1:length(EngineSpeedArray(:,1))
    for c = 1:length(EngineSpeedArray(1,:))
        VolumetricEfficiencyData_frac(r,c) =  0.9;
   end
end

NumberOfCylinders = 6;

% Exhaust Temperature and Residual Mass

[EngineSpeedArray,DesiredFuelArray] = meshgrid(EngineSpeedBreakpoints_rpm, ...
    DesiredFuelBreakpoints_mg);

for r = 1:length(EngineSpeedArray(:,1))
    for c = 1:length(EngineSpeedArray(1,:))
        ExhaustTemperatureData_degC(r,c) = min(800, 400 + 2*DesiredFuelArray(r,c));
        ResidualCorrectionFactorData(r,c) = 1;
   end
end

% figure(2)
% surf(EngineSpeedArray,DesiredFuelArray,ExhaustTemperatureData_degC)
% colorbar
% xlabel('Engine Speed (rpm)')
% ylabel('Desired Fuel (mg)')
% zlabel('Exhaust Temperature (deg C)')
% xlim([600 3000])
% ylim([0 500])

% Flow Function

DischargeCoefficient_frac = 0.7;
PressureRatioBreakpoints_frac = 0:0.01:1;
PR_min = (2 / (k_trb + 1))^ ...
    (k_trb / (k_trb - 1));
PR = max(PR_min,PressureRatioBreakpoints_frac);
FlowFunctionData_frac = PR.^(1/k_trb) .* ...
    sqrt( 2 * k_trb / (k_trb - 1) .* ...
    (1 - PR.^((k_trb - 1) / k_trb)) );

% figure(3)
% plot(PressureRatioBreakpoints_frac,FlowFunctionData_frac)
% xlabel('Pressure Ratio (%/100)')
% ylabel('Throttle Flow Function')
% xlim([0 1])
% ylim([0 1])
% grid on

% Wastegate

WastegateWideOpenArea_mm2 = pi/4*20^2;
EGRValveWideOpenArea_mm2 = pi/4*50^2;

WastegateAreaData_pct = [0 100];
WastegateAreaGagePressBreakpoints_kPa = [250 300];

% Charge Air Cooler

CACFlowRestrictionArea_mm2 = 2000;
CACDischargeCoefficient_frac = 0.7;
CACEffectivenessData_frac = [1 0.95 0.92 0.885 0.87 0.85];
CACEffectivenessMassFlowBP_gps = [0 50 100 200 300 500];

% EGR Cooler

EGRCFlowRestrictionArea_mm2 = 500;
EGRCDischargeCoefficient_frac = 0.7;
EGRCEffectivenessData_frac = [1 0.95 0.92 0.885 0.87 0.85];
EGRCEffectivenessMassFlowBP_gps = [0 50 100 200 300 500];
