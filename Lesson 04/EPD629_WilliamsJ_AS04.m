% J.Williams
% University of Wisconsin-Madison
% EPD629: Powertrain Systems and Controls
% Assignement 04
% 2022-02-18

close all; clear; clc;

s = tf('s') % Define 's' for symbolic transfer function format

% Calibration Parameters
Kt = 0.063; % N*m/A
Kb = Kt; % V/(rad/s)
L = 1e-3; % Henry (Inductance)
R = 0.35; % Ohms
J = 50e-6; % kg*m^2
b = 0.001; % N*m/(rad/s)
K = 0.1; % N*m/rad
MotorTorque = 1; % Magnitude of input torque step (N*m)
ArmatureVoltage = 10; % Magnitude of input voltage step (V)

% Part 1: Throttle with damping, without spring, output angular velocity

G1 = 1/(J*s + b) % Transfer function
figure(1)
step(MotorTorque*G1)
hold on

A1 = -b/J
B1 = 1/J
C1 = 1
D1 = 0
Throttle1_ss = ss(A1, B1, C1, D1) % State space
step(MotorTorque*Throttle1_ss)

% Plot system poles using pzmap
figure(2)
pzmap(Throttle1_ss) % This should be identical to pzmap(G1)

