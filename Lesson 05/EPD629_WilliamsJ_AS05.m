% J.Williams
% University of Wisconsin-Madison
% EPD629: Powertrain Systems and Controls
% Assignement 05
% 2022-02-26

close all;
clear; clc;
 
s = tf('s')  % Define s for symbolic transfer function format
 
% Calibration parameters
Kt = 0.063;  % Nm/A
Kb = Kt;  % V/(rad/s)
L = 1e-3;  % H
R = 0.35;  % Ohm
J = 50e-6;  % kg*m^2
b = 0.001;  % Nm/(rad/s)
K = 0.1;  % Nm/rad
ArmatureVoltage = 10;  % Magnitude of input voltage step (V)
Kprop = 3;  % Controller proportional gain
Kint = 1.5;  % Controller integral gain

G1 = Kt/R % Simplified from Kt/(L*s+R) since assume happens quickly so can ignore inductor
G2 = 1/(J*s)
G2 = feedback(G2, b)
G2 = feedback(G2/s, K)
Gp = feedback(G1*G2, Kb*s) % Plant transfer function
figure(1)
step(ArmatureVoltage*Gp)
figure(2)
pzmap(Gp)


% Part 2: 
% Error for Step Input with proportional control
Gc = Kprop % Controller transfer function
Gre = feedback(1, Gc*Gp) % Error transfer function for reference input
figure(3)
step(Gre)
figure(4)
pzmap(Gre)
erss = dcgain(Gre) % Steady-state error for unit step reference input

% Part 3:
% Error for Torque Step Disturbance with Proportional Control
Gc = Kprop % Controller transfer function
Gde = -R/Kt*feedback(Gp,Gc) % Error transfer function for disturbance
figure(5)
step(Gde)
figure(6)
pzmap(Gde)
edss = dcgain(Gde) % Steady state error for unit step disturbance

% Part 4:
% Error for step input with integral control
Gc = Kint/s % Controller transfer function
Gre = feedback(1, Gc*Gp) % Error transfer function for reference input
figure(7)
step(Gre)
figure(8)
pzmap(Gre)
erss = dcgain(Gre) % Steady-state error for unit step reference input

% Part 5:
% Error for step disturbance with integral control
Gc = Kint/s % Controller transfer function
Gde = -R/Kt*feedback(Gp,Gc) % Error transfer function for disturbance
figure(9)
step(Gde)
figure(10)
pzmap(Gde)
edss = dcgain(Gde) % Steady state error for unit step disturbance

