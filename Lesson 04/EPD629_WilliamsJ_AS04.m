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
% T - b*omega = J*(d/dt[omega])
% Laplace transform is then:
% T(s) - b*omega(s) = J*s*omega(s), which simplifies to:
% T(s) = (J*s + b)*omega(s)
% This can then give the following transfer function:
% omega(s)/T(s) = 1/(J*s + b)

G1 = 1/(J*s + b)    % Transfer function % G1 = omega(s)/T(s)
figure(1)
step(MotorTorque*G1)
hold on

% Recall first equation, T - b*omega = J*(d/dt[omega]) 
% re-order to have as function of rate of change of speed:
% d/dt[omega] = (T - b*omega)/J
% d/dt[omega] = (1/J)*T - (b/J)*omega, then reorder to have omega first:
% d/dt[omega] = (-b/J)*omega + (1/J)*T

A1 = -b/J
B1 = 1/J
C1 = 1
D1 = 0
Throttle1_ss = ss(A1, B1, C1, D1) % State space
step(MotorTorque*Throttle1_ss)

% Plot system poles using pzmap
figure(2)
pzmap(Throttle1_ss) % This should be identical to pzmap(G1)


% Part 2: Throttle with damping, without spring, output angular position
% d/dt[theta] = omega
% Laplace transform is then:
% s*theta(s) = omega(s)
% theta(s)/omega(s) = 1/s
% substituting back into original equation [omega(s)/T(s) = 1/(J*s + b)]:
% s*theta(s)/T(s) = 1/(J*s + b)
% theta(s)/T(s) = 1/(s*(J*s + b))
% theta(s)/T(s) = 1/(J*s^2 + b*s)

G2 = 1/(J*s^2+b*s+0)    % Transfer function % G1 = omega(s)/T(s)
figure(3)
step(MotorTorque*G2)
hold on

A2 = [0, 1; 0, -b/J]
B2 = [0; 1/J]
C2 = [1, 0]
D2 = 0
Throttle2_ss = ss(A2, B2, C2, D2) % State space for part 2
step(MotorTorque*Throttle2_ss)
figure(4)
pzmap(Throttle2_ss)




% Part 3: Throttle without motor, with spring, torque input, angular
% position output
G3 = 1/(J*s^2+b*s+K)    % Transfer function
figure(5)
step(MotorTorque*G3)
hold on
A3 = [0, 1; -K/J, -b/J]
B3 = [0; 1/J]
C3 = [1, 0]
D3 = 0
Throttle3_ss = ss(A3, B3, C3, D3) % State space for part 3
step(MotorTorque*Throttle3_ss)
figure(6)
pzmap(Throttle3_ss)


% Part 4: With motor, spring, voltage input, angular position output
% Motor torque proportional to current
% T = Kt*i % Motor torque
G4 = Kt / ((J*s^2 + b*s + K)*(L*s + R) + (Kb*Kt*s))
figure(7)
step(ArmatureVoltage*G4)
hold on
A4 = [0, 1, 0; -K/J, -b/J, Kt/J; 0, -Kb/L, -R/L]
B4 = [0; 0; 1/L]
C4 = [1, 0, 0]
D4 = 0
Throttle4_ss = ss(A4, B4, C4, D4) % State space for part 4
step(MotorTorque*Throttle4_ss)
figure(8)
pzmap(Throttle4_ss)
