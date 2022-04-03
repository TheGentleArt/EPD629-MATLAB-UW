% J.Williams
% University of Wisconsin-Madison
% EPD629: Powertrain Systems and Controls
% Assignement 09
% 2022-04-02

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assignment:
% This homework assignment is to create root locus plots for an electronic
% throttle system and determine gains for a PI and a PID controller.
%  
% Part 1: Model Initialization File --- PI Control with Zero at -5rad/s (5 points)
%
% The remainder of the assignment is done via simulink.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that this was assignment had a large copy-paste from the homework
% assignment instructions...


clear;
close all;
clc;
 
syms s L R J b K Kt Kb
 
G1 = Kt/(L*s+R)
G2 = 1/(J*s)
G2 = G2/(1+G2*b)
G2 = G2/s
G2 = G2/(1+G2*K)
Gp = G1*G2/(1 + G1*G2*Kb*s)
simplify(Gp)
 
%%
clear
 
s = tf('s')  % Define s for symbolic transfer function format
 
% Calibration parameters
 
Kt = 0.063;  % Nm/A
Kb = Kt;  % V/(rad/s)
% L = 1e-3 was replaced by L = 1e-4 to keep poles real for this assignment
L = 1e-4;  % H  
R = 0.35;  % Ohm
J = 50e-6;  % kg*m^2
b = .001;  % Nm/(rad/s)
K = 0.1;  % Nm/rad
 
% Plant Model Transfer Function
 
G1 = Kt/(L*s+R)
G2 = 1/(J*s)
G2 = feedback(G2,b)
G2 = G2/s
G2 = feedback(G2,K)
Gp = feedback(G1*G2,Kb*s)
 
% Plant poles
 
p = pole(Gp)
p1 = p(1) % Gp pole 1
p2 = p(2) % Gp pole 2
p3 = p(3) % Gp pole 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following creates the root locus:

disp('Part 1: PI Control with Zero at -5 rad/s')

figure(1)
Gc_r1 = (s+5)/s % Gc for root locus with overall gain equal to one (xfer func)
rlocus(Gc_r1*Gp)
% This will be used to find the input for the gain used below...

% Plot the closed loop poles and zeros, as well as the step response:

Gain = 3.98  % Gain from breakaway point in root locus
Gc = Gain*Gc_r1 % Controller transfer function
[NUM, DEN] = tfdata(Gc);
Kp1 = NUM{1}(1)  % Kp1 for Simulink
Ki1 = NUM{1}(2)  % Ki1 for Simulink
T = feedback(Gc*Gp, 1) % Gc*Gp is open loop transfer function, T is closed loop (1 means unity gain)

figure(2)
pzmap(T)

figure(3)
step(T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Part 2: PI Control with Zero at -12 rad/s')
figure(4)
Gc_r2 = (s+12)/s
rlocus(Gc_r2*Gp)
Gain = 4.23  % Gain from breakaway point in root locus
Gc = Gain*Gc_r2 % Controller transfer function
[NUM, DEN] = tfdata(Gc);
Kp2 = NUM{1}(1)  % Kp1 for Simulink
Ki2 = NUM{1}(2)  % Ki1 for Simulink
T = feedback(Gc*Gp, 1) % Gc*Gp is open loop transfer function, T is closed loop (1 means unity gain)
figure(5)
pzmap(T)
figure(6)
step(T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Par 3: PI Control with Zero at Gp Pole 3 (-8.3691 rad/s)')
figure(7)
Gc_r3 = (s-p3)/s
rlocus(Gc_r3*Gp)
Gain = 4.1  % Gain from breakaway point in root locus
Gc = Gain*Gc_r3 % Controller transfer function
[NUM, DEN] = tfdata(Gc);
Kp3 = NUM{1}(1)  % Kp1 for Simulink
Ki3 = NUM{1}(2)  % Ki1 for Simulink
T = feedback(Gc*Gp, 1) % Gc*Gp is open loop transfer function, T is closed loop (1 means unity gain)
figure(8)
pzmap(T)
figure(9)
step(T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Part 4: PI Control with Zero at Gp Pole 3  w/ 0.7<damp_ratio<0.8')
figure(10)
Gc_r4 = (s-p3)/s
rlocus(Gc_r4*Gp)
Gain = 7.08  % Gain from 0.75 damping ratio in root locus plot
Gc = Gain*Gc_r4 % Controller transfer function
[NUM, DEN] = tfdata(Gc);
Kp4 = NUM{1}(1)  % Kp1 for Simulink
Ki4 = NUM{1}(2)  % Ki1 for Simulink
T = feedback(Gc*Gp, 1) % Gc*Gp is open loop transfer function, T is closed loop (1 means unity gain)
figure(11)
pzmap(T)
figure(12)
step(T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Part 5: PID Control with Zero at Gp Pole 2 and Gp Pole 3')
figure(13)
Tau = 0.001; % Filter time constant
Gc_r5 = (s-p2)*(s-p3)/(s*(Tau*s+1)) % Gc for root locus with overall gain=1
rlocus(Gc_r5*Gp)
Gain = 0.0551  % Gain from 0.75 damping ratio in root locus plot
Gc = Gain*Gc_r5 % Controller transfer function
[NUM, DEN] = tfdata(Gc);
Kd5 = NUM{1}(1)  % Kp1 for Simulink
Kp5 = NUM{1}(2)  % Kp1 for Simulink
Ki5 = NUM{1}(3)  % Ki1 for Simulink
T = feedback(Gc*Gp, 1) % Gc*Gp is open loop transfer function, T is closed loop (1 means unity gain)
figure(14)
pzmap(T)
figure(15)
step(T)


% Continue in Simulink...
