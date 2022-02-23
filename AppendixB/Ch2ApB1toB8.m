% J.Williams
% University of Wisconsin-Madison
% EPD629: Powertrain Systems and Controls
% Nise: Control Systems Engineering, Appendix B, Ch2
% 2022-02-23

% https://bcs.wiley.com/he-bcs/Books?action=mininav&bcsId=11568&itemId=1119474221&assetId=473619&resourceId=45902&newwindow=true
clear; clc;

% Ch2apB1
'(ch2apB1)'              % Display label.
' How are you?'          % Display string.
-3.96                    % Display scalar number -3.96.
-4 + 7i                  % Display complex number -4+7i.
-5-6j                    % Display complex number -5-6 j.
(-4+7i)+(-5-6i)          % Add two complex numbers and display sum.
(-4+7j)*(-5-6j)          % Multiply two complex numbers and display product.
M=5                      % Assign 5 to M and display.
N=6                      % Assign 6 to N and display.
P=M+N                    % Assign M+N to P and display.
Q=3+4j                   % Define complex number, Q.
MagQ=abs(Q)              % Find magnitude of Q.
ThetaQ=(180/pi)*angle(Q) % Find the angle of Q in degrees.
'Press space bar to move on...'
pause                    % Have user press key to skip to next portion?

% Ch2apB2
% Polynomials in s can be represented as row vectors containing the coefficients.
% Thus (P_1 = s^3 + 7*s^2 - 3*s + 23) can be represented by the vector shown below with elements
% separated by a space or comma.
'(ch2apB2)'             % Display label
P1 = [1, 7, -3, 23]     % Store polynomial 's^3 + 7*s^2 - 3*s + 23' as P1 and display result.
'Press space bar to move on...'
pause

% Ch2apB2
% Running the previous statements causes MATLAB to display the results.
% Ending the command with a semicolon suppresses the display. Typing an expression
% without a left-hand assignment and without a semicolon causes the expression to be
% evaluated and the result displayed. Enter P2 in the MATLAB Command Window after
% execution.
'(ch2apB3)'         % Display label.
P2 = [3, 5, 7, 8]   % Assign 3s^3 + 5s^2 +7s + 8 to P2 without displaying.
3*5                 % Evaluate 3*5 and display result.
'Press space bar to move on...'
pause

% Ch2apB4
%An F(s) in factored form can be represented in polynomial form. Thus
% P_3 = (s+2)*(s+5)*(s+6) can be transformed into a polynomial using poly (V), where
% V is a row vector containing the roots of the polynomial and poly(V) forms the coefficients
% of the polynomial.
'(ch2apB4)'             % Display label.
P3 = poly([-2, -5, -6]) % Store polynomial (s+2)(s+5)(s+6) as P3 and display the coefficients.
                        % I believe Poly inputs are the roots of the function.
                        % This means that the roots of the polynomial are
                        % known to be at (s+2)=0, (s+5)=0, and (s+6)=0, so
                        % at [-2, -5, -6].
                        % P3 is then returned as the coefficient of the
                        % polynomial that represents a solution to these
                        % roots. 
                        % In this case, the polynomial would be:
                        % s^3 + 13*s^2 + 52*s + 60
'Press space bar to move on...'
pause

% Ch2apB5
% We can find roots of polynomials using the roots (V) command. The
% roots are returned as a column vector. For example, find the roots of
% 5*s^4 + 7*s^3 + 9*s^2 + 3*s + 2 = 0.
'(ch2apB5)'             % Display label.
P4 = [5, 7, 9, 3, 2]    % Store polynomial 5*s^4 + 7*s^3 + 9*s^2 + 3*s + 2 as P4
rootsP4 = roots(P4)     % Find the roots of the polynomial
'Press space bar to move on...'
pause

% Ch2apB6
% Polynomials can be multiplied together using the conv(a,b) command
% (standing for convolve). Thus, P_5 = (s^3 + 7*s^2 + 10*s + 9)*(s^4 + 3*s^3 + 6*s^2 + 2*s + 1) is
% generated as follows:
'(ch2apB6)'     % Display label
P5 = conv([1, 7, 10, 9], [1, 3, 6, 2, 1]) % Convolve the two polynomials (same as multiplying the two polynomials by each other)
'Press space bar to move on...'
pause

% Ch2apB7
% The partial-fraction expansion for F s b s a s can be found using
% the [K, p, k ]= residue (b, a) command (K = residue; p = roots of denominator;
% k = direct quotient, which is found by dividing polynomials prior to performing a partial-
% fraction expansion). We expand F(s) = 7*s^2+9s+12)/(s*(s+7)(s^2+10s+100)) as an
% example. Using the results from MATLAB yields: F(s) =
% (0.2554-0.3382i)/(s+5.000-8.6603i)+0.2554+0.3382i)/(s+5.000+8.6603i)-(0.5280/(s+7)+0.0171/s)
'(ch2apB6)'                             % Display label
numf = [7, 9, 12];                      % Define the numerator of F(s)
denf = conv(poly([0, -7]), [1, 10, 100]) % Define denominator of F(s)
[K,p,k] = residue(numf,denf)            % Find residues and assign to K;
                                        % Find roots of denominator and
                                        % assign to p;
                                        % Find constant and assign to k.
'Press space bar to move on...'
pause

% Ch2apB8 (Example 2.3)
'(ch2apB8) Example 2.3'         % Display label.
numy = 32                       % Numerator of laplace transform
roots_deny = roots([1, 12, 32, 0]) % Roots of polynomial (coefficients) after multiplying both sides by s
deny = poly(roots_deny)         % Denominator of laplace transform
[r, p, k] = residue(numy, deny) % Calculate residues, poles, and direct quotient
                            % Residue is the K terms for equation (2.17) in
                            % the book. This is the coefficient of the
                            % exponential formula once an inverse laplace
                            % transform has taken place.
                            % The poles are the 'a' values of the
                            % exponential (e^(a*t)).
'Press space bar to move on...'
pause
