% J. Williams
% Nise: Control Systems Engineering Ch2 "Try It" Problems
% 2022-02-24

% TryIt 2.1
F = zpk([1, 2], [-1, -2, -2], 2) % Zero-Poles-Gain model
                            % Used to define a transfer function
                            % First vector are the roots of numerator
                            % Second vector is the roots of the denominator
                            % Last scalar is the gain or what the numerator
                            % is multiplied by

% TryIt 2.2
numf = 2;
denf = poly([-1, -2, -2]) % would create a polynomial with given roots
[K, p, k] = residue(numf, denf) % Partial fraction expansion

% TryIt 2.3
F = tf([3], [1, 2, 5, 0]) % Transfer function
                          % tf(numerator/denominator)
                          % Enter the transfer function

% TryIt 2.4
syms s
f = ilaplace(3/(s*(s^2+2*s+5))) %Inverse laplace transform of transfer function defined from TryIt 2.3
pretty(f) % print f in a pretty format

% TryIt 2.5
%%%% did not verify this with the book
numf = 3 % numerator
denf = [1, 2, 5, 0] % denominator
[K, p, f] = residue(numf,denf)





