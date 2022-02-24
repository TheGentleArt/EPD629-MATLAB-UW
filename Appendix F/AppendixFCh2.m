% J. Williams
% Nise: Control Systems Engineering, 8th ed
% Appendix F Problems
% 2022-02-24

% Student Companion Site:
% https://bcs.wiley.com/he-bcs/Books?action=index&itemId=1119474221&bcsId=11568

% 'syms' + 'var' defines the var as a symbolic variable
% IE 'syms s' defines 's' as a symbolic variable

% Transfer functions can be written using symbolics.
% For instance, G=(s+1)∗(s+2)/[(s 2+3∗s+10)∗(s+4)] would
% replace the three statements, numg=poly([ 1 2]), deng=conv([1 3 10),
% [1 4]), and G=tf(numg,deng).

% In order to use symbolic calculation, one must first define the symbolic
% objects. IE 's' usually for laplace transforms, or 't' usually for time.
% Defining done by the 'syms' command. IE 'syms s' defines 's' as a
% symbolic object.
% The need to define symbolic variables only extends to the input
% variables. For instance, using the inverse laplace transform of some
% function with regards to a symbolic 's' will result in matlab creating a
% symbolic 't' by default.


% Ch2ApF1 (chapter 2 appendix F problem 1)
'(ch2apF1)'             % Display label.
syms s                  % Construct symbolic object for Laplace variable 's'.
'Inverse Laplace transform' % Display label.
F=2/[(s+1)*(s+2)^2];    % Define F(s) form case 2 example.
'F(s) from case 2'      % Display label.
pretty (F)              % Pretty print F(s)
f=ilaplace(F);          % Find inverse Laplace transform.
'f(t) for case 2'       % Display label.
pretty(f)               % Pretty print f(t) for Case 2.
F=3/[s*(s^2+2*s+5)];    % Define F(s) from Case 3 example.
'F(s) for Case 3'       % Display label.
pretty(F)               % Pretty print F(s) for Case 3.
f=ilaplace(F);          % Find inverse Laplace transform.
'f(t) for Case 3'       % Display label.
pretty(f)               % Pretty print f(t) for Case 3.
pause

% Ch2ApF2
% This example finds Laplace transforms of time functions using
% 'laplace(f)' command. This command yields F(s) in partial fractions.
% 
% Pretty print is not the only command useful for symbolic toolbox.
% Others include:
% collect(F) : collect common coefficient terms of F
% expand(F): expands product of factors of F
% factor(F): factors F
% simple(F): finds simplest for of F with the least number of terms
% simplify(F): simplifies F
% vpa(expression, places): stands for variable precision arithmetic,
% converts fractional symbolic terms to decimal with a set number of
% decimal places. 

'(ch2apF2)'             % Display label.
syms t                  % Construct symbolic object for time variable 't'.
'Laplace transform'     % Display label.
'f(t) from Case 2'      % Display label.
f=2*exp(-t)-2*t*exp(-2*t)-2*exp(-2*t); % Define f(t) from Case 2 example.
pretty(f)               % pretty print f(t) from case 2 example.
'F(s) for Case 2'       % Display label.
F = laplace(f);         % Find the laplace transform
pretty(F)               % Pretty print partial fractions of F(s) for case 2.
F = simplify(F);        % combine partial fractions
pretty(F)               % pretty print simplified laplace transform

'f(t) for Case 3'       % Display label.
f= 3/5 - 3/5*exp(-t)*[cos(2*t) + (1/2)*sin(2*t)]; % define f(t) for case 3
pretty(f)               % pretty print case 3 f(t) example
'F(s) for Case 3 - Symbolic Fractions' % Display label
F = laplace(f);         % Find laplace transform
pretty(F)               % Pretty print partial fraction of F(s) for case 3
'F(s) for Case 3 - Decimal representation' % Display label
F = vpa(F,3);           % Convert symbolic numerical fraction to 3-place decimal for F(s)
pretty(F)               % Pretty print decimal representatino of F(s) for case 3
'F(s) for case 3 - simplified' % Display label.
F = simplify(F);        % Combine partial fractions.
pretty(F)               % Pretty print combined partial fractions.
pause

