function [quot, remn, remd] = h_fraction(num,den)

% This function returns three values, the quotient and the remainder numerator and denominator
% where num and den are assumed to be integers
%
% num/den = quot+remn/remd  

quot = floor(num/den);
fact = gcd(num,den);
remd = den/fact;
remn = num/fact;
remn = mod(remn,remd);
