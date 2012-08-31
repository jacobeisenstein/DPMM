% Copyright (C) 2007 Jacob Eisenstein, jacobe at csail mit edu;
% distributable under GPL, see README.txt

function [h, hprime] = logalphapdfNI(alpha, k, n)
% Generate the log of the pdf for the alpha variable.
% similar to Rasmussen (2000) eq. 15, but with non-informative priors

try
%    alpha
%    log(alpha)
%    gammaln(alpha)
    h = k * log(alpha) + gammaln(alpha) - gammaln(n+alpha);
catch
    disp('logalphapdfNI fails');
end
hprime = k ./ alpha + psi(alpha) - psi(n+alpha);
