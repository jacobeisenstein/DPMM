function [h, hprime] = logalphapdf(alpha, k, n)
% Generate the log of the pdf for the alpha variable in rasmussen's
% infinite gmm, eq 15.

% Copyright (C) 2005 Michael Mandel, mim at ee columbia edu;
% distributable under GPL, see README.txt

h = (k-3/2)*log(alpha) - 1./(2*alpha) + gammaln(alpha) - ...
    gammaln(n+alpha);
hprime = (k - 3/2)./alpha + 1./(2*alpha.^2) + psi(alpha) - ...
    psi(n+alpha);
