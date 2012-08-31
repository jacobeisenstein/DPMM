function samples = ars(logpdf, pdfargs, N, xi, support)

% Perform adaptive rejection sampling as described in gilks & wild
% '92, and wild & gilks 93.  The PDF must be log-concave.  Draw N
% samples from the pdf passed in as a function handle to its log.  The
% log could be offset by an additive constant, corresponding to an
% unnormalized distribution.  
%
% The pdf function should have prototype [h, hprime] = logpdf(x,
% pdfargs{:}), where x could be a vector of points, h is the value of
% the log pdf and hprime is its derivative.
% 
% The xi argument to this function is a number of points to initially
% evaluate the pdf at, which must be on either side of the
% distribution's mode.  And the support is a 2-vector specifying the
% support of the pdf, defaults to [-inf inf].
%
% This function does not use the lower squeezing bound because it
% is optimized for generating a small number of samples each call.

% Copyright (C) 2005 Michael Mandel, mim at ee columbia edu;
% distributable under GPL, see README.txt

samples = [];

% Don't need to approximate the curve too well, all the sorting and
% whatnot gets expensive
Nxmax = 50;

if(nargin < 5) support = [-inf inf]; end

x = sort(xi);
[h, hprime] = feval(logpdf, x, pdfargs{:});
if(~isfinite(h(1))) 
  x
  h
  hprime
  logpdf
  pdfargs{:}
  size(pdfargs)
  det(pdfargs{2})
  error('h not finite'); 
end

if(support(1) == 0)
  % Cheat!  Get closer and closer to 0 as needed
  while(hprime(1) < 0)
    xt = x(1)/2;
    [ht,hpt] = feval(logpdf, xt, pdfargs{:});
    [x,z,h,hprime,hu,sc,cu] = insert(x,xt,h,ht,hprime,hpt,support);
  end

  while(hprime(end) > 0)
    xt = x(end)*2;
    [ht,hpt] = feval(logpdf, xt, pdfargs{:});
    [x,z,h,hprime,hu,sc,cu] = insert(x,xt,h,ht,hprime,hpt,support);
  end
end

if(hprime(1) < 0 || hprime(end) > 0)
  % If the lower bound isn't 0, can't help it (for now)
  error(['Starting points ' num2str(x) ' do not enclose the' ...
	' mode']);
end


% Avoid under/overflow errors. the envelope and pdf are only
% proporitional to the true pdf, so we can choose any constant
% of proportionality.
offset = max(h);
h = h-offset;

[x,z,h,hprime,hu,sc,cu] = insert(x,[], h,[], hprime,[], support);

Nsamp = 0;
while Nsamp < N
  % Draw 2 random numbers in [0,1]
  u = rand(1,2);
  
  % Find the largest z such that sc(z) < u
  idx = find(sc/cu < u(1));  
  idx = idx(end);
  
  % Figure out the x in that segment that u corresponds to
  xt = x(idx) + (-h(idx) + log(hprime(idx)*(cu*u(1) - sc(idx)) + ...
      exp(hu(idx)))) / hprime(idx);
  [ht,hpt] = feval(logpdf, xt, pdfargs{:});
  ht = ht-offset;
  
  % Figure out what h_u(xt) is a dumb way, uses assumption that the
  % log pdf is concave
  hut = min(hprime.*(xt - x) + h);

  % Decide whether to keep the sample
  if(u(2) < exp(ht - hut))
    Nsamp = Nsamp+1;
    samples(Nsamp) = xt;
  else
% $$$     fprintf('.');
  end

  % Update vectors if necessary
  if(length(x) < Nxmax)
    [x,z,h,hprime,hu,sc,cu] = insert(x,xt,h,ht,hprime,hpt,support);
  end
end



function [x, z, h, hprime, hu, sc, cu] = ...
    insert(x, xnew, h, hnew, hprime, hprimenew, support)
% Insert xnew into x and update all other vectors to reflect the
% new point's addition.

[x,order] = sort([x xnew]);
h = [h hnew]; h = h(order);
hprime = [hprime hprimenew]; hprime = hprime(order);

z = [support(1) x(1:end-1)+(-diff(h)+hprime(2:end).*diff(x)) ./ ...
      diff(hprime) support(end)];
hu = [hprime(1) hprime] .* (z - [x(1) x]) + [h(1) h];

% $$$ plot(z, hu);

sc = [0 cumsum(diff(exp(hu)) ./ hprime)];
cu = sc(end);

