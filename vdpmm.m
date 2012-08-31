function [params gammas assign] = vdpmm(data, varargin)
%function [params assign] = vdpmm(data, varargin)
%
%This is an implementation of variational inference in multivariate gaussian
%dirichlet process mixture models.  It's based on Blei and Jordan 2005 and
%also W.D. Penny "Variational Bayes for d-dimensional Gaussian Mixture
%Models."  The notation (variable names) are largely from Penny.

[K infinite verbose maxits minits eps] = ...
    process_options(varargin,'K',50,'infinite',1,'verbose',1,...
    'maxits',500,'minits',10,'eps',.01);

[params gammas] = vdpmm_init(data,K);
numits = 2;
score = -inf;
score_change = inf;
while (numits < maxits) && (numits < minits || score_change > 1e-4)
    score_change = score;
    [params(numits) kldiv] = vdpmm_maximize(data,params(numits-1),gammas);
    [params(numits) gammas loglike] = vdpmm_expectation(data,params(numits),infinite);
    prev_ll = params(numits).ll;
    numits = numits+1;
    score = kldiv + loglike;
    score_change = (score - score_change) / abs(score);
    if (verbose), disp(sprintf('Iteration: %i\tscore: %f\tdelta: %f',numits,score,score_change)); end
end
[ig assign] = max(gammas');    
end
