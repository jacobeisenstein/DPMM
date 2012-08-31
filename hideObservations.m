% Copyright (C) 2007 Jacob Eisenstein: jacobe at mit dot edu
% distributable under GPL, see README.txt

function [params oldparams] = hideObservations(params, new_class, data)
%function [params oldparams] = hideObservations(params, new_class, data)
%
%updates params when a given set of observations are hidden.  this is
%necessary for Gibbs sampling

oldparams = params;

%update the suff stats
%see Sudderth's thesis, page 46

old_count = params.counts(new_class) + params.kappa;
old_sum = params.sums(new_class,:);
%params.SSE(:,:,new_class) = params.SSE(:,:,new_class) + old_sum' * old_sum / old_count;
params.cholSSE(:,:,new_class) = cholupdate(params.cholSSE(:,:,new_class),old_sum' / sqrt(old_count));

%must add iteratively because of the cholesky update function
for i = 1:size(data,1)
    params.counts(new_class) = params.counts(new_class) - 1;
    params.sums(new_class,:) = params.sums(new_class,:) - data(i,:);
 %   params.SSE(:,:,new_class) = params.SSE(:,:,new_class) - data(i,:)' * data(i,:);
    try
        params.cholSSE(:,:,new_class) = cholupdate(params.cholSSE(:,:,new_class),data(i,:)','-');
    catch
        disp('foo');
    end
end

new_count = params.counts(new_class) + params.kappa;
%params.SSE(:,:,new_class) = params.SSE(:,:,new_class) - new_sum' * new_sum / new_count;
params.cholSSE(:,:,new_class) = cholupdate(params.cholSSE(:,:,new_class),params.sums(new_class,:)' / sqrt(new_count),'-');
end