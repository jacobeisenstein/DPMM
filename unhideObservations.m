% Copyright (C) 2007 Jacob Eisenstein, jacobe at csail mit edu;
% distributable under GPL, see README.txt

function params = unhideObservations(params, new_class, data)
%unhide some observations, updating suff. stats and counts
%see Sudderth's thesis, page 46

old_count = params.kappa + params.counts(new_class);

%for debugging only
%params.SSE(:,:,new_class) = params.SSE(:,:,new_class) + old_mean' * old_mean * old_count;
params.cholSSE(:,:,new_class) = cholupdate(params.cholSSE(:,:,new_class),params.sums(new_class,:)' / sqrt(old_count));

%must add iteratively because of the cholesky update function
params.sums(new_class,:) = params.sums(new_class,:) + sum(data,1);
params.counts(new_class) = params.counts(new_class) + size(data,1);
for i = 1:size(data,1)
 %   params.SSE(:,:,new_class) = params.SSE(:,:,new_class) + data(i,:)' * data(i,:);
    params.cholSSE(:,:,new_class) = cholupdate(params.cholSSE(:,:,new_class),data(i,:)');    
end
new_count = params.kappa + params.counts(new_class);
%params.SSE(:,:,new_class) = params.SSE(:,:,new_class) - new_mean' * new_mean * new_count;
params.cholSSE(:,:,new_class) = cholupdate(params.cholSSE(:,:,new_class),params.sums(new_class,:)' / sqrt(new_count),'-');
%params = handleRemovedClasses(params);
end