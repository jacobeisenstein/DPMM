% Copyright (C) 2007 Jacob Eisenstein: jacobe at mit dot edu
% distributable under GPL, see README.txt

function params = handleRemovedClasses(params)
%params = handleRemovedClasses(params)
if (~isfield(params,'num_fixed') || params.num_fixed == inf)
idxs = find(params(end).counts == 0);
for ctr = idxs
    %reduce all state numbers that are greater than ctr
    params.classes = params.classes - (params.classes >= ctr);
    idxs2 = [1:(ctr-1) (ctr+1):params.num_classes];
    params.counts = params.counts(idxs2);
    params.sums = params.sums(idxs2,:);
    params.cholSSE = params.cholSSE(:,:,idxs2);
    %        params.SSE = params.SSE(:,:,idxs);
    params.num_classes = params.num_classes-1;
end
end
end