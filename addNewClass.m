% Copyright (C) 2007 Jacob Eisenstein: jacobe at mit dot edu
% distributable under GPL, see README.txt

function params = addNewClass(params)
%function params = addNewClass(params)
%adds a new, empty class to the dpmm

    newclassidx = params.num_classes+1;
    params.num_classes = newclassidx;
    params.counts(newclassidx) = 0;
    params.sums(newclassidx,:) = params.kappa * params.initmean;
    params.cholSSE(:,:,newclassidx) = chol(params.nu * params.initcov);
    %params.SSE(:,:,newclassidx) = params.nu * params.initcov;
end
