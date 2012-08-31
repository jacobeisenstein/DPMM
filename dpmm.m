% Copyright (C) 2007 Jacob Eisenstein: jacobe at mit dot edu
% distributable under GPL, see README.txt

function params = dpmm(data, num_its, params)
%function params = dpmm(data, num_its, params)
%standard dirichlet process mixture model, with gaussian observations
%"rao-blackwellised" from, which does not store explicit means or covs

[T dimension] = size(data);

%some stats
debug = false;
allmean = mean(data,1);
allcov = cov(data);

if (~exist('params','var'))
    params(1).alpha = T / 50; %1 / wishrnd(1,1);
    params(1).kappa = .1; %T / 1000; %a pseudo count on the mean
    params(1).nu = 6; %a pseudo-count on the covariance
    params(1).initmean = allmean;
    params(1).initcov = allcov / 10;
    params(1).num_classes = 0;
    params(1).counts = 0;
    params(1).sums = [];
    params(1).cholSSE = [];
    params(1).classes = ones(T,1);
  %  params(1).SSE = [];
    params(1) = addNewClass(params(1));
    params(1) = unhideObservations(params(1), 1, data);
    if debug, if ~checkParams (params(1), data), disp('no check'); end, end
end

start_it = 1+size(params,2);

for it = start_it:(start_it+num_its-1)
    params(it) = params(it-1);
    disp (strcat (sprintf('%i [%i]: ',it,params(it).num_classes), ...
        sprintf(' %i',params(it).counts)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GIBBS SAMPLING %%%%%%%%%%%%%%%%%%%%%%%
    t_order = randperm(T);
    for sctr = 1:T        
        %t = sctr+1;%for debugging
        t = t_order(sctr);

        old_class = params(it).classes(t);

        params(it) = hideObservations(params(it),old_class,data(t,:));
        params(it) = handleRemovedClasses(params(it));
        if debug, if ~checkParams(params(it),data,t), disp('no check at hide'); end, end
        
        %these are the probabilities that we will sample from
        %note we add one to include the chance of adding a new class
        log_p_obs = -inf * ones(params(it).num_classes+1,1);

        p_prior = [];
        
        params(it) = addNewClass(params(it));  %it will be removed if nothing is put in it
        if debug, if ~checkParams(params(it),data,t), disp('no check at add class'); end, end
        
        kappabar = params(it).counts + params(it).kappa;
        nubar = params(it).counts + params(it).nu;
        factor = (kappabar + 1) ./ (kappabar .* (nubar - dimension - 1));
        p_prior = params(it).counts + params(it).alpha * (params(it).counts == 0);
        for i = 1:params(it).num_classes
%            if (params(it).counts(i) == 0), p_prior(i) = params(it).alpha; 
%            else p_prior(i) = params(it).counts(i); end
            try
                %integrating over the parameters of a
                %normal-inverse-Wishart yields student-t.  
                %this can be approximated by a "moment-matched" Gaussian, 
                %see sudderth thesis p 47
                %kappabar = params(it).counts(i) + params(it).kappa;
                %nubar = params(it).counts(i) + params(it).nu;
                %factor = (kappabar + 1) / (kappabar * (nubar - dimension - 1));
                log_p_obs(i) = normpdfln(data(t,:)', ...
                    params(it).sums(i,:)' / kappabar(i),...
                    sqrt(factor(i))*params(it).cholSSE(:,:,i));
            catch
                disp('mvnpdf throws error');
            end
        end
        %lightspeed sample normalizes automatically
        classprobs = p_prior'.*exp(log_p_obs-max(log_p_obs));
        try
            new_class = sample(classprobs);
            if (params(it).counts(new_class) == 0)
%                disp('adding a guy');
            end
            params(it).classes(t) = new_class;
        catch
            disp('could not sample');
        end
        if debug, if ~checkParams(params(it),data,t), disp('no check at sample'); end, end

        params(it) = unhideObservations(params(it),new_class,data(t,:));
        if debug, if ~checkParams(params(it),data), disp('no check at hide'); end, end
        
    end

    %%%%%%%%%%%%%%%%%%%% PARAMETER UPDATES %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %alpha is the "pseudo-count" for new classes.  it is estimated using ARS
    k = params(it).num_classes;
    n = T; 
    
    %can show that derivative is guaranteed to be positive / negative at
    %these points
    deriv_up = 2 / (n - k + 3/2);
    deriv_down = k * n / (n - k + 1);
    
    %this is the version with a conjugate inverse gamma prior on alpha, as
    %in Rasmussen 2000
    params(it).alpha = ars(@logalphapdf, {k, n}, 1, [deriv_up deriv_down], [deriv_up inf]);

    %this is the version with a totally non-informative prior
    %params(it).alpha = ars(@logalphapdfNI, {k, n}, 1, [deriv_up deriv_down], [deriv_up inf]);
end
end

%checks a set of parameters to see if they are self-consistent.
%for debugging
function [total c_basic c_count c_sum] = checkParams(params,data,exclude)
    if exist('exclude','var')
        c_basic = min(params.classes([1:exclude-1 exclude+1:end]) > 0);
    else
        c_basic = min(params.classes > 0);
    end
    c_count = 1;
    c_sum = 1;
    for i = 1:params.num_classes
        statedata = data(find(params.classes == i),:);
        err_amount = params.sums(i,:) - sum(statedata) - params.kappa * params.initmean;
        statecount = size(statedata,1);
        if exist('exclude','var')
            if i == params.classes(exclude) 
                err_amount = err_amount - data(exclude,:); 
                statecount = statecount - 1;
            end
        end
        if (statecount ~= params.counts(i)), c_count = 0; end
        if (sum(err_amount) > .01), c_sum = 0; end
    end        
    total = c_basic * c_count * c_sum;
end