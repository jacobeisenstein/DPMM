function [ params gammas] = vdpmm_init( data, K)
%function [ params gammas] = vdpmm_init( data, K)

[num dim] = size(data);

%initializing variational parameters
gammas = rand(num,K);
gammas = gammas ./ repmat(sum(gammas,2),1,K);
params(1).eq_alpha = 1;
params(1).beta = zeros(K,1);
params(1).a = zeros(K,1);
params(1).mean = zeros(K,dim);
params(1).B = zeros(dim,dim,K);
params(1).g = zeros(K,2);
params(1).ll = -inf;
end

