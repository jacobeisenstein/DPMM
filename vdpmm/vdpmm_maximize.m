function [params kldiv] = vdpmm_maximize(data,params,gammas)
%priors
D = size(data,2);
K = numel(params.a);
a0 = D;
beta0 = 1;
mean0 = mean(data);
B0 = .1 * D * cov(data);

%convenience variables first
Ns = sum(gammas,1) + 1e-10;
mus = zeros(K,D);
sigs = zeros(D,D,K);
mus = gammas' * data ./ repmat(Ns',1,D);
%ag2 = shiftdim(repmat(sqrt(gammas),1,1,D),2);
for i = 1:K
    %            mus(i,:) = sum(repmat(gammas(:,i),1,D).*data) / Ns(i);
    diff0 = data - repmat(mus(i,:),size(data,1),1);
    diff1 = repmat(sqrt(gammas(:,i)),1,D) .* diff0;
    %diff2 = ag2(:,:,i)' .* diff0;
    sigs(:,:,i) = diff1' * diff1;
end

%now the estimates for the variational parameters
params.g(:,1) = 1 + sum(gammas);
%g_{s,2} = Eq[alpha] + sum_n sum_{j=s+1} gamma_j^n
params.g(:,2) = params.eq_alpha + ...
    flipdim(cumsum(flipdim(sum(gammas),2)),2) - sum(gammas);
params.beta = Ns + beta0;
params.a = Ns + a0;
params.mean = (repmat(Ns',1,D) .* mus + beta0 * repmat(mean0,size(Ns,2),1)) ./ repmat(Ns' + beta0,1,D);
for i = 1:K
    diff = mus(i,:) - mean0;
    params.B(:,:,i) = sigs(:,:,i) + Ns(i) * beta0 * diff * diff' / (Ns(i)+beta0) + B0;
end
kldiv = 0; %%%TODO
end
