function [params gammas loglike log_gamma_tilde] = vdpmm_expectation(data,params,infinite)
K = size(params.a,2);
eq_log_Vs = zeros(K,1);
eq_log_1_Vs = zeros(K,1);
log_gamma_tilde = zeros(size(data,1),K);
for i = 1:K
    eq_log_Vs(i) = digamma(params.g(i,1)) - digamma(params.g(i,1)+params.g(i,2));
    eq_log_1_Vs(i) = digamma(params.g(i,2)) - digamma(params.g(i,1)+params.g(i,2));
    log_V_prob(i) = eq_log_Vs(i) + sum(eq_log_1_Vs(1:(i-1)));
    
    eq_log_pi(i) = digamma(params.g(i,1)) - digamma(sum(params.g(:,1)));
    
    pob(:,i) = normwish(data,params.mean(i,:),params.beta(i),params.a(i),params.B(:,:,i));
    %disp(num2str(max(abs(pob - data_specific - normalizer))));
    if (infinite)
        log_gamma_tilde(:,i) = log_V_prob(i) + pob(:,i);
    else
        log_gamma_tilde(:,i) = eq_log_pi(i) + pob(:,i);
    end
end
gammas = exp(log_gamma_tilde);
gammas = gammas ./ repmat(sum(gammas,2),1,K);
%params.ll = sum(sum(gammas .* log_gamma_tilde));
loglike = sum(sum(gammas .* pob));
%this should involve some KL-divergences

%reestimating eq_alpha
h1 = 1; h2 = 1;
w1 = h1 + K - 1;
w2 = h2 - sum(eq_log_1_Vs(1:(end-1)));
params.eq_alpha = w1 / w2;
end