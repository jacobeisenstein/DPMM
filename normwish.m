function logprob = normwish(x,mu,lambda,a,B)
%function logprob = normwish(x,mu,B,lambda,a)
%\int p(m,R) log p(x | m, R) d m d R
% mu ~ N(mu, lambda R)
% R ~ Wish(a,B)
%log probability of p(x|m,v)p(m|mu,aB)p(v|B,beta) or something...
[N D] = size(x);
Gamma_bar = .5 * a * inv(B);
d = x - squeeze(repmat(mu,N,1));
logprob = -.5 * D * log(2*pi) - .5 * logdet(.5*B) + ...
    .5 * sum(digamma(.5 * (a + 1 - (1:D)))) - .5 * D / lambda- ...
    sum((d * Gamma_bar).*d,2);
end

%   Precision = 0.5*hp_posterior.inv_B(:,:,c)*hp_posterior.eta(c);
%   E_log_p_of_x = - 0.5*D*log(pi) - 0.5*detln(hp_posterior.B{c}) ...
%       + 0.5*psi_sum(c) - 0.5*D/(hp_posterior.xi(c));
%     d = data.given_data - repmat(hp_posterior.m(:,c),1,N);
%     E_log_p_of_x = - sum(d.*(Precision*d),1) + E_log_p_of_x; % 1*N