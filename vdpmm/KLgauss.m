function kl = KLgauss(mu1,mu2,cov1,cov2)
mudiff = mu2 - mu1;
if (size(cov1,1)==1)
    kl = .5 * (sum(log(cov2)-log(cov1)) + sum(cov1 ./ cov2) + (mudiff ./ cov2) * mudiff' - numel(mu1));
else
    invcov = inv(cov2);
    kl = .5 * (logdet(cov2)-logdet(cov1) + trace(cov1*invcov) + mudiff * invcov * mudiff' - numel(mu1));
end