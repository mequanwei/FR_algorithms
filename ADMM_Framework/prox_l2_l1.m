function [z] = prox_l2_l1(u, rho, opts)
%PROX_L2_L1 rho/2 \|z - u \|_2^2  +  g(z)
%   where g(z) = lambda (  alpha \|z_g\|_2 + (1-alpha) \|z\|_1 )

lambda = opts.lambda/rho;
alpha = opts.alpha;
indx = (u > lambda*(alpha - 1)) & (u < lambda*(1 - alpha));
z = ((alpha -1)*sign(u) + u/lambda) /(alpha + 1/lambda);
z(indx) = 0;


end

