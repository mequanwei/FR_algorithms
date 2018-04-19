function [z] = prox_l1(u, rho, opts)
%PROX_L1 min_z rho/2 \|z - u \|_2^2  +  g(z)
%    where g(z) = lambda  \|z\|_1

lambda = opts.lambda/rho;
z = max(0,u-lambda)+min(0,u+lambda);

end

