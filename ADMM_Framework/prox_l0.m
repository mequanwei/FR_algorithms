function [z] = prox_l0(u, rho ,opts)
%PROX_L0 min_z rho/2 \|z - u \|_2^2  +  g(z)
%   where g(z) = lambda  \|z\|_0


lambda = 2 * opts.lambda/rho;
z = zeros(length(u),1);
z (abs(u) > sqrt(lambda))= u(abs(u) > sqrt(lambda));

end

