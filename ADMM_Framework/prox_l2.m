function [z ] = prox_l2(u, rho, opts)
%PROX_L2 min_z rho/2 \|z - u\|_2^2  +  g(z)
%        where g(z) = lambda \|x\|_2^2
%   
if ~exist('opts', 'var')
    opts = [];
end

lambda = rho / opts.lambda;
z = u / (1+2*lambda);

end

