function [z] = prox_l21(u, rho, opts)
%PROX_L21 rho/2 \|z - u \|_2^2  +  g(z)
%   where g(z) = lambda  \sum_g \|z_g\|_2

lambda = opts.lambda/rho;
g = opts.g;
z = zeros(length(u),1);


for i = 1:length(g)
    if norm(u(g{i})) > lambda
        z(g{i}) = u(g{i}) -  lambda* (u(g{i})) /norm(u(g{i}));
    end 
end

end

