function [z] = prox_l21_l1(u, rho, opts)
%PROX_L21_L2 rho/2 \|z - u \|_2^2  +  g(z)
%   where g(z) = lambda_1  \sum_g \|z_g\|_2 + lambda_2 \|z\|_1

g = opts.g;
z = zeros(length(u),1);
paramLam1 =  opts.lambda1/rho;
paramLam2 =  opts.lambda2/rho;

hSoftThresholdL1 = @(vX, paramLambda) sign(vX) .* max(abs(vX) - paramLambda, 0);
hSoftThresholdL2 = @(vX, paramLambda) vX .* (1 - (paramLambda / (max(norm(vX, 2), paramLambda))));

for i = 1:length(g)
    ug = u(g{i});
    z(g{i}) = hSoftThresholdL2(hSoftThresholdL1(ug, paramLam1), paramLam2);
end


end

