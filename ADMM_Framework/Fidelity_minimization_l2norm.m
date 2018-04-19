function [x] = Fidelity_minimization_l2norm(D, y, rho, z, opts)
%FIDELITY_MINIMIZATION for L2 norm

% solving min_x 1/2 J(x) + 1/2 rho*\|x-z\|_2^2:
%               where J(x) = \|Dx - y\|_2^2 
%               opts = []   no any other variables
%%
[~,n] = size(D);
x = (D'*D + rho *eye(n))\(D'*y + rho * z);

end

