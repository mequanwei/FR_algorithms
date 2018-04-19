function [x] = Fidelity_minimization_l1norm(D, y, rho, z, opts)
%FIDELITY_MINIMIZATION_L1NORM 
% solving min_x 1/2 J(x) + 1/2 rho*\|x-z\|_2^2:
%               where J(x) = \|Dx - y\|_1 
%               opts.max_iter:   -   max iter time

[~,n] = size(D);
I = eye(n);
%x = (D'*D + rho * I)\(D'*y + rho * z);
if ~exist('opts','var')
    opts = [];
end

max_iter = 1;
if isfield(opts,'max_iter');            max_iter = opts.max_iter;      end


x = rand(n,1);
for i = 1:max_iter
    B = diag(1./abs(D*x - y));
    x = (D'*B*D + rho*I)\(D'*B*y + rho*z);
end


%cost = 0;
%cost = norm(D*x - y, 1);

end

