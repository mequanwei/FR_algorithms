function [cost,aux_variable] = welsch(e,opts)
%WELSCH caculate variables for welsh HQ solution
%   input:      e          -         d_Gi*1 variables
%               opts       -         stureture for opts including:
%
%                      HQ_form   -    form of HQ, 'add' or 'mul'
%                      sigma     -    sigma for huber,if not
%                                   specified, we sigma = 0.
%                      com_cost  _    if true, return cost
%

len_e = length(e);
aux_variable = zeros(len_e,1);
sigma2 = 0.25;
com_cost = true;


if ~exist('opts','var')
    opts = [];
end

if isfield(opts, 'sigma');        sigma2 = (opts.sigma)^2;        end
if isfield(opts,'is_debug');      com_cost = opts.is_debug;       end
if isfield(opts, 'HQ_form')
    if strcmp(opts.HQ_form,'add')
        aux_variable = e - e .* exp(- (e.^2 / sigma2));
    elseif strcmp(opts.HQ_form,'mul')
        aux_variable = exp(- (e.^2 ./ sigma2));
    else
        assert(true,'HQ_form must be "add" or "mul".');
    end
end

if  com_cost
    cost = sum (exp(- e.^2/0.25));
else
    cost = 0;
end
end

