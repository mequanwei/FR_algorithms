function [cost,aux_variable] = l1l2(e,opts)
%L1L2  caculate the variables of l1l2 HQ solution
%   input:      e          -      d_Gi*1 variables
%               opts       -      stureture for opts including:
%                            HQ_form   -    form of HQ, 'add' or 'mul'
%                            is_debug  -    0 or 1, if true, calulate cost

len_e = length(e);
const = 1e-10;
com_cost = true;

aux_variable = zeros(len_e,1);
if ~exist('opts','var')
    opts = [];
end
if isfield(opts,'is_debug');        com_cost = opts.is_debug;       end

if isfield(opts, 'HQ_form')
    if strcmp(opts.HQ_form,'add')
        aux_variable = e - e./(sqrt(e.^2+const));
    elseif strcmp(opts.HQ_form,'mul')
        aux_variable = (e.^2 + const).^(-0.5);
    else
        assert(true,'HQ_form must be "add" or "mul".');
    end
end

if com_cost
    cost = sum(sqrt(e.^2+const));
else
    cost = 0;
end

end

