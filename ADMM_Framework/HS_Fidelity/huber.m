function [cost,aux_variable] = huber(e,opts)
%HUBER caculate variables for huber HQ solution
%   input:      e          -         d_Gi*1 variables
%               opts       -         stureture for opts including:
%
%                      HQ_form   -    form of HQ, 'add' or 'mul'
%                      threshold -    thresholding for huber,if not
%                                   specified, we set (max|e|-min|e|)/3
%                      is_debug  -    if true, return cost;

len_e = length(e);
aux_variable = zeros(len_e,1);
s = aux_variable;
com_cost = true;

if ~exist('opts','var')
    opts = [];
end

if isfield(opts,'is_debug');       com_cost = opts.is_debug;      end

if isfield(opts,'threshold')
    th = opts.threshold;
else
    res = abs(e);
    th = (max(res) - min(res))*0.33 + min(res);
end

if isfield(opts, 'HQ_form')
    if strcmp(opts.HQ_form,'add')
        indx = (res>th);
        aux_variable(indx) = e(indx)-th*sign(e(indx));
    elseif strcmp(opts.HQ_form, 'mul')
        aux_variable(e<=th) = 1;
        aux_variable(e>th)  = th./(e(e>th));
    else
        assert(true,'HQ_form must be "add" or "mul".');
    end
end

if com_cost
    s(abs(e)<=th) = 0;
    s(e<-th) = -1;
    s(e>th) = 1;
    V = diag(1 - s.^2);
    cost = (1/2)* (e'*V*e) + th * (e'-(1/2)*th')*s;
else
    cost = 0;
end

end

