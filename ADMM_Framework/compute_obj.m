function [obj_val] = compute_obj(tr, y, x, opts)
%COMPUTE_OBJ
%                   


J_name = 'l2norm';
g_name = 'l2';
is_debug = 0;

if ~exist('opts','var')
    opts = [];
end


if isfield(opts,'J_name');      J_name = opts.J_name;           end
if isfield(opts,'J_opts');      J_opts = opts.J_opts;           end
if isfield(opts,'g_name');      g_name = opts.g_name;           end
if isfield(opts,'g_opts');      g_opts = opts.g_opts;           end
if isfield(opts,'is_debug');    is_debug = opts.is_debug;       end


e = tr*x -  y;
compute_fide = str2func(J_name);
compute_reg = str2func(['reg_',g_name]);

if ~exist('J_opts','var');      J_opts = [];            end
%if ~exist('g_opts','var');      g_opts = [];            end
J_val = compute_fide(e, J_opts);


if is_debug 
    g_val = compute_reg(x, g_opts);
    obj_val = J_val + g_val;
else
    obj_val = J_val;
end


end

