function [x] = admm_main(D,y,opts)
%ADMM_MAIN min_x J(x) + g(x) s.t. z = x
% 
%%
% 
% * input:
%           D   -   d*n dictionary
%           y   -   d*1 query
%           opts-   structure which field is:
%                   opts.J_name  -  name of J()
%                   opts.J_opts  -  structure for J()
%                                   details refer to specific J()
%                   opts.g_name  -  name of g()
%                   opts.g_opts  -  structure for g()
%                                   details refer to specific g()
%                   opts.rel     -  relative tolreance
%                   opts.abs     -  absolute toleracne 
%                   opts.rho     -  stepsize for dual variable updating in ADMM
%                   opts.mu     -   mu>=1, ratio used to increase rho
% * output:
%           x   -   n*1 coefficent
%%

eta = 1;

rho = 0.08; 
max_iter = 200;
eps_abs = 5e-4;
eps_rel = 5e-4;
max_rho = 1e10;
mu = 1;
is_debug = 0;
J_name = 'l2norm';
J_opts = [];
g_name = 'l2';
g_opts.lambda = 0.1;

if ~exist('opts','var')
    opts = [];
end

if isfield(opts, 'J_name');                J_name = opts.J_name;             end
if isfield(opts, 'J_opts');                J_opts = opts.J_opts;             end
if isfield(opts, 'g_name');                g_name = opts.g_name;             end
if isfield(opts, 'g_opts');                g_opts = opts.g_opts;             end
if isfield(opts, 'rho');                   rho = opts.rho;                   end
if isfield(opts, 'max_iter');              max_iter = opts.max_iter;         end
if isfield(opts, 'eps_abs');               eps_abs = opts.eps_abs;           end
if isfield(opts, 'eps_rel');               eps_rel = opts.eps_rel;           end
if isfield(opts, 'max_rho');               max_rho = opts.max_rho;           end
if isfield(opts, 'mu');                    mu = opts.mu;                     end
if isfield(opts, 'is_debug');              is_debug = opts.is_debug;         end


[d,n] = size(D);
%I = eye(n);
u = ones(n,1)*eta/rho;
z = pinv((D'*D)) * (D'*y);
Fidelity_minimization = str2func(['Fidelity_minimization_',J_name]);
prox_op = str2func(['prox_',g_name]);

if is_debug
    debug_opts.J_name = J_name;
    debug_opts.J_opts = J_opts;
    debug_opts.g_name = g_name;
    debug_opts.g_opts = g_opts;
    debug_opts.is_debug = is_debug;
end

%%
for i = 1:max_iter
    % -----------------updating x------------------
    % min_x 1/2 J(x) + 1/2 rho*\|x-(z-u)\|_2^2:
    x = Fidelity_minimization(D, y, rho, z-u, J_opts);
    
    % -----------------updating z------------------
    % min_z rho/2 \|z - (x+u) \|_2^2  +  g(z)
    z_old = z;
    z = prox_op(x+u, rho, g_opts);
    
    % -----------------updating u------------------
    u = u + x - z;
    
    % -------------stoping criterion---------------
    s_dual = rho*(z_old-z);
    r_pri = x - z;
    eps_pri = sqrt(d) * eps_abs + eps_rel * max([norm(x,2), norm(z,2)]);
    eps_dual = sqrt(n) * eps_abs  + eps_rel * norm(u*rho,2);
    if (norm(r_pri,2) <= eps_pri && norm(s_dual,2) <= eps_dual)
         break;
    end
    rho = min(rho*mu,max_rho);
    
    if is_debug
        dis(i) = compute_obj(D, y, x, debug_opts);
    end
end



end

