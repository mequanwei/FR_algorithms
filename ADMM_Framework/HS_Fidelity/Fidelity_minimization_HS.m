function [x] = Fidelity_minimization_HS(D, y, rho, z, opts)
%FIDELITY_MINIMIZATION_HS 
%    solving min_x 1/2 J(x) + 1/2 rho*\|x-z\|_2^2:
%               where J(x) = \|Dx - y\|_1 
%input:         opts  -  structure for HS fidelity including
%                        func_in_name:   -      string, name of inner func 
%                        threshold_in    -      thresholdolding of inner func         
%                        func_out_name:  -      string, name of out func
%                        threshold_out:  -      thresholdolding of out func
%                        g:              -      set of groups 
%                        is_debug:       -      0 or 1, if true, return
%                                              cost_tatol value
%                        max_iter:       -      max iteration time;


%% initialize variables 
[d,n] = size(D);

if ~exist('opts','var')
    opts = [];
end
%default setting
max_iter = 1;
func_in_name =  'huber';
func_out_name = 'huber';
g{1} =  boolean(zeros(d,1));
is_debug = 0;

% initialization
if isfield(opts,'max_iter');           max_iter = opts.max_iter;        end
if isfield(opts,'func_in_name');       func_in_name=opts.func_in_name;  end
if isfield(opts,'threshold_in');       threshold_in = opts.threshold_in;end
if isfield(opts,'func_out_name');      func_out_name=opts.func_out_name;end
if isfield(opts,'threshold_out');      threshold_out = opts.threshold_out;end
if isfield(opts,'g');                  g = opts.g;                      end
if isfield(opts,'is_debug');           is_debug = opts.is_debug;        end

size_g = length(g);
cost_in = zeros(size_g,1);
nomr_w = zeros(size_g,1);

for i = 1:size_g
    p{i} = zeros(size(y(g{i})));
end


% initialize opts for estimator function
func_in = str2func(func_in_name);
func_out = str2func(func_out_name);
opts_func_in.HQ_form = 'add';
opts_func_in.is_debug = true;               %NOTE: always be true
opts_func_out.HQ_form = 'mul';
opts_func_out.is_debug = false;

if exist('threshold_in','var');     opts_func_in.threshold_in = threshold_in;     end
if exist('threshold_out','var');    opts_func_out.threshold_out = threshold_out;  end
if isfield('is_debug','var');    opts_func_out.is_debug = is_debug;      end       

I = eye(n);
x = (D'*D + rho * I)\(D'*y + rho * z);
%x = rand(n,1);

for iter = 1:max_iter
    for i = 1:size_g
        D_Gi = D(g{i},:);
        y_Gi = y(g{i},:);
        e_Gi = D_Gi*x - y_Gi;
        [cost_in(i),p{i}] = func_in(e_Gi,opts_func_in);
        nomr_w(i) = sum(g{i})/d;
    end

    [~,w] = func_out(cost_in, opts_func_out);
    w = w .* nomr_w(i); 
    temp_inv = zeros(n,n);
    temp_right = zeros(n,1);
    for i = 1:size_g
        temp_inv = temp_inv + w(i)*(D(g{i},:)'*D(g{i},:));
        temp_right = temp_right + w(i)*D(g{i},:)'*(y(g{i})-p{i});
    end
    x = pinv(temp_inv + rho*I) * (temp_right + rho*z);
end

end

