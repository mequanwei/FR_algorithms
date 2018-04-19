function [cost] = HS(e,opts)
%HS caculate cost for HS
% input:            e       -      d*1, error
%                   opts    -      structure for HS, includes:
%                            func_in_name    -   name of inside function
%                            func_out_name   -   name of outside function
%                            threshold_in    -   inside function threshold
%                            threshold_out   -   outside function threshold
% output:           cost    -      1*1, total cost

d = length(e);

%default setting 
g{1} =  boolean(ones(d,1));
func_in_name = 'huber';
func_out_name = 'huber';
opts_in.debug = true;
opts_out.debug = true;


if ~exist('opts','var')
   opts = []; 
end
if isfield(opts,'func_out_name');    func_out_name = opts.func_out_name;          end
if isfield(opts,'func_in_name');     func_in_name = opts.func_in_name;            end
if isfield(opts,'threshold_out');    opts_out.threshold_out = opts.threshold_out; end
if isfield(opts,'threshold_in');     opts_in.threshold_in = opts.threshold_in;    end
if isfield(opts,'g');                g = opts.g;                                  end
 

func_in = str2func(func_in_name);
func_out = str2func(func_out_name);

len_g = length(g);
e_in = zeros(len_g,1);

for i = 1:len_g
    e_in(i) =  func_in(e(g{i}),opts_in);
end
cost = func_out(e_in,opts_out);

end

