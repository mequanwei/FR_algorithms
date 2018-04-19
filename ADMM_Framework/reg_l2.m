function [cost] = reg_l2(x,g_opts)
%REG_L1 cost = lambda \|x\|_2
%  
lambda = g_opts.lambda;
cost = lambda * norm(x,2);

end
