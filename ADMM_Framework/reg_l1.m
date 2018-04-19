function [cost] = reg_l1(x,g_opts)
%REG_L1 cost = lambda \|x\|_1
%  
lambda = g_opts.lambda;
cost = lambda * norm(x,1);

end

