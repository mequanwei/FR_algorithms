function [cost] = reg_l2_l1(x,g_opts)
%REG_L1 cost = lambda ( alpha \|x\|_2 + (1-alpha) \|x\|_1 )
%  
lambda = g_opts.lambda;
alpha = g_opts.alpha;
cost = lambda * ( alpha * norm(x,2) + (1-alpha)* norm(x,1));

end