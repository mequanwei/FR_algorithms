function [cost] = reg_l21_l1(x,g_opts)
%REG_L1 cost = lambda1 sum \|x_g\|_2 + lambda2 \|x\|_1
%  

lambda1 = g_opts.lambda1;
lambda2 = g_opts.lambda2;

g = g_opts.g;
cost = 0;
for i = 1:length(g)
    cost = cost + norm(x(g{i}));
end
cost = lambda1*cost + lambda2*norm(x,1);

end
