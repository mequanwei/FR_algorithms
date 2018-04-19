function [cost] = reg_l21(x,g_opts)
%REG_L1 cost = lambda sum \|x_g\|_2
%  

lambda = g_opts.lambda;
g = g_opts.g;
cost = 0;
for i = 1:length(g)
    cost = cost + norm(x(g{i}));
end
cost = lambda*cost;

end

