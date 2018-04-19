addpath(genpath(cd));


%% SRC
%min_x \|Dx - y\|_1 + lambda \|x\|_1
opts_l1_fid_l1_reg.g_name = 'l1';
opts_l1_fid_l1_reg.g_opts.lambda = 0.1;
opts_l1_fid_l1_reg.J_name = 'l1norm';

%min_x \|Dx - y\|_2 + lambda \|x\|_1
opts_l2_fid_l1_reg.g_name = 'l1';
opts_l2_fid_l1_reg.J_name = 'l2norm';
opts_l2_fid_l1_reg.g_opts.lambda = 0.1;



%% different regularizers
%min_x \|Dx - y\|_2 + lambda \|x\|_0
opts_l2_fid_l0_reg.g_name = 'l0';
opts_l2_fid_l0_reg.J_name = 'l2norm';
opts_l2_fid_l0_reg.g_opts.lambda = 0.1;
opts_l2_fid_l0_reg.max_iter = 2000;

%min_x \|Dx - y\|_2 + lambda \|x\|_2
opts_l2_fid_l2_reg.g_name = 'l2';
opts_l2_fid_l2_reg.J_name = 'l2norm';
opts_l2_fid_l2_reg.g_opts.lambda = 0.1;

%min_x \|Dx - y\|_2 + lambda \|x\|_21
opts_l2_fid_l21_reg.g_name = 'l21';
opts_l2_fid_l21_reg.J_name = 'l2norm';
opts_l2_fid_l21_reg.g_opts.lambda = 0.1;
opts_l2_fid_l21_reg.g_opts.g = Gen_indset_labeled_data(s1_tra);

%min_x \|Dx - y\|_2 + lambda1 \|x\|_21 + lambda2 \|x\|_1

opts_l2_fid_l21_l1_reg.g_name = 'l21_l1';
opts_l2_fid_l21_l1_reg.J_name = 'l2norm';
opts_l2_fid_l21_l1_reg.g_opts.lambda1 = 0.05;
opts_l2_fid_l21_l1_reg.g_opts.lambda2 = 0.05;
opts_l2_fid_l21_l1_reg.g_opts.g = Gen_indset_labeled_data(s1_tra);

%min_x \|Dx_y\|_2 + lambda( alpha\|x\|_21 + (1-lambda) \|x\|_1)
opts_l2_fid_l2_l1_reg.g_name = 'l2_l1';
opts_l2_fid_l2_l1_reg.J_name = 'l2norm';
opts_l2_fid_l2_l1_reg.g_opts.lambda = 0.1;
opts_l2_fid_l2_l1_reg.g_opts.alpha = 0.5;

%% HS estimator
%min_x  l1l2( h(D_Gi x - y_Gi))  + lambda \|x\|_1
%   where h() is a huber estimator, l1l2() is 

opts_HS_fid_l1_reg.J_name = 'HS';
opts_HS_fid_l1_reg.J_opts.func_in_name = 'huber';
opts_HS_fid_l1_reg.J_opts.threshold_in = 0.005;
opts_HS_fid_l1_reg.J_opts.func_out_name = 'l1l2';
opts_HS_fid_l1_reg.is_debug = true;
%opts_HS_fid_l1_reg.g_opts.g = g;
opts_HS_fid_l1_reg.g_name = 'l1';
opts_HS_fid_l1_reg.g_opts.lambda = 0.1;
opts_HS_fid_l1_reg.rho = 0.08;
opts_HS_fid_l1_reg.mu = 1;
[x_HS_fid_l1_reg] = admm_main(D,y,opts_HS_fid_l1_reg);

%min_x  l1l2( h(D_Gi x - y_Gi))  + sum_g lambda1 \|x_g\|_21 + lambda2\|x\|_1
%   where h() is a huber estimator
opts_HS_fid_l21_l1_reg.J_name = 'HS';
opts_HS_fid_l21_l1_reg.J_opts.func_in_name = 'huber';
opts_HS_fid_l21_l1_reg.J_opts.threshold_in = 0.005;
opts_HS_fid_l21_l1_reg.J_opts.func_out_name = 'l1l2';
opts_HS_fid_l21_l1_reg.is_debug = true;
opts_HS_fid_l21_l1_reg.J_opts.g = g;
opts_HS_fid_l21_l1_reg.g_opts.g = Gen_indset_labeled_data(s1_tra);
opts_HS_fid_l21_l1_reg.g_name = 'l21_l1';
opts_HS_fid_l21_l1_reg.g_opts.lambda1 = 0.05;
opts_HS_fid_l21_l1_reg.g_opts.lambda2 = 0.05;
opts_HS_fid_l21_l1_reg.rho = 0.08;
opts_HS_fid_l21_l1_reg.mu = 1;



