close all;
clear all;
clc;
dat_ad = [cd '\database\'];
addpath([cd '\l1_ls_matlab\']); % this code utilies l1_ls toolbox.
addpath([cd '\utilities\']);

% load subset 1,2,and 3 of Extended Yale B, the image size is 96*84.
% each image has unit l2-norm energy
load([dat_ad 'EY_Subset1_DAT.mat']); 
load([dat_ad 'EY_Subset2_DAT.mat']);
load([dat_ad 'EY_Subset3_DAT.mat']);

D          = [subset1_data subset2_data];
D_labels   = [subset1_labels subset2_labels];
Test_DAT   = subset3_data;
testlabels = subset3_labels;
clear subset1_data subset2_data subset3_data
Test_DAT   = Test_DAT;
testlabels = testlabels;

classids   = unique(D_labels);
classnum   = length(classids);

im_h       = 96;
im_w       = 84;

[disc_set,disc_value,Mean_Image]=Eigenface_f(D,260);
disc_value = sqrt((disc_value));
mean_x       =    255*(Mean_Image+0.001*disc_set*disc_value');

block_l    =   0.4;% input the occlusion level by you
height     =   floor(sqrt(im_h*im_w*block_l));
width      =   height;
num_c      =   size(Test_DAT,2);
% width_rand = rand(num_c,1);
% height_rand = rand(num_c,1);
% save rand_w_h width_rand height_rand;
load rand_w_h;    % the random position created by the author
w_a = 1;w_b=im_w-width+1;
r_w = w_a + (w_b-w_a).*width_rand;
h_a = 1;h_b=im_h-height+1;
r_h = h_a + (h_b-h_a).*height_rand;

sigma      =  1e-2;
ID         =  [];
lambda     =  0.005;

nit        =  1;
nIter      =  10;  % for simple, we just do 10 iterations.
BETA_A     =  [8];
MEDIAN_A   =  [0.7];

D          =  D./ repmat(sqrt(sum(D.*D)),[size(D,1) 1]);
beta_a     = BETA_A(1);
median_a   = MEDIAN_A(1);
   
for pro_i = 1:size(Test_DAT,2)
    
    I    =   reshape(Test_DAT(:,pro_i),[96 84]);
    [J]  =   Random_Block_Occlu(uint8(255.*I),r_h(pro_i),...
         r_w(pro_i),height,width); imshow(uint8(J));
     y = double(J(:));

     residual           =   (y-mean_x).^2;
     residual_sort = sort(residual);
     iter = residual_sort(ceil(median_a*length(residual))); beta = beta_a/iter; 
     w = 1./(1+1./exp(-beta*(residual-iter)));
     % w        =   sqrt(w);

     norm_y_D = norm(y);
     y = y./norm(y);

     for nit = 1: nIter

     tem_w = w./max(w);
     index_w = find(tem_w>1e-3);
     % remove the pixels with very small weight
     W_y = w(index_w).*y(index_w);
     W_D = repmat(w(index_w),[1 size(D,2)]).*D(index_w,:);
     
     temp_s = l1_ls(W_D,W_y,lambda*norm(W_y),sigma);
     
     residual = norm_y_D^2.*(y-D*temp_s).^2;
     residual_sort = sort(residual);
     iter = residual_sort(ceil(median_a*length(residual))); beta = beta_a/iter; 
     w = 1./(1+1./exp(-beta*(residual-iter)));
     % w        =   sqrt(w);
     end

    for class = 1:classnum
    s = temp_s (D_labels == class);
    z1 = w.*(y - D(:,D_labels == class)*s);
    gap1(class) = z1(:)'*z1(:);
    end
    
    index = find(gap1==min(gap1));
    ID = [ID index(1)];
end

reco_rate = sum(ID == testlabels)/length(testlabels);