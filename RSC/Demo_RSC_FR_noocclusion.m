% an example on AR database
close all;
clear all;
clc;

dat_ad     =    [cd '\database\'];
addpath([cd '\utilities\']);
load([dat_ad 'AR_DR_DAT']);
% you need create this mat by yourself using AR database.
% the mat 'AR_DAT_DAT'contains four dats:
% NewTrain_DAT: the training data in Session 1
% NewTest_DAT:  the testing data in Session 2
% trainlabels:  the training data's labels
% testlabels:   the testing data's labels
% the image size is 60*43, without any preprocessing

PROC_NUM          =   700;
NewTrain_DAT      =   NewTrain_DAT(:,1:PROC_NUM);
D_labels          =   trainlabels(1:PROC_NUM);
NewTest_DAT       =   NewTest_DAT(:,1:PROC_NUM-1);
testlabels        =   testlabels(1:PROC_NUM-1);
classids          =   unique(D_labels);
classnum          =   length(classids);

nIter             =   2;    % usually only need 2 iterations
residual          =   [];
im_h              =   60;
im_w              =   43;
eigen_num         =   300;  % the dimensionality of eigenface
lambda            =   0.001;
MEDC              =   [0.8];
BETC              =   [8];

[disc_set,disc_value,Mean_Image]=Eigenface_f(NewTrain_DAT,200);
disc_value       =   sqrt((disc_value));
mean_x           =   Mean_Image+0.001*disc_set*disc_value';

[disc_set,disc_value,Mean_Image]=Eigenface_f(NewTrain_DAT,eigen_num);

 ori_D             =  NewTrain_DAT;
 ori_D             =  ori_D./ repmat(sqrt(sum(ori_D.*ori_D)),[size(ori_D,1) 1]);

 median_c   =  MEDC(1);
 beta_c     =  BETC(1);
 ID         =  [];
        
for index_pro =1:size(NewTest_DAT,2)
    
        residual           =   (NewTest_DAT(:,index_pro)-mean_x).^2;
        residual_sort      =   sort(residual);
        iter               =   residual_sort(ceil(median_c*length(residual))); 
        beta               =   beta_c/iter; 
        w                  =   1./(1+1./exp(-beta*(residual-iter)));
        W                  =   diag(w);
     
        for nit = 1: nIter

        Ori_Train_DAT      =  W* NewTrain_DAT;
        Ori_Test_DAT       =  W* NewTest_DAT(:,index_pro);
        Ori_Train_DAT      =  Ori_Train_DAT./ repmat(sqrt(sum(Ori_Train_DAT.*Ori_Train_DAT)),[size(Ori_Train_DAT,1) 1]);
        Ori_Test_DAT       =  Ori_Test_DAT./ norm(Ori_Test_DAT);
        % [disc_set,disc_value,Mean_Image]=Eigenface_f(Ori_Train_DAT,eigen_num);

        Train_DAT         =  disc_set'*Ori_Train_DAT;
        Test_DAT          =  disc_set'*Ori_Test_DAT;

        D                 =  Train_DAT;
        D                 =  D./ repmat(sqrt(sum(D.*D)),[size(D,1) 1]);
        y                 =  Test_DAT;
        y                 =  y./norm(y);

        [x,status]        =  l1_ls(D,y,lambda);

        ori_y             =  NewTest_DAT(:,index_pro);
        norm_y            =  norm(ori_y,2);
        ori_y             =  ori_y./norm(ori_y,2);


        residual           =   norm_y^2*(ori_y-ori_D*x).^2;
        residual_sort      =   sort(residual);
        iter               =   residual_sort(ceil(median_c*length(residual))); 
        beta               =   beta_c/iter; 
        w                  =   1./(1+1./exp(-beta*(residual-iter)));
        W                  =   diag(w);
        end
        
        gap1 = [];
        for class = 1:classnum
         temp_s  =  x (D_labels == class);
         z1      =  y-D(:,D_labels == class)*temp_s;
         gap1(class) = z1(:)'*z1(:);
        end
        index = find(gap1==min(gap1));
        ID(index_pro) = index(1);

end

reco_rate = sum(ID == testlabels')/length(testlabels);



