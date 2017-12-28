function [res_correct_rate] =TPTSR_multi(A,class_num,sample_num,M,refine_time)

% Xu Yong 2011
%M=20;
mu=0.01;
test_num=252;        %experiment times
train_num=5;         %samples used for trainning in each class
test_sample_num=5;   %samples used in tests in each class
%refine_time=3;
[n_A, m_A]=size(A);
oneshot_average_correct_rate1=zeros(1,test_num);
test_table1=combntns(1:10,5);
for test=1:test_num
    disp(['test number=  ',num2str(test)])
    %to set matrix train_A train samples
    % and matrix test_A test samples
    tmp_index_test=zeros(1,test_sample_num*class_num);
    tmp_index_train=zeros(1,train_num*class_num);
    
    %method 1 randomly pick test samples
    %     for tmp_i=1:class_num
    %         tmp=randperm(sample_num);
    %         tmp_index_test(1+(tmp_i-1)*test_sample_num:tmp_i*test_sample_num)=tmp(1:test_sample_num)+(tmp_i-1)*sample_num;
    %         tmp_index_train(1+(tmp_i-1)*train_num:tmp_i*train_num)=tmp(test_sample_num+1:sample_num)+(tmp_i-1)*sample_num;
    %     end
    %   method 2 pick the test sample manually
    train_sample_index=test_table1( test,:);%[1,3,5,7,9,11,13,15];%1:5;
    train_num=size(train_sample_index,2);%samples used for trainning in each class
    test_sample_index=test_table1( 252+1-test,:);%[2,4,6,8,10,12,14];%6:10;
    test_sample_num=size(test_sample_index,2);%samples used in tests in each class
    for tmp_i=1:class_num
        tmp_index_test(1+(tmp_i-1)*test_sample_num:tmp_i*test_sample_num)= test_sample_index+(tmp_i-1)*sample_num;
        tmp_index_train(1+(tmp_i-1)*train_num:tmp_i*train_num)=train_sample_index+(tmp_i-1)*sample_num;
    end
    %method 3 "leave-one-out"
    %      tmp_index_test=test:sample_num:class_num*sample_num;
    %    tmp_index_train=setxor([1:class_num*sample_num], tmp_index_test) ;
    %====================================
    train_A=A( tmp_index_train,:) ;
    test_A=A( tmp_index_test,:) ;
    counter1=0;
    %===========================
    num_of_test_sample=test_sample_num*class_num;
    %query image testing
    for test_sample_No=1:num_of_test_sample
        test_sample=test_A(test_sample_No,1:(m_A-1))';
        train_A_old= train_A;
        neighbor_num=num_of_test_sample;
%         index_rec=1:200;

%         for k=1:refine_time+1
%            clear index_new train_A_new
%            neighbor_num=M+30*(refine_time-k+1);
%           [ index_new,alfa ] =m_neighbors( train_A_old, test_sample', neighbor_num);
%           index_new=sort(index_new);
%           
% %           
% %           index_rec=index_rec(index_new);
% %           figure,plot(index_rec,alfa(index_new))
% %   
%           
%           train_A_new=train_A_old(index_new,:); 
%           clear train_A_old
%           train_A_old=train_A_new;
%           clear train_A_new
%         end
        
        
        [class_index,first_label]=unique(train_A_old(:,m_A)','first');
        [class_index,last_label] =unique(train_A_old(:,m_A)','last');
        class_num_sec=size(class_index,2);
        
        alfa=inv(train_A_old(:,1:(m_A-1))*train_A_old(:,1:(m_A-1))'+mu*eye(neighbor_num))*train_A_old(:,1:(m_A-1))*test_sample;
        residual=zeros(1,class_num_sec);
        for i=1:class_num_sec  
            index_class_i=first_label(i):last_label(i);
            residual(i)=norm(test_sample-train_A_old(index_class_i,1:(m_A-1))'*alfa(index_class_i));
        end
        [min_res,index_min_res]=min(residual);
        d1=class_index(index_min_res);
        counter1=counter1+(d1==test_A(test_sample_No,m_A));
    end
    oneshot_average_correct_rate1(test)=counter1/(test_sample_num*class_num);
    
end
average_correct_rate1=mean(oneshot_average_correct_rate1);average_correct_rate1_std=std(oneshot_average_correct_rate1);
disp(['mean correct rate=  ', num2str(average_correct_rate1),'  std=',num2str(average_correct_rate1_std) ])
res_correct_rate=zeros(1,2);
res_correct_rate(1)=average_correct_rate1;
res_correct_rate(2)=average_correct_rate1_std;



