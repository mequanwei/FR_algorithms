function [ label_recorder distance_recorder,res_correct_rate] = FR_two_stage_subspace2(A,class_num,sample_num)
%Method2
%find M neghbors
M=120;
test_num=1;        %experiment times
train_num=7;       %samples used for trainning in each class
test_sample_num=19;%samples used in tests in each class
label_recorder=zeros(3,test_sample_num*class_num);
distance_recorder=zeros(3,test_sample_num*class_num);
[n_A, m_A]=size(A);
oneshot_average_correct_rate1=zeros(1,test_num);
oneshot_average_correct_rate2=zeros(1,test_num);
oneshot_average_correct_rate3=zeros(1,test_num);
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
        train_sample_index= [1,7,8,9,36,37,38];%[1,3,5,7,9,11,13,15];%1:5;
        train_num=size(train_sample_index,2);%samples used for trainning in each class
        test_sample_index=[4,27:35,56:64];%[2,4,6,8,10,12,14];%6:10;
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
    counter1=0; counter2=0; counter3=0;
    %===========================
    num_of_test_test_sample=test_sample_num*class_num;
     %query image testing
    for test_sample_No=1:num_of_test_test_sample
        test_sample=test_A(test_sample_No,1:(m_A-1))';
         [ index ] =m_neighbors( train_A, test_sample', M);
    
    index=sort(index);
    train_A_sec=train_A(index,:);
    [class_index,first_label]=unique(train_A_sec(:,m_A)','first');
    [class_index,last_label] =unique(train_A_sec(:,m_A)','last');
    class_num_sec=size(class_index,2);
    d1_vector=zeros(1,class_num_sec);
    d2_vector=d1_vector;
    %+++++++++Nearest Subspace++++++++++++++++++
    for i=1:class_num_sec
        dic_A=train_A_sec(first_label(i):last_label(i),1:m_A-1)';
        d1_vector(i)=norm(test_sample-dic_A*inv(dic_A'*dic_A)*dic_A'*test_sample);   %  +0.01*eye(size(dic_A,2))   
        clear dic_A
    end
    %+++++++++Farthest Subspace++++++++++++++++++
       for i=1:class_num_sec
        dic_A=train_A_sec(:,1:m_A-1)';
        dic_A(:,first_label(i):last_label(i))=[];
        d2_vector(i)=norm(test_sample-dic_A*inv(dic_A'*dic_A)*dic_A'*test_sample);    %    +0.01*eye(size(dic_A,2))
        clear dic_A
       end
        tmp_min_d1=min(d1_vector);    
        d1=class_index(find(d1_vector==tmp_min_d1));
        tmp_max_d2=max(d2_vector);
        d2=class_index(find(d2_vector==tmp_max_d2));
        d3_vector=d1_vector./d2_vector;
        tmp_min_d3=min(d3_vector);
        d3=class_index(find(d3_vector==tmp_min_d3));
        
        counter1=counter1+(d1==test_A(test_sample_No,m_A));
        counter2=counter2+(d2==test_A(test_sample_No,m_A));
        counter3=counter3+(d3==test_A(test_sample_No,m_A));
        label_recorder(:,test_sample_No)=[d1,d2,d3];
        distance_recorder(:,test_sample_No)=[tmp_min_d1,tmp_max_d2,tmp_min_d3];

    end
    oneshot_average_correct_rate1(test)=counter1/(test_sample_num*class_num);
    oneshot_average_correct_rate2(test)=counter2/(test_sample_num*class_num);
    oneshot_average_correct_rate3(test)=counter3/(test_sample_num*class_num);
    
end
average_correct_rate1=mean(oneshot_average_correct_rate1);average_correct_rate1_std=std(oneshot_average_correct_rate1);
disp(['correct rate for method 1:  ', num2str(average_correct_rate1),'  std=',num2str(average_correct_rate1_std) ])
average_correct_rate2=mean(oneshot_average_correct_rate2);average_correct_rate2_std=std(oneshot_average_correct_rate2);
disp(['correct rate for method 2:  ', num2str(average_correct_rate2),'  std=',num2str(average_correct_rate2_std)])
average_correct_rate3=mean(oneshot_average_correct_rate3);average_correct_rate3_std=std(oneshot_average_correct_rate3);
disp(['correct rate for method 3:  ', num2str(average_correct_rate3),'  std=',num2str(average_correct_rate3_std)])
res_correct_rate(:,1)=[average_correct_rate1,average_correct_rate2,average_correct_rate3]';
res_correct_rate(:,2)=[average_correct_rate1_std,average_correct_rate2_std,average_correct_rate3_std]';


