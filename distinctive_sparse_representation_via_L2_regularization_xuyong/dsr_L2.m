function [ accur,record0_class ] = dsr_L2( TRA,Y,gama)
%DSR_L2 此处显示有关此函数的摘要
%   此处显示详细说明
     %gama=1e-3;
     %class_num=40;
     %sample_num=10; 
     [tra_num,~] = size(TRA(:,1:end-1));                
     [tes_num,~] = size(Y(:,1:end-1));
     
     
    % test_num_per=tes_num/class_num; 
     %train_num_per=tra_num/class_num;
     
     errors=0;
     ex_data =TRA(:,1:end-1); 
     ex_data1=TRA(:,1:end-1)'; %训练数据
     data =Y(:,1:end-1); 
     data1=Y(:,1:end-1)';        %测试数据
     Label_test = Y(:,end)';
     clas = unique(TRA(:,end));
     for i = 1:size(clas,1)
        indx{i} = (clas(i) == TRA(:,end));
        train_num_per(i) = sum(indx{i});
        indx{i} = indx{i}';
     end
     
     class_num = size(clas,1);
     
     [size1,~]=size(ex_data1);
     M=eye(tra_num);
     for i=1:class_num
         xi = ex_data1(:,indx{i});
         %xi=ex_data1(:,(i-1)*train_num_per+1:i*train_num_per);
         %M((i-1)*train_num_per+1:i*train_num_per,(i-1)*train_num_per+1:i*train_num_per)=xi'*xi;
         M(indx{i},indx{i})=xi'*xi;
     end
     X=ex_data*ex_data';
     T=inv((1+2*gama)*X+2*gama*class_num*M);

     for i=1:tes_num
         y(:,1)=data1(:,i);
         solution=T*ex_data*y;

         contribution0=zeros(size1,class_num);

         for kk=1:class_num
             for hh=1:train_num_per(kk)
                 contribution0(:,kk)=solution((kk-1)*train_num_per(kk)+hh)*ex_data1(:,(kk-1)*train_num_per(kk)+hh)+contribution0(:,kk);
             end
         end
         for kk=1:class_num
             wucha0(kk)=norm(y-contribution0(:,kk));
         end
         [~,record00]=min(wucha0);
         record0_class(i)=clas(record00);

         if record0_class(i)~=Label_test(i)
             errors=errors+1;
         end
     end
     errors_ratio=errors/tes_num;
     accur=1-errors_ratio;
     
end

