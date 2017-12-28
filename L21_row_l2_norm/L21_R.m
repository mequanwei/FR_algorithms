function [ rate ] = L21_R( train_data,test_data,row,col)
% L21 做保真项，行的形式，L2正则
%

[~,n] = size(train_data);
class = unique(train_data(:,n));%类
clas_tra = train_data(:,n);     %训练类标签
clas_tes = test_data(:,n);      %测试类标签
A_cell = cell(1,size(class,1)); %A_cell{i} 指subjecti所有训练数据

for k = 1:size(class,1)
    ind2 = (clas_tra == class(k));
    trak = train_data(ind2,:);
    A = cat(3,reshape(trak(1,1:n-1),row,col),reshape(trak(2,1:n-1),row,col));
    for it = 3:size(trak)
        A = cat(3,A,reshape(trak(3,1:n-1),row,col));
    end
    A_cell{k} = A;
end
right = 0;

for i = 1:size(class,1)
    ind = (clas_tes==class(i));
    tesi = test_data(ind,:);                    %第i类的Y
    class(i)
    for j = 1:size(tesi,1)
        Y = reshape(tesi(j,1:n - 1),row,col);   %第i类的第j个Y
        DIS = zeros(size(class,1),1); 
        for k = 1:size(class,1)
            A = A_cell{k};                                          
            [X,Wl] = train(A,Y,row);            %迭代           
            DIS(k) = Get_dis_row(X,A,Y);
        end
        [~,no] = min(DIS);
        if class(no) == class(i)
            right = right+1;
        end
    end    
end
rate = right /size(test_data,1);
end


function [X,Wl] = train(A,Y,row)
Wl = zeros(1,row) + 1/row;
[~,~,N] = size(A);
T = ones(N,1)*inf;
X = zeros(N,1);
i = 1;
while norm(T-X)>0.0001||i<11 
    i = i +1;
    T = X;
    X = Get_x_row(A,Y,Wl);
    [Wl] = Get_w_row(A,Y,X');
end
end

function [ X ] = Get_x_row( A,Y,W)
%get_x_row
%   此处显示详细说明
[m,n,pi] = size(A);
inv_sum_l = zeros(pi);
D_l = zeros(pi,n);
sum_l = zeros(pi,1);
for l = 1:m
    D = A(l,:,:); %D_l ：pi*m
    for i = 1:pi  % m*n的图像
        D_l(i,:) = D(:,:,i); 
    end
    inv_sum_l = inv_sum_l + W(l)*(D_l*D_l'); 
    sum_l = sum_l + W(l)*D_l * Y(l,:)';
end
X = pinv(inv_sum_l)*sum_l;

end


function [ W ] = Get_w_row( A,Y,X )
%GET_WI 此处显示有关此函数的摘要
%   此处显示详细说明

sigma = 0.001;
[m,n,pi] = size(A);
W = zeros(m,1);
D_l = zeros(pi,n);
for l = 1:m
    D = A(l,:,:); 
    for i = 1:pi  %m*n的图像
        D_l(i,:) = D(:,:,i);  %取D_l的值，行特征 pi*n维
    end
    W(l) = ((Y(l,:)-X*D_l)*(Y(l,:)-X*D_l)' + sigma)^(-1/2);
end

end




