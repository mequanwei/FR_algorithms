function [ acc,label ] = CESR_( tra,tes )
%CESR_ 此处显示有关此函数的摘要
%   此处显示详细说明
tra = tra';
tes = tes';
acc = 0;
A = tra(1:end-1,:);
clas_tra = tra(end,:);
clas = unique(clas_tra);

[~,tes_num]=size(tes); 
for i = 1:size(clas,2)
    ind{i} = (clas_tra==clas(i));
end
for i =1: tes_num
    label = tes(end,i);
    y = tes(1:end-1,i);
    x = CESR(A,y);
    cur_l = tes(end,i);
    for j = 1:size(clas,2)
        idx = ind{j};
        rual(j) = norm(A(:,idx)*x(idx)-y);
    end
    [~,no] = min(rual);
    predi(i) = clas(no);
    if clas(no) == label
        acc = acc+1;
    end
    %disp(['pecent:',num2str(i/tes_num)]);
end
acc = acc/size(tes,2);
end

