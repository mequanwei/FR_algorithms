function [ acc,predi ] = l1_ls_fr( tra,tes,lambda )
% l1_ls ÈËÁ³Ê¶±ð

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
    x = l1_ls(A,y,lambda);
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
    disp(['pecent:',num2str(i/tes_num)]);
end
acc = acc/tes_num;
end

