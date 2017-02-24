function [ acc,label ] = CRC( tra,tes,lambda)
%CRC 此处显示有关此函数的摘要
%   此处显示详细说明
tr_l = tra(:,end);
te_l = tes(:,end);
tr = tra(:,1:end-1);
te = tes(:,1:end-1);
[N,~] =size(tr);
[num,~] = size(te);
clas = unique(tr_l);
ind = boolean(zeros(N,size(clas,1)));

for i = 1:size(clas,1)
    ind(:,i) = (tr_l == clas(i));    
end



label = zeros(num,1);
for i = 1:num
    y = te(i,:);
    p = pinv(tr*tr' + eye(N)*lambda)*(tr*y');
    for j = 1:size(clas,1)
        p_ = p(ind(:,j));
        tr_ = tr(ind(:,j),:);
        dis(j) = norm(y-(p_'*tr_),2)/norm(p_,2);
    end
    [~,No] = min(dis);
    label(i) = clas(No);
    disp(num2str(i/num));
end
    acc = sum(label == te_l)/num;
 
end

