function [ acc, label ] = NMRL2_FR( tra, tes, lambda, m,n)
%NMR_ classification
%   image size : m * n

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
    p = NMRL2(tr, y, lambda, m, n);
    
    for j = 1:size(clas,1)
        p_ = p(ind(:,j));
        tr_ = tr(ind(:,j), :);
        E = reshape(y-(p_'*tr_),m,n);
        dis(j) = sum(svd(E));
    end
    [~,No] = min(dis);
    label(i) = clas(No);
    disp(num2str(i/num));
end
    acc = sum(label == te_l)/num;

end

