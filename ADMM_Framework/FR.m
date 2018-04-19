function [acc] = FR(tra,tes,opts)
%FR  a general framework for linear representation based face recognition
%   此处显示详细说明

tr_l = tra(:,end); % train label: N*1
te_l = tes(:,end); % test label: num*1

tr = tra(:,1:end-1)';  % train data: d*N
te = tes(:,1:end-1)';  % test data: d*num
[~,N] =size(tr); 
[~,num] = size(te);
clas = unique(tr_l);   % class name
ind = boolean(zeros(N,size(clas,1)));
label = zeros(num,1);
dis = zeros(length(clas),1);

for i = 1:size(clas,1)
    ind(:,i) = (tr_l == clas(i)); % index the coressponding class train samples
end

for i = 1:num
    y = te(:,i);
    [x] = admm_main(tr,y,opts);
    opts.is_debug = 0;
    for j = 1:size(clas,1)
        x_ = x(ind(:,j));
        tr_ = tr(:,ind(:,j));
        dis(j) = compute_obj(tr_, y, x_, opts);
    end
    [~,No] = min(dis);
    
    label(i) = clas(No);
    disp(['processing: ',num2str(i),'/',num2str(num)]);
end

acc = sum(label == te_l')/num;
end


