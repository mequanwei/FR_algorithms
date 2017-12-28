function [ res ] = eigenface( data,para )
% each column of data is a face vector
% para decides the kept covariance
x=data;
n=size(data,1);
%PCA==========
k1=0;
[V,D] = eig(cov(x'));
eigenvalue=diag(D);
eigenvalue=eigenvalue/sum(eigenvalue);
eigenvalue=sort(eigenvalue,'descend');
contribution=0;

while contribution<para
    k1=k1+1;
    contribution=eigenvalue(k1)+contribution;
end
res=V(:,(n-k1+1):n);



end

