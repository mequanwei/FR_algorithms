function [ new_eigenvalue ] = eigenfeature_regularization( vector )

%vector:eigen values are sorted in descending order.
m=size(vector,2);
mu=1;
new_eigenvalue=vector;
if mod(m,2)==0
    median=(vector(m/2)+vector(m/2+1))/2;
else
    median=vector((m+1)/2);
end
d=median-vector(m);
index=find(vector<median+mu*d);
lamda_m=vector(index(1)-1);
lamda_1=vector(1);
alfa=lamda_1*lamda_m*(m-1)/(lamda_1-lamda_m);
beta=(m*lamda_m-lamda_1)/(lamda_1-lamda_m);
for i=(index(1)-1):m
   new_eigenvalue(i)=alfa/(i+beta);
   
end

end

