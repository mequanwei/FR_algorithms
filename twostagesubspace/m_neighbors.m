function [ index, alfa, res] =m_neighbors( A, x, M)
%To find M neighbors of x by collabration presentation (Xu Yong 2011Transaction) 
%x=A'*alfa
mu=0.01;
[n_A,m_A]=size(A);
alfa=inv(A(:,1:(m_A-1))*A(:,1:(m_A-1))'+mu*eye(n_A))*A(:,1:(m_A-1))*x';
res=zeros(1,n_A);
% defined by Xu : M neighbors
% for i=1:n_A
%    res(i)=norm(x-alfa(i)*A(i,1:(m_A-1))) ;
% end

% first M significant samples
for i=1:n_A
   res(i)=-alfa(i) ;
end
  [sort_res,index_all]=sortrows(res');
%[sort_res,index_all]=sortrows(abs(alfa));
 index=index_all(1:M);%index_all(n_A+1-M:n_A)'
% index_all=find(alfa>0);
%  index=index_all;


