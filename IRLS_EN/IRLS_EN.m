function [alpha] = IRLS_EN(tra, y, n_class, lambda1, lambda2)
%n_class = 100,lambda1 = 0.1,lambda2 = 0.2;
[N,~] = size(tra);
n_per = int32(N/n_class);


A = tra(:,1:end-1);
A = A';
Label = tra(:,end);
alpha = (A'*A+0.1*eye(N))\(A'*y);
figure, bar(alpha);

for i = 1:20
   W1 = lambda1*(1./abs(alpha));
   %W1 = diag(W1);
   
   W2 = zeros(N,1);   
   for j = 1:n_class
       indx = boolean(zeros(N,1));
       indx(1+n_per*(j-1):n_per*j) = 1;
       idx_alpha = alpha(indx);
       w_cur = 1./norm(idx_alpha,2);
       W2(1+n_per*(j-1):n_per*j) = w_cur;
   end
   W2 = lambda2*W2;
   W = diag((W1+W2));
   
   alpha = (A'*A+W) \ (A'*y);
   ee = 0;
   for j = 1:n_class
       indx = boolean(zeros(1,N));
       indx(1+n_per*(j-1):n_per*j) = 1;
       idx_alpha = alpha(indx);
       ee = ee + norm(idx_alpha,2);
   end
   
   eng(i) = norm(y-A*alpha).^2 + lambda1*norm(alpha,1) + ee;
end

plot(eng);
figure, bar(alpha);

end