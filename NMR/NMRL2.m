function [x] = NMRL2(Dict, B, lambda, m , n)

% Dict:        N * d
% X:           N * 1
% B:           d * 1

Dict = Dict';
B = B';
[d,N] = size(Dict);
x = rand(N,1);

miu = 0.1; % see as step 
tao = 9;  % Threshold for singluar values

eps_abs = 0.005; 
eps_rel = 0.005;

Y = Dict*x-B;
Z = rand(d,1);


for i  = 1:100
%----------- x^{k+1} -------------

g = B + Y - (1/miu) * Z;
x = pinv(Dict'*Dict + (lambda/miu)*eye(N)) * Dict'*g;

%----------- Y^{k+1} -------------

Q = reshape((Dict*x - B + (1/miu)*Z),m,n);
M_Y = SVSO(tao,  Q);
Y_k = Y;
Y = reshape(M_Y, m*n, 1);

%----------- Z^{k+1} -------------
Z = Z + miu * (Dict*x -  Y - B);



M_Y = reshape(Y,m,n);
M_A_x = reshape(Dict*x - Y- B, m, n);
M_Z = reshape(Z, m, n);
dis(i) = sum(svd(M_Y,'econ')) + 0.5 * lambda * (norm(x,2))^2 + trace(M_Z'*M_A_x) + 0.5 * miu * norm(M_A_x, 'fro');



r_pri = Dict*x - Y - B;
s_dual = miu * Dict' * (Y_k - Y);  
eps_pri = sqrt(m*n) * eps_abs + eps_rel * max([norm(Dict*x,2), norm(Y,2), norm(B,2)]);
eps_dual = sqrt(n) * eps_abs  + eps_rel * norm(Dict'*Z,2);

if norm(r_pri,2) <= eps_pri && norm(s_dual,2) <= eps_dual
    break;
end

end 


plot(dis);
figure;
bar(x);

end

%Singular value shrinkage operator
function [Y] = SVSO(tao,Q)

[U,Sigma, V] = svd(Q,'econ');
Sigma = diag(Sigma);
%r = min(size(Q);
Sigma = diag(max(0,(Sigma-tao)));
Y = U*Sigma*V';

end