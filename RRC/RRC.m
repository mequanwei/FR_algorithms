function [ acc,label ] = RRC( tra,tes,m,n )
% input:	m: row size, 
%			n: column size



X = tra(:,1:end-1);
tr_l = tra(:,end);
train_or = X';
train = tra(:,1:end-1)';
test = tes(:,1:end-1)';
gttrain = tr_l';
gttest = tes(:,end)';
classes = unique(gttrain);
fr.dim.m = m;
fr.dim.n = n;
fr.Tor = train_or;
fr.T = train;
fr.alg = 'f-irc'; % 'f-irc', 'f-lr-irc'
fr.epsilon_3 = 0.01;
fr.kappa = 100;
% PARAMETERS ADMM//////////
fr.lambda_star = 0.01; % use 0.01 or 0.05
fr.rho1 = 1;
fr.rho2 = 0.1;
fr.epsilon_1 = 0.01;
fr.epsilon_2 = 0.1;
R = inv(fr.T'*fr.T + eye(size(fr.T'*fr.T,2)).*(fr.rho2/fr.rho1));
fr.Pinv = R;

for ii=1:size(test,2)
    
    % Regression
    [y,T,x] = runFIRC( fr,test(:,ii) );
    
    % Identification scheme
    for class  =  1:size(classes,2)
        s           =   x (gttrain == classes(class));
        residuals(class) = norm(y - T(:,gttrain == classes(class))*s, 2);
    end
    [val,identity] = min(residuals);
    label(ii) = classes(identity);
    
end

s = sum(label == tes(:,end)');
[N,~] = size(tes);
acc = s/N;

end

