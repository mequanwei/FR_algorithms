function mean_rate = PCANN(train, train_label, test, test_label)
%use pca to project the data to a low subspace and use NN to do
%classification

%preset_k=100;

[eigvector, eigvalue] = PCA(train');
total_eigvalue=sum(eigvalue);

%whether use variance to select k
per_variance_kept=0.95;
variance_kept=per_variance_kept*total_eigvalue;
tmp_variance=0;
k=0;
while tmp_variance<variance_kept
    k=k+1;
    tmp_variance=tmp_variance+eigvalue(k);    
end
disp(['k= ', num2str(k)])
%_____________________________ 
pro_vector=eigvector(:,1:k);
new_train=pro_vector'*train;
new_test=pro_vector'*test;
mean_rate=NN(new_train, train_label,new_test, test_label);

end

