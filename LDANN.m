function mean_rate = LDANN(train, train_label, test, test_label)
%use pca to project the data to a low subspace and use NN to do
%classification
%options.Fisherface=1;
options.Regu=1;
[eigvector,~] = LDA(train_label,options,train');
new_train=eigvector'*train;
new_test=eigvector'*test;
mean_rate=NN(new_train, train_label,new_test, test_label);

end

