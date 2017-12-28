function mean_rate = NN(train, train_label, test, test_label)
train_num=size(train,2);
test_num=size(test,2);
counter1=0;

for i=1:test_num
    dif=test(:,i)*ones(1,train_num)-train;
    dist=sum(dif.^2);
    [~,d1]=min(dist);
    counter1=counter1+(train_label(d1)==test_label(i));
end
mean_rate=counter1/test_num;


end

