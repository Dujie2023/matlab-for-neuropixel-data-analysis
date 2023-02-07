function acc = cal_svmacc_npdata(X, Y)
%input data
%X: firing rate in a timewindow for neurons recorded in the same session(ntrials*nunits)
%Y: label variable, (ntrials*1)

svm = fitcsvm(X,Y,'Standardize',true,'KernelFunction','linear');
cv = crossval(svm,'Kfold',10);
[label,~,~] = kfoldPredict(cv);
acc = sum(label==Y)/length(Y);
end