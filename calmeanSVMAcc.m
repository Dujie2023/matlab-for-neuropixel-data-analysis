function [acc] = calmeanSVMAcc(X,Y)
%cal svm acc, but for X in which nunits>10, 
%X: firing rate in a timewindow for neurons recorded in the same session(ntrials*nunits)
%Y: label variable, (ntrials*1)

nunits = size(X,2);
numunit_use = 10;%the number of unnits will be used to train SVM classifier

if nunits-numunit_use > 0

    repeattimes = ceil(nunits/numunit_use);
    acc_all = nan(1,repeattimes);

    for irepeat = 1:repeattimes
        [unituse_idx, ~] = crossvalind('LeaveMOut',nunits,nunits-numunit_use);
        acc_all(irepeat) = cal_svmacc_npdata(X(:,unituse_idx),Y);    
    end
    acc = mean(acc_all);

elseif  nunits-numunit_use == 0
    acc = cal_svmacc_npdata(X,Y);    
    
end

end

