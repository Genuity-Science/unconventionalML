% Helper function to predict outputs that will use to calculate metrics for R.
% function printPredForR(sols,traindatas,testdatas,filename). 
%
%   Parameters:
%   -------------
%   sols: Rx1 cell array of solutions, where R is the number of resamples. Each
%       cell of sols should be kx1 cell array, where k is the number of splits.
%       E.g., sols{1} could be a 3x1 cell array, where we used 3-fold cross
%       validation.
%   
%   traindatas: Rx1 cell array of resampled train data. The train data
%       that sols was trained on 
%   
%   testdatas: Rx1 cell array of resampled train data. The test data that 
%       generate predictions on.
%
%   filename: string. The filename to save the results. Saves the following:
%       y_trains: a matrix of the actual train labels. NxR matrix, where N is 
%       the size of the training data, R is the number of resamplings of the data
%   
%       y_tests: a matrix of the real test labels. N_test x R matrix, where 
%       N_test is the size of the test data, R is the number of resamplings of
%       the data. 
%
%       y_pred_trains: a matrix of the predicted probability of being in class 1.
%       Of size NxR. To get class: y_pred_trains >=0.5
%
%       y_pred_tests: a matrix of the predicted probability of being in class 1.
%       Of size NxR. To get class: y_pred_tests >= 0.5
%


function printPredForR(sols,traindatas,testdatas,filename,type)

if nargin < 5
    type = 'log';
end

sig = @(x)1./(1+exp(-x));

for m = 1 : length(traindatas)
    trdata = traindatas{m};
    tstdata = testdatas{m};
    y_trains(:,m) = trdata(:,1);
    y_tests(:,m) = tstdata(:,1);
    tmpsol = mean(cell2mat(sols{m}));
    if strcmp(type,'log')
        y_pred_trains(:,m) = sig(trdata(:,2:end)*tmpsol');
        y_pred_tests(:,m) = sig(tstdata(:,2:end)*tmpsol');
    else
        [~,~,y_pred_trains{m}] = getMultinomialAcc(tmpsol,traindatas{m});
        [~,~,y_pred_tests{m}] = getMultinomialAcc(tmpsol,testdatas{m});
    end
end

save(filename,'y_*');
