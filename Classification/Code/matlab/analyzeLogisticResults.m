function [r, testperf] = analyzeLogisticResults(sols,traindata,testdata,varargin)
% wrapper function for analysis. Assumes that did cross-validation with k splits.
% Some terminology: 'test data' refers to the held out test data. The train data
% was further split into a 'train split' and a 'test split'. Classical post-
% processing is done on the train data as follows for each split:
%   1. Sort the solutions by their Ising energy on train split (not the entire
%       train data)
%   2. Examine higher energy solutions one at a time and see performance on 
%       a specified performance metric. 
%   3. Keep best performing solution 
% Allows for selection of whether to examine performance on the train split or
% test split. Returns two structs: one of measure of performance on the train 
% data (not split), and the other of performance on the test data (not split). 
%
%   Inputs:
%   --------
%   sols: cell array of solutions, where the number of rows is the number of
%       splits, and the number of columns is the number of penalty strengths
%       (i.e., size(sols,2) == length(lambdas), lambdas defined below)
%
%   traindata: MxN numeric array of training data, where M is the number of
%       training samples, and N is the number of features
%
%   testdata: numeric array of test data
%
%   Optional parameters:
%
%   n_splits: the number of folds used in cross-validation
%
%   fn: a function that maps the raw predictions to output probabilties. By
%       default it is the sigmoid, but if wanted to treat as a linear regression
%       instead of logistic regression, could be identity
%
%   lambdas : a 1xL array of penalty strengths, where L is the number of 
%       candidate regularization strengths that were tried.
%
%   metric: a character indicating the metric. Can either by 'auc', 'acc' or,
%       'log'. Default: 'auc'
%
%   nSols: the number of solutions to examine. Default: 1000. Sorts solutions
%       and includes higher energy solutions one at a time to see better
%       whether improves performance.
%
%   uniqueFlag: whether to only consider unique solutions. By default it is
%       false; i.e., it can kind of weight solutions by the frequency of
%       occurrence
%
%   biasFlag: whether to manually calculate the best bias. Default: true
%
%   isingFlag: whether problem treated as Ising or QUBO or old way. Default:
%       'ising'. Can also be 'qubo' or 'old'.
%
%   iterFlag: whether to use iterative solutions or simple averaging lowest 
%       nSols solutions. Default: true (iterative solutions)
%
%   postProcess: a string to say how to do post process. Can be 'test','train'. 
%       Default: train
%
%   order: the order of the classes. If [1,2] means 0 is the negative class and
%       1 is the positive class. If [2,1] means 0 is the positive class, 1 is 
%       negative class
%
%   Returns:
%   ---------
%   r: a struct array. Contains various measures of performance on the training
%       data. Has the following fields:        
%       TrainItSols: n_splits x length(lambdas) cell array (see Input parameters
%           for n_splits). The solution with best performance on each 
%           training split, found by performing classical post-processing. 
%
%       TrainItObjf: n_splits x length(lambdas) cell array. The measure of 
%           the selected metric (see Input parameters) at each iteration.
%           Each cell is of size 1 x nSols.
%       
%       TrainItEns: n_splits x length(lambdas) cell array. The energy of each
%           TrainItSols on each training split.
%       
%       TestItSols: n_splits x length(lambdas) cell array (see Input parameters
%           for n_splits). The solution with best performance on each testing 
%           split, found by performing classical post-processing.
%
%       TestItObjf: n_splits x length(lambdas) cell array. The measure of the
%           selected metric (see Input parameters) at each iteration on the 
%           test split. Each cell is of size 1 x nSols.
%       
%       TestItEns: n_splits x length(lambdas) cell array. The energy of each
%           TestItSols on each testing split.
%
%       TrainItSolOnTestBacc: The balanced accuracy on the test split (not 
%           the test data) of the TrainItSol.
%       
%       TestItSolOnTestBacc: balanced accuracy on the test split (not the 
%           test data) of TestItSol.
%
%       trainbiases: the bias that gave the best performance on train split. 
%           Mostly unused.
%   
%       testbiases: the bias that gave the best performance on the test split.
%           Mostly unused.
%
%       best_idx: the index of the best performing value of lambda.
%   
%       params: a struct array of the input parameters
%
%       trainmeansolsAcc: accuracy of the solution averaged across train 
%           or test split (depending on whether postProcess is 'train' or 
%           'test') on the entire train data.
%
%       trainmeansolsBacc: bal. accuracy of the solution averaged across train 
%           or test split (depending on whether postProcess is 'train' or 
%           'test') on the entire train data.
%
%       trainmeansolsF1: F1 of the solution averaged across train 
%           or test split (depending on whether postProcess is 'train' or 
%           'test') on the entire train data.
%
%   testperf: a struct array of performance on the train data. Uses the iterative
%       solution found in the post-processing approach. If postProcess is 'train'
%       then uses the solutions that performed best on the train split. If 
%       postprocess is 'test' then uses the solutions that performed best on the 
%       test split. Has the following fields:
%
%       AUCs: array of size n_splits x length(lambdas). The AUC on the test data
%           of each iterative solution. 
%
%       meansolsAUC: array of size 1 x length(lambdas). The AUC on the test data
%           found by taking the element-wise average of the n_splits iterative 
%           solutions.
%
%       Accs: array of size n_splits x length(lambdas). The accuracy on the test 
%           data of each iterative solution. 
%
%       meansolsAcc: array of size 1 x length(lambdas). The accuracy on the test
%           data found by taking the element-wise average of the n_splits 
%           iterative solutions.
%
%       Baccs: array of size n_splits x length(lambdas). The balanced accuracy 
%           on the test data of each iterative solution. 
%
%       meansolsBacc: array of size 1 x length(lambdas). The balanced accuracy
%           on the test data found by taking the element-wise average of the 
%           n_splits iterative solutions.
%
%       F1: array of size n_splits x length(lambdas). The F1 score on the test
%           data of each iterative solution. 
%
%       meansolsF1: array of size 1 x length(lambdas). The F1 score on the test
%           data found by taking the element-wise average of the n_splits 
%           iterative solutions.
%
%       meansolsLogloss: array of size 1 x length(lambdas). The log loss on the
%           test data found by taking the element-wise average of the n_splits 
%           iterative solutions.
%
%       bias: the bias used to evaluate performance measures. Unused, set to 0.
%
%       best_idx: the best index of the lambas, based on train data. Note that
%           it might not give the highest performance on the test data
%
%       best_meansolsBacc: scalar. Balanced accuracy on the test data of the 
%           averaged solution based on best performing \lambda on the train
%           data. Averaged solution refers to taking the element-wise average 
%           across the iterative solutions found during the post-processing on
%           training. 
%
%       best_meansolsAcc: scalar. Accuracy on the test data of the 
%           averaged solution based on best performing \lambda on the train
%           data. Averaged solution refers to taking the element-wise average 
%           across the iterative solutions found during the post-processing on
%           training. 
%
%       best_meansolsAUC: scalar. AUC on the testdata of the 
%           averaged solution based on best performing \lambda on the train
%           data. Averaged solution refers to taking the element-wise average 
%           across the iterative solutions found during the post-processing on
%           training. 
%
%       best_meansolsF1: scalar. F1 score on the test data of the 
%           averaged solution based on best performing \lambda on the train
%           data. Averaged solution refers to taking the element-wise average 
%           across the iterative solutions found during the post-processing on
%           training. 


p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
isfn = @(f) isa(f,'function_handle');
addParameter(p,'biasFlag',true,@islogical)
addParameter(p,'n_splits',3,validScalarPosNum)
addParameter(p,'lambdas',[0 2.^(-3:3)],@isnumeric)
addParameter(p,'uniqueFlag',false,@islogical)
addParameter(p,'iterFlag',true,@islogical)
addParameter(p,'metric','auc',@ischar)
addParameter(p,'nSols',1000,validScalarPosNum);
addParameter(p,'fn',@(x)1./(1+exp(-x)),isfn);
addParameter(p,'type','ising',@ischar);
addParameter(p,'verbose',false,@islogical);
addParameter(p,'postProcess','test',@ischar);
addParameter(p,'evalmetric','bacc',@ischar)
addParameter(p,'order',[1,2],@isnumeric);

parse(p,varargin{:});
params = p.Results;

r = getFoldItSols(sols,traindata,params);
%if params.bestOnTestFold
%    trainitobjf = -r.TrainItSolOnTestBacc;
%    testitobjf = -r.TestItSolOnTestBacc;
%else
    trainitobjf = cellfun(@(x) x(end), r.TrainItObjf);
    testitobjf = cellfun(@(x) x(end), r.TestItObjf);
%end
[besttrainobjf,best_train_ii] = min(nanmean(trainitobjf));
[besttestobjf,best_test_ii] = min(nanmean(testitobjf));

if strcmp(params.postProcess,'test')
    meansol = mean(cell2mat(r.TestItSols(:,best_test_ii)),1);
    sols = r.TestItSols;
    best_idx = best_test_ii;
elseif strcmp(params.postProcess,'train')
    meansol = mean(cell2mat(r.TrainItSols(:,best_train_ii)),1);
    sols = r.TrainItSols;
    best_idx = best_train_ii;
%elseif strcmp(param.postProcess,'hybrid')
%    [~,idx] = max([r.TrainItSolOnTestAcc],[],2
%
end

if params.biasFlag
    if strcmp(params.postProcess,'test')
        bias = mean(r.testbiases(:,best_idx));
    else
        bias = mean(r.trainbiases(:,best_idx));
    end
else
    bias = 0;
end
[~,r.trainmeansolacc,r.trainmeansolbacc,r.trainmeansolF1] = predClasses(meansol,traindata,bias,params.fn);

r.best_idx = best_idx;
r.params = params;
testperf = getPerf(sols,testdata,bias,params.fn,true,params.order);
testperf.meansolsLogloss = calcObjF(meansol,testdata,params.fn,'log');
testperf.bias = bias;
testperf.best_idx = best_idx;
testperf.best_meansolsBacc = testperf.meansolsBacc(best_idx);
testperf.best_meansolsAcc = testperf.meansolsAcc(best_idx);
testperf.best_meansolsAUC = testperf.meansolsAUC(best_idx);
testperf.best_meansolsF1 = testperf.meansolsF1(best_idx);
testperf.params = params;
if params.verbose
    disp(sprintf(['meansolsbacc: ' repmat('%1.2f ',1,length(params.lambdas)) ' best idx: %d\n'], testperf.meansolsBacc,best_idx))
end

end

function r = getFoldItSols(sols,traindata,params)
% Get iterative solution for each fold. Considers solutions one at a time and
% includes them only if results in increase of metric.
% r = getFoldItSols(sols,train_data,k,fn,lambdas,uniqueFlag,metric,nSols)
%   sols: cell array of solutions, where the number of rows is the number of
%       folds, and the number of columns is the number of penalty strengths
%       (i.e., size(sols,2) == length(lambdas), lambdas defined below)
%   train_data: MxN numeric array of training data, where M is the number of
%       training samples, and N is the number of features
%   k: the number of folds used in cross-validation
%   fn: a function that maps the raw predictions to output probabilties. By
%       default it is the sigmoid, but if wanted to treat as a linear regression
%       instead of logistic regression, could be identity
%   lambdas : a 1xL array of penalty strengths
%   uniqueFlag: whether to only consider unique solutions. By default it is
%       false; i.e., it can kind of weight solutions by the frequency of
%       occurrence
%   metric: a character indicating the metric. Can either by 'auc','acc', or 
%   'log'. Default: 'auc'
%   nSols: the number of solutions to examine. Default: 10000

if ~isa(sols{1},'double')
    sols = cellfun(@(x) double(x),sols,'uniformoutput',false);
end

if strcmp(params.type,'qubo')
    sols = cellfun(@(x) double(x==1),sols,'uniformoutput',false);
end
TrainItSols = cell(size(sols));
TrainItObjf = cell(size(sols));
TrainItSolOnTestBacc = zeros(size(sols));
TestItSolOnTestBacc = zeros(size(sols));
TestItSolonTestAUC = zeros(size(sols));
trainbiases = zeros(size(sols));
testbiases = zeros(size(sols));
N = size(traindata,1);
for n = 1 : size(sols,1)
    if params.verbose
        disp(['Fold ' num2str(n)])
    end
    if params.n_splits > 1
        test_idx = floor((n-1)*N/params.n_splits)+1:floor(n*N/params.n_splits);
        train_idx = setdiff(1:N,test_idx);
    else
        train_idx = 1:N;
    end
    [h,J] = generateLogistic(traindata(train_idx,:),params.type);
    [h_t,J_t] = generateLogistic(traindata(test_idx,:),params.type);
    for m = 1 : size(sols,2)
        lambda = params.lambdas(m);
        ens = calcEns(sols{n,m},h+lambda/2,J);
        [~,sidx] = sort(ens);
        tmpsol = sols{n,m}(sidx,:);
        if params.uniqueFlag
            tmpsol = unique(tmpsol,'rows','stable');
        end
        [TrainItSols{n,m},TrainItObjf{n,m}] = getIterSols(tmpsol,traindata(train_idx,:),lambda,params);
        TrainItEns{n,m} = calcEns(TrainItSols{n,m},h+lambda/2,J);
        if params.n_splits > 1
            [TestItSols{n,m},TestItObjf{n,m}] = getIterSols(tmpsol,traindata(test_idx,:),lambda,params);
            TestItEns{n,m} = calcEns(TestItSols{n,m},h_t+lambda/2,J_t);
        else
            TestItSols{n,m} = [];
            TestItObjf{n,m} = Inf;
            TestItEns{n,m} = nan;
        end
%        TrainItSolOnTestPerf(n,m) = calcObjF(TrainItSols{n,m},traindata(test_idx,:),params.fn,params.metric,lambda);
%        TrainItSolOnTestBacc(n,m) = -calcObjF(TrainItSols{n,m},traindata(test_idx,:),params.fn,'bacc',lambda);
        if params.biasFlag
            [trainbiases(n,m),~,TrainItSolOnTestBacc(n,m)] = findBias(TrainItSols{n,m},traindata(test_idx,:),'bacc');
            [testbiases(n,m),~,TestItSolOnTestBacc(n,m)] = findBias(TestItSols{n,m},traindata(test_idx,:),'bacc');
        else
            [~,~,TrainItSolOnTestBacc(n,m)] = predClasses(TrainItSols{n,m},traindata(test_idx,:));
            [~,~,TestItSolOnTestBacc(n,m)] = predClasses(TestItSols{n,m},traindata(test_idx,:));
            TestItSolonTestAUC(n,m) = calcObjF(TestItSols{n,m},traindata(test_idx,:),params.fn,'auc',lambda);
            trainbiases(n,m) = 0;
            testbiases(n,m) = 0;
        end

    end
end


r.TrainItSols = TrainItSols;
r.TrainItObjf = TrainItObjf;
r.TrainItEns = TrainItEns;
r.TestItSols = TestItSols;
r.TestItObjf = TestItObjf;
r.TestItEns = TestItEns;
r.TrainItSolOnTestBacc = TrainItSolOnTestBacc;
r.TestItSolOnTestBacc = TestItSolOnTestBacc;
r.trainbiases = trainbiases;
r.testbiases = testbiases;
end

function [itSols,objf] = getIterSols(sols,data,lambda,params)
% gets the iteratively best solution (i.e., only averages over solutions
% that reduce the objective function.
%
% [itSols,objf] = getIterSols(sols,data,orders,lambda,params)

L =  min(size(sols,1),params.nSols);
objf = zeros(1,L);
sol = sols(1,:);
sol(sol==3) = 0 ;
count = 2;
if strcmp(params.metric,'auc') && nnz(data(:,1)) == 0
    params.metric = 'bacc';
    disp('Not enough unique classes in data, switching to balanced accuracy...')
end
objf(1) = calcObjF(sol,data,params.fn,params.metric,lambda,params.order);
for n = 2 : L
    if params.iterFlag
        dwtmp = sols(n,:);
        tmpsol = sol + (dwtmp-sol)/count;
        
        tmp = calcObjF(tmpsol,data,params.fn,params.metric,lambda,params.order);
        if tmp < objf(n-1)
            objf(n) = tmp;
            sol = tmpsol;
            count = count+1;
        else
            objf(n) = objf(n-1);
        end
    else
        tmpsol = mean(sols(1:n,:),1);
        objf(n) = calcObjF(tmpsol,data,params.fn,params.metric,lambda,params.order);
        if objf(n) < min(objf(1:n-1))
            sol = tmpsol;
        end
%        sol = tmpsol;
    end
end
itSols = sol;

end
