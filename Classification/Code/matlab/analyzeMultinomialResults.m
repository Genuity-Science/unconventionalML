function [r,s] = analyzeMultinomialResults(sols,traindata,testdata,varargin)
% wrapper function for analysis
% [r,s] = analyzeMultinomialResults(sols,traindata,testdata)
%
%   sols: cell array of solutions, where the number of rows is the number of
%       folds, and the number of columns is the number of penalty strengths
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
%   lambdas : a 1xL array of penalty strengths
%
%   objfFlag: a character indicating the metric. False is auc, true is 
%       logistic loss
%
%   nSols: the number of solutions to examine. Default: 10000
%
%   ksplit: whether do a k shuffle and split. If true, then the total
%       number of splits is n_splits; i.e., the training split size is the 
%       the size of total training data divided by k (note that for typical
%       cross-validation, N/k is used for testing and the rest is used for
%       training). The test size split is set to 1/4 of the total training
%       data.
%
%   uniqueFlag: whether to only consider unique solutions. By default it is
%       false; i.e., it can kind of weight solutions by the frequency of
%       occurrence
%
%   testFlag: whether to automatically take the best performing results on
%       the test split
%
%   biasFlag: whether to manually calculate the best bias
%
%   r: a struct array. Contains various measures of performance on the training
%       data. Has the following fields:        
%
%       TrainItSols: n_splits x length(lambdas) cell array (see Input parameters
%           for n_splits). The solution with best performance on each 
%           training split, found by performing classical post-processing. 
%
%       TrainItObjf: n_splits x length(lambdas) cell array. The measure of 
%           the selected metric (see Input parameters) at each iteration.
%           Each cell is of size 1 x nSols.
%       
%       TestItSols: n_splits x length(lambdas) cell array (see Input parameters
%           for n_splits). The solution with best performance on each testing 
%           split, found by performing classical post-processing.
%
%       TestItObjf: n_splits x length(lambdas) cell array. The measure of the
%           selected metric (see Input parameters) at each iteration on the 
%           test split. Each cell is of size 1 x nSols.
%
%   s: struct of a few metrics on the test data. Has the following fields:
%
%       best_cal_acc: scalar. Accuracy of the element-wise averaged solution 
%           on the test data. Average is taking across the iterative solution
%           found on each train split. 
%
%       accs: 1xL numeric array. Accuracy of the averaged solution. 
%
%       best_idx: index of best performing value of \lambda. 
%
p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addParameter(p,'biasFlag',true,@islogical)
addParameter(p,'n_splits',3,validScalarPosNum)
addParameter(p,'lambdas',[0 2.^(-3:3)],@isnumeric)
addParameter(p,'uniqueFlag',false,@islogical)
addParameter(p,'testFlag',true,@islogical)
addParameter(p,'metric','acc',@ischar)
addParameter(p,'nSols',10000,validScalarPosNum);
addParameter(p,'type','ising',@ischar)
addParameter(p,'iterFlag',true,@islogical);
addParameter(p,'nClasses',6,validScalarPosNum)
parse(p,varargin{:});
params = p.Results;

r = getMultinomialFoldItSols(sols,traindata,params);

if params.testFlag
    itobjf = cellfun(@(x) x(end),r.TestItObjf);
    itsols = r.TestItSols;
else
    itobjf = cellfun(@(x) x(end), r.TrainItObjf);
    itsols = r.TrainItSols;
end

[bestitobjf,best_idx] = max(mean(itobjf));
meansol = mean(cell2mat(itsols(:,best_idx)));

accs = zeros(1,size(itsols,2));

for n = 1 : size(itsols,2)
    accs(n) = getMultinomialAcc(mean(cell2mat(itsols(:,n))),testdata);
end

s.best_cal_acc = getMultinomialAcc(meansol,testdata);
s.accs = accs;
s.best_idx = best_idx;

disp(sprintf(['best_cal_acc: ' repmat('%1.2f ',1,length(params.lambdas)) ' best idx: %d\n'], s.best_cal_acc,best_idx))
end

function r = getMultinomialFoldItSols(sols,traindata,params)
% assumes that number of folds are divided evenly

TrainItSols = cell(size(sols));
TrainItObjf = cell(size(sols));
TestItSols = cell(size(sols));
TestItObjf = cell(size(sols));
N = size(traindata,1);
K = size(sols,1);
if ~isa(sols{1},'double')
    sols = cellfun(@(x) double(x),sols,'uniformoutput',false);
end
for n = 1 : params.n_splits
    test_idx = floor((n-1)*N/K)+1:floor(n*N/K);
    train_idx = setdiff(1:N,test_idx);
    [h,J] = generateMultinomial(traindata(train_idx,:),params.type,params.nClasses);
    [h_t,J_t] = generateMultinomial(traindata(test_idx,:),params.type,params.nClasses);
    for m = 1 : size(sols,2)
        ens = calcEns(sols{n,m},h+params.lambdas(m)/2,J);
        [~,sidx] = sort(ens);
        tmpsol = sols{n,m}(sidx,:);
        if params.uniqueFlag
            tmpsol = unique(tmpsol,'rows','stable');
        end
        [TrainItSols{n,m},TrainItObjf{n,m}] = getMultiIterSols(tmpsol,traindata(train_idx,:),params);
        test_ens = calcEns(sols{n,m},h_t+params.lambdas(m)/2,J_t);
        [~,sidx] = sort(ens);
        tmpsol = sols{n,m}(sidx,:);
        if params.uniqueFlag
            tmpsol = unique(tmpsol,'rows','stable');
        end
        [TestItSols{n,m},TestItObjf{n,m}] = getMultiIterSols(tmpsol,traindata(test_idx,:),params);
%        disp(['Finished fold ' num2str(n) ' and lambda number ' num2str(m)])
    end
end

r.TrainItSols = TrainItSols;
r.TrainItObjf = TrainItObjf;
r.TestItSols = TestItSols;
r.TestItObjf = TestItObjf;

end

function [itSols,objf] = getMultiIterSols(sols,data,params)
% gets the iteratively best solution (i.e., only averages over solutions
% that reduce the objective function 
% [itSols,objf] = getIterSols(sols,data,order,N,labels)
labels = unique(data(:,1)); % !!!CHECK AND MAKE SURE THIS IS OK
if length(labels) ~= params.nClasses
    labels = 0:(params.nClasses-1);
end

L =  min(size(sols,1),params.nSols);
objf = zeros(1,L);
sol = sols(1,:);
sol(sol==3) = 0 ;
count = 2;

if strcmp(params.metric,'log')
    [~,~,pred_probs] = getMultinomialAcc(sol,data);
    objf(1) = calcLL(pred_probs,data(:,1),labels);
elseif strcmp(params.metric,'acc')
    objf(1) = getMultinomialAcc(sol,data,labels);
else 
    error('unrecognized metric')
end
for n = 2 : L
    if params.iterFlag
        dwtmp = sols(n,:);
        dwtmp(dwtmp==3) = 0 ;
        tmpsol = sol + (dwtmp-sol)/count;
        if strcmp(params.metric,'log')
            [~,~,pred_probs] = getMultinomialAcc(tmpsol,data);
            tmp = calcLL(pred_probs,data(:,1),labels);
        elseif strcmp(params.metric,'acc')
            tmp = getMultinomialAcc(tmpsol,data,labels);
        end
        if tmp > objf(n-1)
            objf(n) = tmp;
            sol = tmpsol;
            count = count+1;
        else
            objf(n) = objf(n-1);
        end
    else
        tmpsol = mean(sols(1:n,:),1);
%        objf(n) = calcObjF(tmpsol,data,params.fn,params.metric,lambda);
        if strcmp(params.metric,'log')
            [~,~,pred_probs] = getMultinomialAcc(tmpsol,data);
            tmp = calcLL(pred_probs,data(:,1),labels);
        elseif strcmp(params.metric,'acc')
            tmp = getMultinomialAcc(tmpsol,data,labels);
        end
        objf(n) = tmp;
        if objf(n) > max(objf(1:n-1))
            sol = tmpsol;
        end
    end
end
itSols = sol;

end

function ll = calcLL(pred_probs,real_labels,unique_labels)

ll = 0;
for n = 1 : length(unique_labels)
    ll = ll+sum(log(pred_probs(real_labels==unique_labels(n),n)));
end

end
