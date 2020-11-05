function out = getPerf(sols,data,bias,fn,biasFnFlag,order)
% Evaluate some metrics on dataset
% out = getPerf(sols,data,bias,fn,biasFnFlag)
%
%   Parameters
%   ----------
%   sols: cell array of weights. Originally intended to be used on test dataset
%       with k-fold cross validation and L different regularization penalties; 
%       i.e., sols should be a kxL cell array. However, can be used for a 
%       single solution, though should then be a 1x1 cell array.
%
%   data: data matrix. First column assumed to be labels 
%   
%   bias: (optional) numeric, value of bias. Default: 0
%
%   fn: (optional) function handle, f, to take of inner product of sols and row
%       in data matrix. Default: sigmoid. 
%
%   biasFnFlag: (optional) bool. Default: true. If true, add bias inside fn; if
%       false, add bias after take fn; i.e., if true, y_pred = fn(w*x+bias);
%       if false, y_pred = fn(sols*x)+bias, where x is a row in data matrix,
%       w is the vector of weights
%
%   order: (optional). Default: [1,2]. The order of the class labels. If [1,2]
%       then means 0 is negative class, 1 is positive class. If [2,1] means
%       that 0 is positive class, 1 is negative class. Needed to define terms 
%       in the confusion matrix.
%
%   Returns
%   -------
%   out: a struct with the following fields:
%
%       AUCs: numeric array. The AUC for each cell in sols. Same size as sols
%       
%       meansolsAUC: numeric array. Takes the average of the solutions (down 
%           the column) and evaluate the AUC of the averaged solution. E.g., if
%           sols is a 3x7 cell array, meansols will be a numeric array of size 
%           1x7.
%
%       Accs: numeric array. The accuracy for each cell in sols. Same size as
%           sols.
%
%       meansolsAcc: numeric array. Takes the average of the solutions (down 
%           the column) and evaluate the accuracy of the averaged solution. 
%           E.g., if sols is a 3x7 cell array, meansols will be a numeric array
%           of size 1x7.
%
%       Baccs: numeric array. The balanced accuracy for each cell in sols. Same
%           size as sols.
%   
%       meansolsBacc: numeric array. Takes the average of the solutions (down     
%           the column) and evaluate the balanced accuracy of the averaged 
%           solution. E.g., if sols is a 3x7 cell array, meansols will be a 
%           numeric array of size 1x7.
 
if nargin < 3
    bias =0;
end
if nargin < 4
    fn =  @(x)1./(1+exp(-x));
end
if nargin < 5
    biasFnFlag = true; % true means add to raw probs, false means add to prob.
end
if nargin < 6
    order = [1,2];
end
AUC = zeros(size(sols));
meansolsAUC = zeros(1,size(sols,2));
acc = zeros(size(sols));
meansolsAcc = zeros(1,size(sols,2));
bacc = zeros(size(sols));
meansolsBacc = zeros(1,size(sols,2));
F1 = zeros(size(sols));
meansolsF1 = zeros(1,size(sols,2));

for n = 1 : size(sols,1)
    for m = 1 : size(sols,2)
        AUC(n,m) = -calcObjF(sols{n,m},data,fn,'auc');
        [~,acc(n,m),bacc(n,m),f1(n,m)] = predClasses(sols{n,m},data,bias,fn,biasFnFlag,order);
    end
end

for m = 1 : size(sols,2)
    meansol = mean(cell2mat(sols(:,m)),1);
    meansolsAUC(m) = -calcObjF(meansol,data,fn,'auc');
    [~,meansolsAcc(m),meansolsBacc(m),meansolsF1(m)] = predClasses(meansol,data,bias,fn,biasFnFlag,order);
end

out.AUCs = AUC;
out.meansolsAUC = meansolsAUC;
out.Accs = acc;
out.meansolsAcc = meansolsAcc;
out.Baccs = bacc;
out.meansolsBacc = meansolsBacc;
out.F1 = f1;
out.meansolsF1 = meansolsF1;

end
