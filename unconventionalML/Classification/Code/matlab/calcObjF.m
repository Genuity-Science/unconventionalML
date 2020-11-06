function f = calcObjF(weights,data,func,metric,lambda,order)
% f = calcObjF(weights,data,func,metric,lambda,order)
% Helper function to calculate the desired metric. Used for a single set of 
%   weights at a time. Can be used with minimization, so lower means better.
%
%   weights: the matrix of weights (or solutions) to be used. 
%   
%   data: the training matrix which weights were trained on
%   
%   func: func handle. how to transform predictions. Can be identity (in which case
%       probably using as linear regression), or sigmoid (in which case 
%       probably using as logistic regression). Default: identity
%
%   metric: string describing which metric to use (could probably be changed 
%       to a function). Default: 'auc', or, 'log' (log-loss), 'acc', 'bacc'.
%       Raises error if unrecognized string
%
%   lambda: the strength of the regularization
%
%   order: numeric. Default: [1,2]. If [1,2] means that 0 is negative label and
%       1 is positive. If [2,1] means that 0 is positive and 1 is negative.
%
%   Returns:
%   --------
%   f: the performance measure on the input data
if nargin < 3
    func = @(x) x;
end

if nargin < 4
    metric = 'auc';
end

if nargin < 5
    lambda = 0;
end

if nargin < 6
    order = [1,2];
end
act = data(:,1);
o = func(data(:,2:end)*weights');
if strcmp(metric,'auc')
    try
        [~,~,~,AUC] = perfcurve(act,o,1);
        f = -AUC;
    catch
        f = nan;
    end

elseif strcmp(metric,'log')
    % Labels come in as {0,1}, but for this to be accurate should be {-1,1}
    act(act==0) = -1;
    f = sum(log(1+exp(-o.*act))) + lambda*weights*weights';
elseif strcmp(metric,'acc')
    [~,f] = predClasses(weights,data,0,func);
    f = -f;
elseif strcmp(metric,'bacc')
    [~,~,f] = predClasses(weights,data,0,func);
    f = -f;
elseif strcmp(metric,'f1')
    [~,~,~,f] = predClasses(weights,data,0,func,order);
    f = -f;
else
    error('unrecognized metric')
end

