
function [best_bias,acc,bacc] = findBias(sol,data,metric,biasFnFlag,biasRange)
% Unused. Calculates bias when add to argument of sigmoid. Note that this is kind of a 
%   hack and didn't work that generalize that well.
% [best_bias,acc,bacc] = findBias(sol,data,metric,biasFnFlag,biasRange)
%   sol: the vector of weights
%   data: data matrix
%   metric: string which metric can use. Default: 'acc', can be 'bacc'. 
%       otherwise raises error
%   biasFnFlag: bool. If True means add bias within the sigmoid, otherwise add
%       after sigmoid. Default: true
%   biasRange: vector of biases to try. Default: -2:0.01:2

if nargin < 3
    metric = 'acc';
end
if nargin < 4
    biasFnFlag = true;
end
if nargin < 5
    biasRange = -2:0.01:2;
end
    sig = @(x)1./(1+exp(-x));
    accs = zeros(length(biasRange),1);
    baccs = zeros(length(biasRange),1);
    for n = 1 : length(biasRange) 
        [~,accs(n),baccs(n)] = predClasses(sol,data,biasRange(n),sig,biasFnFlag);
    end
    if strcmp(metric,'acc')
        [~,idx] = max(accs);
    elseif strcmp(metric,'bacc')
        [~,idx] = max(baccs);
    else
        error('unknown metric')
    end
    best_bias = biasRange(idx);
    acc = accs(idx);
    bacc = baccs(idx);
end


