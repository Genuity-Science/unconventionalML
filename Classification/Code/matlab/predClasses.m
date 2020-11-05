function [pred_classes,acc,bacc,f1] = predClasses(sol,data,bias,fn,biasFnFlag,order,possible_labels)
% predicts classes given solution and data for binary classification.
% [pred_classes,acc,bacc] = predClasses(sol,data,bias,fn,biasFnFlag)
%
%   Parameters
%   -----------
%   sol: vector of weights
%
%   data: data for which want to predict classes (and perhaps acc, bacc). First
%       column is assumed to be classes. Classes assumed to be in {0,1}.
%
%   bias: (optional) numeric. Default: 0. 
%
%   fn: (optional) function handle. How to transform inner-product of solution
%       with weights. Default: sigmoid function
%   
%   biasFnFlag: bool. Default: true. If true, add bias inside fn; if 
%       false, add bias after take fn
%
%   order: numeric. Default: [1,2]. If [1,2] means that 0 is negative label and
%       1 is positive. If [2,1] means that 0 is positive and 1 is negative.
%
%   possible_labels: numeric. Default: [0, 1]. Means that possible labels are 0 
%       or 1. 
%
%   Returns:
%   --------
%   pred_classes: vector of predicted classes labels 
%   
%   acc: numeric, accuracy on dataset 
%
%   bacc: numeric, balanced accuracy on dataset


if nargin < 3
    bias = 0;
end
if nargin < 4
    fn = @(x)1./(1+exp(-x));
end
if nargin < 5
    biasFnFlag = true; % true means that add the bias inside the sigmoid; false means after sigmoid 
end
if nargin < 6
    order = [1,2];
end
if nargin < 7
    possible_labels = [0 1];
end
if biasFnFlag 
    probs = fn(data(:,2:end)*sol'+bias);
else 
    probs = fn(data(:,2:end)*sol') + bias;
end

pred_classes = probs>=0.5;
acc = nnz(pred_classes == data(:,1))/length(pred_classes);
cm = confusionmat(data(:,1),double(pred_classes),'order',possible_labels);
cm = cm(order,order);
bacc = (cm(1,1)/sum(cm(1,:)) + cm(2,2)/sum(cm(2,:)))*0.5;
pr = cm(2,2)/(cm(2,2)+cm(2,1));
re = cm(2,2)/(cm(2,2)+cm(1,2));
f1 = 2*pr*re/(pr+re);
