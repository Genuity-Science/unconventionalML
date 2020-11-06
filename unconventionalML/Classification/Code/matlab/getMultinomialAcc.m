function [acc,pred_classes,pred_probs] = getMultinomialAcc(sol,data,labels)
% Returns the accuracy and optionally predicted classes for a given dataset and 
% set of labels.
%
% [acc, pred_classes,pred_probs] = getMultinomialAcc(sol,data,labels)
% 
%   sol: vector of weights
%
%   data: matrix of data on which to evaluate solution (not necessarily 
%       training). First column is assumed to be labels
%   
%   labels: (optional) vector of possible class labels. If not provided, uses
%       the unique values from the data matrix
%
% Returns:
%   acc: accuracy of dataset
%
%   pred_classes: vector of predicted classes, length is the number of rows of
%       data
%   
%   pred_probs: vector of raw probabilities for each class. 

if nargin < 3
    labels  = unique(data(:,1));
end

pred_probs = multinomial_probs(data(:,2:end),sol);
pred_classes = zeros(size(pred_probs,1),1);
maxm = max(pred_probs,[],2);
for n = 1 : size(pred_classes,1)
    idx = find(pred_probs(n,:)==maxm(n)); 
    if length(idx) > 1 % if multiple classes with same probability
        tmp = idx(randsample(length(idx),1)); 
    else
        tmp = idx;
    end
    pred_classes(n) = labels(tmp);
end

acc = nnz(pred_classes == data(:,1))/length(pred_classes);

end

function pred_class_probs = multinomial_probs(X,weights)
% Returns the probabilities for a multinomial of K classes. 
%
%   X is a matrix of training data of size NxD
%   weights is a vector of length D*(K-1) (the Kth class is calculated such that
%       the probabilities all sum to 1)

exptmpout = exp(X*reshape(weights,size(X,2),[]));
pred_class_probs = zeros(size(X,1),size(exptmpout,2)+1);
for n = 1:size(exptmpout,2)
    pred_class_probs(:,n) = exptmpout(:,n)./(1+sum(exptmpout,2));
end

pred_class_probs(:,end) = 1./(1+sum(exptmpout,2));

end
