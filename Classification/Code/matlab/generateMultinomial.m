function [h,J] = generateMultinomial(data,type,K)
% generates h and J for approximate multinomial loss.
% [h,J] = generateMultinomial(data,type)
%   data: matrix of training data. Assumes first column is label. 
%   type: string. Default: 'ising'. Otherwise: 'qubo'. Treat the spins as
%       either QUBO or Ising.

if nargin < 2
    type = 'ising';
end

class_labels = data(:,1);
D = data(:,2:end);
P = size(D,2);
labels = unique(class_labels);
%K = length(labels);
if nargin < 3
    K = length(labels);
end
hprime = sum(D,1)'/K;
X = D'*D;
Jp = (K-1)/(2*K*K)*X;
Jpp = -Jp/(K-1); 

tmp = triu(ones(K-1),1);
Q = kron(tmp,Jpp) + kron(eye(K-1),Jp);

bs = zeros((K-1)*P,1);
for n = 1 : K-1
    label = labels(n);
    bs((n-1)*P+1:n*P) = -sum(D(class_labels==label,:),1);
end
if strcmp(type,'ising')
    h = repmat(hprime,K-1,1) + bs;
    J = triu(Q,1)*2;
else
    Q = Q + diag(h);
    [h,J] = quboToIsing(Q);
end


