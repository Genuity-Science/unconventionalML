function [h,J] = generateLogistic(data,type)
% generates h and J for approximate logistic loss.
% [h,J] = generateLogistic(data,type)
%   data: matrix of training data. Assumes first column is label. 
%   type: string. Default: 'ising'. Otherwise: 'qubo'. Treat the spins as
%       either QUBO or Ising.

if nargin < 2
    type = 'ising';
end
y = data(:,1);
D = data(:,2:end);

lh = D'*(0.5-y);
Q = 1/8*D'*D;

if strcmp(type,'ising')
    h = lh;
    J = triu(Q,1)*2;
elseif strcmp(type,'qubo')
    Q = Q + diag(lh);
    [h,J] = quboToIsing(Q);
elseif strcmp(type,'old')
    h = lh + diag(Q);
    J = triu(Q,1)*2;
else
    error('unrecognized type')
end




