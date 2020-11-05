function printSACV(data,filestem,varargin)
% Prints SA input files
% function printSACV(data,filestem,varargin)
%
%   Parameters 
%   ----------
%   data: matrix from which to generate instances
%
%   filestem: string, base name of SA instances
%
%   Optional parameters
%   ------------------- 
%   Parameter type: default 'log'. Options: 'log','multi'
%
%   Parameter isingType: string, default 'ising'. Options: 'qubo', 'old'
%
%   Parameter n_folds: the number of folds. Default: 3
%
%   Parameter inSfx: A file suffix. Default: ''
%
%   Parameter lambdas: an array of regularizer strengths. Should be used with
%       l_str
%
%   Parameter l_str: a cell array of strings. Should be same length as lambdas.

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addParameter(p,'type','log',@ischar)
addParameter(p,'isingType','ising',@ischar)
addParameter(p,'n_folds',3,validScalarPosNum)
addParameter(p,'inSfx','',@ischar)
addParameter(p,'lambdas',[0 2.^(-3:3)],@isnumeric)
addParameter(p,'l_string',{'0','1d8','1d4','1d2','1','2','4','8'},@iscellstr);
parse(p,varargin{:});
params = p.Results;

assert(length(params.lambdas) == length(params.l_string))

if strcmp(params.type,'log')
    fn = @generateLogistic;
elseif strcmp(params.type,'multi')
    fn = @generateMultinomial;
else 
    error('unknown class type')
end

N = size(data,1);

for n = 1 : params.n_folds
    if params.n_folds ~= 1
        test_idx = floor((n-1)*(N/params.n_folds))+1:floor(n*N/params.n_folds);
        train_idx = setdiff(1:N,test_idx);
    else
        train_idx = 1 : N;
    end
    [h,J] = fn(data(train_idx,:),params.isingType);
    for m = 1 : length(params.lambdas)
        filename = [filestem '_Fold' num2str(n) 'L' params.l_string{m} params.inSfx];
        printProblem(h+params.lambdas(m)/2,J,filename)
    end
end
