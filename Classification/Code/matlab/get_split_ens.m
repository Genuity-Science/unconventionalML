function [ens,idxs,test_ens] = get_split_ens(sols,data,varargin)
% compares the energies of a kx1 cell array of dw_sols and sa_sols
% Energies calculated from data and split among k folds
% ens = compare_energies(sols,data,varargin)
%   sols: a kx1 cell array of solutions.
%   data: an NxM data matrix, where N is the number of samples, M is the
%       number of features
%   n_splits: the number of splits of the data for cross-validation
%   lambda : optional parameter. Default: 0
%   type: optional parameter. By default equal to 'ising'. Can also be 'qubo'
% 
%   Returns:
%   --------
%
%   ens: kx1 cell array of energies. The energies for each split (sorted).
%
%   idxs: kx1 cell array of indices. Sort order of the solutions; i.e.,
%       sols{i}(idxs{i},:) will return solutions in sorted order for i = 1 : K
%
%   test_ens: kx1 cell array. The energies on the test split.

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
addParameter(p,'n_splits',3,validScalarPosNum)
addParameter(p,'lambda',0,@isnumeric)
addParameter(p,'type','ising',@ischar)

parse(p,varargin{:});
params = p.Results;

assert(length(sols)==params.n_splits)
if ~isa(sols{1},'double')
    sols = cellfun(@(x) double(x),sols,'uniformoutput',false);
end
N = size(data,1);
if params.n_splits == 1 
    train_idx = 1:N;
    [h,J] = generateLogistic(data(train_idx,:),params.type);
    [ens,idxs] = sort(calcEns(sols{1},h+params.lambda/2,J));
else
    for n = 1 : params.n_splits
        test_idx = floor((n-1)*N/params.n_splits)+1:floor(n*N/params.n_splits);
        train_idx = setdiff(1:N,test_idx);
        [h,J] = generateLogistic(data(train_idx,:),params.type);
        [ens{n,1},idxs{n,1}] = sort(calcEns(sols{n},h+params.lambda/2,J));
        [h_t,J_t] = generateLogistic(data(test_idx,:),params.type);
        test_ens{n,1} = calcEns(sols{n},h_t+params.lambda/2,J_t);
    end
end


