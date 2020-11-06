function [train_idxs,test_idxs] = get_split_idxs(N,k)
% Returns k training and test splits of dataset of length N.
% [train_split,test_split] = get_split_idxs(N,k)
%
%   Parameters:
%   ----------
%
%   N: scalar. The size of the training data
%
%   k: scalar. The number of splits to use. Typically use 3. 
%
%   Returns:
%   --------
%
%   train_idxs: kx1 cell array of the train split indices
%
%   test_idxs: kx1 cell array of the test split indices

for n = 1 : k
    test_idxs{n} = floor((n-1)*N/k)+1:floor(n*N/k);
    train_idxs{n} = setdiff(1:N,test_idxs{n});
end


