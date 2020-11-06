function out = getSACVSols(filestem,k,inSfx,outPfx,l_str)
% Returns struct array of outputs from SA runs with fields sols, ens, probs
%
% out = getSACVSols(filestem,k,inSfx,outPfx)
%
%   filestem: string, base name of SA files
%   
%   k: numeric, number of folds
%
%   inSfx: string, suffix of input files
%
%   OutPfx: string, prefix of output files
%
%   Returns:
%   --------
%   out: a struct array with following fields:
%
%       sols: the matrix of solutions. NxM, where N is the number of repetitions
%           and M is the number of features
%   
%       ens: the energies of each solution in sols. Of size Nx1
%
%       probs: the probability of each solution. 
%
%   Example:
%   --------
%
%   out = getSACVSols('~/Dropbox-Work/Wuxi/SA/bootstrap_resamples/lumAB/Inst1','.txt','out.dat',{'0','1d4','1','4'});
%   Assumes that the files are saved in directory ~/Dropbox-Work/Wuxi/SA/bootstrap_resamples/lumAB
if nargin < 2
    k = 3;
end
if nargin < 3
    inSfx = '';
end
if nargin < 4
    outPfx = '';
end    
if nargin < 5
    l_str = {'0','1d8','1d4','1d2','1','2','4','8'};
end
% the line below is hard-coded but should depend on the strengh of
% regularization used 

%
out = struct();
sols = cell(k,length(l_str));
ens = cell(k,length(l_str));
probs = cell(k,length(l_str));
for n = 1 : k
    for m = 1 : length(l_str)
        filename = [filestem '_Fold' num2str(n) 'L' l_str{m}];
        fIn = [filename inSfx];
        fOut = [filename outPfx];
        [out(n,m).sols,out(n,m).ens,out(n,m).probs] = readSASolution(fIn,fOut);
    end
end
