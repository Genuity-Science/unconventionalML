function [SAsols,SAens,SAenCounts] = readSASolution(SAInput,SAOutput,rawFlag,N)
% Returns the energies, solutions, and from SA instance file and output file
%[SAsols,SAen,SAenCounts] = readSASolution(SAInput,SAOutput.rawFlag,N)
%
%   Parameters
%   ----------
%
%   SAInput: char. The filename of the SA input file
%
%   SAOutput: char. The filename of the output SA file.
%
%   rawFlag: boolean (optional). If rawFlag is true, then return all the 
%       solutions; i.e., not as a histogram. Default: true
%
%   N: int. The number of solutions to return. Default: total number of repetitions
%
%   Returns
%   -------
%   SAsols: N x M array of solutions, where N is specified as an input
%       and M is the number of spins. 
%
%   SAens: N x 1 array of energies. One energy for each spin configuration
%
%   SAencounts: N x 1 array of counts of each solution. If rawFlag is true
%       each entry will be 1. If rawFlag is false, then returns the counts for
%       each solution 
%

if nargin < 3
    rawFlag = true;
end

fid = fopen(SAOutput,'r+');

fseek(fid, 0, 'eof');
fileSize = ftell(fid);
frewind(fid);
%# Read the whole file.
data = fread(fid, fileSize, 'uint8');
%# Count number of line-feeds and increase by one.
len = sum(data == 10);
frewind(fid);


c=textscan(fid,'%f%f%f%s\n',len/2);
SAens = c{1};
SAenCounts = c{2};

if nargin < 4
    N = sum(SAenCounts);
end
sols=textscan(fid,'%s','delimiter','\n');
fclose(fid);
sols=sols{1};
nRQu = length(str2num(sols{1}));  % number of returned qubits
tmpSols = zeros(len/2,nRQu);

for k = 1 : len/2
    tmpSols(k,:) = str2num(sols{k});
end
tmpSols(:,nRQu) = [];

M = dlmread(SAInput);
nQu = M(1);
SAsols = ones(len/2,nQu)*3;

M(1,:) = [];
spins = M(:,1:2);
vals = M(:,3);
tmp = spins(vals==0,:);
missing = tmp(tmp(:,1)==tmp(:,2),1);
preorder = reshape(spins',2,[]);
inputorder = unique(preorder,'stable');
%missing = setdiff(1:nQu,inputorder);
%if (find(inputorder==0))  % now assumes that always starts with 0-based, otherwise doesn't work
    inputorder = inputorder+1;
%end
[s,order] = sort(inputorder);

SAsols(:,s) = tmpSols(:,order);
SAsols(:,missing) = 0; 

if rawFlag
    tmpsols = zeros(N,size(SAsols,2));
    tmpens = zeros(N,1);
    count = 1;
    for n = 1 : length(SAens)
        end_idx = count+SAenCounts(n)-1;
        tmpsols(count:end_idx,:) = repmat(SAsols(n,:),SAenCounts(n),1);
        tmpens(count:end_idx,:) = SAens(n);
        count = count+SAenCounts(n);
    end
    SAenCounts = ones(length(tmpens),1);
    SAsols = tmpsols;
    SAens = tmpens;
end
