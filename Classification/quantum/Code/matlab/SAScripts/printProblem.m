% Script to print problem in 3-column format, where first column is the first spin
% the second column is the second spin, and the third column is the value of the 
% coupling. If the first spin and the second spin are the same, the value is the 
% value of the local field. 
% 
%   printProblem(h,J,outFile,lambda,nFlag)

function printProblem(h,J,outFile,lambda,nFlag)

rng('shuffle')

if nargin < 3
    error('Must give name of output file');
end
if nargin < 4
    lambda = 0;
end
if nargin < 5
    nFlag = false;
end

fid = fopen(outFile,'w+');

h = h + lambda;
nQu = length(h);
J = J + J';
J = triu(J);
mat = J+diag(h);
if nFlag
    mat = mat + rand(size(mat))*0.05*max(max(abs(mat)));
end
fprintf(fid,'%d\n',nQu);
for j=1:nQu
    for k = j:nQu
        if ~mat(j,k) == 0
            fprintf(fid,'%d %d %.12f\n',j-1,k-1,mat(j,k));
        end
    end
end
fclose(fid);

