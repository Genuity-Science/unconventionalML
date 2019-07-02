function ens = calcEns(sols,h,J)
% function ens = calcEns(sols,h,J)
% Calculates the energies of sols given h, J
%   sols: matrix of solutions (possible warning: beware of whether solutions
%       are in  {0,1} or {-1,1}. Most cases should probably be {-1,1}.
%   h: vector of local fields
%   J: matrix of couplings

ens = zeros(size(sols,1),1);

for n = 1 : size(sols,1)
    ens(n) = sols(n,:)*J*sols(n,:)' + sols(n,:)*h;
end
