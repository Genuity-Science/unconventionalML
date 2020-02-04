% Generates solutions looking at the local fields only (trivial to solve)
% return N solutions based on energy.
% sols = getSolsBasedonFields(data)
% 
function sols = getSolsBasedonFields(data)

[h,J] = generateLogistic(data);
sols = -sign(h); 
