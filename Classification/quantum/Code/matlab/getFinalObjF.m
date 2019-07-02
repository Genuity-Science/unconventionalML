function mean_objf = getFinalObjF(results,field_name)
% Helper function to quickly get the final value of a particular objective 
% function for a 1-D struct array of results. Assumes that the field name
% is a cell array 
% mean_objf = getFinalObjF(results,field_name)
mean_objf = zeros(length(results),1);
for n = 1 : length(results)
    mean_objf(n) = mean(cellfun(@(x) x(end),results(n).(field_name)));
end
