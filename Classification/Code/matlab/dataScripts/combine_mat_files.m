function combine_mat_files(input_name_str,output_name,N) 

% Combine .mat files with multiple cuts. Did not figure out how to output
% single .mat file with cell array in R, so have to do it as a separate script.
% 
%   Inputs:
%   -------
%   input_name_str: str, with format specification. To be used with sprintf
%
%   output_name: str, the name of the mat file to save after combining files
%
%   N: int, the number of datasets to combine 
%
%   Example:
%   --------
%   combine_mat_files('split_lists_frac_0.25_dir/resample_%d_data.mat','frac_0.25_data_resamples.mat',50)

traindatas = cell(N,1);
testdatas = cell(N,1);

for n = 1 : N
    disp(n)
    fname = sprintf(input_name_str,n);
    load(fname)
    traindatas{n} = traindata;
    testdatas{n} = testdata;
    try
        exptestdatas{n,1} = exptest;
        valdatas{n,1} = valdata;
    catch

    end
    clearvars *data
end

save(output_name,'*datas')
