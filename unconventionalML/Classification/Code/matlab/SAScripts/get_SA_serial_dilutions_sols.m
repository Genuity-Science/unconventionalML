function get_SA_serial_dilutions_sols(fracs,N,base_dir,dir_name,outPfx,l_str)
%   Script to get lumAB serial dilutions from SA instance files.
%   get_SA_serial_dilutions_sols(fracs,N,base_dir,dir_name,outPfx,l_str)
%   base_dir = '~/Dropbox-Work/Wuxi/SA/luadlusc_splits' the directory where SA files are
%   dir_name: the directory where save the SA mat file
%       e.g.: '~/Dropbox-Work/Wuxi/Results/luadlusc_splits' 

if nargin < 5
    outPfx = '_nr-1000_nswps-1000_b0-0d01_b1-0d03_out.dat';
end
if nargin < 6
    l_str = {'0'};
end
for n = 1:length(fracs)
    frac = num2str(fracs(n));
    disp(frac)
    for m = 1 : N
        filestem = [base_dir 'frac_' frac '_Inst' num2str(m)];
        inSfx = '.txt';
        out = getSACVSols(filestem,3,inSfx,outPfx,l_str);
        sols{m} = reshape({out.sols}',3,[]);
    end
    save_name = ['frac_' frac '_SA' outPfx];
    save_name = strrep(save_name,'_out.dat','.mat');
    disp(save_name)
    save([dir_name save_name],'sols')
    clearvars sols 
end
