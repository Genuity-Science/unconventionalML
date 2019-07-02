% Script to get SA bootstrap solutions from SA instance files.

datasets = {'brcaMatchedTN','ERpn','kirckirp','luadlusc','lumAB','lumAB_gene'};
%l_str = {'0','1d8','1d4','1d2','1','2','4','8'};
l_str = {'0'};
% Suffix for SA instance files. nr is number of repetitions, nswps is num. of 
%   sweeps, b0 is starting beta, b1 is final beta
%b1s = [0.03 0.1 0.3 1];
b1s = [0.03];
dir_name = '~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/';
for p = 1 : length(b1s)
    sfx = ['_nr-1000_nswps-1000_b0-0d01_b1-' num2str(b1s(p))];
    sfx = strrep(sfx,'.','d');
    sfx = [sfx '_out.dat'];

    for n = 1 : 1 %6 : length(datasets)    
        base_dir =[ '~/Dropbox-Work/Wuxi/SA/bootstrap_resamples/' datasets{n} '/'];
        for m = 1 : 100
            out = getSACVSols([base_dir 'Inst' num2str(m)],3,'.txt',sfx,l_str);
            sols{m} = {out.sols}';
        end
        save_name = [datasets{n} '_SA_pc44_' sfx '_sols'];
        save_name = strrep(save_name,'out.dat','');
        save_name = strrep(save_name,'-','_');
        save_name = strrep(save_name,'__','_');
        save_name = strrep(save_name,'.','d');
        eval([save_name '= sols;']); 
        clearvars sols 
    end
end
%save([dir_name 'bootstrap_resamples_SA_sols.mat'], '*sols','-append')
