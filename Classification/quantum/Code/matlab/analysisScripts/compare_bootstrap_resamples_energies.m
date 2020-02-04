% Wrapper file (kind of) for comparing the solutions returned by SA and DW
% for bootstrap resamples. Names of .mat files with results are hard-coded and
% can be changed as needed. 

base_dir = '~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/';
data_dir = '~/Dropbox-Work/Wuxi/Data/';
datasets = {'brcaMatchedTN','ERpn','kirckirp','luadlusc','lumAB'};

% SA results
load([base_dir 'bootstrap_resamples_SA_sols.mat'])

for n = 1 : length(datasets)
    d = datasets{n};
    disp(d)
    load([data_dir d '_bootstrap_resamples.mat']) % datas
    
    % DW results
    load([base_dir d '_bootstrap_resamples/cinit8.0_ri_at_5_nr_1000_out.mat']) 
    % SA solution that will use
    eval(['sa_sols = ' d '_SA_nr_1000_nswps_1000_b0_0d1_b1_3_out_sols;'])
    for m = 1 : 100
        dw_sols = cellfun(@(x) x{1},out{m},'uniformoutput',false);
        [dw_ens,sa_ens] = compare_energies(dw_sols(:,1),sa_sols{m}(:,1),traindatas{m});
        all_dw_ens{n}(m,:) = cellfun(@min,dw_ens)';
        all_sa_ens{n}(m,:) = cellfun(@min,sa_ens)';
    end
end
