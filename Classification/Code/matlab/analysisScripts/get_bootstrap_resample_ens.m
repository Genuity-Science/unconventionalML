
datasets = {'brcaMatchedTN','ERpn','kirckirp','luadlusc','lumAB'};

for n = 1:5
    d = datasets{n};
    disp(d)

    % define path to data directories and load 
    dir_name = ['~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/' d '_bootstrap_resamples/'];
    load(['~/Dropbox-Work/Wuxi/Data/' d '_bootstrap_resamples.mat'])
    load([dir_name 'cinit8.0_ri_at_5_nr_1000_out.mat']);
    saname = [d '_SA_nr_1000_nswps_1000_b0_0d1_b1_3_out_sols'];
    load('~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/bootstrap_resamples_SA_sols',saname)
    eval(['sa_sols=' saname ';']);
    for m = 1 : 100
        sols = cellfun(@(x) x{1},out{m},'uniformoutput',false); 
        dw_ens{m} = get_split_ens(sols,traindatas{m});
        sa_ens{m} = get_split_ens(sa_sols{m},traindatas{m});
    end
    save([dir_name 'ens'], 'dw_ens', 'sa_ens');
end

