% Wrapper file (kind of) for comparing the solutions returned by SA and DW
% for lumAB serial dilutions

base_dir = '~/Dropbox-Work/Wuxi/Results/lumAB_splits/';
data_dir = '~/Dropbox-Work/Wuxi/Data/lumAB_splits/';
fracs = [0.18 0.2 0.3:0.05:0.95];
sfx = '_nr-1000_nswps-1000_b0-0d1_b1-3';
load_names = {'cinit1_ri_sols.mat','cinit8.0_ri_sols.mat','cinit8.0_ri_sols.mat','cinit8.0_ri_at_5_nr_1000_out.mat','cinit1_ri_sols.mat'};

for n = 1 : length(fracs)
    frac = num2str(fracs(n));
    disp(frac)
    fname = ['frac_' frac '_SA' sfx '.mat'];
    load([data_dir 'frac_' frac '_data_resamples.mat'])
    load([base_dir fname])
    sa_sols = sols;
    load([base_dir  'lumAB_splits_frac_' frac '_cinit_3_ri_out_at_5_nr_1000.mat'])
    for m = 1 : 50
        dw_sols = cellfun(@(x) x{1},out{n},'uniformoutput',false);
        [dw_ens,sa_ens] = compare_energies(dw_sols(:,1),sa_sols{m}(:,1),traindatas{m});
        all_dw_ens{n}(m,:) = cellfun(@min,dw_ens)';
        all_sa_ens{n}(m,:) = cellfun(@min,sa_ens)';
    end
end
