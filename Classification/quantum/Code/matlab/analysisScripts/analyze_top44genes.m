% Analysis script for bootstrap resamples. Bootstrap resamples refer to taking
% all the training data and repeatedly selecting 80% for training and 20% for 
% testing. Selected certain parameters for simplicity (see 
% analyze_lumAB_serial_dilutions.m for more parameters that can try).  

clear *results* *testperf*
% list of dataset names
datasets = {'brcaMatchedTN','ERpn','kirckirp','luadlusc','lumAB'};
sig = @(x) 1./(1+exp(-x));
for n = 1 : length(datasets)
    d = datasets{n};
    dir_name = ['~/Dropbox-Work/Wuxi/Results/top44genes/'];
    load(['~/Dropbox-Work/Wuxi/Data/' d '_top44genes_data.mat'])
    % main loop: go through resamples of data
    for p= 1 : 3
        tmp = rand(1000,size(traindata,2)-1);
        tmp(tmp<=0.5) = -1;
        tmp(tmp>0.5) = 1;
        rand_sols_10000{p,1} = tmp;
    end
    [rand_ens_10000,rand_idxs] = get_split_ens(rand_sols_10000,traindata);
    for p = 1 : 3 
        rand_sols_10000{p} = rand_sols_10000{p}(rand_idxs{p},:);
    end
%    [rand_results,rand_testperf] = analyzeLogisticResults(rand_sols,...
%        traindata, testdata, 'uniqueFlag', true, 'iterFlag', false, ...
%        'nSols', 20, 'lambdas',[0],'postProcess','test','biasFlag',false);
%    save([dir_name d '_results'],'rand_results','rand_testperf')
    save([dir_name d '_rand_out.mat'],'-append','rand_sols_10000','rand_ens_10000')
end
