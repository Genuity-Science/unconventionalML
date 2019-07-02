% Analysis script for bootstrap resamples. Bootstrap resamples refer to taking
% all the training data and repeatedly selecting 80% for training and 20% for 
% testing. Selected certain parameters for simplicity (see 
% analyze_lumAB_serial_dilutions.m for more parameters that can try).  

clear *results* *testperf*
% list of dataset names
datasets = {'brcaMatchedTN','ERpn','kirckirp','luadlusc','lumAB'};
sig = @(x) 1./(1+exp(-x));
for n = 4 : length(datasets)
    d = datasets{n};
    dir_name = ['~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/' d ...
        '_bootstrap_resamples/'];
    load(['~/Dropbox-Work/Wuxi/Data/' d '_pc65_bootstrap_resamples.mat'])
    % main loop: go through resamples of data
    for m = 1 : length(traindatas)
        trdata = traindatas{m};
        tstdata = testdatas{m};
        for p= 1 : 3
            tmp = rand(10000,65);
            tmp(tmp<=0.5) = -1;
            tmp(tmp>0.5) = 1;
            sols{p,1} = tmp;
        end
        [rand_ens_pc65_10000{m},rand_idx] = get_split_ens(sols,traindatas{m});
         % don't sort so that can take simulate taking smaller number of reads
%        for p = 1 : 3 
%            sols{p} = sols{p}(rand_idx{p},:); 
%        end
        rand_sols_pc65_10000{m} = sols;
%        [rand_results(m),rand_testperf(m)] = analyzeLogisticResults(sols,...
%            trdata, tstdata, 'uniqueFlag', true, 'iterFlag', true, ...
%            'nSols', 20, 'lambdas',[0],'postProcess','test','biasFlag',false);
%        save([dir_name 'results'],'-append','rand_results','rand_testperf')
        save([dir_name 'rand_out.mat'],'-append','rand_sols_pc65_10000','rand_ens_pc65_10000')
    end
end
