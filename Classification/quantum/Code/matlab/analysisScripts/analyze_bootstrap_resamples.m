% Analysis script for bootstrap resamples. Bootstrap resamples refer to taking
% all the training data and repeatedly selecting 80% for training and 20% for 
% testing. Selected certain parameters for simplicity (see 
% analyze_lumAB_serial_dilutions.m for more parameters that can try).  

clear *results* *testperf*
% list of dataset names
datasets = {'brcaMatchedTN','ERpn','kirckirp','luadlusc','lumAB','lumAB_gene'};
sig = @(x) 1./(1+exp(-x));
nTotSols = [1 5 20 50 100 1000];
nSols = [1 5 10 20 50 100 1000];
cinits = [0.5 1 3 8 16 32];
for n = 6: length(datasets)
    d = datasets{n};
    disp(d)

    % define path to data directories and load 
    dir_name = ['~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/' d '_bootstrap_resamples/'];
    load(['~/Dropbox-Work/Wuxi/Data/' d '_bootstrap_resamples.mat'])

    % the line below can be changed based on the file using. In practice seemed
    % like cinit8.0 worked the best

    for k = 1 : length(cinits)
        cinit = cinits(k);
        load([dir_name 'cinit' num2str(cinit,'%1.1f') '_ri_at_5_nr_1000_out_2.mat'])
        for p = 1 : length(nTotSols);
            for ii = 1 : length(nSols)
                if nSols(ii) > nTotSols(p)
                    continue
                end
                % main loop: go through resamples of data
                rname = ['cinit' num2str(cinit) '_at_5_nsols_' num2str(nSols(ii)) ...
                    '_ntotsols_'  num2str(nTotSols(p)) '_results_simple'];
                tname =  ['cinit' num2str(cinit) '_at_5_nsols_' num2str(nSols(ii)) ...
                    '_ntotsols_' num2str(nTotSols(p)) '_testperf_simple'];
                rname = strrep(rname,'.','d');
                tname = strrep(tname,'.','d');
                s=sprintf('Cinit: %1.1f, nTotSols: %d, nSols: %d',cinit,nTotSols(p),nSols(ii));
                disp(s);
                for m = 1 : length(traindatas)
                    trdata = traindatas{m};
                    tstdata = testdatas{m};
                    sols = cellfun(@(x) x{1},out{m},'uniformoutput',false);
                    sols = cellfun(@(x) x(1:nTotSols(p),1:44),sols,'uniformoutput',false);
                    [results(m),testperf(m)] = analyzeLogisticResults(sols,...
                        trdata, tstdata, 'uniqueFlag', true, 'iterFlag', false, ...
                        'nSols', nSols(ii), 'lambdas',[0],'postProcess','test',...
                        'biasFlag',false,'metric','bacc');
                end
                eval([rname '= results;'])
                eval([tname '= testperf;'])
                try
                    save([dir_name 'results'],'-append',rname,tname)
                catch
                    save([dir_name 'results'],rname,tname);
                end
            end
            clearvars y_*
        end
    end
end
