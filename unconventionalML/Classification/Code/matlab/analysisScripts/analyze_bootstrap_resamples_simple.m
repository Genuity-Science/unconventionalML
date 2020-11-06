function [results,testperf] = analyze_bootstrap_resamples(d,out,saveFlag,cinit)
% Analysis script for bootstrap resamples. Bootstrap resamples refer to taking
% all the training data and repeatedly selecting 80% for training and 20% for 
% testing. Selected certain parameters for simplicity (see 
% analyze_lumAB_serial_dilutions.m for more parameters that can try).  
if nargin < 4
    cinit = 8;
end
sig = @(x) 1./(1+exp(-x));
%nSols = [20 50 100 1000];
%nTotSols = [1 5 20 50 100 1000];
%nSols = [1 5 10 20 50 100 1000];
%cinits = [0.5 1 3 8 16 32];
nTotSols = 1000;
nSols = 20;
%cinits = [8];
%cinits = [16 32];
% define path to data directories and load 
dir_name = ['~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/' d '_bootstrap_resamples/'];
load(['~/Dropbox-Work/Wuxi/Data/' d '_bootstrap_resamples.mat'])
pcs = 44;
% the line below can be changed based on the file using. In practice seemed
% like cinit8.0 worked the best

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
            trdata = traindatas{m}(:,1:pcs+1);
            tstdata = testdatas{m}(:,1:pcs+1);
            sols = cellfun(@(x) x{1},out{m},'uniformoutput',false);
            sols = cellfun(@(x) x(1:nTotSols(p),1:pcs),sols,'uniformoutput',false);
            [results(m),testperf(m)] = analyzeLogisticResults(sols,...
                trdata, tstdata, 'uniqueFlag', true, 'iterFlag', false, ...
                'nSols', nSols(ii), 'lambdas',[0],'postProcess','test',...
                'biasFlag',false,'metric','bacc');
             y_trains(:,m) = trdata(:,1);
             y_tests(:,m) = tstdata(:,1);
             tmpsol = mean(cell2mat(results(m).TestItSols(:,1)));
             y_pred_trains(:,m) = sig(trdata(:,2:end)*tmpsol');
             y_pred_tests(:,m) = sig(tstdata(:,2:end)*tmpsol');
        end
        eval([rname '= results;'])
        eval([tname '= testperf;'])
        if saveFlag
            try
                save([dir_name 'results'],'-append',rname,tname)
            catch
                save([dir_name 'results'],rname,tname);
            end
        end
    end
    % save in a format for plotting with R
    if saveFlag
        save([dir_name 'pred_for_R.mat'],'y_trains','y_tests','y_pred_tests','y_pred_trains');
    end
end


