function [results,testperf] = analyze_bootstrap_resamples_SA(d,sols,saveFlag)
% Analysis script for SA bootstrap resamples. Bootstrap resamples refer to 
% taking all the training data and repeatedly selecting 80% for training and 
% 20% for testing. Selected certain parameters for simplicity (see 
% analyze_lumAB_serial_dilutions.m for more parameters that can try).  

%datasets = {'brcaMatchedTN','ERpn','kirckirp','luadlusc','lumAB','lumAB_gene'};
sig = @(x) 1./(1+exp(-x));

% string suffix that used for naming SA files

nTotSols = 1000; %[1 5 20 50 100 1000];
nSols = 20; %[1 5 10 20 50 100 1000];

%pcs = [24 44 84 117];
pcs = [44];
%b1s = [0.03 0.1 0.3 1 3];
b1s = [0.03];
%sa_mat_name = '~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/bootstrap_resamples_SA_sols.mat';
% previously saved solutions 
dir_name = ['~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/' d ...
            '_bootstrap_resamples/'];
load(['~/Dropbox-Work/Wuxi/Data/' d '_bootstrap_resamples.mat'])
all_traindatas = traindatas;
all_testdatas = testdatas;
% eval(['sols = ' d '_SA_sols;']);
for k = 1 : length(pcs)
    traindatas = cellfun(@(x) x(:,1:pcs(k)+1),all_traindatas,'uniformoutput',false);
    testdatas = cellfun(@(x) x(:,1:pcs(k)+1),all_testdatas,'uniformoutput',false);
    for ii = 1 : length(b1s)
        b1 = num2str(b1s(ii));
        for ij = 1 : length(nTotSols)
            for p = 1 : length(nSols)
                if nSols(p) > nTotSols(ij)
                    continue
                end
                rname = ['SA_nr_1000_nswps_1000_pc' num2str(pcs(k)) '_nsols_' ...
                    num2str(nSols(p)) '_nts_' num2str(nTotSols(ij)) '_b1_' b1 '_results'];
                tname =  ['SA_nr_1000_nswps_1000_pc' num2str(pcs(k)) '_nsols_' ...
                    num2str(nSols(p)) '_nts_' num2str(nTotSols(ij)) '_b1_' b1 '_testperf'];
                rname = strrep(rname,'.','d');
                tname = strrep(tname,'.','d');
                s=sprintf('b1: %s, nTotSols: %d, nSols: %d',b1,nTotSols(ij),nSols(p));
                disp(s);
                for m = 1 : length(traindatas)
                    max_N = min(1000,min(cellfun(@(x) size(x,1),sols{m})));
                    tmpsols = cellfun(@(x) x(randperm(max_N),:),sols{m},'uniformoutput',false);
                    trdata = traindatas{m};
                    tstdata = testdatas{m};
                    tmp_sols = cellfun(@(x) x(1:min(nTotSols(ij),max_N),:),sols{m},'uniformoutput',false);
                    [results(m),testperf(m)] = ...
                        analyzeLogisticResults(tmp_sols,trdata, tstdata, 'uniqueFlag', true,... 
                        'iterFlag', false, 'nSols',nSols(p), 'lambdas',[0],'postProcess','test',...
                        'biasFlag',false);
            
                    % output actual classes and predicted classes for plotting with R
                    y_trains(:,m) = trdata(:,1);
                    y_tests(:,m) = tstdata(:,1);
                    tmpsol = mean(cell2mat(results(m).TestItSols));
                    y_pred_trains(:,m) = sig(trdata(:,2:end)*tmpsol');
                    y_pred_tests(:,m) = sig(tstdata(:,2:end)*tmpsol');
                end
        %        save([dir_name 'SA_pred_for_R2.mat'],'y_trains','y_tests',...
        %             'y_pred_tests','y_pred_trains')            
                eval([rname '= results;'])
                eval([tname '= testperf;'])
                if saveFlag 
                    try
                        save([dir_name 'results'],'-append',rname,tname)
                    catch
                        save([dir_name 'results'],rname,tname)
                    end
                end
            end
        end
    end
end
