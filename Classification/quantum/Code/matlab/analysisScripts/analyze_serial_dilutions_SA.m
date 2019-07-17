function [SA_results,SA_testperf] = analyze_serial_dilutions_SA(fracs,d,saveFlag,sfx)

% Analysis script for lumAB serial dilutions. Serial dilutions refer to 
% starting with a certain cut of the training data and successively decreasing
% the amount of training data. Selected certain parameters for simplicity (see 
% analyze_lumAB_serial_dilutions.m for more parameters that can try).  

sig = @(x) 1./(1+exp(-x));
%lambdas = [0 2.^[-3:3]];
%lambdas = [0 0.25 1 4];

% suffix for SA instance names
dir_name = ['~/Dropbox-Work/Wuxi/Results/' d '_splits/'];
nTotSols = 1000; %[1 5 20 50 100];
nSols = 20; %[1 5 10 20 50 100];
b1s = [0.03 ];
if strcmp(d,'6cancer')
    multiFlag = true;
else
    multiFlag = false;
end
for k = 1 : length(b1s)
    b1 = num2str(b1s(k));
    b1 = strrep(b1,'.','d');
    for p = 1 : length(nTotSols);
        for ii = 1 : length(nSols)
            if nSols(ii) > nTotSols(p)
                continue
            end
            s=sprintf('b1: %s, nTotSols: %d, nSols: %d',b1,nTotSols(p),nSols(ii));
            disp(s)
            for n = 1 : length(fracs)
                frac = num2str(fracs(n));
                fname = ['frac_' frac '_SA' sfx '.mat'];
                load([dir_name fname])
                disp(fname)
                % load data file
                %   traindata is cell array of different cuts of training data for a given
                %       fraction of the data.
                %   valdata is the original test split and remains unchanged between 
                %       different resampling(confusing, but that's how the collaborators 
                %       originally defined it)
                %   testdata is the portion of the original training cut left over after 
                %       taking a portion of the original training cut for training
                %   exptestdata is the original test data plus whatever fraction of data is
                %       not used for training in the training cut
                load(['~/Dropbox-Work/Wuxi/Data/' d '_splits/frac_' frac '_data_resamples.mat'])
                for m = 1 : length(traindatas)
                    disp(sprintf('Frac: %1.2f Iteration: %d',fracs(n),m))
                    trdata = traindatas{m};
                    tstdata = testdatas{m};
                    valdata = valdatas{m};
                    exptstdata = exptestdatas{m};
                    % to save time only analyzing first column, which is lambda=0
                    tmp_sols = cellfun(@(x) x(randperm(1000,nTotSols(p)),:),...
                        sols{m}(:,1),'uniformoutput',false);
                    if multiFlag 
                        [r,t] = analyzeMultinomialResults(tmp_sols,trdata,...
                           valdata, 'uniqueFlag', false, 'iterFlag', true, ...
                           'nSols', nSols(ii), 'lambdas',[0],'metric','acc',...
                           'biasFlag',false);
                    else
                    [r,t] = analyzeLogisticResults(tmp_sols, ...
                            trdata,valdata, 'uniqueFlag', false, 'iterFlag', true, ...
                            'nSols', nSols(ii), 'lambdas',[0],'postProcess','test',...
                            'biasFlag',false);
                    end
                    r.frac = fracs(n);
                    t.frac = fracs(n);
                    SA_results(n,m) = r;
                    SA_testperf(n,m) = t;
                end
            end
            s_sfx = strrep(sfx,'1000','1k');
            sar = ['SA' s_sfx '_nsols_' num2str(nSols(ii)) '_ntsols_' num2str(nTotSols(p)) '_results'];
            sat = ['SA' s_sfx '_nsols_' num2str(nSols(ii)) '_ntsols_' num2str(nTotSols(p)) '_testperf'];
            sar = strrep(sar,'-','_');
            sat = strrep(sat,'-','_');
            eval([sar '=SA_results;' sat '=SA_testperf;']);
            if saveFlag
                try 
                    save([dir_name 'SA_pca_results'],'-append',sar,sat)
                catch
                    save([dir_name 'SA_pca_results'],sar,sat)
                end
            end
        end
    end

    clearvars *datas y_*
end
