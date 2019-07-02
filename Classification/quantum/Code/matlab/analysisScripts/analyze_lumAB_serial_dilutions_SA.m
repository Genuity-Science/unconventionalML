% Analysis script for lumAB serial dilutions. Serial dilutions refer to 
% starting with a certain cut of the training data and successively decreasing
% the amount of training data. Selected certain parameters for simplicity (see 
% analyze_lumAB_serial_dilutions.m for more parameters that can try).  

sig = @(x) 1./(1+exp(-x));
%fracs = [0.18 0.2 0.3:0.05:0.95];
fracs = 0.25;
%lambdas = [0 2.^[-3:3]];
%lambdas = [0 0.25 1 4];

% suffix for SA instance names
dir_name = '~/Dropbox-Work/Wuxi/Results/lumAB_splits/';
%nTotSols = [1 5 20 50 100];
%nSols = [1 5 10 20 50 100];
nTotSols=1000;
nSols = 20;
b1s = [0.03 ];
for k = 1 : length(b1s)
    b1 = num2str(b1s(k));
    b1 = strrep(b1,'.','d');
    if k == 1
        sfx = ['_nr-1000_nswps-1000_b0-0d01_b1-' b1];
    else
        sfx2 = ['_nr-1000_nswps-1000_b0-0d01_b1-' b1];
        sfx = ['_nr-1000_nswps-1000_b0-0d1_b1-' b1];
    end
    for p = 1 : length(nTotSols);
        for ii = 1 : length(nSols)
            if nSols(ii) > nTotSols(p)
                continue
            end
            s=sprintf('b1: %s, nTotSols: %d, nSols: %d',b1,nTotSols(p),nSols(ii));
            disp(s)
            for n = 1 : length(fracs)
                frac = num2str(fracs(n));
                try
                    fname = ['frac_' frac '_SA' sfx '.mat'];
                    load([dir_name fname])
                catch
                    fname = ['frac_' frac '_SA' sfx2 '.mat'];
                    sfx = sfx2;
                    load([dir_name fname])
                end
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
                load(['~/Dropbox-Work/Wuxi/Data/lumAB_splits/frac_' frac '_data_resamples.mat'])
                for m = 1 : length(traindatas)
                    trdata = traindatas2{m};
                    tstdata = testdatas2{m};
                    valdata = valdatas2{m};
                    exptstdata = exptestdatas2{m};
                    % to save time only analyzing first column, which is lambda=0
                    tmp_sols = cellfun(@(x) x(randperm(1000,nTotSols(p)),:),sols{m}(:,1),'uniformoutput',false);
                    [SA_results(m,n),SA_testperf(m,n)] = analyzeLogisticResults(tmp_sols, ...
                            trdata,valdata, 'uniqueFlag', false, 'iterFlag', true, ...
                            'nSols', nSols(ii), 'lambdas',[0],'postProcess','test',...
                            'biasFlag',false);
                end
            end
            s_sfx = strrep(sfx,'1000','1k');
            sar = ['SA' s_sfx '_nsols_' num2str(nSols(ii)) '_ntsols_' num2str(nTotSols(p)) '_results'];
            sat = ['SA' s_sfx '_nsols_' num2str(nSols(ii)) '_ntsols_' num2str(nTotSols(p)) '_testperf'];
            sar = strrep(sar,'-','_');
            sat = strrep(sat,'-','_');
            eval([sar '=SA_results;' sat '=SA_testperf;']);
            save([dir_name 'pca_results'],'-append',sar,sat)
        end
    end

    % save predictions that will read and plot using R 
%    save([dir_name 'preds_for_R/SA_frac_' frac '_pred_for_R2.mat'],'y_*');
    clearvars *datas y_*
end
