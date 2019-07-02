% Analysis script for lumAB serial dilutions. Serial dilutions refer to 
% starting with a certain cut of the training data and successively decreasing
% the amount of training data. Selected certain parameters for simplicity (see 
% analyze_lumAB_serial_dilutions.m for more parameters that can try).  

clearvars *results* *testperf* 

% the fraction of the original training data to use
fracs = [0.18 0.2:0.05:0.95];
nTotSols =1000; %[5 50]; % [20 100 1000];
nSols = 20;%[1 5 10 20 50];% 5 10 20 100 1000];
for p = 1 : length(nTotSols);
    for ii = 1 : length(nSols)
        if nSols(ii) > nTotSols(p)
            continue
        end
        % main loop: go through resamples of data
        s=sprintf('nTotSols: %d, nSols: %d',nTotSols(p),nSols(ii));
        disp(s);
        for n = 1 : length(fracs)
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
            load(['~/Dropbox-Work/Wuxi/Data/lumAB_splits/frac_' ...
                    num2str(fracs(n)) '_data_resamples.mat']) 
        
%            load('~/Dropbox-Work/Wuxi/Results/lumAB_splits/rand_sols','rand_pc44_nsols5_sols');
%            rand_sols = rand_pc44_nsols5_sols;
            % load result file, output in variable named out
            % out is a cell array of cells. Each cell is 1x3, with the first
            % cell being the solutions, the second the coupling strength, and the 
            % third the fraction of unbroken solutions. See runDW() in DW_utils.py
            
            for m = 1 : length(traindatas2)
                for k = 1 : 3
                    tmp = rand(1000,44);
                    tmp(tmp<=0.5) = -1;
                    tmp(tmp>0.5) = 1;
                    rand_sols{m,n}{k,1} = tmp;
                end
                disp(sprintf('Frac: %1.2f Iteration: %d',fracs(n),m))
                tmp_rand_sols = cellfun(@(x) x(1:nTotSols(p),:), rand_sols{m,n},'uniformoutput',false);
                [r,t] = analyzeLogisticResults(tmp_rand_sols,traindatas2{m},valdatas2{m},...
                        'nSols',nSols(ii),'uniqueFlag',true,'lambdas',[0],'iterFlag',true,...
                        'postProcess','test','biasFlag',false);
                r.frac = fracs(n);
                t.frac = fracs(n);
                results(n,m) = r;
                testperf(n,m) = t;
            end
            clear *datas*
        end

        rand_pc44_ntotsols_20_sols = rand_sols;
        save('~/Dropbox-Work/Wuxi/Results/lumAB_splits/rand_sols','rand_pc44_ntotsols_20_sols')   

        rname = ['rand_pc44_nsols_' num2str(nSols(ii)) '_ntotsols' num2str(nTotSols(p)) '_results'];
        tname = ['rand_pc44_nsols_' num2str(nSols(ii)) '_ntotsols' num2str(nTotSols(p)) '_testperf'];
        eval([rname '= results;'])
        eval([tname '=testperf;']);
    save('~/Dropbox-Work/Wuxi/Results/lumAB_splits/rand_pca_results',rname,tname,'-append')
    end
end

