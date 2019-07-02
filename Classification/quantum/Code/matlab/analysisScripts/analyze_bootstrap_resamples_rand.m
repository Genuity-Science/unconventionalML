% Analysis script for bootstrap resamples. Bootstrap resamples refer to taking
% all the training data and repeatedly selecting 80% for training and 20% for 
% testing. Selected certain parameters for simplicity (see 
% analyze_lumAB_serial_dilutions.m for more parameters that can try).  

clear *results* *testperf*
% list of dataset names
datasets = {'brcaMatchedTN','ERpn','kirckirp','luadlusc','lumAB','lumAB_gene'};
sig = @(x) 1./(1+exp(-x));
% total number of solutions to try
nTotSols = [1 5 20 50 100 1000];
% number of solutions to include in the iterative procedure
nSols = [1 5 10 20 50 100 1000];
%nSols = [1000];
%pcs = [24 44 65 84 117];
pcs = [44];
for n = 6 : length(datasets)
    d = datasets{n};
    disp(d)
    dir_name = ['~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/' d ...
        '_bootstrap_resamples/'];
    load(['~/Dropbox-Work/Wuxi/Data/' d '_bootstrap_resamples.mat'])
    all_traindatas = traindatas2;
    all_testdatas = testdatas2; 
    for ii = 1 : length(pcs)
        traindatas = cellfun(@(x) x(:,1:pcs(ii)+1),all_traindatas,'uniformoutput',false);
        testdatas = cellfun(@(x) x(:,1:pcs(ii)+1),all_testdatas,'uniformoutput',false);
        sname = ['rand_sols_pc' num2str(pcs(ii))];
        try
            load([dir_name 'rand_out'],sname);
            eval(['sols =' sname ';']);
        catch
            % uncomment below if want to regenerate random solutions 
            disp('Generating random solutions')
            for m = 1 : 100
                for p= 1 : 3
                    tmp = rand(1000,pcs(ii));
                    tmp(tmp<=0.5) = -1;
                    tmp(tmp>0.5) = 1;
                    sols{m}{p,1} = tmp;
                end
            end
            eval([sname '=sols;']);
            try
                save([dir_name 'rand_out.mat'],'-append',sname)
            catch
                save([dir_name 'rand_out.mat'],sname)
            end
        end
        % main loop: go through resamples of data
        for k = 1 : length(nTotSols)
            for ij = 1 : length(nSols)

                if nSols(ij) > nTotSols(k)
                    continue
                end
                s=sprintf('nTotSols: %d, nSols: %d',nTotSols(k),nSols(ij));
                disp(s);
                for m = 1 : length(traindatas)
                    trdata = traindatas{m};
                    tstdata = testdatas{m};
                    s = cellfun(@(x) x(1:nTotSols(k),:),sols{m},'uniformoutput',false);
                    [r(m),t(m)] = analyzeLogisticResults(s,...
                        trdata, tstdata, 'uniqueFlag', true, 'iterFlag', false, ...
                        'nSols', nSols(ij), 'lambdas',[0],'postProcess','test','biasFlag',false);
                end
                rname = ['rand_pc' num2str(pcs(ii)) '_nsols_' num2str(nSols(ij)) ...
                    '_nTotSols_' num2str(nTotSols(k)) '_results'];
                tname = ['rand_pc' num2str(pcs(ii)) '_nsols_' num2str(nSols(ij)) ...
                    '_nTotSols_' num2str(nTotSols(k)) '_testperf'];
                eval([rname '=r;' tname '=t;'])
                save([dir_name 'results'],'-append',rname,tname)
            end
        end
    end
end
