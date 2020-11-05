% this function uses a sliding window of the energies to determine whether
% there is an "optimal" energy range where the balanced accuracy is greatest.
% As a bonus also gets the log-loss versus the accuracy and compares. Should
% think about whether use iterative or average. SHould I also see how the bacc
% results do on the logloss and vice versa? I think so...For now just look 
% at optimize for bacc and take both the logloss and the bacc

datasets = {'brcaMatchedTN','ERpn','kirckirp','luadlusc','lumAB'};
sig = @(x) 1./(1+exp(-x));
b1s = [0.03 0.1 0.3 1];
for n = 2 : length(datasets)
    d = datasets{n};
    disp(d)
    load(['~/Dropbox-Work/Wuxi/Data/' d '_bootstrap_resamples.mat'])
    base_dir = ['~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/' d '_bootstrap_resamples/'];
    sa_b1_loss = struct();
    for k = 1 : length(b1s)
        sa_sname = [d '_SA_pc44_nr_1000_nswps_1000_b0_0d1_b1_' num2str(b1s(k)) '_sols'];
        sa_sname = strrep(sa_sname,'.','d');
        load('~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/bootstrap_resamples_SA_sols',sa_sname);
    %% sa solutions come sorted
        eval(['sa_sols = ' sa_sname ';'])
        clearvars(sa_sname)
        for m = 1 : 100
            disp(m);
            trdata = traindatas{m};
            tstdata = testdatas{m};
            [h,J] = generateLogistic(trdata,'ising');
            N_SA_sols = min(cellfun(@length,sa_sols{m}));
            for ii = 1 : N_SA_sols-19
                tmp_sa_sols = cellfun(@(x) x(ii:ii+19,:),sa_sols{m},'uniformoutput',false);
                [sa_r,sa_t] = analyzeLogisticResults(tmp_sa_sols,trdata, tstdata, 'uniqueFlag', true, 'iterFlag', false, 'nSols', 20, 'lambdas',[0],'postProcess','test','biasFlag',false);
                sa_b1_loss(m,k).avg_bacc(ii) = sa_t.best_meansolsBacc;
                sa_b1_loss(m,k).avg_logloss(ii) = sa_t.meansolsLogloss;
                sa_b1_loss(m,k).avg_en(ii) = calcEns(mean(cell2mat(sa_r.TestItSols)),h,J);
                sa_b1_loss(m,k).avg_auc(ii) = sa_t.meansolsAUC;
            end
        end
    end
    try
        save([base_dir 'testmetrics_vs_slidingenergies.mat'],'-append','sa_b1_loss')
    catch
        save([base_dir 'testmetrics_vs_slidingenergies.mat'],'sa_b1_loss')
    end
end
