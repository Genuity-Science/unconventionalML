% this script uses a sliding window of the energies to determine whether
% there is an "optimal" energy range where the balanced accuracy is greatest.
% As a bonus also gets the log-loss versus the accuracy and compares. 

sig = @(x) 1./(1+exp(-x));
datasets = {'brcaMatchedTN','ERpn','kirckirp','luadlusc','lumAB'};
for n =  1 : length(datasets) 
    d = datasets{n};
    disp(d)
    load('~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/bootstrap_resamples_SA_sols',[d '*_3_sols']);
    % sa solutions come sorted
    eval(['sa_sols = ' d '_SA_pc44_nr_1000_nswps_1000_b0_0d1_b1_3_sols;'])
    base_dir = ['~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/' d '_bootstrap_resamples/'];
    load([base_dir 'cinit8.0_ri_at_5_nr_1000_out.mat']);
    load(['~/Dropbox-Work/Wuxi/Data/' d '_bootstrap_resamples']);
    dw_loss = struct();
    sa_loss = struct();
    rand_loss = struct();
    for m = 1 : 100
        disp(m);
        trdata = traindatas{m};
        tstdata = testdatas{m};
        [h,J] = generateLogistic(trdata,'ising');
        for p = 1 : 3
            tmp = rand(1000,44);
            tmp(tmp<=0.5) = -1;
            tmp(tmp>0.5) = 1;
            rand_sols{m}{p,1} = tmp;
        end
        dw_sols = cellfun(@(x) double(x{1}),out{m},'uniformoutput',false); 
        [rand_ens{m},rand_idx] = get_split_ens(rand_sols{m},trdata);
        [~,dw_idx] = get_split_ens(dw_sols,trdata);
        % sort solutions
        for p = 1 : 3
            dw_sols{p} = dw_sols{p}(dw_idx{p},:);
            rand_sols{m}{p} = rand_sols{m}{p}(rand_idx{p},:);
        end
        % go through solutions 20 at a time 
        for ii = 1 : 981
            tmp_dw_sols = cellfun(@(x) x(ii:ii+19,:),dw_sols,'uniformoutput',false);
            tmp_rand_sols = cellfun(@(x) x(ii:ii+19,:),rand_sols{m},'uniformoutput',false);
            [dw_r,dw_t] = analyzeLogisticResults(tmp_dw_sols, trdata,trdata, 'uniqueFlag', true, 'iterFlag', true, 'nSols', 20, 'lambdas',[0],'postProcess','test','biasFlag',false);
            [rand_r,rand_t] = analyzeLogisticResults(tmp_rand_sols, trdata,trdata, 'uniqueFlag', true, 'iterFlag', true, 'nSols', 20, 'lambdas',[0],'postProcess','test','biasFlag',false);
            dw_loss(m).it_bacc(ii) = dw_t.best_meansolsBacc;
            dw_loss(m).it_logloss(ii) = dw_t.meansolsLogloss;
            dw_loss(m).it_en(ii) = calcEns(mean(cell2mat(dw_r.TestItSols)),h,J);
            rand_loss(m).it_bacc(ii) = rand_t.best_meansolsBacc;
            rand_loss(m).it_logloss(ii) = rand_t.meansolsLogloss;
            rand_loss(m).it_en(ii) = calcEns(mean(cell2mat(rand_r.TestItSols)),h,J);
           
            [dw_r,dw_t] = analyzeLogisticResults(tmp_dw_sols, trdata,trdata, 'uniqueFlag', true, 'iterFlag', false, 'nSols', 20, 'lambdas',[0],'postProcess','test','biasFlag',false);
            [rand_r,rand_t] = analyzeLogisticResults(tmp_rand_sols,trdata, trdata, 'uniqueFlag', true, 'iterFlag', false, 'nSols', 20, 'lambdas',[0],'postProcess','test','biasFlag',false);
            dw_loss(m).avg_bacc(ii) = dw_t.best_meansolsBacc;
            dw_loss(m).avg_logloss(ii) = dw_t.meansolsLogloss;
            dw_loss(m).avg_en(ii) = calcEns(mean(cell2mat(dw_r.TestItSols)),h,J);
            dw_loss(m).avg_auc(ii) = dw_t.best_meansolsAUC;
            rand_loss(m).avg_bacc(ii) = rand_t.best_meansolsBacc;
            rand_loss(m).avg_logloss(ii) = rand_t.meansolsLogloss;
            rand_loss(m).avg_en(ii) = calcEns(mean(cell2mat(rand_r.TestItSols)),h,J);
            rand_loss(m).avg_auc(ii) = rand_t.meansolsAUC;
        end
        N_SA_sols = min(cellfun(@length,sa_sols{m}));
        for ii = 1 : N_SA_sols-19
            tmp_sa_sols = cellfun(@(x) x(ii:ii+19,:),sa_sols{m},'uniformoutput',false);
            [sa_r,sa_t] = analyzeLogisticResults(tmp_sa_sols,trdata, tstdata, 'uniqueFlag', true, 'iterFlag', true, 'nSols', 20, 'lambdas',[0],'postProcess','test','biasFlag',false);
            sa_loss(m).it_bacc(ii) = sa_t.best_meansolsBacc;
            sa_loss(m).it_logloss(ii) = sa_t.meansolsLogloss;
            sa_loss(m).it_en(ii) = calcEns(mean(cell2mat(sa_r.TestItSols)),h,J);
            [sa_r,sa_t] = analyzeLogisticResults(tmp_sa_sols,trdata, trdata, 'uniqueFlag', true, 'iterFlag', false, 'nSols', 20, 'lambdas',[0],'postProcess','test','biasFlag',false);
            sa_loss(m).avg_bacc(ii) = sa_t.best_meansolsBacc;
            sa_loss(m).avg_logloss(ii) = sa_t.meansolsLogloss;
            sa_loss(m).avg_en(ii) = calcEns(mean(cell2mat(sa_r.TestItSols)),h,J);
            sa_loss(m).avg_auc(ii) = sa_t.meansolsAUC;
        end
    end
    dw_train_loss = dw_loss;
    sa_train_loss = sa_loss;
    rand_train_loss = rand_loss;
    save([base_dir 'testmetrics_vs_slidingenergies.mat'],'-append','dw_train_loss','sa_train_loss','rand_train_loss'); 
end
