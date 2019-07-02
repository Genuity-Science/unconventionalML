% Find a variety of test metrics when include an increasing number of
% solutions. Not sure why I'm doing this? Anyway, Daniel thinks it's
% relevant for some reason.

datasets = {'brcaMatchedTN','ERpn','kirckirp','luadlusc','lumAB'};
sig = @(x) 1./(1+exp(-x));
for n = 1 : length(datasets)
    d = datasets{n};
    disp(d)
    % sa solutions come sorted
    sa_sols_name = [d '_SA_pc44_nr_1000_nswps_1000_b0_0d1_b1_3_sols'];
    load('~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/bootstrap_resamples_SA_sols',sa_sols_name);
    eval(['sa_sols = ' sa_sols_name ';']);
    base_dir = ['~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/' d '_bootstrap_resamples/'];
    load([base_dir 'cinit8.0_ri_at_5_nr_1000_out.mat']);
    load(['~/Dropbox-Work/Wuxi/Data/' d '_bootstrap_resamples']);
    dw_loss = struct();
    sa_loss = struct();
    rand_loss = struct();
    load([base_dir 'rand_out'],'rand_sols_pc44');
    rand_sols = rand_sols_pc44;
    nSols = [1 10 20 50 100 500 1000];
    for m = 1 : 100
        disp(m);
        trdata = traindatas{m};
        tstdata = testdatas{m};
        [h,J] = generateLogistic(trdata,'ising');
        [h_t,J_t] = generateLogistic(tstdata,'ising');
        dw_sols = cellfun(@(x) double(x{1}),out{m},'uniformoutput',false);
        [~,dw_idx] = get_split_ens(dw_sols,trdata);
        % sort solutions
        [rand_ens,rand_idx] = get_split_ens(rand_sols{m},trdata);
        for p = 1 : 3
            dw_sols{p} = dw_sols{p}(dw_idx{p},:);
            tmp_rand_sols{p,1} = rand_sols{m}{p}(rand_idx{p},:);
        end
%        [dw_r,dw_t] = analyzeLogisticResults(dw_sols, trdata,tstdata, 'uniqueFlag', true, 'iterFlag', true, 'nSols', nSols(ii), 'lambdas',[0],'postProcess','test','biasFlag',false,'metric','bacc');
%        [rand_r,rand_t] = analyzeLogisticResults(rand_sols{m}, trdata,tstdata, 'uniqueFlag', true, 'iterFlag', true, 'nSols', nSols(ii), 'lambdas',[0],'postProcess','test','biasFlag',false,'metric','bacc');
%        dw_loss(m).it_bacc(ii) = dw_t.best_meansolsBacc;
%        dw_loss(m).it_auc(ii) = dw_t.best_meansolsAUC;
%        dw_loss(m).it_logloss(ii) = dw_t.meansolsLogloss;
%        dw_loss(m).it_en(ii) = calcEns(mean(cell2mat(dw_r.TestItSols)),h,J);
%        rand_loss(m).it_bacc(ii) = rand_t.best_meansolsBacc;
%        rand_loss(m).it_logloss(ii) = rand_t.meansolsLogloss;
%        rand_loss(m).it_en(ii) = calcEns(mean(cell2mat(rand_r.TestItSols)),h,J);
%        rand_loss(m).it_auc(ii) = rand_t.meansolsAUC;
       
        [dw_r,dw_t] = analyzeLogisticResults(dw_sols, trdata,tstdata, 'uniqueFlag', true, 'iterFlag', false, 'nSols', nSols(end), 'lambdas',[0],'postProcess','test','biasFlag',false);
        [rand_r,rand_t] = analyzeLogisticResults(tmp_rand_sols,trdata, tstdata, 'uniqueFlag', true, 'iterFlag', false, 'nSols', nSols(end), 'lambdas',[0],'postProcess','test','biasFlag',false);
        [sa_r,sa_t] = analyzeLogisticResults(sa_sols{m},trdata, tstdata, 'uniqueFlag', true, 'iterFlag', false, 'nSols', nSols(end), 'lambdas',[0],'postProcess','test','biasFlag',false);
        for ii = 1 : length(nSols)
            for p = 1 : 3
                [~,idx] = min(dw_r.TestItObjf{p}(1:nSols(ii)));
                dw_tmpsols{p,1} = mean(dw_sols{p}(1:idx,:),1);
                [~,idx] = min(sa_r.TestItObjf{p}(1:min(nSols(ii),sa_r.TestItObjf{p})));
                sa_tmpsols{p,1} = mean(sa_sols{m}{p}(1:idx,:),1);
                [~,idx] = min(rand_r.TestItObjf{p}(1:nSols(ii)));
                rand_tmpsols{p,1} = mean(tmp_rand_sols{p}(1:idx,:),1);
            end
            dw_t = getPerf(dw_tmpsols,tstdata,0,sig,true);
            dw_t.meansolsLogloss = calcObjF(mean(cell2mat(dw_tmpsols)),tstdata,sig,'log');
            dw_loss(m).avg_bacc(ii) = dw_t.meansolsBacc;
            dw_loss(m).avg_auc(ii) = dw_t.meansolsAUC;
            dw_loss(m).avg_logloss(ii) = dw_t.meansolsLogloss;
%            dw_loss(m).avg_en(ii) = calcEns(mean(cell2mat(dw_r.TestItSols)),h,J);
            dw_loss(m).avg_en(ii) = calcEns(mean(cell2mat(dw_tmpsols)),h_t,J_t);
            rand_t = getPerf(rand_tmpsols,tstdata,0,sig,true);
            rand_t.meansolsLogloss = calcObjF(mean(cell2mat(rand_tmpsols)),tstdata,sig,'log');
            rand_loss(m).avg_bacc(ii) = rand_t.meansolsBacc;
            rand_loss(m).avg_logloss(ii) = rand_t.meansolsLogloss;
%            rand_loss(m).avg_en(ii) = calcEns(mean(cell2mat(rand_r.TestItSols)),h,J);
            rand_loss(m).avg_en(ii) = calcEns(mean(cell2mat(rand_tmpsols)),h_t,J_t);
            rand_loss(m).avg_auc(ii) = rand_t.meansolsAUC;
            
    %        [sa_r,sa_t] = analyzeLogisticResults(sa_sols{m},trdata, tstdata, 'uniqueFlag', true, 'iterFlag', true, 'nSols', nSols(ii), 'lambdas',[0],'postProcess','test','biasFlag',false);
    %        sa_loss(m).it_bacc(ii) = sa_t.best_meansolsBacc;
    %        sa_loss(m).it_auc(ii) = sa_t.best_meansolsAUC;
    %        sa_loss(m).it_logloss(ii) = sa_t.meansolsLogloss;
    %        sa_loss(m).it_en(ii) = calcEns(mean(cell2mat(sa_r.TestItSols)),h,J);
            sa_t = getPerf(sa_tmpsols,tstdata,0,sig,true);
            sa_t.meansolsLogloss = calcObjF(mean(cell2mat(sa_tmpsols)),tstdata,sig,'log');
            sa_loss(m).avg_bacc(ii) = sa_t.meansolsBacc;
            sa_loss(m).avg_logloss(ii) = sa_t.meansolsLogloss;
%            sa_loss(m).avg_en(ii) = calcEns(mean(cell2mat(sa_r.TestItSols)),h,J);
            sa_loss(m).avg_en(ii) = calcEns(mean(cell2mat(sa_tmpsols)),h_t,J_t);
            sa_loss(m).avg_auc(ii) = sa_t.meansolsAUC;
        end
    end
    save([base_dir 'testmetrics_vsnsols.mat'],'dw_loss','sa_loss','rand_loss'); 
%    save([base_dir 'testmetrics_vsnsols.mat'],'-append','rand_loss')
end
