% Script to analyze solutions for field method. 

datasets = {'brcaMatchedTN','ERpn','kirckirp','luadlusc','lumAB','lumAB_gene'};
sig = @(x) 1./(1+exp(-x));
for n = 6: 6
    d = datasets{n};
    disp(d)
    dir_name = ['~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/' d ...
                '_bootstrap_resamples/'];

    % load data             
    try
        load(['~/Dropbox-Work/Wuxi/Data/' d '_pc118_bootstrap_resamples.mat'])
    catch
        load(['~/Dropbox-Work/Wuxi/Data/' d '_bootstrap_resamples.mat'])
        traindatas = traindatas2;
        testdatas = testdatas2;
    end

    train_splits = get_split_idxs(size(traindatas{1},1),3);
    % loop over 100 resamples 
    for m = 1 : 100
        trdata = traindatas{m}(:,1:45);
        tstdata = testdatas{m}(:,1:45);
        for k = 1 : 3
            [h,J] = generateLogistic(trdata(train_splits{k},:));
            cl_sols{m}{k,1} = -sign(h)'; 
        end
        [field_r(m),field_t(m)] = analyzeLogisticResults(cl_sols{m},trdata, tstdata, 'uniqueFlag', true,... 
                    'iterFlag', false, 'nSols', 1, 'lambdas',[0],'postProcess','test',...
                    'biasFlag',false);
        y_trains(:,m) = trdata(:,1);
        y_tests(:,m) = tstdata(:,1);
        tmpsol = mean(cell2mat(cl_r(m).TestItSols(:,1)));
        y_pred_trains(:,m) = sig(trdata(:,2:end)*tmpsol');
        y_pred_tests(:,m) = sig(tstdata(:,2:end)*tmpsol');
    end
    save([dir_name 'results'],'-append','field_r','field_t')
    save([dir_name 'field_pred_for_R.mat'],'y_trains','y_tests','y_pred_tests','y_pred_trains');
    clear y_*
end
