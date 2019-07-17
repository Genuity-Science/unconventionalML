% Output serial dilution predictions for R. In case changed or 
% want specific result. See analyze_serial_dilutions.m
function output_serial_dilutions_preds_for_R(d,r,fname);
% E.g: output_serial_dilutions_preds_for_R('lumAB',rand_pc44_nsols20_results,'~/Dropbox-Work/Wuxi/Results/lumAB_splits/preds_for_R/rand_frac_%s_nsols20_preds_for_R.mat')

sig = @(x) 1./(1+exp(-x));

for n = 1 : size(r,1)
    frac = r(n,1).frac;
    load(['~/Dropbox-Work/Wuxi/Data/' d '_splits/frac_' num2str(frac) '_data_resamples.mat']);
    for m = 1 : length(traindatas) 
        y_trains(:,m) = traindatas{m}(:,1);
        y_exptests(:,m) = exptestdatas{m}(:,1);
        y_valids(:,m) = valdatas{m}(:,1);
        y_tests(:,m) = testdatas{m}(:,1);

        meansol = mean(cell2mat(r(n,m).TestItSols));
        if strcmp(d,'6cancer')
            [~,~,y_pred_trains{m,1}] = getMultinomialAcc(meansol,traindatas{m},0:5);
            [~,~,y_pred_tests{m,1}] = getMultinomialAcc(meansol,testdatas{m},0:5);
            [~,~,y_pred_exptests{m,1}] = getMultinomialAcc(meansol,exptestdatas{m},0:5);
            [~,~,y_pred_vals{m,1}] = getMultinomialAcc(meansol,valdatas{m},0:5);
        else
            y_pred_trains(:,m) = sig(traindatas{m}(:,2:end)*meansol');
            y_pred_exptests(:,m) = sig(exptestdatas{m}(:,2:end)*meansol');
            y_pred_tests(:,m) = sig(testdatas{m}(:,2:end)*meansol');
            y_pred_vals(:,m) = sig(valdatas{m}(:,2:end)*meansol');
        end
    end
    save(sprintf(fname,num2str(frac)),'y_*')
    clearvars y_*
end
