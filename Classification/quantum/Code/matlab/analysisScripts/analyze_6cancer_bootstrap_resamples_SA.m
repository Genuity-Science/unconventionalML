function [r,t] = analyze_multinomial_bootstrap_resamples_SA(sols,saveFlag)

% Analysis script for multinomial. Try many different evaluation metrics (auc, 
% balanced accuracy, log loss),ways of combining solutions (iterative, 
% average), number of solutions to combine (20,100). See also
% analyze_lumAB_serial_dilutions_simple.m for a script with one parameter
% selected, based on trying multiple 
load('~/Dropbox-Work/Wuxi/Data/6cancer_bootstrap_resamples.mat')
for m = 1 : length(sols)
    tmp_sols = cellfun(@(x) x(1:1000,:),sols{m},'uniformoutput',false);
    [r(m),t(m)] = analyzeMultinomialResults(tmp_sols,traindatas{m},testdatas{m},...
        'metric','acc','nSols',20,'uniqueFlag',true,'lambdas',0,...
        'iterFlag',true);
    [~,~,y_pred_trains{m}] = getMultinomialAcc(mean(cell2mat(r(m).TestItSols)),traindatas{m});
    [~,~,y_pred_tests{m}] = getMultinomialAcc(mean(cell2mat(r(m).TestItSols)),testdatas{m});
    y_trains(:,m) = traindatas{m}(:,1);
    y_tests(:,m) = testdatas{m}(:,1);
end
base_dir = '~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/6cancer_bootstrap_resamples/';
rname = ['SA_nr_1000_nswps_1000_b0_0d01_b1_0d03_nsols20_nts1000_results'];
tname = ['SA_nr_1000_nswps_1000_b0_0d01_b1_0d03_nsols20_nts1000_testperf'];
eval([rname '=r;']);
eval([tname '=t;']);
if saveFlag
    try
        save([base_dir 'results.mat'],rname,tname,'-append')
    catch
        save([base_dir 'results.mat'],rname,tname)
    end
    save([base_dir 'sa_nsols20_ntotsols1000_pred_for_R'],'y_*')
end
