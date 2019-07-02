% Analysis script for multinomial. Try many different evaluation metrics (auc, 
% balanced accuracy, log loss),ways of combining solutions (iterative, 
% average), number of solutions to combine (20,100). See also
% analyze_lumAB_serial_dilutions_simple.m for a script with one parameter
% selected, based on trying multiple 

clearvars results testperf 
%cinits = [1 5 8 10 30];
cinits = [3 8 12 30];
load('~/Dropbox-Work/Wuxi/Data/6cancer_bootstrap_resamples.mat')
base_dir = '~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/6cancer_bootstrap_resamples/';
for n = 1 : length(cinits)
    load([base_dir sprintf('cinit%1.1f_ri_sols.mat',cinits(n))]);
    for m = 1 : length(out)
        sols = cellfun(@(x) x{1},out{m},'uniformoutput',false);
%        sols = cellfun(@(x) x(1:5,:),sols,'uniformoutput',false);
        [r(m),t(m)] = analyzeMultinomialResults(sols,traindatas{m},testdatas{m},...
            'metric','acc','nSols',20,'uniqueFlag',true,'lambdas',0,...
            'iterFlag',false,'testFlag',true);
        [~,~,y_pred_trains{m}] = getMultinomialAcc(mean(cell2mat(r(m).TestItSols)),traindatas{m});
        [~,~,y_pred_tests{m}] = getMultinomialAcc(mean(cell2mat(r(m).TestItSols)),testdatas{m});
        y_trains(:,m) = traindatas{m}(:,1);
        y_tests(:,m) = testdatas{m}(:,1);
    end
    rname = ['cinit' num2str(cinits(n)) '_nsols20_ntotsols1000_results'];
    tname = ['cinit' num2str(cinits(n)) '_nsols20_ntotsols1000_testperf'];
    eval([rname '=r;']);
    eval([tname '=t;']);
    save(['~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/6cancer_bootstrap_resamples/results.mat'],rname,tname,'-append')
end
