% analysis script for random solutions for 6 cancer dataset.
base_dir = '~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/6cancer_bootstrap_resamples/';
load('~/Dropbox-Work/Wuxi/Data/6cancer_bootstrap_resamples.mat')
for n = 1 : 100
    for m = 1 : 3 
        tmp = rand(1000,65); 
        tmp(tmp<=0.5) = -1; 
        tmp(tmp>0.5) = 1; 
        sols{m,1} = tmp; 
    end
    [results(n),testperf(n)] = analyzeMultinomialResults(sols,traindatas{n},testdatas{n},...
                    'metric','acc','nSols',20,'uniqueFlag',true,'lambdas',0,...
                                'iterFlag',false);
    [~,~,y_pred_trains{n}] = getMultinomialAcc(mean(cell2mat(results(n).TestItSols)),traindatas{n});
    [~,~,y_pred_tests{n}] = getMultinomialAcc(mean(cell2mat(results(n).TestItSols)),testdatas{n});
    y_trains(:,n) = traindatas{n}(:,1);
    y_tests(:,n) = testdatas{n}(:,1);
end
rname = 'rand_nsols20_ntotsols1000_results';
tname = 'rand_nsols20_ntotsols1000_testperf';
eval([rname '=results;'])
eval([tname '=testperf;'])

save([base_dir 'results.mat'],rname,tname,'-append')
save([base_dir 'rand_nsols20_ntotsols1000_pred_for_R'],'y_*')
