% analysis script for random solutions for 6 cancer dataset.
base_dir = '~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/6cancer_bootstrap_resamples/';
load('~/Dropbox-Work/Wuxi/Data/6cancer_bootstrap_resamples.mat')
train_splits = get_split_idxs(size(traindatas{1},1),3);
lambdas = 
for n = 1 : 100
    trdata = traindatas{n};
    for m = 1 : length(lambdas)
        for k = 1 : 3
            [h,J] = generateMultinomial(trdata(train_splits{k},:));
            field_sols{n}{k,1} = -sign(h+lambdas(m)/2)'; 
        end
        [results(n,m),testperf(n,m)] = analyzeMultinomialResults(field_sols{n},traindatas{n},testdatas{n},...
                        'metric','acc','nSols',1,'uniqueFlag',true,'lambdas',lambdas(m),...
                                    'iterFlag',false);
        [~,~,y_pred_trains{n}] = getMultinomialAcc(mean(cell2mat(results(n).TestItSols)),traindatas{n});
        [~,~,y_pred_tests{n}] = getMultinomialAcc(mean(cell2mat(results(n).TestItSols)),testdatas{n});
        y_trains(:,n) = traindatas{n}(:,1);
        y_tests(:,n) = testdatas{n}(:,1);
    end
end
rname = 'field_results';
tname = 'field_testperf';
eval([rname '=results;'])
eval([tname '=testperf;'])

save([base_dir 'results.mat'],rname,tname,'-append')
i%save([base_dir 'field_pred_for_R'],'y_*')
