% Analysis script for multinomial. Try many different evaluation metrics (auc, 
% balanced accuracy, log loss),ways of combining solutions (iterative, 
% average), number of solutions to combine (20,100). See also
% analyze_lumAB_serial_dilutions_simple.m for a script with one parameter
% selected, based on trying multiple 

clearvars results testperf 
cinits = [1 5 8 10 30]
metrics = {'acc','log'};
nsols = [20, 100,1000];
cmb_sols = [true,false];
cmb_names = {'iter','avg'};
n_metrics = length(metrics);
n_nsols = length(nsols);


for n = 1 : length(cinits)
    load(['~/Dropbox-Work/Wuxi/Data/6cancer_boostrap_data_resamples.mat'])
    eval(['out = cinit' num2str(cinits(n)) '_ri_out']);
    for m = 1 : length(out)
        sols = cellfun(@(x) x{1},out{m},'uniformoutput',false);
        sprintf('Frac: %1.2f Iteration: %d',fracs(n),m)
        for k = 1 : n_metrics*n_nsols*2
            [ii,ij,ik] = ind2sub([n_metrics n_nsols 2],k);
            [r(k),t(k)] = analyzeMultinomialResults(sols,traindatas{m},...
                exptestdatas{m},'metric',metrics{ii},'nSols',nsols(ij),...
                'uniqueFlag',true,'lambdas',0,'iterFlag',cmb_sols(ik));
        end
        [~,idx] = max([r.trainmeansolacc]);
        best_r = r(idx);
        [ii,ij,ik] = ind2sub([n_metrics,n_nsols,2],idx);
        best_r.hyperparams = sprintf('metric: %s,%s %d',metrics{ii},cmb_names{ik},nsols(ij));
        disp(best_r.hyperparams)
        best_results(m,n) = best_r;
        best_testperf(m,n) = t(idx);
    end
    clear *datas*
end
