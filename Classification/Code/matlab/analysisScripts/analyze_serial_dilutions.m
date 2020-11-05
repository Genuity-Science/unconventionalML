% Analysis script for lumAB serial dilutions. Serial dilutions refer to 
% starting with a certain cut of the training data and successively decreasing
% the amount of training data. Try many different evaluation metrics (auc, 
% balanced accuracy, log loss),ways of combining solutions (iterative, 
% average), number of solutions to combine (20,100). See also
% analyze_lumAB_serial_dilutions_simple.m for a script with one parameter
% selected, based on trying multiple 

clearvars results testperf 
% the fraction of the original training data to use
fracs = 0.1:0.05:0.95;
% define search space over evaluation metrics
metrics = {'auc','bacc','log'}; %auc, balanced accuracy, log loss
nsols = [20, 100]; % number of solutions to try to combine
cmb_sols = [true,false]; % true means iterative, false means average
cmb_names = {'iter','avg'};

n_metrics = length(metrics);
n_nsols = length(nsols);
dir_name = '~/Dropbox-Work/Wuxi/Results/lumAB_splits/'
% loop through fractions
if strcmp(d,'6cancer')
    multiFlag = true;
else
    multiFlag = false;
end
for n = 1 : length(fracs)
    % load data file
    %   traindata is cell array of different cuts of training data for a given
    %       fraction of the data.
    %   valdata is the original test split and remains unchanged between 
    %       different resampling(confusing, but that's how the collaborators 
    %       originally defined it)
    %   testdata is the portion of the original training cut left over after 
    %       taking a portion of the original training cut for training
    %   exptestdata is the original test data plus whatever fraction of data is
    %       not used for training in the training cut
    load(['~/Dropbox-Work/Wuxi/Data/lumAB_splits/frac_' ...
        num2str(fracs(n)) '_data_resamples.mat'])

    % load result file. Output in variable called out
    % out is a cell array of cells. Each cell is 1x3, with the first
    % cell being the solutions, the second the coupling strength, and the 
    % third the fraction of unbroken solutions. See runDW() in DW_utils.py
    load([dir_name 'lumAB_splits_frac_' num2str(fracs(n)) ...
        '_cinit_8.0_ri_out_at_5_nr_1000.mat'])

    for m = 1 : length(out)
        sols = cellfun(@(x) x{1}(:,1:44),out{m},'uniformoutput',false);
        sprintf('Frac: %1.2f Iteration: %d',fracs(n),m)
        trdata = traindatas2{m};
        valdata = valdatas2{m};
        % flatten for loops so code looks nicer (maybe faster too?)
        for k = 1 : n_metrics*n_nsols*2
            [ii,ij,ik] = ind2sub([n_metrics n_nsols 2],k);
            % *datas2 is a scaling of the data that gave best results. 
            if multiFlag
                [r,t] = analyzeMultinomialResults(sols,trdata,...
                    valdata, 'uniqueFlag', false, 'iterFlag', true, ...
                    'nSols', 1, 'lambdas',[0],'metric','acc',...
                    'biasFlag',false);
            else
                [r(k),t(k)] = analyzeLogisticResults(sols,trdata,valdata,...
                'metric',metrics{ii},'nSols',nsols(ij),'uniqueFlag',true,...
                'lambdas',0,'iterFlag',cmb_sols(ik),'postProcess','test');
            end
        end

        % find set of parameters that perform best on balanced accuracy (could
        % select other metric)
        [~,idx] = max(mean([r.TrainItSolOnTestBacc]));
        best_r = r(idx);
        [ii,ij,ik] = ind2sub([n_metrics,n_nsols,2],idx);
        best_r.hyperparams = sprintf('metric: %s,%s %d',metrics{ii},cmb_names{ik},nsols(ij));
        disp(['bacc: ' num2str(t(idx).best_meansolsBacc) ' ' best_r.hyperparams])
        best_results(m,n) = best_r;
        best_testperf(m,n) = t(idx);
    end
    clear *datas*
    pca_best_results = best_results;
    pca_best_testperf = best_testperf;
    % save
    save([dir_name 'pca_results'],'pca_best_results', 'pca_best_testperf', '-append')
end
