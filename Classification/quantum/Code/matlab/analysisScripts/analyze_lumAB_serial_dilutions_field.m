% Analysis script for lumAB serial dilutions. Serial dilutions refer to 
% starting with a certain cut of the training data and successively decreasing
% the amount of training data. Selected certain parameters for simplicity (see 
% analyze_lumAB_serial_dilutions.m for more parameters that can try).  

clearvars *results* *testperf* 

% the fraction of the original training data to use
fracs = [0.18 0.2:0.05:0.95];
%fracs = 0.25;
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

    train_splits = get_split_idxs(size(traindatas2{1},1),3);
    % loop over 100 resamples 
    % load result file, output in variable named out
    % out is a cell array of cells. Each cell is 1x3, with the first
    % cell being the solutions, the second the coupling strength, and the 
    % third the fraction of unbroken solutions. See runDW() in DW_utils.py
    
    for m = 1 : length(traindatas2)
        trdata = traindatas2{m};
        valdata = valdatas2{m};

        for k = 1 : 3
            [h,J] = generateLogistic(trdata(train_splits{k},:));
            field_sols{m}{k,1} = -sign(h)'; 
        end
        disp(sprintf('Frac: %1.2f Iteration: %d',fracs(n),m))
        [r,t] = analyzeLogisticResults(field_sols{m},trdata,valdata,...
                'nSols',1,'uniqueFlag',true,'lambdas',[0],'iterFlag',true,...
                'postProcess','test','biasFlag',false);
        r.frac = fracs(n);
        t.frac = fracs(n);
        results(n,m) = r;
        testperf(n,m) = t;
    end
    
    clear *datas*
end

%sname = 'rand_pc44_nsols5_ntotsols_20_sols';
rname = 'field_pc44_results';
tname = 'field_pc44_testperf';
eval([rname '= results;'])
eval([tname '=testperf;']);
%eval([sname '=rand_sols;']);
save('~/Dropbox-Work/Wuxi/Results/lumAB_splits/field_pca_results',rname,tname,'-append')
