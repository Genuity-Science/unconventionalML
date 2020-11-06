function [results,testperf] = analyze_serial_dilutions_simple(fracs,d,saveFlag)
% Analysis script for lumAB serial dilutions. Serial dilutions refer to 
% starting with a certain cut of the training data and successively decreasing
% the amount of training data. Selected certain parameters for simplicity (see 
% analyze_lumAB_serial_dilutions.m for more parameters that can try).  

% the fraction of the original training data to use
%fracs = 0.25;
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
    load(['~/Dropbox-Work/Wuxi/Data/' d '_splits/frac_' ...
            num2str(fracs(n)) '_data_resamples.mat']) 

    train_splits = get_split_idxs(size(traindatas{1},1),3);
    % loop over 100 resamples 
    % load result file, output in variable named out
    % out is a cell array of cells. Each cell is 1x3, with the first
    % cell being the solutions, the second the coupling strength, and the 
    % third the fraction of unbroken solutions. See runDW() in DW_utils.py
    
    for m = 1 : length(traindatas)
        trdata = traindatas{m};
        valdata = valdatas{m};
        tic 
        for k = 1 : 3
            if multiFlag
                [h,J] = generateMultinomial(trdata(train_splits{k},:));
            else
                [h,J] = generateLogistic(trdata(train_splits{k},:));
            end
            field_sols{m}{k,1} = -sign(h)'; 
        end
        field_time = toc;
        disp(sprintf('Frac: %1.2f Iteration: %d',fracs(n),m))
        if multiFlag 
            [r,t] = analyzeMultinomialResults(field_sols{m},traindatas{m},...
               valdatas{m}, 'uniqueFlag', false, 'iterFlag', true, ...
               'nSols', 1, 'lambdas',[0],'metric','acc',...
               'biasFlag',false);
        else 
            [r,t] = analyzeLogisticResults(field_sols{m},trdata,valdata,...
                'nSols',1,'uniqueFlag',true,'lambdas',[0],'iterFlag',true,...
                'postProcess','test','biasFlag',false);
        end 
        r.frac = fracs(n);
        r.time = field_time;
        t.frac = fracs(n);
        results(n,m) = r;
        testperf(n,m) = t;
    end
    
    clear *datas*
end

rname = 'field_results';
tname = 'field_testperf';
eval([rname '= results;'])
eval([tname '=testperf;']);

if saveFlag
    save(['~/Dropbox-Work/Wuxi/Results/' d '_splits/field_sols'],'field_sols');
    
    try
        save(['~/Dropbox-Work/Wuxi/Results/' d '_splits/field_results'],rname,tname,'-append')
    catch
        save(['~/Dropbox-Work/Wuxi/Results/' d '_splits/field_results'],rname,tname)
    end
end
