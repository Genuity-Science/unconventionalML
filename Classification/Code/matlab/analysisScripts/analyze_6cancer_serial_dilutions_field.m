function [field_r,field_t] = analyze_6cancer_serial_dilutions_field(fracs,d,saveFlag)

% Analysis script for lumAB serial dilutions. Serial dilutions refer to 
% starting with a certain cut of the training data and successively decreasing
% the amount of training data. Selected certain parameters for simplicity (see 
% analyze_lumAB_serial_dilutions.m for more parameters that can try).  


% suffix for SA instance names
dir_name = ['~/Dropbox-Work/Wuxi/Results/' d '_splits/'];

for n = 1 : length(fracs)
    frac = num2str(fracs(n));
    load(['~/Dropbox-Work/Wuxi/Data/' d '_splits/frac_' frac '_data_resamples.mat'])
    train_splits = get_split_idxs(size(traindatas{1},1),3);
    lambdas = 0;
    for m = 1 : length(traindatas)
        for k = 1 : 3
            [h,J] = generateMultinomial(traindatas{m}(train_splits{k},:));
            field_sols{m}{k,1} = -sign(h)'; 
        end
        [r,t] = analyzeMultinomialResults(field_sols{m},traindatas{m},valdatas{m},...
                        'metric','acc','nSols',1,'uniqueFlag',true,'lambdas',0,...
                                    'iterFlag',false);
        r.frac = fracs(n);
        t.frac = fracs(n);
        field_r(n,m) = r;
        field_t(n,m) = t;
    end
    rname = 'field_results';
    tname = 'field_testperf';
    eval([rname '=field_r;'])
    eval([tname '=field_t;'])
    
    if saveFlag
        try 
            save([dir_name 'field_results'],'-append',rname,tname)
        catch
            save([dir_name 'field_results'],rname,tname)
        end
    end
end
