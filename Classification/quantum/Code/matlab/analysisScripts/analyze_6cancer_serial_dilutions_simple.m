function [results,testperf] = analyze_6cancer_serial_dilutions_simple(fracs,d,saveFlag)
% Analysis script for multinomial. Try many different evaluation metrics (auc, 
% balanced accuracy, log loss),ways of combining solutions (iterative, 
% average), number of solutions to combine (20,100). See also
% analyze_lumAB_serial_dilutions_simple.m for a script with one parameter
% selected, based on trying multiple 

if nargin < 3
    saveFlag = false;
end
clearvars results testperf 
cinits = 30;%[1 5 8 10 30];
dir_name = ['~/Dropbox-Work/Wuxi/Results/' d '_splits/'];
for k = 1 : length(cinits)
    cinit = cinits(k);
    for n = 1 : length(fracs)
        load(['~/Dropbox-Work/Wuxi/Data/' d '_splits/frac_' ...
            num2str(fracs(n)) '_data_resamples.mat']) 
        load([dir_name d '_splits_frac_' num2str(fracs(n)) ...
            '_cinit_' num2str(cinit,'%1.1f') '_ri_out_at_5_nr_1000.mat']) 
        for m = 1 : length(out)
            sols = cellfun(@(x) x{1},out{m},'uniformoutput',false);
            [r,t] = analyzeMultinomialResults(sols,traindatas{m},valdatas{m},...
                'metric','acc','nSols',20,'uniqueFlag',true,'lambdas',0,...
                'iterFlag',false,'testFlag',true);
            r.frac = fracs(n);
            t.frac = fracs(n);
            results(n,m) = r;
            testperf(n,m) = t;
        end
    end
    rname = ['cinit' num2str(cinits(k)) '_nsols20_ntotsols1000_results'];
    tname = ['cinit' num2str(cinits(k)) '_nsols20_ntotsols1000_testperf'];
    eval([rname '=results;']);
    eval([tname '=testperf;']);
    if saveFlag
        try
            save([dir_name 'DW_pca_results'],rname,tname,'-append')
        catch
            save([dir_name 'DW_pca_results'],rname,tname)
        end
    end
end
