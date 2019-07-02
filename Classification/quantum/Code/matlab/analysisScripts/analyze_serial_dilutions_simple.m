function [results,testperf] = analyze_serial_dilutions_simple(fracs,d,saveFlag,cinits)

% Analysis script for lumAB serial dilutions. Serial dilutions refer to 
% starting with a certain cut of the training data and successively decreasing
% the amount of training data. Selected certain parameters for simplicity (see 
% analyze_lumAB_serial_dilutions.m for more parameters that can try).  

if nargin < 3
    saveFlag = false;
end
if nargin < 4
    cinits = [1 3 8];
end

% the fraction of the original training data to use
nTotSols = [1000 ];
nSols = [1];% 5 10 20 50 100];
dir_name = ['~/Dropbox-Work/Wuxi/Results/' d '_splits/'];

if strcmp(d,'6cancer')
    multiFlag = true;
else
    multiFlag = false;
end
%nsols = 20;

    % l ad data file
    %   traindata is cell array of different cuts of training data for a given
    %       fraction of the data.
    %   valdata is the original test split and remains unchanged between 
    %       different resampling(confusing, but that's how the collaborators 
    %       originally defined it)
    %   testdata is the portion of the original training cut left over after 
    %       taking a portion of the original training cut for training
    %   exptestdata is the original test data plus whatever fraction of data is
    %       not used for training in the training cut

    % load result file, output in variable named out
    % out is a cell array of cells. Each cell is 1x3, with the first
    % cell being the solutions, the second the coupling strength, and the 
    % third the fraction of unbroken solutions. See runDW() in DW_utils.py
for k = 1 : length(cinits)
    cinit = cinits(k);
    for p = 1 : length(nTotSols);
        for ii = 1 : length(nSols)
            if nSols(ii) > nTotSols(p)
                continue
            end
            % main loop: go through resamples of data
            s=sprintf('Cinit: %1.1f, nTotSols: %d, nSols: %d',cinit,nTotSols(p),nSols(ii));
            disp(s);
%            rname = ['cinit' num2str(cinit) '_ntotsols_'  num2str(nTotSols(p)) '_minsol_results'];
%            tname =  ['cinit' num2str(cinit) '_ntotsols_' num2str(nTotSols(p)) '_minsol_testperf'];
            rname = ['cinit' num2str(cinit) '_at_5_nsols_' num2str(nSols(ii)) ...
                '_ntotsols_'  num2str(nTotSols(p)) '_results'];
            tname =  ['cinit' num2str(cinit) '_at_5_nsols_' num2str(nSols(ii)) ...
                '_ntotsols_' num2str(nTotSols(p)) '_testperf'];
            rname = strrep(rname,'.','d');
            tname = strrep(tname,'.','d');
            for n = 1 : length(fracs)
                load(['~/Dropbox-Work/Wuxi/Data/' d '_splits/frac_' ...
                    num2str(fracs(n)) '_data_resamples.mat']) 
                load([dir_name d '_splits_frac_' num2str(fracs(n)) ...
                    '_cinit_' num2str(cinit,'%1.1f') '_ri_out_at_5_nr_1000.mat']) 

                for m = 1 : length(out)
                    sols = cellfun(@(x) x{1}(1:nTotSols(p),:),out{m},'uniformoutput',false);
                    trdata = traindatas{m};
                    valdata = valdatas{m};
                    disp(sprintf('Frac: %1.2f Iteration: %d',fracs(n),m))
                    if multiFlag 
                        [r,t] = analyzeMultinomialResults(sols,trdata,...
                           valdata, 'uniqueFlag', false, 'iterFlag', true, ...
                           'nSols', nSols(ii), 'lambdas',[0],'metric','acc',...
                           'biasFlag',false);
                    else
                        [r,t] = analyzeLogisticResults(sols,trdata,valdata,...
                            'nSols',nSols(ii),'uniqueFlag',true,'lambdas',[0],'iterFlag',false,...
                            'postProcess','test','biasFlag',false);
                    end
                    r.frac = fracs(n);
                    t.frac = fracs(n);
                    results(n,m) = r;
                    testperf(n,m) = t;
                end
                
            end
            eval([rname '= results;']);
            eval([tname '= testperf;'])
            if saveFlag
                try
                    save([dir_name 'DW_pca_results'],rname,tname,'-append')
                catch
                    save([dir_name 'DW_pca_results'],rname,tname)
                end
            end
        end
    end
    clear *datas*
end

%rname = ['cinit8_at_5_ntotsols_' num2str(ntotsols) '_minsol_results'];
%tname = ['cinit8_at_5_ntotsols_' num2str(ntotsols) '_minsol_testperf'];

