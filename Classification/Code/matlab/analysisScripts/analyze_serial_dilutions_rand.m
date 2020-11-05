function [results,testperf] = analyze_serial_dilutions_rand(fracs,d,saveFlag)

% Analysis script for lumAB serial dilutions. Serial dilutions refer to 
% starting with a certain cut of the training data and successively decreasing
% the amount of training data. Selected certain parameters for simplicity (see 
% analyze_lumAB_serial_dilutions.m for more parameters that can try).  

% the fraction of the original training data to use
nTotSols =1000; %[5 50]; % [20 100 1000];
nMaxSols = 1000;
nSols = 20;%[1 5 10 20 50];% 5 10 20 100 1000];
if strcmp(d,'6cancer')
    multiFlag = true;
    Nfeats = 65;
else
    multiFlag = false;
    Nfeats = 44;
end
for p = 1 : length(nTotSols);
    for ii = 1 : length(nSols)
        if nSols(ii) > nTotSols(p)
            continue
        end
        % main loop: go through resamples of data
        s=sprintf('nTotSols: %d, nSols: %d',nTotSols(p),nSols(ii));
        disp(s);
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
        
%            load('~/Dropbox-Work/Wuxi/Results/lumAB_splits/rand_sols','rand_pc44_nsols5_sols');
%            rand_sols = rand_pc44_nsols5_sols;
            % load result file, output in variable named out
            % out is a cell array of cells. Each cell is 1x3, with the first
            % cell being the solutions, the second the coupling strength, and the 
            % third the fraction of unbroken solutions. See runDW() in DW_utils.py
            
            for m = 1 : length(traindatas)
                tic
                for k = 1 : 3
                    tmp = rand(nMaxSols,Nfeats);;
                    tmp(tmp<=0.5) = -1;
                    tmp(tmp>0.5) = 1;
                    rand_sols{m,n}{k,1} = tmp;
                end
                ti = toc;
                trdata = traindatas{m};
                valdata = valdatas{m};
                disp(sprintf('Frac: %1.2f Iteration: %d',fracs(n),m))

                tmp_sols = cellfun(@(x) x(1:nTotSols(p),:), ...
                    rand_sols{m,n},'uniformoutput',false);
                    if multiFlag 
                        [r,t] = analyzeMultinomialResults(tmp_sols,trdata,...
                           valdata, 'uniqueFlag', false, 'iterFlag', true, ...
                           'nSols',nSols(ii), 'lambdas',[0],'metric','acc',...
                           'biasFlag',false);
                    else
                trdata = traindatas{m};
                trdata(:,2:end) = trdata(:,2:end)*1;
                [r,t] = analyzeLogisticResults(tmp_sols,trdata,...
                        valdatas{m}, 'uniqueFlag', false, 'iterFlag', true, ...
                        'nSols', nSols(ii), 'lambdas',[0],'postProcess','test',...
                        'biasFlag',false,'metric','bacc');
                    end
                r.frac = fracs(n);
                t.frac = fracs(n);
                r.time = ti;
                results(n,m) = r;
                testperf(n,m) = t;
            end
            clear *datas*
        end

        if saveFlag
            sname =  ['rand_ntotsols_' num2str(nTotSols(p)) '_sols'];
            eval([sname '=rand_sols;'])
            save(['~/Dropbox-Work/Wuxi/Results/' d '_splits/rand_sols'],sname)

            rname = ['rand_nsols' num2str(nSols(ii)) '_ntotsols' ...
                num2str(nTotSols(p)) '_results'];
            tname = ['rand_nsols' num2str(nSols(ii)) '_ntotsols' ...
                num2str(nTotSols(p)) '_testperf'];
            eval([rname '= results;'])
            eval([tname '=testperf;']);
            try
                save(['~/Dropbox-Work/Wuxi/Results/' d '_splits/rand_results'],...
                    rname,tname,'-append')
            catch
                save(['~/Dropbox-Work/Wuxi/Results/' d '_splits/rand_results'],...
                    rname,tname)
            end
        end
    end
end

