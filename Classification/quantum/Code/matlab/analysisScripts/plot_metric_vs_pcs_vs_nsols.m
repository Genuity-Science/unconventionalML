pcs = [24 44 65 84 117];
nsols = [20 50 100 1000];

datasets = {'brcaMatchedTN','ERpn','kirckirp','luadlusc','lumAB'};
for d = 1 : length(datasets)
    load(['~/Dropbox-Work/Wuxi/Results/bootstrap_resamples/' datasets{d} ...
            '_bootstrap_resamples/results'],'SA*testperf*','rand*testperf')
    for n = 1 : length(pcs)
        for m = 1 : length(nsols)
            try
                eval(['rt = rand_pc' num2str(pcs(n)) '_nsols_' ...
                        num2str(nsols(m)) '_testperf;'])
            catch
                eval(['rt = rand_pc' num2str(pcs(n)+1) '_nsols_' ...
                        num2str(nsols(m)) '_testperf;'])
            end
            eval(['st = SA_nr_1000_nswps_1000_pc' num2str(pcs(n)) '_nsols_' ...
                    num2str(nsols(m)) '_testperf_simple;'])
            SA_metric{d}(n,m) = mean([st.meansolsBacc]);
            rand_metric{d}(n,m) = mean([rt.meansolsBacc]);
            std_diff_metric{d}(n,m) = std([st.meansolsBacc]-[rt.meansolsBacc]);
        end
    end
    clear *testperf*
end

figure
for n = 1 : 5 
    subplot(2,3,n); 
    surf(nsols,pcs,SA_metric{n}-rand_metric{n}); 
    set(gca,'xscale','log'); 
    title(datasets{n});
    xlabel('Number of solutions'); 
    ylabel('Number of pcs')
    zlabel('SA bal acc - Rand bal acc')
end

figure; 
for m = 1 : 5 
    subplot(2,3,m); 
    hold on; 
    for n = 1 : 5 
        errorbar(nsols,SA_metric{m}(n,:)-rand_metric{m}(n,:),std_diff_metric{m}(n,:)/10); 
    end
    set(gca,'xscale','log');
    xlim([10 1300])
    xlabel('Number of solutions')
    ylabel('SA bal acc - Rand bal acc')
    title(datasets{m})
end
legend(strsplit(num2str(pcs)))
%surf(pcs,nsols,SA_metric)
%surf(pcs,nsols,rand_metric)
