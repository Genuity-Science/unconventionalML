% script to print SA files for bootstrap resamples. 

datasets = {'brcaMatchedTN','ERpn','kirckirp','luadlusc','lumAB'};
%pcs = [25 85 118];
pcs = [44];
for n = 2 : length(datasets)
    d = datasets{n};
    load(['~/Dropbox-Work/Wuxi/Data/' d '_pc118_bootstrap_resamples.mat'],'traindatas')
    all_traindatas = traindatas;
    for k = 1 : length(pcs)
        traindatas = cellfun(@(x) x(:,1:pcs(k)+1),all_traindatas,'uniformoutput',false);
        for m = 1 : length(traindatas)
            printSACV(traindatas{m},sprintf('~/Dropbox-Work/Wuxi/SA/bootstrap_resamples/%s/Inst%d_pc%d',d,m,pcs(k)),'lambdas',0,'l_string',{'0'},'inSfx','.txt')
            disp(['Printing ' d ' Instance ' num2str(m)])
        end
    end
end
