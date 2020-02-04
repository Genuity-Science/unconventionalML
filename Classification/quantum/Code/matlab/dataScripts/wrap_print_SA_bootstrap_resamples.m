% script to print SA files for bootstrap resamples. 

datasets = {'brcaMatchedTN','ERpn','kirckirp','luadlusc','lumAB','lumAB_gene','6cancer'};
%pcs = [25 85 118];


for n = [6] % 5 : length(datasets)
    if n ==7 
        pcs = 13;
        type = 'multi';
    else
        pcs = [44];
        type = 'log';
    end
    d = datasets{n};
    load(['~/Dropbox-Work/Wuxi/Data/' d '_bootstrap_resamples.mat'],'traindatas')
%    d = 'lumAB_diffexp_gene';
    all_traindatas = traindatas;
    for k = 1 : length(pcs)
        traindatas = cellfun(@(x) x(:,1:pcs(k)+1),all_traindatas,'uniformoutput',false);
        for m = 1 : length(traindatas)
            printSACV(traindatas{m},sprintf('~/Dropbox-Work/Wuxi/SA/bootstrap_resamples/%s/Inst%d_pc%d',...
                      d,m,pcs(k)),'lambdas',0,'l_string',{'0'},'inSfx','.txt','type',type)
            disp(['Printing ' d ' Instance ' num2str(m)])
        end
    end
end
