% script to print SA files for bootstrap resamples. 

%fracs = [0.02:0.02:0.2 0.25:0.05:0.95];
%fracs = [0.06 0.1:0.05:0.95];
fracs = [0.05:0.05:0.95];
%fracs = [0.10 0.18 0.2:0.05:0.95];
%d = 'ERpn_splits'
%d = 'lumAB_splits'
d = '6cancer_splits';
for n = 1 : length(fracs)
    load(['~/Dropbox-Work/Wuxi/Data/' d '/frac_' num2str(fracs(n)) '_data_resamples.mat'],'traindatas')
    for m = 1 : length(traindatas)
        fbase = ['~/Dropbox-Work/Wuxi/SA/' d '/frac_' num2str(fracs(n)) '_Inst' num2str(m)];
        if contains(d,'6cancer')
            type = 'multi';
        else
            type = 'log';
        end
        printSACV(traindatas{m},fbase,'lambdas',[0],'l_string',{'0'},'inSfx','.txt','type',type)
        disp(['Printing ' d ' Instance ' num2str(m)])
    end
end
