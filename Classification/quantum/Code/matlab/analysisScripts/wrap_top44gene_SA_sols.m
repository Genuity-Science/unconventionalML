datasets = {'brcaMatchedTN','ERpn','kirckirp','luadlusc','lumAB'};
Outsfxs = {'','_nSweeps_10000_b1-1','_nSweeps_1000_b1-1'};
for n = 1 : length(datasets)
    d = datasets{n};
    load(['~/Dropbox-Work/Wuxi/Data/' d '_top44genes_data.mat'])
    for m = 1 : length(Outsfxs)
        out = getSACVSols(['~/Dropbox-Work/Wuxi/SA/top44genes/' d '_top44genesFold'],3,'.txt',Outsfxs{m},{'0'});
        sols = {out.sols}';
        base_name = [d '_SA' strrep(Outsfxs{m},'-','_')];
        eval([base_name '_sols' '=sols;'])
        rname = [base_name '_results'];
        tname = [base_name '_testperf'];
        [r,t] = analyzeLogisticResults(sols,traindata, testdata, ...
            'uniqueFlag', true, 'iterFlag', false, 'nSols', 20, 'lambdas',[0],...
            'postProcess','test','biasFlag',false);
        eval([rname '=r;'])
        eval([tname '=t;'])
    end
end
