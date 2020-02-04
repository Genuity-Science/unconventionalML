"""
Script to run bootstrap resamples for binomial datasets.
"""

# import packages, add path for DW_utils
import sys, os
sys.path.insert(1,os.path.join('/home/richard/Dropbox-Work/Wuxi/Code/python'))
import pickle,gzip, numpy as np, scipy.io as sio, DW_utils 

# load embedding
embedding = pickle.load(gzip.open('/home/richard/Dropbox-Work/Wuxi/Embeddings/ISIK44.pkl.gz'))
datasets = ['ERpn','luadlusc','lumAB']#['brcaMatchedTN','ERpn','kirckirp','luadlusc','lumAB']
# try different coupling strengths
#cinits = [0.5, 1.0, 3.0, 8.0, 16.0, 32.0]
cinits = [8.0]
n_pcs = 44
for dataset in datasets:
    print dataset
    base_dir = '/home/richard/Dropbox-Work/Wuxi/Results/bootstrap_resamples/' + dataset + '_bootstrap_resamples/'
    data_file = '/home/richard/Dropbox-Work/Wuxi/Data/' + dataset + '_bootstrap_resamples.mat'
    traindatas = sio.loadmat(data_file)['traindatas'].flatten()

    for cinit in cinits : 
        outlist_new = np.empty(100,dtype=object)
        for n in range(100):
            tdata = traindatas[n][:,:n_pcs+1]
            out = DW_utils.runCV(tdata,[0],embedding,solver_name = 'ISI',num_reads=1000,num_gauges=10,coupling_init=cinit,stop_point=-0.01,method='vote',annealing_time=5)
            outlist_new[n] = out
        sio.savemat(base_dir + 'cinit' +str(cinit) + '_DW_ri_at_5_nr_1000_out.mat',{'out':outlist_new},do_compression=True)
