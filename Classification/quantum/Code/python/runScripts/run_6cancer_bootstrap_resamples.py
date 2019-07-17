"""
Script to run bootstrap resamples with 6cancer dataset. 
"""

# load necessary packages
import pickle,gzip, numpy as np, scipy.io as sio, os, path

# insert to path to find DW_utils
sys.path.insert(1,os.path.join('/home/richard/Dropbox-Work/Wuxi/Code/python/'))
import DW_utils 

# load embedding
embedding = pickle.load(gzip.open('/home/richard/Dropbox-Work/Wuxi/Embeddings/DW2KQ_EmbedK65.pkl.gz'))
dataset = '6cancer'
    
# array of coupling strengths to try
#cinits = [3.0, 8.0, 30.0]
cinits = [8.0, 12.0, 30.0]

# load the training datas for resamples
base_dir = '/home/richard/Dropbox-Work/Wuxi/Results/bootstrap_resamples/' + dataset + '_bootstrap_resamples/'
traindatas = sio.loadmat('/home/richard/Dropbox-Work/Wuxi/Data/6cancer_bootstrap_resamples.mat')['traindatas'].flatten()
for cinit in cinits :
    outlist_new = np.empty(100,dtype=object)
    for n in range(100):
        tdata = traindatas[n]
        out = DW_utils.runMultiCV(tdata,[0],embedding,solver_name = 'DW',
                num_reads=10000,num_gauges=10,coupling_init=cinit,
                stop_point=0.001,method='vote')
        outlist_new[n] = out
    sio.savemat(base_dir + 'cinit' + str(cinit) + '_ri_sols.mat',{'out':outlist_new},do_compression=True)
