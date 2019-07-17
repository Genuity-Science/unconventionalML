"""
Script to run incremental decrease with dataset given by filename f (see line 22)
"""

# import packages
import gzip, pickle,sys,os
import numpy as np
import scipy.io as sio
file_path = '/home/richard/Dropbox-Work/Wuxi/'
sys.path.insert(1,os.path.join(file_path + 'Code/python'))
import DW_utils

# can try different lambda
#lambdas = [0,0.5,2]
embedding = pickle.load(gzip.open(file_path + '/Embeddings/ISIK44.pkl.gz'))

# define the list of fractions (may differ based on the dataset)
fracs = np.r_[0.06, np.arange(0.1,1.0,0.05)]
frac_str = np.array(str.split(str(fracs).strip('[]')))

cinits = [12.0, 16.0]
f = 'ERpn'

for cinit in cinits :
    for frac in frac_str:

        data_file = file_path + 'Data/' + f + '_splits/frac_' + frac + '_data_resamples.mat'
        print data_file
        savename = file_path + 'Results/' + f + '_splits/' +f + '_splits_frac_' + frac + '_cinit_' + str(cinit) + '_ri_out_at_5_nr_1000.mat'
        if os.path.isfile(savename):
            continue
        traindatas = sio.loadmat(data_file)['traindatas'].flatten()

        # initialize np array to hold outputs
        outlist = np.empty(len(traindatas),dtype=object)
        for n in range(len(traindatas)):
            out = DW_utils.runCV(traindatas[n], [0], embedding, solver_name='ISI', 
                    num_reads=1000, num_gauges=10,coupling_init=cinit,
                    stop_point=-0.01,method='vote',annealing_time=5)
            outlist[n] = out

        # save 
        sio.savemat(savename,{'out':outlist},do_compression=True)
