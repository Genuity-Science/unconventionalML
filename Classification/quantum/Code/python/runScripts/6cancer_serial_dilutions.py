"""
Script to run incremental decrease with 6cancer dataset. Must be run on at least a DW2000Q 
"""

# import packages
import gzip, pickle,sys,os
import numpy as np
import scipy.io as sio
file_path = '/home/richard/Dropbox-Work/Wuxi/'
sys.path.insert(1,os.path.join(file_path + 'Code/python'))
import DW_utils

# can try different lambdas
#lambdas = [0,0.5,2]
embedding = pickle.load(gzip.open(file_path + 'Embeddings/DW2KQ_2_1_EmbedK65.pkl.gz'))
f = '6cancer'

# get the list of fractions 
fracs = np.r_[np.arange(0.02,0.21,0.02),np.arange(0.05,0.96,0.05)]
frac_str = np.array(str.split(str(fracs).strip('[]')))

# define list of cinits to try. The higher the value of cinit, the higher the noise 
# (shrinks rest of problem parameters)
cinits = [8.0, 30.0]

for cinit in cinits :
    for frac in frac_str:

        data_file = file_path + 'Data/' + f + '_splits/frac_' + frac + '_data_resamples.mat'
        print data_file
        savename = file_path + 'Results/' + f + '_splits/' +f + '_splits_frac_' + frac + '_cinit_' + str(cinit) + '_ri_out_at_5_nr_1000.mat'
        traindatas = sio.loadmat(data_file)['traindatas'].flatten()
        
        # don't rerun if already saved
        if os.path.isfile(savename):
            continue

        # initialize np array to hold outputs
        outlist = np.empty(len(traindatas),dtype=object)
        for n in range(len(traindatas)):
            out = DW_utils.runMultiCV(traindatas[n], [0], embedding, 
                    solver_name='DW', num_reads=1000, num_gauges=10,
                    coupling_init=cinit,stop_point=-0.01,method='vote',
                    annealing_time=5)

            outlist[n] = out

        # save 
        sio.savemat(savename,{'out':outlist},do_compression=True)
