'''
Author: Richard Li
Collection of various utilities to run DW with python package. For unknown
reasons, Matlab would not return all the solutions requested for some choices
of parameters. Python implementation acted as expected.
'''


import numpy as np
from dwave_sapi2.remote import RemoteConnection
from dwave_sapi2.core import async_solve_ising, await_completion #solve_ising
from dwave_sapi2.embedding import embed_problem, unembed_answer
from dwave_sapi2.util import get_hardware_adjacency, qubo_to_ising
from chimera_embedding import processor
from sklearn.model_selection import KFold

import time, pickle, gzip

def connectSolver(solver_name):
    '''
    Connects to solver.
    
    Parameters
    ----------
    
    solver_name: str. Name of solver. Can be 'NASA','ISI',or 'DW'. 
    
    Returns: 
    --------
    
    solver: DW SAPI solver object
    '''
    # connect to solver 
    if solver_name == 'NASA':
        url = 'https://qfe.nas.nasa.gov/sapi';
        token = '' 
        remote_connection = RemoteConnection(url, token)
        solver = remote_connection.get_solver('C16')
    elif solver_name == 'ISI':
        url = 'https://usci.qcc.isi.edu/sapi';
        token = ''
        remote_connection = RemoteConnection(url, token)
        solver = remote_connection.get_solver('DW2X')
    elif solver_name == 'DW':
        url = 'https://cloud.dwavesys.com/sapi'
        token = '' 
        remote_connection = RemoteConnection(url, token)
        solver = remote_connection.get_solver('DW_2000Q_2_1');
    else :
        NameError('Unrecognized solver name')

    return solver

def runDW(h, J, embedding, stop_point=0.25, num_reads=1000, 
        coupling_init=1.0, coupling_increment=0.1, min_solver_calls=1, 
        max_solver_calls=1000, method='vote', last=True, num_gauges=1, 
        solver_name='NASA', annealing_time=20):

    ''' 
    Submits an instance to DW. 
    
    Parameters
    -----
    h : list, a list of fields
    
    J : a dictionary, where keys are a tuple corresponding to the coupling
    
    embedding : a list of lists. Can use DW sapi to generate
    
    stop_point :float, default: 0.25. 
        Stop increasing coupling strength when returns at least this fraction of
        solutions are unbroken.
        
    num_reads: int, default: 1000. 
        The number of reads.
    
    coupling_init: float, default: 1.0. 
        The initial value of coupling, the value of the ferromagnetic coupling 
        between physical qubits. If number of unbroken of solutions is not at 
        least stop_point, then the magnitude of coupling is incremented by 
        coupling_increment. Note however, that the though we specify 
        coupling_init as positive, the coupling is negative. For example, 
        Suppose coupling_init=1.0, coupling_increment (defined below) is 0.1, 
        and stop_point = 0.25. The initial physical ferromagnetic coupling 
        strength will be -1.0. If stop_point isn't reached, coupling is 
        incremented by 0.1, or in other words, the new chain strength is -1.1. 
        coupling is incremented by coupling_increment until stop_point is 
        reached.
        
    coupling_increment: float, default: 0.1. 
        Increment of coupling strength,
    
    min_solver_calls: int, default: 1. 
        The minimum number of solver calls.
    
    max_solver_calls: int, default: 1000. 
        The maximum number of solver calls.
    
    method: str, 'minimize_energy', 'vote', or 'discard', default: 'minimize_energy'
        How to deal with broken chains. 'minimize_energy' uses the energy
        minimization decoding. 'vote' uses majority vote decoding. 'discard' 
        discard broken chains.
    
    last: bool, default: True
        If True, return the last num_reads solutions. If False, return the first
        num_reads solutions.
        
    num_gauges: int, default: 1
        Number of gauge transformations.
        
    solver_name: str, 'NASA', 'ISI', or 'DW', default: 'NASA'
        Which solver to use. 'NASA' uses NASA's DW2000Q. 'ISI' uses ISI's 
        DW2X. 'DW' uses DW's DW2000Q.
        
    Returns
    -------
    A tuple of sols, c, ratio
    
    sols: numpy ndarray, shape = [num_reads, num_spins]
        Solutions where each row is a set of spins (num_spins dependent on 
        solver)
        
    c: float 
        The final value of the coupling strength used

    ratio: float
        The final fraction of unbroken solutions returned at coupling_strength c
    
    '''
    
    meths = ['discard','vote','minimize_energy']
    assert(method in meths)

    solver = connectSolver(solver_name)

    A = get_hardware_adjacency(solver)
    # embed problem
    (h0, j0, jc, new_emb) = embed_problem(h, J, embedding, A)
    
    # scale problem 
    maxjh = max(max(np.abs(h0)),max(np.abs(j0.values())))
    h0 = [el/maxjh for el in h0]
    j0 = {ij:v/maxjh for ij,v in zip(j0.keys(),j0.values())}    

    ratio = 0
    sols = []
    ncalls = 0
    l = coupling_init
    print coupling_init

    jc = dict.fromkeys(jc,-l)
    emb_j = j0.copy()
    emb_j.update(jc)

    kwargs = {'num_reads':num_reads, 'num_spin_reversal_transforms':num_gauges, 'answer_mode':'raw', 
        'annealing_time':annealing_time}
    # iteratively increase ferromagentic strength until returns a certain ratio
    # of solutions where there are no broken chains
    while ratio <= stop_point and ncalls < max_solver_calls :
        jc = dict.fromkeys(jc,-l)
        emb_j = j0.copy()
        emb_j.update(jc)
        if solver_name == 'ISI':
            _check_wait()
        problem = async_solve_ising(solver, h0, emb_j, **kwargs);
        await_completion([problem],1,50000)
        answer = problem.result()
        result = unembed_answer(answer['solutions'],new_emb,broken_chains='discard',h=h,j=J)
        sols += result
        nSols = len(result)
        l = l + coupling_increment
        ratio = nSols/float(len(answer['solutions']))
        ncalls = ncalls+1;
        print ratio

    # Don't remember why do this. Maybe for good measure
    if method == 'discard' :
        nAdd = int(1/ratio)+1
    else:
        nAdd = 1

    l = l - coupling_increment

    for n in range(nAdd):
        if solver_name == 'ISI':
            _check_wait()
        problem = async_solve_ising(solver, h0, emb_j, **kwargs);
        await_completion([problem],1,50000)
        answer = problem.result()
        result = unembed_answer(answer['solutions'],new_emb,broken_chains=method,h=h,j=J)
        sols += result

    if len(sols) < num_reads:
        # try one more time
        problem = async_solve_ising(solver, h0, emb_j, **kwargs);
        await_completion([problem],1,50000)
        answer = problem.result()
        result = unembed_answer(answer['solutions'],new_emb,broken_chains=method,h=h,j=J)
        sols += result
         
    if last:
        if len(sols) >= num_reads :
            sols = sols[-num_reads:]
        else:
            print 'WARNING! DID NOT COLLECT ENOUGH READS...continuing'
    else :
        sols = sols[:num_reads];
    
    return (np.array(sols,dtype=np.int8),l,ratio)

def runDW_batch(h, J, embedding, stop_point=0.25, num_reads=1000, 
        coupling_init=1.0, coupling_increment=0.1, min_solver_calls=1, 
        max_solver_calls=1000, method='vote', last=True, num_gauges=1, 
        solver_name='NASA', returnProblems=True):
    ''' 
    Submits an instance to DW as a batch. Note that when used, sometimes 
    the solutions are markedly different than when use runDW (no batch). 
    Generally using run_DW() seems to be a better idea
    
    Parameters
    -----
    h : list of lists, with each list is a list of fields
    
    J : a list of dictionary, where keys are a tuple corresponding to the
        coupling. Should be the same length as h.
    
    embedding : a list of lists. Can use DW sapi to generate
    
    stop_point :float, default: 0.25.
        Stop increasing coupling strength when returns at least this fraction of
        solutions are unbroken.
        
    num_reads: int, default: 1000. 
        The number of reads.
    
    coupling_init: float, default: 1.0. 
        The initial value of coupling, the value of the ferromagnetic coupling 
        between physical qubits. If number of unbroken of solutions is not at 
        least stop_point, then the magnitude of coupling is incremented by 
        coupling_increment. Note however, that the though we specify 
        coupling_init as positive, the coupling is negative. For example, 
        Suppose coupling_init=1.0, coupling_increment (defined below) is 0.1, 
        and stop_point = 0.25. The initial physical ferromagnetic coupling 
        strength will be -1.0. If stop_point isn't reached, coupling is 
        incremented by 0.1, or in other words, the new chain strength is -1.1. 
        coupling is incremented by coupling_increment until stop_point is 
        reached.
        
    coupling_increment: float, default: 0.1. 
        Increment of coupling strength,
    
    min_solver_calls: int, default: 1. 
        The minimum number of solver calls.
    
    max_solver_calls: int, default: 1000. 
        The maximum number of solver calls.
    
    method: str, 'minimize_energy', 'vote', or 'discard', default: 'minimize_energy'
        How to deal with broken chains. 'minimize_energy' uses the energy
        minimization decoding. 'vote' uses majority vote decoding. 'discard' 
        discard broken chains.
    
    last: bool, default: True
        If True, return the last num_reads solutions. If False, return the first
        num_reads solutions.
        
    num_gauges: int, default: 1
        Number of gauge transformations.
        
    solver_name: str, 'NASA', 'ISI', or 'DW', default: 'NASA'
        Which solver to use. 'NASA' uses NASA's DW2000Q. 'ISI' uses ISI's 
        DW2X. 'DW' uses DW's DW2000Q.

    returnProblems: bool
        Determines what it returns. If True, return problems, new_emb. If False
        return solutions only
        
    Returns
    -------

    if returnProblems is True, returns problems, new_emb (to be used with 
        get_async_sols)
        problems: list 
            list of problems from async_solve_ising
        new_emb: list
            list of embeddings returned from embed_problem

    if returnProblems is False, returns solutions
        sols: np array
            Array of solutions 
    
    '''
    
    meths = ['discard','vote','minimize_energy']
    assert(method in meths)

    if solver_name == 'NASA':
        url = 'https://qfe.nas.nasa.gov/sapi';
        token = 'NASA-870f7ee194d029923ad8f9cd063de357ba53b838';
        remote_connection = RemoteConnection(url, token)
        solver = remote_connection.get_solver('C16')
    elif solver_name == 'ISI':
        url = 'https://usci.qcc.isi.edu/sapi';
        token = 'QUCB-089028555cb44b4f3da34cd4c6dd4a73ec859bc8';
        remote_connection = RemoteConnection(url, token)
        solver = remote_connection.get_solver('DW2X')
    elif solver_name == 'DW':
        url = 'https://cloud.dwavesys.com/sapi'
        token = 'usc-171bafd63a1b07635fd696db283ad4c28b820d14'
        remote_connection = RemoteConnection(url, token)
        solver = remote_connection.get_solver('DW_2000Q_2_1');
    else :
        NameError('Unrecognized solver name')

    A = get_hardware_adjacency(solver)

    h0 = []
    j0 = []
    jc = []
    new_emb = []
    for n in range(len(h)) :
        (h0t,j0t,jct,new_embt) = embed_problem(h[0],J[0],embedding,A)
        maxjh = max(max(np.abs(h0t)),max(np.abs(j0t.values())))
        h0t = [el/maxjh for el in h0t]
        j0t = {ij:v/maxjh for ij,v in zip(j0t.keys(),j0t.values())}
        h0.append(h0t)
        j0.append(j0t)
        jc.append(jct)
        new_emb.append(new_embt)
    ncalls = 0
    sols = np.empty(len(h0),dtype=object)
    if isinstance(coupling_init,list) :
        l = coupling_init
    else :
        l = [coupling_init]*len(h)
    print np.unique(l)


    kwargs = {'num_reads':num_reads, 'num_spin_reversal_transforms':num_gauges, 'answer_mode':'raw'}
    problem = []
    for n in range(len(h0)) :
        jct = dict.fromkeys(jc[n],-l[n]);
        emb_j = j0[n].copy()
        emb_j.update(jct)
        if solver_name == 'ISI':
            _check_wait()
        problem.append(async_solve_ising(solver, h0[n], emb_j, **kwargs));
    await_completion(problem,len(h),50000)
    
    if returnProblems :
        return problem,new_emb
        
    for n in range(len(h0)) :
        answer = problem[n].result()
        sols[n] = np.array(unembed_answer(answer['solutions'],new_emb[n],broken_chains=method,h=h[n],j=J[n]),dtype=np.int8)
#    return problem,new_emb
    return np.array(sols)
    
def get_async_sols(problems,hs,Js,new_embs,method='vote') :
    """
    Returns solution for async problem.
    
    Parameters
    ----------
    problems: list
        list of problem handles returned by async_solve_ising

    hs: list of lists
        local fields

    Js: list of dicts
        problem couplings

    new_embs: list
        embeddings returned from embed_problem

    method: str
        decoding method. Can be 'minimize_energy','vote', or 'discard'

    Returns
    -------

    """
    sols = np.empty(len(hs),dtype=object)
    for n in range(len(hs)): 
        answer = problems[n].result()
        sols[n] = np.array(unembed_answer(answer['solutions'],new_embs[n],broken_chains=method,h=hs[n],j=Js[n]),dtype=np.int8)

    return sols

def _check_wait() :

    '''
    Internal function for use with ISI's DW. Check whether should wait to 
    submit problem, depending on ISI's schedule.
    '''

    times = np.array([(10800, 17970), (30600, 37770), (47700, 54870), (64800, 68370),
             (80100, 87270)])

    t = list(time.localtime())
    if t[3] == 0 and t[4] < 14:
        t[3] = 24
    current_time = np.dot(t[3:6],[3600,60,1]) 
    ar = np.sum(times < current_time,axis=1)
    idx = np.argmax(ar==0)

    if ar[idx-1] != 1 :
        wait_time = times[idx][0] - current_time
        print ("Pausing for %d seconds..." % wait_time)
        time.sleep(wait_time)
        
def generate_qubo(data):

    '''
    Generates a QUBO as analog to a linear regression.

    Parameters
    ----------
    data: np array
        Training matrix. First column should be response variable

    Returns
    -------
    A tuple of Q,h,J
    Q: dict
        A dictionary corresponding to the QUBO matrix
    h: list
        local fields
    J: dict
        the non-zero couplings

    '''
    D = data[:,1:]
    y = data[:,0]
    X = np.dot(D.T,D)
    l = 2*np.dot(D.T,y)
    Q = X - np.diag(l)
    h = list(np.diag(Q))
    J = np.triu(Q,1)
    J_dict = { index:2*v for index, v in np.ndenumerate(J) if not np.isclose(v,0)}
    return (dict(np.ndenumerate(Q)),h,J_dict)
    
def q_to_i(Q) :
    # wrapper function for DW's qubo_to_ising
    # h,J,offset = q_to_i(Q)
    return qubo_to_ising(Q)
    
def generate_logistic_params(data,isingFlag=True):

    '''
    Generates parameters corresponding if treat problem as logistic regression

    Parameters
    ----------
    data: np array
        Training matrix. First column should be response variable
    
    isingFlag: bool
        True if return parameters as Ising, False if return as QUBO

    Returns
    ------

    tuple of h, J_dict,J

    h: np array
        local fields
    J_dict: dict
        couplings (DW Sapi needs iterable)
    J: np array
        matrix of couplings

    '''
    D = data[:,1:]
    y = data[:,0]
    lh = np.dot(0.5-y,D)
    Q = 1./8*np.dot(D.T,D)
    if isingFlag:
        h = lh
        J = np.triu(Q,1)
        J_dict = { index:2*v for index, v in np.ndenumerate(J) if not np.isclose(v,0)}
        return h,J_dict,J*2
    else :
        Q = Q + np.diag(lh)
        h,J_dict,offset = q_to_i(dict(np.ndenumerate(Q)))
        J = np.zeros((len(h),len(h)))
        for key, value in J_dict.iteritems():
            J[key] = value
        return np.array(h),J_dict,J

def generate_multinomial(data,isingFlag=True) : 
    
    '''
    Generates the h's and J's for a multinomial regression problem. Given K
    classes and P features, will return h's and J's for P*(K-1) weights.
   
    Parameters
    ---------- 
    data: a numpy matrix. 
        Training data. Assumes the first column is the classes

    isingFlag: bool
        True if return parameters as Ising, False if return as QUBO

    Returns
    -------
    h, J_dict, J

    h: np array
        local fields
    J_dict: dict
        couplings (DW Sapi needs iterable)
    J: np array
        matrix of couplings

    '''
    
    class_labels = data[:,0]
    D = data[:,1:]
    P = D.shape[1]
    labels = np.unique(class_labels);
    K = len(labels)

    hprime = np.sum(D,axis=0)/float(K)
    X = np.dot(D.T,D)
    Jp = (K-1.)/(2*K*K)*X
    Jpp = -Jp/(K-1.) 

    tmp = np.triu(np.ones((K-1,K-1)),1)
    Q = np.kron(tmp,Jpp) + np.kron(np.eye(K-1),Jp) 
    
    bs = np.zeros((K-1)*P)
    for n in range(K-1) : 
        label = labels[n]
        bs[n*P:(n+1)*P] = -np.sum(D[np.where(class_labels==label)],axis=0)
    
    lh = np.tile(hprime,K-1) + bs
    if isingFlag:
        h = lh
        J = np.triu(Q,1)
        J_dict = { index:2*v for index, v in np.ndenumerate(J) if not np.isclose(v,0)}
        return h,J_dict,J*2
    else :
        Q = Q + np.diag(lh)
        h,J_dict,offset = q_to_i(dict(np.ndenumerate(Q)))
        J = np.zeros((len(h),len(h)))
        for key, value in J_dict.iteritems():
            J[key] = value
        return np.array(h),J_dict,J

def runCV(data,lambdas,embedding,k=3,logFlag=True,kSplit=False,**kwargs):
    '''
    Runs a k-fold CV on dataset. 

    Parameters
    ----------
    data: np matrix
        Training data. First column is labels
    lambdas: list (or array)
        Regularization strengths to use 
    embedding: list of lists
        Embedding for DW
    k: int
        number of folds. By default k=3
    logFlag: bool
        if True, run as logistic regression (i.e., use generate_logistic_params)
        if False, run as regression.
        By default logFlag=True
    kSplit: bool
        if True, use 1/4 of data for testing, 1/k for training
        if False, use 1/k for training, remainder for testing
    **kwargs: dict
        Parameters for DW (num_reads,annealing_time,etc.; see DW documentation
        for list of parameters)

    Returns
    -------
    out: np array, dtype=object
        Of size k, len(lambdas). Array of outputs from runDW (tuple of 
        solutions, ferromagnetic coupling strength, and fraction of unbroken 
        solutions)
    
    '''
    
    if kSplit:
        N = data.shape[0]
        order = np.arange(N)
        kflist = []
        len_test = N/4
        for n in range(k) :
            start_idx = n*N/k
            end_idx = (n+1)*N/k
            train_idx = order[start_idx:(n+1)*N/k]
            new_order = np.r_[order[end_idx:],order[:end_idx]]
            test_idx = new_order[0:len_test]
            kflist.append((train_idx,test_idx))
    else:
        kf = KFold(n_splits=k)
        kflist = list(kf.split(data));
    
    out = np.empty((k,len(lambdas)),dtype=object)   
    count = 0
   
    for train_idx, test_idx in kflist:
    #for train_idx, test_idx in fivesplitlist[3:5] :
        train_data = data[train_idx]
        if logFlag:
            h,J,_ = generate_logistic_params(train_data)  
        else:
            Q,hi,Ji = generate_qubo(train_data)
            h,J,_ = q_to_i(Q)
        for m in range(len(lambdas)):
            # basically run with the input coupling strength, don't get to stopping ratio
           out[count,m] = runDW(h+lambdas[m]/2,J,embedding,**kwargs)
        count+=1

    return out

def getLogAsyncInstances(datas,k=3,lambdas=[0],isingFlag=True) :
    '''
    Generate list of instances to be used with runDW_batch. Kind of deprecated.

    Parameters
    ----------
    datas: list
        A list of training datas (e.g., multiple binomial comparisons for
        for different cancer types)
    k: int
        The number of splits to use
    lambdas: list
        Strength of regularization
    isingFlag: bool 
        If True, treat problem as Ising
        If False, treat problem as QUBO

    Returns
    -------
    hs: list
        list of local fields. Length is k x len(datas)
    Js: list
        list of dicts. 
    '''

    if k > 1:
        kf = KFold(n_splits=k)
    hs = []
    Js = []
    for n in range(len(datas)):
        data = datas[n]
        if k > 1:
            kflist = list(kf.split(data))
        else :
            kflist = [[range(data.shape[0]),[]]]
        for train_idx,test_idx in kflist:
            train_data = data[train_idx]
            h,J,_ = generate_logistic_params(train_data,isingFlag)
            for l in lambdas :
                hs.append(h+l/2)
                Js.append(J)
    return (hs,Js)

def runLogisticRepeatedKFold(data,lambdas,embedding,frac,Nrep=50,**kwargs) :
    """
    Does a repeated k-fold CV on dataset.

    Parameters
    ----------
    data: np array
        Training data. Assumes first column is the labels

    lambdas: list
        Strength of regularization.

    embedding: list of lists
        For use with DW SAPI

    frac: float
        fraction of data to use for training split

    Nrep: int
        Number of repetitions to do. Default Nrep=50

    **kwargs: dict
        Parameters for DW (num_reads,annealing_time,etc.; see DW documentation
        for list of parameters)
        
    Returns
    -------
    out: np array, dtype= object
        Of size Nrep, len(lambdas). Array of output from runDW
    """
    
    N = data.shape[0]
    L = int(N*frac)
    out = np.empty((Nrep,len(lambdas)),dtype=object)
    for n in range(Nrep):
        rp = np.random.permutation(N)
        train_idx = rp[0:L]
        h,J,_ = generate_logistic_params(data[train_idx])
        for m in range(len(lambdas)):
            tmp = runDW(h+lambdas[m]/2,J,embedding,**kwargs)
            out[n,m] = tmp + (train_idx,)
        print ("Finished rep %d of %d" % (n,Nrep))
        
    return out 

def runMultiCV(data,lambdas,embedding,k=3,kSplit=False,isingFlag=True,cinit=None,**kwargs):
    '''
    Runs a k-fold CV on dataset for multinomial regression.

    Parameters
    ----------
    data: np array
        Training data. Assumes first column is the labels

    lambdas: list
        Strength of regularization.

    embedding: list of lists
        For use with DW SAPI
    
    k: int
        number of folds. By default k=3
    
    kSplit: bool  
        if True, use 1/4 of data for testing, 1/k for training
        if False, use 1/k for training, remainder for testing

    isingFlag: bool
        True if return parameters as Ising, False if return as QUBO

    cinit: list or None
        list of initial coupling strengths to use. If None, use 1.0. This is
        a little convoluted and could be improved, because cinit may also
        appear in **kwargs.

    **kwargs: dict
        Parameters for DW (num_reads,annealing_time,etc.; see DW documentation
        for list of parameters)

    Returns
    -------
    out: np array, dtype= object
        Of size k, len(lambdas). Array of output from runDW
    '''
    
    if kSplit:
        N = data.shape[0]
        order = np.arange(N)
        kflist = []
        len_test = N/4
        for n in range(k) :
            start_idx = n*N/k
            end_idx = (n+1)*N/k
            train_idx = order[start_idx:(n+1)*N/k]
            new_order = np.r_[order[end_idx:],order[:end_idx]]
            test_idx = new_order[0:len_test]
            kflist.append((train_idx,test_idx))
    else:
        kf = KFold(n_splits=k)
        kflist = list(kf.split(data));
    
    out = np.empty((k,len(lambdas)),dtype=object)   
    count = 0
    if isinstance(cinit,(int,float)) :
        cinits = [cinit]*k
    elif cinit is None:
        if 'coupling_init' in kwargs:
            cinits = [kwargs['coupling_init']]*k
        else :
            cinits = [1.0]*k
    else :
        cinits = cinit  
    for train_idx, test_idx in kflist:
    #for train_idx, test_idx in fivesplitlist[3:5] :
        train_data = data[train_idx]
        h,J,_ = generate_multinomial(train_data,isingFlag)
        kwargs['coupling_init'] = cinits[count]
        for m in range(len(lambdas)):
            # basically run with the input coupling strength, don't get to stopping ratio
            try :
               out[count,m] = runDW(h+lambdas[m]/2,J,embedding,**kwargs)
            except :
                pickle.dump(out,open('tmpout.pkl','w'))
                raise 
        print 'Completed run {}'.format(count)
        count+=1

    return out    

def heuristicEmbed(solver_name,K) :
    '''
    Finds heuristic embedding for solver_name. 
    
    Parameters:
    ----------
    solver_name: str. Name of solver to connect to. Can be 'NASA','DW', or 'ISI'.
    
    K: the length of embedding. If K is an int, generates a complete graph 
        embedding. Otherwise, assumes K is an iterable with couplings specified.
    '''

    solver = connectSolver(solver_name) 
    hardware_adj = get_hardware_adjacency(solver)

    # generate adjacency for complete graph
    if isinstance(K,int) :
        J = np.triu(np.ones((K,K)))
        J_it = [index for index, v in np.ndenumerate(J) if not np.isclose(v,0)]
    else:
        J_it = K
    embedding = find_embedding(J_int,hardware_adj)
    
    return embedding
    
def balancedEmbed(solver_name,K):
    '''
    Finds balanced embedding for solver_name. 
    
    Parameters:
    ----------
    solver_name: str. Name of solver to connect to. Can be 'NASA','DW', or 'ISI'.
    
    K: the length of embedding. Assumes a complete graph embedding
    '''
    solver = connectSolver(solver_name) 
    hardware_adj = get_hardware_adjacency(solver)
    if solver_name == 'ISI':
        M = 12
        N = 12
    else:
        M = 16
        N = 16
        
    embedder = processor(hardware_adj, M=M, N=N, L=4)
    embedding = embedder.tightestNativeClique(K)
    
    return embedding
