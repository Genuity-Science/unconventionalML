'''
Author: Richard Li

Collection of functions to analyze DW results in python. I originally thought 
this would run faster in Python than in Matlab, but it turns out that my
Matlab implementation was faster. Keeping just in case. For documentation look
in corresponding Matlab files.

'''

import numpy as np
#from import sklearn.metrics average_precision_score, matthews_corrcoef,roc_auc_score
from DW_utils import generate_multinomial
def getIterSols(sols, data, metric, sign, pred=lambda x:x, N=20, methodFlag=0):
    '''
    Given an np array of sols, finds the iterative solutions based on the metric.

    Parameters
    ----------

    sols: array 
        solutions returned from DW.
    data: array 
        the training data matrix from which sols was found
    metric: function
        the metric to be evaluated. of the form metric(y_act,y_pred)
    sign : int
        whether positive or negative. Depends on the metric. By default this 
        function minimizes the objective function, so if goal is to maximize 
        the auc, for example, should be -1. 
    methodFlag: int
        how to treat the inputs.
        methodFlag = 0 means treat everything in {-1,1}
        methodFlag = 1 means treat all solutions as {0,1}
    '''
    assert(sign==-1 or sign==1)
    sols[sols==3] = 0
    sol = sols[0]
    
    if methodFlag == 1 :
        sol[sol<0] = 0
        
    count = 2.0
    n_sols = sols.shape[0]
    L = min(N,n_sols)
    objf = np.empty(L)
    y_act = data[:,0]
    objf[0] = metric(y_act,pred(np.dot(data[:,1:],sol)))*sign
    
    for n in range(1,L):
        if methodFlag == 1 :
            sols[n][sols[n]<0] = 0
        tmp_sol = sol + (sols[n] - sol)/count
        tmp = metric(y_act,pred(np.dot(data[:,1:],tmp_sol)))*sign
        if tmp < objf[n-1]:
            objf[n] = tmp
            sol = tmp_sol
            count += 1
        else:
            objf[n] = objf[n-1]
    return sol,objf

def getCVIterSols(sol_list,data,folds,metric,sign,pred=lambda x:x, N=20, methodFlag=0) :
    '''
    Assumes sol_list is same length as folds. Return a list of the same thing,
    where each item in list has 4 items: the training itersol, training objf,
    testing itersol, test obj
    '''

    results = []
    for n,sols in enumerate(sol_list):
        trainout = getIterSols(sols,data[folds[n][0]],metric,sign,pred=pred,
            N=N,methodFlag=methodFlag)
        testout = getIterSols(sols,data[folds[n][1]],metric,sign,pred=pred,
            N=N,methodFlag=methodFlag)
            
        results.append(trainout+testout)

    return results
    
def calcEns(sols,h,J) :
    h = np.array(h)
    J_mat = np.zeros((len(h),len(h)))
    for key in J.keys():
        J_mat[key] = J[key]
    
    ens = [np.dot(sol,h)+np.dot(sol,np.dot(sol,J_mat)) for sol in sols]
    return np.array(ens)
    
    
def convertToClassifer(train_data,skipFirstColumn=True,asIsing=False) :
    # for now converts to classification problem based on means. Returns an
    # accuracy rate which shows the performance of the classifiers
    
    y = train_data[:,0]
    feat_mat = train_data[:,1:]

    if skipFirstColumn:
        feat_mat = feat_mat[:,1:]
        
    new_feat_mat = np.copy(feat_mat)
    
    t_idx = np.where(y==1)[0]
    f_idx = np.where(y==0)[0]
    feat_means = [(np.mean(feat_mat[t_idx,n]),np.mean(feat_mat[f_idx,n])) 
        for n in range(feat_mat.shape[1])]

    acc_rate = np.empty(feat_mat.shape[1])
    for i,val in enumerate(feat_means):
        m = np.mean(val)
        tmp = _sigmoid((feat_mat[:,i]-m)*np.sign(val[0]-val[1]))
        if asIsing:
            new_feat_mat[:,i] = tmp*2 - 1
        else:
            new_feat_mat[:,i] = tmp
            
        tmp[tmp<0.5] = 0
        tmp[tmp>=0.5]= 1
        acc_rate[i] = np.equal(tmp,train_data[:,0]).sum()/float(train_data.shape[0])


    if skipFirstColumn:
        new_feat_mat = np.hstack((np.reshape(train_data[:,1],(-1,1)),new_feat_mat))
        
    new_train_data = np.hstack((np.reshape(y,(-1,1)),new_feat_mat))
    
    return new_train_data,acc_rate
    
def _sigmoid(x):
    return 1/(1+np.exp(-x))
    
    
def takeLinearCombinations(data, startColumn=2) :
    '''
    startColumn means the column where features start. Default is 2, assuming
    that column 0 is the response value, column 1 is the bias term
    
    Returns the linear combination of features, as well as a matrix of the signs
    that were used to generate linear combinations (would probably need this and
    apply same transformation to test to do fair comparison). Each column in the
    matrix corresponds to the linear transformation used
    '''

    M = data.shape[1] - startColumn
    signs = np.empty((M,M))
    new_data = np.copy(data)
    for m in range(M) :
        a = np.random.random(data.shape[1]-startColumn)
        a[a>=0.5] = 1
        a[a<0.5] = -1
        new_data[:,startColumn+m] = np.sum(lumABData[:,2:]*a,axis=1)
        signs[:,m] = a
        
    return new_data,signs
    
def processFeatherInput(data,N=45,insertBias=True):
    '''
    Assumes that response value is binary for now
    '''
    new_data = data.iloc[:,0:N].values
#    y_vals = np.unique(new_data[:,0])
    if insertBias:
        new_data = np.concatenate((new_data[:,[0]], np.ones((len(new_data),1)),new_data[:,1:]),axis=1)
    new_data[:,0] = new_data[:,0] == new_data[0,0] #y_vals[0]
    new_data = new_data.astype(float)
    
    return new_data

def processFeatherClasses(class_labels, classes=None) :
    '''
    Takes a vector (or similar) of classes and returns a vector of categorical 
    integers. Used when original data comes in Feather format.
    
    Parameters 
    -----
    
    list_classes: a list of the classes
    
    classes: a list of the unique classes. For example, if want to use the 
        same order of classes found in the training data set as in testing 
    '''

    if classes is None:
        classes = np.unique(class_labels)
        
    cat_classes = np.zeros(len(class_labels))
    for ii,cl in enumerate(classes) :
        cat_classes[np.where(class_labels==cl)] = ii
    
    return cat_classes 
    
def thresholdProbs(y,th=0.5):
    new_y = np.copy(y)
    new_y[y>=th] = 1
    new_y[y<=th] = 0
    return new_y 
    
def timeRandsample(sols,data,order,N=20) :

    L =  np.min((sols.shape[0],N))
    objf = np.zeros(L)
    sol = sols[0]
    sol[sol==3] = 0
    count = 2
    objf[0] = get_acc(sol,data,order);
    for n in range(1,L):
        dwtmp = sols[n]
        dwtmp[dwtmp==3] = 0 ;
        tmpsol = sol + (dwtmp-sol)/float(count)
        tmp = get_acc(tmpsol,data,order)
        if tmp > objf[n-1] :
            objf[n] = tmp
            sol = tmpsol
            count += 1
        else :
            objf[n] = objf[n-1]
    return sol,objf

def get_acc(sol,data,order) :

    pred_probs = multinomial_probs(data[order,2:],sol);
    pred_classes = np.zeros(pred_probs.shape[0]);
    maxm = np.max(pred_probs,axis=1)
#    class_range = np.array(range(pred_probs.shape[1]))
    for n in range(len(pred_classes)):
        pred = pred_probs[n]
        
#        pred_classes[n] = np.random.choice(class_range[pred==maxm[n]])
        pred_classes[n] = np.random.choice(np.flatnonzero(pred==maxm[n]))
#       choices = np.flatnonzero(pred==maxm[n])
#        if len(choices) > 1:
#            pred_classes[n] = np.random.choice(choices)
#        else :
#            pred_classes[n] = choices[0]
    acc = np.count_nonzero(pred_classes == data[order,0])/len(pred_classes)

    return acc

def multinomial_probs(X,weights) :

    exptmpout = np.exp(np.dot(X,np.reshape(weights,(X.shape[1],-1),order='F'))); # check this later to make sure gives same result
    pred_class_probs = np.zeros((X.shape[0],exptmpout.shape[1]+1))
    denom = 1 + np.sum(exptmpout,axis=1)
    pred_class_probs[:,0:exptmpout.shape[1]] = exptmpout/denom[:,None]
    pred_class_probs[:,-1] = 1./denom
    
    return pred_class_probs
