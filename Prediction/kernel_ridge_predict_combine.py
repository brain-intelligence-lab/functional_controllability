import numpy as np
import warnings
import scipy.io as sio
from scipy.stats import pearsonr
from sklearn.model_selection import cross_val_predict
from sklearn.kernel_ridge import KernelRidge

K = 10  
run_session = '1LR'
mykernel = 'cosine' # cosine  linear
alpha_list = [0.1,1,2,5,10,15,20,50,80,100,150,200,300,1000,10000,100000,1000000]
default_folder = 'save_measurement/'

warnings.filterwarnings('ignore')

sortedFeat = ['PicSeq_Unadj','CardSort_Unadj','Flanker_Unadj','PMAT24_A_CR','ReadEng_Unadj',\
              'PicVocab_Unadj','ProcSpeed_Unadj','DDisc_AUC_40K','VSPLOT_TC','SCPT_SEN','SCPT_SPEC',\
              'IWRD_TOT','ListSort_Unadj','MMSE_Score','PSQI_Score','Endurance_Unadj','GaitSpeed_Comp',\
              'Dexterity_Unadj','Strength_Unadj','Odor_Unadj','PainInterf_Tscore','Taste_Unadj','Mars_Final',\
              'Emotion_Task_Face_Acc','Language_Task_Math_Avg_Difficulty_Level',\
              'Language_Task_Story_Avg_Difficulty_Level','Relational_Task_Acc','Social_Task_Perc_Random',\
              'Social_Task_Perc_TOM','WM_Task_Acc','NEOFAC_A','NEOFAC_O','NEOFAC_C',\
              'NEOFAC_N','NEOFAC_E','ER40_CR','ER40ANG','ER40FEAR','ER40HAP','ER40NOE',\
              'ER40SAD','AngAffect_Unadj','AngHostil_Unadj','AngAggr_Unadj','FearAffect_Unadj',\
              'FearSomat_Unadj','Sadness_Unadj','LifeSatisf_Unadj','MeanPurp_Unadj','PosAffect_Unadj',\
              'Friendship_Unadj','Loneliness_Unadj','PercHostil_Unadj','PercReject_Unadj','EmotSupp_Unadj',\
              'InstruSupp_Unadj','PercStress_Unadj','SelfEff_Unadj']
              
def normalization(feature):
    [N,M] = feature.shape
    
    for i in range(N):
        avg = feature[i,:].mean()
        stt = feature[i,:].std()
        feature[i,:] = (feature[i,:]-avg)/stt
    return feature
    
def norm_task(task):
    [N,M] = task.shape 
    for i in range(M):
        task[:,i] = (task[:,i]-task[:,i].mean())/task[:,i].std()
    return task

def load_tasks():
    ls = 'tasks_'+run_session+'.mat'
    behavior = sio.loadmat(ls)
    tasks = behavior['behaviour_data']
    return tasks
    
def load_feature_control(session,item):
    load_name = default_folder+'Control_'+session+'_all'+run_session+'.mat'
    mat_system = sio.loadmat(load_name)
    feature_name = item+'M_'+session
    feature = mat_system[feature_name]
    feature = normalization(feature)
    return feature
    
def load_feature_graph(session,item):
    load_name = default_folder+'Graph_'+session+'_all'+run_session+'.mat'
    mat_system = sio.loadmat(load_name)
    feature_name = item+'_'+session
    feature = mat_system[feature_name]
    feature = normalization(feature)
    return feature
    
def cross_score(feature,task):
    [N,sub_num] = feature.shape
    data_X = feature.T
    data_y = task
    model = KernelRidge(alpha=alpha_list, kernel=mykernel)
    predicted = cross_val_predict(model, data_X, data_y, cv=K, n_jobs = K)
    M = task.shape[1]
    m_r = np.zeros(M)
    for i in range(M):
        try:
            m_r[i] = pearsonr(predicted[:,i],data_y[:,i])[0]
        except:
            pass # avoid NaN
        print(' -> ',sortedFeat[i],'corr:', m_r[i])
    return m_r

if __name__ == '__main__':
    session = ['FC','AR']
    
    job1 = 'control'
    job2 = 'graph'
    
    items1 = ['AC','MC','EV']
    items2 = ['Node','efficiency','Pcoef']
    
    tasks = load_tasks() # load all 58 tasks
    tasks = norm_task(tasks) # normalization

    
    max_correlation_combine = np.zeros([36,len(sortedFeat)])
    
    count = 0
    for i1 in session:
        for i2 in session:
            for j1 in items1:
                for j2 in items2:
                    feature_control = load_feature_control(i1, j1)
                    feature_graph = load_feature_graph(i2,j2)
                    feature_combine = np.row_stack((feature_control, feature_graph))
                    cor = cross_score(feature_combine,tasks)
                    max_correlation_combine[count,:] = cor 
                    count += 1
    max_correlation_combine = max_correlation_combine.max(0)
    sio.savemat('max_correlation_combine.mat',{'max_correlation_combine':max_correlation_combine})   
                    
        
        