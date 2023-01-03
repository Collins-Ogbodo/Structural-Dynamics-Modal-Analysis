#Import python module
import os,pickle,zipfile
import pandas as pd
import numpy as np

#Import and process FRF data
iters = [10]
reps = [2]
test_series = "DS_RLE"

#empty list for individual sensor plot
sensor_frf_lists = {}
sensor_frf_freq_lists = {}
sensor_frf_mean = {}

for iter_ in iters:
    
    for rep in reps: 
        test_rep = rep
        test_name = test_series + '_' + str(iter_) + '_' + str(test_rep)
        zip_file_path = "Test_Data/" + test_series + "/"+ test_series \
            + "_" + str(iter_)+ '/' + test_series + "_" \
                + str(iter_) + "_" + str(rep) + '.zip'
        #directory for dumping the extraction file
        ext_dump_dir = 'Extraction_Dump'
        if not os.path.exists(ext_dump_dir):
             os.makedirs(ext_dump_dir)
        #unzip files 
        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            zip_ref.extractall(ext_dump_dir)

        #create a list of files
        list_files = os.listdir(ext_dump_dir)
        list_files.remove('meta_data.json')
        
        #Remove FRC which is an empty list
        list_files.remove('FRC.pickle')
        
        #create an empty list for frf and coh and their frequencies for each sensor
        if not sensor_frf_lists:
            sensor_frf_lists = {key:[] for key in list_files}
            sensor_frf_freq_lists = {key:[] for key in list_files}

        #create combine all repetition for each sensor
        for sensor in list_files:
            #open file
            infile = open(os.path.join(ext_dump_dir,sensor), 'rb')
            sensor_ = pickle.load(infile)
            infile.close()
            
            all_data = sensor_['data']
            #extract individual data
            frf= all_data[-1]
            
            #zip all data with corresponding frequencies
            frf_data, frf_freq = frf['data_y'],frf['data_x']
            
            #Append the zipped data in a list for each sensors
            sensor_frf_lists[sensor].append(frf_data)
            sensor_frf_freq_lists[sensor].append(frf_freq) 

            
    # Empty dic. for all iteration of individual sensor
    if not sensor_frf_mean:
        sensor_frf_mean = {key:[] for key in list_files}
        sensor_frf_freq_mean = {key:[] for key in list_files}
        
    for sensors in sensor_frf_lists:
        freqs_frf = sensor_frf_freq_lists[sensors]
        frf = sensor_frf_lists[sensors]
        
        #find average of all iterations
        frf_average = [sum(frf_list) / len(frf_list) for frf_list in zip(*sensor_frf_lists[sensors])]
        frf_freq_average = [sum(frf_freq_list) / len(frf_freq_list) for frf_freq_list in zip(*sensor_frf_freq_lists[sensors])]
        
        #Store average for this sensor in the empty dictionary
        sensor_frf_mean[sensors] = frf_average
        sensor_frf_freq_mean[sensors] = frf_freq_average

#%%

#Convert average frf for all sensor to dataframe
data_frf = pd.DataFrame(sensor_frf_mean)   
data_frf = data_frf.to_numpy()

#Convert average frequency for all sensor to dataframe
data_freq = pd.DataFrame(sensor_frf_freq_mean)   
data_freq = data_freq.to_numpy()
freq = np.array(data_freq[:,0])

#%%
#Order of degree of freedom
N = 5
#Normalisation of frequency
freq = (freq-min(freq))/(max(freq)-min(freq))

m = 2*N-1  #number of numerator polynomial terms
n = 2*N    #number of denominator polynomial terms

#Number of modes
num_mode = 2
#Power of fitting polynomial
max_k = 2* num_mode

#Generating Orthorgonal Polynomial
def orthogonal(frf, freq, max_k, thetha_phi):
    #This function generates the orthogonal polynomial
    #and returns the transformation matric and Polynomials
    
    #Calculating the weighting function at the i_th frequency
    if thetha_phi == 'num':
        q = np.ones(len(freq),len(freq))
    elif thetha_phi == 'denum':
        q = np.array(abs(data_frf)**2)
    #Calculating Polynomials
    R_minus_1 = np.zeros(len(freq),len(freq))
    R_0 = np.ones(len(freq),len(freq))*(1/np.sqrt(2*sum(q)))
    R=[R_minus_1, R_0]
    
    #computing the weighing function of the Rational Fraction
    for i in range(max_k):
        v_k = 2*sum(freq*R[i+1]*R[i]*q)
        s_k = freq*R[i+1]-v_k*R[i]
        d_k = np.sqrt(2*sum((s_k**2)*q))
        R_ = s_k/d_k
        R.append(R_)
    #Orthogonal Matrix
    R = R[1:-1]
    
    #Complex complex part of polynomial
    j_k = [np.sqrt(-1)]*max_k
    P = R*j_k
    
    return P
matrix_phi = orthogonal(frf, freq, m, 'num')
matrix_theta = orthogonal(frf, freq, n, 'denum')

    
        
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
