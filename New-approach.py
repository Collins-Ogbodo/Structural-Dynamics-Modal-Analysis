import os,pickle,zipfile

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
        
        #Data dictionary label
        data_labels = [label[:-7] for label in list_files]
        data_labels = [label.replace("-", "_" ) for label in data_labels]
        
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
        sensor_frf_mean = {key:[] for key in data_labels}
        sensor_frf_freq_mean = {key:[] for key in data_labels}
        
    for sensors in sensor_frf_lists:
        freqs_frf = sensor_frf_freq_lists[sensors]
        frf = sensor_frf_lists[sensors]
        
        #find average of all iterations
        frf_average = [sum(frf_list) / len(frf_list) for frf_list in zip(*sensor_frf_lists[sensors])]
        frf_freq_average = [sum(frf_freq_list) / len(frf_freq_list) for frf_freq_list in zip(*sensor_frf_freq_lists[sensors])]
        sensor = sensors[:-7]
        sensor = sensor.replace('-','_')
        #Store average for this sensor in the empty dictionary
        sensor_frf_mean[sensor] = frf_average
        sensor_frf_freq_mean[sensor] = frf_freq_average

#%%
import numpy as np
import scipy as sp
def rational_polynomial_method(frf, fs, degree=2):
    """
    Estimate the natural frequency, damping factor, modal constant, and mode shape
    of a system from its complex frequency response function (FRF).
    
    Parameters
    ----------
    frf : array-like
        The complex frequency response function of the system.
    fs : float
        The sampling frequency of the FRF.
    degree : int, optional
        The degree of the rational polynomial approximation. Default is 2.
        
    Returns
    -------
    omega : float
        The natural frequency of the system.
    zeta : float
        The damping factor of the system.
    k : float
        The modal constant of the system.
    mode_shape : array-like
        The mode shape of the system.
    """
    # Convert the FRF to a rational polynomial
    r, p, k = sp.signal.residue(frf, np.ones(degree+1), tol=1e-6)
    
    # Find the poles of the rational polynomial
    poles = p[np.abs(k) > 1e-6]
    
    # Find the natural frequency and damping factor of the system
    omega = np.abs(poles[0])
    zeta = -np.real(poles[0]) / omega
    
    # Find the modal constant
    k = np.abs(k[0])
    
    # Find the mode shape
    mode_shape = r[0] / k
    
    return omega, zeta, k, mode_shape

#%%

def rational_polynomial(frf, n_modes=6):
    """
    Estimate natural frequencies, damping ratios, modal constants, and mode shapes
    using the rational polynomial method.
    
    Parameters:
    -----------
    frf : array_like
        Complex frequency response function.
    n_modes : int, optional
        Number of modes to estimate. Default is 1.
        
    Returns:
    --------
    freqs : ndarray
        Estimated natural frequencies.
    damping : ndarray
        Estimated damping ratios.
    constants : ndarray
        Estimated modal constants.
    shapes : ndarray
        Estimated mode shapes.
    """
    # Compute the poles and residues of the FRF using partial fraction expansion
    poles, residues = np.polynomial.polynomial.polypow(frf, m=n_modes)
    
    # Extract the natural frequencies and damping ratios from the poles
    freqs = np.abs(poles)
    damping = -np.real(poles) / freqs
    
    # Compute the modal constants and mode shapes
    constants = np.abs(residues)
    shapes = np.angle(residues)
    
    return freqs, damping, constants, shapes

#%%
import pandas as pd
#Convert average frf for all sensor to dataframe
data_frf = pd.DataFrame(sensor_frf_mean)   
data_frf = data_frf.to_numpy()

#%%

dof = 2
m= 2*dof-1; 
n= 2*dof;
Freq = {}
Freq['Freq'] = sensor_frf_freq_mean['EXH']
sensor_para = []
#df.columns[2:]

frf = data_frf[:,0]
fs = np.array(Freq['Freq'])
wn, zeta, phi, c = rational_polynomial_method(frf, fs)
modal_para = {}
modal_para['wn'] = wn
modal_para['zeta'] = zeta
modal_para['phi'] = phi
modal_para['c'] = c
sensor_para.append(modal_para)
#%%
dof = 2
m= 2*dof-1; 
n= 2*dof;
Freq = {}
Freq['Freq'] = sensor_frf_freq_mean['EXH']
sensor_para = []
#df.columns[2:]
for sensors in data_frf.columns[0:]:
    frf = np.array(sensors)
    fs = np.array(Freq['Freq'])
    wn, zeta, phi, c = rational_polynomial_method(frf, fs, n, m)
    modal_para = {}
    modal_para['wn'] = wn
    modal_para['zeta'] = zeta
    modal_para['phi'] = phi
    modal_para['c'] = c
    sensor_para.append(modal_para)
    
    
