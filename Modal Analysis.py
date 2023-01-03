import os,pickle,zipfile
import pandas as pd
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import signal
import pyoma as oma

#%%
#Plot for the Sensor data for each repeat and iteration
iters = [1]
reps = [1, 2]
test_series = "BR_AR"

#empty list for individual sensor plot
sensor_Autospec_lists = {}
sensor_Autospec_freq_lists = {}
sensor_Autospec_mean = {}

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
        list_files.pop(1)
        
        #create an empty list for frf and coh and their frequencies for each sensor
        if not sensor_Autospec_lists:
            sensor_Autospec_lists = {key:[] for key in list_files}
            sensor_Autospec_freq_lists = {key:[] for key in list_files}
        
        #create combine all repetition for each sensor
        for sensor in list_files:
            #open file
            infile = open(os.path.join(ext_dump_dir,sensor), 'rb')
            sensor_ = pickle.load(infile)
            infile.close()
            
            all_data = sensor_['data']
            Autospec= all_data[1]
            
            #zip all data with corresponding frequencies
            Autospec_data, Autospec_freq = Autospec['data_y'],Autospec['data_x']
            
            #Append the zipped data in a list for each sensors
            sensor_Autospec_lists[sensor].append(Autospec_data) 
            sensor_Autospec_freq_lists[sensor].append(Autospec_freq) 

    # Empty dic. for all iteration of individual sensor
    if not sensor_Autospec_mean:
        sensor_Autospec_mean = {key:[] for key in list_files}
        sensor_Autospec_freq_mean = {key:[] for key in list_files}
        
    for sensors in sensor_Autospec_lists:
        freqs_Autospec = sensor_Autospec_freq_lists[sensors]
        Autospec = sensor_Autospec_lists[sensors]
 
        #find average of all iterations
        Autospec_average = [sum(Autospec_list) / len(Autospec_list) for Autospec_list in zip(*sensor_Autospec_lists[sensors])]
        Autospec_freq_average = [sum(Autospec_freq_list) / len(Autospec_freq_list) for Autospec_freq_list in zip(*sensor_Autospec_freq_lists[sensors])]
        
        #Store average for this sensor in the empty dictionary
        sensor_Autospec_mean[sensors] = Autospec_average
        sensor_Autospec_freq_mean[sensors] = Autospec_freq_average
    
    #Empty dictionary for the PSD
    PSD = {}
    
    #Division function
    def zero_div(x, y):
        if y == 0.0:
            return 0.0
        else:
            return x/y
     #Extract the real part of the autopower spectrum   
    for sensor in sensor_Autospec_mean.keys():
        Autospec_mean_dumb = []
        for i in sensor_Autospec_mean[sensor]:
            real_ = i.real
            Autospec_mean_dumb.append(real_)
        sensor_Autospec_mean[sensor] = Autospec_mean_dumb
        
        #Power spectrum density PSD
        PSD[sensor] = [zero_div(PSDs,freq) for PSDs , freq in zip(sensor_Autospec_mean[sensor], sensor_Autospec_freq_mean[sensors])]
        
f =PSD
#%%
#Convert average frf data to dataframe format
data = pd.DataFrame(PSD)
data = data.to_numpy()
# Sampling frequency
fs = 512 # [Hz] Sampling Frequency
df = 0.0625
# Run FDD
FDD = oma.FDDsvp(data, fs, df)
fdd = FDD
#%%

# Define list/array with the peaks identified from the plot
FreQ = [153] # identified peaks

# Extract the modal properties
# extracting modal properties using standard FDD
Res_FDD = oma.FDDmodEX(FreQ, FDD[1])
# extracting modal properties using Enhanced-FDD
Res_EFDD = oma.EFDDmodEX(FreQ, FDD[1], method='EFDD')
# extracting modal properties using FSDD with additional input
# arguments to customize the return of function, e.g. return plot
Res_FSDD = oma.EFDDmodEX(FreQ, FDD[1], method='FSDD', npmax = 35, MAClim=0.95, plot=True)


