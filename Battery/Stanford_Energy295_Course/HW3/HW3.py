
import numpy as np  
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt   

# ----------------------------------------------------------------------
#   Main
# ---------------------------------------------------------------------- 
def main():  
    
    
    capacity_test_filename = 'Files/INR21700_M50T_T23_OCV_W8.xlsx'
    hppt_test_filename     = 'Files/INR21700_M50T_T23_HPPC_N0_W8.xlsx' #  hybrid pulse power characterization test filename
    validation_filename    = 'Files/INR21700_M50T_T23_UDDS_W8.xlsx'  
    
    battery = dict({'Manufacturer'    : 'LG Chem',
                   'Rated_Capacity'   : 4.85,
                   'Nominal_Voltage'  : 3.63,
                   'Maximum_Voltage'  : 4.2,
                   'Minimum_Voltage'  : 2.5,
                   'Cathode_Chemistry': 'NMC',
                   'Anode_Chemistry'  : 'Graphite'})
    
    capacity_test_raw_data , battery =  compute_nominal_capacity_and_OCV_vs_SOC(capacity_test_filename,battery)
    parameter_estimation(hppt_test_filename,battery,capacity_raw_test_data)
    return 


    
def compute_nominal_capacity_and_OCV_vs_SOC(capacity_test_filename,battery):
    
    # -------------------------------------------------------
    # read data  
    # -------------------------------------------------------
    
    # charging is positive current , discharging is negative current
    raw_data = pd.read_excel(capacity_test_filename)
    
    # -------------------------------------------------------
    # computations 
    # -------------------------------------------------------
    
    # integrate to current to get Ah
    Q     = sp.integrate.cumtrapz(y = raw_data['Current(A)'], x=raw_data['Test_Time(s)']) 
    Q_nom = Q[-1]
    
    # compute SOC as a function of voltage
    SOC   = Q/Q_nom
    
    # append nominal capacity to dictionary
    battery['Nominal_Capacity'] = Q_nom/3600
    
    # append SOC onto dataframe
    raw_data['SOC']             = SOC
    raw_data['Charge(Ah)']     = Q/3600 
    
    # -------------------------------------------------------
    # plots  
    # -------------------------------------------------------
    fig = plt.figure()        
    fig.set_size_inches(10,10)  
    
    # Voltage vs Time  
    axis_1 = fig.add_subplot(2,2,1)  
    axis_1.plot(raw_data['Test_Time(s)'] ,raw_data['Voltage(V)'] , 'k-')
    axis_1.set_xlabel('Time (s)')
    axis_1.set_ylabel('Voltage (OCV)')    
    
    # Current vs Time 
    axis_2 = fig.add_subplot(2,2,2)  
    axis_2.plot(raw_data['Test_Time(s)'] ,raw_data['Current(A)'] , 'k-')
    axis_2.set_xlabel('Time (s)')
    axis_2.set_ylabel('Current (A)')   
    
    # Charge vs Time 
    axis_3 = fig.add_subplot(2,2,3)  
    axis_3.plot(raw_data['Test_Time(s)'] ,raw_data['Charge(Ah)'] , 'k-')
    axis_3.set_xlabel('Time (s)')
    axis_3.set_ylabel('Charge (Ah)')  
    
    # SOC vs V 
    axis_4 = fig.add_subplot(2,2,4)  
    axis_4.plot(raw_data['SOC'] ,raw_data['Voltage (OCV)'] , 'k-')
    axis_4.set_xlabel('SOC')
    axis_4.set_ylabel('Voltage (OCV)')   
    
    return  raw_data , battery

def parameter_estimation(hppt_test_filename,battery,capacity_raw_test_data):
    # read data  
    # charging is positive current , discharging is negative current
    hppt_test_raw_data = pd.read_excel(hppt_test_filename)
    
    return 


if __name__ == '__main__': 
    main()    
    plt.show()   