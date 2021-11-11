
import numpy as np  
from scipy.integrate    import  cumtrapz
import pandas as pd
from scipy.optimize import differential_evolution 
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
    parameter_estimation(hppt_test_filename,battery,capacity_test_raw_data)
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
    Q     = np.zeros(len(raw_data['Current(A)']))
    Q[1:] = cumtrapz(y = -raw_data['Current(A)'], x=raw_data['Test_Time(s)']) 
    Q_nom = Q[-1]
    
    # compute SOC as a function of voltage
    SOC = Q/Q_nom   
    # appen nominal capacity to dictionary
    battery['Nominal_Capacity'] = Q_nom/3600
    
    # append SOC onto dataframe
    raw_data['SOC']             = SOC[::-1]
    raw_data['Charge(Ah)']      = Q/3600 
    
    # -------------------------------------------------------
    # plots  
    # -------------------------------------------------------
    fig = plt.figure()        
    fig.set_size_inches(10,10)  
    
    # Voltage vs Time  
    axis_1 = fig.add_subplot(2,2,1)  
    axis_1.plot(raw_data['Test_Time(s)']/3600 ,raw_data['Voltage(V)'] , 'k-')
    axis_1.set_xlabel('Time (hrs)')
    axis_1.set_ylabel('Voltage (OCV)')    
    
    # Current vs Time 
    axis_2 = fig.add_subplot(2,2,2)  
    axis_2.plot(raw_data['Test_Time(s)']/3600 ,-raw_data['Current(A)'] , 'k-')
    axis_2.set_xlabel('Time (hrs)')
    axis_2.set_ylim(0,0.3)
    axis_2.set_ylabel('Current (A)')   
    
    # Charge vs Time 
    axis_3 = fig.add_subplot(2,2,3)  
    axis_3.plot(raw_data['Test_Time(s)']/3600 ,raw_data['Charge(Ah)'] , 'k-')
    axis_3.set_xlabel('Time (hrs)')
    axis_3.set_ylabel('Charge (Ah)')  
    
    # SOC vs V 
    axis_4 = fig.add_subplot(2,2,4)  
    axis_4.plot(raw_data['SOC'] ,raw_data['Voltage(V)'] , 'k-')
    axis_4.set_xlabel('SOC')
    axis_4.set_ylabel('Voltage (OCV)')   
    
    return  raw_data , battery

def parameter_estimation(hppt_test_filename,battery,capacity_raw_test_data):
    # read data  
    # charging is positive current , discharging is negative current
    hppt_test_raw_data = pd.read_excel(hppt_test_filename)
    
    time             = np.array(hppt_test_raw_data['Test_Time(s)'])
    voltage          = np.array(hppt_test_raw_data['Voltage(V)'])
    current          = np.array(hppt_test_raw_data['Current(A)'])
    cell_capacity    = battery['Nominal_Capacity']  

    
    R0_discharge_guess ,R1_discharge_guess ,C1_discharge_guess, R2_discharge_guess, \
        R2_discharge_guess, R0_charge_guess,R1_charge_guess,R2_charge_guess , R2_charge_guess\
        = graph_estimation(time,voltage,current,cell_capacity)
    
    
    
 
 
 
 
    
    number_of_pulses = 9
     
    # find first index instance where charging occures i.e.  current > 0
    charging_start_idx = np.where(current>0)[0][0]
    
    # find where current cuts off again i.e. end of charging 
    charging_end_idx   = np.where(current[charging_start_idx:] ==0)[0][0] + charging_start_idx
    
    
    pre_pulse_buffer     = 5
    post_pulse_buffer    = 35    
    trim_idx             = charging_end_idx
    SOC_0                = 1
     
    
    R0_discharge_guess = np.zeros(number_of_pulses)
    R1_discharge_guess = np.zeros(number_of_pulses)
    C1_discharge_guess = np.zeros(number_of_pulses)
    R2_discharge_guess = np.zeros(number_of_pulses)
    C2_discharge_guess = np.zeros(number_of_pulses)
    R0_charge_guess   = np.zeros(number_of_pulses)
    R1_charge_guess   = np.zeros(number_of_pulses)
    C1_charge_guess   = np.zeros(number_of_pulses)
    R2_charge_guess   = np.zeros(number_of_pulses)
    C2_charge_guess   = np.zeros(number_of_pulses)
    
    # trim data 
    time_trimmed_1       = time[trim_idx:]
    voltage_trimmed_1    = voltage[trim_idx:]
    current_trimmed_1    = current[trim_idx:] 
    
    for i in range(number_of_pulses):
        
        # -----------------------
        # DISCHARGE
        # -----------------------        
        # find indices  where discharge first occures 
        discharge_idx_start  = np.where(current_trimmed_1<0)[0][0]
        discharge_idx_end    = np.where(current_trimmed_1[discharge_idx_start:]==0)[0][0] + discharge_idx_start
        
        # compute SOC 
        discharge_time = time_trimmed_1[discharge_idx_start:discharge_idx_end] - time_trimmed_1[discharge_idx_start]
        SOC            = np.zeros(len(discharge_time))
        SOC[0]         = SOC_0
        SOC[1:]        = SOC_0 + cumtrapz(current_trimmed_1[discharge_idx_start:discharge_idx_end],discharge_time)/(cell_capacity*3600)
        
        # trim data 
        time_trimmed_2       = time_trimmed_1[discharge_idx_end:]
        voltage_trimmed_2    = voltage_trimmed_1[discharge_idx_end:]
        current_trimmed_2    = current_trimmed_1[discharge_idx_end:]
        
        # -----------------------
        # DISCHARGE PULSE 
        # -----------------------
        # find indices where discharge pulse occures 
        discharge_pulse_idx_start  = np.where(current_trimmed_2<0)[0][0]
        discharge_pulse_idx_end    = np.where(current_trimmed_2[discharge_pulse_idx_start:]==0)[0][0] + discharge_pulse_idx_start
        I0_discharge_pulse               = -np.mean(current_trimmed_2[discharge_pulse_idx_start:discharge_pulse_idx_end])
        V4_discharge_pulse               = voltage_trimmed_2[discharge_pulse_idx_start-pre_pulse_buffer]       
        V3_discharge_pulse               = voltage_trimmed_2[discharge_pulse_idx_end+post_pulse_buffer]
        V2_discharge_pulse               = voltage_trimmed_2[discharge_pulse_idx_end+1]
        delta_V_discharge_pulse          = V3_discharge_pulse-V2_discharge_pulse
        V_tau_discharge_pulse            = 0.63*delta_V_discharge_pulse + V2_discharge_pulse 
        tau_CT_idx_discharge_pulse       = np.argmin(abs(voltage_trimmed_2[discharge_pulse_idx_end:(discharge_pulse_idx_end+post_pulse_buffer)]-V_tau_discharge_pulse))
        tau_CT_discharge_pulse           = time_trimmed_2[discharge_pulse_idx_end:discharge_pulse_idx_end+post_pulse_buffer][tau_CT_idx_discharge_pulse] - time_trimmed_2[discharge_pulse_idx_end]
        V1_discharge_pulse              = voltage_trimmed_2[discharge_pulse_idx_end-1]
        
        # find where time constant occures   
        R0_discharge_pulse_guess         = abs(V4_discharge_pulse -V2_discharge_pulse)/I0_discharge_pulse
        R1_discharge_pulse_guess         = V4_discharge_pulse-V1_discharge_pulse-R0_discharge_pulse_guess
        C1_discharge_pulse_guess         = tau_CT_discharge_pulse/R1_discharge_pulse_guess
        R2_discharge_pulse_guess         = R1_discharge_pulse_guess*0.1
        C2_discharge_pulse_guess         = C1_discharge_pulse_guess*0.1
        
        discharge_pulse_voltage    = voltage_trimmed_2[discharge_pulse_idx_start:discharge_pulse_idx_end]
        discharge_pulse_time       = time_trimmed_2[discharge_pulse_idx_start:discharge_pulse_idx_end] - time_trimmed_1[discharge_pulse_idx_start]
        discharge_pulse_SOC        = np.zeros(len(discharge_pulse_time))
        discharge_pulse_SOC[0]     = SOC_0
        discharge_pulse_SOC[1:]    = SOC_0 + cumtrapz(current_trimmed_2[discharge_pulse_idx_start:discharge_pulse_idx_end],discharge_pulse_time)/(cell_capacity*3600)
        
        ECM_guess_parameters       = dict({'R0_guess' : R0_discharge_pulse_guess,
                                           'R1_guess' : R1_discharge_pulse_guess,
                                           'C1_guess' : C1_discharge_pulse_guess,
                                           'R2_guess' : R2_discharge_pulse_guess,
                                           'C2_guess' : C2_discharge_pulse_guess})
        
        # find optimal parameters 
        ECM_optimal_parameters      = optimize_ECM_parameters(battery,discharge_pulse_SOC,discharge_pulse_voltage,ECM_guess_parameters)  
        
        # store optimal parameters
        
        
        
        # -----------------------
        # CHARGE PULSE         
        # -----------------------     

        # find indices where charge pulse occures 
        charge_pulse_idx_start  = np.where(current_trimmed_2>0)[0][0]
        charge_pulse_idx_end    = np.where(current_trimmed_2[charge_pulse_idx_start:]==0)[0][0] + charge_pulse_idx_start
        I0_charge_pulse               = -np.mean(current_trimmed_2[charge_pulse_idx_start:charge_pulse_idx_end])
        V4_charge_pulse               = voltage_trimmed_2[charge_pulse_idx_start-pre_pulse_buffer]       
        V3_charge_pulse               = voltage_trimmed_2[charge_pulse_idx_end+post_pulse_buffer]
        V2_charge_pulse               = voltage_trimmed_2[charge_pulse_idx_end+1]
        delta_V_charge_pulse          = V3_charge_pulse-V2_charge_pulse
        V_tau_charge_pulse            = 0.63*delta_V_charge_pulse + V2_charge_pulse 
        tau_CT_idx_charge_pulse       = np.argmin(abs(voltage_trimmed_2[charge_pulse_idx_end:(charge_pulse_idx_end+post_pulse_buffer)]-V_tau_charge_pulse))
        tau_CT_charge_pulse           = time_trimmed_2[charge_pulse_idx_end:charge_pulse_idx_end+post_pulse_buffer][tau_CT_idx_charge_pulse] - time_trimmed_2[charge_pulse_idx_end]
        V1_charge_pulse              = voltage_trimmed_2[charge_pulse_idx_end-1]
        
        # find where time constant occures   
        R0_charge_pulse_guess         = abs(V4_charge_pulse -V2_charge_pulse)/I0_charge_pulse
        R1_charge_pulse_guess         = V4_charge_pulse-V1_charge_pulse-R0_charge_pulse_guess
        C1_charge_pulse_guess         = tau_CT_charge_pulse/R1_charge_pulse_guess
        R2_charge_pulse_guess         = R1_charge_pulse_guess*0.1
        C2_charge_pulse_guess         = C1_charge_pulse_guess*0.1
        
        charge_pulse_voltage    = voltage_trimmed_2[charge_pulse_idx_start:charge_pulse_idx_end]
        charge_pulse_time       = time_trimmed_2[charge_pulse_idx_start:charge_pulse_idx_end] - time_trimmed_1[charge_pulse_idx_start]
        charge_pulse_SOC        = np.zeros(len(charge_pulse_time))
        charge_pulse_SOC[0]     = SOC_0
        charge_pulse_SOC[1:]    = SOC_0 + cumtrapz(current_trimmed_2[charge_pulse_idx_start:charge_pulse_idx_end],charge_pulse_time)/(cell_capacity*3600)
        
        ECM_guess_parameters       = dict({'R0_guess' : R0_charge_pulse_guess,
                                           'R1_guess' : R1_charge_pulse_guess,
                                           'C1_guess' : C1_charge_pulse_guess,
                                           'R2_guess' : R2_charge_pulse_guess,
                                           'C2_guess' : C2_charge_pulse_guess})
        
        # find optimal parameters 
        ECM_optimal_parameters      = optimize_ECM_parameters(battery,charge_pulse_SOC,charge_pulse_voltage,ECM_guess_parameters)  
        
        
        # ------------------------------------------------
        # Trim Data for next discharge-charge pulse
        # ------------------------------------------------
        trim_idx             = np.maximum(discharge_pulse_idx_end,charge_pulse_idx_end)
        time_trimmed_1       = time_trimmed_2[trim_idx:]
        voltage_trimmed_1    = voltage_trimmed_2[trim_idx:]
        current_trimmed_1    = current_trimmed_2[trim_idx:]   
        
        SOC_0 = SOC[-1]
        
    return  R0_discharge_guess ,R1_discharge_guess ,C1_discharge_guess, R2_discharge_guess, R2_discharge_guess, R0_charge_guess,R1_charge_guess,R2_charge_guess , R2_charge_guess\

def optimize_ECM_parameters(battery,pulse_SOC,pulse_voltage,ECM_guess_parameters):
    
    # differential evolution 
    bounds = [ ( ECM_guess_parameters['R0_guess']*0.9 , ECM_guess_parameters['R0_guess']*1.1  ),
               ( ECM_guess_parameters['R1_guess']*0.9 , ECM_guess_parameters['R1_guess']*1.1 ),
               ( ECM_guess_parameters['C1_guess']*0.9 , ECM_guess_parameters['C1_guess']*1.1  ),
               ( ECM_guess_parameters['R2_guess']*0.9 , ECM_guess_parameters['R2_guess']*1.1   ),
               ( ECM_guess_parameters['C2_guess']*0.9 , ECM_guess_parameters['C2_guess']*1.1  )]
    
    results = differential_evolution(randels_model(args = pulse_SOC,pulse_voltage), bounds)
    
    ECM_opt_parameters = dict({'R0_opt' : results[0],
                               'R1_opt' : results[1],
                               'C1_opt' : results[2],
                               'R2_opt' : results[3],
                               'C2_opt' : results[4]})  
    
    return ECM_opt_parameters 

def randels_model(x, SOC, voltage): 
    R0 = x[0]  
    R1 = x[1]  
    C1 = x[2]  
    R2 = x[3]  
    C2 = x[4]   
     
    V_measured - V_model
    return 


if __name__ == '__main__': 
    main()    
    plt.show()   