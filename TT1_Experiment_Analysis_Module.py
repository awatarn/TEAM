# TT1_Experiment_Analysis
# Created by Apiwat Wisitsorasak on 29 Sep 2023
# - Update many part: PrintSummary, Compute Bt, NeMax, IpMax [1 Oct 2023]

import numpy as np
import matplotlib.pyplot as plt
from icecream import ic
import pandas as pd
from cycler import cycler
from scipy.interpolate import interp1d

class TT1Discharge:
    # global variables
    _aa = 0.20              # minor radius [m]
    _bb = 0.25              # radius of vacuum chamber [m]
    _RR = 0.65              # major radius [m]
    _TF_coil_turn = 40      # number of turns in each TF coil
    _TF_coil_rin = 0.33     # inner radius of TF coil [m]
    _TF_coil_rout = 0.44    # outer radius of TF coil [m]
    _TF_coil_r = (0.33 + 0.44) / 2   # average radius of TF coil [m]

    def __init__(self, ShotNo, Path_Output):
        self.ShotNo = ShotNo
        self.Path = '%s/%s'%(Path_Output, ShotNo)
        self.exp = {}
        self.DischargeTime = -1         # Discharge time interval
        self.DischargeTime0 = -1        # beginning time of the discharge, measure from IP
        self.DischargeTime1 = -1        # final time of the discharge, measure from IP
        self.Bt0Max = -1                # maximum of Bt at center (95 percentile)
        self.NeBar = -1                 # average values of the density (Ne)
        self.NeMax = -1                 # maximum value of Ne (95 percentile)
        self.IpMax = -1                 # maximum value of the plasma current (95 percentile)
        self.Date = 'None'              # Date and time of experiment
        print('Shot Number: %d'%(self.ShotNo))
        print('Path: ', self.Path)

    def Print_DischargeSummary(self, OutputMode=0):
        # Print summary of the discharge
        # OutputMode = 0 (separate by a comma, for Excel)
        # OutputMode = 1 (pretty print-out)
        print('-----------------------------------------------------------')
        print('                            Summary                        ')
        print('-----------------------------------------------------------')
        item = 'Shot Number'
        value = self.ShotNo
        if OutputMode == 0:
            print('{0}, {1}'.format(item, value))
        elif OutputMode == 1:
            print('{0:<40}, {1}'.format(item, value))

        item = 'Date and time'
        value = self.Date
        if OutputMode == 0:
            print('{0}, {1}'.format(item, value))
        elif OutputMode == 1:
            print('{0:<40}, {1}'.format(item, value))  

        item = 'Discharge time [ms]'
        value = self.DischargeTime
        if OutputMode == 0:
            print('{0}, {1}'.format(item, value))
        elif OutputMode == 1:
            print('{0:<40}, {1:.4}'.format(item, value))            

        item = 'Discharge time (begin) [ms]'
        value = self.DischargeTime0
        if OutputMode == 0:
            print('{0}, {1}'.format(item, value))
        elif OutputMode == 1:
            print('{0:<40}, {1:.4}'.format(item, value)) 

        item = 'Discharge time (end) [ms]'
        value = self.DischargeTime1
        if OutputMode == 0:
            print('{0}, {1}'.format(item, value))
        elif OutputMode == 1:
            print('{0:<40}, {1:.4}'.format(item, value))  

        item = 'Bt0 (max, average > 95%) [T]'
        value = self.Bt0Max
        if OutputMode == 0:
            print('{0}, {1}'.format(item, value))
        elif OutputMode == 1:
            print('{0:<40}, {1:.4}'.format(item, value))  

        item = 'Ip (max, average > 95%) [kA]'
        value = self.IpMax/1000
        if OutputMode == 0:
            print('{0}, {1}'.format(item, value))
        elif OutputMode == 1:
            print('{0:<40}, {1:.4}'.format(item, value))             

        item = 'Ne (max, average > 95%) [m**-3]'
        value = self.NeMax
        if OutputMode == 0:
            print('{0}, {1}'.format(item, value))
        elif OutputMode == 1:
            print('{0:<40}, {1:.4}'.format(item, value))

        item = 'Ne (average) [m**-3]'
        value = self.NeBar
        if OutputMode == 0:
            print('{0}, {1}'.format(item, value))
        elif OutputMode == 1:
            print('{0:<40}, {1:.4}'.format(item, value))

    def Compute_DischargeTime(self, var_i='IP1', ThresholdPercentage=5.0):
        # Compute the discharge time from IP
        # ThresholdPercentage = minimum threshold value of the current [percentage]

        ThresholdValue = (ThresholdPercentage/100)*np.max(self.exp[var_i][var_i])
        
        # dt = self.exp[var_i]['t'][10] - self.exp[var_i]['t'][9]
        # DischargeTime = dt * np.sum(self.exp[var_i][var_i] > ThresholdValue)
        # self.DischargeTime = DischargeTime

        t0_index = np.where(self.exp[var_i][var_i] > ThresholdValue)[0][0]
        t1_index = np.where(self.exp[var_i][var_i] > ThresholdValue)[0][-1]

        DischargeTime = self.exp[var_i]['t'][t1_index] - self.exp[var_i]['t'][t0_index]
        self.DischargeTime = DischargeTime

        self.DischargeTime0  = self.exp[var_i]['t'][t0_index]
        self.DischargeTime1  = self.exp[var_i]['t'][t1_index]

    def Compute_Ne(self, var_i):
        # convert HCN signal to density, 5V = 1.7e19 m^-3
        key = 'NE%d'%(int(var_i[3]))
        self.exp[key] = self.exp[var_i].__copy__()
        self.exp[key].rename(columns={var_i:key}, inplace=True)
        self.exp[key][key] = self.exp[key][key] * (1.7e19 / 5.0)

        # Compute the average density, must be done after Compute_DischargeTime()
        if self.DischargeTime > 0:
            t0_index = np.where(self.exp[key]['t'] > self.DischargeTime0)[0][0]
            t1_index = np.where(self.exp[key]['t'] > self.DischargeTime1)[0][0]
            self.NeBar = np.average(self.exp[key][key][t0_index:t1_index])

            NeMax = np.percentile(self.exp[key][key], 95)
            temp_index = np.where(self.exp[key][key] > NeMax)[0]
            self.NeMax = np.average(self.exp[key][key][temp_index])
        else:
            print('[Warning] Cannot compute the average density since DischargeTime < 0')
            self.NeBar = -1
    
    def Compute_Bt0(self, var_i):
        # convert IT1 or IT2 to Bt (center)
        # From the calculation of Bt in Shimizu's code (based on Biot-Savart law),
        # if current = 0.26e6/40 = , Bt (center) = 1.2820298031567656.
        # Bt is proportional to Ip, Therefore, Bt = 1.2820298031567656 / 6500 * IT1
        
        key = 'BT0'
        self.exp[key] = self.exp[var_i].__copy__()
        self.exp[key].rename(columns={var_i:key}, inplace=True)
        self.exp[key][key] = np.abs(self.exp[key][key] * (1.2820298031567656 / 6500.0))

        Bt095 = np.percentile(self.exp[key][key], 95)
        temp_index = np.where(self.exp[key][key] > Bt095)[0]
        self.Bt0Max = np.average(self.exp[key][key][temp_index])

    def Get_unit(self, var_i):
        # Get unit of a variable
        
        if 'IP' in var_i:
            unit = 'A'
        elif 'HA' in var_i:
            unit = 'V'
        elif 'HA' in var_i:
            unit = 'V'
        elif 'HCN' in var_i:
            unit = 'V'
        elif 'BT' in var_i:
            unit = 'T'
        elif 'NE' in var_i:
            unit = 'm**-3'
        elif 'VP' in var_i:
            unit = 'V'
        elif 'IT' in var_i:
            unit = 'A'
        elif 'HA' in var_i:
            unit = 'V'
        elif 'DIA' in var_i:
            unit = 'Wb'
        elif 'GP' in var_i:
            unit = 'V'
        else:
            unit = 'A.U.'
        
        return unit

    
    def Read_0D_TN(self, variable, long_display=True):        
        if type(variable) != list:            
            var = [variable]
        else:
            var = variable
        
        for var_i in var:            
            print('Loading variable: %10s'%(var_i), end='    ')
            path_data = '%s/%s.txt'%(self.Path, var_i)
            if long_display:
                print('Path to data: ', path_data)

            # Get Unit and date of experiment, just the first variable
            if var_i == var[0]:
                fp_i = open(path_data, 'r')
                j = 0
                while j < 10:
                    line = fp_i.readline()
                    if "Create" in line:
                        i0 = line.find('=')
                        self.Date = line[i0+2:-2]
                        continue
                    j = j + 1
    
            # Load data
            colnames = ['t', var_i]
            df = pd.read_csv(path_data, delim_whitespace=True, header=None, names=colnames, skiprows=8)
            self.exp[var_i] = df.copy()

            if 'HCN' in var_i: # convert to density
                self.Compute_Ne(var_i)
            
            if ('IT1' in var_i) or ('IT2' in var_i): # Compute the toroidal field
                self.Compute_Bt0(var_i)
            
            if 'IP' in var_i:
                self.Compute_DischargeTime(var_i, ThresholdPercentage=5.0)

                # Compute IpMax
                Ip95 = np.percentile(self.exp[var_i][var_i], 95)
                temp_index = np.where(self.exp[var_i][var_i] > Ip95)[0]
                self.IpMax = np.average(self.exp[var_i][var_i][temp_index])

            
    





def PlotSingleColumn(shots, var, **kwargs):
    # Plot experimental results of TT-1 discharge
    # - Add interpolation, PrintSummary, Unit [1 Oct 2023]    

    
    # default keyword arguments
    defaultKwargs = {'lw':1, 'color':'k', 'tinit':0, 'tfinal':1000,
                     'grid':True, 'width':8, 'height':1.2, 'savefig':False,
                     'unit':True, 'Interpolation':False, 'InterpolationPeriod':1,
                     'ShowSummary':True}
    kwargs = { **defaultKwargs, **kwargs }
    color=['r', 'g', 'b', 'c', 'm', 'y', 'k', 'r', 'g', 'b', 'c', 'm', 'y', 'k']

    if type(shots) != list:
        # if the input is not a list, convert it to a list
        shots = [shots]

    TotalDischarge = len(shots)
    TotalVar = len(var)

    print("Total number of discharges: %d"%(TotalDischarge))
    for shot in shots:
        print('%d, '%(shot.ShotNo), end='')
    print(' ')
    print("Total variables to plot: %d"%(TotalVar))
    for i in var:
        print('%s, '%(i), end='')
    print(' ')

    fig, ax = plt.subplots(TotalVar, 1, sharex=True)

    # set figure size
    fig.set_figheight(kwargs['height']*TotalVar)
    fig.set_figwidth(kwargs['width'])
    
    for j in range(len(shots)): # Loop of different shot numbers
        shot = shots[j]
        # print(shot.ShotNo)

        for i in range(TotalVar): # Loop of variables            

            tinit_index = np.where(shot.exp[var[i]]['t'] >= kwargs['tinit'])[0][0]        
            if len( np.where(shot.exp[var[i]]['t'] >= kwargs['tfinal'])[0]) == 0:
                tfinal_index = -1
            else:
                tfinal_index = np.where(shot.exp[var[i]]['t'] >= kwargs['tfinal'])[0][0]

            x = shot.exp[var[i]]['t'][tinit_index:tfinal_index]
            y = shot.exp[var[i]][var[i]][tinit_index:tfinal_index]

            if kwargs['Interpolation']: # interpolate the data with given InterpolationPeriod
                if tinit_index == 0: #  might cause an error in the interp1d
                    tinit_index = 1
                if tfinal_index == -1: #  might cause an error in the interp1d
                    tfinal_index = len(shot.exp[var[i]]['t'])-1

                x_interp = np.arange(shot.exp[var[i]]['t'][tinit_index], shot.exp[var[i]]['t'][tfinal_index], kwargs['InterpolationPeriod'])
                interp_func = interp1d(x, y)
                x = x_interp
                y = interp_func(x)

            # Main plot
            ax[i].plot(x, y, linewidth=kwargs['lw'], label='%d'%(shot.ShotNo), color=color[j])

            # ylim
            if 'GP' in var[i]:
                ax[i].set_ylim(-1, 5)

            # Unit
            if kwargs['unit']:
                ax[i].set_ylabel('%s [%s]'%(var[i], shot.Get_unit(var[i])))
            else:
                ax[i].set_ylabel(var[i])

            # Grid      
            ax[i].grid(kwargs['grid'])

            if (i==0) and (j == len(shots)-1):
                ax[i].legend(frameon=False)

            if (i == 0) and (kwargs['ShowSummary']):
                shot.Print_DischargeSummary(OutputMode=1)

            



    plt.xlabel('time [ms]')

    if kwargs['savefig']:
        suffix = "_".join([str(i.ShotNo) for i in shots])
        fn = 'DischargePlots_%s.png'%(suffix)
        print("Save figure to: %s"%(fn))
        fig.savefig(fn, dpi=200)

    plt.show()

    return fig