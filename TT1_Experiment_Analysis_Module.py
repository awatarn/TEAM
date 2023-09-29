# TT1_Experiment_Analysis
# Created by Apiwat Wisitsorasak

import numpy as np
import matplotlib.pyplot as plt
from icecream import ic
import pandas as pd
from cycler import cycler

class TT1Discharge:
    def __init__(self, ShotNo, Path_Output):
        self.ShotNo = ShotNo
        self.Path = '%s/%s'%(Path_Output, ShotNo)
        self.exp = {}
        print('Shot Number: %d'%(self.ShotNo))
        print('Path: ', self.Path)
    
    def Read_0D_TN(self, variable, long_display=True):        
        if type(variable) != list:            
            var = [variable]
        else:
            var = variable
        
        for var_i in var:
            print('Loading variable: %s'%(var_i))
            path_data = '%s/%s.txt'%(self.Path, var_i)
            if long_display:
                print('Path to data: ', path_data)
    
            colnames = ['t', var_i]
            df = pd.read_csv(path_data, delim_whitespace=True, header=None, names=colnames, skiprows=8)
            self.exp[var_i] = df.copy()

            if 'HCN' in var_i: # convert HCN signal to density, 5V = 1.7e19 m^-3
                key = 'NE%d'%(int(var_i[3]))
                self.exp[key] = df.copy()
                self.exp[key].rename(columns={var_i:key}, inplace=True)
                self.exp[key][key] = self.exp[key][key] * (1.7e19 / 5.0)


def PlotSingleColumn(shots, var, **kwargs):
    
    # default keyword arguments
    defaultKwargs = {'lw':1, 'color':'k', 'tinit':0, 'tfinal':1000,
                     'grid':True, 'width':8, 'height':1.2, 'savefig':False}
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
        print(shot.ShotNo)

        for i in range(TotalVar): # Loop of variables            

            tinit_index = np.where(shot.exp[var[i]]['t'] >= kwargs['tinit'])[0][0]        
            if len( np.where(shot.exp[var[i]]['t'] >= kwargs['tfinal'])[0]) == 0:
                tfinal_index = -1
            else:
                tfinal_index = np.where(shot.exp[var[i]]['t'] >= kwargs['tfinal'])[0][0]

            x = shot.exp[var[i]]['t'][tinit_index:tfinal_index]
            y = shot.exp[var[i]][var[i]][tinit_index:tfinal_index]
            ax[i].plot(x, y, linewidth=kwargs['lw'], label='%d'%(shot.ShotNo), color=color[j])
            ax[i].set_ylabel(var[i])
                        
            ax[i].grid(kwargs['grid'])

            if (i==0) and (j == len(shots)-1):
                ax[i].legend(frameon=False)



    plt.xlabel('time [ms]')

    if kwargs['savefig']:
        suffix = "_".join([str(i.ShotNo) for i in shots])
        fn = 'DischargePlots_%s.png'%(suffix)
        print("Save figure to: %s"%(fn))
        fig.savefig(fn, dpi=200)

    plt.show()

    return fig