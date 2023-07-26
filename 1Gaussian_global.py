import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lmfit import Parameters, minimize, report_fit
import math


def gauss2(x1, amp1, cen1, sigma1):          # This is the equation we're fitting to
    """Gaussian lineshape."""
    return amp1 * np.exp(-(x-cen1)**2 / (sigma1**2))


def gauss1(x, amp, cen, sigma):          # This is the equation we're fitting to
    """Gaussian lineshape."""
    return amp * np.exp(-(x-cen)**2 / (sigma**2))



def gauss_dataset(params, i, x):
    """Calculate Gaussian lineshape from parameters for data set."""
    amp1 = params['amp1_%i' % (i+1)]
    cen1 = params['cen1_%i' % (i+1)]
    sig1 = params['sig1_%i' % (i+1)]

    return gauss2(x, amp1, cen1, sig1)


def gauss_plot(params, i):
    """Calculate Gaussian lineshape from parameters for data set."""
    amp1 = params['amp1_%i' % (i+1)]
    cen1 = params['cen1_%i' % (i+1)]
    sig1 = params['sig1_%i' % (i+1)]
    
    return gauss1(x2, amp1, cen1, sig1)



def gauss2pl(x1, amp1, cen1, sigma1):          # This is the equation we're fitting to
    """Gaussian lineshape."""
    return amp1 * np.exp(-(x2-cen1)**2 / (sigma1**2))


def objective(params, x, data):
    """Calculate total residual for fits of Gaussians to several data sets."""
    ndata, _ = data.shape
    resid = 0.0*data[:]

    # make residual per data set
    for i in range(ndata):
        resid[i, :] = data[i, :] - gauss_dataset(params, i, x)

    # now flatten this to a 1D array, as minimize() needs
    return resid.flatten()


x2 = np.linspace(0.005, 1, 200)  
x = np.linspace(0.05, 1, 20)       # Make the 
file=r"/Users/Mathew/Documents/Current analysis/20230718_lysate/Timepoints/All_FRET_Small.csv"
datas = pd.read_table(file,header=0)
trans_datas=np.transpose(datas)

trans_datas.to_csv(file + '_' + 'trans.csv', sep = '\t')

data=trans_datas.to_numpy()
data=np.delete(data, 0, 0)


fit_params = Parameters()
for iy, y in enumerate(data):
    fit_params.add('amp1_%i' % (iy+1), value=50, min=0.0, max=100000)
    fit_params.add('cen1_%i' % (iy+1), value=0.35, min=0.2, max=0.5)
    fit_params.add('sig1_%i' % (iy+1), value=0.1, min=0.09, max=0.5)

    
    

for iy in range(2,len(data)+1):
    fit_params['cen1_%i' % iy].expr = 'cen1_1'
    fit_params['sig1_%i' % iy].expr = 'sig1_1'

    
out = minimize(objective, fit_params, args=(x, data))
report_fit(out.params)

font = {'family' : 'normal',
       
        'size'   : 20}

plt.rc('font', **font)

for i in range(len(data)):
    y_fit = gauss_plot(out.params, i)
    # plt.plot(x, data[i, :], 'o', color='#949494',markeredgecolor='black')
    plt.bar(x,data[i, :], align='center',width=0.05,facecolor='white',edgecolor='#949494')
    plt.plot(x2, y_fit, '-',color='#00008b')
    # plt.plot(x2, y_fit2, '-',color='#f7a813')
    # plt.plot(x2, y_fit, '-',color='#00008b')
    plt.xlabel('FRET Efficiency')
    plt.ylabel('Number of events')
    # plt.ylim((0,0.02))


    plt.savefig(file + '_' +str(i)+'.pdf')
    plt.show()


cen=out.params['cen1_1']
wid1=out.params['sig1_1']


print(cen,wid1)

low=[]
low_error=[]


for i in range(len(data)):
    amp1_name='amp1_'+str(i+1)
    amp1=out.params[amp1_name].value
    amp1_error=out.params[amp1_name].stderr
    
    width1_name='sig1_'+str(i+1)
    width1=out.params[width1_name].value
    width1_error=out.params[width1_name].stderr
    
    int1=width1*amp1*math.sqrt(math.pi)/0.05
    int1_error=math.sqrt(math.pi*((width1**2)*(amp1_error**2)+(amp1**2)*(width1_error**2)))/0.05
    low.append(int1)
    low_error.append(int1_error)
    
   
  
    
output=pd.DataFrame()

output['low']=low
output['low_error']=low_error


output.to_csv(file + '_' + 'Output.csv', sep = '\t')