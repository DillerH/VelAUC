# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 15:00:13 2024
First Derivative 
@author: andrew.hale736
"""
import pandas as pd
import matplotlib.pyplot as plt
import scipy.signal as ss
from numpy import log
from scipy import stats
import scipy.constants as sc

gCrit = {3:1.15, 4:1.48, 5:1.71, 6:1.89, 7:2.02, 8:2.13,\
         9:2.21, 10:2.29, 11:2.34, 12:2.41, 13:2.46, 14:2.51,\
             15:2.55, 16:2.59, 17:2.62, 18:2.65, 19:2.68, 20:2.71,\
                 21:2.73, 22:2.76, 23:2.78, 24:2.80, 25:2.82, 26:2.84,\
                     27:2.86, 28:2.88, 29:2.89, 30:2.91, 31:2.92, 32:2.94,\
                         33:2.95, 34:2.97, 35:2.98, 36:2.99, 37:3.00, 38:3.01,\
                             39:3.03, 40:3.04, 41:3.05, 42:3.06, 43:3.07, 44:3.08,\
                                 45:3.09, 46:3.10, 47:3.11, 48:3.12, 49:3.13, 50:3.13,\
                                     51:3.14, 52:3.15, 53:3.15, 54:3.16, 55:3.17, 56:3.18,\
                                         57:3.18, 58:3.19, 59:3.195, 60:3.20, 61:3.20, 62:3.21,\
                                             63:3.21, 64:3.22, 64:3.22, 65:3.23, 66:3.23, 67:3.24,\
                                                 68:3.24, 69:3.25, 70:3.26, 71:3.26, 72:3.27, 73:3.27,\
                                                     74:3.28, 75:3.28, 76:3.29, 77:3.29, 78:3.30, 79:3.30,\
                                                         80:3.31}

"""
make it calculate D, MW, and s
make some plots like are in the book
separate loop for group plot goes through data a second time
"""
def compileName(fileName):
    strlist = []
    fileNameList = []
    wlist = []
    w2tlist = []
    timelist = []
    templist = []
    Klist = []
    wconversion = 0.10471975512
    
    hnd=open(fileName,'r')
    while True:
        buff = hnd.readline()
        if not buff:
            break
        try:
            strlist.append(buff.strip())
        except ValueError:
            pass
    hnd.close()
    i = 9
    while i < len(strlist):
        templist = strlist[i].split()
        fileNameList.append(templist[9])
        i+=3
    i = 10
    while i < len(strlist):
        templist = strlist[i].split()
        wlist.append(wconversion * int(templist[1].replace(",","")))
        timelist.append(int(templist[3].replace(",","")) - 1157)
        Klist.append(273.15 + float(templist[5].replace(",","")))
        w2tlist.append(float(templist[7]))
        i+=3
    inputlist = list(zip(fileNameList, timelist, Klist, wlist, w2tlist))
    df_out = pd.DataFrame(inputlist, columns=['FileName','time','temp','w','w2t'])
    df_out['xbar'] = 0 
    df_out['xmen'] = 0
    df_out['dx'] = 0
    df_out['lnxbar'] = 0 
    return df_out

def compileData(fileName):
    strlist = []
    rlist = []
    abslist = []
    templist = []
    hnd=open(fileName,'r')
    while True:
        buff=hnd.readline()
        if not buff:
            break
        try:
            strlist.append(buff.strip())
        except ValueError:
            pass
    hnd.close()
    for i in range(0,3):
        strlist.pop(0)
    for i in range(0,len(strlist)):
        templist = strlist[i].split()
        rlist.append(float(templist[0])-6)
        abslist.append(float(templist[1]))
    ziplist = list(zip(rlist,abslist))
    df_out = pd.DataFrame(ziplist,columns=['r','Abs'])
    return df_out

def cleanData(df_in, reflist):
    df_in['Abs'] = ss.savgol_filter(df_in['Abs'],7,3) ##ss filtering
    
    df_cut = df_in.query('Abs < 0') #think about making this just this line
    indexMax = df_cut.index[-1] + 3
    reflist.append(df_in.at[indexMax-2,'r'])
    df_out = df_in.truncate(before=indexMax)
    return df_out

def dervData(dfI):
    midpoint = (dfI['r']+dfI['r'].shift(1))/2 #Thanks PJ for this line it takes the data frame x1 +x2 then divides by two plots into a series
    dAbs = dfI['Abs'].diff()/dfI['r'].diff()
    inlist = list(zip(midpoint, dAbs))
    dfdx = pd.DataFrame(inlist, columns=['r','d1'])
    dfdx = dfdx.dropna()
    return dfdx

def findxBar(df_in):
    max_d1_index = df_in['d1'].idxmax()#find dA/dr maximum
    xbar = df_in.at[max_d1_index,'r']#declares xbar to be the
    return xbar

def grpPlot(df_Vel):
    plt.plot(df_Vel['r'],df_Vel['Abs'])
    plt.xlabel('radius(cm)',fontsize=15)
    plt.ylabel(r'$A_{280}$',fontsize=15)
    plt.title('Velocity Group Plot',fontsize=20)
    plt.savefig('GroupPlot.pdf')
    return True

def pltVelPlot(df_Vel, dfdx, time, velPlot=True, dervPlot=True):
    if velPlot == True and dervPlot == False:
        plt.figure(figsize=(5,4))
        plt.plot(df_Vel['r'],df_Vel['Abs'],'g-')
        plt.xlabel('radius(cm)',fontsize=15)
        plt.ylabel(r'$A_{280}$',fontsize=15)
        plt.title('Velocity Plot at {} min'.format(int(time/60)),fontsize=20)
    
    if velPlot == False and dervPlot == True:
        plt.figure(figsize=(5,4))
        plt.plot(dfdx['r'],dfdx['d1'],'r-')
        plt.xlabel('radius(cm)',fontsize=15)
        plt.ylabel(r'$A_{280}$',fontsize=15)
        plt.title('Velocity Plot at {} min'.format(int(time/60)),fontsize=20)

    if velPlot == True and dervPlot == True:
        fig, ax1 = plt.subplots()
        color = 'tab:red'
        ax1.set_xlabel('radius(cm)', fontsize=15)
        ax1.set_ylabel(r'$\frac{dA_{280}}{dr}$',color=color,fontsize=15)
        ax1.plot(dfdx['r'],dfdx['d1'],color=color)
        ax1.tick_params(axis='both',direction='in')
        
        ax2 = ax1.twinx()
        color = 'tab:green'
        ax2.set_ylabel(r'$A_{280}$',fontsize=15,color=color)
        ax2.plot(df_Vel['r'],df_Vel['Abs'], color=color,linewidth=3)
        ax2.set_title('Velocity and Derivative Plot {} min'.format(int(time/60)),fontsize=20)
        plt.show()

def GrubsTest(df_ref,strY):
    """ so this outliar test essentially works by stating
    if a point is x standard deviations(as according to the dictionary)
     if above x(GCrit) then value is outliar and to find Gvalue
     
    """
    while True:
        GCrit = gCrit[len(df_ref)]
        indexlist = []
        ylist = []
        Glist = []
        seriesX = df_ref['time']
        seriesY = df_ref[strY]
        res1 = stats.linregress(seriesX, seriesY)
        ylist = [res1.slope*i + res1.intercept for i in seriesX]
        ycalc = pd.Series(ylist)
        seriesDelta = abs(ycalc - seriesY)
        mu = seriesDelta.mean()
        sigma = seriesDelta.std()
        Glist = [(abs(i - mu)/sigma) for i in seriesDelta]
        df_ref['GValue'] = Glist
        df_query = df_ref.query('GValue > @GCrit')
        if len(df_query) == 0:
            break
        indexlist = df_query.index.tolist()
        for i in indexlist:
            df_ref = df_ref.drop(index = i)
        df_ref = df_ref.reset_index(drop=True)
        df_ref = df_ref.drop(columns=['GValue'])
    return df_ref

def getlogSlope(df_ref,plotregress=False):#using equation 16
    res = stats.linregress(df_ref['time'],df_ref['lnxbar'])
    if plotregress == True:
        linefit = [i*res.slope + res.intercept for i in df_ref['time']]
        plt.plot(df_ref['time'],df_ref['lnxbar'],'bo')
        plt.plot(df_ref['time'],linefit,'m-')
        plt.xlabel('time(s)',fontsize=15)
        plt.ylabel(r'$\ln{(\bar{x})}$',fontsize=15)
        plt.title('Natural log of the boundary vs time',fontsize=20)
        plt.tight_layout()
        plt.savefig('dlnxdt.pdf')
        plt.show()
    return res.slope

def getdxSlope(df_ref,plotregress=False):#using equation 2
    res = stats.linregress(df_ref['time'],df_ref['xbar'])
    if plotregress == True:
        linefit = [i*res.slope + res.intercept for i in df_ref['time']]
        plt.plot(df_ref['time'],df_ref['xbar'],'ro')
        plt.plot(df_ref['time'],linefit,linestyle='solid',color='tab:brown')
        plt.xlabel('time(s)',fontsize=15)
        plt.ylabel(r'$\bar{x}$',fontsize=15)
        plt.title('Boundary position as a function of time',fontsize=20)
        plt.savefig('dxdt.pdf')
        plt.show()
    return res.slope

def eqnsixteen(dln, w):
    w2 = w**2
    s = dln/w2
    S = s * 10**13
    return s,S

def eqntwentysix():
    pass

def eqntwo(df_ref, slope, plotregress=False):
    #s = dxbar/dt/xbar/w2
    w = df_ref.at[0,'w']
    w2 = w**2
    df_ref['s_2'] = slope / w2 / df_ref['xbar'] * (10**13)
    if plotregress == True:
        res = stats.linregress(df_ref['time'],df_ref['s_2'])
        linefit = [i*res.slope + res.intercept for i in df_ref['time']]
        plt.plot(df_ref['time'],df_ref['s_2'],'go')
        plt.plot(df_ref['time'],linefit,linestyle='solid',color='tab:purple')
        plt.xlabel('time(s)',fontsize=15)
        plt.ylabel('s (sec)',fontsize=15)
        plt.title('s plot as a function of time',fontsize=20)
        plt.savefig('splottime.pdf')
        plt.show()
        
    return df_ref, res.intercept

def MW(s):
    eta = 0.01002#poise units water at 20
    hydration = 0.5
    rho = .99823 #g/mL water at 20
    s = s *(10**-13) #s
    vbar = 0.7301 #mL/g
    perin = 1
    
    
    term1 = 6 *sc.pi * sc.N_A * s * eta
    term2 = (3*vbar / 4/ sc.pi / sc.N_A)**(1/3)
    term3 = (1 + (hydration/(vbar*rho)))**(1/3)
    top = term1 * term2 *term3 * perin
    bottom = (1 - (vbar*rho))
    solution = (top / bottom)**(3/2)
    return solution

def main():
    #1
    df_master = compileName('CELL2ABS.LOG')
    xbarlist = []
    xmenlist = []
    bool1 = False
    #2
    for i in range(0,4): #30 makes all plots linear
        for j in range(0,10):
            if i == 0 and j == 0:
                pass
            else:
                fileName = '000{}{}.RA2'.format(i,j)
                indexno = int('{}{}'.format(i,j))
                dfVel = compileData(fileName)
                dfVel = cleanData(dfVel, xmenlist)
                dfd1 = dervData(dfVel)
                xbarlist.append(findxBar(dfd1))
                #pltVelPlot(dfVel, dfd1,df_master.loc[indexno,'time']) #individual plots
                #bool1 = grpPlot(dfVel) #Group plots does not work with pltVelPlot
    if bool1 == True:
        #save fig
        plt.show()
    #3
    xmenseries = pd.Series(xmenlist)
    xbarseries = pd.Series(xbarlist)
    df_master['xmen'] = xmenseries
    df_master['xbar'] = xbarseries
    df_master = df_master.dropna()
    df_master['dx'] = df_master['xbar']-df_master['xmen']
    df_master['lnxbar'] = log(df_master['xbar'])
    df_master = GrubsTest(df_master, 'xbar')
    dlnxdt = getlogSlope(df_master, plotregress=True)
    dxbardt = getdxSlope(df_master, plotregress=True)
    s,S = eqnsixteen(dlnxdt,w = df_master.at[0,'w'])
    df_master,s_0 = eqntwo(df_master, dxbardt,plotregress=True)
    MW1 = MW(S)
    MW2 = MW(s_0)
    print('eqn 16: S={}Sv, MW={}MDa,\n eqn2: S={}Sv, MW={}MDa'.format(S,MW1/10**6,s_0,MW2/10**6))
    

#M^(2/3)*(1-vrho)/(N6pi eta) (3v/4piN)^(1/3)/(1+deltaH/vrho)^(1/3) /(f/f0)
#vbar calc from amino acid sequence
#perrin factor french paper
#delta H hydration draft we never got published



if __name__ == "__main__":
    main()

