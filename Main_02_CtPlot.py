import os
import glob
import pandas as pd
# import welib
import fastlib
import weio
import numpy as np
try:
    from pybra import clean_exceptions
except:
    pass
print(weio.__file__)


def spanwiseAD(tsAvg,vr_bar,rho,R,nB,suffix=''):
    nr=len(vr_bar)
    Columns     = [('r/R_[-]', vr_bar)]
    Columns.append(fastlib.extractSpanTS(tsAvg,nr,'B1N{:d}Alpha_[deg]','B1Alpha_[deg]'))
    Columns.append(fastlib.extractSpanTS(tsAvg,nr,'B1N{:d}AxInd_[-]'  ,'B1AxInd_[-]'  ))
    Columns.append(fastlib.extractSpanTS(tsAvg,nr,'B1N{:d}Cl_[-]'     ,'B1AxCl_[-]'   ))
    Columns.append(fastlib.extractSpanTS(tsAvg,nr,'B1N{:d}Cd_[-]'     ,'B1AxCd_[-]'   ))
    Columns.append(fastlib.extractSpanTS(tsAvg,nr,'B1N{:d}Vrel_[m/s]' ,'B1Vrel_[m/s]' ))
    Columns.append(fastlib.extractSpanTS(tsAvg,nr,'B1N{:d}Fx_[N/m]'   ,'B1Fx_[N/m]'   ))
    if Columns[-1][1] is not None:
        r=vr_bar*R
        Fx =Columns[-1][1]
        U0=tsAvg['Wind1VelX_[m/s]']
        print(U0)
        Ct=nB*Fx/(0.5 * rho * 2 * U0**2 * np.pi * r)
        Ct[vr_bar<0.01] = 0
        Columns.append(('B1Ct_[-]', Ct))
    Columns.append(fastlib.extractSpanTS(tsAvg,nr,'B1N{:d}Fy_[N/m]'   ,'B1Fy_[N/m]'   ))



    data     = np.column_stack([c for _,c in Columns if c is not None])
    ColNames = [n for n,_ in Columns if n is not None]
    if len(ColNames)<=0:
        print('No spanwise aero data from sim {}'.s)
    else:
        dfRad = pd.DataFrame(data= data, columns = ColNames)
        dfRad.to_csv(os.path.join(Outdir,'Parametric_Spanwise_Aero'+suffix+'.csv'),sep='\t',index=False)



def spanwisePostPro(FST_In,suffix=None):
    # --- Init
    fst = weio.FASTInputDeck(FST_In)
    ED  = fst.ED
    AD  = fst.Aero
    R =  ED ['TipRad']


    r_FST_struct = fastlib.ED_BldGag(ED)
    r_FST_aero = fastlib.AD_BldGag(AD,AD.Bld1) + ED['HubRad']

    df    = weio.read(FST_In.replace('.fst','.outb')).toDataFrame()
    dfAvg = fastlib.averageDF(df,avgMethod = 'constantwindow',avgParam=5)
    Ct=dfAvg['RtAeroCt_[-]'].values


    if suffix is None:
        suffix='Ct{:.2f}'.format(Ct[0])

    spanwiseAD(dfAvg.iloc[0], r_FST_aero/R, AD['AirDens'], R, 3, suffix)




# --- Main Parameters
# SimDir='OpenFAST_Parametric/'
# Outdir = './'
# 
# fstfiles=glob.glob(SimDir + '*.fst')
# print(fstfiles)
# for FST_In in fstfiles:
#     spanwisePostPro(FST_In)



# --- Main Parameters
FST_In = 'OpenFAST_AD/Main_Onshore_OF2.fst'
Outdir = './'
spanwisePostPro(FST_In,'PrescCt')

