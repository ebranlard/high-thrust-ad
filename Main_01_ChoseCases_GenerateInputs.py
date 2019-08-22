import numpy as np
import os
import weio
import fastlib



def ParametricCT(CT):
    # --- Parameters for this script
    ref_dir          = 'OpenFAST_AD/'   # Folder where the fast input files are located (will be copied)
    work_dir         = 'Parametric_Ct_CFD/'     # Output folder (will be created)
    main_file        = 'Main_Onshore_OF2.fst'  # Main file in ref_dir, used as a template
    FAST_EXE         = 'OpenFAST2_x64s_ebra.exe' # Location of a FAST exe (and dll)
    # --- Defining the parametric study  (list of dictionnaries with keys as FAST parameters)
    Tmax   = 600
    BaseDict = {'FAST|TMax': Tmax, 'FAST|DT': 0.01, 'FAST|DT_Out': 0.1}
    PARAMS=[]
    for Ct in CT:
        p=BaseDict.copy()
        p['__name__']       = 'CT_{:03.2f}'.format(Ct)
        p['FAST|CompInflow']       = 2
        p['AeroFile|WakeMod']      = 0
        p['AeroFile|PrescribedAD'] = True
        p['AeroFile|PrescribedCt'] = Ct
        PARAMS.append(p)
    # --- Generating all files in a workdir
    fastfiles=fastlib.templateReplace(ref_dir,PARAMS,workdir=work_dir,RemoveRefSubFiles=True,main_file=main_file)

    # --- Creating a batch script just in case
    fastlib.writeBatch(os.path.join(work_dir,'_RUN_ALL.bat'),fastfiles,fastExe=FAST_EXE)
    # --- Running the simulations
    #fastlib.run_fastfiles(fastfiles,fastExe=FAST_EXE,parallel=True,ShowOutputs=False,nCores=2)

def ParametricPitch(Pitch,BEM):
    # --- Parameters for this script
    ref_dir          = 'OpenFAST/'   # Folder where the fast input files are located (will be copied)
    if BEM:
        work_dir         = 'Parametric_Pitch_BEM/'     # Output folder (will be created)
    else:
        work_dir         = 'Parametric_Pitch_CFD/'     # Output folder (will be created)
    main_file        = 'Main_Onshore_OF2.fst'  # Main file in ref_dir, used as a template
    FAST_EXE         = 'OpenFAST2_x64s_ebra.exe' # Location of a FAST exe (and dll)

    # --- Defining the parametric study  (list of dictionnaries with keys as FAST parameters)
    if BEM:
        Tmax   = 10
    else:
        Tmax   = 600

    BaseDict = {'FAST|TMax': Tmax, 'FAST|DT': 0.01, 'FAST|DT_Out': 0.1}
    PARAMS=[]
    for pitch in Pitch:
        p=BaseDict.copy()
        p['__name__']       = 'Pitch_{:04.1f}'.format(pitch)
        p['EDFile|BlPitch(1)']     = pitch
        p['EDFile|BlPitch(2)']     = pitch
        p['EDFile|BlPitch(3)']     = pitch
        PARAMS.append(p)
    # --- Generating all files in a workdir
    fastfiles=fastlib.templateReplace(ref_dir,PARAMS,workdir=work_dir,RemoveRefSubFiles=True,main_file=main_file)

    if BEM:
        # --- Creating a batch script just in case
        fastlib.writeBatch(os.path.join(work_dir,'_RUN_ALL.bat'),fastfiles,fastExe=FAST_EXE)
        # --- Running the simulations
        fastlib.run_fastfiles(fastfiles,fastExe=FAST_EXE,parallel=True,ShowOutputs=False,nCores=2)

        # --- Simple Postprocessing
        outFiles = [os.path.splitext(f)[0]+'.outb' for f in fastfiles]
        #avg_results = fastlib.averagePostPro(outFiles,avgMethod='periods',avgParam=1, ColMap = {'WS_[m/s]':'Wind1VelX_[m/s]'},ColSort='WS_[m/s]')
        avg_results = fastlib.averagePostPro(outFiles,avgMethod='constantwindow',avgParam=Tmax/2, ColMap = {'WS_[m/s]':'Wind1VelX_[m/s]'},ColSort='BldPitch1_[deg]')
        avg_results.to_csv(os.path.join(work_dir,'ParametricPitch.csv'),sep='\t',index=False)
        return avg_results

if __name__=='__main__':
    #Lambda = 15
    #R      = 63.5
    #U0     = 10
    #RPM    =  Lambda*U0/R * 60 / (2*np.pi)
    #print(RPM)


    df    = weio.read('ParametricStudyPitch.csv').toDataFrame()
    df = df.sort_values('RtAeroCt_[-]')
    ct    = df['RtAeroCt_[-]'].values
    pitch = df['BldPitch1_[deg]'].values

    CT=np.arange(0.1,2,0.2)
    PITCH = np.interp(CT,ct,pitch)

    #import matplotlib.pyplot as plt
    #plt.plot(pitch,ct)
    #plt.plot(PITCH,CT,'o')
    #plt.show()
    ParametricCT(CT)
    ParametricPitch(PITCH,BEM=False)

    #ParametricPitch(PITCH,BEM=True)


