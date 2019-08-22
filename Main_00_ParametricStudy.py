import numpy as np
import os
import fastlib

def ParametricExample():
    """ Example to run a set of OpenFAST simulations (parametric study)

    This script uses a reference directory (`ref_dir`) which contains a reference input file (.fst)
    1) The reference directory is copied to a working directory (`work_dir`).
    2) All the fast input files are generated in this directory based on a list of dictionaries (`PARAMS`).
    For each dictionary in this list:
       - The keys are "path" to a input parameter, e.g. `EDFile|RotSpeed`  or `FAST|TMax`.
         These should correspond to the variables used in the FAST inputs files.
       - The values are the values corresponding to this parameter
    For instance:
         PARAMS[0]['EDFile|RotSpeed']       = 5
         PARAMS[0]['InflowFile|HWindSpeed'] = 10

    3) The simulations are run, successively distributed on `nCores` CPUs.
    4) The output files are read, and averaged based on a method (e.g. average over a set of periods,
        see averagePostPro in fastlib for the different averaging methods).
       A pandas DataFrame is returned

    """
    # --- Parameters for this script
    ref_dir          = 'OpenFAST/'   # Folder where the fast input files are located (will be copied)
    work_dir         = '_ParametricStudy/'     # Output folder (will be created)
    main_file        = 'Main_Onshore_OF2.fst'  # Main file in ref_dir, used as a template
    FAST_EXE         = 'OpenFAST2_x64s_ebra.exe' # Location of a FAST exe (and dll)

    # --- Defining the parametric study  (list of dictionnaries with keys as FAST parameters)
    #Lambda = 15
    #R      = 63.5
    #U0     = 10
    #RPM    =  Lambda*U0/R * 60 / (2*np.pi)
    #print(RPM)
    Tmax   = 10
    PITCH = np.linspace(-7,6.4,30)

    BaseDict = {'FAST|TMax': Tmax, 'FAST|DT': 0.01, 'FAST|DT_Out': 0.1}
    PARAMS=[]
    for pitch in PITCH:
        p=BaseDict.copy()
        p['__name__']       = 'Pitch_{:04.1f}'.format(pitch)
        p['EDFile|BlPitch(1)']     = pitch
        p['EDFile|BlPitch(2)']     = pitch
        p['EDFile|BlPitch(3)']     = pitch
        PARAMS.append(p)
    # --- Generating all files in a workdir
    fastfiles=fastlib.templateReplace(ref_dir,PARAMS,workdir=work_dir,RemoveRefSubFiles=True,main_file=main_file)
    print(fastfiles)

    # --- Creating a batch script just in case
    fastlib.writeBatch(os.path.join(work_dir,'_RUN_ALL.bat'),fastfiles,fastExe=FAST_EXE)
    # --- Running the simulations
    fastlib.run_fastfiles(fastfiles,fastExe=FAST_EXE,parallel=True,ShowOutputs=False,nCores=2)

    # --- Simple Postprocessing
    outFiles = [os.path.splitext(f)[0]+'.outb' for f in fastfiles]
    #avg_results = fastlib.averagePostPro(outFiles,avgMethod='periods',avgParam=1, ColMap = {'WS_[m/s]':'Wind1VelX_[m/s]'},ColSort='WS_[m/s]')
    avg_results = fastlib.averagePostPro(outFiles,avgMethod='constantwindow',avgParam=Tmax/2, ColMap = {'WS_[m/s]':'Wind1VelX_[m/s]'},ColSort='BldPitch1_[deg]')
    avg_results.to_csv('ParametricStudyPitch.csv',sep='\t',index=False)
    print(avg_results['RtAeroCt_[-]'])

    return avg_results



if __name__=='__main__':
    ParametricExample()
