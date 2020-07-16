from __future__ import print_function
import numpy as np
import os
import weio
import fastlib
import fileinput
try:
    from pybra.clean_exceptions import *
except:
    pass



def ParametricCT(CT,CTname,BEM=False):
    # --- Parameters for this script
    ref_dir          = 'OpenFAST_AD/'   # Folder where the fast input files are located (will be copied)
    if BEM:
        work_dir         = 'Parametric_Ct_BEM/'     # Output folder (will be created)
    else:
        work_dir         = 'Parametric_Ct_CFD/'     # Output folder (will be created)
    main_file        = 'Main_Onshore_OF2.fst'  # Main file in ref_dir, used as a template
    FAST_EXE         = 'openfast_x64s_AD.exe' # Location of a FAST exe (and dll)
    print('>>> Parametric CT' + work_dir)
    # --- Defining the parametric study  (list of dictionnaries with keys as FAST parameters)
    Tmax   = 1000
    BaseDict = {'FAST|TMax': Tmax, 'FAST|DT': 0.01, 'FAST|DT_Out': 0.1}
    PARAMS=[]
    for Ct,Ctname in zip(CT,CTname):
        p=BaseDict.copy()
        p['__name__']       = 'CT_{:03.2f}'.format(Ctname)
        if BEM:
            p['FAST|CompInflow']       = 1
            p['AeroFile|WakeMod']      = 0
        else:
            p['FAST|CompInflow']       = 2
            p['AeroFile|WakeMod']      = 0
        p['AeroFile|PrescribedAD'] = True
        p['AeroFile|PrescribedCt'] = Ct
        PARAMS.append(p)
    # --- Generating all files in a workdir
    fastfiles=fastlib.templateReplace(PARAMS,ref_dir,workdir=work_dir,RemoveRefSubFiles=True,main_file=main_file, oneSimPerDir=True)

    # --- Creating a batch script just in case
    #fastlib.writeBatch(os.path.join(work_dir,'_RUN_ALL.bat'),fastfiles,fastExe=FAST_EXE)
    # --- Running the simulations
    if BEM:
        fastlib.run_fastfiles(fastfiles,fastExe=FAST_EXE,parallel=True,ShowOutputs=False,nCores=2)

    if not BEM:
        # --- Replacing in Nalu input file
        for f in fastfiles:
            parent   = os.path.dirname(f)
            basename = os.path.basename(f)
            nalu_file = os.path.join(parent, 'alm_simulation.yaml')
            print(nalu_file)

            # Python 2/3 compatible, otherwise use "with"
            file = fileinput.FileInput(nalu_file, inplace=True, backup='.bak')
            for line in file:
                print(line.replace('XXX.fst' , basename), end='')
            file.close()


def ParametricPitch(Pitch,CT,BEM):
    # --- Parameters for this script
    ref_dir          = 'OpenFAST/'   # Folder where the fast input files are located (will be copied)
    work_dir         = 'Parametric_Pitch_CFD/'     # Output folder (will be created)
    main_file        = 'Main_Onshore_OF2.fst'  # Main file in ref_dir, used as a template
    FAST_EXE         = 'OpenFAST2_x64s_ebra.exe' # Location of a FAST exe (and dll)
    Tmax             = 1000
    if BEM:
        Tmax   = 10
        work_dir         = '_BEM/Parametric_Pitch_BEM/'     # Output folder (will be created)

    # --- Defining the parametric study  (list of dictionnaries with keys as FAST parameters)
    print('>>> Parametric Pitch' + work_dir)

    BaseDict = {'FAST|TMax': Tmax, 'FAST|DT': 0.01, 'FAST|DT_Out': 0.1}
    PARAMS=[]
    for pitch,Ct in zip(Pitch,CT):
        p=BaseDict.copy()
        #p['__name__']       = 'Pitch_{:04.1f}'.format(pitch)
        print('Ct',Ct, 'Pitch',pitch)
        p['__name__']       = 'CT_{:03.2f}'.format(Ct)
        p['EDFile|BlPitch(1)']     = pitch
        p['EDFile|BlPitch(2)']     = pitch
        p['EDFile|BlPitch(3)']     = pitch
        if not BEM:
            p['FAST|CompInflow']       = 2
            p['AeroFile|WakeMod']      = 0
        PARAMS.append(p)

    oneSimPerDir=not BEM # CFD wants one sim per dir

    # --- Generating all files in a workdir
    fastfiles=fastlib.templateReplace(PARAMS,ref_dir,workdir=work_dir,RemoveRefSubFiles=True,main_file=main_file, oneSimPerDir=oneSimPerDir)

    # --- Creating a batch script just in case
    #fastlib.writeBatch(os.path.join(work_dir,'_RUN_ALL.bat'),fastfiles,fastExe=FAST_EXE)
    if BEM:
        # --- Running the simulations
        fastlib.run_fastfiles(fastfiles,fastExe=FAST_EXE,parallel=True,ShowOutputs=False,nCores=2)

        # --- Simple Postprocessing
        outFiles = [os.path.splitext(f)[0]+'.outb' for f in fastfiles]
        #avg_results = fastlib.averagePostPro(outFiles,avgMethod='periods',avgParam=1, ColMap = {'WS_[m/s]':'Wind1VelX_[m/s]'},ColSort='WS_[m/s]')
        avg_results = fastlib.averagePostPro(outFiles,avgMethod='constantwindow',avgParam=Tmax/2, ColMap = {'WS_[m/s]':'Wind1VelX_[m/s]'},ColSort='BldPitch1_[deg]')
        avg_results.to_csv(os.path.join(work_dir,'ParametricPitch.csv'),sep='\t',index=False)
        return avg_results
    else:
        # --- Replacing in Nalu input file
        for f in fastfiles:
            parent   = os.path.dirname(f)
            basename = os.path.basename(f)
            nalu_file = os.path.join(parent, 'alm_simulation.yaml')
            print(nalu_file)

            file= fileinput.FileInput(nalu_file, inplace=True, backup='.bak')
            for line in file:
                print(line.replace('XXX.fst' , basename), end='')
            file.close()

def ParametricOmega(Omega,CT,BEM):
    # --- Parameters for this script
    ref_dir          = 'OpenFAST/'   # Folder where the fast input files are located (will be copied)
    work_dir         = 'Parametric_RPM_CFD/'     # Output folder (will be created)
    main_file        = 'Main_Onshore_OF2.fst'  # Main file in ref_dir, used as a template
    FAST_EXE         = 'OpenFAST2_x64s_ebra.exe' # Location of a FAST exe (and dll)
    Tmax             = 1000
    if BEM:
        Tmax   = 10
        work_dir         = '_BEM/Parametric_RPM_BEM/'     # Output folder (will be created)

    # --- Defining the parametric study  (list of dictionnaries with keys as FAST parameters)
    print('>>> Parametric RPM' + work_dir)

    BaseDict = {'FAST|TMax': Tmax, 'FAST|DT': 0.01, 'FAST|DT_Out': 0.1}
    PARAMS=[]
    for rpm,Ct in zip(RPM,CT):
        p=BaseDict.copy()
        #p['__name__']       = 'Pitch_{:04.1f}'.format(pitch)
        print('Ct',Ct, 'RPM',rpm)
        p['__name__']       = 'CT_{:03.2f}'.format(Ct)
        p['EDFile|RotSpeed'] = rpm
        if not BEM:
            p['FAST|CompInflow']       = 2
            p['AeroFile|WakeMod']      = 0
        PARAMS.append(p)

    oneSimPerDir=not BEM # CFD wants one sim per dir

    # --- Generating all files in a workdir
    fastfiles=fastlib.templateReplace(PARAMS,ref_dir,workdir=work_dir,RemoveRefSubFiles=True,main_file=main_file, oneSimPerDir=oneSimPerDir)

    # --- Creating a batch script just in case
    #fastlib.writeBatch(os.path.join(work_dir,'_RUN_ALL.bat'),fastfiles,fastExe=FAST_EXE)
    if BEM:
        # --- Running the simulations
        fastlib.run_fastfiles(fastfiles,fastExe=FAST_EXE,parallel=True,ShowOutputs=False,nCores=2)

        # --- Simple Postprocessing
        outFiles = [os.path.splitext(f)[0]+'.outb' for f in fastfiles]
        #avg_results = fastlib.averagePostPro(outFiles,avgMethod='periods',avgParam=1, ColMap = {'WS_[m/s]':'Wind1VelX_[m/s]'},ColSort='WS_[m/s]')
        avg_results = fastlib.averagePostPro(outFiles,avgMethod='constantwindow',avgParam=Tmax/2, ColMap = {'WS_[m/s]':'Wind1VelX_[m/s]'},ColSort='RotSpeed_[rpm]')
        avg_results.to_csv(os.path.join(work_dir,'ParametricRPM.csv'),sep='\t',index=False)
        return avg_results
    else:
        # --- Replacing in Nalu input file
        for f in fastfiles:
            parent   = os.path.dirname(f)
            basename = os.path.basename(f)
            nalu_file = os.path.join(parent, 'alm_simulation.yaml')
            print(nalu_file)

            with fileinput.FileInput(nalu_file, inplace=True, backup='.bak') as file:
                for line in file:
                    print(line.replace('XXX.fst' , basename), end='')



if __name__=='__main__':
    CT   = np.arange(0.1,2,0.2)

    # --- Parametric rpm
    df    = weio.read('OMEGA_CT.csv').toDataFrame()
    df = df.sort_values('RtAeroCt_[-]')
    ct    = df['RtAeroCt_[-]'].values
    rpm = df['RotSpeed_[rpm]'].values

    RPM = np.interp(CT,ct,rpm)
    ParametricOmega(RPM,CT,BEM=False)

    
    # --- Parametric pitch
    df    = weio.read('ParametricStudyPitch.csv').toDataFrame()
    df = df.sort_values('RtAeroCt_[-]')
    ct    = df['RtAeroCt_[-]'].values
    pitch = df['BldPitch1_[deg]'].values
    PITCH = np.interp(CT,ct,pitch)

    ParametricPitch(PITCH,CT,BEM=False)
    #   ParametricPitch(PITCH,CT,BEM=True)

    # --- Parametric CT
    CT_for_AD=CT*(1+0.6**2/3)
    ParametricCT(CT_for_AD, CT, BEM=False)
    #ParametricCT(CT_for_AD, CT, BEM=True)


