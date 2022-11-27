#!/home/xing/Downloads/software/anaconda3/bin/python3
import time
import numpy as np
from phonopy import load
from phonopy.spectrum.dynamic_structure_factor import atomic_form_factor_WK1995
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import colors
from matplotlib.colors import ListedColormap
import scipy
from scipy import ndimage
from matplotlib import ticker
from constants import THz2meV
import scipy.linalg

plt.rcParams.update({'font.size':24})
start = time.time()

#### START OF USER INPUTS #####
primitive_cell = [[1,0,0],[0,1,0],[0,0,1]]
supercell = [2,2,2]
### Use either FORCE_SETS or FORCE_CONSTANTS ### 
forcesets_file = None
forceconstants_file = 'FORCE_CONSTANTS' 
## Q and E inputs ##
Q_start = np.array([1,3,0]) # Q bands are in respect of primitive cell if you give primitive_cell
Q_end = np.array([3,1,0])
Q_steps = 40
Temperature = 300
# Neutron coherent scattering length can be found at https://www.ncnr.nist.gov/resources/n-lengths/
coh_scatter_length ={'Na': 3.63, 'Cl': 9.5770}
E_min = 0
E_max = 25
E_step = 0.02
xlabels = '[2+H, 2-H, 0] (r.l.u.)'
savefigas = 'temp.png'
colormax = 2
colormin = 0
e_resolution = 50 ## unit of E_step
num_ybins = 6 ## number of yticks 
## save txt ##
savetxt_as = 'temp_sqe.txt'

#### END OF USER INPUTS ####

## The script below are mainly copied from Phonopy website
## But I changed something to make it output Q and SQE for plotting
def get_AFF_func(f_params):
    def func(symbol, Q):
        return atomic_form_factor_WK1995(Q, f_params[symbol])
    return func

def run(phonon,
        Qpoints,
        temperature,
        atomic_form_factor_func=None,
        scattering_lengths=None):
    # Transformation to the Q-points in reciprocal primitive basis vectors
    Q_prim = np.dot(Qpoints, phonon.primitive_matrix)
    # Q_prim must be passed to the phonopy dynamical structure factor code.
    phonon.run_dynamic_structure_factor(
        Q_prim,
        temperature,
        atomic_form_factor_func=atomic_form_factor_func,
        scattering_lengths=scattering_lengths,
        freq_min=1e-3)
    dsf = phonon.dynamic_structure_factor
    q_cartesian = np.dot(dsf.qpoints,
                         np.linalg.inv(phonon.primitive.get_cell()).T)
    distances = np.sqrt((q_cartesian ** 2).sum(axis=1))

    #Comment out original prints. They write what are what
    #print("# [1] Distance from Gamma point,")
    #print("# [2-4] Q-points in cubic reciprocal space, ")
    #print("# [5-8] 4 band frequencies in meV (becaues of degeneracy), ")
    #print("# [9-12] 4 dynamic structure factors.")
    #print("# For degenerate bands, dynamic structure factors are summed.")
    #print("")

    Qpoints = np.array(Qpoints)
    #print (dsf.frequencies.shape)
    #print (dsf.dynamic_structure_factors.shape)
    # Use as iterator
    #for i in range(len(dsf.frequencies.shape[0])):
    #    SandE = []
    #    for j in range(len(dsf.frequencies.shape[1])):
    #        SandE.append([dsf.frequencies[i,j],Qpoints[i],dsf.dynamic_structure_factors[i,j]])

    SandE = np.array([dsf.frequencies,dsf.dynamic_structure_factors])
    return SandE


if __name__ == '__main__':
    phonon = load(supercell_matrix=supercell,
                  primitive_matrix=primitive_cell,
                  unitcell_filename="POSCAR",
                  force_sets_filename=forcesets_file,
                  force_constants_filename=forceconstants_file
                  )

    q_start = Q_start
    q_end = Q_end
    band = []
    Qpoints = []

    for i in range(Q_steps+1):
        band.append(list(q_start+(q_end-q_start)/Q_steps*i))
    Qpoints = band
    ### To avoid DW singularity at Q=[0,0,0]
    Qpoints_abs = scipy.linalg.norm(Qpoints,axis=1)
    zero_index = np.where(Qpoints_abs==0)[0]
    if zero_index.shape[0] !=0:
        Qpoints[zero_index[0]] = np.array([1e-6,1e-6,1e-6])
    #print (Qpoints)

    # Mesh sampling phonon calculation is needed for Debye-Waller factor.
    # This must be done with is_mesh_symmetry=False and with_eigenvectors=True.
    mesh = [11, 11, 11]
    phonon.run_mesh(mesh,
                    is_mesh_symmetry=False, # symmetry must be off
                    with_eigenvectors=True) # eigenvectors must be true
    temperature = Temperature

    IXS = False

    if IXS:
        # For IXS, atomic form factor is needed and given as a function as
        # a parameter.
        # D. Waasmaier and A. Kirfel, Acta Cryst. A51, 416 (1995)
        # f(Q) = \sum_i a_i \exp((-b_i Q^2) + c
        # Q is in angstron^-1
        # a1, b1, a2, b2, a3, b3, a4, b4, a5, b5, c
        f_params = {'Na': [3.148690, 2.594987, 4.073989, 6.046925,
                           0.767888, 0.070139, 0.995612, 14.1226457,
                           0.968249, 0.217037, 0.045300],  # 1+
                    'Cl': [1.061802, 0.144727, 7.139886, 1.171795,
                           6.524271, 19.467656, 2.355626, 60.320301,
                           35.829404, 0.000436, -34.916604]}  # 1-
        AFF_func = get_AFF_func(f_params)
        run(phonon,
            Qpoints,
            temperature,
            atomic_form_factor_func=AFF_func)
    else:
        # For INS, scattering length has to be given.
        # The following values is obtained at (Coh b)
        # https://www.nist.gov/ncnr/neutron-scattering-lengths-list
        output = run(phonon,
                     Qpoints,
                     temperature,
                     scattering_lengths=coh_scatter_length) # At this point, we need to input coh scattering length manually. But it's good for isotopes.
        
    ## output has shape as (2,len(Qpoints),branches), the [0,:,:] is for frequency and the [1,:,:] is for SQE; The frequency is in THz unit
    for i in range(len(Qpoints)):
        output[0,i,:] *= THz2meV

    ## Resolution for energy and Q ##
    def Gauss(evalues,qvalues,ecenter,qcenter,sigmae,sigmaq):
        Ne = 1.0/np.sqrt(2*np.pi*sigmae**2)
        Nq = 1.0/np.sqrt(2*np.pi*sigmaq**2)
        delta_q = qvalues-qcenter
        q_square = np.dot(delta_q,delta_q)
        print (evalues.shape,ecenter.shape,q_square.shape)
        Res = Ne*np.exp(-(evalues-ecenter)**2/(2*sigmae**2))*Nq*np.exp(-(q_square)/(2*sigmaq**2))
        return Res

    ## Bin the SQE ##
    MinimumEnergy = E_min
    MaximumEnergy = E_max
    EnergyStep = E_step
    deltae = EnergyStep 
    ne = int(1 + (MaximumEnergy - MinimumEnergy)/EnergyStep)
    Evec = [MinimumEnergy,EnergyStep,ne, MaximumEnergy]
    evalues = np.arange(E_min,ne*deltae+deltae,deltae)    

    nql = Q_steps+1
    BinnedSQE=np.zeros((int(nql),int(ne)))

    ntotal = 0

    # Energy binning
    #for ih in range(nql): # qpoints
    #    for j in range(len(output[0,0,:])): # branche
    #        EIndex=int(round(output[0][ih][j]/deltae))
    #        if EIndex < ne :
    #            BinnedSQE[ih,EIndex] += output[1][ih][j]
                #ntotal=ntotal+1
                #if (ntotal%1000==0):#print status only every 1000 points
                #    print (ntotal,' of ',nqtot,' ',ntotal*100.0/nqtot,' percent')
    
    # Energy binning with convolution
    for ih in range(nql): # qpoints
        for j in range(len(output[0,0,:])): # branch
            EIndex=int(round(output[0][ih][j]/deltae))
            if EIndex < ne :
                BinnedSQE[ih,EIndex] += output[1][ih][j]
        BinnedSQE[ih,:] = ndimage.gaussian_filter(BinnedSQE[ih,:],e_resolution)

    mapName = 'viridis'
    zeromask = np.where(BinnedSQE == 0.0)
    BinnedSQE[zeromask] = 1e-6
    np.savetxt(savetxt_as,BinnedSQE.T)

    plotVals = np.log10(BinnedSQE.T)
    #####
    #plotVals = pl.log10(BinnedSQE.T)
    norm=colors.Normalize(np.min(plotVals),np.max(plotVals))
    figure, ax = plt.subplots(figsize=(10,8))

    im = ax.imshow(plotVals,
               #extent=(xmin,xmax,ymin,ymax),
               aspect='auto',
               interpolation='bicubic',
               #norm=norm,
               vmin = colormin,
               vmax = colormax,
               cmap=mapName)

    xpos = np.linspace(0,nql-1,6)
    q_range = q_end-q_start
    index = np.nonzero(q_range)
    x_max = q_range[index[0][0]]
    #xticks = np.linspace(0,x_max,6)
    #fmt = lambda x: "{:.2f}".format(x) # the function to only keep two digits
    #ax.set_xticks(xpos)
    #ax.set_xticklabels([fmt(i) for i in xticks])
    ax.set_xticks([0,nql-1])
    ax.set_xticklabels([q_start,q_end])
 
    ybin_num = num_ybins
    ypos = np.linspace(0,ne,ybin_num)
    yticks = np.linspace(MinimumEnergy,MaximumEnergy,ybin_num)
    fmt = lambda x: "{:.1f}".format(x) # the function to only keep two digits
    ax.set_yticks(ypos)
    ax.set_yticklabels([fmt(i) for i in yticks])

    ax.set_xlim(0,nql-1)
    ax.set_ylim(0,ne)

    ax.set_xlabel(xlabels)
    ax.set_ylabel('Energy (meV)')

    cbaxes = figure.add_axes([0.90, 0.23, 0.03, 0.6])
    cb2 = figure.colorbar(im,cax=cbaxes,orientation = 'vertical')
    tick_locator = ticker.MaxNLocator(nbins=int((colormax-colormin)//1))
    cb2.locator = tick_locator
    cb2.update_ticks()

    plt.subplots_adjust(left=0.15,right = 0.85,bottom = 0.15, top = 0.9, wspace = 0.15,hspace=0.1)
    plt.savefig(savefigas)

    end3 = time.time()
    print ('Total time: %f seconds'%(end3 - start))
