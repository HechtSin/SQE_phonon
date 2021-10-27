
#Simphonies Input file

SystemName = 'STO' #String containing chemical formula


ScatteringType='neutron'  #neutron or x-ray


FirstPrinciplesCalculationType = 'VASP'   #Only VASP is supported at this point


#File names as they are given in the working directory
ForceConstants = 'FORCE_CONSTANTS'  #
StructureFile = 'POSCAR'

#parallelization options
ParallelCalculation = False
NumProcessors = 4

# Force Constant Optimization
Optimization = False

ResolutionConvolution= True   #True/False: convolve calculated SQE with estimated instrumental energy resolution
                   #implementation needed


InstrumentName='ARCS'  #Can be either ARCS, CNCS, HERIX, Generic; Used for estimated energy resolution convolution
                   #if more accurate resolution is desired, please specify values:
                   #FC, for fermi chopper frequency (INS only)
                   # Implementation needed here
                   #Generic is Gaussian Width
                   #Add parameter for Q resolution
                   
                   

ResolutionType='instrument' #Can be either polynomial (provided polynomial coefficients by InstrumentScientist
                            #'constant' - a constant resolution in energy transfer (used for IXS)
                            #'polynomial' - uses fourth order polynomial to calculate resolution
                            #'instrument' - calculates the estimated energy resolution using a windsor approximation and instrument parameters
                            
IncidentEnergy = 40  #Needed only when ResolutionType is set to Instrument. In this case energy max should not exceed incident energy, lest the values are meaningless.
ElasticFWHM = 1  #elastic line energy resolution, FWHM
QResolution = 0.2 #constant values only supported at this time

PolynomialCoef = [ 1.02188553396572,  
                  -0.0594951491635698,  
                  0.000791281838645456,  
                  2.58258569481579e-05,   
                  -1.64952077361938e-07  ]  #coeffiecient required for polynomial fitting of resolution function
                            
#if primitive cell is different from conventional cell
#user needs to be consistent with how the cells are defined.

PrimitiveVectors = [[1,0,0],[0,1,0],[0,0,1]]  #Vectors describing primitive cell
#PrimitiveVectors = [[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]]  #Vectors describing primitive cell
SuperCellVectors = [[2,0,0],[0,2,0],[0,0,2]] #vectors describing supercell modfication
#SuperCellVectors = [[0,1,1],[1,0,1],[1,1,0]] #vectors describing supercell modfication


#If your system is hexagonal, you will need to provide vectors to transform to the 
#rhombohedral system for phonopy.  Otherwise, leave RhombohedralModification False.

RhombohedralModification= False  #If calculation was performed with a rhombohedral unit cell, set this flag to TRUE and
                                #also specify Rhombohedral vectors

RhombohedralVectors = [[1.0/3.0,-1.0/3.0,-1.0/3.0],
                       [-2.0/3.0,-1.0/3.0,-1.0/3.0],
                       [1.0/3.0,2.0/3.0,-1.0/3.0]]  #A list of vectors the primitive rhombohedral cell
                                                    #requires RhombohedralModification to
                                                    #to be true This should likely be generalized so
                                                    #that one can make any arbitrary primitive axis translation.

Temperature = 300.0  #Kelvin

#Calculation type can be 'sqe' or 'chi'
CalculationType = 'sqe'

CalculateFullVolume = False

CalculateDirection = True
NumberOfDirections = 1  #Number of directions to calculate.  


 #DirectionSettings: for each desired calculation input a dictionary of values with the following fields.  Must be consistent with
#NumberOfDirections defined above.
#Each calculated direction needs:
#'coordinates' Coordinate system specifying the direction; can be for arbitrary direction but all vectors must be perpendicular
#'namestring' : name for the output files; convention is to use abbreviated system name and calculated direction
#'energyMin' : minimum energy for the calculation; best if something non-zero (ie, off Bragg peak)
#'energyMax' : maximum energy for the calculation.
#'energyBin' : Step size in energy;

#'origin': start value of changing index for Q slice
#'step': step size for changing index
#'final': final position for changing index
#'changeIndex': numerical index of point to change (indexed from 0)
#'indexShift': arbitrary offset from coordinate system for calculation
#'SliceIntegration' : True to integrate the slice along the other two directions; False otherwise
#'IntegrationWidth' : Amount to integrate in the perpendicular direction; gives start, stop, and step size for each perpendicular direction.

ReciprocalVolumeCalcSettings = [{'coordinates':[[1,0,0],[0,1,0],[0,0,1]],
                      'namestring':'Si_FullVolume',
                      'titleString' : 'Si_FullVolume ',
                      'xLabel' : 'H00 (rlu)',
                      'energyMin':0.5,
                      'energyMax':100.0,
                      'energyBin':0.5,
                      'QxMin': 0.0,
                      'QxStep':0.02,
                      'QxMax': 1.00,
                      'QyMin': 0.0,
                      'QyStep': 0.02,
                      'QyMax': 1.0,
                      'QzMin':0.0,
                      'QzStep': 0.01,
                      'QzMax': 0.10,
                      'changeIndex': 0,
                      'indexShift': [0,0,0],}
                                ]

#DirectionSettings = [{'coordinates': [[1,-1,1],[0,0,1],[1,-1,0]],
#DirectionSettings = [{'coordinates': [[1,1,1],[-1,-1,2],[1,-1,0]],
DirectionSettings = [{'coordinates': [[1,0,0],[0,1,0],[0,0,1]],
                      'namestring':'4K0',
                      'titleString' : 'blank',
                      'xLabel' : '4 K 0',
                      'energyMin':0,
                      'energyMax':25,
                      'energyBin':0.01,
                      'rluOrigin': 0,
                      'rluStop':1,
                      'rluStep':0.01,
                      'changeDirection': 1,
                      'rluOffset': [4,0,0],
                      'SliceIntegration':False,
                      'IntegrationWidth':[[-0.01,0.01,0.005],[-0.01,0.01,0.005]]
                      }]
# ,
#                       {'coordinates': [[1,1,0],[0,0,1],[-1,1,0]],
#                       'namestring':'Si_000_220',
#                       'titleString' : 'Si HH0',
#                       'xLabel' : 'HH0 (rlu)',
#                       'energyMin':0.5,
#                       'energyMax':100.0,
#                       'energyBin':0.5,
#                       'origin': 0.0,
#                       'step':0.01,
#                       'final': 2.0,
#                       'changeIndex': 1,
#                       'indexShift': [0,0,0],
#                       'SliceIntegration':False,
#                       'IntegrationWidth':[[-0.2,0.2,0.05],[-0.2,0.2,0.05]]}]



#Flags for making plots.  Plot Limits is color scaling limits for logarithmic intensity scale
Plots = True
dispersion = True
PlotLimits = [-1,2]

OutputFileTypes = [] #single value or lists of values for file outputs 

                                    #values for calculating eigenvector projections
                                    #in beta, leave set to False for now.
EigenvectorWeightMap = False
EigenvectorProjectionMaps = False
CalculateDos = False
