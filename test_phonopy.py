import phonopy
import phonopy.structure.brillouin_zone
import numpy as np
from phonopy.structure.brillouin_zone import get_qpoints_in_Brillouin_zone
primitive_vectors = np.array([[1,0,0],[0,1,0],[0,0,1]])
qpoints = [[3,0.5,0],[3,0.2,0]]
a = get_qpoints_in_Brillouin_zone(primitive_vectors,qpoints)

print a
