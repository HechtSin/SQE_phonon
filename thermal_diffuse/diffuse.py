import numpy as np
import matplotlib.pyplot as plt
import os

for i in np.arange(0,4.1,0.1):
    q_start = "Q_start = np.array([%.2f,0,0]) \n"%(i)
    q_end = "Q_end = np.array([%.2f,0,4]) \n"%(i)
    f = open("SQE_gauss.py")
    print ("a = ")
    print (i)
    f_out = open("temp.py","wt")
    for line in f:
        if 'Q_start =' in line:
            f_out.write(q_start)
        elif 'Q_end =' in line:
            f_out.write(q_end)
        else:
            f_out.write(line)
    f.close()
    f_out.close()
    os.system("python temp.py")
    os.system("mv temp_sqe.txt sqe_%.2f.txt" % (i))
