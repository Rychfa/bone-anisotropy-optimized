import numpy as np

# Get voxel size
ldim = 3
# Output voxelmodel in VTK format for paraview
with open("/home/jarunanp/Documents/myCourses/19_FS_How_to_Write_Fast_Numerical_Code/projects/bonemap/ground_truth.txt", "r") as g :
    data = g.readlines()

# Open file to output
with open('fabForVTK.dat','w') as f:

    m1 = 0.0
    m2 = 0.0
    m3 = 0.0
    count = 0 
    for line in data:
        values = line.strip().split(',')
        # read in
        i = int(values[0])
        j = int(values[1])
        k = int(values[2])

        eval0 = float(values[3])
        eval1 = float(values[4])
        eval2 = float(values[5])

        evec00 = float(values[6])
        evec01 = float(values[7])
        evec02 = float(values[8])
        evec10 = float(values[9])
        evec11 = float(values[10])
        evec12 = float(values[11])
        evec20 = float(values[12])
        evec21 = float(values[13])
        evec22 = float(values[14])
        # scale the fabric tensor
        trace = eval0+eval1+eval2
        if (trace > 0):
            x = i * ldim
            y = j * ldim
            z = k * ldim

            count += 1
            #scale_factor = 3.0/trace
            eval0 = eval0
            eval1 = eval1
            eval2 = eval2
            #
            eval_min = min(eval0, eval1, eval2)
            #
            if (eval0 == eval_min):
                m1 = evec00
                m2 = evec01
                m3 = evec02
            elif (eval1 == eval_min):
                m1 = evec10
                m2 = evec11
                m3 = evec12
            elif (eval2 == eval_min):
                m1 = evec20
                m2 = evec21
                m3 = evec22

            if (m1+m2+m3 > 0.0):
               f.write('{:d} {:d} {:d} '.format(x, y, z))
               f.write('{:f} {:f} {:f}\n'.format(m1, m2, m3))
        
