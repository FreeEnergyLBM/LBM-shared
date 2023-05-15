import os, math, re, sys
import struct
import numpy as np
import matplotlib.pyplot as plt

HeaderFile = open("data/Header.mat", 'rb')

LX=struct.unpack('=i', HeaderFile.read(4))[0]

LY=struct.unpack('=i', HeaderFile.read(4))[0]

LZ=struct.unpack('=i', HeaderFile.read(4))[0]

ndim=struct.unpack('=i', HeaderFile.read(4))[0]

t_zero = 0
tstart = 0

tend = struct.unpack('=i', HeaderFile.read(4))[0]
tinc = struct.unpack('=i', HeaderFile.read(4))[0]

slicepos=0

sliceaxis=2
if LY==1:
    sliceaxis=1
elif LX==1:
    sliceaxis=0

print(tend)
outDirName = "figures"
os.system("mkdir %s"%outDirName)

for t in range(tstart,tend+1,tinc):
    print("t=%s"%t)
    t_file =t+t_zero

    file_name = "data/"+"OrderParameter_t%li.mat"%t_file

    File = open(file_name, 'rb')
    file_name = "data/"+"Velocity_t%li.mat"%t_file

    File2 = open(file_name, 'rb')
    print(file_name)

    def coord_k(k, LY, LZ):
        """From a k value, determines its xk, yk, and zk."""    

        xk = math.floor(k/(LY*LZ))
        yk = math.floor((k - xk*LZ*LY)/LZ)
        zk = k - xk*LZ*LY - yk*LZ
        return xk, yk, zk

    NLatt=LX*LY*LZ
    rho = np.zeros((LX,LY,LZ))
    v = np.zeros((LX,LY,LZ,ndim))

    for k in range(0,NLatt,1):
        (xk,yk,zk) = coord_k(k,LY,LZ)
        rho[xk,yk,zk] = struct.unpack('=d', File.read(8))[0]

        for i in range(ndim):
            v[xk,yk,zk,i] = struct.unpack('=d', File2.read(8))[0]
            #print(ndim)

    
    File.close()
    
    fig,ax=plt.subplots(1,2,figsize=(12,12))

    output = "%s/component_plot_%012d.png"%(outDirName,t)

    ax[0].imshow(rho.take(indices=slicepos,axis=sliceaxis).transpose(),interpolation='nearest',origin='upper')
    #ax.imshow((v.take(indices=0,axis=3).take(indices=slicepos,axis=sliceaxis)),interpolation='nearest',origin='upper')

    step=1
    X,Y=np.meshgrid(np.linspace(0,LX-1,int((LX)/step)),np.linspace(0,LY-1,int((LY)/step)))

    ax[0].quiver(X.T,Y.T,v[:,:,:,0].take(indices=slicepos,axis=sliceaxis),-v[:,:,:,3-sliceaxis].take(indices=slicepos,axis=sliceaxis),scale=0.05)
    ax[1].plot(v[:,:,:,0].take(indices=slicepos,axis=sliceaxis)[int(LX/2),:])
    plt.savefig(output, dpi=200, format='png')
    plt.close(fig)