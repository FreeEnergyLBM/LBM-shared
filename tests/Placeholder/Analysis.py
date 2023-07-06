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
os.system("mkdir -p %s"%outDirName)

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
    rho0 = np.zeros((LX,LY,LZ))
    rho = np.zeros((LX,LY,LZ))
    rho2 = np.zeros((LX,LY,LZ))
    rho4 = np.zeros((LX,LY,LZ))
    v = np.zeros((LX,LY,LZ,ndim))

    for k in range(0,NLatt,1):
        (xk,yk,zk) = coord_k(k,LY,LZ)
        #rho0[xk,yk,zk] = struct.unpack('=d', File.read(8))[0]
        rho[xk,yk,zk] = struct.unpack('=d', File.read(8))[0]
        
        #rho2[xk,yk,zk] = struct.unpack('=d', File3.read(8))[0]
        #rho4[xk,yk,zk] = struct.unpack('=d', File4.read(8))[0]
        for i in range(ndim):
            v[xk,yk,zk,i] = struct.unpack('=d', File2.read(8))[0]
            #print(ndim)

    
    File.close()
    
    fig,ax=plt.subplots(1,1,figsize=(6,6))

    output = "%s/component_plot_%012d.png"%(outDirName,t)
    #rho3=2*0.01*(rho2-0.2)*(rho2-1)*(2*rho2-0.2-1)-0.0128*rho4
    im=ax.imshow(np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose(),interpolation='nearest',origin='upper')
    #ax.imshow((v.take(indices=0,axis=3).take(indices=slicepos,axis=sliceaxis)),interpolation='nearest',origin='upper')
    #print(np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose()[70,70])
    #ax.scatter(70,70)
    step=1
    X,Y=np.meshgrid(np.linspace(0,LX-1,int((LX)/step)),np.linspace(0,LY-1,int((LY)/step)))

    ax.quiver(X.T,Y.T,np.flip(-v[:,:,:,0].take(indices=slicepos,axis=sliceaxis)),np.flip(v[:,:,:,3-sliceaxis].take(indices=slicepos,axis=sliceaxis)),width=0.005,headwidth=6.5,headlength=2.5)
    fig.colorbar(im)
    #ax.scatter(49,49)
    plt.savefig(output, dpi=200, format='png')
    plt.close(fig)
    print(np.amax(v))
    print(np.sum(rho))


plt.figure()
plt.plot(v[5,:,0,0])
plt.savefig("test.png", dpi=200, format='png')