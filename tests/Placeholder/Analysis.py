import os, math, re, sys
import struct
import numpy as np
from matplotlib import pyplot as plt

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
v = np.zeros((LX,LY,LZ,ndim))
for t in range(tstart,tend+1,tinc):
    print("t=%s"%t)
    t_file =t+t_zero

    #file_name = "data/"+"MassSink_t%li.mat"%t_file
    #file_name = "data/"+"Pressure_t%li.mat"%t_file
    #file_name = "data/"+"Humidity_t%li.mat"%t_file
    file_name = "data/"+"OrderParameter_t%li.mat"%t_file
    #file_name = "data/"+"ChemicalPotential_t%li.mat"%t_file
    #file_name = "data/"+"BoundaryLabels_t%li.mat"%t_file

    File = open(file_name, 'rb')
    
    file_name = "data/"+"Velocity_t%li.mat"%t_file

    File2 = open(file_name, 'rb')

    #file_name = "data/"+"GradientHumidity_t%li.mat"%t_file

    #File3 = open(file_name, 'rb')

    #file_name = "data/"+"Density_t%li.mat"%t_file

    #File3 = open(file_name, 'rb')
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

    dat=File.read()
    rho = np.ndarray((LX,LY,LZ),'=d',dat,0,(8*LZ*LY,8*LZ,8))

    dat=File2.read()
    v = np.ndarray((LX,LY,LZ,ndim),'=d',dat,0,(ndim*8*LZ*LY,ndim*8*LZ,ndim*8,8))

    File.close()
    File2.close()
    
    fig,ax=plt.subplots(1,1,figsize=(6,6))

    output = "%s/component_plot_%012d.png"%(outDirName,t)
    #rho3=2*0.01*(rho2-0.2)*(rho2-1)*(2*rho2-0.2-1)-0.0128*rho4
    rgbv = np.zeros((LY,LX))
    rgbv[:,:] = np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose()
    #rgbv[:,:,1] = np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose()
    #rgbv[:,:,2] = np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose()
    
    #im=ax.imshow(np.flip(rho0.take(indices=slicepos,axis=sliceaxis)).transpose(),interpolation='nearest',origin='upper')

    im=ax.imshow(rgbv,interpolation='nearest',origin='upper')
    #im=ax.imshow(np.sqrt((gh.take(indices=0,axis=3).take(indices=slicepos,axis=sliceaxis))**2+(gh.take(indices=1,axis=3).take(indices=slicepos,axis=sliceaxis))**2),interpolation='nearest',origin='upper')
    #im=ax.imshow(np.sqrt((gh.take(indices=0,axis=2).take(indices=slicepos,axis=sliceaxis))**2),interpolation='nearest',origin='upper')
    #im=ax.imshow(np.sqrt((v.take(indices=0,axis=3).take(indices=slicepos,axis=sliceaxis))**2+(v.take(indices=1,axis=3).take(indices=slicepos,axis=sliceaxis))**2),interpolation='nearest',origin='upper')
    #print(np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose()[70,70])
    #ax.scatter(70,70)
    step=1
    X,Y=np.meshgrid(np.linspace(0,LX-1,int((LX)/step)),np.linspace(0,LY-1,int((LY)/step)))
    #print(np.sum(np.sqrt((gh.take(indices=0,axis=3).take(indices=slicepos,axis=sliceaxis))**2+(gh.take(indices=1,axis=3).take(indices=slicepos,axis=sliceaxis))**2)))
    ax.quiver(X.T,Y.T,np.flip(v[:,:,:,0].take(indices=slicepos,axis=sliceaxis)),np.flip(-v[:,:,:,1].take(indices=slicepos,axis=sliceaxis)),width=0.001,headwidth=2.5,headlength=1.5)
    print(np.sum(rho),1/6/(49)*np.log(1/(1-0.5)))
    fig.colorbar(im)
    #ax.scatter(49,49)
    plt.savefig(output, dpi=200, format='png')
    plt.close(fig)
    plt.figure()
    plt.plot(rho[:,int(LY/2),0])
    plt.savefig("test_%012d.png"%(t), dpi=200, format='png')

    print(np.amax(v[int(LX/2),:,0]))



plt.figure()
plt.plot(v[int(LX/2),:,0])
plt.savefig("test.png", dpi=200, format='png')