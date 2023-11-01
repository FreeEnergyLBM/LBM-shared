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
tstart = 85000

tend = struct.unpack('=i', HeaderFile.read(4))[0]
tinc = struct.unpack('=i', HeaderFile.read(4))[0]

slicepos=0

sliceaxis=2
if LY==1:
    sliceaxis=1
elif LX==1:
    sliceaxis=0

height = np.array([])

densityG = 1
vapourFrac = 0.5
diffusion = 0.008
densityL = 1
H = 98

def mlfit(x):
    return (1-(x/np.amax(x))**2)**(1/6-0.5)

print(tend)
outDirName = "figures"
os.system("mkdir -p %s"%outDirName)
v = np.zeros((LX,LY,LZ,ndim))
for t in range(tstart,tend+1,tinc):
    print("t=%s"%t)
    t_file =t+t_zero

    file_name = "data/"+"MassSink_t%li.mat"%t_file
    #file_name = "data/"+"Pressure_t%li.mat"%t_file
    #file_name = "data/"+"Density_t%li.mat"%t_file
    #file_name = "data/"+"Humidity_t%li.mat"%t_file
    #file_name = "data/"+"OrderParameter_t%li.mat"%t_file
    #file_name = "data/"+"ChemicalPotential_t%li.mat"%t_file
    #file_name = "data/"+"BoundaryLabels_t%li.mat"%t_file

    File = open(file_name, 'rb')

    file_name = "data/"+"MassSink_t%li.mat"%t_file
    #file_name = "data/"+"Pressure_t%li.mat"%t_file
    #file_name = "data/"+"Density_t%li.mat"%t_file
    #file_name = "data/"+"Humidity_t%li.mat"%t_file
    #file_name = "data/"+"OrderParameter_t%li.mat"%t_file
    #file_name = "data/"+"ChemicalPotential_t%li.mat"%t_file
    #file_name = "data/"+"BoundaryLabels_t%li.mat"%t_file

    File0 = open(file_name, 'rb')
    
    file_name = "data/"+"Velocity_t%li.mat"%t_file

    File2 = open(file_name, 'rb')

    file_name = "data/"+"GradientHumidity_t%li.mat"%t_file

    File3 = open(file_name, 'rb')

    file_name = "data/"+"OrderParameter_t%li.mat"%t_file

    File4 = open(file_name, 'rb')

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
    rho2 = np.zeros((LX,LY,LZ))
    v = np.zeros((LX,LY,LZ,ndim))

    dat=File.read()
    rho = np.ndarray((LX,LY),'=d',dat,0,(8*LY,8))

    dat=File0.read()
    rho2 = np.ndarray((LX,LY),'=d',dat,0,(8*LY,8))
    #rho = np.ndarray((LX,LY),'=i',dat,0,(4*LY,4))

    dat=File2.read()
    v = np.ndarray((LX,LY,ndim),'=d',dat,0,(ndim*8*LY,ndim*8,8))

    dat=File3.read()
    gh = np.ndarray((LX,LY,ndim),'=d',dat,0,(ndim*8*LY,ndim*8,8))

    dat=File4.read()
    c = np.ndarray((LX,LY),'=d',dat,0,(8*LY,8))

    File.close()
    File2.close()
    File3.close()
    File4.close()
    
    fig,ax=plt.subplots(1,1,figsize=(6,6))

    output = "%s/component_plot_%012d.png"%(outDirName,t)
    #rho3=2*0.01*(rho2-0.2)*(rho2-1)*(2*rho2-0.2-1)-0.0128*rho4
    rgbv = np.zeros((LY,LX))
    rgbv[:,:] = np.flip(rho).transpose()
    #rgbv[:,:] = np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose()
    #rgbv[:,:,1] = np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose()
    #rgbv[:,:,2] = np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose()
    
    #im=ax.imshow(np.flip(rho0.take(indices=slicepos,axis=sliceaxis)).transpose(),interpolation='nearest',origin='upper')
    i2 = np.argmax(np.logical_and((c[0,:] <= 0.5),(c[0,:] > 0.1)))

    #print("sum",np.sum([0.5*(liquid[0,i+1]-liquid[0,i-1]) for i in range(i2,H-2)]))
    i1 = i2 - 1
    c1_1 = c[0,i1]
    c1_2 = c[0,i2]
    h = i1 + (c1_1 - 0.5) / (c1_1 - c1_2) + 1
    print(h)
    im=ax.imshow(rgbv,interpolation='nearest',origin='upper')
    #im=ax.imshow(np.sqrt((gh.take(indices=0,axis=2).take(indices=slicepos,axis=sliceaxis))**2+(gh.take(indices=1,axis=2).take(indices=slicepos,axis=sliceaxis))**2),interpolation='nearest',origin='upper')
    #im=ax.imshow(np.sqrt((gh.take(indices=0,axis=2))**2+(gh.take(indices=1,axis=2))**2),interpolation='nearest',origin='upper')
    #im=ax.imshow(np.sqrt((gh.take(indices=0,axis=2).take(indices=slicepos,axis=sliceaxis))**2),interpolation='nearest',origin='upper')
    #im=ax.imshow(np.sqrt((v.take(indices=0,axis=3).take(indices=slicepos,axis=sliceaxis))**2+(v.take(indices=1,axis=3).take(indices=slicepos,axis=sliceaxis))**2),interpolation='nearest',origin='upper')
    #print(np.flip(rho.take(indices=slicepos,axis=sliceaxis)).transpose()[70,70])
    #ax.scatter(70,70)
    step=1
    X,Y=np.meshgrid(np.linspace(0,LX-1,int((LX)/step)),np.linspace(0,LY-1,int((LY)/step)))
    #print(np.sum(np.sqrt((gh.take(indices=0,axis=3).take(indices=slicepos,axis=sliceaxis))**2+(gh.take(indices=1,axis=3).take(indices=slicepos,axis=sliceaxis))**2)))
    #ax.quiver(X.T,Y.T,np.flip(v[:,:,:,0].take(indices=slicepos,axis=sliceaxis)),np.flip(-v[:,:,:,1].take(indices=slicepos,axis=sliceaxis)),width=0.001,headwidth=2.5,headlength=1.5)
    print(np.sum(rho[int(LX/2),:])/(0.002/(100-h)*np.log(1/(1-0.1))))
    fig.colorbar(im)
    #ax.scatter(49,49)
    plt.savefig(output, dpi=200, format='png')
    plt.close(fig)
    #plt.figure()
    #plt.plot(rho[int(LX/2),:])
    #plt.plot(rho[:,int(LY/2),0])
    #plt.savefig("test_%012d.png"%(t), dpi=200, format='png')
    height = np.append(height, h)
    print(gh[int(LX/2),51,1])

    x=np.where(np.logical_and(rho[:,:]>=0.45, rho[:,:]<=0.55))[0]
    y=np.where(np.logical_and(rho[:,:]>=0.45, rho[:,:]<=0.55))[1]

    fig,ax=plt.subplots(1,1,figsize=(6,6))
    plt.scatter(x,rho2[x,y]/np.abs(np.amin(rho2)))
    try:
        plt.plot(np.split(x,2)[1],mlfit(np.split(x,2)[1])/np.abs(np.amin(mlfit(np.split(x,2)[1]))))
    except:
        pass
    output = "%s/ml_plot_%012d.png"%(outDirName,t)
    plt.savefig(output, dpi=200, format='png')
    plt.close(fig)



plt.figure()
plt.plot(rho[int(LX/2),:])
plt.savefig("test.png", dpi=200, format='png')


plt.figure()
plt.xlabel("Time (lattice Units)")
plt.ylabel("Height (lattice Units)")
# plt.ylim([0,1.1])

t = tinc*np.linspace(0,height.size-1,height.size)

plt.scatter(t,height[0:],label="Data")

def height_func(t, h0, k):

    a = np.array(2*k*t + (H-h0)**2)
    a[a<0] = 0
    return H - np.sqrt(a)

def height_func2(t, h0, k):
    dt = np.diff(t, prepend=[0])
    intGradVapour = np.cumsum(gradVapour * dt)
    return h0 - k * intGradVapour

#popt, _ = curve_fit(height_func, t, height, p0=[70,-1e-3])
#print(popt)
# print(height_func(t, popt[0], popt[1]))
#plt.plot(t, height_func(t, popt[0], popt[1]), 'r-', label="Fit")

k = (densityG*1/(1-vapourFrac)) *vapourFrac * diffusion/densityL

h0 = height[0]
print(h0, k)
#print(popt[1]/k)
height_an = height_func(t, h0, k)
#height_an = height_func2(t, h0, 0.5*k/vapourFrac)
#height_an2 = height_func2(t, h0, k/vapourFrac)
plt.plot(t, height_an, 'k-', label="Theory")
#plt.plot(humidity[0,:])
#plt.plot(masssink[0,:]/np.amax(masssink[0,:]))
#plt.plot(0.5-0.5*np.tanh(0.5*(np.linspace(0,humidity[0,:].size-1,humidity[0,:].size)-h)))
#plt.plot(liquid[0,:])
#plt.plot(np.linspace(height[-1],399,height.size),np.amax(humidity[0,:])*(1-np.linspace(0,1,height.size)))
#plt.scatter(height[-1],0)
#plt.plot(t, height_an2, 'r-', label="Matching mass loss")
#plt.plot([0.1483, 0.16367, 0.16935, 0.16904, 0.16935])
print((height[0]-height[-1])/(height_an[0]-height_an[-1]))
plt.legend()
plt.savefig("massloss.png", dpi=500, format='png')
