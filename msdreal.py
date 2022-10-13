import MDAnalysis as md
import math
import matplotlib.pyplot as plt
import numpy as np
import csv


startconfiguration="start.xyz"
filename="pore.xyz"
steps=200000
MSD=True


g=globals()
velocity=[]
correlation2=[]
correlation0=[]
correlation1=[]
correlationxy=[]

correlation_length=20000
shift=1
correlation_list=[]
window=steps-correlation_length
correlation_tau=[]

def converttolist():
    correlation_list=[]
    traj=md.Universe(startconfiguration,filename)
    for ks in traj.trajectory[0:steps]:
        correlation_list.append(traj.atoms.center_of_mass()[1])
#    print(correlation_list)
    return correlation_list

def msdcenterofmass2(correlation_list,correlation_length,shift,steps):
#    correlation_list,correlation_z=[], []
    window=steps-correlation_length
#    correlation_z=0
    print(window)
#    tau=1
    print(correlation_length)
    for i in range(1,correlation_length,1):
        correlation_z=0
        for tau in range(1,window,shift):
#            print("tau",tau)
#            print("i:",i)
#            print(correlation_list)
#            print(correlation_list[0])
            correlation_z=correlation_z+(math.pow(correlation_list[tau+i]-correlation_list[tau],2))
#            print(correlation_z)
        correlation_tau.append(correlation_z/(window/shift))
    return correlation_tau
#    return correlation_tau

#print(type(traj.atoms.center_of_mass()[1]))#-traj.atoms.center_of_mass()[1][tau+shift])   

def writefiles(correlation_z,correlation_length):
    x=np.arange(1,correlation_length,1)
    with open("diffusion.dat","w") as outfile:
        writer=csv.writer(outfile, delimiter="\t")
        writer.writerows(zip(x,correlation_z))

def writecoefficient(correlation_z,correlation_length):
    with open("diff_coeff.dat","w") as outfile2:
        writer=csv.writer(outfile2,delimiter=" ")
        writer.writerow("diff_z")
        diffusion=((correlation_tau[-1] - correlation_tau[10000])/(correlation_length/2))*10E-05
        writer.writerow(str(diffusion))

if MSD==True:
    correlation_list=converttolist()
    correlation_tau=msdcenterofmass2(correlation_list,correlation_length,shift,steps)#,steps,correlation_length,shift)
    x=np.arange(1,correlation_length,1)
#    print(x)
    plt.plot(x,correlation_tau)
    plt.savefig("diffusion_real.png")
    plt.show()
    diffusion=((correlation_tau[-1] - correlation_tau[10000])/(correlation_length-10000))*10E-05
    print(diffusion)
    writefiles(correlation_tau,correlation_length)
    writecoefficient(correlation_tau,correlation_length)
#def msdcenterofmass2():
#    for i in range(3):
#    traj=md.Universe(startconfiguration,filename)
#        for ks in traj.trajectory[0:steps-correlation_length]
#            current_step=traj.atoms.center_of_mass(){i]
#    for tau in range(0,steps-correlation_length,shift):
#        current_step=
#             print(type(traj.atoms.center_of_mass()[1]))#-traj.atoms.center_of_mass()[1][tau+shift])   
                
                
                
                
                
                #        print(traj.atoms.center_of_mass())




def msdcenterofmass():
    for i in range (3):
        dodo=False
        print(i)
        correlation3=0
        traj=md.Universe(startconfiguration,filename)
        for ks in traj.trajectory[0:steps]:
#        print(traj.atoms.center_of_mass()[i])
            current_step=traj.atoms.center_of_mass()[i]
#        last_step=traj.atoms.center_of_mass()[i]
            if dodo==True:
                velocity=(last_step-current_step)
                correlation3=correlation3+math.sqrt(velocity*velocity)    
                g['correlation{0}'.format(i)].append(correlation3)            
            last_step=traj.atoms.center_of_mass()[i]
            dodo=True
    return correlation0,correlation2,correlation1



def combinexydiffusion(correlation_x,correlation_y):
#    print(correlation_x,correlation_y)
    correlation_xy=[]
    for i, s in zip(correlation_x, correlation_y):
        correlation_xy.append((i+s)/2)
    return correlation_xy

def writefiles(correlation_z,correlation_length):
    x=np.arange(1,correlation_length,1)
    with open("diffusion.dat","w") as outfile:
        writer=csv.writer(outfile, delimiter="\t")
        writer.writerows(zip(x,correlation_z))

def writecoefficient(correlation_x,correlation_y,correlation_z,steps):
    with open("diff_coeff.dat","w") as outfile2:
        writer=csv.writer(outfile2,delimiter=" ")
        writer.writerow("diff_x/y"+"\t"+"diff_z")
        writer.writerow(str((calcmsd(correlation2,steps)+calcmsd(correlation0,steps))/2)+"\t"+str(calcmsd(correlation1,steps)))
        writer.writerow(str(calcdiffratio(correlation0,correlation2,correlation1,steps)))


def plotmsd(correlation_xy,correlation_z,steps):
    x=np.arange(1,steps,1)
    plt.plot(x,correlation_z,label="z")
    plt.plot(x,correlation_xy,label="x/y")
#    plt.plot(x,correlation_x,label="x")
    plt.legend(loc="upper left")
    plt.savefig("diffusion.png")
    plt.show()

def calcmsd(correlation,steps):
    diffusion_coeff=(((correlation[-1]-correlation[10000])/(steps-10000))*2E-05)
    return diffusion_coeff


def calcdiffratio(correlation0,correlation2,correlation1,steps):
    ratio=calcmsd(correlation1,steps)/((calcmsd(correlation0,steps)+calcmsd(correlation2,steps))/2)
    return ratio


#if MSD==True:
#    converttolist()
#    msdcenterofmass2()
#    plotmsd(correlationxy,correlation1,steps)
#    print("y/x:",(calcmsd(correlation2,steps)+calcmsd(correlation0,steps))/2,"z:",calcmsd(correlation1,steps))
#    print("ratio:",calcdiffratio(correlation0,correlation2,correlation1,steps))
#    correlationxy=combinexydiffusion(correlation0,correlation2)
#    print(correlationxy)
#    writefiles(correlationxy,correlation1,steps)
#    writecoefficient(correlation2,correlation0,correlation2,steps)
#    plotmsd(correlationxy,correlation1,steps)
