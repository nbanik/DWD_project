import glob
import pickle
import numpy as np
import matplotlib
matplotlib.use('agg')
from numpy.polynomial import Polynomial
from scipy import ndimage, signal, interpolate, integrate,optimize
from galpy.orbit import Orbit
from galpy.potential import MWPotential2014, turn_physical_off, MiyamotoNagaiPotential, plotDensities,evaluateDensities,SpiralArmsPotential,vcirc
from galpy.util import bovy_conversion, save_pickles, bovy_coords, bovy_plot
import pal5_util
import SCFbar_util_new
import seaborn as sns
import astropy.units as u
from galpy import potential
from galpy.potential import DehnenSmoothWrapperPotential as DehnenWrap
from optparse import OptionParser

Ac,As=SCFbar_util_new.compute_Acos_Asin()

Mbar=10**10.

vo=220.
ro=8.

def get_options():
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    
        
                         
    parser.add_option("--npoints",dest='nr',default=10,
                      type='int',
                      help="no of points")
                      
    parser.add_option("--tage",dest='tage',default=-9.,
                      type='float',
                      help="how far back in the past")
                      
                        
    parser.add_option("--ind",dest='ind',default=0,
                      help="ind")
                      
    
    return parser
    
parser= get_options()
options,args= parser.parse_args()

nr=options.nr
ind=options.ind
tage=options.tage


barpot,nobarpot=SCFbar_util_new.MWPotentialSCFbar_grow(Mbar,Acos=Ac,Asin=As,pat_speed=39.,fin_phi_deg=27.,t_on=-5.,tgrow=2,tstream=10.)

def rho_bulge(r,r1=8.,alpha=1.8,rc=1.9,Mbulge=0.5*10**10.):
    
    def integrand(r):
        return 4.*np.pi*(r**2.)*(r1/r)**(alpha) * np.exp(-(r/rc)**2.)
    
    norm=Mbulge/(integrate.quad(integrand,0.,np.inf)[0])
         
    return (norm)*(r1/r)**(alpha) * np.exp(-(r/rc)**2.)
    
def cdf_bulge(r):
    def integrand(x):
        return 4.*np.pi*(x**2.)*rho_bulge(x)
    
    norm = integrate.quad(integrand,0.,np.inf)[0]
    out=integrate.quad(integrand,0.,r)[0]
    
    return out/norm
    
def spherical_to_cylindrical(r,theta,phi):
    R = r*np.sin(phi)
    z = r*np.cos(phi)
         
    return (R,z,theta)


rr=np.linspace(0.001,10.,1000)    
cdfrho=[cdf_bulge(rr[ii]) for ii in range(len(rr))]
icdf= interpolate.InterpolatedUnivariateSpline(cdfrho,rr,k=1)

rand_r=np.empty(nr)

for jj in range(nr):
    rand_r[jj] = icdf(np.random.uniform())
    
rand_theta=np.random.uniform(0.,2.*np.pi,nr)
rand_phi = np.arccos(1.- 2.*np.random.uniform(0.,1.,nr))

R,z,theta=spherical_to_cylindrical(rand_r,rand_theta,rand_phi)

vR=np.zeros(nr)
vz=np.zeros(nr)
vT=vcirc(nobarpot,R/8.)

t_age=np.linspace(tage,0.,1001)/bovy_conversion.time_in_Gyr(vo,ro)  
coord=[]
orbits=[]
    
for ii in range(nr):
        coord.append([R[ii]/8.,vR[ii],vT[ii],z[ii]/8.,vz[ii],theta[ii]])
        orbits.append(Orbit(coord[ii])) 
        orbits[ii].turn_physical_off()
        orbits[ii].integrate(t_age,barpot)
        
tout = np.linspace(tage,0.,10)/bovy_conversion.time_in_Gyr(vo,ro) 

for kk in tout : 
        tt = round(kk*bovy_conversion.time_in_Gyr(vo,ro),1)        
        fo=open('sampled_bulge/N{}_sample_bulge_9Gyr_timestep{}Gyr_{}.dat'.format(nr,tt,ind),'w')
        fo.write('# t   x   y   z   vx   vy   vz' + '\n')
        
        for i in range(nr):
            x1=orbits[i].x(kk)
            y1=orbits[i].y(kk)
            z1=orbits[i].z(kk)
            vx1=orbits[i].vx(kk)
            vy1=orbits[i].vy(kk)
            vz1=orbits[i].vz(kk)
            fo.write(str(kk) + "  " +str(x1) + "  " + str(y1) + "  " + str(z1) + "  " + str(vx1) + "  " + str(vy1) + "  " + str(vz1) + "\n")
            
        fo.close()



