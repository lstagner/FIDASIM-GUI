#!/usr/bin/env python

from scipy.io import netcdf
import matplotlib.pyplot as plt
import numpy as np
import os

def get_data(file,vars=[]):
    try:
        f=netcdf.netcdf_file(file,'r')
    except IOError:
        print('Error: Cannot open file'+file)
        return 0

    if vars == []: vars=f.variables.keys()
    d=dict((v,f.variables[v].data) for v in vars)
    f.close()
    return d

def get_dimensions(file,vars=[]):
    try:
        f=netcdf.netcdf_file(file,'r')
    except IOError:
        print('Cannot open file: '+file)
        return 0

    if vars == []: vars=f.variables.keys()
    d=dict((k,f.variables[k].dimensions) for k in vars)
    f.close()
    return d

def plot_npa(dir):

   runid=os.path.basename(os.path.normpath(dir))
   inputs=get_data(dir+runid+'_inputs.cdf')
 
   ##Get Grid
   x_grid=inputs['x_grid']
   y_grid=inputs['y_grid']
   z_grid=inputs['z_grid']

   ##Get viewing chords
   chan_id=inputs['chan_id']
   xlens=inputs['xlens']
   ylens=inputs['ylens']
   zlens=inputs['zlens']
   xlos=inputs['xlos']
   ylos=inputs['ylos']
   zlos=inputs['zlos']

   ##Get NPA data
   npa=get_data(dir+runid+'_npa.cdf')
   ipos=npa['ipos']
   flux=npa['flux']
   energy=npa['energy']
   counts=npa['counts']

   ##Get NPA weight flux
   npa_weight=get_data(dir+runid+'_npa_weight_function.cdf')
   wflux=npa_weight['flux']
   wenergy=npa_weight['energy']

   ##Get Neutral Density
   neut=get_data(dir+runid+'_neutrals.cdf')
   dens=neut['fdens'].sum(0).sum(0)+neut['hdens'].sum(0).sum(0)+neut['tdens'].sum(0).sum(0)+neut['halodens'].sum(0).sum(0)
 
   ##Plot chords overplotted with neutral density and birth positions
   fig, ax = plt.subplots()
   ax.plot(x_grid[0,:,:],y_grid[0,:,:],'k,');
   ax.contour(x_grid[0,:,:],y_grid[0,:,:],dens,20);
   cnt=0
   for i in range(len(chan_id)):
       if chan_id[i] == 0: continue
       ax.plot([xlos[i],xlens[i]],[ylos[i],ylens[i]],label='Chan: '+str(cnt+1))
       ax.plot(ipos[cnt,0,0:counts[cnt]],ipos[cnt,1,0:counts[cnt]],'k,',alpha=.3)       
       cnt=cnt+1

   ax.set_xlim(min(x_grid[0,0,:]),max(x_grid[0,0,:]))
   ax.set_ylim(min(y_grid[0,:,0]),max(y_grid[0,:,0]))
   ax.legend()
   ax.set_xlabel('x [cm]')
   ax.set_ylabel('y [cm]')

   ##Plot MC and WF flux
   fig2,ax2=plt.subplots(nrows=len(counts),sharex=True)
   for i in range(len(counts)):
       ax2[i].step(energy,flux[i,:],label='Chan: '+str(i+1)+' MC')
       if sum(wflux[i,:]) > 0: ax2[i].plot(wenergy,wflux[i,:],label='Chan: '+str(i+1)+' WF')
       ax2[i].set_ylabel('Energy Flux')
       ax2[i].legend()

   ax2[-1].set_xlabel('Energy [keV]')
   plt.show()

   ##delete variables to save memory
   del inputs
   del npa,neut,dens
   del npa_weight,energy,counts,flux,ipos
   del x_grid,y_grid,z_grid

def plot_fida_weights(dir):
   from matplotlib.widgets import Slider,RadioButtons
   
   runid=os.path.basename(os.path.normpath(dir))
   wght=get_data(dir+runid+'_fida_weight_function.cdf')
   energy=wght['energy']
   lam=wght['lambda']
   pitch=wght['pitch']
   rad=wght['radius']
   wfunct=np.ma.array(wght['wfunct'])
   w = wfunct == 0
   wfunct[w]=np.ma.masked

   ##Defaults
   wl=655.0
   chan=0
   print('Chan: '+str(chan))
   print('Wavelength: '+str(lam[np.argmin(np.abs(lam-wl))])+' nm')

   fig, ax = plt.subplots()
   plt.subplots_adjust(left=0.25, bottom=0.25)
   c=ax.contourf(energy,pitch,wfunct[chan,:,:,np.argmin(np.abs(lam-wl))])
   ax.set_xlabel('Energy [keV]')
   ax.set_ylabel('Pitch')

   axwl  = plt.axes([0.25, 0.15, 0.65, 0.03])
   axch = plt.axes([0.25, 0.1, 0.65, 0.03])
   swl = Slider(axwl, 'Lambda', min(lam), max(lam), valinit=wl)
   sch = Slider(axch, 'Chan', 0, len(rad), valinit=0)
   def update(val):
       wl2=swl.val
       cn=int(round(sch.val))
       print('Chan: '+str(cn))
       print('Wavelength: '+str(lam[np.argmin(np.abs(lam-wl2))])+' nm')
       ax.cla()
       ax.contourf(energy,pitch,wfunct[cn,:,:,np.argmin(np.abs(lam-wl2))],30)
       ax.set_xlabel('Energy [keV]')
       ax.set_ylabel('Pitch')
       plt.draw()

   swl.on_changed(update)
   sch.on_changed(update)

   plt.show()
