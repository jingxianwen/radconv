#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function, absolute_import, unicode_literals
import numpy as np
import matplotlib.pyplot as pl
import matplotlib.animation as animation
from matplotlib.figure import Figure
import subprocess
import os, sys, traceback
PATH = os.path.dirname(os.path.abspath(__file__))

class array(np.ndarray):
  '''
  A custom subclass of numpy ndarray with some extra
  attributes.
  
  '''
  
  def __new__(cls, input_array, unit = None, name = None, var = None):
    # Input array is an already formed ndarray instance
    # We first cast to be our class type
    obj = np.asarray(input_array).view(cls)
    # add the new attribute to the created instance
    obj.unit = unit
    obj.name = name
    obj.var = var
    # Finally, we must return the newly created object:
    return obj

  def __array_finalize__(self, obj):
    # see InfoArray.__array_finalize__ for comments
    if obj is None: return
    self.var = getattr(obj, 'var', None)
    self.unit = getattr(obj, 'unit', None)
    self.name = getattr(obj, 'name', None)
  
  def __array_wrap__(self, out_arr, context = None):
    # Call the parent
    return np.ndarray.__array_wrap__(self, out_arr, context)
  
  @property
  def label(self):
    if self.unit is not None:
      return '%s (%s)' % (self.name, self.unit)
    else:
      return '%s' % (self.name)

class RadConv(object):
  '''
  
  '''
  
  def __init__(self, **kwargs):
    '''
    
    '''
    
    # Initialize params
    self.reset()
    
    # Update with user values 
    for k in kwargs.keys():
      if k.lower() not in self.__dict__.keys():
        kwargs.pop(k, None)
    self.__dict__.update(kwargs)
    
    # Separate params into radiation and radconv params (for the namelists)
    self._rkeys =  ['solar_constant', 'del_sol', 'del_sw', 
                    'atm_abs', 'sw_diff', 'linear_tau', 'albedo', 
                    'window', 'wv_exponent', 'solar_exponent', 'wv_tau', 
                    'reference_slp', 'p_broad']
    self._rckeys = ['nlayer', 'niter', 'delta_t', 't_init', 'q_init', 'ocean_depth', 
                    'surf_k', 't_eps', 'q_eps', 'gridpower',
                    'moistconv', 'thin', 'conv_iter', 'max_dt']

  def __iter__(self):
    '''
    
    '''
    
    for name in self._arrays.keys():
      yield getattr(self, name)
  
  def __getitem__(self, name):
    '''
    
    '''
    
    for k, v in self._arrays.items():
      if name.lower() == v[0].lower():
        return getattr(self, k)
    raise ValueError('No array named `%s`.' % name)
  
  def reset(self):
    '''
    
    '''
    
    # Radiation parameters
    self.solar_constant =     1360.
    self.del_sol =            1.4
    self.del_sw =             0.0
    self.atm_abs =            0.0
    self.sw_diff =            0.0
    self.linear_tau =         0.1
    self.albedo =             0.3
    self.window =             0.0
    self.wv_exponent =        4.0
    self.solar_exponent =     4.0
    self.wv_tau =             1.0
    self.reference_slp =      1.e5
    self.p_broad =            False

    # Radiative/convective parameters
    self.gridpower =          1.
    self.nlayer =             50
    self.niter =              10000
    self.delta_t =            1.
    self.t_init =             280.
    self.t_eps =              0.01
    self.q_init =             0.1
    self.q_eps =              0.1
    self.ocean_depth =        0.001
    self.surf_k =             1.
    self.moistconv =          True
    self.thin =               10
    self.conv_iter =          20
    self.max_dt =             5.
    
    # List of all arrays
    self._arrays = {
                    't':            ('Temperature', 'K'), 
                    'p':            ('Pressure', 'bar'), 
                    'q':            ('Q', None),
                    'convergence':  ('Convergence', None),
                    'time':         ('Time', 's'),
                    'dtrans':       ('dtrans', None),
                    'up':           ('Upward flux', None),
                    'down':         ('Downward flux', None),
                    'net':          ('Net flux', None),
                    'solar_down':   ('Downward Solar flux', None),
                    'flux_rad':     ('Radiative flux', None),
                    'flux_sw':      ('Shortwave flux', None),
                    'b':            ('b', None),
                    'lw_tau':       ('Longwave Optical Depth', None),
                    'tdt_rad':      ('tdt_rad', None),
                    'tdt_sw':       ('tdt_sw', None),
                    'rain':         ('Rain Profile', None),
                    't_ad':         ('Dry Adiabat', 'K'),
                    'ts':           ('Surface temperature', 'K'),
                    'ps':           ('Surface pressure', 'bar'),
                    'ph':           ('Half Pressure', 'bar'),
                    't_moistad':    ('Moist adiabat', 'K')
                   }
    
    # Initialize them
    for key, val in self._arrays.items():
      setattr(self, key, array(np.empty(0), var = key, name = val[0], unit = val[1]))
    
  def run(self):
    '''
  
    '''
        
    # Create the namelist files
    with open(os.path.join(PATH, 'radconv.nml'), 'w') as f:
      radconv_nml = ', '.join(['%s = %s' % (k,v) for k,v in self.__dict__.items() if k.lower() in self._rckeys])
      radconv_nml.replace('True', '.true.')
      radconv_nml.replace('False', '.false.')
      print('&radconv_nml %s/' % radconv_nml, file = f)
    with open(os.path.join(PATH, 'radiation.nml'), 'w') as f:
      radiation_nml = ', '.join(['%s = %s' % (k,v) for k,v in self.__dict__.items() if k.lower() in self._rkeys])
      radiation_nml.replace('True', '.true.')
      radiation_nml.replace('False', '.false.')
      print('&radiation_nml %s/' % radiation_nml, file = f)
  
    # Call the fortran code
    print("Running...")
    subprocess.call([os.path.join(PATH, 'radconv.out')], cwd = PATH)
    print("Done!")

    self.load()

    return 
  
  def load(self):
    '''
    
    '''
    
    # Load the results
    for arr in self:
      
      try:
        file = np.loadtxt(os.path.join(PATH, 'output', '%s.dat' % arr.var))
      except:
        import pdb; pdb.set_trace()
      setattr(self, arr.var, array(file, 
            var = arr.var, name = arr.name, unit = arr.unit))
  
  def animate(self, *arrs):
    '''
    
    '''
    
    # Get array
    if len(arrs) == 0:
      arrs = [self.t]
    
    curve = [None for a in arrs]
    convlabel = [None for a in arrs]
    
    # Plot the animation
    if len(arrs) == 1:
      width, height = 8, 6
    elif len(arrs) == 2:
      width, height = 12, 5
    else:
      width, height = 15, 5
    fig, ax = pl.subplots(1, len(arrs), figsize = (width,height))
    ax = np.atleast_1d(ax)
    iterlabel = pl.suptitle('Iteration 0', fontsize = 16)
    
    for i, arr in enumerate(arrs):
      convlabel[i] = ax[i].annotate('$\Delta$ = %.2f' % np.inf, xy = (0.95, 0.965), ha = 'right', va = 'top', 
                                    xycoords = 'axes fraction', fontsize = 12)
      # Profiles
      if arr.shape[1] == self.p.shape[1]:
        ax[i].plot(arr[0], self.p[0], color = 'b', ls = '--', alpha = 0.25)
        curve[i], = ax[i].plot(arr[0], self.p[0], color = 'b')
      else:
        ax[i].plot(arr[0], self.ph[0], color = 'b', ls = '--', alpha = 0.25)
        curve[i], = ax[i].plot(arr[0], self.ph[0], color = 'b')
    
      # Special case for temperature
      if arr.var == 't':
        # Show the surface
        surft = ax[i].text(self.ts[0], self.ps[0], 'S', ha = 'center', va = 'center', fontsize = 8, color = 'k')
        # Show the adiabat
        adiabat, = ax[i].plot(self.t_ad[0], self.p[0], color = 'r', alpha = 0.25)
        # Show the moist adiabat
        moistad, = ax[i].plot(self.t_moistad[0], self.p[0], color = 'g', alpha = 0.25) 
        label1 = ax[i].text(self.t_ad[0,self.nlayer // 2], self.p[0, self.nlayer // 2], 
                            'Dry adiabat', fontsize = 8, color = 'r', ha = 'right', va = 'top',
                            alpha = 0.5, clip_on = True)
        label2 = ax[i].text(self.t_moistad[0,self.nlayer // 3], self.p[0,self.nlayer // 3], 
                            'Moist adiabat', fontsize = 8, color = 'g', ha = 'right', va = 'top',
                            alpha = 0.5, clip_on = True)
      ax[i].set_xlabel(arr.label, fontsize = 14)
      ax[i].set_ylabel('Pressure (bar)', fontsize = 14)

    def updatefig(i):
      iterlabel.set_text('Iteration %d' % i)
      
      for j, arr in enumerate(arrs):
        convlabel[j].set_text(r'$\Delta$ = %.2f' % self.convergence[i])
        curve[j].set_xdata(arr[i])
        if arr.shape[1] == self.p.shape[1]:
          curve[j].set_ydata(self.p[i])
        else:
          curve[j].set_ydata(self.ph[i])
        if arr.var == 't':
          surft.set_position((self.ts[i],self.ps[i]))
          adiabat.set_xdata(self.t_ad[i])
          adiabat.set_ydata(self.p[i])
          moistad.set_xdata(self.t_moistad[i])
          moistad.set_ydata(self.p[i])
          label1.set_position((self.t_ad[i,self.nlayer // 2], self.p[i,self.nlayer // 2]))
          label2.set_position((self.t_moistad[i,self.nlayer // 3], self.p[i,self.nlayer // 3]))
      
      return curve
      
    anim = animation.FuncAnimation(fig, updatefig, frames = len(self.p), interval=50, repeat = True)
    
    for i, arr in enumerate(arrs):
      ax[i].invert_yaxis()
      ax[i].set_ylim(1.05,-0.05)
      if arr.var == 't':
        pad = max(0.75, 0.75 * (arr[1:].max() - arr[1:].min()))
      else:
        pad = max(0.1, 0.1 * (arr[1:].max() - arr[1:].min()))
      xmin = arr[1:].min() - pad
      xmax = arr[1:].max() + pad    
      ax[i].set_xlim(xmin, xmax)
    
    pl.show()