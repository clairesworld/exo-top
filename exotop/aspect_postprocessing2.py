from vtk import *
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import scipy
import matplotlib.pyplot as plt
from scipy.integrate import trapz
from scipy.ndimage.interpolation import zoom
from scipy.fftpack import dct, idct
from scipy import ndimage
import xml.etree.ElementTree
import os
import h5py
import re
import pandas as pd
from scipy.signal import periodogram
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable

rasterized=True

def unique_rows(a):
    # Given a numpy array, return another numpy array with only the unique rows
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

def unique_rows_indices(a):
    # Given a numpy array, return another numpy array with only the unique rows
    a = np.ascontiguousarray(a)
    unique_a, indices = np.unique(a.view([('', a.dtype)]*a.shape[1]), return_index=True)
    return indices
    #return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

def reduce_dims(a, transpose=True): 
    # output array from 3D to 2D, if simulation is 2D
    if transpose:
        return a[:,:,0].T
    else:
        return a[:,:,0]    
    
def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def max_slope(x,y, which='max', plot=False, tol=1):
    # find maximum gradient or minimum gradient - only works for 1D
    grad = np.diff(y) / np.diff(x)
    if which=='max':
        g = grad
        gmax = np.max(g) 
    elif which=='min':
        g = grad
        gmax = np.min(g)
    elif which=='zero':
        g = abs(grad)
        gmax = np.min(abs(g)) 
    else:
        print('max, min, or zero?')
    i_max = np.nonzero(abs(g == gmax))[0][0] # would add one to take right hand value
    x_grad_max = x[i_max]
    y_grad_max = y[i_max]
    x_grad_max0 = x[np.nonzero(g == gmax)[0][0]-tol]
    y_grad_max0 = y[np.nonzero(g == gmax)[0][0]-tol]
    x_grad_max1 = x[np.nonzero(g == gmax)[0][0]+tol]
    y_grad_max1 = y[np.nonzero(g == gmax)[0][0]+tol]
    m = (y_grad_max1-y_grad_max0)/(x_grad_max1-x_grad_max0)
#     m = grad_max
    if plot:
        plt.figure(figsize=(4,4))
        plt.plot(x, y, c='k', label='input array')
        plt.scatter(x_grad_max, y_grad_max, c='g')
        b = y_grad_max - m*x_grad_max
        plt.plot(x, m*x + b, c='g', ls='--', lw=0.5, label='slope extremus')
        plt.legend(frameon=False)
    return m, i_max

class Aspect_Data():
    def __init__(self, directory):
        self.directory = directory
        
        xdmf_filename = directory + "solution.xdmf"
        pvd_filename = directory + "solution.pvd"

        xdmf_check = os.path.isfile(xdmf_filename)
        pvd_check = os.path.isfile(pvd_filename)

        if xdmf_check:
            self.format='hdf5'
        elif pvd_check:
            self.format='vtu'
        else:
            raise Exception('No Aspect Data found in folder '+directory)
        print("Aspect data found in", self.format, "format")
        self.mesh_file=''
        
        self.read_parameters()
        self.print_summary()

    def read_mesh(self, n):
        if not(hasattr(self, 'snames')):
            self.get_solution_filenames()
        
        #if self.format=='hdf5':
        #    mesh_filename = self.directory+"mesh-00000.h5"
        #if self.format=='vtu':
        #    mesh_filename = self.directory+"solution-09999.pvtu"
        mesh_filename = self.directory+self.mnames[n]
        self.mesh_file= self.mnames[n]
        
        print("Reading mesh from", mesh_filename)
        
        if self.format=='hdf5':
            mesh_file = h5py.File(mesh_filename,"r")
            nodes_data = mesh_file['nodes']
            coords = nodes_data[:,:]
            
        if self.format=='vtu':
            reader = vtk.vtkXMLPUnstructuredGridReader()
            reader.SetFileName(mesh_filename)
            reader.Update()
            my_vtk_array = reader.GetOutput().GetPointData().GetArray("T")
            nodes_vtk_array= reader.GetOutput().GetPoints().GetData()
            coords = vtk_to_numpy(nodes_vtk_array)
    
        # Just collect unique rows (duplication due to parallel files)
        coords = coords.round(decimals=14)
        #idx = unique_rows_indices(coords)
        u, idx = np.unique(coords, axis=0, return_index = True)
        
        ind = np.lexsort((coords[idx,2], coords[idx,1], coords[idx,0]))
        self.indices = idx[ind]
        self.coords = coords[self.indices, :]
        
        xm = self.coords[:,0]
        ym = self.coords[:,1]
        zm = self.coords[:,2]

        xu = np.unique(xm)
        yu = np.unique(ym)
        zu = np.unique(zm)

        self.x, self.y, self.z = xu, yu, zu
                
        nx, ny, nz = len(self.x), len(self.y), len(self.z)
        
        self.array_shape = (nx, ny, nz)
        xm.shape = self.array_shape
        ym.shape = self.array_shape
        zm.shape = self.array_shape

        self.xm = xm
        self.ym = ym
        self.zm = zm
        
        return self.x, self.y, self.z
    
    def get_solution_filenames(self):
        if self.format == 'hdf5':
            root_filename = self.directory+"solution.xdmf"
            
        if self.format == 'vtu':
            root_filename = self.directory+"solution.pvd"
            
        print("Reading filelist from", root_filename)
        
        snames=[]
        mnames=[]
        times=[]
        
        if self.format == 'hdf5':            
            root = xml.etree.ElementTree.parse(root_filename).getroot()

            for neighbor in root.iter('Time'):
                times.append(float(neighbor.attrib['Value']))
            
            for neighbor in root.iter('DataItem'):
                text = neighbor.text
                text = text.strip()
                if text[-1]=='T':
                    snames.append(text[:-3])
                if text[-5:]=='nodes':
                    mnames.append(text[:-7])

        if self.format == 'vtu':
            root = xml.etree.ElementTree.parse(root_filename).getroot()

            for child in root:
                for atype in child.findall('DataSet'):
                    snames.append(atype.get('file'))
                    mnames.append(atype.get('file'))
                    times.append(float(atype.get('timestep')))
        
        self.times = times
        self.snames = snames
        self.mnames = mnames

        return times, snames

    def final_step(self):
        str_f = self.get_solution_filenames()[1][-1]
        return int(re.search(r'\d+', str_f).group(0))
        
    def read_temperature(self, n):
        if not(hasattr(self, 'snames')):
            self.get_solution_filenames()
            
        if self.mesh_file != self.mnames[n]:
            self.read_mesh(n)
        
        filename = self.directory + self.snames[n]
        print("Reading temperature from", filename)
        
        if self.format == 'hdf5':
            f = h5py.File(filename,"r")
            T_data=f['T']
            T = T_data[:,0]
            
        if self.format == 'vtu':
            reader = vtk.vtkXMLPUnstructuredGridReader()
            reader.SetFileName(filename)
            reader.Update()
            my_vtk_array = reader.GetOutput().GetPointData().GetArray("T")
            T = vtk_to_numpy(my_vtk_array)
#             print('loaded shape', np.shape(T))
            
        T = T[self.indices]
        T.shape = self.array_shape
        
#         print('final shape', np.shape(T))
        return self.x, self.y, self.z, T

    def read_velocity(self, n):
        if not(hasattr(self, 'snames')):
            self.get_solution_filenames()
            
        if self.mesh_file != self.mnames[n]:
            self.read_mesh(n)
        
        filename = self.directory + self.snames[n]
        print("Reading velocity from", filename)
        
        if self.format == 'hdf5':
            f = h5py.File(filename,"r")
            vel_data=f['velocity']
            u = vel_data[:,0]
            v = vel_data[:,1]
            w = vel_data[:,2]
            
        if self.format == 'vtu':
            reader = vtk.vtkXMLPUnstructuredGridReader()
            reader.SetFileName(filename)
            reader.Update()
            my_vtk_array = reader.GetOutput().GetPointData().GetArray("velocity")
            vel_data = vtk_to_numpy(my_vtk_array)
            u = vel_data[:,0]
            v = vel_data[:,1]
            w = vel_data[:,2]
            
        u = u[self.indices]
        u.shape = self.array_shape
        
        v = v[self.indices]
        v.shape = self.array_shape
        
        w = w[self.indices]
        w.shape = self.array_shape

        return self.x, self.y, self.z, u, v, w
    
    def read_statistics(self):
        filename = self.directory + "statistics"
        print("Reading statistics from", filename)
        
#         df = pd.read_csv(filename, header=None,skiprows=None,
#                          index_col=False, delimiter=r"\s+", engine='python', comment='#')
#         self.stats_timestep = np.array(df.iloc[:, 0])
#         self.stats_time = np.array(df.iloc[:, 1])
#         self.stats_rms_velocity = np.array(df.iloc[:, 10])
#         self.stats_average_T = np.array(df.iloc[:, 13])
#         self.stats_heatflux_left = np.array(df.iloc[:, 16])
#         self.stats_heatflux_right = np.array(df.iloc[:, 17])
#         self.stats_heatflux_bottom = np.array(df.iloc[:, 18])
#         self.stats_heatflux_top = np.array(df.iloc[:, 19])
        
        data = np.genfromtxt(filename, skip_header=26)
        all_data = np.genfromtxt(filename, skip_header=26, dtype=None)
        self.stats_timestep = np.array([d[0] for d in all_data])
        self.stats_time = data[:,1]
        self.stats_rms_velocity = data[:,10]
        self.stats_average_T = data[:,13]
        self.stats_heatflux_left = data[:,16]
        self.stats_heatflux_right = data[:,17]
        self.stats_heatflux_bottom = data[:,18]
        self.stats_heatflux_top = data[:,19]
#         self.stats_sol_files = [d[-1][-14:] for d in all_data]
    
    def read_parameters(self):
        # parse the parameters file into a python dictionary
        filename = self.directory + "parameters.prm"
        print("Reading parameters from", filename)
    
        parameter_file = open(filename).readlines()

        d = {}
        keys = []

        def nested_set(dic, keys, value):
            for k in keys[:-1]:
                dic = dic[k]
            dic[keys[-1]] = value

        def num(s):
            try:
                return int(s)
            except ValueError:
                try: 
                    return float(s)
                except ValueError:
                    return s

        for line in parameter_file :
            line = line.partition('#')[0]  # remove comments
            line = line.strip()
            if line:
                if re.match("subsection(.*)", line):
                    sub_name = line[11:]
                    sub_name = sub_name.strip()
                    keys.append(sub_name)
                    nested_set(d, keys, {})

                if re.match("end(.*)", line):
                    keys.pop()            
                    
                if re.match("set(.*)", line):
                    i = line.index('=')
                    key = line[4:i].strip()
                    value = line[i+1:].strip()
                    keys.append(key)
                    nested_set(d, keys, num(value))
                    keys.pop()

        self.parameters = d
        
    def rayleigh(self):
        # Calculate the Rayleigh number from the input parameters
        p = self.parameters
        T0 = p['Boundary temperature model']['Box']['Bottom temperature']
        T1 = p['Boundary temperature model']['Box']['Top temperature']
        g = p['Gravity model']['Vertical']['Magnitude']
        cp = p['Material model']['Simple model']['Reference specific heat']
        alpha = p['Material model']['Simple model']['Thermal expansion coefficient']
        k = p['Material model']['Simple model']['Thermal conductivity']
        rho = p['Material model']['Simple model']['Reference density']
        eta = p['Material model']['Simple model']['Viscosity']
        Z = p['Geometry model']['Box']['Z extent']

        deltaT = T0-T1
        kappa = k/(rho*cp)
        Ra = (rho*g*alpha*deltaT * (Z**3))/(eta*kappa)
        
        self.Ra = Ra
        return Ra
    
    def Nu(self):
        # Nusselt number with no internal heating
        return self.stats_heatflux_top/self.parameters['Geometry model']['Box']['X extent']

    def lid_thickness(self, u, v, tol=1, cut=False, plot=True): # 2D only
        x = self.x
        y = self.y
        # get horizontal average of vertical velocity
        u = reduce_dims(u)
        v = reduce_dims(v)
        mag = np.sqrt(u**2 + v**2)
        mag_int = np.trapz(mag, x=x)
        mag_av = mag_int/(np.shape(v)[1]-1)
        if plot:
            self.plot_profile(mag, xlabel='velocity magnitude')
        if cut:
            # take upper half but need to inspect each case individually
            mag_av = mag_av[int(len(mag_av)/2):]
            y = y[int(len(y)/2):]

            
        #    
            
        # maximum gradient of averaged velocity profile
        grad = np.diff(mag_av) / np.diff(y)
        grad_max = np.min(grad) # actually want the minimum because you expect a negative slope
        i_max = np.nonzero(grad == grad_max)[0][0] # would add one to take right hand value
        x_grad_max = mag_av[i_max]
        y_grad_max = y[i_max]
        if plot:
            plt.scatter(x_grad_max, y_grad_max, c='k')

        # intersection of this tangent with y-axis
        x_vel = np.linspace(0, np.max(mag_av))
        x_grad_max0 = mag_av[np.nonzero(grad == grad_max)[0][0]-tol]
        y_grad_max0 = y[np.nonzero(grad == grad_max)[0][0]-tol]
        x_grad_max1 = mag_av[np.nonzero(grad == grad_max)[0][0]+tol]
        y_grad_max1 = y[np.nonzero(grad == grad_max)[0][0]+tol]
        m1 = (y_grad_max1-y_grad_max0)/(x_grad_max1-x_grad_max0)
    #     print('m should be',m1)
        
        #
        
#         m1, idx = max_slope(y, mag_av, maximum=False, plot=plot, tol=tol) # want most negative slope
        
        b = y_grad_max - m1*x_grad_max
        y_tan = m1*x_vel + b
        if plot:
            plt.plot(x_vel, y_tan, c='g', ls='--')

        self.D_l = b # store in object
        return b
    
    def lid_base_temperature(self, T, y_L):
        x = self.x
        y = self.y
        T = reduce_dims(T)
        T_av = np.trapz(T, x=x)/(np.shape(T)[1]-1)
        # find T at y_L
#         print('index of y_L:', find_nearest_idx(y, y_L))
        T_L = T_av[find_nearest_idx(y, y_L)]
        return T_L, y_L
        
    def internal_temperature(self, T):
        x = self.x
        y = self.y
        T = reduce_dims(T)
        T_av = np.trapz(T, x=x)/(np.shape(T)[1]-1)

        # find inflection point for max core temperature
      
        
        return T_max, y_max
        
    def dT_rh(self, T_i, T_L):
        # rheological temperature scale
        return -(T_L - T_i)
    
    def delta_rh(self, y_L, y_T_i):
        # upper thermal boundary layer thickness
        return y_L - y_T_i
    
    def vbcs(self):
        # Determine velocity boundary conditions from the input parameters
        try:
            zero_vel = self.parameters["Model settings"]["Zero velocity boundary indicators"]
        except:
            # Aspect 2.0 new syntax
            zero_vel = self.parameters["Boundary velocity model"]["Zero velocity boundary indicators"]

        top_bc = "free"
        bottom_bc = "free"
        if zero_vel.find("bottom") != -1:
            bottom_bc = "rigid"
        if zero_vel.find("top") != -1:
            top_bc = "rigid"
            
        return (top_bc, bottom_bc)
    
    def print_summary(self):
        p = self.parameters
        X = p['Geometry model']['Box']['X extent']
        Y = p['Geometry model']['Box']['Y extent']
        Z = p['Geometry model']['Box']['Z extent']

        print("Simulation in box", X, "x", Y, "x", Z, ", Rayleigh number = ", self.rayleigh())
        
        vbcs = self.vbcs()
        print("Boundaries are", vbcs[0] + "-" + vbcs[1])

    def read_dynamic_topography(self, timestep):
        filename = self.directory+"dynamic_topography." + str(timestep)
        scalefactor = 1.0/self.parameters['Material model']['Simple model']['Thermal expansion coefficient']
        
        data = np.genfromtxt(filename, delimiter=' ')

        ind = np.lexsort((data[:,2], data[:, 1], data[:, 0]))
        data = data[ind]

        x, y, z, topo = data[:,0], data[:,1], data[:,2], data[:,3]

        xu, yu, zu = np.unique(x), np.unique(y), np.unique(z)
        nx, ny, nz = len(xu), len(yu), len(zu)
        topo.shape = (nx, ny)
        topo = topo * scalefactor
        
        return xu, yu, topo
    
    def plot_profile(self, s, xlabel=None):
        # s is a 2D array
        x = self.x
        y = self.y
        fig, ax = plt.subplots(figsize=(4,4))
        try:
            a = np.trapz(s, x=x)/(np.shape(s)[1]-1)
            ax.plot(a, y, c='k')
        except ValueError:
            s = reduce_dims(s)
            a = np.trapz(s, x=x)/(np.shape(s)[1]-1)
            ax.plot(a, y, c='k')
        ax.set_xlabel(xlabel, fontsize=18)
        ax.set_ylabel('depth', fontsize=18)
        ax.set_ylim(y.min(), y.max())
        return fig, ax
        
    def plot_mesh(self, s, vlabel=None, cmap='coolwarm'):
        fig, ax = plt.subplots(figsize=(8,4))
        ax.set_xlabel('x', fontsize=18)
        ax.set_ylabel('y', fontsize=18)
#         im = ax.imshow(reduce_dims(s), origin='lower', cmap=cmap)
        im = ax.pcolormesh(self.x, self.y, reduce_dims(s), cmap=cmap)
        ax.axis('equal')
        divider = make_axes_locatable(ax)
        cax = divider.new_vertical(size="20%", pad=0.5, pack_start=True)
        fig.add_axes(cax)
        fig.colorbar(im, cax=cax, orientation="horizontal")
        cax.set_xlabel(vlabel, fontsize=18)       
        return fig, ax
        
    # end class

