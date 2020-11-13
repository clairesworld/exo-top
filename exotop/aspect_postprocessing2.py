import sys
src_paths = ['/usr/lib/python36.zip', 'usr/lib/python3.6', '/usr/lib/python3.6/lib-dynload', '/usr/local/lib/python3.6/dist-packages', '/usr/lib/python3/dist-packages']
for s in src_paths:
    sys.path.insert(0, s)
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapz
import xml.etree.ElementTree
import os
import h5py
import re
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

def horizontal_mean(A, x, axis=None):
    if axis is None:
        axis = np.argmax(np.shape(A))  # assume longest axis (probably box width)
    int_x = trapz(A, x=x, axis=axis)
    int_x = int_x / ((max(x)-min(x)))
    if A.ndim == 2:
        return int_x.T
    elif A.ndim == 3:
        return int_x.T[0]

class Aspect_Data():
    def __init__(self, directory, read_parameters=True, read_statistics=False, verbose=True):
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
        if verbose:
            print("Aspect data found in", self.format, "format")
        self.mesh_file=''
        
        if read_parameters:
            self.read_parameters(verbose=verbose)
        if read_statistics:
            self.read_statistics(verbose=verbose)
        if verbose:
            self.print_summary()

    def read_mesh(self, n, verbose=True):
        if not(hasattr(self, 'snames')):
            self.get_solution_filenames(verbose=verbose)
        
        #if self.format=='hdf5':
        #    mesh_filename = self.directory+"mesh-00000.h5"
        #if self.format=='vtu':
        #    mesh_filename = self.directory+"solution-09999.pvtu"
        mesh_filename = self.directory+self.mnames[n]
        self.mesh_file= self.mnames[n]
        
        if verbose:
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
    
    def get_solution_filenames(self, verbose=True):
        if self.format == 'hdf5':
            root_filename = self.directory+"solution.xdmf"
            
        if self.format == 'vtu':
            root_filename = self.directory+"solution.pvd"
        
        if verbose:
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
        
    def read_temperature(self, n, verbose=True):
        if not(hasattr(self, 'snames')):
            self.get_solution_filenames()
            
        if self.mesh_file != self.mnames[n]:
            self.read_mesh(n, verbose=verbose)
        
        filename = self.directory + self.snames[n]
        if verbose:
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

    def read_velocity(self, n, verbose=True):
        if not(hasattr(self, 'snames')):
            self.get_solution_filenames()
            
        if self.mesh_file != self.mnames[n]:
            self.read_mesh(n, verbose=verbose)
        
        filename = self.directory + self.snames[n]
        if verbose:
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
    
    def read_statistics(self, skip_header=26, verbose=True):
        filename = self.directory + "statistics"
        if verbose:
            print("Reading statistics from", filename)
        
        data = np.genfromtxt(filename, skip_header=skip_header)
        all_data = np.genfromtxt(filename, skip_header=26, dtype=None)
        self.stats_timestep = np.array([d[0] for d in all_data])
        self.stats_time = data[:,1]
        self.stats_rms_velocity = data[:,10]
        self.stats_average_T = data[:,13]
        self.stats_heatflux_left = data[:,16]
        self.stats_heatflux_right = data[:,17]
        self.stats_heatflux_bottom = data[:,18]
        self.stats_heatflux_top = data[:,19]
        self.stats_average_viscosity = data[:,22]

    def read_stats_sol_files(self, col_vis=20, skip_header=26):
        filename = self.directory + "statistics"
        all_data = np.genfromtxt(filename, skip_header=skip_header, dtype=None)
        last = 0
        files = np.zeros(len(all_data), dtype=np.int8)
        #  find last instance that's not ""
        for n, d in enumerate(all_data):
            s = d[col_vis][-14:]
            if not s.decode()=='""':
                last = int(re.search(r'\d+', s.decode()).group(0))
            files[n] = last
        self.sol_files = files
        return files

    def find_time_at_sol(self, n=None, sol_files=None, return_indices=True, i_vis=20, skip_header=26):
        # input solution file and get first timestep - for n or all solutions
        # n is the number in the solution filename
        if sol_files is None:
            sol_files = self.read_stats_sol_files(col_vis=i_vis, skip_header=skip_header)
        u, indices = np.unique(sol_files, return_index=True)
        if return_indices:
            if n is None:
                return indices
            else:
                return indices[n]
        else:
            time = self.stats_time
            if n is None:
                return time[indices]
            else:
                return time[indices[n]]
        
    def read_parameters(self, verbose=True):
        # parse the parameters file into a python dictionary
        filename = self.directory + "parameters.prm"
        if verbose:
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
    
    def Ra_1(self):
        p = self.parameters
        Ra_0 = p['Material model']['Nondimensional model']['Ra']
        deltaeta = p['Material model']['Nondimensional model']['Viscosity temperature prefactor']
        deltaeta = np.round(np.exp(deltaeta))
        return Ra_0*deltaeta
    
    def Nu(self, k=1):
        # Nusselt number with no internal heating
        p = self.parameters
        F = self.stats_heatflux_top/p['Geometry model']['Box']['X extent']
        d = p['Geometry model']['Box']['Y extent']
        dT = p['Boundary temperature model']['Box']['Bottom temperature'] - p['Boundary temperature model']['Box']['Top temperature']
        Nu = d*F/(k*dT)
        self.Nu = Nu
        return Nu

    def ubl_thickness(self, n=None, T_l=None, T_i=None, k=1, **kwargs):
        # get upper boundary layer thickness required to produce surface heat flux (F=Nu) between T_i and T_l
        # corresponds to rheological sublayer delta_rh in Moresi & Solomatov 2000
        if T_i is None:
            T_i = self.internal_temperature(self, n=n, **kwargs)
        if T_l is None:
            T_l = self.lid_base_temperature(self, n=n, **kwargs)
        ts = self.find_time_at_sol(n)
        F = self.stats_heatflux_top[ts]/self.parameters['Geometry model']['Box']['X extent']
        return k*(T_i - T_l)/F

    def delta_0(self, delta_rh=None, delta_L=None, **kwargs):
        # total lithosphere thickness, from MS95
        if delta_rh is None:
            delta_rh = self.ubl_thickness(**kwargs)
        if delta_L is None:
            delta_L = self.lid_thickness(**kwargs)
        return delta_rh + delta_L

    def delta_Ti_MS95(self, n, Nu=None, T=None, tol=1e-5, **kwargs):
        # corresponds to the stagnant lid thickness + rheological sublayer (upper tbl)
        def delta_0(Ti, Nu):
            return Ti/Nu
        def delta_1(Ti, Nu):
            return (1-Ti)/Nu
        def T_i(d0, d1, T, x, z):  # area-averaged temperature between thermal boundary layers
            Ix = trapz(T, x=x)
            z_d0 = find_nearest_idx(z, 1 - d0)
            z_d1 = find_nearest_idx(z, d1)
            Iz = trapz(Ix[z_d0:z_d1], x=z[z_d0:z_d1])
            return 1/(1 - d0 - d1)*Iz

        x = self.x
        y = self.y
        if Nu is None:
            Nu = self.Nu(**kwargs)
        if T is None:
            _, _, _, T = self.read_temperature(n, verbose=verbose)
        Ti = np.mean(horizontal_mean(T, x))  # initial guess
        diff = 10
        while diff >= tol:
            Ti_guess = Ti
            d0 = delta_0(Ti_guess, Nu)
            d1 = delta_1(Ti_guess, Nu)
            Ti = T_i(d0, d1, T, x, z)
            diff = abs(Ti - Ti_guess)
        return Ti, d0

    def lid_thickness(self, u=None, v=None, n=None, tol=1, cut=False, plot=False, cutdiv=2, **kwargs):
        # stagnant lid thickness (mechanical boundary layer) from Moresi & Solomatov 2000
        x = self.x
        y = self.y
        if (u is None) or (v is None):
            _, _, _, u, v, _ = self.read_velocity(n, verbose=False)
        # get horizontal average of vertical velocity
        mag = np.sqrt(u**2 + v**2)
        mag_av = horizontal_mean(mag, x)
        if plot:
            self.plot_profile(mag, xlabel='velocity magnitude')
        if not cut:
            cutdiv=1
        b=1.1
        while b>1:
            # take upper half by defualt (cutdiv=2) but need to inspect each case individually
            mag_avprime = mag_av[int(len(mag_av)/cutdiv):]
            yprime = y[int(len(y)/cutdiv):]

            # maximum gradient of averaged velocity profile
            grad = np.diff(mag_avprime, axis=0) / np.diff(yprime)
            grad_max = np.min(grad) # actually want the minimum because you expect a negative slope
            i_max = np.nonzero(grad == grad_max)[0][0] # would add one to take right hand value
            x_grad_max = mag_avprime[i_max]
            y_grad_max = yprime[i_max]
            if plot:
                plt.scatter(x_grad_max, y_grad_max, c='k', label='max grad: ({:04.1f}),({:04.1f})'.format(x_grad_max, 
                                                                                                          y_grad_max))
                plt.axhline(y=np.min(yprime), alpha=0.2, c='k', ls='--')

            # intersection of this tangent with y-axis
            x_vel = np.linspace(0, np.max(mag_avprime))
            x_grad_max0 = mag_avprime[np.nonzero(grad == grad_max)[0][0]-tol]
            y_grad_max0 = yprime[np.nonzero(grad == grad_max)[0][0]-tol]
            x_grad_max1 = mag_avprime[np.nonzero(grad == grad_max)[0][0]+tol]
            y_grad_max1 = yprime[np.nonzero(grad == grad_max)[0][0]+tol]
            m1 = (y_grad_max1-y_grad_max0)/(x_grad_max1-x_grad_max0)
            b = y_grad_max - m1*x_grad_max
            if b>1: # if this doesn't work it's probably because lid base is below 50% depth, need to recut profile
                print('\n recutting')
                cutdiv = cutdiv+0.2
        y_tan = m1*x_vel + b
        if plot:  # overplot
            plt.plot(x_vel, y_tan, c='g', ls='--', label='tangent to max gradient')
            plt.legend()
        return b
    
    def lid_base_temperature(self, n=None, T=None, T_av=None, delta_L=None, u=None, v=None, cut=False, plot=False,
                             verbose=False, **kwargs):
        x = self.x
        y = self.y
        if T_av is None:
            if T is None:
                _, _, _, T = self.read_temperature(n, verbose=verbose, **kwargs)
            T_av = horizontal_mean(T, x)
        if (delta_L is None):
            if (u is None) or (v is None):
                _, _, _, u, v, _ = self.read_velocity(n, verbose=verbose, **kwargs)
            delta_L = self.lid_thickness(u=u, v=v, cut=cut, plot=plot, **kwargs)
        # find T at delta_L
        T_l = T_av[find_nearest_idx(y, delta_L)]
        return T_l
        
    def internal_temperature(self, n=None, T=None, T_av=None, plot=False, return_coords=False, **kwargs):
        # almost-isothermal temperature of core of convecting cell
        # note: MS2000 define this as the maximal horizontally-averaged temperature in the layer
        x = self.x
        y = self.y
        if T_av is None:
            if T is None:
                _, _, _, T = self.read_temperature(n, verbose=verbose)
            T_av = horizontal_mean(T, x)

        # find inflection point for max core temperature
        z = y
        f_prime = np.gradient(T_av) # differential approximation
        idx = np.where(np.diff(np.sign(f_prime)))[0] # Find the inflection point.
        y_infections = z[idx]
        T_inflections = T_av[idx]
        
        if plot:
            print ('inflection point', y_infections) 
            fig, ax = plt.subplots (figsize = (7, 7))
            ax.plot (z, T_av, 'bo-', ms = 2)
            ax.plot (y_infections, T_inflections, 'ro', ms = 5)
            ax.set_xlabel('depth')
            ax.set_ylabel('T')

        if return_coords:
            return T_inflections[-1], y_infections[-1]
        else:
            return T_inflections[-1]
        
    def dT_rh(self, T_l=None, T_i=None, **kwargs):
        # rheological temperature scale (temperature drop in unstable part of lid), e.g. eq. (A7) in SM2000
        if T_i is None:
            T_i = self.internal_temperature(self, **kwargs)
        if T_l is None:
            T_l = self.lid_base_temperature(self, **kwargs)
        return -(T_l - T_i)
    
    def T_components(self, n=None, T=None, T_i=None, T_l=None, delta_rh=None, delta_L=None, u=None, v=None, cut=False, plot=False,
                     verbose=False, **kwargs):
        # return RHS of h' \propto (dT_rh/dT_m)*(delta_u/d_m)

        x = self.x
        y = self.y
        if T is None:
            _, _, _, T = self.read_temperature(n, verbose=verbose)
        T_av = horizontal_mean(T, x)
        p = self.parameters
        d_m = p['Geometry model']['Box']['Y extent']
        dT_m = p['Boundary temperature model']['Box']['Bottom temperature'] - p['Boundary temperature model']['Box']['Top temperature']
        if delta_L is None:
            delta_L = self.lid_thickness(u=u, v=v, cut=cut, plot=plot)
        if T_i is None:
            T_i = self.internal_temperature(T_av=T_av, **kwargs)
        if T_l is None:
            T_l = self.lid_base_temperature(T_av=T_av, delta_L=delta_L, cut=cut, **kwargs)
        if delta_rh is None:
            delta_rh = self.ubl_thickness(n=n, T_l=T_l, T_i=T_i, **kwargs)
        delta_0 = self.delta_0(delta_rh=delta_rh, delta_L=delta_L)  # mechanical boundary layer MS95
        dT_rh = self.dT_rh(T_l=T_l, T_i=T_i)
        self.T_params = {'dT_rh':dT_rh, 'dT_m':dT_m, 'delta_rh':delta_rh, 'd_m':d_m, 'delta_L':delta_L, 'T_l':T_l, 'T_i':T_i, 'delta_0':delta_0, 'T_av':T_av, 'y':y}
        return self.T_params
    
    def surface_mobility(self, n=None, delta_0=None, delta_rh=None, delta_l=None, u=None, **kwargs):
        # stagnant lid criterion S <<1 from Moresi & Solomatov 2000
        if u is None:
            _, _, _, u, v, _ = self.read_velocity(n, verbose=False)
        u_0 = horizontal_mean(u, self.x)[-1]
        if delta_0 is None:
            if delta_rh is None:
                delta_rh = ubl_thickness(self, n=n, **kwargs)
            if delta_l is None:
                delta_l = lid_thickness(self, n=n, u=u, **kwargs)
            delta_0 = delta_rh + delta_l
        return delta_0**2 * u_0

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
    
    def plot_profile(self, s, n=None, xlabel='', ylabel='depth', fig=None, ax=None, **plotkwargs):
        # s is a 2D or 1D array
        try:
            x = self.x
            y = self.y
        except AttributeError:
            if n is None:
                n = self.final_step()
            self.read_mesh(n=n, verbose=False)
            x = self.x
            y = self.y
        if fig is None:
            fig, ax = plt.subplots(figsize=(4,4))
        if s.ndim == 2:  # s is not horizontally- averaged yet
            s = horizontal_mean(s, x)
        elif s.ndim > 2:  # still in 3D default shape
            s = reduce_dims(s)
            s = horizontal_mean(s, x)
        ax.plot(s, y, **plotkwargs)
        ax.set_xlabel(xlabel, fontsize=18)
        ax.set_ylabel(ylabel, fontsize=18)
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
        
    def write_ascii(self, A=None, fname='ascii', ext='.txt', path='', n=None, default_field='T', **kwargs):
        # format ascii file for initialising future ASPECT runs. A is 2D array
        if n is None:
            n = self.final_step()

        if A is not None:
            if self.mesh_file != self.mnames[n]:
                self.read_mesh(n, verbose=verbose)
            x = self.x
            y = self.y
        else:  # load default_field
            if default_field == 'T':
                x, y, _, A = self.read_temperature(n, verbose=False)

        nx = len(x)
        ny = len(y)
        header = 'POINTS: {:d} {:d}'.format(nx, ny) + '\nColumns: x y temperature'
        fpath = path+fname+ext

        xv, yv = np.meshgrid(x, y)
        out = np.zeros((nx*ny, 3))
        # out = np.vstack((xv, yv, np.zeros_like(xv)))
        A = reduce_dims(A)
        print('xv, yv, A, out', np.shape(xv), np.shape(yv), np.shape(A), np.shape(out))

        row = 0
        for jj in range(ny):
            for ii in range(nx):
                out[row, 0] = x[ii]
                out[row, 1] = y[jj]
                out[row, 2] = A[jj, ii]
                row = row + 1
                # out[row, 2] = A[jj, ii]

        np.savetxt(fpath, out, delimiter=" ", fmt="%s", header=header)
        print('writing to ascii file', fpath)

    def hello(self):
        print('                                    i am here')
    # end class