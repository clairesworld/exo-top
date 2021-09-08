import numpy as np
import matplotlib.pyplot as plt
import pyshtools
import pprint
pp = pprint.PrettyPrinter(indent=4)


def read_hoggard(filename='NGS-2015-07-01303-s13.csv', file_path='models/', to_km=False, lmax=30):
    sh = np.loadtxt(file_path+filename, delimiter=",")
    cilm = np.zeros((2, lmax+1, lmax+1))
    for i in range(sh.shape[0]):
        l = int(sh[i,0])
        m = int(sh[i,1])
        coeff = sh[i,2]
        if to_km:
            if m>0:
                cilm[1, l, m] = coeff/1e3  # 1e3 for converting from m to km
            else:
                cilm[0, l, -m] = coeff/1e3

    # Hoggard uses "fully normalized" coefficients with the Condon-Shortley phase convention (like seismologists)
    clm= pyshtools.SHCoeffs.from_array(cilm,lmax=30,normalization="ortho",csphase=-1)
    
    return cilm, clm


def do(filename='NGS-2015-07-01303-s13.csv', fig_path='figs/', save=False, file_path='models/'):
    """ run John F. Rudge's script to plot dynamic topography from Hoggard (2016)
    filename : str
        Hoggard et al. (2016) file for dynamic topography
    fig_path: str
        Location to save figures
    save : boolean
        Save figures
    """
    lmax_data = 30   # spherical harmonic degree of Hoggard model
    lmax_plot = 120  # spherical harmonic degree to use for plotting purposes

    cilm, clm = read_hoggard(filename, file_path, to_km=True, lmax=lmax_data)

#     # Mimic Figure 5 -- plot of l2 norm spectrum (not really power, out by 4pi)
#     fig1, ax1 = clm.plot_spectrum(unit="per_l",xscale='lin',yscale='log',convention="l2norm")
#     plt.xlim([0,30])
#     plt.ylabel("L2 norm (km$^2$)")
#     plt.xlabel("Degree, $l$")
#     plt.ylim(np.array([1e-2, 1e1]))
#     if save:
#         plt.savefig(fig_path+'hoggard_l2norm.pdf',bbox_inches="tight")

#     plt.figure()
#     spectrum = clm.spectrum(unit='per_lm')
#     R = 6371.0 # km
#     l = clm.degrees()
#     plt.loglog(l, 4.0*np.pi*R*R*spectrum)
#     plt.xlim(-0.5+2.0*np.pi*R/5000.0, -0.5+2.0*np.pi*R/200.0)
#     plt.ylim(1e1,1e6)
#     plt.xlabel("Degree, $l$")
#     plt.ylabel("Power spectral density (km$^2$ km$^2$)")
#     if save:
#         plt.savefig(fig_path+'hoggard_2d_PSD.pdf',bbox_inches="tight")

    # Expand onto a regular lat/lon grid for plotting
    topo = clm.expand(lmax=lmax_plot)
    data = topo.data
    lats = topo.lats()
    lons = topo.lons()

    # Aid plotting by repeating the 0 degree longitude as 360 degree longitude
    lons = np.hstack([lons,np.array([360.0])])
    v = data[:,0]
    v=v.reshape((v.shape[0],1))
    data = np.hstack([data, v])

    # Cartopy plot
    import cartopy.crs as ccrs

    data_crs = ccrs.PlateCarree()
    proj_crs = ccrs.Mollweide(central_longitude=22.5)

    fig=plt.figure(figsize=(12,7))
    ax = plt.axes(projection=proj_crs)
    ax.set_global()
    ax.coastlines()
    V = np.arange(-2000,2200,200)/1e3
    cf= ax.contourf(lons, lats, data,cmap='RdBu_r',levels=V, transform=data_crs, extend="both")
    ct = ax.contour(lons,lats,data, levels=V, colors='black', linewidths = 0.5, linestyles = 'solid', transform=data_crs)
    cbar = plt.colorbar(cf, orientation = 'horizontal', label='Dynamic topography (km)', fraction = 0.07)

    if save:
        plt.savefig(fig_path+'hoggard_dyntop.pdf',bbox_inches="tight")
     
    
    # histogram of lat/lon data
    fig=plt.figure()
    xv, yv = np.meshgrid(lons, lats, sparse=False, indexing='ij')
    proj_points = proj_crs.transform_points(x=xv, y=yv, z=data.T, src_crs=data_crs)
    plt.hist(proj_points[:,:,2].flatten(), bins=20, histtype='step', color='k')
    plt.xlabel('Dynamic topography (km)')
    plt.title('RMS = %.2f km'%np.sqrt(np.mean(proj_points[:,:,2].flatten()**2)))
    if save:
        plt.savefig(fig_path+'hoggard_dyntop_hist.pdf',bbox_inches="tight")

    # Produce cartopy plot of just degrees 0-3
    # Expand onto a regular lat/lon grid for plotting
    clm.coeffs[:,4:,:]=0.0
    topo = clm.expand(lmax=lmax_plot)
    data = topo.data
    lats = topo.lats()
    lons = topo.lons()
    
    # histogram of lat/lon data
    fig=plt.figure()
    xv, yv = np.meshgrid(lons, lats, sparse=False, indexing='ij')
    proj_points = proj_crs.transform_points(x=xv, y=yv, z=data.T, src_crs=data_crs)
    plt.hist(proj_points[:,:,2].flatten(), bins=20, histtype='step', color='k')
    plt.xlabel('Dynamic topography (km)')
    plt.title('RMS = %.2f km'%np.sqrt(np.mean(proj_points[:,:,2].flatten()**2)))
    if save:
        plt.savefig(fig_path+'hoggard_dyntop_low_degree_hist.pdf',bbox_inches="tight")

    # Aid plotting by repeating the 0 degree longitude as 360 degree longitude
    lons = np.hstack([lons,np.array([360.0])])
    v = data[:,0]
    v=v.reshape((v.shape[0],1))
    data = np.hstack([data, v])

    # Cartopy plot

    fig=plt.figure(figsize=(12,7))
    ax = plt.axes(projection=proj_crs)
    ax.set_global()
    ax.coastlines()
    V = np.arange(-2000,2200,200)/1e3
    cf= ax.contourf(lons, lats, data,cmap='RdBu_r',levels=V, transform=data_crs, extend="both")
    ct = ax.contour(lons,lats,data, levels=V, colors='black', linewidths = 0.5, linestyles = 'solid', transform=data_crs)
    cbar = plt.colorbar(cf, orientation = 'horizontal', label='Dynamic topography (km)', fraction = 0.07)

    if save:
        plt.savefig(fig_path+'hoggard_dyntop_low_degree.pdf',bbox_inches="tight")


    plt.show()
