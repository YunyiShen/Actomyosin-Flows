import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import matplotlib.gridspec as gridspec
from matplotlib import colors
from scipy.ndimage import zoom
import os
from tqdm import tqdm
from matplotlib.animation import FuncAnimation, PillowWriter
from celluloid import Camera

from matplotlib.colors import LinearSegmentedColormap

cmap = LinearSegmentedColormap.from_list(
    "my_cmap",
    ["cyan", "magenta"]  # start → end color
)


visc_range = [4000, 20000]
for drag_range in [
                           #[200, 500],
                           #[200, 1000], 
                           #[200, 3000],
                           #[500, 3000],  
                           
                           
                            #[2000, 2000],
                            #[3000, 3000],
                            #[5000, 5000],
                            #[2000, 10000],
                            #[3000, 10000],
                            #[5000, 10000]
                            
                            
                            #[20000, 20000],
                            #[30000, 30000],
                            #[50000, 50000],
                            #[20000, 100000],
                            #[30000, 100000],
                            #[50000, 100000]
                            
                            [0,0]
                            
                            #[5, 20],
                            #[50, 200],
                            #[4000, 20000],
                            #[8000, 40000]
                            
                           #[500, 5000],
                           #[1000, 5000],
                           
                           ]:

    stress_max = 1000.0
    cell_radius = 0.5
    dt = 0.004
    dx =0.01098901098901099 #0.012345679012345678
    N = 91#81
    tmax = 600
    chi_thr = 0.2

    res_file = f"./simulations/modelcell2D_Stokes_maxstress{stress_max}_drag{drag_range[0]}-{drag_range[1]}_size{cell_radius}_visc{visc_range[0]}-{visc_range[1]}_dt{dt}_dx{dx}_tmax{tmax}.npz"
    if not os.path.exists(res_file):
        print(f"{res_file} does not exist")
        continue
    
    simulation = np.load(res_file)
    #breakpoint()
    u, v, p, stress_ext_save, t = simulation['u'], simulation['v'], simulation['p'], simulation['stress_ext']/1000, simulation['t']
    chi = simulation['chi']
    u[:,chi > 0.5*chi_thr] = np.nan
    v[:, chi > 0.5*chi_thr] = np.nan

    X, Y, N = simulation['x'], simulation['y'], N
    #stress_ext_save[:, chi > 0.2] = np.nan

    #u[chi > 0.8] = np.nan
    nx, ny = N, N
    n_frame = len(t)
    n_plot_time_series = 60
    u_ts = u[::(n_frame//n_plot_time_series)]
    v_ts = v[::(n_frame//n_plot_time_series)]
    stress_ext_ts = stress_ext_save[::(n_frame//n_plot_time_series)]
    #breakpoint()
    t_ts = t[::(n_frame//n_plot_time_series)]
    
    # fig to hold frames
    fig = plt.figure(figsize=(3*3.5, 9.0))
    gs = gridspec.GridSpec(3, 5, width_ratios=[1, 1, 0.03, 1, 0.03])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    cax1 = fig.add_subplot(gs[0, 2])
    ax2 = fig.add_subplot(gs[0, 3])
    ax3 = fig.add_subplot(gs[1, 1])
    cax2 = fig.add_subplot(gs[1, 2])
    ax4 = fig.add_subplot(gs[1, 3])
    ax5 = fig.add_subplot(gs[2, 1])
    cax3 = fig.add_subplot(gs[0, 4])
    cax4 = fig.add_subplot(gs[2, 2])
    cax5 = fig.add_subplot(gs[1, 4])
    ax6 = fig.add_subplot(gs[2, 3])
    cax6 = fig.add_subplot(gs[2, 4])
    camera = Camera(fig)
    for tt in tqdm(range(u_ts.shape[0])):
        U = u_ts[tt]
        V = v_ts[tt]
        X2 = X.reshape((nx, ny))
        Y2 = Y.reshape((nx, ny))
        U2 = U.reshape((nx, ny))
        V2 = V.reshape((nx, ny))
        stress_ext = stress_ext_ts[tt].reshape((nx, ny))
        vel_size = np.linalg.norm(np.stack((U2, V2), axis = 0), axis = 0)
    
        u_size = np.abs(U2)
        v_size = np.abs(V2)
    
        #breakpoint()
        #### three panel plot ###

        


        # interpolate a bit
        stressfine  = zoom(stress_ext, 6, order=1)   
        chifine = zoom(chi.reshape(N, N), 6, order=1)
        vel_sizef = zoom(vel_size, 6, order = 1)

        Xf = zoom(X2, 6, order=1)
        Yf = zoom(Y2, 6, order=1)
    
        u_sizef = zoom(u_size, 6, order = 1)
        v_sizef = zoom(v_size, 6, order = 1)
    
        u_f = zoom(U2, 6, order = 1)
        v_f = zoom(V2, 6, order = 1)

        mask = (chifine > 0.1) 
    
        stress_masked = np.ma.array(stressfine, mask=mask)
        vel_masked = np.ma.array(vel_sizef, mask = mask)
        u_sizef_masked = np.ma.array(u_sizef, mask = mask)
        v_sizef_masked = np.ma.array(v_sizef, mask = mask)
        u_f_masked = np.ma.array(u_f, mask = mask)
        v_f_masked = np.ma.array(v_f, mask = mask)
    
    
    

        cf = ax0.contourf(Xf, Yf, stress_masked, 
                  corner_mask=True, antialiased=True,
            levels=np.linspace(0, stress_max/1000,  21), 
            cmap='viridis')

        ax0.set_xticks([])
        ax0.set_yticks([])
        ax0.set_axis_off()
        ax0.set_aspect('equal')
        ax0.set_title("stress")


        ## overlay 
        
        cf = ax1.contourf(Xf, Yf, stress_masked, 
                  corner_mask=True, antialiased=True,
                  vmin=0,
                  vmax=stress_max/1000,
            levels=np.linspace(0, stress_max/1000,  21), 
            cmap='viridis')
        #breakpoint()

        ax1.quiver(
            X2[::4, ::4], Y2[::4, ::4],
            U2[::4, ::4], V2[::4, ::4],
            scale=2e-3,
            width=0.01,
            color='black'
        )
        ax1.set_xticks([])
        ax1.set_yticks([])
        ax1.set_axis_off()

        ax1.set_aspect('equal')
        ax1.set_title("overlay")  


        import matplotlib





        
        #cbar = fig.colorbar(cf, cax=cax1, extend="both")
        norm = matplotlib.colors.Normalize(vmin=0,
                  vmax=stress_max/1000)
        sm = matplotlib.cm.ScalarMappable(norm=norm, cmap="viridis")
        sm.set_array([])  # required by colorbar
        cbar = fig.colorbar(sm, cax=cax1, ticks=np.linspace(0, stress_max/1000, 6))
        cbar.set_label('Contractile stress by actomyosin (Pa)') 

    
        ## vel 
        
        cf = ax2.contourf(Xf, Yf, vel_masked * 1000 * 60, 
                  corner_mask=True, antialiased=True,
                  vmin=0,
                  vmax=np.nanmax(vel_masked * 1000 * 60),
            levels=np.linspace(0, np.nanmax(vel_masked * 1000 * 60),  21), 
            #linewidths=0,
            cmap=cmap)
    
        ax2.set_xticks([])
        ax2.set_yticks([])
        ax2.set_axis_off()

        ax2.set_aspect('equal')
        ax2.set_title("velocity")
    
        
        #cbar = fig.colorbar(cf, cax=cax3, extend="both")
        norm = matplotlib.colors.Normalize(vmin=0,
                  vmax=np.nanmax(vel_masked * 1000 * 60))
        sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])  # required by colorbar
        cbar = fig.colorbar(sm, cax=cax3, ticks=np.linspace(0, np.nanmax(vel_masked * 1000 * 60), 6))
        cbar.set_label('Velocity (µm/min)') 
    
    
        ########################
        ###### x and y #########
        ########################
        ## vel 
        
        cf = ax3.contourf(Xf, Yf, u_sizef_masked * 1000 * 60, 
                  corner_mask=True, antialiased=True,
                  vmin=0,
                  vmax=np.nanmax(vel_masked * 1000 * 60),
            levels=np.linspace(0, np.nanmax(vel_masked * 1000 * 60),  21), 
            #linewidths=0,
            cmap=cmap)
    
        ax3.set_xticks([])
        ax3.set_yticks([])
        ax3.set_axis_off()

        ax3.set_aspect('equal')
        ax3.set_title("x velocity, size")
    
        
        #cbar = fig.colorbar(cf, cax=cax2, extend="both")
        norm = matplotlib.colors.Normalize(vmin=0,
                  vmax=np.nanmax(vel_masked * 1000 * 60))
        sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])  # required by colorbar
        cbar = fig.colorbar(sm, cax=cax2, ticks=np.linspace(0, np.nanmax(vel_masked * 1000 * 60), 6))
        cbar.set_label('Velocity (µm/min)') 
    
        
        cf = ax4.contourf(Xf, Yf, v_sizef_masked * 1000 * 60, 
                  corner_mask=True, antialiased=True,
                  vmin=0,
                  vmax=np.nanmax(vel_masked * 1000 * 60),
            levels=np.linspace(0, np.nanmax(vel_masked * 1000 * 60),  21), 
            #linewidths=0,
            cmap=cmap)
    
        ax4.set_xticks([])
        ax4.set_yticks([])
        ax4.set_axis_off()

        ax4.set_aspect('equal')
        ax4.set_title("y velocity, size")
    
        
        #cbar = fig.colorbar(cf, cax=cax5, extend="both")
        norm = matplotlib.colors.Normalize(vmin=0,
                  vmax=np.nanmax(vel_masked * 1000 * 60))
        sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])  # required by colorbar
        cbar = fig.colorbar(sm, cax=cax5, ticks=np.linspace(0, np.nanmax(vel_masked * 1000 * 60), 6))
        cbar.set_label('Velocity (µm/min)')
    
        #ax0.set_xticks([])
        #ax0.set_yticks([])
    
        ########### with sign ##########
    
    
        
        cf = ax5.contourf(Xf, Yf, u_f_masked * 1000 * 60, 
                  corner_mask=True, antialiased=True,
                  vmin=-np.nanmax(vel_masked * 1000 * 60),
                  vmax=np.nanmax(vel_masked * 1000 * 60),
            levels=np.linspace(-np.nanmax(vel_masked * 1000 * 60), 
                               np.nanmax(vel_masked * 1000 * 60),  21), 
            #linewidths=0,
            cmap=cmap)
    
        ax5.set_xticks([])
        ax5.set_yticks([])
        ax5.set_axis_off()

        ax5.set_aspect('equal')
        ax5.set_title("x velocity")
    
        
        #cbar = fig.colorbar(cf, cax=cax4, extend="both")
        norm = matplotlib.colors.Normalize(vmin=-np.nanmax(vel_masked * 1000 * 60),
                  vmax=np.nanmax(vel_masked * 1000 * 60))
        sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])  # required by colorbar
        cbar = fig.colorbar(sm, cax=cax4, ticks=np.linspace(-np.nanmax(vel_masked * 1000 * 60), 
                                                            np.nanmax(vel_masked * 1000 * 60), 6))
        cbar.set_label('Velocity (µm/min)') 
    
        
        cf = ax6.contourf(Xf, Yf, v_f_masked * 1000 * 60, 
                  corner_mask=True, antialiased=True,
                  vmin=-np.nanmax(vel_masked * 1000 * 60),
                  vmax=np.nanmax(vel_masked * 1000 * 60),
            levels=np.linspace(-np.nanmax(vel_masked * 1000 * 60), 
                               np.nanmax(vel_masked * 1000 * 60),  21), 
            #linewidths=0,
            cmap=cmap)
    
        ax6.set_xticks([])
        ax6.set_yticks([])
        ax6.set_axis_off()

        ax6.set_aspect('equal')
        ax6.set_title("y velocity")
    
        
        #cbar = fig.colorbar(cf, cax=cax6, extend="both")
        norm = matplotlib.colors.Normalize(vmin=-np.nanmax(vel_masked * 1000 * 60),
                  vmax=np.nanmax(vel_masked * 1000 * 60))
        sm = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])  # required by colorbar
        cbar = fig.colorbar(sm, cax=cax6, ticks=np.linspace(-np.nanmax(vel_masked * 1000 * 60), np.nanmax(vel_masked * 1000 * 60), 6))
        cbar.set_label('Velocity (µm/min)')
    
    
    

        fig.tight_layout()
        #fig.savefig(f"./Fig_6/modelcell2D_Stokes_maxstress{stress_max}_drag{drag_range[0]}-{drag_range[1]}_size{cell_radius}_visc{visc_range[0]}-{visc_range[1]}_dt{dt}_dx{dx}_N{N}_tmax{tmax}.pdf")
        camera.snap()
    animation = camera.animate()
    animation.save(f"./Fig_6/modelcell2D_Stokes_maxstress{stress_max}_drag{drag_range[0]}-{drag_range[1]}_size{cell_radius}_visc{visc_range[0]}-{visc_range[1]}_dt{dt}_dx{dx}_N{N}_tmax{tmax}.mp4", writer='ffmpeg', fps=20)