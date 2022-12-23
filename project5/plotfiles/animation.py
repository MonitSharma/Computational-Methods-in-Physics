import numpy as np
from numpy.ma import masked_array
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
import pyarma as pa
import sys

#Input arguments that determines what is being plotted.
filename = sys.argv[1]
z_axis_label = sys.argv[2]
slits_overlay = int(sys.argv[3])
number_of_snaps = int(sys.argv[4])

snaps = []
for i in range(number_of_snaps):
    snaps.append(float(sys.argv[i+5]))

# Set up a 2D xy grid
h = 0.005
x_points = np.arange(0, 1+h, h)
y_points = np.arange(0, 1+h, h)
x, y = np.meshgrid(x_points, y_points, sparse=True)



# Fill z_data_list with f(x,y,t)
z_data_list = []
snapshot_index_list = []
c = 0

if z_axis_label == 'Real(u)':   #if we want to plot the real part
    A = pa.cx_cube() #Create pa.mat object (just as arma::mat in C++)
    A.load("./datafiles/"+str(filename)) #Load the content of the matrix you saved into your Python program.
    # A function for a Gaussian that is travelling
    # in the x direction and broadening as time passes


    # Array of time points (Dynamically allocates T)
    dt = 0.000025
    t_points = np.arange(0, dt*(np.shape(A)[0]), dt)

    for t in t_points:
        z_data = np.rot90(np.array(np.real(A[pa.single_slice, c])))
        c += 1
        z_data_list.append(z_data)

        #Finds if the timestap matches a snapshot an appends index
        if any(t== t_snap for t_snap in snaps):
            snapshot_index_list.append(int(c-1))
    filename = filename.split('.t')[0]+'_real'

elif z_axis_label == 'Imag(u)': #if we want to plot the imaginary part
    A = pa.cx_cube() #Create pa.mat object (just as arma::mat in C++)
    A.load("./datafiles/"+str(filename)) #Load the content of the matrix you saved into your Python program.
    # A function for a Gaussian that is travelling
    # in the x direction and broadening as time passes


    # Array of time points (Dynamically allocates T)
    dt = 0.000025
    t_points = np.arange(0, dt*(np.shape(A)[0]), dt)

    for t in t_points:
        z_data = np.rot90(np.array(np.imag(A[pa.single_slice, c])))
        c += 1
        z_data_list.append(z_data)

        #Finds if the timestap matches a snapshot an appends index
        if any(t== t_snap for t_snap in snaps):
            snapshot_index_list.append(int(c-1))
    filename = filename.split('.t')[0]+'_imag'

else: #if the we want to plot p
    A = pa.cx_cube() #Create pa.mat object (just as arma::mat in C++)
    A.load("./datafiles/"+str(filename)) #Load the content of the matrix you saved into your Python program.
    # A function for a Gaussian that is travelling
    # in the x direction and broadening as time passes


    # Array of time points (Dynamically allocates T)
    dt = 0.000025
    t_points = np.arange(0, dt*(np.shape(A)[0]), dt)

    for t in t_points:
        z_data = np.rot90(np.array(A[pa.single_slice, c]))
        c += 1

        out = np.real(np.multiply(np.conj(z_data), z_data))
        z_data_list.append(out)
        
            #Finds if the timestap matches a snapshot an appends index
        if any(t== t_snap for t_snap in snaps):
            snapshot_index_list.append(int(c-1))
    filename = filename.split('.t')[0]

#Fetches the different potentials in case overlay is requested by input arguments
V_0 = pa.mat()
V_1 = pa.mat()
V_2 = pa.mat()
V_3 = pa.mat()
V_0.load("./datafiles/box.dat")
V_1.load("./datafiles/box_single_slit.dat")
V_2.load("./datafiles/box_double_slit.dat")
V_3.load("./datafiles/box_triple_slit.dat")
V_0 = np.rot90(np.array(V_0))
V_1 = np.rot90(np.array(V_1))
V_2 = np.rot90(np.array(V_2))
V_3 = np.rot90(np.array(V_3))

V = [V_0,V_1,V_2,V_3]




#
# Now the list z_data_list contains a series of "frames" of z(x,y,t),
# where each frame can be plotted as a 2D image using imshow. Let's
# animate it!
#

# Some settings
fontsize = 12
t_min = t_points[0]
x_min, x_max = x_points[0], x_points[-1]
y_min, y_max = y_points[0], y_points[-1]

# Create figure
fig = plt.figure()
ax = plt.gca()




#-----Setting up plot format------
norm = matplotlib.cm.colors.Normalize(vmin=np.min(z_data_list[0]), vmax=np.max(z_data_list[0]))


V[slits_overlay] = np.ma.masked_where(V[slits_overlay] > 10.00, V[slits_overlay])
pb = ax.imshow(V[slits_overlay],extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("Greys_r"), norm=norm)
img = ax.imshow(z_data_list[0],extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)

# Axis labels
plt.xlabel("x [Units of distance /1]", fontsize=fontsize)
plt.ylabel("y [Units of distance /1]", fontsize=fontsize)
plt.xticks(fontsize=fontsize)
plt.yticks(fontsize=fontsize)

# Add a colourbar
cbar = fig.colorbar(img, ax=ax)
cbar.set_label(z_axis_label, fontsize=fontsize)
cbar.ax.tick_params(labelsize=fontsize)

# Add a text element showing the time
time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white",
                    horizontalalignment="right", verticalalignment="top", fontsize=fontsize)
#-----Setting up plot format(end)------



#img = ax.imshow(z_data_list[0], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)

# Function that takes care of updating the z data and other things for each frame
def animation(i):
    # Normalize the colour scale to the current frame
    norm = matplotlib.cm.colors.Normalize(vmin=np.min(z_data_list[i]), vmax=np.max(z_data_list[i]))
    img.set_norm(norm)

    # Update z data(Sums with V only to make the ovarlay for the boundary conditions)
    img.set_data(np.add(z_data_list[i], V[slits_overlay]))

    # Update the time label
    current_time = t_min + i * dt
    time_txt.set_text("t = {:.3e}".format(current_time))

    return img


# Use matplotlib.animation.FuncAnimation to put it all together
anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(z_data_list), 2), repeat=False, blit=0)


# Save the animation
anim.save('./figures/'+filename+'_animation.gif', writer="ffmpeg", fps=15)

#Saves snapshots if arguments say so
for ind in snapshot_index_list:
    norm = matplotlib.cm.colors.Normalize(vmin=np.min(z_data_list[ind]), vmax=np.max(z_data_list[ind]))
    V[slits_overlay] = np.ma.masked_where(V[slits_overlay] > 10.00, V[slits_overlay])
    pb = ax.imshow(V[slits_overlay],extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("Greys_r"), norm=norm)
    img = ax.imshow(z_data_list[ind],extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)
    img.set_norm(norm)
    img.set_data(np.add(z_data_list[ind], V[slits_overlay]))
    time_txt.set_text("t = {:.3e}".format(t_points[ind]))
    plt.savefig('./figures/'+filename+'_time'+str(t_points[ind])+'.pdf')
