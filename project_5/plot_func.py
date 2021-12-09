import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pyarma as pa

#
# Now the list z_data_list contains a series of "frames" of z(x,y,t), 
# where each frame can be plotted as a 2D image using imshow. Let's
# animate it!
#

def animate(z_data_list, name):
    # Create figure
    fig = plt.figure()
    ax = plt.gca()

    # Create a colour scale normalization according to the max z value in the first frame
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[0]))

    # Plot the first frame
    img = ax.imshow(z_data_list[0], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)

    # Axis labels
    plt.xlabel("x", fontsize=fontsize)
    plt.ylabel("y", fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.yticks(fontsize=fontsize)

    # Add a colourbar
    cbar = fig.colorbar(img, ax=ax)
    cbar.set_label("z(x,y,t)", fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)

    # Add a text element showing the time
    time_txt = plt.text(0.95, 0.95, "t = {:.3e}".format(t_min), color="white", 
                        horizontalalignment="right", verticalalignment="top", fontsize=fontsize)

    # Function that takes care of updating the z data and other things for each frame
    def animation(i):
        # Normalize the colour scale to the current frame?
        norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[i]))
        img.set_norm(norm)

        # Update z data
        img.set_data(z_data_list[i])

        # Update the time label
        current_time = t_min + i * dt
        time_txt.set_text("t = {:.3e}".format(current_time))

        return img

    # Use matplotlib.animation.FuncAnimation to put it all together
    anim = FuncAnimation(fig, animation, interval=1, frames=np.arange(0, len(z_data_list), 2), repeat=False, blit=0)

    # Run the animation!
    plt.show()

    # # Save the animation
    anim.save('./animation_'+name+'.mp4', writer="ffmpeg", bitrate=-1, fps=30)

# Some settings
print("1. Compare probability deviation")
print("2. Create double slit animation and figure")
choice = int(input("Input your choice: "))
if choice == 1:
    filenames = ['double_slit_7', 'double_slit_7_2']
    label_names = ['No potential', 'Two slits']
    t = np.linspace(0,0.008, int(0.008/(2.5e-5))+1)
    fig, ax = plt.subplots()
    
    for i, j in zip(filenames, label_names):
        d = pa.cx_cube()
        d.load(i)
        d = np.array(d)
        d = np.swapaxes(d,0,1)
        z_d_l = (d*np.conjugate(d)).real
        p = 1-np.sum(z_d_l, axis=(1,2))
        ax.plot(t,abs(p), label = j)
    ax.set_ylabel('Absolute value of probability deviation')
    ax.set_xlabel('Dimensionless time [1]')
    ax.set_yscale('log')
    ax.legend(), ax.grid(), plt.savefig('Time.pdf'), plt.clf()

    for filename in filenames: 
        d = pa.cx_cube()
        d.load(filename)
        d = np.array(d)
        d = np.swapaxes(d,0,1)
        z_d_l = (d*np.conjugate(d)).real
        animate(z_d_l, name = filename)

elif choice == 2:
    filename = "double_slit_8" 
    data = pa.cx_cube()
    data.load(filename)
    data = np.array(data)
    # M, N, M
    t = np.linspace(0, 0.008, int(0.008/2.5e-5) +1)

    data = np.swapaxes(data,0,1)
    z_data_list = (data*np.conjugate(data)).real
    p = 1 - np.sum(z_data_list, axis=(1,2))
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.plot(t, abs(p))
    ax.set_ylabel('Absolute value of probability deviation')
    ax.set_xlabel('Dimensionless time [1]')
    ax.grid(); plt.show()
    #print(np.shape(z_data_list))

    fontsize = 12
    t_min = 0
    x_min, x_max = 0, 1
    y_min, y_max = 0, 1
    dt = 2.5e-5


    i_list = [0, 50, -1]
    fig, ax = plt.subplots(ncols=len(i_list), figsize=(10,5))
    k = 0

    for i in i_list:
        norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.max(z_data_list[i]))
        img = ax[k].imshow(z_data_list[i], extent=[x_min,x_max,y_min,y_max], cmap=plt.get_cmap("viridis"), norm=norm)
        ax[k].set_xlabel("x", fontsize=fontsize)
        ax[k].set_ylabel("y", fontsize=fontsize)
        #ax[k].set_xticks(fontsize=fontsize)
        #ax[k].set_yticks(fontsize=fontsize)
        cbar = fig.colorbar(img, ax=ax[k],fraction=0.046, pad=0.04)
        cbar.set_label("z(x,y,t)", fontsize=fontsize)
        cbar.ax.tick_params(labelsize=fontsize)
        time_txt = ax[k].text(0.95, 0.95, "t = {:.3e}".format(t[i]), color="white", 
                        horizontalalignment="right", verticalalignment="top", fontsize=fontsize)
        k += 1
    plt.tight_layout()
    plt.savefig("snapshots_"+filename)
