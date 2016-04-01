import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as anm


#plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'
plt.close('all')
data = np.loadtxt('solar_system.dat')
data2 = data[:,15:30]

fig = plt.figure()
ax = p3.Axes3D(fig)
ax.set_xlim3d([np.min(data2[:,0::3]), np.max(data2[:,0::3])])
ax.set_xlabel('X')
ax.set_ylim3d([np.min(data2[:,1::3]), np.max(data2[:,1::3])])
ax.set_ylabel('Y')
ax.set_zlim3d([np.min(data2[:,2::3]), np.max(data2[:,2::3])])
ax.set_zlabel('Z')

# choose a different color for each trajectory
colors = plt.cm.jet(np.linspace(0, 1, np.size(data2[0,:])/3))

# set up lines and points
lines = sum([ax.plot([], [], [], '-', c=c)
             for c in colors], [])
pts = sum([ax.plot([], [], [], 'o', c=c)
           for c in colors], [])
               
ax.view_init(30, 0)

data3 = np.reshape(data2,(np.size(data2[0,:])/3,np.size(data2[:,0]),3))
n = 0
for i in np.arange(0,int(np.size(data2[0,:])/3),1):
    data3[i,:,0:3] = data2[:,i+n:i+n+3]
    n = n + 2

def init():
    for line, pt in zip(lines, pts):
        line.set_data([], [])
        line.set_3d_properties([])

        pt.set_data([], [])
        pt.set_3d_properties([])
    return pts + lines,



def animate(i):
    # we'll step two time-steps per frame.  This leads to nice results.
    i = (1 * i) % data3.shape[1]

    for line, pt, xi in zip(lines, pts, data3):
        x, y, z = xi[:i,0:3].T
        line.set_data(x, y)
        line.set_3d_properties(z)

        pt.set_data(x[-1:], y[-1:])
        pt.set_3d_properties(z[-1:])

    #ax.view_init(30, 0.3 * i)
    fig.canvas.draw()
    return pts + lines
    
anim = anm.FuncAnimation(fig, animate, init_func=init,
                              frames=int(np.size(data2[:,0])), interval=1, blit=True)

#writer = anm.writers['ffmpeg'](fps=30)
#anim.save('inner_sol_sys.mp4', fps=15, extra_args=['-vcodec', 'libx264'])#, 'ffmpeg_file', fps=15, extra_args=['-vcodec', 'libx264']
