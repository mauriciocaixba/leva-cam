"""
The reciprocating radial roller follower of a plate
cam is to rise 2 in with simple harmonic motion
in 180◦ of cam rotation and return with simple
harmonic motion in the remaining 180◦. If the roller
radius is 0.375 in and the prime-circle radius is 2 in,
construct the displacement diagram, the pitch curve,
and the cam profile for clockwise cam rotation.
"""
#%% Libraries
from DiskCamMechanismLibrary import PDCamRollerFollower
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

# %% Simple Harmonic Motion
def SimpleHarmonicMotion(th,L):
    y = 0.5*L*(1-np.cos(th))
    yp = 0.5*L*np.sin(th)
    ypp = 0.5*L*np.cos(th)
    return y,yp,ypp
    
# %% problem data
L=2
Rprime=2 #prime radius circle
rd=0.375 #roller radius
Rdrill=3/16 # drill bit radius (cam center)
eccentricity = 0.0
FollowerAng = np.pi/2 # Angular position of the follower in radians
theta = np.linspace(0,1,500)*2*np.pi # angular sweep from zero to 2 pi radians
# calculate displacement, velocity, acceleration
y,yp,ypp = SimpleHarmonicMotion(theta,L) 

# Group data in dictionary, for other parameters consult the documentation of DiskCamMechanismLibrary
CamData={'theta':theta,
         'y':y,
         'yp':yp,
         'ypp':ypp,
         'Rbase':Rprime,
         'Rhole':Rdrill,
         'epsilon':eccentricity,
         'FollowerAng':FollowerAng,
         'Followerwidth': 4/16,
         'turn_direction':'clockwise',
         'Rroller':rd
        }

#%% Calculating the Cam Profile
Cam=PDCamRollerFollower(**CamData)

#%% Motion diagram
figMD=plt.figure()
Cam.PlotMotionDiagram(figMD)

#%% Plot the cam profile
figPCam=plt.figure()
Cam.PlotCamRollerFollower(figPCam)

#%% Cam animation
fig, ax=plt.subplots()
ax.set_axis_off()
init_func=Cam.initAnim(ax),
dpi=100
width = 1920/dpi
hight = 1080/dpi
fig.set_size_inches(width,hight)

anim3 = FuncAnimation(fig, Cam, frames=np.arange(1000),
                    interval=100, blit=False)
plt.show()

#%% Saving the cam animation to a file
writer = animation.writers['ffmpeg'](fps=30)
anim3.save('Cam01.mp4',writer=writer,dpi=dpi)
