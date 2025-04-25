"""
PLATE CAM WITH RECIPROCATING FLAT-FACE FOLLOWER

The reciprocating flat-face follower of a plate
cam is to rise 2 in with simple harmonic motion
in 180◦ of cam rotation and return with simple
harmonic motion in the remaining 180◦. The prime-circle
radius is 1.5 in, and the cam rotates counterclockwise.
Construct the displacement diagram and the cam
profile, offsetting the follower stem by 0.75 in in
the direction that reduces the bending of the follower
during rise.
"""
#%% Libraries
from DiskCamMechanismLibrary import PDCamFlatFaceFollower
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
Rbase=1.5 #prime radius circle
Rbroca=3/16 # drill bit radius (cam center)
eccentricity = 0.75
# Angular position of the follower in radians
FollowerAng = np.pi/2
# angular sweep from zero to 2 pi radians
theta = np.linspace(0,1,500)*2*np.pi
# calculate displacement, velocity, acceleration
y,yp,ypp = SimpleHarmonicMotion(theta,L)

# Group data in dictionary, for other parameters consult the documentation of DiskCamMechanismLibrary
CamData={'theta':theta,
         'y':y,
         'yp':yp,
         'ypp':ypp,
         'Rbase':Rbase,
         'Rhole':Rbroca,
         'epsilon':eccentricity,
         'FollowerAng':FollowerAng,
         'Followerwidth': 4/16,
         'turn_direction':'anti-clockwise',
        }

#%% Calculating the Cam Profile
Cam=PDCamFlatFaceFollower(**CamData)

#%% Motion diagram
figMD=plt.figure()
Cam.PlotMotionDiagram(figMD)

#%% Plot the cam profile
figPCam=plt.figure()
Cam.PlotCamFlatFollower(figPCam)

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
anim3.save('mp4/Cam602.mp4',writer=writer,dpi=dpi)