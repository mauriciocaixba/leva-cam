"""
El seguidor de movimiento alternativo, radial y de rodillo,
de una leva de placa debe subir 2 pulg con movimiento armónico
simple en 180° de rotación de la leva, y retornar con movimiento
armónico simple en los 180° restantes. Si el radio del rodillo es de
0.375 pulg y el del círculo primario es de 2 pulg constrúyase
el diagrama de desplazamientos, la curva de paso y el perfil de la
leva para una rotación de ésta en el mismo sentido que el movimiento
de las manecillas del reloj.

"""
#%% Librerias
from DiskCamMechanismLibrary import PDCamRollerFollower
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation

# %% Movimiento armonico simple
def MovArmonicoSimple(th,L):
    y = 0.5*L*(1-np.cos(th))
    yp = 0.5*L*np.sin(th)
    ypp = 0.5*L*np.cos(th)
    return y,yp,ypp
    
# %% Datos del problema
L=2
Rbase=2 #radio primario
rd=0.375 #radio rodillo
Rbroca=3/16 # Radio de la broca (centro de la leva)
excentricidad = 0.0
posAngularSeguidor = np.pi/2 # posicion angular del seguidor en radianes
theta = np.linspace(0,1,500)*2*np.pi # barrido angular de cero a 2pi radianes
# calcular desplazamiento, velocidad, aceleracion
y,yp,ypp = MovArmonicoSimple(theta,L) 

# Agrupar datos en diccionario, para otros parametros consultar la documentacion de DiskCamMechanismLibrary
CamData={'theta':theta,
         'y':y,
         'yp':yp,
         'ypp':ypp,
         'Rbase':Rbase,
         'Rhole':Rbroca,
         'epsilon':excentricidad,
         'FollowerAng':posAngularSeguidor,
         'Followerwidth': 4/16,
         'turn_direction':'clockwise',
         'Rroller':rd
        }

#%% Calcular el perfil de la Leva
Leva=PDCamRollerFollower(**CamData)

#%% Diagrama de movimiento
figMD=plt.figure()
Leva.PlotMotionDiagram(figMD)

#%% Graficar el perfil de la leva
figPCam=plt.figure()
Leva.PlotCamRollerFollower(figPCam)

#%% Animación de la leva
fig, ax=plt.subplots()
ax.set_axis_off()
init_func=Leva.initAnim(ax),
dpi=100
width = 1920/dpi
hight = 1080/dpi
fig.set_size_inches(width,hight)

anim3 = FuncAnimation(fig, Leva, frames=np.arange(1000),
                    interval=100, blit=False)
plt.show()

#%% Guardar la animación de la leva en un archivo
writer = animation.writers['ffmpeg'](fps=30)
anim3.save('Leva01.mp4',writer=writer,dpi=dpi)
