#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DiskCamMechanismLibrary
=====

This is a collection of python classes to produce the planar disk cam profile
from a given set of kinematic curves (motion diagrams). It is intended to serve as an educational tool for undergraduate students.


classes
--------------------
    
1. PDCamKnifeFollower(**DataCam)
    Planar Disk Cam with Translating Knife-Edge Translating Follower (Centered and eccentric)
    
2. PDCamRollerFollower(**DataCam)
    Planar Disk Cam with Translating Roller Follower (Centered and eccentric)
    
3. PDCamFlatFaceFollower(**DataCam)
    Planar Disk Cam with Translating Flat-faced Follower (Centered and eccentric)
    
4. PDCamOscillatingRollerFollower(**DataCam)
    Planar Disk Cam with Roller Oscilatory Follower
    

For further explanation, please see the corresponding documentation.

Note:
--------
NumPy and matplotlib have to be installed already on a python 3.7 system.


@author: Mauricio Caixba-Sanchez
"""

class PDCamKnifeFollower:
    """ 
    PDCamKnifeFollower(**CamData)

    
    Description
    ------------
    
    This is a design tool for Planar Disk Cam with Knife-Edge Follower
        (Centered an eccentric)
        
    The main functionalities are:
        1. It gives a plot with the Motion Diagram.
        2. It gives the cam's profile data in rectangular coordinates.
        3. The former data is ploted in a matplotlib Figure.
        4. Also there is a set of methods which help to animate cam's rotation.
    
    Parameters
    -----------
    
    The cam linkage parameters have to be provided via a python dictionary.
    
    theta : data-type
        1D NumPy array with N equal-spaced values between (0,2pi)
        This variable is the x-axis of the displacement diagram.
                
    y : data-type
        (1D NumPy array) Knife-Edge Follower displacement as function of theta.
        
    yp : data-type
        (1D NumPy array) Knife-Edge Follower velocity as function of theta.
        
    ypp : data-type
        (1D NumPy array) Knife-Edge Follower aceleration as function of theta.
        
    Rbase : float
        Base Circle Radius.      

       
        Note:                
            y, yp and ypp have to be the same length.
       

    Optional Parameters:
    --------------------
        
    epsilon : float
        Knife Follower eccentricity (positive for right position and negative otherwise)
        
    FollowerAng : float
        Knife Follower's axis angle (radians)
        
    Followerwidth : float
        Knife-Edge Follower width
        
    Followerhight : float
        Knife-Edge Follower hight
        
    Rhole : float
        Hole drill radius. For make-up purpose.
        
    turn_direction : str
        Direction of rotation ('clockwise' or 'anti-clockwise')

    CamProfileColor : str or color array
        Color of the Cam Profile according to matplotlib's color format
        
    KnifeFollowerColor : str or color array
        Color of the Knife-Edge according to matplotlib's color format
        
    BlockGuideColor : str or color array
        Color of the Block Guide according to matplotlib's color format
        
    CenterCamBgColor : str or color array
        Bakground color of the center Cam according to matplotlib's color format
        
    CenterCamFgColor : str or color array
        Foreground color of the center Cam according to matplotlib's color format
        
    
    Attributes
    -----------
    Xp : data-type
        1D NumPy array which contents a sequence of x-rectangular coordinates
        that describes the cam profile.
        
    Yp : data-type
        1D NumPy array which contents a sequence of y-rectangular coordinates
        that describes the cam profile.
        
    xFollower : data-type
        1D NumPy array which contents a sequence of x-rectangular coordinates
        that describes the Follower geometry at the initial condition.
        
    yFollower : data-type
        1D NumPy array which contents a sequence of y-rectangular coordinates
        that describes the Follower geometry at the initial condition.
        
    Xblock1: data-type
        1D NumPy array which contents a sequence of x-rectangular coordinates
        that partially describes the Follower's Guide.
    
    Yblock1 : data-type
        1D NumPy array which contents a sequence of y-rectangular coordinates
        that partially describes the Follower's Guide.
    
    Xblock2 : data-type
        1D NumPy array which contents a sequence of x-rectangular coordinates
        that partially describes the Follower's Guide.
    
    Yblock2 : data-type
        1D NumPy array which contents a sequence of y-rectangular coordinates
        that partially describes the Follower's Guide.
    
    
    Methods
    --------
                    
    PlotMotionDiagram() 
        It plots the position, velocity an acceleration of alternating follower
        as function of the angular displacement.
                         
    CamForKnifeFollower() 
        It calculates the rectangular coordinates of the cam profile according
        to the given data.
   
    """
    from numpy import sqrt
    
    def __init__(self,**kwargs):
        """
        numpy and matplotlib have to be installed on the system
        """
        from numpy import pi
        """
        Default arguments (optional arguments)
        """
        self.epsilon=0
        self.FollowerAng=pi/2
        self.Followerwidth=None
        self.Followerhight=None
        self.Rhole=None
        self.turn_direction='clockwise'
        self.CamProfileColor='goldenrod'
        self.KnifeFollowerColor='brown'
        self.BlockGuideColor='darkslategray'
        self.CenterCamBgColor=[0.1,0.1,0.4]
        self.CenterCamFgColor='yellowgreen'        
        """
        Updating input arguments
        """
        self.__dict__.update(**kwargs)         
        """
        Calculating the Cam Profile and initial data for Knife-Edge Follower
        """
        self.setparam()        
        self.CamForKnifeFollower()
        self.xFollower,self.yFollower=self.KnifeFollower(self.Xp[0],self.Yp[0])               
        self.Xblock1, self.Yblock1, self.Xblock2, self.Yblock2=self.FollowerGuide()
        
        
    def setparam(self):
        from numpy import sqrt
        if self.Followerwidth==None:
            self.Followerwidth=max(self.y)/8
        if self.Followerhight==None:
            self.Followerhight=3*max(self.y)
        if self.Rhole==None:
            self.Rhole=0.15*self.Rbase
        self.Lmin=0.5*max(self.y)   
        self.dmax=sqrt(self.Rbase**2-self.epsilon**2)+1.1*self.Followerwidth+self.Lmin+max(self.y)
        self.dmin=sqrt(self.Rbase**2-self.epsilon**2)+1.1*self.Followerwidth
        
    def PlotMotionDiagram(self,fig):
        """PlotMotionDiagram(fig)
        
        Description:
        ------------
        It plots the Motion diagram associated with the Cam instance
        
        Parameters:
        -----------
            
            fig : FigureClass
                matplotlib figure instance
                
                How to get it example
                
                >>> from matplotlib.pyplot import figure
                
                >>> fig=figure()
                
        Return:
        --------
        
            axes : matplotlib.axes._subplots.AxesSubplot object array
                matplotlib axes array where the diagrams where plotted.
                More elements can be added by the user to such axes.
        """        
        axes=fig.subplots(3,1, sharex=True)    
        axes[0].plot(self.theta,self.y)
        axes[0].set_title('Displacement Diagram')
        axes[0].grid(True)
        axes[1].plot(self.theta,self.yp)
        axes[1].set_title('Velocity Diagram')
        axes[1].grid(True)
        axes[2].plot(self.theta,self.ypp)
        axes[2].set_title('Acceleration Diagram')
        axes[2].set(xlabel='Angle [radians]')
        axes[2].grid(True)              
        
        return axes  
        
    def CamForKnifeFollower(self):
        """CamForKnifeFollower()
        
        Description
        ------------
        
        It computes the rectangular coordinates of the Cam Profile for
        a Knife-Edge Follower.
        
        Returns
        --------
        
        Xp : data-type
            1D NumPy array which contents a sequence of x-rectangular coordinates
            that describes the cam profile.
        
        Yp : data-type
            1D NumPy array which contents a sequence of y-rectangular coordinates
            that describes the cam profile.        
        """
        #Note: This method is called within the __init__method and 
        #its values are stored in the Xp and Yp as attributes.
        from numpy import exp
        epsilon=self.epsilon
        direction=self.turn_direction
        Rbase=self.Rbase
        theta,y=self.theta,self.y
        beta=self.FollowerAng        
        Q=(((Rbase**2-epsilon**2)**0.5+y)-1j*epsilon)
        if direction=='clockwise':            
            R=Q*exp(1j*(beta+theta))
        elif direction=='anti-clockwise':
            R=Q*exp(1j*(beta-theta))
        else:
            print('Please specify a valid direction of rotation')
        Xp=R.real
        Yp=R.imag
        self.Xp,self.Yp=Xp,Yp
        self.Q=Q # Complex distance From de Cam center to the Knife follower tip
        return Xp,Yp
    
    def PressureAngle(self):
        """
        PressureAngle()
        
        Description
        ------------
        It computes the pressure angle between the follower axis and the cam profile
        
        Return
        ______
        
        phi : data-type
            1D NumPy array which contents the pressure angle in radians for each
            cam's angular displacement.        
        """
        from numpy import absolute, arccos        
        Q=self.Q
        cosphi=Q.real/absolute(Q)
        phi=arccos(cosphi)        
        return phi
    
    def PlotPressureAngle(self,Axis):
        """
        PlotPressureAngle(Axis)
        
        Description
        ------------
        
        It plots the cam's pressure angle curve.
        
        Parameters
        -----------
        
        Axis : matplotlib.axes._subplots.AxesSubplot object
            matplotlib axis where the pressure angle curve will be plotted.

        Returns:
            
            line : matplotlib.lines.Line2D
                Gives a Line2D instance with x and y data en sequences xdata, ydata.        
        """        
        from numpy import pi
        phi=self.PressureAngle()*180/pi
        line,=Axis.plot(self.theta,phi,color='blue')
        Axis.set_xlabel(r'$\theta$ [rad]')
        Axis.set_ylabel(r'Angle Pressure [degrees]')
        Axis.set_title('Angular displacement vs Angle Pressure')
        return line
        
    def KnifeFollower(self,xtip,ytip,angle=None):
        """
        KnifeFollower(xtip,ytip,angle=None)
        
        Description
        ------------
        It gives the rectangular coordinates to render a Matplotlib PolyFill object
        to represent de oscillating knife follower shape. This is an internal method 
        for use of other internal methods.
        
        Parameters
        ----------
        xtip : float
            Tip Follower's x-coordinate position
            
        ytip : float
            Tip Follower's y-coordinate position
            
        angle : float
            Angular position of the Follower, default value (None) have been set to the
            original Follower axis's angle (Default or given by the user).
            
        Returns
        --------
        (data-type,data-type) 
            Tuple with two 1D NumPy array which contents the (x,y)-coordinates of the follower shape at
            a given position.
        """
        from numpy import  exp, array    
        wth=self.Followerwidth
        ht=self.Followerhight
        if angle==None:
            angle=self.FollowerAng        
        xFollower0=array([0,0.866*wth,0.866*wth,1.1*wth,1.1*wth,ht,ht,1.1*wth,1.1*wth,0.866*wth,0.866*wth,0])
        yFollower0=array([0,0.5*wth,wth,wth,0.5*wth,0.5*wth,-0.5*wth,-0.5*wth,-wth,-wth,-0.5*wth,0])
        Follower=(xFollower0+1j*yFollower0)*exp(1j*angle)+(xtip+1j*ytip)
        return Follower.real, Follower.imag
    
    def FollowerGuide(self,angle=None):
        """
        FollowerGuide(angle=None)
        
        Description
        ------------
        It gives the rectangular coordinates to render a Matplotlib PolyFill object
        to represent a simple guide shape. This is an internal method for use 
        of other internal methods.
        
        Parameter
        ----------
        
        angle : float
            Angular position of the axis Follower, default value (None) have been set to the
            original Follower axis's angle (Default or given by the user).
        
        Returns
        ---------------------------
        Xblock1: data-type
        1D NumPy array which contents a sequence of x-rectangular coordinates
        that partially describes the Follower's Guide.
    
        Yblock1 : data-type
            1D NumPy array which contents a sequence of y-rectangular coordinates
            that partially describes the Follower's Guide.
    
        Xblock2 : data-type
            1D NumPy array which contents a sequence of x-rectangular coordinates
            that partially describes the Follower's Guide.
    
        Yblock2 : data-type
            1D NumPy array which contents a sequence of y-rectangular coordinates
            that partially describes the Follower's Guide.
        
        """
        from numpy import exp, array
        if self.Followerwidth==None:
            self.Followerwidth=max(self.y)/8            
        wth=self.Followerwidth
        dz=self.dmax
        if angle==None:
            angle=self.FollowerAng        
        xb=array([dz,1.4*dz,1.4*dz,dz])
        yb=array([0.5*wth,0.5*wth,2.5*wth,2.5*wth])
        Zblock1=(xb+1j*yb-1j*self.epsilon)*exp(1j*angle)
        Zblock2=(xb-1j*yb-1j*self.epsilon)*exp(1j*angle)
        Xblock1=Zblock1.real
        Yblock1=Zblock1.imag
        Xblock2=Zblock2.real
        Yblock2=Zblock2.imag
        return Xblock1, Yblock1, Xblock2, Yblock2
    
    def Linkagecoord(self,width,x1,y1,x2,y2,angle=None):
        """
        Internal method for computing spring's geometry parts
        """
        from numpy import pi, piecewise, linspace, cos, sin, absolute, exp
        eps=self.epsilon
        if angle==None:
            gamm=self.FollowerAng
        else:
            gamm=angle
        q0=x1+1j*y1
        q1=x2+1j*y2
        height=absolute(q1-q0)
        betas=linspace(0,2*pi,50)
        Xreal=0.5*width*cos(betas)
        Ximag=piecewise(betas,[betas<pi,betas>=pi],[lambda x: 
            (0.5*width*sin(x)+height),
                    lambda x: 0.5*width*sin(x)])
        X=Xreal+1j*Ximag
        W=1j*(q1-q0)*X/absolute(q1-q0)+q1
        Linkage=(W-1j*eps)*exp(1j*gamm)        
        return Linkage
    
    def SpringBlock(self, dmax, dmin, Lmin,N):
        """
        Internal method for computing spring geometry
        """
        from numpy import linspace, ones_like, arange, concatenate
        wth=self.Followerwidth
        r=0.85*Lmin/(2*N)
        L=dmax-dmin
        p=(L-2*r)/(N-0.5)
        x1=linspace(dmin+r,dmax-0.5*p-r,N)
        y1=(0.75*wth)*ones_like(x1)
        x2=linspace(dmin+r+0.5*p,dmax-r,N)
        y2=-(0.75*wth)*ones_like(x2)

        for kk in arange(N+1):
            
            if kk==0:
                linkfw=self.Linkagecoord(2*r,x1[kk],y1[kk],x2[kk],y2[kk])
                linkfw=linkfw.reshape(len(linkfw),1)
                linksfwchain=linkfw
            elif kk<N:
                linkfw=self.Linkagecoord(2*r,x1[kk],y1[kk],x2[kk],y2[kk])
                linkfw=linkfw.reshape(len(linkfw),1)
                linksfwchain=concatenate((linksfwchain,linkfw),axis=1)
                
            if kk==0:
                linkbk=self.Linkagecoord(2*r,x1[kk],y1[kk],x1[kk],y2[kk])
                linkbk=linkbk.reshape(len(linkbk),1)
                linksbkchain=linkbk
            elif kk==N:
                linkbk=self.Linkagecoord(2*r,x2[kk-1],y1[kk-1],x2[kk-1],y2[kk-1])
                linkbk=linkbk.reshape(len(linkbk),1)
                linksbkchain=concatenate((linksbkchain,linkbk),axis=1)
            else:
                linkbk=self.Linkagecoord(2*r,x1[kk],y1[kk],x2[kk-1],y2[kk-1])
                linkbk=linkbk.reshape(len(linkbk),1)
                linksbkchain=concatenate((linksbkchain,linkbk),axis=1)
            
        return linksfwchain, linksbkchain
        
        
    def PlotCamKnifeFollower(self,fig,detailed=True):
        """ 
        PlotCamKnifeFollower(fig,detailed=True)
        
        Description:
        ------------
        This plots the Cam Profile with the Knife-Edge Follower and
        a decorative spring.
        
        Parameters:
        -----------
            
            fig : FigureClass
                matplotlib figure instance
                
                How to get it example
                
                >>> from matplotlib.pyplot import figure
                
                >>> fig=figure()
                
            detailed : bool
                If True, all the decorative elements are depicted, otherwise only
                the cam profile is rendered with a center mark.
                
        Return:
        --------
        
            axis : matplotlib.axes._subplots.AxesSubplot object
                matplotlib axis where the forms where plotted.
                More elements can be added by the user to such axis.        
       """
        from numpy import cos,sin, arange
        Rbase,Rh,theta=self.Rbase,self.Rhole,self.theta
        xp,yp=self.Xp, self.Yp
        xf,yf=self.xFollower, self.yFollower
        xb1,yb1=self.Xblock1, self.Yblock1
        xb2,yb2=self.Xblock2, self.Yblock2
                        
        dmin=self.dmin+self.y[0]
        Lmin=self.Lmin
        dmax=self.dmax   
        linksfwchain, linksbkchain=self.SpringBlock(dmax, dmin, Lmin,12)        

        axis=fig.add_subplot(111)
        axis.fill(xp,yp,color=self.CamProfileColor) 
        axis.set_aspect('equal')
        if detailed==True:
            axis.plot(Rbase*cos(theta),Rbase*sin(theta),'w--')
            axis.fill(Rh*cos(theta),Rh*sin(theta),'w',linewidth=3)
            for column in arange(linksbkchain.shape[1]):
                axis.fill(linksbkchain[:,column].real,linksbkchain[:,column].imag,'k')            
            axis.fill(xf,yf,color=self.KnifeFollowerColor)
            for column in arange(linksfwchain.shape[1]):
                axis.fill(linksfwchain[:,column].real,linksfwchain[:,column].imag,'gray')          
            axis.plot(0,0,'k+',linewidth=3)
            axis.fill(xb1,yb1, color=self.BlockGuideColor)
            axis.fill(xb2,yb2, color=self.BlockGuideColor)
        else:
            axis.plot(0,0,'w+',linewidth=3)
            
        return axis
            
        
    def initAnim(self, ax):
        """
        initAnim(ax)
        
        Description
        ------------
        Method to initialize the Animation Parameters.
        
        Parameters:
        --------
        
            ax : matplotlib.axes._subplots.AxesSubplot object
                matplotlib axis where the initial geometries will be plotted.        
        """
        from numpy import sin, cos, sqrt, arange, linspace, concatenate, exp, pi
        
        dmin=self.dmin+self.y[0]
        Lmin=self.Lmin
        dmax=self.dmax    
        linksfgchain, linksbkchain=self.SpringBlock(dmax, dmin, Lmin,12)
        self.axAnimation=ax
        colorblk=self.BlockGuideColor
        Xblock1=self.Xblock1
        Yblock1=self.Yblock1
        Xblock2=self.Xblock2
        Yblock2=self.Yblock2
        # Render Follower's Guide Block
        self.axAnimation.fill(Xblock1,Yblock1,color=colorblk)
        self.axAnimation.fill(Xblock2,Yblock2,color=colorblk)        
        # Render Cam Profile
        self.poly1,=self.axAnimation.fill(self.Xp,self.Yp,color=self.CamProfileColor)
        # Render spring's rear view
        self.springbk=[]
        for col in arange(linksbkchain.shape[1]):
                l1,=self.axAnimation.fill(linksbkchain[:,col].real,linksbkchain[:,col].imag,'k')
                self.springbk.append(l1)
        # Render Knife-Edge Follower
        self.poly2,=self.axAnimation.fill(self.xFollower,self.yFollower,
                                          color=self.KnifeFollowerColor)
        # Render spring's front view
        self.springfg=[]
        for col in arange(linksfgchain.shape[1]):
                l1,=self.axAnimation.fill(linksfgchain[:,col].real,linksfgchain[:,col].imag,'gray')
                self.springfg.append(l1)
        self.axAnimation.fill(self.Rhole*cos(self.theta),self.Rhole*sin(self.theta),
                              color=self.CenterCamBgColor)
        
        # Calculating and render the center mark on the cam profile
        alp1=linspace(45,135,10)*pi/180
        alp2=linspace(-45,-135,10)*pi/180
        self.alp=concatenate((alp1,alp2))
        arc=0.8*self.Rhole*exp(1j*self.alp)        
        self.artc,=self.axAnimation.fill(arc.real,arc.imag, color=self.CenterCamFgColor)
        
        # Setting the Axes limits according to the dimensions of the main parts of the animation
        LimCam=max(sqrt((self.Xp**2+self.Yp**2)))
        LimBlk1x=max(self.Xblock1)
        LimBlk2x=max(self.Xblock2)
        LimBlk1y=max(self.Yblock1)
        LimBlk2y=max(self.Yblock2)
        Limxp=1.05*max([LimCam,LimBlk1x,LimBlk2x])
        Limxn=1.05*min([-LimCam,LimBlk1x,LimBlk2x])
        Limyp=1.05*max([LimCam,LimBlk1y,LimBlk2y])
        Limyn=1.05*min([-LimCam,LimBlk1y,LimBlk2y])        
        self.axAnimation.set_xlim(Limxn,Limxp)
        self.axAnimation.set_ylim(Limyn,Limyp)
        self.axAnimation.set_aspect('equal')
        return self.poly1, self.poly2
        
    def __call__(self,k):
        """Method to update Animation data
        """
        from numpy import mod, exp, transpose, array, arange
        
        i=5*k
        ii=mod(i,len(self.theta))
   
        # Calculating New center mark and the cam profile data
        arc=0.8*self.Rhole*exp(1j*self.alp)
        
        if self.turn_direction=='clockwise':
            z=(self.Xp+1j*self.Yp)*exp(-1j*self.theta[ii])            
            xdata=z.real
            ydata=z.imag
            arc=arc*exp(-1j*self.theta[ii]) 
        elif self.turn_direction=='anti-clockwise':
            z=(self.Xp+1j*self.Yp)*exp(1j*self.theta[ii])            
            xdata=z.real
            ydata=z.imag
            arc=arc*exp(1j*self.theta[ii])
        dataxy=array([xdata,ydata]).transpose()
        
        # Updating center mark data
        arc_data=array([arc.real,arc.imag]).transpose()               
        self.artc.set_xy(arc_data)
        
        # Updating Cam profile data
        self.poly1.set_xy(dataxy)
        
        # Updating Knife-Edge Follower
        dis=self.y[ii]*exp(1j*self.FollowerAng)
        xkna=self.xFollower+dis.real
        ykna=self.yFollower+dis.imag       
        
        dataknife=array([xkna,ykna]).transpose()      
        self.poly2.set_xy(dataknife)
        
        # Calculating and Updating Spring Data
        dmin=self.dmin+self.y[ii]
        Lmin=self.Lmin
        dmax=self.dmax    
        linksfgchain, linksbkchain=self.SpringBlock(dmax, dmin, Lmin,12)
        
        for link,jj in zip(self.springbk,arange(len(self.springbk))):
            link.set_xy(transpose(array([linksbkchain[:,jj].real,linksbkchain[:,jj].imag])))
        
        for link,jj in zip(self.springfg,arange(len(self.springfg))):
            link.set_xy(transpose(array([linksfgchain[:,jj].real,linksfgchain[:,jj].imag])))        
            
        return self.poly1, self.poly2

class PDCamRollerFollower:
    """ 
    PDCamRollerFollower(**CamData)
    
    Description
    ------------
    
    This is a design tool for Planar Disk Cam with Alternating Roller Follower
        (Centered an eccentric)
        
    The main functionalities are:
        1. It gives a plot with the Motion Diagram.
        2. It gives the cam's profile data in rectangular coordinates.
        3. The former data is ploted in a matplotlib Figure.
        4. Also there is a set of methods which help to animate cam's rotation.
    
    Parameters
    -----------
    
    The cam linkage parameters have to be provided via a python dictionary.
    
    theta : data-type
        1D NumPy array with N equal-spaced values between (0,2pi)
        This variable is the x-axis of the displacement diagram.
                
    y : data-type
        (1D NumPy array) Knife-Edge Follower displacement as function of theta.
        
    yp : data-type
        (1D NumPy array) Knife-Edge Follower velocity as function of theta.
        
    ypp : data-type
        (1D NumPy array) Knife-Edge Follower aceleration as function of theta.
        
    Rbase : float
        Base circle radius.
        
    Rroller : float
        Roller circle radius.

       
        Note:                
            y, yp and ypp have to be the same length.
       

    Optional Parameters:
    --------------------
        
    epsilon : float
        Alternating Roller Follower eccentricity 
        (positive for right position and negative otherwise)
        
    FollowerAng : float
        Follower's axis angle (radians)
        
    Followerwidth : float
        Alternating Roller Follower width
        
    Followerhight : float
        Alternating Roller Follower hight
        
    Rhole : float
        Hole drill radius. For make-up purpose.
        
    turn_direction : str
        Direction of rotation ('clockwise' or 'anti-clockwise')

    CamProfileColor : str or color array
        Color of the Cam Profile according to matplotlib's color format
        
    RollerFollowerColor : str or color array
        Color of the Follower according to matplotlib's color format
        
    RollerColor : str or color array
        Color of the roller according to matplotlib's color format
        
    BlockGuideColor : str or color array
        Color of the Block Guide according to matplotlib's color format
        
    CenterCamBgColor : str or color array
        Bakground color of the center Cam according to matplotlib's color format
        
    CenterCamFgColor : str or color array
        Foreground color of the center Cam according to matplotlib's color format
        
    
    Attributes
    -----------
    Xr : data-type
            1D NumPy array which contents a sequence of x-rectangular coordinates
            that describes the primitive trayectory (roller center path).
        
    Yr : data-type
            1D NumPy array which contents a sequence of y-rectangular coordinates
            that describes the primitive trayectory (roller center path).
    
    Xp : data-type
        1D NumPy array which contents a sequence of x-rectangular coordinates
        that describes the cam profile.
        
    Yp : data-type
        1D NumPy array which contents a sequence of y-rectangular coordinates
        that describes the cam profile.
        
    xFollower : data-type
        1D NumPy array which contents a sequence of x-rectangular coordinates
        that describes the Follower geometry at the initial condition.
        
    yFollower : data-type
        1D NumPy array which contents a sequence of y-rectangular coordinates
        that describes the Follower geometry at the initial condition.
        
    Xblock1: data-type
        1D NumPy array which contents a sequence of x-rectangular coordinates
        that partially describes the Follower's Guide.
    
    Yblock1 : data-type
        1D NumPy array which contents a sequence of y-rectangular coordinates
        that partially describes the Follower's Guide.
    
    Xblock2 : data-type
        1D NumPy array which contents a sequence of x-rectangular coordinates
        that partially describes the Follower's Guide.
    
    Yblock2 : data-type
        1D NumPy array which contents a sequence of y-rectangular coordinates
        that partially describes the Follower's Guide.
    
    
    Methods
    --------
                    
    PlotMotionDiagram() 
        It plots the position, velocity an acceleration of alternating follower
        as function of the angular displacement.
     
    CamForRollerFollower()
   
    """
    from numpy import sqrt
    
    def __init__(self,**kwargs):
        from numpy import pi
        """
        numpy and matplotlib have to be installed on the system
        """
        self.epsilon=0
        self.FollowerAng=pi/2
        self.Followerwidth=None
        self.Followerhight=None
        self.turn_direction='clockwise'
        self.CamProfileColor='goldenrod'
        self.RollerFollowerColor='brown'
        self.RollerColor=[0.1,0.1,0.4]
        self.BlockGuideColor='darkslategray'
        self.CenterCamBgColor=[0.1,0.1,0.4]
        self.CenterCamFgColor='yellowgreen'        
        """
        Updating Mandatory arguments
        """
        self.__dict__.update(**kwargs)         
        """
        Calculating the Cam Profile and initial data
        """ 
        
        self.updateparam()
        self.CamForRollerFollower()
        self.xFollower,self.yFollower=self.RollerFollower(self.Xr[0],self.Yr[0])
        self.Xblock1, self.Yblock1, self.Xblock2, self.Yblock2=self.FollowerGuide()

    def updateparam(self):
        from numpy import sqrt
        if self.Followerwidth==None:
            self.Followerwidth=max(self.y)/8
        if self.Followerhight==None:
            self.Followerhight=3*max(self.y)            
       
        self.Lmin=0.5*max(self.y)   
        self.dmax=(sqrt(self.Rbase**2-self.epsilon**2)
        +2.1*self.Rroller + 0.24*self.Followerwidth+self.Lmin+max(self.y))
        self.dmin=sqrt(self.Rbase**2-self.epsilon**2)+2.1*self.Rroller +0.24*self.Followerwidth
        
        
    def PlotMotionDiagram(self,fig):
        """PlotMotionDiagram(fig)
        
        Description:
        ------------
        It plots the Motion diagram associated with the Cam instance
        
        Parameters:
        -----------
            
            fig : FigureClass
                matplotlib figure instance
                
                How to get it example
                
                >>> from matplotlib.pyplot import figure
                
                >>> fig=figure()
                
        Return:
        --------
        
            axes : matplotlib.axes._subplots.AxesSubplot object array
                matplotlib axes array where the diagrams where plotted.
                More elements can be added by the user to such axes.
        """       
        axes=fig.subplots(3,1, sharex=True)    
        axes[0].plot(self.theta,self.y)
        axes[0].set_title('Displacement Diagram')
        axes[0].grid(True)
        axes[1].plot(self.theta,self.yp)
        axes[1].set_title('Velocity Diagram')
        axes[1].grid(True)
        axes[2].plot(self.theta,self.ypp)
        axes[2].set_title('Acceleration Diagram')
        axes[2].set(xlabel='Angle [radians]')
        axes[2].grid(True)              
        
        return axes  
   
    def CamForRollerFollower(self):
        """CamForRollerFollower()
        
        Description
        ------------
        
        It computes the rectangular coordinates of the Cam Profile for
        an Alternating Roller Follower.
        
        Returns
        --------
        
        Xr : data-type
            1D NumPy array which contents a sequence of x-rectangular coordinates
            that describes the primitive trayectory (roller center path).
        
        Yr : data-type
            1D NumPy array which contents a sequence of y-rectangular coordinates
            that describes the primitive trayectory (roller center path).
            
        Xp : data-type
            1D NumPy array which contents a sequence of x-rectangular coordinates
            that describes the cam profile.
        
        Yp : data-type
            1D NumPy array which contents a sequence of y-rectangular coordinates
            that describes the cam profile.
        """

        #Note: This method is called within the __init__ method and 
        #its values are stored in the Xp and Yp as attributes.

        from numpy import exp, arctan2, pi
        epsilon=self.epsilon
        direction=self.turn_direction
        Rbase=self.Rbase
        Rroller=self.Rroller
        theta,y=self.theta,self.y
        yp=self.yp
        beta=self.FollowerAng        
        q=((Rbase+Rroller)**2-epsilon**2)**0.5+y-1j*epsilon
        
        if direction=='clockwise':
            q2=yp+1j*q
            phi=arctan2(q2.imag,q2.real)
            A=Rroller*exp(1j*(phi+pi/2))
            p=q+A
            P=p*exp(1j*(beta+theta))            
            R=q*exp(1j*(beta+theta))
            Xr=R.real
            Yr=R.imag
            Xp=P.real
            Yp=P.imag
        elif direction=='anti-clockwise':
            q2=yp-1j*q
            phi=arctan2(q2.imag,q2.real)
            A=Rroller*exp(1j*(phi-pi/2))
            p=q+A
            P=p*exp(1j*(beta-theta))
            R=q*exp(1j*(beta-theta))            
            Xr=R.real
            Yr=R.imag
            Xp=P.real
            Yp=P.imag
        else:
            print('Please specify a valid direction of rotation')
        self.Xr,self.Yr,self.Xp,self.Yp=Xr,Yr,Xp,Yp
        self.Q=q # Complex distance From de Cam center to the Roller center
        self.P=p # Complex distance From de Cam center to the Roller tangent
        return self.Xr,self.Yr,self.Xp,self.Yp
    
    def PressureAngle(self):
        """
        PressureAngle()
        
        Description
        ------------
        It computes the pressure angle between the follower axis and the cam profile
        
        Return
        ______
        
        phi : data-type
            1D NumPy array which contents the pressure angle in radians for each
            cam's angular displacement.        
        """
        from numpy import absolute, arccos        
        P=self.P
        cosphi=P.real/absolute(P)
        phi=arccos(cosphi)        
        return phi
    
    def PlotPressureAngle(self,Axis):
        """
        PlotPressureAngle(Axis)
        
        Description
        ------------
        
        It plots the cam's pressure angle curve.
        
        Parameters
        -----------
        
        Axis : matplotlib.axes._subplots.AxesSubplot object
            matplotlib axis where the pressure angle curve will be plotted.

        Returns:
            
            line : matplotlib.lines.Line2D
                Gives a Line2D instance with x and y data en sequences xdata, ydata.        
        """
        from numpy import pi
        phi=self.PressureAngle()*180/pi
        line,=Axis.plot(self.theta,phi,color='blue')
        Axis.set_xlabel(r'$\theta$ [rad]')
        Axis.set_ylabel(r'Angle Pressure [degrees]')
        Axis.set_title('Angular displacement vs Angle Pressure')
        return line           

    def RollerFollower(self,xtip,ytip,angle=None):
        """
        RollerFollower(xtip,ytip,angle=None)
        
        Description
        ------------
        It gives the rectangular coordinates to render a Matplotlib PolyFill object
        to represent de alternating roller follower shape. This is an internal method 
        for use of other internal methods.
        
        Parameters
        ----------
        xtip : float
            Tip Follower's x-coordinate position
            
        ytip : float
            Tip Follower's y-coordinate position
            
        angle : float
            Angular position of the Follower, default value (None) have been set to the
            original Follower axis's angle (Default or given by the user).
            
        Returns
        --------
        (data-type,data-type) 
            Tuple with two 1D NumPy array which contents the (x,y)-coordinates
            of the follower shape at a given position.
        """

        from numpy import  exp, array,linspace, pi, concatenate
        wth=self.Followerwidth
        ht=self.Followerhight
        Rr=self.Rroller
        if angle==None:
            angle=self.FollowerAng
        xp1=array([1.1*Rr,1.1*Rr,1.1*Rr+0.24*wth,1.1*Rr+0.24*wth,ht,
                   ht,1.1*Rr+0.24*wth,1.1*Rr+0.24*wth,1.1*Rr,1.1*Rr])
        yp1=array([0.5*wth,wth,wth,0.5*wth,0.5*wth,
                   -0.5*wth,-0.5*wth,-wth,-wth,-0.5*wth])        
        th=linspace(270,90,20)*pi/180
        scir=0.5*wth*exp(1j*th)        
        xrollerF0=concatenate((xp1,scir.real,array([1.1*Rr])))
        yrollerF0=concatenate((yp1,scir.imag,array([0.5*wth])))       
        Follower=(xrollerF0+1j*yrollerF0)*exp(1j*angle)+(xtip+1j*ytip)
        xFollower=Follower.real
        yFollower=Follower.imag
        return xFollower, yFollower
    
    def FollowerGuide(self,angle=None):
        """
        FollowerGuide(angle=None)
        
        Description
        ------------
        It gives the rectangular coordinates to render a Matplotlib PolyFill object
        to represent a simple guide shape. This is an internal method for use 
        of other internal methods.
        
        Parameter
        ----------
        
        angle : float
            Angular position of the axis Follower, default value (None) have been set to the
            original Follower axis's angle (Default or given by the user).
        
        Returns
        ---------------------------
        Xblock1: data-type
            1D NumPy array which contents a sequence of x-rectangular coordinates
            that partially describes the Follower's Guide.
            
        Xblock1: data-type
            1D NumPy array which contents a sequence of x-rectangular coordinates
            that partially describes the Follower's Guide.
    
        Yblock1 : data-type
            1D NumPy array which contents a sequence of y-rectangular coordinates
            that partially describes the Follower's Guide.
    
        Xblock2 : data-type
            1D NumPy array which contents a sequence of x-rectangular coordinates
            that partially describes the Follower's Guide.
        
        """
        from numpy import exp, array
          
        wth=self.Followerwidth
        dz=self.dmax
        if angle==None:
            angle=self.FollowerAng        
        xb=array([dz,1.4*dz,1.4*dz,dz])
        yb=array([0.5*wth,0.5*wth,2.5*wth,2.5*wth])
        Zblock1=(xb+1j*yb-1j*self.epsilon)*exp(1j*angle)
        Zblock2=(xb-1j*yb-1j*self.epsilon)*exp(1j*angle)
        Xblock1=Zblock1.real
        Yblock1=Zblock1.imag
        Xblock2=Zblock2.real
        Yblock2=Zblock2.imag
        return Xblock1, Yblock1, Xblock2, Yblock2
    
    def Linkagecoord(self,width,x1,y1,x2,y2,angle=None):
        """
        Internal method for computing spring's geometry parts
        """
        from numpy import pi, piecewise, linspace, cos, sin, absolute, exp
        eps=self.epsilon
        if angle==None:
            gamm=self.FollowerAng
        else:
            gamm=angle
        q0=x1+1j*y1
        q1=x2+1j*y2
        height=absolute(q1-q0)
        betas=linspace(0,2*pi,50)
        Xreal=0.5*width*cos(betas)
        Ximag=piecewise(betas,[betas<pi,betas>=pi],[lambda x: 
            (0.5*width*sin(x)+height),
                    lambda x: 0.5*width*sin(x)])
        X=Xreal+1j*Ximag
        W=1j*(q1-q0)*X/absolute(q1-q0)+q1
        Linkage=(W-1j*eps)*exp(1j*gamm)       
        
        return Linkage
    
    def SpringBlock(self, dmax, dmin, Lmin,N):
        """
        Internal method for computing spring geometry
        """
        from numpy import linspace, ones_like, arange, concatenate
        if self.Followerwidth==None:
            self.Followerwidth=max(self.y)/8
        wth=self.Followerwidth
        r=0.85*Lmin/(2*N)
        L=dmax-dmin
        p=(L-2*r)/(N-0.5)
        x1=linspace(dmin+r,dmax-0.5*p-r,N)
        y1=(0.75*wth)*ones_like(x1)
        x2=linspace(dmin+r+0.5*p,dmax-r,N)
        y2=-(0.75*wth)*ones_like(x2)

        for kk in arange(N+1):
            
            if kk==0:
                linkfw=self.Linkagecoord(2*r,x1[kk],y1[kk],x2[kk],y2[kk])
                linkfw=linkfw.reshape(len(linkfw),1)
                linksfwchain=linkfw
            elif kk<N:
                linkfw=self.Linkagecoord(2*r,x1[kk],y1[kk],x2[kk],y2[kk])
                linkfw=linkfw.reshape(len(linkfw),1)
                linksfwchain=concatenate((linksfwchain,linkfw),axis=1)
                
            if kk==0:
                linkbk=self.Linkagecoord(2*r,x1[kk],y1[kk],x1[kk],y2[kk])
                linkbk=linkbk.reshape(len(linkbk),1)
                linksbkchain=linkbk
            elif kk==N:
                linkbk=self.Linkagecoord(2*r,x2[kk-1],y1[kk-1],x2[kk-1],y2[kk-1])
                linkbk=linkbk.reshape(len(linkbk),1)
                linksbkchain=concatenate((linksbkchain,linkbk),axis=1)
            else:
                linkbk=self.Linkagecoord(2*r,x1[kk],y1[kk],x2[kk-1],y2[kk-1])
                linkbk=linkbk.reshape(len(linkbk),1)
                linksbkchain=concatenate((linksbkchain,linkbk),axis=1)
            
        return linksfwchain, linksbkchain
    
    def Circle(self,radius,Xcenter,Ycenter):
        from numpy import exp, linspace, pi
        th=linspace(0,360,50)*pi/180
        Circ=radius*exp(1j*th)+Xcenter+1j*Ycenter
        return Circ
        
        
    def PlotCamRollerFollower(self,fig,detailed=True):
        """ 
        PlotCamRollerFollower(fig,detailed=True)
        
        Description:
        ------------
        This plots the Cam Profile with the Roller Follower and
        a decorative spring.
        
        Parameters:
        -----------
            
            fig : FigureClass
                matplotlib figure instance
                
                How to get it example
                
                >>> from matplotlib.pyplot import figure
                
                >>> fig=figure()
                
            detailed : bool
                If True, all the decorative elements are depicted, otherwise only
                the cam profile is rendered with a center mark.
                
        Return:
        --------
        
            axis : matplotlib.axes._subplots.AxesSubplot object
                matplotlib axis where the forms where plotted.
                More elements can be added by the user to such axis.        
       """

        from numpy import cos,sin, arange
        Rbase,Rh,theta=self.Rbase,self.Rhole,self.theta
        xp,yp=self.Xp, self.Yp
        xr,yr=self.Xr, self.Yr
        xf,yf=self.xFollower, self.yFollower
        xb1,yb1=self.Xblock1, self.Yblock1
        xb2,yb2=self.Xblock2, self.Yblock2
        
        CRoller=self.Circle(self.Rroller,xr[0],yr[0])                        
        dmin=self.dmin+self.y[0]
        Lmin=self.Lmin
        dmax=self.dmax   
        linksfwchain, linksbkchain=self.SpringBlock(dmax, dmin, Lmin,12)        

        axis=fig.add_subplot(111)
        axis.fill(xp,yp,color=self.CamProfileColor) 
        axis.set_aspect('equal')
        if detailed==True:
            axis.plot(xr,yr,'b--')
            axis.plot(Rbase*cos(theta),Rbase*sin(theta),'w--')
            axis.fill(Rh*cos(theta),Rh*sin(theta),'w',linewidth=3)
            axis.fill(CRoller.real,CRoller.imag,color=self.RollerColor)
            for column in arange(linksbkchain.shape[1]):
                    axis.fill(linksbkchain[:,column].real,linksbkchain[:,column].imag,'k')            
            axis.fill(xf,yf,color=self.RollerFollowerColor)
            for column in arange(linksfwchain.shape[1]):
                axis.fill(linksfwchain[:,column].real,linksfwchain[:,column].imag,'gray')          
            axis.plot(0,0,'k+',linewidth=3)
            axis.fill(xb1,yb1, color=self.BlockGuideColor)
            axis.fill(xb2,yb2, color=self.BlockGuideColor)
        else:
            axis.plot(0,0,'w+',linewidth=3)
        return axis
            
        
    def initAnim(self, ax):
        """
        initAnim(ax)
        
        Description
        ------------
        Method to initialize the Animation Parameters.
        
        Parameters:
        --------
        
            ax : matplotlib.axes._subplots.AxesSubplot object
                matplotlib axis where the initial geometries will be plotted.        
        """
        from numpy import sin, cos, sqrt, arange, linspace, concatenate, exp, pi
        self.axAnimation=ax
        
        # Calculating and Render Follower's Guide Block
        
        self.axAnimation.fill(self.Xblock1,self.Yblock1,color=self.BlockGuideColor)
        self.axAnimation.fill(self.Xblock2,self.Yblock2,color=self.BlockGuideColor)
        
        # Render Cam Profile
        self.poly1,=self.axAnimation.fill(self.Xp,self.Yp,color=self.CamProfileColor)
        
        # Calculating and render the roller
        Roller=self.Circle(self.Rroller,self.Xr[0], self.Yr[0])
        xdatRoll=Roller.real
        ydatRoll=Roller.imag
        self.rollpol,=self.axAnimation.fill(xdatRoll,ydatRoll,
                                            color=self.RollerColor)
        
        # Calculating Spring group coordinates
        dmin=self.dmin+self.y[0]
        Lmin=self.Lmin
        dmax=self.dmax    
        linksfgchain, linksbkchain=self.SpringBlock(dmax, dmin, Lmin,12)
        
        # Render spring's rear view        
        self.springbk=[]
        for column in arange(linksbkchain.shape[1]):
                l1,=self.axAnimation.fill(linksbkchain[:,column].real,linksbkchain[:,column].imag,'k')
                self.springbk.append(l1)
        # Render Roller Follower
        self.poly2,=self.axAnimation.fill(self.xFollower,self.yFollower,
                                          color=self.RollerFollowerColor)
        # Render spring's front view
        
        self.springfg=[]
        for column in arange(linksfgchain.shape[1]):
                l1,=self.axAnimation.fill(linksfgchain[:,column].real,linksfgchain[:,column].imag,
                                          'gray')
                self.springfg.append(l1)
        
        # Calculating and render the center mark on the cam profile
        self.axAnimation.fill(self.Rhole*cos(self.theta),self.Rhole*sin(self.theta),
                              color=self.CenterCamBgColor)        
        alp1=linspace(45,135,10)*pi/180
        alp2=linspace(-45,-135,10)*pi/180
        self.alp=concatenate((alp1,alp2))
        arc=0.8*self.Rhole*exp(1j*self.alp)        
        self.artc,=self.axAnimation.fill(arc.real,arc.imag, color=self.CenterCamFgColor)
        
        # Setting the Axes limits according to the dimensions of the main parts of the animation
        LimCam=max(sqrt((self.Xp**2+self.Yp**2)))
        LimBlk1x=max(self.Xblock1)
        LimBlk2x=max(self.Xblock2)
        LimBlk1y=max(self.Yblock1)
        LimBlk2y=max(self.Yblock2)
        Limxp=1.05*max([LimCam,LimBlk1x,LimBlk2x])
        Limxn=1.05*min([-LimCam,LimBlk1x,LimBlk2x])
        Limyp=1.05*max([LimCam,LimBlk1y,LimBlk2y])
        Limyn=1.05*min([-LimCam,LimBlk1y,LimBlk2y])        
        self.axAnimation.set_xlim(Limxn,Limxp)
        self.axAnimation.set_ylim(Limyn,Limyp)
        self.axAnimation.set_aspect('equal')
        
        return self.poly1, self.poly2
        
    def __call__(self,k):
        """Method to update Animation data
        """
        from numpy import mod, exp, transpose, array, arange
        
        i=5*k
        ii=mod(i,len(self.theta))
   
        # Calculating New center mark and the cam profile data
        arc=0.8*self.Rhole*exp(1j*self.alp)
        
        if self.turn_direction=='clockwise':
            z=(self.Xp+1j*self.Yp)*exp(-1j*self.theta[ii])            
            xdata=z.real
            ydata=z.imag
            arc=arc*exp(-1j*self.theta[ii]) 
        elif self.turn_direction=='anti-clockwise':
            z=(self.Xp+1j*self.Yp)*exp(1j*self.theta[ii])            
            xdata=z.real
            ydata=z.imag
            arc=arc*exp(1j*self.theta[ii])
        # Updating Cam profile data
        dataxy=array([xdata,ydata]).transpose()
        self.poly1.set_xy(dataxy)
        
        # Updating center mark data
        arc_data=array([arc.real,arc.imag]).transpose()               
        self.artc.set_xy(arc_data)       
        
        # Updating Roller Follower data
        dis=self.y[ii]*exp(1j*self.FollowerAng)
        xkna=self.xFollower+dis.real
        ykna=self.yFollower+dis.imag
        datafollower=array([xkna,ykna]).transpose()        
        self.poly2.set_xy(datafollower)
        
        # Updating Spring Data
        dmin=self.dmin+self.y[ii]
        Lmin=self.Lmin
        dmax=self.dmax    
        linksfgchain, linksbkchain=self.SpringBlock(dmax, dmin, Lmin,12)
        
        for link,jj in zip(self.springbk,arange(len(self.springbk))):
            link.set_xy(transpose(array([linksbkchain[:,jj].real,linksbkchain[:,jj].imag])))
        
        for link,jj in zip(self.springfg,arange(len(self.springfg))):
            link.set_xy(transpose(array([linksfgchain[:,jj].real,linksfgchain[:,jj].imag])))
        
        # Updating  Roller Data
        
        Roller=self.Circle(self.Rroller,self.Xr[0]+dis.real, self.Yr[0]+dis.imag)
        xdatRoll=Roller.real
        ydatRoll=Roller.imag
        datafollower=array([xdatRoll,ydatRoll]).transpose()
        self.rollpol.set_xy(datafollower)
            
        return self.poly1, self.poly2
    
class PDCamFlatFaceFollower():    
    from numpy import sqrt
    """ 
    PDCamFlatFaceFollower(**CamData)

    
    Description
    ------------
    
    This is a design tool for Planar Disk Cam with Flat-Face Follower
        (Centered an eccentric)
        
    The main functionalities are:
        1. It gives a plot with the Motion Diagram.
        2. It gives the cam's profile data in rectangular coordinates.
        3. The former data is ploted in a matplotlib Figure.
        4. Also there is a set of methods which help to animate cam's rotation.
    
    Parameters
    -----------
    
    The cam linkage parameters have to be provided via a python dictionary.
    
    theta : data-type
        1D NumPy array with N equal-spaced values between (0,2pi)
        This variable is the x-axis of the displacement diagram.
                
    y : data-type
        (1D NumPy array) Flat-Face Follower displacement as function of theta.
        
    yp : data-type
        (1D NumPy array) Flat-Face Follower velocity as function of theta.
        
    ypp : data-type
        (1D NumPy array) Flat-Face Follower aceleration as function of theta.
        
    Rbase : float
        Base Circle Radius.      

       
        Note:                
            y, yp and ypp have to be the same length.
       

    Optional Parameters:
    --------------------
        
    epsilon : float
        Follower eccentricity (positive for right position and negative otherwise)
        
    FollowerAng : float
        Follower's axis angle (radians)
        
    Followerwidth : float
        Follower width
        
    Followerhight : float
        Follower hight
        
    Followertickness : float
        Flat-face tickness
        
    Rhole : float
        Hole drill radius. For make-up purpose.
        
    turn_direction : str
        Direction of rotation ('clockwise' or 'anti-clockwise')

    CamProfileColor : str or color array
        Color of the Cam Profile according to matplotlib's color format
        
    FlatFollowerColor : str or color array
        Color of the Knife-Edge according to matplotlib's color format
        
    BlockGuideColor : str or color array
        Color of the Block Guide according to matplotlib's color format
        
    CenterCamBgColor : str or color array
        Bakground color of the center Cam according to matplotlib's color format
        
    CenterCamFgColor : str or color array
        Foreground color of the center Cam according to matplotlib's color format
        
    
    Attributes
    -----------
    Xp : data-type
        1D NumPy array which contents a sequence of x-rectangular coordinates
        that describes the cam profile.
        
    Yp : data-type
        1D NumPy array which contents a sequence of y-rectangular coordinates
        that describes the cam profile.
        
    xFollower : data-type
        1D NumPy array which contents a sequence of x-rectangular coordinates
        that describes the Follower geometry at the initial condition.
        
    yFollower : data-type
        1D NumPy array which contents a sequence of y-rectangular coordinates
        that describes the Follower geometry at the initial condition.
        
    Xblock1: data-type
        1D NumPy array which contents a sequence of x-rectangular coordinates
        that partially describes the Follower's Guide.
    
    Yblock1 : data-type
        1D NumPy array which contents a sequence of y-rectangular coordinates
        that partially describes the Follower's Guide.
    
    Xblock2 : data-type
        1D NumPy array which contents a sequence of x-rectangular coordinates
        that partially describes the Follower's Guide.
    
    Yblock2 : data-type
        1D NumPy array which contents a sequence of y-rectangular coordinates
        that partially describes the Follower's Guide.
    
    
    Methods
    --------
                    
    PlotMotionDiagram() 
        It plots the position, velocity an acceleration of alternating follower
        as function of the angular displacement.
                         

   
    """
    
    def __init__(self,**kwargs):
        from numpy import pi
        """
        numpy and matplotlib have to be installed on the system
        """
        self.epsilon=0
        self.FollowerAng=pi/2
        self.Followerwidth=None
        self.Followerhight=None
        self.Followertickness=None
        self.turn_direction='clockwise'
        self.CamProfileColor='goldenrod'
        self.FlatFollowerColor='brown'
        self.BlockGuideColor='darkslategray'
        self.CenterCamBgColor=[0.1,0.1,0.4]
        self.CenterCamFgColor='yellowgreen'        
        """
        Updating Mandatory arguments
        """
        self.__dict__.update(**kwargs)         
        """
        Calculating the Cam Profile and initial data
        """
        self.setparam()        
        self.CamForFlatFollower()
        self.xFollower, self.yFollower=self.FlatFollower(self.Xp[0],self.Yp[0])               
        self.Xblock1, self.Yblock1, self.Xblock2, self.Yblock2=self.FollowerGuide()
        
    def setparam(self):
        
        if self.Followerwidth==None:
            self.Followerwidth=max(self.y)/8
        if self.Followerhight==None:
            self.Followerhight=3*max(self.y)
        if self.Followertickness==None:
            self.Followertickness=0.25*self.Followerwidth
        self.Lmin=0.5*max(self.y)   
        self.dmax=self.Rbase+self.Followertickness+self.Lmin+max(self.y)
        self.dmin=self.Rbase+self.Followertickness
        
    def PlotMotionDiagram(self,fig):
        """PlotMotionDiagram(fig)
        
        Description:
        ------------
        It plots the Motion diagram associated with the Cam instance
        
        Parameters:
        -----------
            
            fig : FigureClass
                matplotlib figure instance
                
                How to get it example
                
                >>> from matplotlib.pyplot import figure
                
                >>> fig=figure()
                
        Return:
        --------
        
            axes : matplotlib.axes._subplots.AxesSubplot object array
                matplotlib axes array where the diagrams where plotted.
                More elements can be added by the user to such axes.
        """   
        
        axes=fig.subplots(3,1, sharex=True)    
        axes[0].plot(self.theta,self.y)
        axes[0].set_title('Displacement Diagram')
        axes[0].grid(True)
        axes[1].plot(self.theta,self.yp)
        axes[1].set_title('Velocity Diagram')
        axes[1].grid(True)
        axes[2].plot(self.theta,self.ypp)
        axes[2].set_title('Acceleration Diagram')
        axes[2].set(xlabel='Angle [radians]')
        axes[2].grid(True)              
        
        return axes  
        
    def CamForFlatFollower(self):
        """CamForFlatFollower()
        
        Description
        ------------
        
        It computes the rectangular coordinates of the Cam Profile for
        a Flat-Face Follower.
        
        Returns
        --------
        
        Xp : data-type
            1D NumPy array which contents a sequence of x-rectangular coordinates
            that describes the cam profile.
        
        Yp : data-type
            1D NumPy array which contents a sequence of y-rectangular coordinates
            that describes the cam profile.        
        """
        #Note: This method is called within the __init__method and 
        #its values are stored in the Xp and Yp as attributes.

        from numpy import exp
        direction=self.turn_direction
        Rbase=self.Rbase
        theta,y=self.theta,self.y
        s=self.yp
        beta=self.FollowerAng        
        
        if direction=='clockwise':
            Q=(Rbase+y)+1j*s            
            R=Q*exp(1j*(beta+theta))
            Xp=R.real
            Yp=R.imag
        elif direction=='anti-clockwise':
            Q=(Rbase+y)-1j*s
            R=Q*exp(1j*(beta-theta))
            Xp=R.real
            Yp=R.imag
        else:
            print('Please specify a valid direction of rotation')
        self.Xp,self.Yp=Xp,Yp
        self.Q=Q # Complex distance From de Cam center to the Flat follower tip
        return self.Xp,self.Yp
    
    def PressureAngle(self):
        """
        PressureAngle()
        
        Description
        ------------
        It computes the pressure angle between the follower axis and the cam profile
        
        Return
        ______
        
        phi : data-type
            1D NumPy array which contents the pressure angle in radians for each
            cam's angular displacement.        
        """
        from numpy import absolute, arccos        
        Q=self.Q
        cosphi=Q.real/absolute(Q)
        phi=arccos(cosphi)        
        return phi
    
    def PlotPressureAngle(self,Axis):
        """
        PlotPressureAngle(Axis)
        
        Description
        ------------
        
        It plots the cam's pressure angle curve.
        
        Parameters
        -----------
        
        Axis : matplotlib.axes._subplots.AxesSubplot object
            matplotlib axis where the pressure angle curve will be plotted.

        Returns:
            
            line : matplotlib.lines.Line2D
                Gives a Line2D instance with x and y data en sequences xdata, ydata.        
        """
        from numpy import pi
        phi=self.PressureAngle()*180/pi
        line,=Axis.plot(self.theta,phi,color='blue')
        Axis.set_xlabel(r'$\theta$ [rad]')
        Axis.set_ylabel(r'Angle Pressure [degrees]')
        Axis.set_title('Angular displacement vs Angle Pressure')
        return line
        
    def FlatFollower(self,xtip,ytip,angle=None):
        """
        FlatFollower(xtip,ytip,angle=None)
        
        Description
        ------------
        Gives the rectangular coordinates to render a Matplotlib PolyFill object
        to represent de alternating flat-face follower shape. This is an internal method 
        for use of other internal methods.
        
        Parameters
        ----------
        xtip : float
            Tip Follower's x-coordinate position
            
        ytip : float
            Tip Follower's y-coordinate position
            
        angle : float
            Angular position of the Follower, default value (None) have been set to the
            original Follower axis's angle (Default or given by the user).
            
        Returns
        --------
        (data-type,data-type) 
            Tuple with two 1D NumPy array which contents the (x,y)-coordinates of the follower shape at
            a given position.
        """
        from numpy import  exp, array, linspace, pi, concatenate    
        wth=self.Followerwidth
        ht=self.Followerhight
        thick=self.Followertickness
        ypmax=1.2*max(self.yp)
        ypmin=1.2*min(self.yp)
        if angle==None:
            angle=self.FollowerAng
        th1=linspace(0,-90,10)*pi/180
        th2=linspace(90,0,10)*pi/180
            
        xp1=array([thick,ht,ht,thick])
        xp2=thick*exp(1j*th1).real
        xp3=thick*exp(1j*th2).real
        
        yp1=array([0.5*wth,0.5*wth,-0.5*wth,-0.5*wth])-self.epsilon
        if self.turn_direction=='anti-clockwise':
            yp2=thick*exp(1j*th1).imag-ypmax
            yp3=thick*exp(1j*th2).imag-ypmin
        else:
            yp2=thick*exp(1j*th1).imag+ypmin
            yp3=thick*exp(1j*th2).imag+ypmax
        
        xFollower0=concatenate((xp1,xp2,xp3))
        yFollower0=concatenate((yp1,yp2,yp3))
        
        Follower=(xFollower0+1j*yFollower0)*exp(1j*angle)+(xtip+1j*ytip)
        xFollower=Follower.real
        yFollower=Follower.imag
        return xFollower, yFollower
    
    def FollowerGuide(self,angle=None):
        """
        FollowerGuide(angle=None)
        
        Description
        ------------
        It gives the rectangular coordinates to render a Matplotlib PolyFill object
        to represent a simple guide shape. This is an internal method for use 
        of other internal methods.
        
        Parameter
        ----------
        
        angle : float
            Angular position of the axis Follower, default value (None) have been set to the
            original Follower axis's angle (Default or given by the user).
        
        Returns
        ---------------------------
        Xblock1: data-type
            1D NumPy array which contents a sequence of x-rectangular coordinates
            that partially describes the Follower's Guide.
    
        Yblock1 : data-type
            1D NumPy array which contents a sequence of y-rectangular coordinates
            that partially describes the Follower's Guide.
    
        Xblock2 : data-type
            1D NumPy array which contents a sequence of x-rectangular coordinates
            that partially describes the Follower's Guide.
    
        Yblock2 : data-type
            1D NumPy array which contents a sequence of y-rectangular coordinates
            that partially describes the Follower's Guide.
        
        """
        from numpy import exp, array
        if self.Followerwidth==None:
            self.Followerwidth=max(self.y)/8            
        wth=self.Followerwidth
        dz=self.dmax
        if angle==None:
            angle=self.FollowerAng        
        xb=array([dz,1.4*dz,1.4*dz,dz])
        yb=array([0.5*wth,0.5*wth,2.5*wth,2.5*wth])
        Zblock1=(xb+1j*yb-1j*self.epsilon)*exp(1j*angle)
        Zblock2=(xb-1j*yb-1j*self.epsilon)*exp(1j*angle)
        Xblock1=Zblock1.real
        Yblock1=Zblock1.imag
        Xblock2=Zblock2.real
        Yblock2=Zblock2.imag
        return Xblock1, Yblock1, Xblock2, Yblock2
    
    def Linkagecoord(self,width,x1,y1,x2,y2,angle=None):
        """
        Internal method for computing spring's geometry parts
        """
        from numpy import pi, piecewise, linspace, cos, sin, absolute, exp
        eps=self.epsilon
        if angle==None:
            gamm=self.FollowerAng
        else:
            gamm=angle
        q0=x1+1j*y1
        q1=x2+1j*y2
        height=absolute(q1-q0)
        betas=linspace(0,2*pi,50)
        Xreal=0.5*width*cos(betas)
        Ximag=piecewise(betas,[betas<pi,betas>=pi],[lambda x: 
            (0.5*width*sin(x)+height),
                    lambda x: 0.5*width*sin(x)])
        X=Xreal+1j*Ximag
        W=1j*(q1-q0)*X/absolute(q1-q0)+q1
        Linkage=(W-1j*eps)*exp(1j*gamm)       
        
        return Linkage
    
    def SpringBlock(self, dmax, dmin, Lmin,N):
        """
        Internal method for computing spring geometry
        """
        from numpy import linspace, ones_like, arange, concatenate
        wth=self.Followerwidth
        r=0.85*Lmin/(2*N)
        L=dmax-dmin
        p=(L-2*r)/(N-0.5)
        x1=linspace(dmin+r,dmax-0.5*p-r,N)
        y1=(0.75*wth)*ones_like(x1)
        x2=linspace(dmin+r+0.5*p,dmax-r,N)
        y2=-(0.75*wth)*ones_like(x2)

        for kk in arange(N+1):
            
            if kk==0:
                linkfw=self.Linkagecoord(2*r,x1[kk],y1[kk],x2[kk],y2[kk])
                linkfw=linkfw.reshape(len(linkfw),1)
                linksfwchain=linkfw
            elif kk<N:
                linkfw=self.Linkagecoord(2*r,x1[kk],y1[kk],x2[kk],y2[kk])
                linkfw=linkfw.reshape(len(linkfw),1)
                linksfwchain=concatenate((linksfwchain,linkfw),axis=1)
                
            if kk==0:
                linkbk=self.Linkagecoord(2*r,x1[kk],y1[kk],x1[kk],y2[kk])
                linkbk=linkbk.reshape(len(linkbk),1)
                linksbkchain=linkbk
            elif kk==N:
                linkbk=self.Linkagecoord(2*r,x2[kk-1],y1[kk-1],x2[kk-1],y2[kk-1])
                linkbk=linkbk.reshape(len(linkbk),1)
                linksbkchain=concatenate((linksbkchain,linkbk),axis=1)
            else:
                linkbk=self.Linkagecoord(2*r,x1[kk],y1[kk],x2[kk-1],y2[kk-1])
                linkbk=linkbk.reshape(len(linkbk),1)
                linksbkchain=concatenate((linksbkchain,linkbk),axis=1)
            
        return linksfwchain, linksbkchain
        
        
    def PlotCamFlatFollower(self,fig,detailed=True):
        """ 
        PlotCamFlatFollower(fig,detailed=True)
        
        Description:
        ------------
        This plots the Cam Profile with the Flat-Face Follower and
        a decorative spring.
        
        Parameters:
        -----------
            
            fig : FigureClass
                matplotlib figure instance
                
                How to get it example
                
                >>> from matplotlib.pyplot import figure
                
                >>> fig=figure()
                
            detailed : bool
                If True, all the decorative elements are depicted, otherwise only
                the cam profile is rendered with a center mark.
                
        Return:
        --------
        
            Axis : matplotlib.axes._subplots.AxesSubplot object
                matplotlib axis where the forms where plotted.
                More elements can be added by the user to such axis.        
       """

        from numpy import cos,sin, arange
        Rbase,Rh,theta=self.Rbase,self.Rhole,self.theta
        xp,yp=self.Xp, self.Yp
        xf,yf=self.xFollower, self.yFollower
        xb1,yb1=self.Xblock1, self.Yblock1
        xb2,yb2=self.Xblock2, self.Yblock2
                        
        dmin=self.dmin+self.y[0]
        Lmin=self.Lmin
        dmax=self.dmax   
        linksfwchain, linksbkchain=self.SpringBlock(dmax, dmin, Lmin,12)        

        axis=fig.add_subplot(111)
        axis.fill(xp,yp,color=self.CamProfileColor) 
        axis.set_aspect('equal')
        if detailed==True:
            axis.plot(Rbase*cos(theta),Rbase*sin(theta),'w--')
            axis.fill(Rh*cos(theta),Rh*sin(theta),'w',linewidth=3)
            for column in arange(linksbkchain.shape[1]):
                axis.fill(linksbkchain[:,column].real,linksbkchain[:,column].imag,'k')            
            axis.fill(xf,yf,color=self.FlatFollowerColor)
            for column in arange(linksfwchain.shape[1]):
                axis.fill(linksfwchain[:,column].real,linksfwchain[:,column].imag,'gray')          
            axis.plot(0,0,'k+',linewidth=3)
            axis.fill(xb1,yb1, color=self.BlockGuideColor)
            axis.fill(xb2,yb2, color=self.BlockGuideColor)
        else:
            axis.plot(0,0,'w+',linewidth=3)
            
        
    def initAnim(self,ax):
        """
        initAnim(ax)
        
        Description
        ------------
        Method to initialize the Animation Parameters.
        
        Parameters:
        --------
        
            ax : matplotlib.axes._subplots.AxesSubplot object
                matplotlib axis where the initial geometries will be plotted.        
        """
        from numpy import sin, cos, sqrt, arange, linspace, concatenate, exp, pi
        
        dmin=self.dmin+self.y[0]
        Lmin=self.Lmin
        dmax=self.dmax    
        linksfgchain, linksbkchain=self.SpringBlock(dmax, dmin, Lmin,12)
        self.axAnimation=ax
        colorblk=self.BlockGuideColor
        Xblock1=self.Xblock1
        Yblock1=self.Yblock1
        Xblock2=self.Xblock2
        Yblock2=self.Yblock2
        # Render Follower's Guide Block
        self.axAnimation.fill(Xblock1,Yblock1,color=colorblk)
        self.axAnimation.fill(Xblock2,Yblock2,color=colorblk)        
        # Render Cam Profile
        self.poly1,=self.axAnimation.fill(self.Xp,self.Yp,color=self.CamProfileColor)
        # Render spring's rear view
        self.springbk=[]
        for col in arange(linksbkchain.shape[1]):
                l1,=self.axAnimation.fill(linksbkchain[:,col].real,linksbkchain[:,col].imag,'k')
                self.springbk.append(l1)
        # Render Flat-Face Follower
        self.poly2,=self.axAnimation.fill(self.xFollower,self.yFollower, color=self.FlatFollowerColor)
        # Render spring's front view
        self.springfg=[]
        for col in arange(linksfgchain.shape[1]):
                l1,=self.axAnimation.fill(linksfgchain[:,col].real,linksfgchain[:,col].imag,'gray')
                self.springfg.append(l1)
        self.axAnimation.fill(self.Rhole*cos(self.theta),self.Rhole*sin(self.theta),
                              color=self.CenterCamBgColor)
        
        # Calculating and render the center mark on the cam profile
        alp1=linspace(45,135,10)*pi/180
        alp2=linspace(-45,-135,10)*pi/180
        self.alp=concatenate((alp1,alp2))
        arc=0.8*self.Rhole*exp(1j*self.alp)        
        self.artc,=self.axAnimation.fill(arc.real,arc.imag, color=self.CenterCamFgColor)
        
        # Setting the Axes limits according to the dimensions of the main parts of the animation
        LimCam=max(sqrt((self.Xp**2+self.Yp**2)))
        LimBlk1x=max(self.Xblock1)
        LimBlk2x=max(self.Xblock2)
        LimBlk1y=max(self.Yblock1)
        LimBlk2y=max(self.Yblock2)
        Limxp=1.05*max([LimCam,LimBlk1x,LimBlk2x])
        Limxn=1.05*min([-LimCam,LimBlk1x,LimBlk2x])
        Limyp=1.05*max([LimCam,LimBlk1y,LimBlk2y])
        Limyn=1.05*min([-LimCam,LimBlk1y,LimBlk2y])
        
        self.axAnimation.set_xlim(Limxn,Limxp)
        self.axAnimation.set_ylim(Limyn,Limyp)
        self.axAnimation.set_aspect('equal')
        return self.poly1, self.poly2
        
    def __call__(self,k):
        """Method to update Animation data
        """
        from numpy import mod, exp, transpose, array, arange
        
        i=5*k
        ii=mod(i,len(self.theta))
   
        # Calculating New center mark and the cam profile data
        arc=0.8*self.Rhole*exp(1j*self.alp)
        
        if self.turn_direction=='clockwise':
            z=(self.Xp+1j*self.Yp)*exp(-1j*self.theta[ii])            
            xdata=z.real
            ydata=z.imag
            arc=arc*exp(-1j*self.theta[ii]) 
        elif self.turn_direction=='anti-clockwise':
            z=(self.Xp+1j*self.Yp)*exp(1j*self.theta[ii])            
            xdata=z.real
            ydata=z.imag
            arc=arc*exp(1j*self.theta[ii])
        dataxy=array([xdata,ydata]).transpose()
        
        # Updating center mark data
        arc_data=array([arc.real,arc.imag]).transpose()               
        self.artc.set_xy(arc_data)
        
        # Updating Cam profile data
        self.poly1.set_xy(dataxy)
        
        dis=self.y[ii]*exp(1j*self.FollowerAng)
        xkna=self.xFollower+dis.real
        ykna=self.yFollower+dis.imag
        
        # Updating Flat-Face Follower
        datafollower=array([xkna,ykna]).transpose()      
        self.poly2.set_xy(datafollower)
        
        # Calculating and Updating Spring Data
        dmin=self.dmin+self.y[ii]
        Lmin=self.Lmin
        dmax=self.dmax    
        linksfgchain, linksbkchain=self.SpringBlock(dmax, dmin, Lmin,12)
        
        for link,jj in zip(self.springbk,arange(len(self.springbk))):
            link.set_xy(transpose(array([linksbkchain[:,jj].real,linksbkchain[:,jj].imag])))
        
        for link,jj in zip(self.springfg,arange(len(self.springfg))):
            link.set_xy(transpose(array([linksfgchain[:,jj].real,linksfgchain[:,jj].imag])))        
            
        return self.poly1, self.poly2

class PDCamOscillatingRollerFollower:
    """ 
    PDCamRollerFollower(**CamData)
    
    Description
    ------------
    
    This is a design tool for Planar Disk Cam with Angular Oscillating Roller Follower
        (Centered an eccentric)
        
    The main functionalities are:
        1. It gives a plot with the Motion Diagram.
        2. It gives the cam's profile data in rectangular coordinates.
        3. The former data is ploted in a matplotlib Figure.
        4. Also there is a set of methods which help to animate cam's rotation.
    
    Parameters
    -----------
    
    The cam linkage parameters have to be provided via a python dictionary.
    
    theta : data-type
        1D NumPy array with N equal-spaced values between (0,2pi)
        This variable is the x-axis of the displacement diagram.
                
    phi : data-type
        (1D NumPy array) Follower's angular displacement as function of theta.
        
    phip : data-type
        (1D NumPy array) Follower's angular velocity as function of theta.
        
    phipp : data-type
        (1D NumPy array) Follower's angular aceleration as function of theta.
        
    Rbase : float
        Base circle radius.
        
    Rroller : float
        Roller circle radius.
        
    Rpiv : float
        Radial distance from de cam center to the Follower's pivot position.
        
    Rl : float
        Radial distance from de Follower's pivot position to the roller center.
          
        Note:                
            phi, phip and phipp have to be the same length.
       

    Optional Parameters:
    --------------------        

    PivAngpos : float
        Pivot angular position , default pi/2 (radians)
        
    Pivpos : float
        Location of the pivot ('right' or 'left'), try both options to see the difference.
        
    Followerwidth : float
        Oscillating Roller Follower width        
      
    Rhole : float
        Center mark radius. For make-up purpose.
        
    turn_direction : str
        Direction of rotation ('clockwise' or 'anti-clockwise')

    CamProfileColor : str or color array
        Color of the Cam Profile according to matplotlib's color format
        
    RollerFollowerColor : str or color array
        Color of the Follower according to matplotlib's color format
        
    RollerColor : str or color array
        Color of the roller according to matplotlib's color format
        
    PivotColor : str or color array
        Color of the Pivot block according to matplotlib's color format
        
    CenterCamBgColor : str or color array
        Bakground color of the center Cam according to matplotlib's color format
        
    CenterCamFgColor : str or color array
        Foreground color of the center Cam according to matplotlib's color format
        
    
    Attributes
    -----------
    Xr : data-type
            1D NumPy array which contents a sequence of x-rectangular coordinates
            that describes the primitive trayectory (roller center path).
        
    Yr : data-type
            1D NumPy array which contents a sequence of y-rectangular coordinates
            that describes the primitive trayectory (roller center path).
    
    Xp : data-type
        1D NumPy array which contents a sequence of x-rectangular coordinates
        that describes the cam profile.
        
    Yp : data-type
        1D NumPy array which contents a sequence of y-rectangular coordinates
        that describes the cam profile.
        
    xFollower : data-type
        1D NumPy array which contents a sequence of x-rectangular coordinates
        that describes the Follower geometry at the initial condition.
        
    yFollower : data-type
        1D NumPy array which contents a sequence of y-rectangular coordinates
        that describes the Follower geometry at the initial condition. 
    
    Methods
    --------
                    
    PlotMotionDiagram() 
        It plots the position, velocity an acceleration of alternating follower
        as function of the angular displacement.
     
  
   
    """
    from numpy import sqrt
    
    def __init__(self,**kwargs):
        from numpy import pi, exp
        """
        numpy and matplotlib have to be installed on the system
        """
        self.PivAngpos=pi/2
        self.Pivpos='right'
        self.Followerwidth=None
        self.turn_direction='clockwise'
        self.CamProfileColor='goldenrod'
        self.RollerFollowerColor='brown'
        self.RollerColor='olive'
        self.PivotColor='darkslategray'
        self.CenterCamBgColor=[0.1,0.1,0.4]
        self.CenterCamFgColor='yellowgreen'        
        """
        Updating Input arguments
        """
        self.__dict__.update(**kwargs)         
        """
        Calculating the Cam Profile and initial data
        """ 
        
        self.updateparam()
        self.Xr,self.Yr,self.Xp,self.Yp=self.CamForOscRollerFollower()
        q0=self.Rpiv*exp(1j*self.PivAngpos)
        xpiv=q0.real
        ypiv=q0.imag
        xend=self.Xr[0]
        yend=self.Yr[0]
        self.xFollower, self.yFollower=self.OscRollerFollower(xpiv,ypiv,xend,yend)
        self.xpiv=xpiv
        self.ypiv=ypiv

    def updateparam(self):
        if self.Followerwidth==None:
            self.Followerwidth=self.Rroller     
        
    def PlotMotionDiagram(self,fig):
        """PlotMotionDiagram(fig)
        
        Description:
        ------------
        It plots the Motion diagram associated with the Cam instance
        
        Parameters:
        -----------
            
            fig : FigureClass
                matplotlib figure instance
                
                How to get it example
                
                >>> from matplotlib.pyplot import figure
                
                >>> fig=figure()
                
        Return:
        --------
        
            axes : matplotlib.axes._subplots.AxesSubplot object array
                matplotlib axes array where the diagrams where plotted.
                More elements can be added by the user to such axes.
        """        
        axes=fig.subplots(3,1, sharex=True)    
        axes[0].plot(self.theta,self.phi)
        axes[0].set_title('Angular Displacement Diagram')
        axes[0].grid(True)
        axes[1].plot(self.theta,self.phip)
        axes[1].set_title('Angular Velocity Diagram')
        axes[1].grid(True)
        axes[2].plot(self.theta,self.phipp)
        axes[2].set_title('Angular Acceleration Diagram')
        axes[2].set(xlabel='Angle [radians]')
        axes[2].grid(True)              
        
        return axes  
   
    def CamForOscRollerFollower(self):
        """CamOscRollerFollower()
        
        Description
        ------------
        
        It computes the rectangular coordinates of the Cam Profile for
        a Angular Oscillating Roller Follower.
        
        Returns
        --------
        
        Xp : data-type
            1D NumPy array which contents a sequence of x-rectangular coordinates
            that describes the cam profile.
        
        Yp : data-type
            1D NumPy array which contents a sequence of y-rectangular coordinates
            that describes the cam profile.        
        """

        #Note: This method is called within the __init__method and 
        #its values are stored in the Xp and Yp as attributes.

        from numpy import exp, arctan2, arccos
        phi=self.phi
        phip=self.phip
        th=self.theta
        Rb=self.Rbase
        Rro=self.Rroller
        Rl=self.Rl
        Rpiv=self.Rpiv
        gamma=self.PivAngpos
        Pivpos=self.Pivpos       
        direction=self.turn_direction
        
        alpha0=arccos((Rpiv**2+Rl**2-(Rb+Rro)**2)/(2*Rpiv*Rl))
        
        if direction=='clockwise':
            rot=1
        elif direction=='anti-clockwise':
            rot=-1        
        else:
            print('Please specify a valid direction of rotation')
            print('Default is taken')
            rot=1
            
        if Pivpos=='right':
            pos=1
        elif Pivpos=='left':
            pos=-1
        else:
            print('Please specify a valid position (''left'' or ''right'')')
            print('Default is taken')
            pos=1
            
        Zd=-Rpiv+(Rl+Rl*pos*rot*phip)*exp(1j*(alpha0+phi)*pos)            
        delta=arctan2(Zd.imag,Zd.real)
        
        a=Rro*exp(1j*delta)
        q=Rpiv-Rl*exp(1j*(alpha0+phi)*pos)
        r=q+a
        Q=q*exp(1j*(th*rot+gamma))
        R=r*exp(1j*(th*rot+gamma))
        Xr=Q.real
        Yr=Q.imag
        Xp=R.real
        Yp=R.imag           
        self.q=q
        return Xr,Yr,Xp,Yp
    
    
    
    def OscRollerFollower(self,xpiv,ypiv,xend,yend):
        """
        OscRollerFollower(xtip,ytip,angle=None)
        
        Description
        ------------
        Gives the rectangular coordinates to render a Matplotlib PolyFill object
        to represent de angular alternating roller follower shape. This is an internal method 
        for use of other internal methods.
        
        Parameters
        ----------
        xtip : float
            Tip Follower's x-coordinate position
            
        ytip : float
            Tip Follower's y-coordinate position
            
        angle : float
            Angular position of the Follower, default value (None) have been set to the
            original Follower axis's angle (Default or given by the user).
            
        Returns
        --------
        (data-type,data-type) 
            Tuple with two 1D NumPy array which contents the (x,y)-coordinates of the follower shape at
            a given position.
        """
        from numpy import  pi, piecewise, linspace, cos, sin, absolute
        width=self.Followerwidth        
        q0=xpiv+1j*ypiv        
        q1=xend+1j*yend
        height=absolute(q1-q0)
        betas=linspace(0,2*pi,50)
        Xreal=0.5*width*cos(betas)
        Ximag=piecewise(betas,[betas<pi,betas>=pi],[lambda x: 
            (0.5*width*sin(x)+height),
                    lambda x: 0.5*width*sin(x)])
        X=Xreal+1j*Ximag
        Follower=1j*(q1-q0)*X/height+q1      
        xFollower=Follower.real
        yFollower=Follower.imag
        return xFollower, yFollower
    
    
    def Circle(self,radius,Xcenter,Ycenter):
        from numpy import exp, linspace, pi
        th=linspace(0,360,50)*pi/180
        Circ=radius*exp(1j*th)+Xcenter+1j*Ycenter
        return Circ
        
        
    def PlotCamOscRollerFollower(self,fig,detailed=True):
        """ 
        PlotCamOscRollerFollower(fig,detailed=True)
        
        Description:
        ------------
        This plots the Cam Profile with the Angular Alternating Roller Follower.
        
        Parameters:
        -----------
            
            fig : FigureClass
                matplotlib figure instance
                
                How to get it example
                
                >>> from matplotlib.pyplot import figure
                
                >>> fig=figure()
                
            detailed : bool
                If True, all the decorative elements are depicted, otherwise only
                the cam profile is rendered with a center mark.
                
        Return:
        --------
        
            Axis : matplotlib.axes._subplots.AxesSubplot object
                matplotlib axis where the forms where plotted.
                More elements can be added by the user to such axis.        
       """
        from numpy import cos,sin
        Rbase,Rh,theta=self.Rbase,self.Rhole,self.theta
        xp,yp=self.Xp, self.Yp
        xr,yr=self.Xr, self.Yr
        xf,yf=self.xFollower, self.yFollower
        
        CRoller=self.Circle(self.Rroller,xr[0],yr[0])                        
        axis=fig.add_subplot(111)
        axis.fill(xp,yp,color=self.CamProfileColor) 
        axis.set_aspect('equal')
        if detailed==True:
            axis.plot(xr,yr,'b--')
            axis.plot(Rbase*cos(theta),Rbase*sin(theta),'w--')
            axis.fill(Rh*cos(theta),Rh*sin(theta),'w',linewidth=3)
            axis.fill(CRoller.real,CRoller.imag,color=self.RollerColor)
            axis.fill(xf,yf,color=self.RollerFollowerColor)
            axis.plot(0,0,'k+',linewidth=3)
        else:
            axis.plot(0,0,'w+',linewidth=3)
            
        return axis
            
        
    def initAnim(self, ax):
        """
        initAnim(ax)
        
        Description
        ------------
        Method to initialize the Animation Parameters.
        
        Parameters:
        --------
        
            ax : matplotlib.axes._subplots.AxesSubplot object
                matplotlib axis where the initial geometries will be plotted.        
        """
        from numpy import sin, cos, sqrt,linspace, concatenate, exp, pi
        self.axAnimation=ax
        
        # Render Cam Profile
        self.poly1,=self.axAnimation.fill(self.Xp,self.Yp,color=self.CamProfileColor)
        
        # Calculating and render the roller
        Roller=self.Circle(self.Rroller,self.Xr[0], self.Yr[0])
        xdatRoll=Roller.real
        ydatRoll=Roller.imag
        self.rollpol,=self.axAnimation.fill(xdatRoll,ydatRoll,
                                            color=self.RollerColor)
        
        # Render Roller Follower
        self.poly2,=self.axAnimation.fill(self.xFollower,self.yFollower,
                                          color=self.RollerFollowerColor)
        
        # Calculating and render the center mark on the cam profile
        self.axAnimation.fill(self.Rhole*cos(self.theta),self.Rhole*sin(self.theta),
                              color=self.CenterCamBgColor)        
        alp1=linspace(45,135,10)*pi/180
        alp2=linspace(-45,-135,10)*pi/180
        self.alp=concatenate((alp1,alp2))
        arc=0.8*self.Rhole*exp(1j*self.alp)        
        self.artc,=self.axAnimation.fill(arc.real,arc.imag, color=self.CenterCamFgColor)
        
        # Setting the Axes limits according to the dimensions of the main parts of the animation
        LimCam=max(sqrt((self.Xp**2+self.Yp**2)))
        LimFllwerx=max(self.xFollower)
        LimFllwery=max(self.yFollower)
        Limxp=1.05*max([LimCam,LimFllwerx])
        Limxn=1.05*min([-LimCam,LimFllwerx])
        Limyp=1.05*max([LimCam,LimFllwery])
        Limyn=1.05*min([-LimCam,LimFllwery])        
        self.axAnimation.set_xlim(Limxn,Limxp)
        self.axAnimation.set_ylim(Limyn,Limyp)
        self.axAnimation.set_aspect('equal')
        
        return self.poly1, self.poly2
        
    def __call__(self,k):
        """Method to update Animation data
        """
        from numpy import mod, exp, array
        
        i=5*k
        ii=mod(i,len(self.theta))
   
        # Calculating New center mark and the cam profile data
        arc=0.8*self.Rhole*exp(1j*self.alp)
        
        if self.turn_direction=='clockwise':
            z=(self.Xp+1j*self.Yp)*exp(-1j*self.theta[ii])            
            xdata=z.real
            ydata=z.imag
            arc=arc*exp(-1j*self.theta[ii]) 
        elif self.turn_direction=='anti-clockwise':
            z=(self.Xp+1j*self.Yp)*exp(1j*self.theta[ii])            
            xdata=z.real
            ydata=z.imag
            arc=arc*exp(1j*self.theta[ii])
        # Updating Cam profile data
        dataxy=array([xdata,ydata]).transpose()
        self.poly1.set_xy(dataxy)
        
        # Updating center mark data
        arc_data=array([arc.real,arc.imag]).transpose()               
        self.artc.set_xy(arc_data)       
        
        # Updating Roller Follower data
        
        xpiv=self.xpiv
        ypiv=self.ypiv        
        qq=self.q[ii]*exp(1j*self.PivAngpos)
        xend=qq.real 
        yend=qq.imag
        xkna,ykna=self.OscRollerFollower(xpiv,ypiv,xend,yend)
        datafollower=array([xkna,ykna]).transpose()        
        self.poly2.set_xy(datafollower)
        
        # Updating  Roller Data
        
        Roller=self.Circle(self.Rroller,xend,yend)
        xdatRoll=Roller.real
        ydatRoll=Roller.imag
        datafollower=array([xdatRoll,ydatRoll]).transpose()
        self.rollpol.set_xy(datafollower)
            
        return self.poly1, self.poly2
