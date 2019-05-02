# -*- coding: utf-8 -*-

global config
#config='14out9'
#config='15out6'
#config='mobie'
config='15out6cam'

from numpy import arange, sin, pi
import matplotlib
from random import randint
from multiprocessing import Process, Queue
from threading import *
import time
from matplotlib.backends.backend_wxagg import NavigationToolbar2GTK3 as NavigationToolbar
matplotlib.use('WXAgg')
from scipy import optimize
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import numpy
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure
from scipy import optimize
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import numpy
import wx
import wx.grid as gridlib
import sys
import copy
from matplotlib.patches import Rectangle
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.anchored_artists import AnchoredAuxTransformBox
from mpl_toolkits.axes_grid.anchored_artists import AnchoredText
from mpl_toolkits.axes_grid.anchored_artists import AnchoredDrawingArea
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredEllipse
fig=plt.figure(1, figsize=(3,3))
ax = plt.subplot(111)
#params = (2, 3)

def on_xlims_change(axes):
    print "updated xlims: ", axes.get_xlim()

def on_ylims_change(axes):
    print "updated ylims: ", axes.get_ylim()



#optimization function 
def print_fun(x, f, accepted):
        print("at minima %.4f accepted %d" % (f, int(accepted)))
def f(p):
	global redo_x,redo_y,redo_coll_dx_x,redo_coll_dy_x,redo_coll_dz_x,redo_coll_dx_y,redo_coll_dy_y,redo_coll_dz_y,grat_tx_x,grat_tx_y,redo_x1,redo_y1,redo_x2,redo_y2
	global blueo_x,blueo_y,blueo_coll_dx_x,blueo_coll_dy_x,blueo_coll_dz_x,blueo_coll_dx_y,blueo_coll_dy_y,blueo_coll_dz_y,grat_tx_x,grat_tx_y,blueo_x1,blueo_y1,blueo_x2,blueo_y2
        global redo_x,redo_y,redo_mirr_dx_x,redo_mirr_dy_x,redo_mirr_dz_x,redo_mirr_dx_y,redo_mirr_dy_y,redo_mirr_dz_y,grat_tx_x,grat_tx_y
	global blueo_x,blueo_y,blueo_dchr_dx_x,blueo_dchr_dy_x,blueo_dchr_dz_x,blueo_dchr_dx_y,blueo_dchr_dy_y,blueo_dchr_dz_y,grat_tx_x,grat_tx_y
        global redo_x,redo_y,redo_gra_dx_x,redo_gra_dy_x,redo_gra_dz_x,redo_gra_dx_y,redo_gra_dy_y,redo_gra_dz_y,grat_tx_x,grat_tx_y
	global blueo_x,blueo_y,blueo_gra_dx_x,blueo_gra_dy_x,blueo_gra_dz_x,blueo_gra_dx_y,blueo_gra_dy_y,blueo_gra_dz_y,grat_tx_x,grat_tx_y
	global active,rxrot2,ryrot2,bxrot2,bxrot2,coll_comp,mirr_comp,dchr_comp,rgra_comp,bgra_comp
	
	#print p[0:6]
	redo_x2= copy.copy(redo_x1)
	redo_y2= copy.copy(redo_y1)
	
	act_coll_dx=p[0] 
	act_coll_dy=p[1] 
	act_coll_dz=p[2] 
	act_coll_tx=p[3] 
	act_coll_ty=p[4] 
	act_coll_tz=p[5] 
	#print act_coll_tx
	act_mirr_dx=p[6] 
	act_mirr_dy=p[7] 
	act_mirr_dz=p[8] 
	act_mirr_tx=p[9] 
	act_mirr_ty=p[10] 
	act_mirr_tz=p[11] 	

	act_rgra_dx=p[12] 
	act_rgra_dy=p[13] 
	act_rgra_dz=p[14] 
	act_rgra_tx=p[15] 
	act_rgra_ty=p[16] 
	act_rgra_tz=p[17] 
	
	act_dchr_dx=p[18] 
	act_dchr_dy=p[19] 
	act_dchr_dz=p[20] 
	act_dchr_tx=p[21] 
	act_dchr_ty=p[22] 
	act_dchr_tz=p[23] 
	
	act_bgra_dx=p[24] 
	act_bgra_dy=p[25] 
	act_bgra_dz=p[26] 
	act_bgra_tx=p[27] 
	act_bgra_ty=p[28] 
	act_bgra_tz=p[29] 
	det_red_rot=p[30]
	det_blue_rot=p[31]
	det_red_shiftx=p[32]
	det_red_shifty=p[33]
	det_blue_shiftx=p[34]
	det_blue_shifty=p[35]		
		
	
	if 0 in active:
	    if 0 in coll_comp:
		redo_x2= redo_x2+act_coll_dx*redo_coll_dx_x 
	    if 1 in coll_comp:
		redo_x2= redo_x2+act_coll_dy*redo_coll_dy_x
	    if 2 in coll_comp:
		redo_x2= redo_x2+act_coll_dz*redo_coll_dz_x
	    if 3 in coll_comp:
		redo_x2= redo_x2+act_coll_tx*redo_coll_tx_x
	    if 4 in coll_comp:
		redo_x2= redo_x2+act_coll_ty*redo_coll_ty_x
	    if 5 in coll_comp:
		redo_x2= redo_x2+act_coll_tz*redo_coll_tz_x										
#		+act_coll_dy*redo_coll_dy_x	+ act_coll_dz*redo_coll_dz_x + act_coll_tx*redo_coll_tx_x +act_coll_ty*redo_coll_ty_x+ act_coll_tz*redo_coll_tz_x 
	if 1 in active:
    	    if 0 in mirr_comp:
		redo_x2= redo_x2+act_mirr_dx*redo_mirr_dx_x 
	    if 1 in mirr_comp:
		redo_x2= redo_x2+act_mirr_dy*redo_mirr_dy_x
	    if 2 in mirr_comp:
		redo_x2= redo_x2+act_mirr_dz*redo_mirr_dz_x
	    if 3 in mirr_comp:
		redo_x2= redo_x2+act_mirr_tx*redo_mirr_tx_x
	    if 4 in mirr_comp:
		redo_x2= redo_x2+act_mirr_ty*redo_mirr_ty_x
	    if 5 in mirr_comp:
		redo_x2= redo_x2+act_mirr_tz*redo_mirr_tz_x										

#	    redo_x2= redo_x2+act_mirr_dx*redo_mirr_dx_x +act_mirr_dy*redo_mirr_dy_x	+ act_mirr_dz*redo_mirr_dz_x + act_mirr_tx*redo_mirr_tx_x +act_mirr_ty*redo_mirr_ty_x+ act_mirr_tz*redo_mirr_tz_x 
	if 3 in active:
	    redo_x2= redo_x2+act_rgra_dx*redo_gra_dx_x +act_rgra_dy*redo_gra_dy_x	+ act_rgra_dz*redo_gra_dz_x + act_rgra_tx*redo_gra_tx_x +act_rgra_ty*redo_gra_ty_x+ act_rgra_tz*redo_gra_tz_x
	if 0 in active:
	    if 0 in coll_comp:
		redo_y2= redo_y2+act_coll_dx*redo_coll_dx_y 
	    if 1 in coll_comp:
		redo_y2= redo_y2+act_coll_dy*redo_coll_dy_y
	    if 2 in coll_comp:
		redo_y2= redo_y2+act_coll_dz*redo_coll_dz_y
	    if 3 in coll_comp:
		redo_y2= redo_y2+act_coll_tx*redo_coll_tx_y
	    if 4 in coll_comp:
		redo_y2= redo_y2+act_coll_ty*redo_coll_ty_y
	    if 5 in coll_comp:
		redo_y2= redo_y2+act_coll_tz*redo_coll_tz_y										
	    
#	    redo_y2= redo_y2+act_coll_dx*redo_coll_dx_y +act_coll_dy*redo_coll_dy_y	+ act_coll_dz*redo_coll_dz_y + act_coll_tx*redo_coll_tx_y +act_coll_ty*redo_coll_ty_y+ act_coll_tz*redo_coll_tz_y 
	if 1 in active:
	    if 0 in mirr_comp:
		redo_y2= redo_y2+act_mirr_dx*redo_mirr_dx_y 
	    if 1 in mirr_comp:
		redo_y2= redo_y2+act_mirr_dy*redo_mirr_dy_y
	    if 2 in mirr_comp:
		redo_y2= redo_y2+act_mirr_dz*redo_mirr_dz_y
	    if 3 in mirr_comp:
		redo_y2= redo_y2+act_mirr_tx*redo_mirr_tx_y
	    if 4 in mirr_comp:
		redo_y2= redo_y2+act_mirr_ty*redo_mirr_ty_y
	    if 5 in mirr_comp:
		redo_y2= redo_y2+act_mirr_tz*redo_mirr_tz_y	
		
#	    redo_y2= redo_y2+act_mirr_dx*redo_mirr_dx_y +act_mirr_dy*redo_mirr_dy_y	+ act_mirr_dz*redo_mirr_dz_y + act_mirr_tx*redo_mirr_tx_y +act_mirr_ty*redo_mirr_ty_y+ act_mirr_tz*redo_mirr_tz_y 
	if 3 in active:
	    redo_y2= redo_y2+act_rgra_dx*redo_gra_dx_y +act_rgra_dy*redo_gra_dy_y	+ act_rgra_dz*redo_gra_dz_y + act_rgra_tx*redo_gra_tx_y +act_rgra_ty*redo_gra_ty_y+ act_rgra_tz*redo_gra_tz_y
	if 5 in active:
	    rxrot2 = redo_x2*np.cos(det_red_rot)-redo_y2*np.sin(det_red_rot)
            ryrot2=  redo_x2*np.sin(det_red_rot)+redo_y2*np.cos(det_red_rot)
	    redo_x2=copy.copy(rxrot2)
	    redo_y2=copy.copy(ryrot2) 

	if 6 in active:
	    redo_x2=redo_x2+det_red_shiftx
	    redo_y2=redo_y2+det_red_shifty
	          
        blueo_x2= copy.copy(blueo_x1)
	blueo_y2= copy.copy(blueo_y1)
	
	
	if 0 in active:
	    if 0 in coll_comp:
		blueo_x2= blueo_x2+act_coll_dx*blueo_coll_dx_x 
	    if 1 in coll_comp:
		blueo_x2= blueo_x2+act_coll_dy*blueo_coll_dy_x
	    if 2 in coll_comp:
		blueo_x2= blueo_x2+act_coll_dz*blueo_coll_dz_x
	    if 3 in coll_comp:
		blueo_x2= blueo_x2+act_coll_tx*blueo_coll_tx_x
	    if 4 in coll_comp:
		blueo_x2= blueo_x2+act_coll_ty*blueo_coll_ty_x
	    if 5 in coll_comp:
		blueo_x2= blueo_x2+act_coll_tz*blueo_coll_tz_x										
	    
#	    blueo_x2= blueo_x2+act_coll_dx*blueo_coll_dx_x +act_coll_dy*blueo_coll_dy_x	+ act_coll_dz*blueo_coll_dz_x + act_coll_tx*blueo_coll_tx_x +act_coll_ty*blueo_coll_ty_x+ act_coll_tz*blueo_coll_tz_x 
	if 2 in active:
	    blueo_x2= blueo_x2+act_dchr_dx*blueo_dchr_dx_x +act_dchr_dy*blueo_dchr_dy_x	+ act_dchr_dz*blueo_dchr_dz_x + act_dchr_tx*blueo_dchr_tx_x +act_dchr_ty*blueo_dchr_ty_x+ act_dchr_tz*blueo_dchr_tz_x 
	if 4 in active:
	    blueo_x2= blueo_x2+act_bgra_dx*blueo_gra_dx_x +act_bgra_dy*blueo_gra_dy_x	+ act_bgra_dz*blueo_gra_dz_x + act_bgra_tx*blueo_gra_tx_x +act_bgra_ty*blueo_gra_ty_x+ act_bgra_tz*blueo_gra_tz_x
	if 0 in active:
	    if 0 in coll_comp:
		blueo_y2= blueo_y2+act_coll_dx*blueo_coll_dx_y 
	    if 1 in coll_comp:
		blueo_y2= blueo_y2+act_coll_dy*blueo_coll_dy_y
	    if 2 in coll_comp:
		blueo_y2= blueo_y2+act_coll_dz*blueo_coll_dz_y
	    if 3 in coll_comp:
		blueo_y2= blueo_y2+act_coll_tx*blueo_coll_tx_y
	    if 4 in coll_comp:
		blueo_y2= blueo_y2+act_coll_ty*blueo_coll_ty_y
	    if 5 in coll_comp:
		blueo_y2= blueo_y2+act_coll_tz*blueo_coll_tz_y	    
#	    blueo_y2= blueo_y2+act_coll_dx*blueo_coll_dx_y +act_coll_dy*blueo_coll_dy_y	+ act_coll_dz*blueo_coll_dz_y + act_coll_tx*blueo_coll_tx_y +act_coll_ty*blueo_coll_ty_y+ act_coll_tz*blueo_coll_tz_y 
	if 2 in active:
	    blueo_y2= blueo_y2+act_dchr_dx*blueo_dchr_dx_y +act_dchr_dy*blueo_dchr_dy_y	+ act_dchr_dz*blueo_dchr_dz_y + act_dchr_tx*blueo_dchr_tx_y +act_dchr_ty*blueo_dchr_ty_y+ act_dchr_tz*blueo_dchr_tz_y 
	if 4 in active:
	    blueo_y2= blueo_y2+act_bgra_dx*blueo_gra_dx_y +act_bgra_dy*blueo_gra_dy_y	+ act_bgra_dz*blueo_gra_dz_y + act_bgra_tx*blueo_gra_tx_y +act_bgra_ty*blueo_gra_ty_y+ act_bgra_tz*blueo_gra_tz_y
	if 5 in active:	
	    bxrot2 = blueo_x2*np.cos(det_blue_rot)-blueo_y2*np.sin(det_blue_rot)
            byrot2=  blueo_x2*np.sin(det_blue_rot)+blueo_y2*np.cos(det_blue_rot)
	    blueo_x2=copy.copy(bxrot2)
	    blueo_y2=copy.copy(byrot2)
	if 6 in active:	
	    blueo_x2=blueo_x2+det_blue_shiftx
	    blueo_y2=blueo_y2+det_blue_shifty        
	
	#print 'hello',det_blue_shiftx
	err=(sum(np.sqrt((redo_x-redo_x2)**2+(redo_y-redo_y2)**2)) + sum(np.sqrt((blueo_x-blueo_x2)**2+(blueo_y-blueo_y2)**2)))#/100#/(len(redo_y2)+len(blueo_y2))
	#print  sum(np.sqrt((redo_x-redo_x2)**2+(redo_y-redo_y2)**2))/len(redo_y2)
	#print  sum(np.sqrt((blueo_x-blueo_x2)**2+(blueo_y-blueo_y2)**2))/len(blueo_y2)
	wx.Yield()
	#a,b = params
        #if err<1:
	    #print redo_x2[1],redo_x[1],redo_y2[1],redo_y[1]
	#print err
	return err

#--StdOut Class--
class MyTakeStep(object):
   def __init__(self, stepsize=0.5):
       self.stepsize = stepsize
   def __call__(self, x):
       self.pp=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
       self.l=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
       self.u=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
       x=copy.copy(self.pp)
       if 0 in active :
	    self.pp[0]=500
	    self.pp[1]=500
	    self.pp[2]=500
	    self.pp[3]=500
	    self.pp[4]=500
	    self.pp[5]=500
	    self.l[0]=-500
	    self.l[1]=-500
	    self.l[2]=-500
	    self.l[3]=-500
	    self.l[4]=-500
	    self.l[5]=-500
	    self.u[0]=500
	    self.u[1]=500
	    self.u[2]=500
	    self.u[3]=500
	    self.u[4]=500
	    self.u[5]=500
       if 1 in active :
	    self.pp[6]=500
	    self.pp[7]=500
	    self.pp[8]=500
	    self.pp[9]=500
	    self.pp[10]=500
	    self.pp[11]=500
	    self.l[6]=-500
	    self.l[7]=-500
	    self.l[8]=-500
	    self.l[9]=-500
	    self.l[10]=-500
	    self.l[11]=-500
	    self.u[6]=500
	    self.u[7]=500
	    self.u[8]=500
	    self.u[9]=500
	    self.u[10]=500
	    self.u[11]=500
       if 3 in active :
	    self.pp[12]=500
	    self.pp[13]=500
	    self.pp[14]=500
	    self.pp[15]=500
	    self.pp[16]=500
	    self.pp[17]=500
	    self.l[12]=-500
	    self.l[13]=-500
	    self.l[14]=-500
	    self.l[15]=-500
	    self.l[16]=-500
	    self.l[17]=-500
	    self.u[12]=500
	    self.u[13]=500
	    self.u[14]=500
	    self.u[15]=500
	    self.u[16]=500
	    self.u[17]=500
       if 2 in active :
	    self.pp[18]=500
	    self.pp[19]=500
	    self.pp[20]=500
	    self.pp[21]=500
	    self.pp[22]=500
	    self.pp[23]=500
	    self.l[18]=-500
	    self.l[19]=-500
	    self.l[20]=-500
	    self.l[21]=-500
	    self.l[22]=-500
	    self.l[23]=-500
	    self.u[18]=500
	    self.u[19]=500
	    self.u[20]=500
	    self.u[21]=500
	    self.u[22]=500
	    self.u[23]=500
       if 4 in active :
	    self.pp[24]=500
	    self.pp[25]=500
	    self.pp[26]=500
	    self.pp[27]=500
	    self.pp[28]=500
	    self.pp[29]=500
	    self.l[24]=-500
	    self.l[25]=-500
	    self.l[26]=-500
	    self.l[27]=-500
	    self.l[28]=-500
	    self.l[29]=-500
	    self.u[24]=500
	    self.u[25]=500
	    self.u[26]=500
	    self.u[27]=500
	    self.u[28]=500
	    self.u[29]=500
       if 5 in active :
	    self.pp[30]=.001
	    self.pp[31]=.001
	    self.u[30]=3.14
	    self.u[31]=3.14
	    self.l[30]=-3.14
	    self.l[31]=-3.14
       if 6 in active :
	    self.pp[32]=.001
	    self.pp[33]=.001
	    self.pp[34]=.001
	    self.pp[35]=.001
    
       if 0 in active :
	   s = self.pp[0]
	   if 0 in coll_comp:
	       x[0] += np.random.uniform(-2.*s, 2.*s)
	   if 1 in coll_comp:
	       x[1] += np.random.uniform(-2.*s, 2.*s)	       
	   if 2 in coll_comp:
	       x[2] += np.random.uniform(-2.*s, 2.*s)
	   if 3 in coll_comp:
	       x[3] += np.random.uniform(-2.*s, 2.*s)
	   if 4 in coll_comp:
	       x[4] += np.random.uniform(-2.*s, 2.*s)
	   if 5 in coll_comp:
	       x[5] += np.random.uniform(-2.*s, 2.*s)	       	       	       	    
       #print x[0:6] 
       if 1 in active :
	   s = self.pp[6]
	   if 0 in mirr_comp:
	       x[6] += np.random.uniform(-2.*s, 2.*s)
	   if 1 in mirr_comp:
	       x[7] += np.random.uniform(-2.*s, 2.*s)	       
	   if 2 in mirr_comp:
	       x[8] += np.random.uniform(-2.*s, 2.*s)
	   if 3 in mirr_comp:
	       x[9] += np.random.uniform(-2.*s, 2.*s)
	   if 4 in mirr_comp:
	       x[10] += np.random.uniform(-2.*s, 2.*s)
	   if 5 in mirr_comp:
	       x[11] += np.random.uniform(-2.*s, 2.*s)	  	    
#        x[6:12] += np.random.uniform(-2.*s, 2.*s,6)
       if 3 in active :
	    s = self.pp[12]
            x[12:18] += np.random.uniform(-2.*s, 2.*s,6)
       if 2 in active :
	    s = self.pp[18]
            x[18:24] += np.random.uniform(-2.*s, 2.*s,6)
       if 4 in active :
	    s = self.pp[24]
            x[24:30] += np.random.uniform(-2.*s, 2.*s,6)
       if 5 in active :
	    s = .5
            x[30:32] += np.random.uniform(-s, s,2)
       if 6 in active :
	    s = .1
            x[32:36] += np.random.uniform(-s, s,4)	    
       return x


class RedirectText(object):
    def __init__(self,aWxTextCtrl):
        self.out=aWxTextCtrl
 
    def write(self,string):
        self.out.WriteText(string)

#Matplotlib Class
class CanvasPanel(wx.Panel):
    def __init__(self, parent):
	global red_x,red_y,blue_x,blue_y
        wx.Panel.__init__(self, parent)
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.toolbar = NavigationToolbar(self.canvas)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
	self.sizer.Add(self.toolbar, 0, wx.EXPAND)
        self.SetSizer(self.sizer)
        self.Fit()
	self.flexscale=1000
	self.resscale=10000	
	self.pathscale=1000		

    def drawblue(self):
	global red_x,red_y,blueo_x,blueo_y
	self.axes.clear()
	self.axes.scatter(blueo_x, blueo_y,color='blue',s=7)

        self.axes.set_ylabel("Dispersion-Direction mm", fontsize=8)
        self.axes.set_xlabel("Slit-Direction mm", fontsize=8)
	self.axes.set_title("Blue FCS Field", fontsize=8)
	self.axes.axis([-150,150,-150,150])
	self.canvas.draw()


    def drawanalyzeblue(self):
        global blue_coll_dx_x_flex,blue_coll_dy_x_flex,blue_coll_dz_x_flex,blue_coll_tx_x_flex,blue_coll_ty_x_flex,blue_coll_tz_x_flex,blue_gra_dx_x_flex,blue_gra_dy_x_flex,blue_gra_dz_x_flex,blue_gra_tx_x_flex,blue_gra_ty_x_flex,blue_gra_tz_x_flex, blue_dchr_dx_x_flex,blue_dchr_dy_x_flex,blue_dchr_dz_x_flex,blue_dchr_tx_x_flex,blue_dchr_ty_x_flex,blue_dchr_tz_x_flex
        global blue_coll_dx_y_flex,blue_coll_dy_y_flex,blue_coll_dz_y_flex,blue_coll_tx_y_flex,blue_coll_ty_y_flex,blue_coll_tz_y_flex,blue_gra_dx_y_flex,blue_gra_dy_y_flex,blue_gra_dz_y_flex,blue_gra_tx_y_flex,blue_gra_ty_y_flex,blue_gra_tz_y_flex, blue_dchr_dx_y_flex,blue_dchr_dy_y_flex,blue_dchr_dz_y_flex,blue_dchr_tx_y_flex,blue_dchr_ty_y_flex,blue_dchr_tz_y_flex
        global blue_coll_dx_x_comp,blue_coll_dy_x_comp,blue_coll_dz_x_comp,blue_coll_tx_x_comp,blue_coll_ty_x_comp,blue_coll_tz_x_comp,blue_det_x_comp,blue_det_y_comp, blue_dchr_dx_x_comp,blue_dchr_dy_x_comp,blue_dchr_dz_x_comp,blue_dchr_tx_x_comp,blue_dchr_ty_x_comp,blue_dchr_tz_x_comp
        global blue_coll_dx_y_comp,blue_coll_dy_y_comp,blue_coll_dz_y_comp,blue_coll_tx_y_comp,blue_coll_ty_y_comp,blue_coll_tz_y_comp,blue_dchr_dx_y_comp,blue_dchr_dy_y_comp,blue_dchr_dz_y_comp,blue_dchr_tx_y_comp,blue_dchr_ty_y_comp,blue_dchr_tz_y_comp
        global active,flexboxval,compboxval,resboxval, red_flex_x,red_flex_y,red_comp_x,red_comp_y,blue_flex_x,blue_flex_y,blue_comp_x,blue_comp_y
	self.axes.clear()
        arrx=[blue_flex_x,blue_comp_x,blue_coll_dx_x_flex,blue_coll_dy_x_flex,blue_coll_dz_x_flex,blue_coll_tx_x_flex,blue_coll_ty_x_flex,blue_coll_tz_x_flex,blue_gra_dx_x_flex,blue_gra_dy_x_flex,blue_gra_dz_x_flex,blue_gra_tx_x_flex,blue_gra_ty_x_flex,blue_gra_tz_x_flex, blue_dchr_dx_x_flex,blue_dchr_dy_x_flex,blue_dchr_dz_x_flex,blue_dchr_tx_x_flex,blue_dchr_ty_x_flex,blue_dchr_tz_x_flex,blue_coll_dx_x_comp,blue_coll_dy_x_comp,blue_coll_dz_x_comp,blue_coll_tx_x_comp,blue_coll_ty_x_comp,blue_coll_tz_x_comp,blue_det_x_comp,blue_dchr_dx_x_comp,blue_dchr_dy_x_comp,blue_dchr_dz_x_comp,blue_dchr_tx_x_comp,blue_dchr_ty_x_comp,blue_dchr_tz_x_comp]
        arry=[blue_flex_y,blue_comp_y,blue_coll_dx_y_flex,blue_coll_dy_y_flex,blue_coll_dz_y_flex,blue_coll_tx_y_flex,blue_coll_ty_y_flex,blue_coll_tz_y_flex,blue_gra_dx_y_flex,blue_gra_dy_y_flex,blue_gra_dz_y_flex,blue_gra_tx_y_flex,blue_gra_ty_y_flex,blue_gra_tz_y_flex, blue_dchr_dx_y_flex,blue_dchr_dy_y_flex,blue_dchr_dz_y_flex,blue_dchr_tx_y_flex,blue_dchr_ty_y_flex,blue_dchr_tz_y_flex,blue_coll_dx_y_comp,blue_coll_dy_y_comp,blue_coll_dz_y_comp,blue_coll_tx_y_comp,blue_coll_ty_y_comp,blue_coll_tz_y_comp,blue_det_y_comp,blue_dchr_dx_y_comp,blue_dchr_dy_y_comp,blue_dchr_dz_y_comp,blue_dchr_tx_y_comp,blue_dchr_ty_y_comp,blue_dchr_tz_y_comp]
	arrxmax=np.abs(np.array(arrx)).max()
	arrymax=np.abs(np.array(arry)).max()
	arrxmax=arrxmax+.1*arrxmax
	arrymax=arrymax+.1*arrymax
        maxx=np.array([arrxmax,arrymax]).max()
        if flexboxval==True:
	    
	    self.axes.arrow(0,0,blue_coll_dx_x_flex,blue_coll_dx_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    self.axes.text(blue_coll_dx_x_flex+.025*maxx,blue_coll_dx_y_flex,'coll_dx',fontsize=12,color = 'red')
	    self.axes.arrow(0,0,blue_coll_dy_x_flex,blue_coll_dy_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    self.axes.text(blue_coll_dy_x_flex+.025*maxx,blue_coll_dy_y_flex,'coll_dy',fontsize=12,color='red')
	    #self.axes.arrow(0,0,blue_coll_dz_x_flex,blue_coll_dz_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color='green',alpha=.5)
	    #self.axes.text(blue_coll_dz_x_flex,blue_coll_dz_y_flex,'coll_dz',fontsize=12,color='green')
	    self.axes.arrow(0,0,blue_coll_tx_x_flex,blue_coll_tx_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    self.axes.text(blue_coll_tx_x_flex+.025*maxx,blue_coll_tx_y_flex,'coll_tx',fontsize=12,color='red')
	    self.axes.arrow(0,0,blue_coll_ty_x_flex,blue_coll_ty_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    self.axes.text(blue_coll_ty_x_flex,blue_coll_ty_y_flex+.025*maxx,'coll_ty',fontsize=12,color='red')
	    #self.axes.arrow(0,0,blue_coll_tz_x_flex,blue_coll_tz_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'purple',alpha=.5)
	    #self.axes.text(blue_coll_tz_x_flex,blue_coll_tz_y_flex,'coll_tz',fontsize=12,color='purple')
	    
	    #self.axes.arrow(0,0,blue_gra_dx_x_flex,blue_gra_dx_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    #self.axes.text(blue_gra_dx_x_flex,blue_gra_dx_y_flex,'gra_dx',fontsize=12,color = 'red')
	    #self.axes.arrow(0,0,blue_gra_dy_x_flex,blue_gra_dy_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    #self.axes.text(blue_gra_dy_x_flex,blue_gra_dy_y_flex,'gra_dy',fontsize=12,color='red')
	    #self.axes.arrow(0,0,blue_gra_dz_x_flex,blue_gra_dz_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color='green',alpha=.5)
	    #self.axes.text(blue_gra_dz_x_flex,blue_gra_dz_y_flex,'gra_dz',fontsize=12,color='green')
	    self.axes.arrow(0,0,blue_gra_tx_x_flex,blue_gra_tx_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    self.axes.text(blue_gra_tx_x_flex+.025*maxx,blue_gra_tx_y_flex,'gra_tx',fontsize=12,color='red')
	    self.axes.arrow(0,0,blue_gra_ty_x_flex,blue_gra_ty_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    self.axes.text(blue_gra_ty_x_flex,blue_gra_ty_y_flex+.025*maxx,'gra_ty',fontsize=12,color='red')
	    self.axes.arrow(0,0,blue_gra_tz_x_flex,blue_gra_tz_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    self.axes.text(blue_gra_tz_x_flex,blue_gra_tz_y_flex+.025*maxx,'gra_tz',fontsize=12,color='red')
	    
	    #self.axes.arrow(0,0,blue_dchr_dx_x_flex,blue_dchr_dx_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    #self.axes.text(blue_dchr_dx_x_flex,blue_dchr_dx_y_flex,'dchr_dx',fontsize=12,color = 'red')
	    #self.axes.arrow(0,0,blue_dchr_dy_x_flex,blue_dchr_dy_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    #self.axes.text(blue_dchr_dy_x_flex,blue_dchr_dy_y_flex,'dchr_dy',fontsize=12,color='red')
	    #self.axes.arrow(0,0,blue_dchr_dz_x_flex,blue_dchr_dz_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color='green',alpha=.5)
	    #self.axes.text(blue_dchr_dz_x_flex,blue_dchr_dz_y_flex,'dchr_dz',fontsize=12,color='green')
	    self.axes.arrow(0,0,blue_dchr_tx_x_flex,blue_dchr_tx_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    self.axes.text(blue_dchr_tx_x_flex+.025*maxx,blue_dchr_tx_y_flex,'dchr_tx',fontsize=12,color='red')
	    self.axes.arrow(0,0,blue_dchr_ty_x_flex,blue_dchr_ty_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    self.axes.text(blue_dchr_ty_x_flex,blue_dchr_ty_y_flex+.025*maxx,'dchr_ty',fontsize=12,color='red')
	    #self.axes.arrow(0,0,blue_dchr_tz_x_flex,blue_dchr_tz_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'purple',alpha=.5)
	    #self.axes.text(blue_dchr_tz_x_flex,blue_dchr_tz_y_flex,'dchr_tz',fontsize=12,color='purple')
	    if resboxval==True:
		self.axes.arrow(0,0,blue_flex_x,blue_flex_y,head_width=.05*maxx, head_length=.1*maxx,width=(.005*maxx), color = 'red',alpha=.5)
		self.axes.text(blue_flex_x,blue_flex_y+.025*maxx,'Flex_Total',fontsize=12,color='red')
		
        if compboxval==True:
	    
	    if 0 in active:
		self.axes.arrow(0,0,blue_coll_dx_x_comp,blue_coll_dx_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'green',alpha=.5)
		self.axes.text(blue_coll_dx_x_comp,blue_coll_dx_y_comp+.025*maxx,'coll_dx',fontsize=12,color = 'green')
		self.axes.arrow(0,0,blue_coll_dy_x_comp,blue_coll_dy_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'green',alpha=.5)
		self.axes.text(blue_coll_dy_x_comp+.025*maxx,blue_coll_dy_y_comp,'coll_dy',fontsize=12,color='green')
	    #self.axes.arrow(0,0,blue_coll_dz_x_comp,blue_coll_dz_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color='green',alpha=.5)
	    #self.axes.text(blue_coll_dz_x_comp,blue_coll_dz_y_comp,'coll_dz',fontsize=12,color='green')
		self.axes.arrow(0,0,blue_coll_tx_x_comp,blue_coll_tx_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'green',alpha=.5)
		self.axes.text(blue_coll_tx_x_comp+.025*maxx,blue_coll_tx_y_comp,'coll_tx',fontsize=12,color='green')
		self.axes.arrow(0,0,blue_coll_ty_x_comp,blue_coll_ty_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'green',alpha=.5)
		self.axes.text(blue_coll_ty_x_comp,blue_coll_ty_y_comp+.025*maxx,'coll_ty',fontsize=12,color='green')
	    #self.axes.arrow(0,0,blue_coll_tz_x_comp,blue_coll_tz_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'purple',alpha=.5)
	    #self.axes.text(blue_coll_tz_x_comp,blue_coll_tz_y_comp,'coll_tz',fontsize=12,color='purple')
	    
	    if 6 in active:
		self.axes.arrow(0,0,blue_det_x_comp,blue_det_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'green',alpha=.5)
		self.axes.text(blue_det_x_comp+.025*maxx,blue_det_y_comp+.025*maxx,'det_shft',fontsize=12,color='green')
    
    
	    
	    #self.axes.arrow(0,0,blue_dchr_dx_x_comp,blue_dchr_dx_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'green',alpha=.5)
	    #self.axes.text(blue_dchr_dx_x_comp,blue_dchr_dx_y_comp,'dchr_dx',fontsize=12,color = 'green')
	    #self.axes.arrow(0,0,blue_dchr_dy_x_comp,blue_dchr_dy_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    #self.axes.text(blue_dchr_dy_x_comp,blue_dchr_dy_y_comp,'dchr_dy',fontsize=12,color='red')
	    #self.axes.arrow(0,0,blue_dchr_dz_x_comp,blue_dchr_dz_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color='green',alpha=.5)
	    #self.axes.text(blue_dchr_dz_x_comp,blue_dchr_dz_y_comp,'dchr_dz',fontsize=12,color='green')
	    if 2 in active:
		self.axes.arrow(0,0,blue_dchr_tx_x_comp,blue_dchr_tx_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'green',alpha=.5)
		self.axes.text(blue_dchr_tx_x_comp+.025*maxx,blue_dchr_tx_y_comp,'dchr_tx',fontsize=12,color='green')
		self.axes.arrow(0,0,blue_dchr_ty_x_comp,blue_dchr_ty_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'green',alpha=.5)
		self.axes.text(blue_dchr_ty_x_comp,blue_dchr_ty_y_comp+.025*maxx,'dchr_ty',fontsize=12,color='green')
	    #self.axes.arrow(0,0,blue_dchr_tz_x_comp,blue_dchr_tz_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'purple',alpha=.5)
	    #self.axes.text(blue_dchr_tz_x_comp,blue_dchr_tz_y_comp,'dchr_tz',fontsize=12,color='purple')

	    if resboxval==True:
		self.axes.arrow(0,0,blue_comp_x,blue_comp_y,head_width=.05*maxx, head_length=.1*maxx,width=(.005*maxx), color = 'green',alpha=.5)
		self.axes.text(blue_comp_x,blue_comp_y+.025*maxx,'Comp_Total',fontsize=12,color='green')
	
	
	
	
	
	
	
	
        self.axes.set_ylabel("Shift on Detector (mm)", fontsize=12)
        self.axes.set_xlabel("Shift on Detector (mm)", fontsize=12)
	self.axes.axhline(0, color='black',linestyle='--')
        self.axes.axvline(0, color='black',linestyle='--')
	self.axes.axis([-maxx,maxx,-maxx,maxx])
	self.canvas.draw()


    def drawanalyzered(self):
        global red_coll_dx_x_flex,red_coll_dy_x_flex,red_coll_dz_x_flex,red_coll_tx_x_flex,red_coll_ty_x_flex,red_coll_tz_x_flex,red_gra_dx_x_flex,red_gra_dy_x_flex,red_gra_dz_x_flex,red_gra_tx_x_flex,red_gra_ty_x_flex,red_gra_tz_x_flex, red_mirr_dx_x_flex,red_mirr_dy_x_flex,red_mirr_dz_x_flex,red_mirr_tx_x_flex,red_mirr_ty_x_flex,red_mirr_tz_x_flex
        global red_coll_dx_y_flex,red_coll_dy_y_flex,red_coll_dz_y_flex,red_coll_tx_y_flex,red_coll_ty_y_flex,red_coll_tz_y_flex,red_gra_dx_y_flex,red_gra_dy_y_flex,red_gra_dz_y_flex,red_gra_tx_y_flex,red_gra_ty_y_flex,red_gra_tz_y_flex, red_mirr_dx_y_flex,red_mirr_dy_y_flex,red_mirr_dz_y_flex,red_mirr_tx_y_flex,red_mirr_ty_y_flex,red_mirr_tz_y_flex
        global red_coll_dx_x_comp,red_coll_dy_x_comp,red_coll_dz_x_comp,red_coll_tx_x_comp,red_coll_ty_x_comp,red_coll_tz_x_comp,red_det_x_comp,red_det_y_comp, red_mirr_dx_x_comp,red_mirr_dy_x_comp,red_mirr_dz_x_comp,red_mirr_tx_x_comp,red_mirr_ty_x_comp,red_mirr_tz_x_comp
        global red_coll_dx_y_comp,red_coll_dy_y_comp,red_coll_dz_y_comp,red_coll_tx_y_comp,red_coll_ty_y_comp,red_coll_tz_y_comp,red_mirr_dx_y_comp,red_mirr_dy_y_comp,red_mirr_dz_y_comp,red_mirr_tx_y_comp,red_mirr_ty_y_comp,red_mirr_tz_y_comp
        global active,flexboxval,compboxval,resboxval, red_flex_x,red_flex_y,red_comp_x,red_comp_y,blue_flex_x,blue_flex_y,blue_comp_x,blue_comp_y

	self.axes.clear()
        arrx=[red_flex_x,red_comp_x,red_coll_dx_x_flex,red_coll_dy_x_flex,red_coll_dz_x_flex,red_coll_tx_x_flex,red_coll_ty_x_flex,red_coll_tz_x_flex,red_gra_dx_x_flex,red_gra_dy_x_flex,red_gra_dz_x_flex,red_gra_tx_x_flex,red_gra_ty_x_flex,red_gra_tz_x_flex, red_mirr_dx_x_flex,red_mirr_dy_x_flex,red_mirr_dz_x_flex,red_mirr_tx_x_flex,red_mirr_ty_x_flex,red_mirr_tz_x_flex,red_coll_dx_x_comp,red_coll_dy_x_comp,red_coll_dz_x_comp,red_coll_tx_x_comp,red_coll_ty_x_comp,red_coll_tz_x_comp,red_det_x_comp,red_mirr_dx_x_comp,red_mirr_dy_x_comp,red_mirr_dz_x_comp,red_mirr_tx_x_comp,red_mirr_ty_x_comp,red_mirr_tz_x_comp]
        arry=[red_flex_y,red_comp_y,red_coll_dx_y_flex,red_coll_dy_y_flex,red_coll_dz_y_flex,red_coll_tx_y_flex,red_coll_ty_y_flex,red_coll_tz_y_flex,red_gra_dx_y_flex,red_gra_dy_y_flex,red_gra_dz_y_flex,red_gra_tx_y_flex,red_gra_ty_y_flex,red_gra_tz_y_flex, red_mirr_dx_y_flex,red_mirr_dy_y_flex,red_mirr_dz_y_flex,red_mirr_tx_y_flex,red_mirr_ty_y_flex,red_mirr_tz_y_flex,red_coll_dx_y_comp,red_coll_dy_y_comp,red_coll_dz_y_comp,red_coll_tx_y_comp,red_coll_ty_y_comp,red_coll_tz_y_comp,red_det_y_comp,red_mirr_dx_y_comp,red_mirr_dy_y_comp,red_mirr_dz_y_comp,red_mirr_tx_y_comp,red_mirr_ty_y_comp,red_mirr_tz_y_comp]
	arrxmax=np.abs(np.array(arrx)).max()
	arrymax=np.abs(np.array(arry)).max()
	arrxmax=arrxmax+.1*arrxmax
	arrymax=arrymax+.1*arrymax
        maxx=np.array([arrxmax,arrymax]).max()
	
	if flexboxval==True:
		
	    self.axes.arrow(0,0,red_coll_dx_x_flex,red_coll_dx_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    self.axes.text(red_coll_dx_x_flex,red_coll_dx_y_flex+.025*maxx,'coll_dx',fontsize=12,color = 'red')
	    self.axes.arrow(0,0,red_coll_dy_x_flex,red_coll_dy_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    self.axes.text(red_coll_dy_x_flex+.025*maxx,red_coll_dy_y_flex,'coll_dy',fontsize=12,color='red')
	    #self.axes.arrow(0,0,red_coll_dz_x_flex,red_coll_dz_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color='green',alpha=.5)
	    #self.axes.text(red_coll_dz_x_flex,red_coll_dz_y_flex,'coll_dz',fontsize=12,color='green')
	    self.axes.arrow(0,0,red_coll_tx_x_flex,red_coll_tx_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    self.axes.text(red_coll_tx_x_flex+.025*maxx,red_coll_tx_y_flex,'coll_tx',fontsize=12,color='red')
	    self.axes.arrow(0,0,red_coll_ty_x_flex,red_coll_ty_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    self.axes.text(red_coll_ty_x_flex,red_coll_ty_y_flex+.025*maxx,'coll_ty',fontsize=12,color='red')
	    #self.axes.arrow(0,0,red_coll_tz_x_flex,red_coll_tz_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'purple',alpha=.5)
	    #self.axes.text(red_coll_tz_x_flex,red_coll_tz_y_flex,'coll_tz',fontsize=12,color='purple')
	    
	    #self.axes.arrow(0,0,red_gra_dx_x_flex,red_gra_dx_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    #self.axes.text(red_gra_dx_x_flex,red_gra_dx_y_flex,'gra_dx',fontsize=12,color = 'red')
	    #self.axes.arrow(0,0,red_gra_dy_x_flex,red_gra_dy_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'blue',alpha=.5)
	    #self.axes.text(red_gra_dy_x_flex,red_gra_dy_y_flex,'gra_dy',fontsize=12,color='blue')
	    #self.axes.arrow(0,0,red_gra_dz_x_flex,red_gra_dz_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color='green',alpha=.5)
	    #self.axes.text(red_gra_dz_x_flex,red_gra_dz_y_flex,'gra_dz',fontsize=12,color='green')
	    self.axes.arrow(0,0,red_gra_tx_x_flex,red_gra_tx_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    self.axes.text(red_gra_tx_x_flex+.025*maxx,red_gra_tx_y_flex,'gra_tx',fontsize=12,color='red')
	    self.axes.arrow(0,0,red_gra_ty_x_flex,red_gra_ty_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    self.axes.text(red_gra_ty_x_flex,red_gra_ty_y_flex+.025*maxx,'gra_ty',fontsize=12,color='red')
	    self.axes.arrow(0,0,red_gra_tz_x_flex,red_gra_tz_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    self.axes.text(red_gra_tz_x_flex,red_gra_tz_y_flex+.025*maxx,'gra_tz',fontsize=12,color='red')
	    
	    #self.axes.arrow(0,0,red_mirr_dx_x_flex,red_mirr_dx_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    #self.axes.text(red_mirr_dx_x_flex,red_mirr_dx_y_flex,'mirr_dx',fontsize=12,color = 'red')
	    #self.axes.arrow(0,0,red_mirr_dy_x_flex,red_mirr_dy_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'blue',alpha=.5)
	    #self.axes.text(red_mirr_dy_x_flex,red_mirr_dy_y_flex,'mirr_dy',fontsize=12,color='blue')
	    #self.axes.arrow(0,0,red_mirr_dz_x_flex,red_mirr_dz_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color='green',alpha=.5)
	    #self.axes.text(red_mirr_dz_x_flex,red_mirr_dz_y_flex,'mirr_dz',fontsize=12,color='green')
	    self.axes.arrow(0,0,red_mirr_tx_x_flex,red_mirr_tx_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    self.axes.text(red_mirr_tx_x_flex+.025*maxx,red_mirr_tx_y_flex,'mirr_tx',fontsize=12,color='red')
	    self.axes.arrow(0,0,red_mirr_ty_x_flex,red_mirr_ty_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'red',alpha=.5)
	    self.axes.text(red_mirr_ty_x_flex,red_mirr_ty_y_flex+.025*maxx,'mirr_ty',fontsize=12,color='red')
	    #self.axes.arrow(0,0,red_mirr_tz_x_flex,red_mirr_tz_y_flex,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'purple',alpha=.5)
	    #self.axes.text(red_mirr_tz_x_flex,red_mirr_tz_y_flex,'mirr_tz',fontsize=12,color='purple')
	    if resboxval==True:
		self.axes.arrow(0,0,red_flex_x,red_flex_y,head_width=.05*maxx, head_length=.1*maxx,width=(.005*maxx), color = 'red',alpha=.5)
		self.axes.text(red_flex_x,red_flex_y+.025*maxx,'Flex_Total',fontsize=12,color='red')
	    
        if compboxval==True:
		
	    if 0 in active:
		self.axes.arrow(0,0,red_coll_dx_x_comp,red_coll_dx_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'green',alpha=.5)
		self.axes.text(red_coll_dx_x_comp,red_coll_dx_y_comp+.025*maxx,'coll_dx',fontsize=12,color = 'green')
		self.axes.arrow(0,0,red_coll_dy_x_comp,red_coll_dy_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'green',alpha=.5)
		self.axes.text(red_coll_dy_x_comp+.025*maxx,red_coll_dy_y_comp,'coll_dy',fontsize=12,color='green')
	    #self.axes.arrow(0,0,red_coll_dz_x_comp,red_coll_dz_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color='green',alpha=.5)
	    #self.axes.text(red_coll_dz_x_comp,red_coll_dz_y_comp,'coll_dz',fontsize=12,color='green')
		self.axes.arrow(0,0,red_coll_tx_x_comp,red_coll_tx_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'green',alpha=.5)
		self.axes.text(red_coll_tx_x_comp+.025*maxx,red_coll_tx_y_comp,'coll_tx',fontsize=12,color='green')
		self.axes.arrow(0,0,red_coll_ty_x_comp,red_coll_ty_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'green',alpha=.5)
		self.axes.text(red_coll_ty_x_comp,red_coll_ty_y_comp+.025*maxx,'coll_ty',fontsize=12,color='green')
	    #self.axes.arrow(0,0,red_coll_tz_x_comp,red_coll_tz_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'purple',alpha=.5)
	    #self.axes.text(red_coll_tz_x_comp,red_coll_tz_y_comp,'coll_tz',fontsize=12,color='purple')
	    
	    if 6 in active:
		self.axes.arrow(0,0,red_det_x_comp,red_det_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'green',alpha=.5)
		self.axes.text(red_det_x_comp+.025*maxx,red_det_y_comp+.025*maxx,'det_shft',fontsize=12,color='green')
    
    
	    
	    #self.axes.arrow(0,0,red_mirr_dx_x_comp,red_mirr_dx_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'green',alpha=.5)
	    #self.axes.text(red_mirr_dx_x_comp,red_mirr_dx_y_comp,'mirr_dx',fontsize=12,color = 'green')
	    #self.axes.arrow(0,0,red_mirr_dy_x_comp,red_mirr_dy_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'blue',alpha=.5)
	    #self.axes.text(red_mirr_dy_x_comp,red_mirr_dy_y_comp,'mirr_dy',fontsize=12,color='blue')
	    #self.axes.arrow(0,0,red_mirr_dz_x_comp,red_mirr_dz_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color='green',alpha=.5)
	    #self.axes.text(red_mirr_dz_x_comp,red_mirr_dz_y_comp,'mirr_dz',fontsize=12,color='green')
	    if 1 in active:
		self.axes.arrow(0,0,red_mirr_tx_x_comp,red_mirr_tx_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'green',alpha=.5)
		self.axes.text(red_mirr_tx_x_comp+.025*maxx,red_mirr_tx_y_comp,'mirr_tx',fontsize=12,color='green')
		self.axes.arrow(0,0,red_mirr_ty_x_comp,red_mirr_ty_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'green',alpha=.5)
		self.axes.text(red_mirr_ty_x_comp,red_mirr_ty_y_comp+.025*maxx,'mirr_ty',fontsize=12,color='green')
	    #self.axes.arrow(0,0,red_mirr_tz_x_comp,red_mirr_tz_y_comp,head_width=.025*maxx, head_length=.05*maxx,width=(.005*maxx), color = 'purple',alpha=.5)
	    #self.axes.text(red_mirr_tz_x_comp,red_mirr_tz_y_comp,'mirr_tz',fontsize=12,color='purple')
	    if resboxval==True:
		self.axes.arrow(0,0,red_comp_x,red_comp_y,head_width=.05*maxx, head_length=.1*maxx,width=(.005*maxx), color = 'green',alpha=.5)
		self.axes.text(red_comp_x,red_comp_y+.025*maxx,'Comp_Total',fontsize=12,color='green')    
	    
	    
	    
	    
	    
	
	
	
        self.axes.set_ylabel("Shift on Detector (mm)", fontsize=12)
        self.axes.set_xlabel("Shift on Detector (mm)", fontsize=12)
        self.axes.axhline(0, color='black',linestyle='--')
        self.axes.axvline(0, color='black',linestyle='--')
	self.axes.axis([-maxx,maxx,-maxx,maxx])
	self.canvas.draw()
	
    def drawred(self):
	global redo_x,red_y,blueo_x,blueo_y
	self.axes.clear()
        self.axes.scatter(redo_x, redo_y,color='red',s=7)
	self.axes.set_ylabel("Dispersion-Direction mm", fontsize=12)
        self.axes.set_xlabel("Slit-Direction mm", fontsize=12)
	self.axes.set_title("Red FCS Field", fontsize=14)
	self.axes.axis([-150,150,-150,150])
	self.canvas.draw()
    
    def drawredfull(self):
	global red_x,red_y
	self.axes.clear()
        self.axes.scatter(red_x, red_y,color='red',s=7)
        self.axes.set_ylabel("Dispersion-Direction mm", fontsize=12)
        self.axes.set_xlabel("Slit-Direction mm", fontsize=12)
	self.axes.set_title("Red Full Field", fontsize=14)
	#at = AnchoredText("Pixel",prop=dict(size=15), frameon=False,loc=2)
	#box.drawing_area.add_artist(at)
	#self.axes.annotate('Pixel Relative Size', xy=(.2, .1),xycoords='axes fraction', xytext=(.3, .94))
	#self.axes.add_artist(at)	
        self.axes.callbacks.connect('xlim_changed', on_xlims_change)
        self.axes.callbacks.connect('ylim_changed', on_ylims_change)
	
	self.axes.axis([-150,150,-150,150])
	self.canvas.draw()
    def drawbluefull(self):
	global blue_x,blue_y
	self.axes.clear()
	self.axes.scatter(blue_x, blue_y,color='blue',s=7)
        self.axes.set_ylabel("Dispersion-Direction mm", fontsize=12)
        self.axes.set_xlabel("Slit-Direction mm", fontsize=12)
	self.axes.set_title("Blue Full Field", fontsize=14)
	self.axes.axis([-150,150,-150,150])	
	self.canvas.draw()
    def drawblueshiftreal(self):
	global red_x1,red_y1,blue_x1,blue_y1
	self.axes.clear()
	self.axes.scatter(blue_x, blue_y,s=7)
	self.axes.hold('True')
	self.axes.scatter(blue_x+(blue_x1-blue_x), blue_y+(blue_y1-blue_y),color='red',s=7)
	for i in range(blue_x.size):
	   self.axes.plot([blue_x[i],blue_x[i]+(blue_x1[i]-blue_x[i])],[blue_y[i],blue_y[i]+(blue_y1[i]-blue_y[i])],color='blue',linewidth=2, alpha=.6) 
	#box = AnchoredAuxTransformBox(ax.transData, loc=2)
        ##el = Ellipse((0,0), width=1, height=1, angle=30) # in data coordinates!
	##box.drawing_area.add_artist(el)
        #self.axes.add_artist(box)
	self.axes.hold('False')
	self.axes.set_ylabel("Dispersion-Direction mm", fontsize=12)
        self.axes.set_xlabel("Slit-Direction mm", fontsize=12)
	self.axes.set_title("Blue flexure Effect", fontsize=14)
        box = AnchoredAuxTransformBox(self.axes.transData, loc=2,frameon=False)
        el = Rectangle((0,0), width=15, height=15,alpha=.5) 
        box.drawing_area.add_artist(el)
	self.axes.add_artist(box)
	self.axes.axis([-150,150,-150,150])
	self.canvas.draw()
	
    def drawredshiftreal(self):
	global red_x1,red_y1,blue_x1,blue_y1
	self.axes.clear()
	self.axes.scatter(red_x, red_y,s=7)
	self.axes.hold('True')
	self.axes.scatter(red_x+(red_x1-red_x), red_y+(red_y1-red_y),color='red',s=7)
	for i in range(red_x.size):
	   self.axes.plot([red_x[i],red_x[i]+(red_x1[i]-red_x[i])],[red_y[i],red_y[i]+(red_y1[i]-red_y[i])],color='red',linewidth=2, alpha=.6) 
	#box = AnchoredAuxTransformBox(ax.transData, loc=2)
        #el = Ellipse((0,0), width=.1, height=.1, angle=30) # in data coordinates!
	#box.drawing_area.add_artist(el)
        #self.axes.add_artist(box)
	self.axes.hold('False')
        self.axes.set_ylabel("Dispersion-Direction mm", fontsize=12)
        self.axes.set_xlabel("Slit-Direction mm", fontsize=12)
	self.axes.set_title("Red flexure Effect", fontsize=14)	
        box = AnchoredAuxTransformBox(self.axes.transData, loc=2,frameon=False)
        el = Rectangle((0,0), width=15, height=15,alpha=.5) 
        box.drawing_area.add_artist(el)
	self.axes.add_artist(box)
	self.axes.axis([-150,150,-150,150])
	self.canvas.draw()

	
    def drawblueshift(self):
	global red_x1,red_y1,blue_x1,blue_y1
	self.axes.clear()
	self.axes.scatter(blue_x, blue_y,color='black',s=7)
	self.axes.hold('True')
	self.axes.scatter(blue_x+(blue_x1-blue_x)*self.flexscale, blue_y+(blue_y1-blue_y)*self.flexscale,color='blue',s=7)
	for i in range(blue_x.size):
	   self.axes.plot([blue_x[i],blue_x[i]+(blue_x1[i]-blue_x[i])*self.flexscale],[blue_y[i],blue_y[i]+(blue_y1[i]-blue_y[i])*self.flexscale],color='green',linewidth=2, alpha=.6) 
	#box = AnchoredAuxTransformBox(ax.transData, loc=2)
        ##el = Ellipse((0,0), width=1, height=1, angle=30) # in data coordinates!
	##box.drawing_area.add_artist(el)
        #self.axes.add_artist(box)
	self.axes.hold('False')
	self.axes.set_ylabel("Dispersion-Direction mm", fontsize=12)
        self.axes.set_xlabel("Slit-Direction mm", fontsize=12)
	self.axes.set_title("Blue flexure Effect", fontsize=14)
        box = AnchoredAuxTransformBox(self.axes.transData, loc=2,frameon=False)
        el = Rectangle((0,0), width=15, height=15,alpha=.5) 
        box.drawing_area.add_artist(el)
	self.axes.add_artist(box)
	self.axes.axis([-150,150,-150,150])
	self.canvas.draw()
	
    def drawredshift(self):
	global red_x1,red_y1,blue_x1,blue_y1
	self.axes.clear()
	self.axes.scatter(red_x, red_y,color='black',s=7)
	self.axes.hold('True')
	self.axes.scatter(red_x+(red_x1-red_x)*self.flexscale, red_y+(red_y1-red_y)*self.flexscale,color='red',s=7)
	for i in range(red_x.size):
	   self.axes.plot([red_x[i],red_x[i]+(red_x1[i]-red_x[i])*self.flexscale],[red_y[i],red_y[i]+(red_y1[i]-red_y[i])*self.flexscale],color='green',linewidth=2, alpha=.6) 
	#box = AnchoredAuxTransformBox(ax.transData, loc=2)
        #el = Ellipse((0,0), width=.1, height=.1, angle=30) # in data coordinates!
	#box.drawing_area.add_artist(el)
        #self.axes.add_artist(box)
	self.axes.hold('False')
        self.axes.set_ylabel("Dispersion-Direction mm", fontsize=12)
        self.axes.set_xlabel("Slit-Direction mm", fontsize=12)
	self.axes.set_title("Red flexure Effect", fontsize=14)	
        box = AnchoredAuxTransformBox(self.axes.transData, loc=2,frameon=False)
        el = Rectangle((0,0), width=(.015*self.flexscale), height=(.015*self.flexscale),alpha=.5) 
        box.drawing_area.add_artist(el)
	self.axes.add_artist(box)
	self.axes.axis([-150,150,-150,150])
	self.canvas.draw()

    def drawredres(self):
	global red_x1,red_y1,red_x2,red_y2
	self.axes.clear()
	self.axes.scatter(red_x, red_y,color='black',s=7)
	self.axes.hold('True')
	self.axes.scatter(red_x+(red_x2-red_x)*self.resscale, red_y+(red_y2-red_y)*self.resscale,color='red',s=7)
	for i in range(red_x.size):
	   self.axes.plot([red_x[i],red_x[i]+(red_x2[i]-red_x[i])*self.resscale],[red_y[i],red_y[i]+(red_y2[i]-red_y[i])*self.resscale],color='green',linewidth=2, alpha=.6) 

	self.axes.hold('False')
        self.axes.set_ylabel("Dispersion-Direction mm", fontsize=12)
        self.axes.set_xlabel("Slit-Direction mm", fontsize=12)
	self.axes.set_title("Red Residuals", fontsize=14)	
        box = AnchoredAuxTransformBox(self.axes.transData, loc=2,frameon=False)
        el = Rectangle((0,0), width=(.015*self.resscale), height=(.015*self.resscale),alpha=.3) 
        box.drawing_area.add_artist(el)
	self.axes.add_artist(box)
	self.axes.axis([-150,150,-150,150])
	self.canvas.draw()

    def drawblueres(self):
	global blue_x1,blue_y1,blue_x2,blue_y2
	self.axes.clear()
	self.axes.scatter(blue_x, blue_y,color='black',s=7)
	self.axes.hold('True')
       # self.axes.scatter(blue_x2, blue_y2,color='blue')
	self.axes.scatter(blue_x+(blue_x2-blue_x)*self.resscale, blue_y+(blue_y2-blue_y)*self.resscale,color='blue',s=7)	
	for i in range(blue_x.size):
	   self.axes.plot([blue_x[i],blue_x[i]+(blue_x2[i]-blue_x[i])*self.resscale],[blue_y[i],blue_y[i]+(blue_y2[i]-blue_y[i])*self.resscale],color='green',linewidth=2, alpha=.6) 
	self.axes.hold('False')
        self.axes.set_ylabel("Dispersion-Direction mm", fontsize=12)
        self.axes.set_xlabel("Slit-Direction mm", fontsize=12)
	self.axes.set_title("Blue Residuals", fontsize=14)	
        box = AnchoredAuxTransformBox(self.axes.transData, loc=2,frameon=False)
        el = Rectangle((0,0), width=(.015*self.resscale), height=(.015*self.resscale),alpha=.3) 
        box.drawing_area.add_artist(el)
	self.axes.add_artist(box)
	self.axes.axis([-150,150,-150,150])	
	self.canvas.draw()
	
	
    def drawredcomp(self):
	global red_x,red_y,red_x1,red_y1,red_x2,red_y2
	self.axes.clear()
	self.axes.scatter(red_x, red_y,color='black',s=7)
	self.axes.hold('True')
	self.axes.scatter(red_x+(red_x1-red_x)*1000, red_y+(red_y1-red_y)*1000,color='grey',s=7)
	self.axes.scatter(red_x+(red_x2-red_x)*1000, red_y+(red_y2-red_y)*1000,color='red',s=7)
	for i in range(red_x.size):
	   self.axes.plot([red_x[i],red_x[i]+(red_x1[i]-red_x[i])*1000],[red_y[i],red_y[i]+(red_y1[i]-red_y[i])*1000],color='green',linewidth=2, alpha=.6)
	   #self.axes.quiver([red_x[i],red_y[i]],red_x2[i]-red_x[i],red_y2[i]-red_y[i])
	for i in range(red_x.size):
	   self.axes.plot([red_x[i]+(red_x1[i]-red_x[i])*1000,red_x[i]+(red_x2[i]-red_x[i])*1000],[red_y[i]+(red_y1[i]-red_y[i])*1000,red_y[i]+(red_y2[i]-red_y[i])*1000],color='green',linewidth=2, alpha=.6)	   
        self.axes.hold('False')
	self.axes.set_ylabel("Dispersion-Direction mm", fontsize=12)
        self.axes.set_xlabel("Slit-Direction mm", fontsize=12)
	self.axes.set_title("Red Compensation Path", fontsize=14)	
        box = AnchoredAuxTransformBox(self.axes.transData, loc=2,frameon=False)
        el = Rectangle((0,0), width=15, height=15,alpha=.5) 
        box.drawing_area.add_artist(el)
	self.axes.add_artist(box)
	self.axes.axis([-150,150,-150,150])
	self.canvas.draw()
   
    def drawbluecomp(self):
	global blue_x,blue_y,blue_x1,blue_y1,blue_x2,blue_y2
	self.axes.clear()
	self.axes.scatter(blue_x, blue_y,color='black',s=7)
	self.axes.hold('True')
	self.axes.scatter(blue_x+(blue_x1-blue_x)*1000, blue_y+(blue_y1-blue_y)*1000,color='grey',s=7)
	self.axes.scatter(blue_x+(blue_x2-blue_x)*1000, blue_y+(blue_y2-blue_y)*1000,color='blue',s=7)
	for i in range(blue_x.size):
	   self.axes.plot([blue_x[i],blue_x[i]+(blue_x1[i]-blue_x[i])*1000],[blue_y[i],blue_y[i]+(blue_y1[i]-blue_y[i])*1000],color='green',linewidth=2, alpha=.6)
	   #self.axes.quiver([blue_x[i],blue_y[i]],blue_x2[i]-blue_x[i],blue_y2[i]-blue_y[i])
	for i in range(blue_x.size):
	   self.axes.plot([blue_x[i]+(blue_x1[i]-blue_x[i])*1000,blue_x[i]+(blue_x2[i]-blue_x[i])*1000],[blue_y[i]+(blue_y1[i]-blue_y[i])*1000,blue_y[i]+(blue_y2[i]-blue_y[i])*1000],color='green',linewidth=2, alpha=.6)	   

	self.axes.hold('False')
        self.axes.set_ylabel("Dispersion-Direction mm", fontsize=12)
        self.axes.set_xlabel("Slit-Direction mm", fontsize=12)
	self.axes.set_title("Blue Compensation Path", fontsize=14)	
        box = AnchoredAuxTransformBox(self.axes.transData, loc=2,frameon=False)
        el = Rectangle((0,0), width=15, height=15,alpha=.5) 
        box.drawing_area.add_artist(el)
	self.axes.add_artist(box)
	self.axes.axis([-150,150,-150,150])
	self.canvas.draw()
	
	
	
class TabOne(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
	self.selectBtn = wx.Button(self, label="Compensate")
        self.selectBtn.Bind(wx.EVT_BUTTON, self.onGetSelection)
	
	self.selectBtn1 = wx.Button(self, label="Show Path")
        self.selectBtn1.Bind(wx.EVT_BUTTON, self.onShowComp)
	
	self.selectBtn2 = wx.Button(self, label="Show Residuals")
        self.selectBtn2.Bind(wx.EVT_BUTTON, self.onShowComp)
	
	self.cb = wx.CheckBox(self, -1,'Use edge field for optimization', (10, 10))
	self.cb.SetValue(False)
	self.quote4 = wx.StaticText(self, label=" Number of Iterations", pos=(20, 30))
	self.selectIter = wx.TextCtrl(self,value='10')
	
	self.quote3 = wx.StaticText(self, label=" Residual.View.Scale:", pos=(20, 30))
	self.selectResScale = wx.TextCtrl(self,value='10000')
	
	lblname = wx.StaticText(self, label="Select Active Elements")
	self.list_choice = ['Collimator','Fold_Mirror', 'Dichroic',  'Red_Channel_Grating', 'Blue_Channel_Grating', 'Detector_Rotation','Detector_Translation']
        pos = (5, 20)
        self.list_param = wx.CheckListBox(self, wx.ID_ANY, pos, wx.DefaultSize,
                                         self.list_choice, style=1)
	self.list_col = ['dx','dy', 'dz','tx','ty','tz']
        pos = (5, 20)
        self.list_coll = wx.CheckListBox(self, wx.ID_ANY, pos, wx.DefaultSize,
                                         self.list_col, style=1)

	
	self.list_mir = ['dx','dy', 'dz','tx','ty','tz']
        pos = (5, 20)
        self.list_mirr = wx.CheckListBox(self, wx.ID_ANY, pos, wx.DefaultSize,
                                         self.list_mir, style=1) 	

	self.list_dch = ['dx','dy', 'dz','tx','ty','tz']
        pos = (5, 20)
        self.list_dchr = wx.CheckListBox(self, wx.ID_ANY, pos, wx.DefaultSize,
                                         self.list_dch, style=1) 					 			 

	self.list_rgr = ['dx','dy', 'dz','tx','ty','tz']
        pos = (5, 20)
        self.list_rgra = wx.CheckListBox(self, wx.ID_ANY, pos, wx.DefaultSize,
                                         self.list_rgr, style=1) 	

	self.list_bgr = ['dx','dy', 'dz','tx','ty','tz']
        pos = (5, 20)
        self.list_bgra = wx.CheckListBox(self, wx.ID_ANY, pos, wx.DefaultSize,
                                         self.list_bgr, style=1) 	
        self.list_coll.SetChecked([0,1,2,3,4,5])
        self.list_mirr.SetChecked([0,1,2,3,4,5])
        self.list_dchr.SetChecked([0,1,2,3,4,5])
        self.list_rgra.SetChecked([0,1,2,3,4,5])
        self.list_bgra.SetChecked([0,1,2,3,4,5])
										 				 
	self.list_mirr.Hide()
#	self.list_coll.Hide()
	self.list_bgra.Hide()
	self.list_rgra.Hide()
	self.list_dchr.Hide()
					 				                                   
	sizer1 = wx.BoxSizer(wx.VERTICAL)
	#sizerhor = wx.BoxSizer(wx.HORIZONTAL)
	sizer1.Add(lblname, 0, wx.ALL|wx.CENTER, 5)
	sizer1.Add(self.list_param, 0, wx.ALL|wx.CENTER, 5)
	sizer1.Add(self.selectBtn, 0, wx.ALL|wx.CENTER, 5)
	
	sizer2 = wx.BoxSizer(wx.VERTICAL)
	
	sizerhor = wx.BoxSizer(wx.HORIZONTAL)
	sizerhor.Add(self.quote3, 0, wx.ALL|wx.CENTER, 5)
	sizerhor.Add(self.selectResScale, 0, wx.ALL|wx.CENTER, 5)
	
	sizerhor1 = wx.BoxSizer(wx.HORIZONTAL)
	sizerhor1.Add(self.quote4, 0, wx.ALL|wx.CENTER, 5)
	sizerhor1.Add(self.selectIter, 0, wx.ALL|wx.CENTER, 5)
	
	sizer2.Add(self.cb, 0, wx.ALL|wx.CENTER, 5)
	sizer2.Add(sizerhor1, 0, wx.ALL|wx.CENTER, 5)
	sizer2.Add(sizerhor, 0, wx.ALL|wx.CENTER, 5)
	sizer2.Add(self.selectBtn1, 0, wx.ALL|wx.CENTER, 5)
	sizer2.Add(self.selectBtn2, 0, wx.ALL|wx.CENTER, 5)
	#sizer.Add(sizerhor, 0, wx.ALL|wx.CENTER, 5)
	
	sizer = wx.BoxSizer(wx.HORIZONTAL)
	sizer.Add(sizer1, 0, wx.ALL|wx.CENTER, 5)
	sizer.Add(self.list_coll, 0, wx.ALL|wx.CENTER, 5)
	sizer.Add(self.list_mirr, 0, wx.ALL|wx.CENTER, 5)
	sizer.Add(self.list_dchr, 0, wx.ALL|wx.CENTER, 5)
	sizer.Add(self.list_rgra, 0, wx.ALL|wx.CENTER, 5)
	sizer.Add(self.list_bgra, 0, wx.ALL|wx.CENTER, 5)
	sizer.Add(sizer2, 1, wx.ALL|wx.CENTER, 5)
        self.SetSizer(sizer)
    def onGetSelection(self, event):
	print 'hello'
    def onShowComp(self, event):
	print 'hello'

#Grid Input for flexure
class Tabtwo(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
	self.currentlySelectedCell = (0, 0)
        self.myGrid = gridlib.Grid(self)
        self.myGrid.CreateGrid(9, 6)
        self.myGrid.SetColLabelValue(0, "dx")
        self.myGrid.SetColLabelValue(1, "dy")
        self.myGrid.SetColLabelValue(2, "dz")
        self.myGrid.SetColLabelValue(3, "tx")
        self.myGrid.SetColLabelValue(4, "ty")
        self.myGrid.SetColLabelValue(5, "tz")
        self.myGrid.SetRowLabelValue(0, "Collimator")
	self.myGrid.SetRowLabelValue(1, "Mirror")
        self.myGrid.SetRowLabelValue(2, "Dichroic")
        self.myGrid.SetRowLabelValue(3, "Red_Grating")
        self.myGrid.SetRowLabelValue(4, "Blue Grating")
	self.myGrid.SetRowLabelValue(5, " ")	
	self.myGrid.SetRowLabelValue(6, "Det Motion")
	self.myGrid.SetRowLabelValue(7, " ")
	self.myGrid.SetRowLabelValue(8, "Det Rotation (deg)")
        self.myGrid.SetRowLabelSize(150)
	#self.myGrid.Bind(gridlib.EVT_GRID_SELECT_CELL, self.onSingleSelect)
        #self.myGrid.Bind(gridlib.EVT_GRID_RANGE_SELECT, self.onDragSelection)
	
        self.quote1 = wx.StaticText(self, label=" dx,dy,dz  [ 1 unit = .001 mm ]   tx,ty,tz [1 unit = .00001 degree]", style=wx.ALIGN_RIGHT)
	self.quote2 = wx.StaticText(self, label=" Range: Dis +/-", pos=(20, 30))
	self.quote5 = wx.StaticText(self, label=" Rot +/-", pos=(20, 30))
	self.quote3 = wx.StaticText(self, label=" Flex.View.Scale:", pos=(20, 30))
	self.myGrid.SetCellValue(5,0,'Red x')
	self.myGrid.SetCellValue(5,1,'Red y')	
	self.myGrid.SetCellValue(5,2,'Blue x')
	self.myGrid.SetCellValue(5,3,'Blue y')		
	self.myGrid.SetCellValue(7,0,'Red')
	self.myGrid.SetCellValue(7,1,'Blue')

	self.myGrid.SetCellBackgroundColour(5,0,(240,240,240))
	self.myGrid.SetCellBackgroundColour(5,1,(240,240,240))
	self.myGrid.SetCellBackgroundColour(5,2,(240,240,240))
	self.myGrid.SetCellBackgroundColour(5,3,(240,240,240))
	
	self.myGrid.SetCellBackgroundColour(7,0,(240,240,240))
	self.myGrid.SetCellBackgroundColour(7,1,(240,240,240))
		
	for i in range(5):
	    for j in range(6):
	       self.myGrid.SetCellValue(i,j,'0')

	for i in range(6,9,2):
	    for j in range(2):
	       self.myGrid.SetCellValue(i,j,'0')
	self.myGrid.SetCellValue(6,2,'0')   
	self.myGrid.SetCellValue(6,3,'0')     
	
	
	
	DropDownList = []
        Options = {"Reset","Random","Collimator","All","Grating"}
        for TextNode in Options:
             DropDownList.append( TextNode )
	self.cmbx=wx.ComboBox(self,value="Select Preset flexure",choices=DropDownList)

	DropDownList = []
        Options = {"15out6prx","14out9prx","15out6cam","MOBIE"}
        for TextNode in Options:
             DropDownList.append( TextNode )
	self.configcmbx=wx.ComboBox(self,value="Select Config:",choices=DropDownList)
	
	self.selectBtn1 = wx.Button(self, label="Show FCS field")
	self.selectBtn4 = wx.Button(self, label="Save")
	self.quote4 = wx.StaticText(self, label=" File Name:", pos=(20, 30))	
	self.selectFile = wx.TextCtrl(self,value='flexure.txt',size=(200,-1))
	self.selectBtn5 = wx.Button(self, label="Load")

	#self.selectBtn1.Bind(wx.EVT_BUTTON, self.onFCSButton)
	
	self.selectBtn3 = wx.Button(self, label="Show Full field")
	#self.selectBtn1.Bind(wx.EVT_BUTTON, self.onFCSButton)
	
	self.selectBtn2 = wx.Button(self, label="Flex")
	self.selectDisRange = wx.TextCtrl(self,value='1000')
	self.selectRotRange = wx.TextCtrl(self,value='1000')
	self.selectFlexScale = wx.TextCtrl(self,value='1000')
	#self.selectBtn2.Bind(wx.EVT_BUTTON, self.onFlexButton)
          
	sizergrid = wx.BoxSizer(wx.VERTICAL)
	sizergrid.Add(self.myGrid, 1, wx.EXPAND)
	sizergrid.Add(self.quote1, 0, wx.EXPAND)
	#sizergrid.Add(self.quote2, 0, wx.EXPAND)
	#sizergrid.Add(self.quote2, 0, wx.EXPAND)
        sizergridhor = wx.BoxSizer(wx.HORIZONTAL)
        #sizerhor.Add(selectBtn, 0, wx.ALL|wx.CENTER, 5)
	sizergridhor.Add(self.configcmbx, 0, wx.ALL|wx.CENTER, 5)
	sizergridhor.Add(self.selectBtn1, 0, wx.ALL|wx.CENTER, 5)
	sizergridhor.Add(self.selectBtn3, 0, wx.ALL|wx.CENTER, 5)
	sizergridhor.Add(self.quote3, 0, wx.ALL|wx.CENTER, 5)
	sizergridhor.Add(self.selectFlexScale, 0, wx.ALL|wx.CENTER, 5)
	
		
	sizergridhor2 = wx.BoxSizer(wx.HORIZONTAL)
	sizergridhor2.Add(self.selectBtn2, 0, wx.ALL|wx.CENTER, 5)
	sizergridhor2.Add(self.cmbx, 0, wx.ALL|wx.CENTER, 5)
	sizergridhor2.Add(self.quote2, 0, wx.ALL|wx.CENTER, 5)	
	sizergridhor2.Add(self.selectDisRange, 0, wx.ALL|wx.CENTER, 5)	
	sizergridhor2.Add(self.quote5, 0, wx.ALL|wx.CENTER, 5)	
	sizergridhor2.Add(self.selectRotRange, 0, wx.ALL|wx.CENTER, 5)		

	sizergridhor3 = wx.BoxSizer(wx.HORIZONTAL)
	sizergridhor3.Add(self.selectBtn5, 0, wx.ALL|wx.CENTER, 5)	
	sizergridhor3.Add(self.selectBtn4, 0, wx.ALL|wx.CENTER, 5)
	sizergridhor3.Add(self.quote4, 0, wx.ALL|wx.CENTER, 5)	
	sizergridhor3.Add(self.selectFile, 0, wx.ALL|wx.CENTER, 5)		

	
	sizergrid.Add(sizergridhor3, 0, wx.ALL|wx.CENTER, 5)
	sizergrid.Add(sizergridhor, 0, wx.ALL|wx.CENTER, 5)
	sizergrid.Add(sizergridhor2, 0, wx.ALL|wx.CENTER, 5)
	self.SetSizer(sizergrid)    


	
class Tabthree(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
	self.currentlySelectedCell = (0, 0)
        self.myGrid = gridlib.Grid(self)
        self.myGrid.CreateGrid(9, 6)
        self.myGrid.SetColLabelValue(0, "dx")
        self.myGrid.SetColLabelValue(1, "dy")
        self.myGrid.SetColLabelValue(2, "dz")
        self.myGrid.SetColLabelValue(3, "tx")
        self.myGrid.SetColLabelValue(4, "ty")
        self.myGrid.SetColLabelValue(5, "tz")
        self.myGrid.SetRowLabelValue(0, "Collimator")
	self.myGrid.SetRowLabelValue(1, "Mirror")
        self.myGrid.SetRowLabelValue(2, "Dichroic")
        self.myGrid.SetRowLabelValue(3, "Red_Grating")
        self.myGrid.SetRowLabelValue(4, "Blue Grating")
	self.myGrid.SetRowLabelValue(5, " ")	
	self.myGrid.SetRowLabelValue(6, "Det Motion")
	self.myGrid.SetRowLabelValue(7, " ")
	self.myGrid.SetRowLabelValue(8, "Det Rotation (deg)")
        self.myGrid.SetRowLabelSize(150)
	#self.myGrid.Bind(gridlib.EVT_GRID_SELECT_CELL, self.onSingleSelect)
        #self.myGrid.Bind(gridlib.EVT_GRID_RANGE_SELECT, self.onDragSelection)
	
        self.quote1 = wx.StaticText(self, label=" dx,dy,dz  [ 1 unit = .001 mm ]   tx,ty,tz [1 unit = .00001 degree]", style=wx.ALIGN_RIGHT)
	self.quote2 = wx.StaticText(self, label=" Range: +/-", pos=(20, 30))
	self.quote3 = wx.StaticText(self, label=" Flex.View.Scale:", pos=(20, 30))
	self.myGrid.SetCellValue(5,0,'Red x')
	self.myGrid.SetCellValue(5,1,'Red y')	
	self.myGrid.SetCellValue(5,2,'Blue x')
	self.myGrid.SetCellValue(5,3,'Blue y')		
	self.myGrid.SetCellValue(7,0,'Red')
	self.myGrid.SetCellValue(7,1,'Blue')

	self.myGrid.SetCellBackgroundColour(5,0,(240,240,240))
	self.myGrid.SetCellBackgroundColour(5,1,(240,240,240))
	self.myGrid.SetCellBackgroundColour(5,2,(240,240,240))
	self.myGrid.SetCellBackgroundColour(5,3,(240,240,240))
	
	self.myGrid.SetCellBackgroundColour(7,0,(240,240,240))
	self.myGrid.SetCellBackgroundColour(7,1,(240,240,240))
		
	for i in range(5):
	    for j in range(6):
	       self.myGrid.SetCellValue(i,j,'0')

	for i in range(6,9,2):
	    for j in range(2):
	       self.myGrid.SetCellValue(i,j,'0')
	self.myGrid.SetCellValue(6,2,'0')   
	self.myGrid.SetCellValue(6,3,'0')     
	
	self.selectBtn4 = wx.Button(self, label="Save")
	self.quote4 = wx.StaticText(self, label=" File Name:", pos=(20, 30))	
	self.selectFile = wx.TextCtrl(self,value='compensation.txt',size=(200,-1))	
	self.quote5 = wx.StaticText(self, label=" Flexure/Compensation Components", pos=(20, 30))	
	self.selectBtn5 = wx.Button(self, label="Analyze")
	self.flexbox=wx.CheckBox(self, label = 'flexure',pos = (10,10))
	self.compbox=wx.CheckBox(self, label = 'compensation',pos = (10,10))
	self.resbox=wx.CheckBox(self, label = 'Resultant',pos = (10,10))	
	sizergrid = wx.BoxSizer(wx.VERTICAL)
	sizergridhor3 = wx.BoxSizer(wx.HORIZONTAL)
	sizergridhor3.Add(self.selectBtn4, 0, wx.ALL|wx.CENTER, 5)
	sizergridhor3.Add(self.quote4, 0, wx.ALL|wx.CENTER, 5)	
	sizergridhor3.Add(self.selectFile, 0, wx.ALL|wx.CENTER, 5)
	sizergridhor4 = wx.BoxSizer(wx.HORIZONTAL)	
	sizergridhor4.Add(self.quote5, 0, wx.ALL|wx.CENTER, 5)		
	sizergridhor4.Add(self.selectBtn5, 0, wx.ALL|wx.CENTER, 5)	
	sizergridhor4.Add(self.flexbox, 0, wx.ALL|wx.CENTER, 5)
	sizergridhor4.Add(self.compbox, 0, wx.ALL|wx.CENTER, 5)	
	sizergridhor4.Add(self.resbox, 0, wx.ALL|wx.CENTER, 5)
	sizergrid.Add(self.myGrid, 1, wx.EXPAND)
	sizergrid.Add(sizergridhor3, 0, wx.ALL|wx.CENTER, 5)
	sizergrid.Add(sizergridhor4, 0, wx.ALL|wx.CENTER, 5)	
	self.SetSizer(sizergrid)    
	
 

 
#def f(p, *params):
	#global red_x,red_y,red_coll_dx_x,red_coll_dy_x,red_coll_dz_x,red_coll_dx_y,red_coll_dy_y,red_coll_dz_y,grat_tx_x,grat_tx_y,red_x1,red_y1,red_x2,red_y2
	#global blue_x,blue_y,blue_coll_dx_x,blue_coll_dy_x,blue_coll_dz_x,blue_coll_dx_y,blue_coll_dy_y,blue_coll_dz_y,grat_tx_x,grat_tx_y,blue_x1,blue_y1,blue_x2,blue_y2
	#global bluegra_x,bluegra_y,bluegra_coll_dx_x,bluegra_coll_dy_x,bluegra_coll_dz_x,bluegra_coll_dx_y,bluegra_coll_dy_y,bluegra_coll_dz_y,grat_tx_x,grat_tx_y,bluegra_x1,bluegra_y1,bluegra_x2,bluegra_y2

	#act_coll_dx_x= p[0]
	#act_coll_dy_x= p[1]
	#act_coll_dz_x= p[2]
	#act_coll_dx_y= p[3]
	#act_coll_dy_y= p[4]
	#act_coll_dz_y= p[5]
	
	
	#det_tz=p[2]*0.0174532925 
	#red_x2=red_x1+ act_coll_dx_x*red_coll_dx_x +act_coll_dy_y*red_coll_dy_y	+ act_coll_dz_y*red_coll_dz_y	
	#red_y2=red_y1+ act_coll_dx_y*red_coll_dx_y +act_coll_dy_y*red_coll_dy_y+act_coll_dz_y*red_coll_dy_y
	#blue_x2=blue_x1+ act_coll_dx_x*blue_coll_dx_x +act_coll_dy_y*blue_coll_dy_y	+ act_coll_dz_y*blue_coll_dz_y	
	#blue_y2=blue_y1+ act_coll_dx_y*blue_coll_dx_y +act_coll_dy_y*blue_coll_dy_y+act_coll_dz_y*blue_coll_dy_y     #print x2.mean(),y2.mean()
	#err=sum(np.sqrt((red_x-red_x2)**2+(red_y-red_y2)**2))+ sum(np.sqrt((blue_x-blue_x2)**2+(blue_y-blue_y2)**2))
	#print err
	#a,b = params
	#return err

class MyPanel(wx.Panel):
    """"""
 
    #----------------------------------------------------------------------
    def __init__(self, parent):
        """Constructor"""
        wx.Panel.__init__(self, parent)
    
	log = wx.TextCtrl(self, wx.ID_ANY,style = wx.TE_MULTILINE|wx.TE_READONLY|wx.HSCROLL)
	log.SetBackgroundColour((0,100,0))
	log.SetForegroundColour((240,230,140))
	redir=RedirectText(log)
	sys.stdout=redir

	nb1 = wx.Notebook(self)
        nb2 = wx.Notebook(self)
        self.tab1 = TabOne(nb1)
	self.tab2 = Tabtwo(nb1)
	self.tab3 = Tabthree(nb1)
        self.panell = CanvasPanel(nb2)
        self.panell2 = CanvasPanel(nb2)
        sizerv = wx.BoxSizer(wx.VERTICAL)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        nb1.AddPage(self.tab2, "Flexure")
        nb1.AddPage(self.tab1, "Compensation")
        nb1.AddPage(self.tab3, "Analyze")	
        nb2.AddPage(self.panell, "Red_Camera")
        nb2.AddPage(self.panell2, "Blue_Camera")
	#self.tab2.myGrid.Bind(gridlib.EVT_GRID_SELECT_CELL, self.onSingleSelect)
        #self.tab2.myGrid.Bind(gridlib.EVT_GRID_RANGE_SELECT, self.onDragSelection)
	self.tab2.selectBtn1.Bind(wx.EVT_BUTTON, self.onFCSButton)
	self.tab2.selectBtn2.Bind(wx.EVT_BUTTON, self.onFlexButton)
	self.tab2.selectBtn3.Bind(wx.EVT_BUTTON, self.onFullButton)
	self.tab1.selectBtn1.Bind(wx.EVT_BUTTON, self.onCompButton)
	self.tab1.selectBtn2.Bind(wx.EVT_BUTTON, self.onResButton)
	self.tab1.list_param.Bind(wx.EVT_LISTBOX, self.onSelectEle)
	self.tab2.cmbx.Bind(wx.EVT_COMBOBOX, self.onCmbxSelection)
	self.tab2.configcmbx.Bind(wx.EVT_COMBOBOX, self.onConfigSelect)	
	self.tab1.selectBtn.Bind(wx.EVT_BUTTON, self.onOpButton)
	self.tab2.selectBtn4.Bind(wx.EVT_BUTTON, self.onFlexSave)
	self.tab2.selectBtn5.Bind(wx.EVT_BUTTON, self.onFlexLoad)
	self.tab3.selectBtn4.Bind(wx.EVT_BUTTON, self.onCompSave)
	self.tab3.selectBtn5.Bind(wx.EVT_BUTTON, self.onAnalyze)	
#	self.tab3.flexbox.Bind(wx.EVT_CHECKBOX,self.onflexbox)
#	self.tab3.compbox.Bind(wx.EVT_CHECKBOX,self.oncompbox)
	sizerv.Add(log, 1, wx.EXPAND)
	sizerv.Add(nb1, 2, wx.EXPAND)
        sizer.Add(sizerv, 1, wx.EXPAND)
	sizer.Add(nb2, 1, wx.EXPAND)
        self.SetSizer(sizer)

    def onSelectEle(self, event):
	ele=event.GetString()	
	if ele=='Collimator':
	    	self.tab1.list_mirr.Hide()
	        self.tab1.list_coll.Show()
	        self.tab1.list_bgra.Hide()
	        self.tab1.list_rgra.Hide()
	        self.tab1.list_dchr.Hide()
		self.tab1.Layout()
	if ele=='Fold_Mirror':
	    	self.tab1.list_mirr.Show()
	        self.tab1.list_coll.Hide()
	        self.tab1.list_bgra.Hide()
	        self.tab1.list_rgra.Hide()
	        self.tab1.list_dchr.Hide()
		self.tab1.Layout()
	if ele=='Dichroic':
	    	self.tab1.list_mirr.Hide()
	        self.tab1.list_coll.Hide()
	        self.tab1.list_bgra.Hide()
	        self.tab1.list_rgra.Hide()
	        self.tab1.list_dchr.Hide()
		self.tab1.Layout()
	if ele=='Red_Channel_Grating':
	    	self.tab1.list_mirr.Hide()
	        self.tab1.list_coll.Hide()
	        self.tab1.list_bgra.Hide()
	        self.tab1.list_rgra.Hide()
	        self.tab1.list_dchr.Hide()
		self.tab1.Layout()
	if ele=='Blue_Channel_Grating':
	    	self.tab1.list_mirr.Hide()
	        self.tab1.list_coll.Hide()
	        self.tab1.list_bgra.Hide()
	        self.tab1.list_rgra.Hide()
	        self.tab1.list_dchr.Hide()	
		self.tab1.Layout()
	if ele=='Detector_Rotation':
	    	self.tab1.list_mirr.Hide()
	        self.tab1.list_coll.Hide()
	        self.tab1.list_bgra.Hide()
	        self.tab1.list_rgra.Hide()
	        self.tab1.list_dchr.Hide()	
		self.tab1.Layout()
	if ele=='Detector_Translation':
	    	self.tab1.list_mirr.Hide()
	        self.tab1.list_coll.Hide()
	        self.tab1.list_bgra.Hide()
	        self.tab1.list_rgra.Hide()
	        self.tab1.list_dchr.Hide()	
		self.tab1.Layout()											
	
		
    def onAnalyze(self, event):	
        global red_coll_dx_x_flex,red_coll_dy_x_flex,red_coll_dz_x_flex,red_coll_tx_x_flex,red_coll_ty_x_flex,red_coll_tz_x_flex,red_gra_dx_x_flex,red_gra_dy_x_flex,red_gra_dz_x_flex,red_gra_tx_x_flex,red_gra_ty_x_flex,red_gra_tz_x_flex, red_mirr_dx_x_flex,red_mirr_dy_x_flex,red_mirr_dz_x_flex,red_mirr_tx_x_flex,red_mirr_ty_x_flex,red_mirr_tz_x_flex
        global red_coll_dx_y_flex,red_coll_dy_y_flex,red_coll_dz_y_flex,red_coll_tx_y_flex,red_coll_ty_y_flex,red_coll_tz_y_flex,red_gra_dx_y_flex,red_gra_dy_y_flex,red_gra_dz_y_flex,red_gra_tx_y_flex,red_gra_ty_y_flex,red_gra_tz_y_flex, red_mirr_dx_y_flex,red_mirr_dy_y_flex,red_mirr_dz_y_flex,red_mirr_tx_y_flex,red_mirr_ty_y_flex,red_mirr_tz_y_flex
        global red_coll_dx_x_comp,red_coll_dy_x_comp,red_coll_dz_x_comp,red_coll_tx_x_comp,red_coll_ty_x_comp,red_coll_tz_x_comp,red_det_x_comp,red_det_y_comp, red_mirr_dx_x_comp,red_mirr_dy_x_comp,red_mirr_dz_x_comp,red_mirr_tx_x_comp,red_mirr_ty_x_comp,red_mirr_tz_x_comp
        global red_coll_dx_y_comp,red_coll_dy_y_comp,red_coll_dz_y_comp,red_coll_tx_y_comp,red_coll_ty_y_comp,red_coll_tz_y_comp,red_mirr_dx_y_comp,red_mirr_dy_y_comp,red_mirr_dz_y_comp,red_mirr_tx_y_comp,red_mirr_ty_y_comp,red_mirr_tz_y_comp

        global blue_coll_dx_x_flex,blue_coll_dy_x_flex,blue_coll_dz_x_flex,blue_coll_tx_x_flex,blue_coll_ty_x_flex,blue_coll_tz_x_flex,blue_gra_dx_x_flex,blue_gra_dy_x_flex,blue_gra_dz_x_flex,blue_gra_tx_x_flex,blue_gra_ty_x_flex,blue_gra_tz_x_flex, blue_dchr_dx_x_flex,blue_dchr_dy_x_flex,blue_dchr_dz_x_flex,blue_dchr_tx_x_flex,blue_dchr_ty_x_flex,blue_dchr_tz_x_flex
        global blue_coll_dx_y_flex,blue_coll_dy_y_flex,blue_coll_dz_y_flex,blue_coll_tx_y_flex,blue_coll_ty_y_flex,blue_coll_tz_y_flex,blue_gra_dx_y_flex,blue_gra_dy_y_flex,blue_gra_dz_y_flex,blue_gra_tx_y_flex,blue_gra_ty_y_flex,blue_gra_tz_y_flex, blue_dchr_dx_y_flex,blue_dchr_dy_y_flex,blue_dchr_dz_y_flex,blue_dchr_tx_y_flex,blue_dchr_ty_y_flex,blue_dchr_tz_y_flex
        global blue_coll_dx_x_comp,blue_coll_dy_x_comp,blue_coll_dz_x_comp,blue_coll_tx_x_comp,blue_coll_ty_x_comp,blue_coll_tz_x_comp,blue_det_x_comp,blue_det_y_comp, blue_dchr_dx_x_comp,blue_dchr_dy_x_comp,blue_dchr_dz_x_comp,blue_dchr_tx_x_comp,blue_dchr_ty_x_comp,blue_dchr_tz_x_comp
        global blue_coll_dx_y_comp,blue_coll_dy_y_comp,blue_coll_dz_y_comp,blue_coll_tx_y_comp,blue_coll_ty_y_comp,blue_coll_tz_y_comp,blue_dchr_dx_y_comp,blue_dchr_dy_y_comp,blue_dchr_dz_y_comp,blue_dchr_tx_y_comp,blue_dchr_ty_y_comp,blue_dchr_tz_y_comp
        global flexboxval,compboxval,resboxval, red_flex_x,red_flex_y,red_comp_x,red_comp_y,blue_flex_x,blue_flex_y,blue_comp_x,blue_comp_y
	
	flexboxval=self.tab3.flexbox.GetValue()	
	compboxval=self.tab3.compbox.GetValue()	
	resboxval=self.tab3.resbox.GetValue()	
	blue_gra_dx_x_mn= blue_gra_dx_x.mean()
	blue_gra_dx_y_mn= blue_gra_dx_y.mean()
	blue_gra_dy_x_mn= blue_gra_dy_x.mean()
	blue_gra_dy_y_mn= blue_gra_dy_y.mean()
	blue_gra_dz_x_mn= blue_gra_dz_x.mean()
	blue_gra_dz_y_mn= blue_gra_dz_y.mean()
	blue_gra_tx_x_mn= blue_gra_tx_x.mean()
	blue_gra_tx_y_mn= blue_gra_tx_y.mean()
	blue_gra_ty_x_mn= blue_gra_ty_x.mean()
	blue_gra_ty_y_mn= blue_gra_ty_y.mean()
	blue_gra_tz_x_mn= blue_gra_tz_x.mean()
	blue_gra_tz_y_mn= blue_gra_tz_y.mean()
	
	red_gra_dx_x_mn= red_gra_dx_x.mean()
	red_gra_dx_y_mn= red_gra_dx_y.mean()
	red_gra_dy_x_mn= red_gra_dy_x.mean()
	red_gra_dy_y_mn= red_gra_dy_y.mean()
	red_gra_dz_x_mn= red_gra_dz_x.mean()
	red_gra_dz_y_mn= red_gra_dz_y.mean()
	red_gra_tx_x_mn= red_gra_tx_x.mean()
	red_gra_tx_y_mn= red_gra_tx_y.mean()
	red_gra_ty_x_mn= red_gra_ty_x.mean()
	red_gra_ty_y_mn= red_gra_ty_y.mean()
	red_gra_tz_x_mn= red_gra_tz_x.mean()
	red_gra_tz_y_mn= red_gra_tz_y.mean()
	
	blue_coll_dx_x_mn= blue_coll_dx_x.mean()
	blue_coll_dx_y_mn= blue_coll_dx_y.mean()
	blue_coll_dy_x_mn= blue_coll_dy_x.mean()
	blue_coll_dy_y_mn= blue_coll_dy_y.mean()
	blue_coll_dz_x_mn= blue_coll_dz_x.mean()
	blue_coll_dz_y_mn= blue_coll_dz_y.mean()
	blue_coll_tx_x_mn= blue_coll_tx_x.mean()
	blue_coll_tx_y_mn= blue_coll_tx_y.mean()
	blue_coll_ty_x_mn= blue_coll_ty_x.mean()
	blue_coll_ty_y_mn= blue_coll_ty_y.mean()
	blue_coll_tz_x_mn= blue_coll_tz_x.mean()
	blue_coll_tz_y_mn= blue_coll_tz_y.mean()
	
	red_coll_dx_x_mn= red_coll_dx_x.mean()
	red_coll_dx_y_mn= red_coll_dx_y.mean()
	red_coll_dy_x_mn= red_coll_dy_x.mean()
	red_coll_dy_y_mn= red_coll_dy_y.mean()
	red_coll_dz_x_mn= red_coll_dz_x.mean()
	red_coll_dz_y_mn= red_coll_dz_y.mean()
	red_coll_tx_x_mn= red_coll_tx_x.mean()
	red_coll_tx_y_mn= red_coll_tx_y.mean()
	red_coll_ty_x_mn= red_coll_ty_x.mean()
	red_coll_ty_y_mn= red_coll_ty_y.mean()
	red_coll_tz_x_mn= red_coll_tz_x.mean()
	red_coll_tz_y_mn= red_coll_tz_y.mean()
	
	blue_dchr_dx_x_mn= blue_dchr_dx_x.mean()
	blue_dchr_dx_y_mn= blue_dchr_dx_y.mean()
	blue_dchr_dy_x_mn= blue_dchr_dy_x.mean()
	blue_dchr_dy_y_mn= blue_dchr_dy_y.mean()
	blue_dchr_dz_x_mn= blue_dchr_dz_x.mean()
	blue_dchr_dz_y_mn= blue_dchr_dz_y.mean()
	blue_dchr_tx_x_mn= blue_dchr_tx_x.mean()
	blue_dchr_tx_y_mn= blue_dchr_tx_y.mean()
	blue_dchr_ty_x_mn= blue_dchr_ty_x.mean()
	blue_dchr_ty_y_mn= blue_dchr_ty_y.mean()
	blue_dchr_tz_x_mn= blue_dchr_tz_x.mean()
	blue_dchr_tz_y_mn= blue_dchr_tz_y.mean()
	
	red_mirr_dx_x_mn= red_mirr_dx_x.mean()
	red_mirr_dx_y_mn= red_mirr_dx_y.mean()
	red_mirr_dy_x_mn= red_mirr_dy_x.mean()
	red_mirr_dy_y_mn= red_mirr_dy_y.mean()
	red_mirr_dz_x_mn= red_mirr_dz_x.mean()
	red_mirr_dz_y_mn= red_mirr_dz_y.mean()
	red_mirr_tx_x_mn= red_mirr_tx_x.mean()
	red_mirr_tx_y_mn= red_mirr_tx_y.mean()
	red_mirr_ty_x_mn= red_mirr_ty_x.mean()
	red_mirr_ty_y_mn= red_mirr_ty_y.mean()
	red_mirr_tz_x_mn= red_mirr_tz_x.mean()
	red_mirr_tz_y_mn= red_mirr_tz_y.mean()

        act_coll_dx = float(self.tab2.myGrid.GetCellValue(0,0))
	act_coll_dy = float(self.tab2.myGrid.GetCellValue(0,1))
	act_coll_dz = float(self.tab2.myGrid.GetCellValue(0,2))
	act_coll_tx = float(self.tab2.myGrid.GetCellValue(0,3))
	act_coll_ty = float(self.tab2.myGrid.GetCellValue(0,4))
	act_coll_tz = float(self.tab2.myGrid.GetCellValue(0,5))
	
	act_mirr_dx = float(self.tab2.myGrid.GetCellValue(1,0))
	act_mirr_dy = float(self.tab2.myGrid.GetCellValue(1,1))
	act_mirr_dz = float(self.tab2.myGrid.GetCellValue(1,2))
	act_mirr_tx = float(self.tab2.myGrid.GetCellValue(1,3))
	act_mirr_ty = float(self.tab2.myGrid.GetCellValue(1,4))
	act_mirr_tz = float(self.tab2.myGrid.GetCellValue(1,5))
	
	act_dchr_dx = float(self.tab2.myGrid.GetCellValue(2,0))
	act_dchr_dy = float(self.tab2.myGrid.GetCellValue(2,1))
	act_dchr_dz = float(self.tab2.myGrid.GetCellValue(2,2))
	act_dchr_tx = float(self.tab2.myGrid.GetCellValue(2,3))
	act_dchr_ty = float(self.tab2.myGrid.GetCellValue(2,4))
	act_dchr_tz = float(self.tab2.myGrid.GetCellValue(2,5))
	
	act_rgra_dx = float(self.tab2.myGrid.GetCellValue(3,0))
	act_rgra_dy = float(self.tab2.myGrid.GetCellValue(3,1))
	act_rgra_dz = float(self.tab2.myGrid.GetCellValue(3,2))
	act_rgra_tx = float(self.tab2.myGrid.GetCellValue(3,3))
	act_rgra_ty = float(self.tab2.myGrid.GetCellValue(3,4))
	act_rgra_tz = float(self.tab2.myGrid.GetCellValue(3,5))
	
	act_bgra_dx = float(self.tab2.myGrid.GetCellValue(4,0))
	act_bgra_dy = float(self.tab2.myGrid.GetCellValue(4,1))
	act_bgra_dz = float(self.tab2.myGrid.GetCellValue(4,2))
	act_bgra_tx = float(self.tab2.myGrid.GetCellValue(4,3))
	act_bgra_ty = float(self.tab2.myGrid.GetCellValue(4,4))
	act_bgra_tz = float(self.tab2.myGrid.GetCellValue(4,5))

	det_red_shiftx = float(self.tab2.myGrid.GetCellValue(6,0))*1e-3
	det_red_shifty = float(self.tab2.myGrid.GetCellValue(6,1))*1e-3
	det_blue_shiftx = float(self.tab2.myGrid.GetCellValue(6,2))*1e-3
	det_blue_shifty = float(self.tab2.myGrid.GetCellValue(6,3))*1e-3			
	det_red_rot = np.radians(float(self.tab2.myGrid.GetCellValue(8,0)))
	det_blue_rot = np.radians(float(self.tab2.myGrid.GetCellValue(8,1)))
	
		
	red_coll_dx_x_flex= act_coll_dx*red_coll_dx_x_mn
	red_coll_dy_x_flex= act_coll_dy*red_coll_dy_x_mn
	red_coll_dz_x_flex= act_coll_dz*red_coll_dz_x_mn
	red_coll_tx_x_flex= act_coll_tx*red_coll_tx_x_mn
	red_coll_ty_x_flex= act_coll_ty*red_coll_ty_x_mn
	red_coll_tz_x_flex= act_coll_tz*red_coll_tz_x_mn
	red_coll_dx_y_flex= act_coll_dx*red_coll_dx_y_mn
	red_coll_dy_y_flex= act_coll_dy*red_coll_dy_y_mn
	red_coll_dz_y_flex= act_coll_dz*red_coll_dz_y_mn
	red_coll_tx_y_flex= act_coll_tx*red_coll_tx_y_mn
	red_coll_ty_y_flex= act_coll_ty*red_coll_ty_y_mn
	red_coll_tz_y_flex= act_coll_tz*red_coll_tz_y_mn
	
	red_gra_dx_x_flex= act_rgra_dx*red_gra_dx_x_mn
	red_gra_dy_x_flex= act_rgra_dy*red_gra_dy_x_mn
	red_gra_dz_x_flex= act_rgra_dz*red_gra_dz_x_mn
	red_gra_tx_x_flex= act_rgra_tx*red_gra_tx_x_mn
	red_gra_ty_x_flex= act_rgra_ty*red_gra_ty_x_mn
	red_gra_tz_x_flex= act_rgra_tz*red_gra_tz_x_mn
	red_gra_dx_y_flex= act_rgra_dx*red_gra_dx_y_mn
	red_gra_dy_y_flex= act_rgra_dy*red_gra_dy_y_mn
	red_gra_dz_y_flex= act_rgra_dz*red_gra_dz_y_mn
	red_gra_tx_y_flex= act_rgra_tx*red_gra_tx_y_mn
	red_gra_ty_y_flex= act_rgra_ty*red_gra_ty_y_mn
	red_gra_tz_y_flex= act_rgra_tz*red_gra_tz_y_mn
	
	red_mirr_dx_x_flex= act_mirr_dx*red_mirr_dx_x_mn
	red_mirr_dy_x_flex= act_mirr_dy*red_mirr_dy_x_mn
	red_mirr_dz_x_flex= act_mirr_dz*red_mirr_dz_x_mn
	red_mirr_tx_x_flex= act_mirr_tx*red_mirr_tx_x_mn
	red_mirr_ty_x_flex= act_mirr_ty*red_mirr_ty_x_mn
	red_mirr_tz_x_flex= act_mirr_tz*red_mirr_tz_x_mn
	red_mirr_dx_y_flex= act_mirr_dx*red_mirr_dx_y_mn
	red_mirr_dy_y_flex= act_mirr_dy*red_mirr_dy_y_mn
	red_mirr_dz_y_flex= act_mirr_dz*red_mirr_dz_y_mn
	red_mirr_tx_y_flex= act_mirr_tx*red_mirr_tx_y_mn
	red_mirr_ty_y_flex= act_mirr_ty*red_mirr_ty_y_mn
	red_mirr_tz_y_flex= act_mirr_tz*red_mirr_tz_y_mn

		
	blue_coll_dx_x_flex= act_coll_dx*blue_coll_dx_x_mn
	blue_coll_dy_x_flex= act_coll_dy*blue_coll_dy_x_mn
	blue_coll_dz_x_flex= act_coll_dz*blue_coll_dz_x_mn
	blue_coll_tx_x_flex= act_coll_tx*blue_coll_tx_x_mn
	blue_coll_ty_x_flex= act_coll_ty*blue_coll_ty_x_mn
	blue_coll_tz_x_flex= act_coll_tz*blue_coll_tz_x_mn
	blue_coll_dx_y_flex= act_coll_dx*blue_coll_dx_y_mn
	blue_coll_dy_y_flex= act_coll_dy*blue_coll_dy_y_mn
	blue_coll_dz_y_flex= act_coll_dz*blue_coll_dz_y_mn
	blue_coll_tx_y_flex= act_coll_tx*blue_coll_tx_y_mn
	blue_coll_ty_y_flex= act_coll_ty*blue_coll_ty_y_mn
	blue_coll_tz_y_flex= act_coll_tz*blue_coll_tz_y_mn
	
	blue_gra_dx_x_flex= act_bgra_dx*blue_gra_dx_x_mn
	blue_gra_dy_x_flex= act_bgra_dy*blue_gra_dy_x_mn
	blue_gra_dz_x_flex= act_bgra_dz*blue_gra_dz_x_mn
	blue_gra_tx_x_flex= act_bgra_tx*blue_gra_tx_x_mn
	blue_gra_ty_x_flex= act_bgra_ty*blue_gra_ty_x_mn
	blue_gra_tz_x_flex= act_bgra_tz*blue_gra_tz_x_mn
	blue_gra_dx_y_flex= act_bgra_dx*blue_gra_dx_y_mn
	blue_gra_dy_y_flex= act_bgra_dy*blue_gra_dy_y_mn
	blue_gra_dz_y_flex= act_bgra_dz*blue_gra_dz_y_mn
	blue_gra_tx_y_flex= act_bgra_tx*blue_gra_tx_y_mn
	blue_gra_ty_y_flex= act_bgra_ty*blue_gra_ty_y_mn
	blue_gra_tz_y_flex= act_bgra_tz*blue_gra_tz_y_mn
	
	blue_dchr_dx_x_flex= act_dchr_dx*blue_dchr_dx_x_mn
	blue_dchr_dy_x_flex= act_dchr_dy*blue_dchr_dy_x_mn
	blue_dchr_dz_x_flex= act_dchr_dz*blue_dchr_dz_x_mn
	blue_dchr_tx_x_flex= act_dchr_tx*blue_dchr_tx_x_mn
	blue_dchr_ty_x_flex= act_dchr_ty*blue_dchr_ty_x_mn
	blue_dchr_tz_x_flex= act_dchr_tz*blue_dchr_tz_x_mn
	blue_dchr_dx_y_flex= act_dchr_dx*blue_dchr_dx_y_mn
	blue_dchr_dy_y_flex= act_dchr_dy*blue_dchr_dy_y_mn
	blue_dchr_dz_y_flex= act_dchr_dz*blue_dchr_dz_y_mn
	blue_dchr_tx_y_flex= act_dchr_tx*blue_dchr_tx_y_mn
	blue_dchr_ty_y_flex= act_dchr_ty*blue_dchr_ty_y_mn
	blue_dchr_tz_y_flex= act_dchr_tz*blue_dchr_tz_y_mn
	
	
        act_coll_dx = float(self.tab3.myGrid.GetCellValue(0,0))
	act_coll_dy = float(self.tab3.myGrid.GetCellValue(0,1))
	act_coll_dz = float(self.tab3.myGrid.GetCellValue(0,2))
	act_coll_tx = float(self.tab3.myGrid.GetCellValue(0,3))
	act_coll_ty = float(self.tab3.myGrid.GetCellValue(0,4))
	act_coll_tz = float(self.tab3.myGrid.GetCellValue(0,5))
	
	act_mirr_dx = float(self.tab3.myGrid.GetCellValue(1,0))
	act_mirr_dy = float(self.tab3.myGrid.GetCellValue(1,1))
	act_mirr_dz = float(self.tab3.myGrid.GetCellValue(1,2))
	act_mirr_tx = float(self.tab3.myGrid.GetCellValue(1,3))
	act_mirr_ty = float(self.tab3.myGrid.GetCellValue(1,4))
	act_mirr_tz = float(self.tab3.myGrid.GetCellValue(1,5))
	
	act_dchr_dx = float(self.tab3.myGrid.GetCellValue(2,0))
	act_dchr_dy = float(self.tab3.myGrid.GetCellValue(2,1))
	act_dchr_dz = float(self.tab3.myGrid.GetCellValue(2,2))
	act_dchr_tx = float(self.tab3.myGrid.GetCellValue(2,3))
	act_dchr_ty = float(self.tab3.myGrid.GetCellValue(2,4))
	act_dchr_tz = float(self.tab3.myGrid.GetCellValue(2,5))
	
	act_rgra_dx = float(self.tab3.myGrid.GetCellValue(3,0))
	act_rgra_dy = float(self.tab3.myGrid.GetCellValue(3,1))
	act_rgra_dz = float(self.tab3.myGrid.GetCellValue(3,2))
	act_rgra_tx = float(self.tab3.myGrid.GetCellValue(3,3))
	act_rgra_ty = float(self.tab3.myGrid.GetCellValue(3,4))
	act_rgra_tz = float(self.tab3.myGrid.GetCellValue(3,5))
	
	act_bgra_dx = float(self.tab3.myGrid.GetCellValue(4,0))
	act_bgra_dy = float(self.tab3.myGrid.GetCellValue(4,1))
	act_bgra_dz = float(self.tab3.myGrid.GetCellValue(4,2))
	act_bgra_tx = float(self.tab3.myGrid.GetCellValue(4,3))
	act_bgra_ty = float(self.tab3.myGrid.GetCellValue(4,4))
	act_bgra_tz = float(self.tab3.myGrid.GetCellValue(4,5))

	det_red_shiftx = float(self.tab3.myGrid.GetCellValue(6,0))*1e-3
	det_red_shifty = float(self.tab3.myGrid.GetCellValue(6,1))*1e-3
	det_blue_shiftx = float(self.tab3.myGrid.GetCellValue(6,2))*1e-3
	det_blue_shifty = float(self.tab3.myGrid.GetCellValue(6,3))*1e-3			
	det_red_rot = np.radians(float(self.tab3.myGrid.GetCellValue(8,0)))
	det_blue_rot = np.radians(float(self.tab3.myGrid.GetCellValue(8,1)))
	
	red_coll_dx_x_comp= act_coll_dx*red_coll_dx_x_mn
	red_coll_dy_x_comp= act_coll_dy*red_coll_dy_x_mn
	red_coll_dz_x_comp= act_coll_dz*red_coll_dz_x_mn
	red_coll_tx_x_comp= act_coll_tx*red_coll_tx_x_mn
	red_coll_ty_x_comp= act_coll_ty*red_coll_ty_x_mn
	red_coll_tz_x_comp= act_coll_tz*red_coll_tz_x_mn
	red_coll_dx_y_comp= act_coll_dx*red_coll_dx_y_mn
	red_coll_dy_y_comp= act_coll_dy*red_coll_dy_y_mn
	red_coll_dz_y_comp= act_coll_dz*red_coll_dz_y_mn
	red_coll_tx_y_comp= act_coll_tx*red_coll_tx_y_mn
	red_coll_ty_y_comp= act_coll_ty*red_coll_ty_y_mn
	red_coll_tz_y_comp= act_coll_tz*red_coll_tz_y_mn
	
	red_det_x_comp= det_red_shiftx
	red_det_y_comp= det_red_shifty

	
	red_mirr_dx_x_comp= act_mirr_dx*red_mirr_dx_x_mn
	red_mirr_dy_x_comp= act_mirr_dy*red_mirr_dy_x_mn
	red_mirr_dz_x_comp= act_mirr_dz*red_mirr_dz_x_mn
	red_mirr_tx_x_comp= act_mirr_tx*red_mirr_tx_x_mn
	red_mirr_ty_x_comp= act_mirr_ty*red_mirr_ty_x_mn
	red_mirr_tz_x_comp= act_mirr_tz*red_mirr_tz_x_mn
	red_mirr_dx_y_comp= act_mirr_dx*red_mirr_dx_y_mn
	red_mirr_dy_y_comp= act_mirr_dy*red_mirr_dy_y_mn
	red_mirr_dz_y_comp= act_mirr_dz*red_mirr_dz_y_mn
	red_mirr_tx_y_comp= act_mirr_tx*red_mirr_tx_y_mn
	red_mirr_ty_y_comp= act_mirr_ty*red_mirr_ty_y_mn
	red_mirr_tz_y_comp= act_mirr_tz*red_mirr_tz_y_mn


	blue_coll_dx_x_comp= act_coll_dx*blue_coll_dx_x_mn
	blue_coll_dy_x_comp= act_coll_dy*blue_coll_dy_x_mn
	blue_coll_dz_x_comp= act_coll_dz*blue_coll_dz_x_mn
	blue_coll_tx_x_comp= act_coll_tx*blue_coll_tx_x_mn
	blue_coll_ty_x_comp= act_coll_ty*blue_coll_ty_x_mn
	blue_coll_tz_x_comp= act_coll_tz*blue_coll_tz_x_mn
	blue_coll_dx_y_comp= act_coll_dx*blue_coll_dx_y_mn
	blue_coll_dy_y_comp= act_coll_dy*blue_coll_dy_y_mn
	blue_coll_dz_y_comp= act_coll_dz*blue_coll_dz_y_mn
	blue_coll_tx_y_comp= act_coll_tx*blue_coll_tx_y_mn
	blue_coll_ty_y_comp= act_coll_ty*blue_coll_ty_y_mn
	blue_coll_tz_y_comp= act_coll_tz*blue_coll_tz_y_mn
	
	blue_det_x_comp= det_blue_shiftx
	blue_det_y_comp= det_blue_shifty

	
	blue_dchr_dx_x_comp= act_dchr_dx*blue_dchr_dx_x_mn
	blue_dchr_dy_x_comp= act_dchr_dy*blue_dchr_dy_x_mn
	blue_dchr_dz_x_comp= act_dchr_dz*blue_dchr_dz_x_mn
	blue_dchr_tx_x_comp= act_dchr_tx*blue_dchr_tx_x_mn
	blue_dchr_ty_x_comp= act_dchr_ty*blue_dchr_ty_x_mn
	blue_dchr_tz_x_comp= act_dchr_tz*blue_dchr_tz_x_mn
	blue_dchr_dx_y_comp= act_dchr_dx*blue_dchr_dx_y_mn
	blue_dchr_dy_y_comp= act_dchr_dy*blue_dchr_dy_y_mn
	blue_dchr_dz_y_comp= act_dchr_dz*blue_dchr_dz_y_mn
	blue_dchr_tx_y_comp= act_dchr_tx*blue_dchr_tx_y_mn
	blue_dchr_ty_y_comp= act_dchr_ty*blue_dchr_ty_y_mn
	blue_dchr_tz_y_comp= act_dchr_tz*blue_dchr_tz_y_mn




        red_flex_x =red_coll_dx_x_flex+red_coll_dy_x_flex+red_coll_dz_x_flex+red_coll_tx_x_flex+red_coll_ty_x_flex+red_coll_tz_x_flex+red_gra_dx_x_flex+red_gra_dy_x_flex+red_gra_dz_x_flex+red_gra_tx_x_flex+red_gra_ty_x_flex+red_gra_tz_x_flex+ red_mirr_dx_x_flex+red_mirr_dy_x_flex+red_mirr_dz_x_flex+red_mirr_tx_x_flex+red_mirr_ty_x_flex+red_mirr_tz_x_flex
        red_flex_y =red_coll_dx_y_flex+red_coll_dy_y_flex+red_coll_dz_y_flex+red_coll_tx_y_flex+red_coll_ty_y_flex+red_coll_tz_y_flex+red_gra_dx_y_flex+red_gra_dy_y_flex+red_gra_dz_y_flex+red_gra_tx_y_flex+red_gra_ty_y_flex+red_gra_tz_y_flex+ red_mirr_dx_y_flex+red_mirr_dy_y_flex+red_mirr_dz_y_flex+red_mirr_tx_y_flex+red_mirr_ty_y_flex+red_mirr_tz_y_flex
        red_comp_x =red_coll_dx_x_comp+red_coll_dy_x_comp+red_coll_dz_x_comp+red_coll_tx_x_comp+red_coll_ty_x_comp+red_coll_tz_x_comp+red_det_x_comp+ red_mirr_dx_x_comp+red_mirr_dy_x_comp+red_mirr_dz_x_comp+red_mirr_tx_x_comp+red_mirr_ty_x_comp+red_mirr_tz_x_comp
	red_comp_y=red_coll_dx_y_comp+red_coll_dy_y_comp+red_coll_dz_y_comp+red_coll_tx_y_comp+red_coll_ty_y_comp+red_coll_tz_y_comp+red_mirr_dx_y_comp+red_mirr_dy_y_comp+red_mirr_dz_y_comp+red_mirr_tx_y_comp+red_mirr_ty_y_comp+red_mirr_tz_y_comp+red_det_y_comp

        blue_flex_x = blue_coll_dx_x_flex+blue_coll_dy_x_flex+blue_coll_dz_x_flex+blue_coll_tx_x_flex+blue_coll_ty_x_flex+blue_coll_tz_x_flex+blue_gra_dx_x_flex+blue_gra_dy_x_flex+blue_gra_dz_x_flex+blue_gra_tx_x_flex+blue_gra_ty_x_flex+blue_gra_tz_x_flex+ blue_dchr_dx_x_flex+blue_dchr_dy_x_flex+blue_dchr_dz_x_flex+blue_dchr_tx_x_flex+blue_dchr_ty_x_flex+blue_dchr_tz_x_flex
        blue_flex_y = blue_coll_dx_y_flex+blue_coll_dy_y_flex+blue_coll_dz_y_flex+blue_coll_tx_y_flex+blue_coll_ty_y_flex+blue_coll_tz_y_flex+blue_gra_dx_y_flex+blue_gra_dy_y_flex+blue_gra_dz_y_flex+blue_gra_tx_y_flex+blue_gra_ty_y_flex+blue_gra_tz_y_flex+ blue_dchr_dx_y_flex+blue_dchr_dy_y_flex+blue_dchr_dz_y_flex+blue_dchr_tx_y_flex+blue_dchr_ty_y_flex+blue_dchr_tz_y_flex
        blue_comp_x = blue_coll_dx_x_comp+blue_coll_dy_x_comp+blue_coll_dz_x_comp+blue_coll_tx_x_comp+blue_coll_ty_x_comp+blue_coll_tz_x_comp+blue_det_x_comp+ blue_dchr_dx_x_comp+blue_dchr_dy_x_comp+blue_dchr_dz_x_comp+blue_dchr_tx_x_comp+blue_dchr_ty_x_comp+blue_dchr_tz_x_comp
        blue_comp_y = blue_coll_dx_y_comp+blue_coll_dy_y_comp+blue_coll_dz_y_comp+blue_coll_tx_y_comp+blue_coll_ty_y_comp+blue_coll_tz_y_comp+blue_dchr_dx_y_comp+blue_dchr_dy_y_comp+blue_dchr_dz_y_comp+blue_dchr_tx_y_comp+blue_dchr_ty_y_comp+blue_dchr_tz_y_comp+blue_det_y_comp







	
	self.panell.drawanalyzered()
	self.panell2.drawanalyzeblue()
	

    def onFlexButton(self, event):
	global redo_x,redo_y,redo_coll_dx_x,redo_coll_dy_x,redo_coll_dz_x,redo_coll_dx_y,redo_coll_dy_y,redo_coll_dz_y,redo_coll_tx_x,redo_coll_ty_x,redo_coll_tz_x,redo_coll_tx_y,redo_coll_ty_y,redo_coll_tz_y,redo_x1,redo_y1,redo_x2,redo_y2
	global blueo_x,blueo_y,blueo_coll_dx_x,blueo_coll_dy_x,blueo_coll_dz_x,blueo_coll_dx_y,blueo_coll_dy_y,blueo_coll_dz_y,blueo_coll_tx_x,blueo_coll_ty_x,blueo_coll_tz_x,blueo_coll_tx_y,blueo_coll_ty_y,blueo_coll_tz_y,blueo_x1,blueo_y1,blueo_x2,blueo_y2
        global redo_x,redo_y,redo_mirr_dx_x,redo_mirr_dy_x,redo_mirr_dz_x,redo_mirr_dx_y,redo_mirr_dy_y,redo_mirr_dz_y,redo_mirr_tx_x,redo_mirr_ty_x,redo_mirr_tz_x,redo_mirr_tx_y,redo_mirr_ty_y,redo_mirr_tz_y
	global blueo_x,blueo_y,blueo_dchr_dx_x,blueo_dchr_dy_x,blueo_dchr_dz_x,blueo_dchr_dx_y,blueo_dchr_dy_y,blueo_dchr_dz_y,blueo_dchr_tx_x,blueo_dchr_ty_x,blueo_dchr_tz_x,blueo_dchr_tx_y,blueo_dchr_ty_y,blueo_dchr_tz_y
        global redo_x,redo_y,redo_gra_dx_x,redo_gra_dy_x,redo_gra_dz_x,redo_gra_dx_y,redo_gra_dy_y,redo_gra_dz_y,redo_gra_tx_x,redo_gra_ty_x,redo_gra_tz_x,redo_gra_tx_y,redo_gra_ty_y,redo_gra_tz_y
	global blueo_x,blueo_y,blueo_gra_dx_x,blueo_gra_dy_x,blueo_gra_dz_x,blueo_gra_dx_y,blueo_gra_dy_y,blueo_gra_dz_y,blueo_gra_tx_x,blueo_gra_ty_x,blueo_gra_tz_x,blueo_gra_tx_y,blueo_gra_ty_y,blueo_gra_tz_y
	
	global red_x,red_y,red_coll_dx_x,red_coll_dy_x,red_coll_dz_x,red_coll_dx_y,red_coll_dy_y,red_coll_dz_y,red_coll_tx_x,red_coll_ty_x,red_coll_tz_x,red_coll_tx_y,red_coll_ty_y,red_coll_tz_y,red_x1,red_y1,red_x2,red_y2
	global blue_x,blue_y,blue_coll_dx_x,blue_coll_dy_x,blue_coll_dz_x,blue_coll_dx_y,blue_coll_dy_y,blue_coll_dz_y,blue_coll_tx_x,blue_coll_ty_x,blue_coll_tz_x,blue_coll_tx_y,blue_coll_ty_y,blue_coll_tz_y,blue_x1,blue_y1,blue_x2,blue_y2
        global red_x,red_y,red_mirr_dx_x,red_mirr_dy_x,red_mirr_dz_x,red_mirr_dx_y,red_mirr_dy_y,red_mirr_dz_y,red_mirr_tx_x,red_mirr_ty_x,red_mirr_tz_x,red_mirr_tx_y,red_mirr_ty_y,red_mirr_tz_y
	global blue_x,blue_y,blue_dchr_dx_x,blue_dchr_dy_x,blue_dchr_dz_x,blue_dchr_dx_y,blue_dchr_dy_y,blue_dchr_dz_y,blue_dchr_tx_x,blue_dchr_ty_x,blue_dchr_tz_x,blue_dchr_tx_y,blue_dchr_ty_y,blue_dchr_tz_y
        global red_x,red_y,red_gra_dx_x,red_gra_dy_x,red_gra_dz_x,red_gra_dx_y,red_gra_dy_y,red_gra_dz_y,red_gra_tx_x,red_gra_ty_x,red_gra_tz_x,red_gra_tx_y,red_gra_ty_y,red_gra_tz_y
	global blue_x,blue_y,blue_gra_dx_x,blue_gra_dy_x,blue_gra_dz_x,blue_gra_dx_y,blue_gra_dy_y,blue_gra_dz_y,blue_gra_tx_x,blue_gra_ty_x,blue_gra_tz_x,blue_gra_tx_y,blue_gra_ty_y,blue_gra_tz_y
	
	global redo1_x,redo1_y,redo1_coll_dx_x,redo1_coll_dy_x,redo1_coll_dz_x,redo1_coll_dx_y,redo1_coll_dy_y,redo1_coll_dz_y,redo1_coll_tx_x,redo1_coll_ty_x,redo1_coll_tz_x,redo1_coll_tx_y,redo1_coll_ty_y,redo1_coll_tz_y,redo1_x1,redo1_y1,redo1_x2,redo1_y2
	global blueo1_x,blueo1_y,blueo1_coll_dx_x,blueo1_coll_dy_x,blueo1_coll_dz_x,blueo1_coll_dx_y,blueo1_coll_dy_y,blueo1_coll_dz_y,blueo1_coll_tx_x,blueo1_coll_ty_x,blueo1_coll_tz_x,blueo1_coll_tx_y,blueo1_coll_ty_y,blueo1_coll_tz_y,blueo1_x1,blueo1_y1,blueo1_x2,blueo1_y2
        global redo1_x,redo1_y,redo1_mirr_dx_x,redo1_mirr_dy_x,redo1_mirr_dz_x,redo1_mirr_dx_y,redo1_mirr_dy_y,redo1_mirr_dz_y,redo1_mirr_tx_x,redo1_mirr_ty_x,redo1_mirr_tz_x,redo1_mirr_tx_y,redo1_mirr_ty_y,redo1_mirr_tz_y
	global blueo1_x,blueo1_y,blueo1_dchr_dx_x,blueo1_dchr_dy_x,blueo1_dchr_dz_x,blueo1_dchr_dx_y,blueo1_dchr_dy_y,blueo1_dchr_dz_y,blueo1_dchr_tx_x,blueo1_dchr_ty_x,blueo1_dchr_tz_x,blueo1_dchr_tx_y,blueo1_dchr_ty_y,blueo1_dchr_tz_y
        global redo1_x,redo1_y,redo1_gra_dx_x,redo1_gra_dy_x,redo1_gra_dz_x,redo1_gra_dx_y,redo1_gra_dy_y,redo1_gra_dz_y,redo1_gra_tx_x,redo1_gra_ty_x,redo1_gra_tz_x,redo1_gra_tx_y,redo1_gra_ty_y,redo1_gra_tz_y
	global blueo1_x,blueo1_y,blueo1_gra_dx_x,blueo1_gra_dy_x,blueo1_gra_dz_x,blueo1_gra_dx_y,blueo1_gra_dy_y,blueo1_gra_dz_y,blueo1_gra_tx_x,blueo1_gra_ty_x,blueo1_gra_tz_x,blueo1_gra_tx_y,blueo1_gra_ty_y,blueo1_gra_tz_y
		


	self.act_coll_dx = float(self.tab2.myGrid.GetCellValue(0,0))
	self.act_coll_dy = float(self.tab2.myGrid.GetCellValue(0,1))
	self.act_coll_dz = float(self.tab2.myGrid.GetCellValue(0,2))
	self.act_coll_tx = float(self.tab2.myGrid.GetCellValue(0,3))
	self.act_coll_ty = float(self.tab2.myGrid.GetCellValue(0,4))
	self.act_coll_tz = float(self.tab2.myGrid.GetCellValue(0,5))
	
	self.act_mirr_dx = float(self.tab2.myGrid.GetCellValue(1,0))
	self.act_mirr_dy = float(self.tab2.myGrid.GetCellValue(1,1))
	self.act_mirr_dz = float(self.tab2.myGrid.GetCellValue(1,2))
	self.act_mirr_tx = float(self.tab2.myGrid.GetCellValue(1,3))
	self.act_mirr_ty = float(self.tab2.myGrid.GetCellValue(1,4))
	self.act_mirr_tz = float(self.tab2.myGrid.GetCellValue(1,5))
	
	self.act_dchr_dx = float(self.tab2.myGrid.GetCellValue(2,0))
	self.act_dchr_dy = float(self.tab2.myGrid.GetCellValue(2,1))
	self.act_dchr_dz = float(self.tab2.myGrid.GetCellValue(2,2))
	self.act_dchr_tx = float(self.tab2.myGrid.GetCellValue(2,3))
	self.act_dchr_ty = float(self.tab2.myGrid.GetCellValue(2,4))
	self.act_dchr_tz = float(self.tab2.myGrid.GetCellValue(2,5))
	
	self.act_rgra_dx = float(self.tab2.myGrid.GetCellValue(3,0))
	self.act_rgra_dy = float(self.tab2.myGrid.GetCellValue(3,1))
	self.act_rgra_dz = float(self.tab2.myGrid.GetCellValue(3,2))
	self.act_rgra_tx = float(self.tab2.myGrid.GetCellValue(3,3))
	self.act_rgra_ty = float(self.tab2.myGrid.GetCellValue(3,4))
	self.act_rgra_tz = float(self.tab2.myGrid.GetCellValue(3,5))
	
	self.act_bgra_dx = float(self.tab2.myGrid.GetCellValue(4,0))
	self.act_bgra_dy = float(self.tab2.myGrid.GetCellValue(4,1))
	self.act_bgra_dz = float(self.tab2.myGrid.GetCellValue(4,2))
	self.act_bgra_tx = float(self.tab2.myGrid.GetCellValue(4,3))
	self.act_bgra_ty = float(self.tab2.myGrid.GetCellValue(4,4))
	self.act_bgra_tz = float(self.tab2.myGrid.GetCellValue(4,5))

	self.det_red_shiftx = float(self.tab2.myGrid.GetCellValue(6,0))*1e-3
	self.det_red_shifty = float(self.tab2.myGrid.GetCellValue(6,1))*1e-3
	self.det_blue_shiftx = float(self.tab2.myGrid.GetCellValue(6,2))*1e-3
	self.det_blue_shifty = float(self.tab2.myGrid.GetCellValue(6,3))*1e-3			
	self.det_red_rot = np.radians(float(self.tab2.myGrid.GetCellValue(8,0)))
	self.det_blue_rot = np.radians(float(self.tab2.myGrid.GetCellValue(8,1)))
	
	
	fullfield=self.tab1.cb.GetValue()
	if fullfield==True:
	    redo_x= copy.copy(redo1_x)
	    redo_y= copy.copy(redo1_y)
	    redo_coll_dx_x= copy.copy(redo1_coll_dx_x)
	    redo_coll_dx_y= copy.copy(redo1_coll_dx_y)
	    redo_coll_dy_x= copy.copy(redo1_coll_dy_x)
	    redo_coll_dy_y= copy.copy(redo1_coll_dy_y)
	    redo_coll_dz_x= copy.copy(redo1_coll_dz_x)
	    redo_coll_dz_y= copy.copy(redo1_coll_dz_y)
	    redo_coll_tx_x= copy.copy(redo1_coll_tx_x)
	    redo_coll_tx_y= copy.copy(redo1_coll_tx_y)
	    redo_coll_ty_x= copy.copy(redo1_coll_ty_x)
	    redo_coll_ty_y= copy.copy(redo1_coll_ty_y)
	    redo_coll_tz_x= copy.copy(redo1_coll_tz_x)
	    redo_coll_tz_y= copy.copy(redo1_coll_tz_y)
	    
	    redo_mirr_dx_x= copy.copy(redo1_mirr_dx_x)
	    redo_mirr_dx_y= copy.copy(redo1_mirr_dx_y)
	    redo_mirr_dy_x= copy.copy(redo1_mirr_dy_x)
	    redo_mirr_dy_y= copy.copy(redo1_mirr_dy_y)
	    redo_mirr_dz_x= copy.copy(redo1_mirr_dz_x)
	    redo_mirr_dz_y= copy.copy(redo1_mirr_dz_y)
	    redo_mirr_tx_x= copy.copy(redo1_mirr_tx_x)
	    redo_mirr_tx_y= copy.copy(redo1_mirr_tx_y)
	    redo_mirr_ty_x= copy.copy(redo1_mirr_ty_x)
	    redo_mirr_ty_y= copy.copy(redo1_mirr_ty_y)
	    redo_mirr_tz_x= copy.copy(redo1_mirr_tz_x)
	    redo_mirr_tz_y= copy.copy(redo1_mirr_tz_y)
	    
	    redo_gra_dx_x= copy.copy(redo1_gra_dx_x)
	    redo_gra_dx_y= copy.copy(redo1_gra_dx_y)
	    redo_gra_dy_x= copy.copy(redo1_gra_dy_x)
	    redo_gra_dy_y= copy.copy(redo1_gra_dy_y)
	    redo_gra_dz_x= copy.copy(redo1_gra_dz_x)
	    redo_gra_dz_y= copy.copy(redo1_gra_dz_y)
	    redo_gra_tx_x= copy.copy(redo1_gra_tx_x)
	    redo_gra_tx_y= copy.copy(redo1_gra_tx_y)
	    redo_gra_ty_x= copy.copy(redo1_gra_ty_x)
	    redo_gra_ty_y= copy.copy(redo1_gra_ty_y)
	    redo_gra_tz_x= copy.copy(redo1_gra_tz_x)
	    redo_gra_tz_y= copy.copy(redo1_gra_tz_y)
	    
	    blueo_x= copy.copy(blueo1_x)
	    blueo_y= copy.copy(blueo1_y)
	    blueo_coll_dx_x= copy.copy(blueo1_coll_dx_x)
	    blueo_coll_dx_y= copy.copy(blueo1_coll_dx_y)
	    blueo_coll_dy_x= copy.copy(blueo1_coll_dy_x)
	    blueo_coll_dy_y= copy.copy(blueo1_coll_dy_y)
	    blueo_coll_dz_x= copy.copy(blueo1_coll_dz_x)
	    blueo_coll_dz_y= copy.copy(blueo1_coll_dz_y)
	    blueo_coll_tx_x= copy.copy(blueo1_coll_tx_x)
	    blueo_coll_tx_y= copy.copy(blueo1_coll_tx_y)
	    blueo_coll_ty_x= copy.copy(blueo1_coll_ty_x)
	    blueo_coll_ty_y= copy.copy(blueo1_coll_ty_y)
	    blueo_coll_tz_x= copy.copy(blueo1_coll_tz_x)
	    blueo_coll_tz_y= copy.copy(blueo1_coll_tz_y)
	    
	    blueo_dchr_dx_x= copy.copy(blueo1_dchr_dx_x)
	    blueo_dchr_dx_y= copy.copy(blueo1_dchr_dx_y)
	    blueo_dchr_dy_x= copy.copy(blueo1_dchr_dy_x)
	    blueo_dchr_dy_y= copy.copy(blueo1_dchr_dy_y)
	    blueo_dchr_dz_x= copy.copy(blueo1_dchr_dz_x)
	    blueo_dchr_dz_y= copy.copy(blueo1_dchr_dz_y)
	    blueo_dchr_tx_x= copy.copy(blueo1_dchr_tx_x)
	    blueo_dchr_tx_y= copy.copy(blueo1_dchr_tx_y)
	    blueo_dchr_ty_x= copy.copy(blueo1_dchr_ty_x)
	    blueo_dchr_ty_y= copy.copy(blueo1_dchr_ty_y)
	    blueo_dchr_tz_x= copy.copy(blueo1_dchr_tz_x)
	    blueo_dchr_tz_y= copy.copy(blueo1_dchr_tz_y)
	    
	    blueo_gra_dx_x= copy.copy(blueo1_gra_dx_x)
	    blueo_gra_dx_y= copy.copy(blueo1_gra_dx_y)
	    blueo_gra_dy_x= copy.copy(blueo1_gra_dy_x)
	    blueo_gra_dy_y= copy.copy(blueo1_gra_dy_y)
	    blueo_gra_dz_x= copy.copy(blueo1_gra_dz_x)
	    blueo_gra_dz_y= copy.copy(blueo1_gra_dz_y)
	    blueo_gra_tx_x= copy.copy(blueo1_gra_tx_x)
	    blueo_gra_tx_y= copy.copy(blueo1_gra_tx_y)
	    blueo_gra_ty_x= copy.copy(blueo1_gra_ty_x)
	    blueo_gra_ty_y= copy.copy(blueo1_gra_ty_y)
	    blueo_gra_tz_x= copy.copy(blueo1_gra_tz_x)
	    blueo_gra_tz_y= copy.copy(blueo1_gra_tz_y)

	if fullfield==False:
	    redo_x= copy.copy(red_x)
	    redo_y= copy.copy(red_y)
	    redo_coll_dx_x= copy.copy(red_coll_dx_x)
	    redo_coll_dx_y= copy.copy(red_coll_dx_y)
	    redo_coll_dy_x= copy.copy(red_coll_dy_x)
	    redo_coll_dy_y= copy.copy(red_coll_dy_y)
	    redo_coll_dz_x= copy.copy(red_coll_dz_x)
	    redo_coll_dz_y= copy.copy(red_coll_dz_y)
	    redo_coll_tx_x= copy.copy(red_coll_tx_x)
	    redo_coll_tx_y= copy.copy(red_coll_tx_y)
	    redo_coll_ty_x= copy.copy(red_coll_ty_x)
	    redo_coll_ty_y= copy.copy(red_coll_ty_y)
	    redo_coll_tz_x= copy.copy(red_coll_tz_x)
	    redo_coll_tz_y= copy.copy(red_coll_tz_y)
	    
	    redo_mirr_dx_x= copy.copy(red_mirr_dx_x)
	    redo_mirr_dx_y= copy.copy(red_mirr_dx_y)
	    redo_mirr_dy_x= copy.copy(red_mirr_dy_x)
	    redo_mirr_dy_y= copy.copy(red_mirr_dy_y)
	    redo_mirr_dz_x= copy.copy(red_mirr_dz_x)
	    redo_mirr_dz_y= copy.copy(red_mirr_dz_y)
	    redo_mirr_tx_x= copy.copy(red_mirr_tx_x)
	    redo_mirr_tx_y= copy.copy(red_mirr_tx_y)
	    redo_mirr_ty_x= copy.copy(red_mirr_ty_x)
	    redo_mirr_ty_y= copy.copy(red_mirr_ty_y)
	    redo_mirr_tz_x= copy.copy(red_mirr_tz_x)
	    redo_mirr_tz_y= copy.copy(red_mirr_tz_y)
	    
	    redo_gra_dx_x= copy.copy(red_gra_dx_x)
	    redo_gra_dx_y= copy.copy(red_gra_dx_y)
	    redo_gra_dy_x= copy.copy(red_gra_dy_x)
	    redo_gra_dy_y= copy.copy(red_gra_dy_y)
	    redo_gra_dz_x= copy.copy(red_gra_dz_x)
	    redo_gra_dz_y= copy.copy(red_gra_dz_y)
	    redo_gra_tx_x= copy.copy(red_gra_tx_x)
	    redo_gra_tx_y= copy.copy(red_gra_tx_y)
	    redo_gra_ty_x= copy.copy(red_gra_ty_x)
	    redo_gra_ty_y= copy.copy(red_gra_ty_y)
	    redo_gra_tz_x= copy.copy(red_gra_tz_x)
	    redo_gra_tz_y= copy.copy(red_gra_tz_y)
	    
	    blueo_x= copy.copy(blue_x)
	    blueo_y= copy.copy(blue_y)
	    blueo_coll_dx_x= copy.copy(blue_coll_dx_x)
	    blueo_coll_dx_y= copy.copy(blue_coll_dx_y)
	    blueo_coll_dy_x= copy.copy(blue_coll_dy_x)
	    blueo_coll_dy_y= copy.copy(blue_coll_dy_y)
	    blueo_coll_dz_x= copy.copy(blue_coll_dz_x)
	    blueo_coll_dz_y= copy.copy(blue_coll_dz_y)
	    blueo_coll_tx_x= copy.copy(blue_coll_tx_x)
	    blueo_coll_tx_y= copy.copy(blue_coll_tx_y)
	    blueo_coll_ty_x= copy.copy(blue_coll_ty_x)
	    blueo_coll_ty_y= copy.copy(blue_coll_ty_y)
	    blueo_coll_tz_x= copy.copy(blue_coll_tz_x)
	    blueo_coll_tz_y= copy.copy(blue_coll_tz_y)
	    
	    blueo_dchr_dx_x= copy.copy(blue_dchr_dx_x)
	    blueo_dchr_dx_y= copy.copy(blue_dchr_dx_y)
	    blueo_dchr_dy_x= copy.copy(blue_dchr_dy_x)
	    blueo_dchr_dy_y= copy.copy(blue_dchr_dy_y)
	    blueo_dchr_dz_x= copy.copy(blue_dchr_dz_x)
	    blueo_dchr_dz_y= copy.copy(blue_dchr_dz_y)
	    blueo_dchr_tx_x= copy.copy(blue_dchr_tx_x)
	    blueo_dchr_tx_y= copy.copy(blue_dchr_tx_y)
	    blueo_dchr_ty_x= copy.copy(blue_dchr_ty_x)
	    blueo_dchr_ty_y= copy.copy(blue_dchr_ty_y)
	    blueo_dchr_tz_x= copy.copy(blue_dchr_tz_x)
	    blueo_dchr_tz_y= copy.copy(blue_dchr_tz_y)
	    
	    blueo_gra_dx_x= copy.copy(blue_gra_dx_x)
	    blueo_gra_dx_y= copy.copy(blue_gra_dx_y)
	    blueo_gra_dy_x= copy.copy(blue_gra_dy_x)
	    blueo_gra_dy_y= copy.copy(blue_gra_dy_y)
	    blueo_gra_dz_x= copy.copy(blue_gra_dz_x)
	    blueo_gra_dz_y= copy.copy(blue_gra_dz_y)
	    blueo_gra_tx_x= copy.copy(blue_gra_tx_x)
	    blueo_gra_tx_y= copy.copy(blue_gra_tx_y)
	    blueo_gra_ty_x= copy.copy(blue_gra_ty_x)
	    blueo_gra_ty_y= copy.copy(blue_gra_ty_y)
	    blueo_gra_tz_x= copy.copy(blue_gra_tz_x)
	    blueo_gra_tz_y= copy.copy(blue_gra_tz_y)
		
	#print self.det_blue_rot
	redo_x1= redo_x+self.act_coll_dx*redo_coll_dx_x +self.act_coll_dy*redo_coll_dy_x	+ self.act_coll_dz*redo_coll_dz_x + self.act_coll_tx*redo_coll_tx_x +self.act_coll_ty*redo_coll_ty_x+ self.act_coll_tz*redo_coll_tz_x + \
		 self.act_mirr_dx*redo_mirr_dx_x +self.act_mirr_dy*redo_mirr_dy_x	+ self.act_mirr_dz*redo_mirr_dz_x + self.act_mirr_tx*redo_mirr_tx_x +self.act_mirr_ty*redo_mirr_ty_x+ self.act_mirr_tz*redo_mirr_tz_x + \
		 self.act_rgra_dx*redo_gra_dx_x +self.act_rgra_dy*redo_gra_dy_x	+ self.act_rgra_dz*redo_gra_dz_x + self.act_rgra_tx*redo_gra_tx_x +self.act_rgra_ty*redo_gra_ty_x+ self.act_rgra_tz*redo_gra_tz_x
	redo_y1= redo_y+self.act_coll_dx*redo_coll_dx_y +self.act_coll_dy*redo_coll_dy_y	+ self.act_coll_dz*redo_coll_dz_y + self.act_coll_tx*redo_coll_tx_y +self.act_coll_ty*redo_coll_ty_y+ self.act_coll_tz*redo_coll_tz_y + \
		 self.act_mirr_dx*redo_mirr_dx_y +self.act_mirr_dy*redo_mirr_dy_y	+ self.act_mirr_dz*redo_mirr_dz_y + self.act_mirr_tx*redo_mirr_tx_y +self.act_mirr_ty*redo_mirr_ty_y+ self.act_mirr_tz*redo_mirr_tz_y + \
		 self.act_rgra_dx*redo_gra_dx_y +self.act_rgra_dy*redo_gra_dy_y	+ self.act_rgra_dz*redo_gra_dz_y + self.act_rgra_tx*redo_gra_tx_y +self.act_rgra_ty*redo_gra_ty_y+ self.act_rgra_tz*redo_gra_tz_y

	blueo_x1= blueo_x+self.act_coll_dx*blueo_coll_dx_x +self.act_coll_dy*blueo_coll_dy_x	+ self.act_coll_dz*blueo_coll_dz_x + self.act_coll_tx*blueo_coll_tx_x +self.act_coll_ty*blueo_coll_ty_x+ self.act_coll_tz*blueo_coll_tz_x + \
		 self.act_dchr_dx*blueo_dchr_dx_x +self.act_dchr_dy*blueo_dchr_dy_x	+ self.act_dchr_dz*blueo_dchr_dz_x + self.act_dchr_tx*blueo_dchr_tx_x +self.act_dchr_ty*blueo_dchr_ty_x+ self.act_dchr_tz*blueo_dchr_tz_x + \
		 self.act_bgra_dx*blueo_gra_dx_x +self.act_bgra_dy*blueo_gra_dy_x	+ self.act_bgra_dz*blueo_gra_dz_x + self.act_bgra_tx*blueo_gra_tx_x +self.act_bgra_ty*blueo_gra_ty_x+ self.act_bgra_tz*blueo_gra_tz_x
	blueo_y1= blueo_y+self.act_coll_dx*blueo_coll_dx_y +self.act_coll_dy*blueo_coll_dy_y	+ self.act_coll_dz*blueo_coll_dz_y + self.act_coll_tx*blueo_coll_tx_y +self.act_coll_ty*blueo_coll_ty_y+ self.act_coll_tz*blueo_coll_tz_y + \
		  self.act_dchr_dx*blueo_dchr_dx_y +self.act_dchr_dy*blueo_dchr_dy_y	+ self.act_dchr_dz*blueo_dchr_dz_y + self.act_dchr_tx*blueo_dchr_tx_y +self.act_dchr_ty*blueo_dchr_ty_y+ self.act_dchr_tz*blueo_dchr_tz_y + \
		  self.act_bgra_dx*blueo_gra_dx_y +self.act_bgra_dy*blueo_gra_dy_y	+ self.act_bgra_dz*blueo_gra_dz_y + self.act_bgra_tx*blueo_gra_tx_y +self.act_bgra_ty*blueo_gra_ty_y+ self.act_bgra_tz*blueo_gra_tz_y

	red_x1= red_x+self.act_coll_dx*red_coll_dx_x +self.act_coll_dy*red_coll_dy_x	+ self.act_coll_dz*red_coll_dz_x + self.act_coll_tx*red_coll_tx_x +self.act_coll_ty*red_coll_ty_x+ self.act_coll_tz*red_coll_tz_x + \
		 self.act_mirr_dx*red_mirr_dx_x +self.act_mirr_dy*red_mirr_dy_x	+ self.act_mirr_dz*red_mirr_dz_x + self.act_mirr_tx*red_mirr_tx_x +self.act_mirr_ty*red_mirr_ty_x+ self.act_mirr_tz*red_mirr_tz_x + \
		 self.act_rgra_dx*red_gra_dx_x +self.act_rgra_dy*red_gra_dy_x	+ self.act_rgra_dz*red_gra_dz_x + self.act_rgra_tx*red_gra_tx_x +self.act_rgra_ty*red_gra_ty_x+ self.act_rgra_tz*red_gra_tz_x
	red_y1= red_y+self.act_coll_dx*red_coll_dx_y +self.act_coll_dy*red_coll_dy_y	+ self.act_coll_dz*red_coll_dz_y + self.act_coll_tx*red_coll_tx_y +self.act_coll_ty*red_coll_ty_y+ self.act_coll_tz*red_coll_tz_y + \
		 self.act_mirr_dx*red_mirr_dx_y +self.act_mirr_dy*red_mirr_dy_y	+ self.act_mirr_dz*red_mirr_dz_y + self.act_mirr_tx*red_mirr_tx_y +self.act_mirr_ty*red_mirr_ty_y+ self.act_mirr_tz*red_mirr_tz_y + \
		 self.act_rgra_dx*red_gra_dx_y +self.act_rgra_dy*red_gra_dy_y	+ self.act_rgra_dz*red_gra_dz_y + self.act_rgra_tx*red_gra_tx_y +self.act_rgra_ty*red_gra_ty_y+ self.act_rgra_tz*red_gra_tz_y

	blue_x1= blue_x+self.act_coll_dx*blue_coll_dx_x +self.act_coll_dy*blue_coll_dy_x	+ self.act_coll_dz*blue_coll_dz_x + self.act_coll_tx*blue_coll_tx_x +self.act_coll_ty*blue_coll_ty_x+ self.act_coll_tz*blue_coll_tz_x + \
		 self.act_dchr_dx*blue_dchr_dx_x +self.act_dchr_dy*blue_dchr_dy_x	+ self.act_dchr_dz*blue_dchr_dz_x + self.act_dchr_tx*blue_dchr_tx_x +self.act_dchr_ty*blue_dchr_ty_x+ self.act_dchr_tz*blue_dchr_tz_x + \
		 self.act_bgra_dx*blue_gra_dx_x +self.act_bgra_dy*blue_gra_dy_x	+ self.act_bgra_dz*blue_gra_dz_x + self.act_bgra_tx*blue_gra_tx_x +self.act_bgra_ty*blue_gra_ty_x+ self.act_bgra_tz*blue_gra_tz_x
	blue_y1= blue_y+self.act_coll_dx*blue_coll_dx_y +self.act_coll_dy*blue_coll_dy_y	+ self.act_coll_dz*blue_coll_dz_y + self.act_coll_tx*blue_coll_tx_y +self.act_coll_ty*blue_coll_ty_y+ self.act_coll_tz*blue_coll_tz_y + \
		  self.act_dchr_dx*blue_dchr_dx_y +self.act_dchr_dy*blue_dchr_dy_y+ self.act_dchr_dz*blue_dchr_dz_y + self.act_dchr_tx*blue_dchr_tx_y +self.act_dchr_ty*blue_dchr_ty_y+ self.act_dchr_tz*blue_dchr_tz_y + \
		  self.act_bgra_dx*blue_gra_dx_y +self.act_bgra_dy*blue_gra_dy_y+ self.act_bgra_dz*blue_gra_dz_y + self.act_bgra_tx*blue_gra_tx_y +self.act_bgra_ty*blue_gra_ty_y+ self.act_bgra_tz*blue_gra_tz_y
        #print 'first',red_x[1:13]
	rxrot2 = red_x1*np.cos(self.det_red_rot)-red_y1*np.sin(self.det_red_rot)
        ryrot2=  red_x1*np.sin(self.det_red_rot)+red_y1*np.cos(self.det_red_rot)
	
	bxrot2 = blue_x1*np.cos(self.det_blue_rot)-blue_y1*np.sin(self.det_blue_rot)
        byrot2=  blue_x1*np.sin(self.det_blue_rot)+blue_y1*np.cos(self.det_blue_rot)
	
	red_x1=rxrot2+self.det_red_shiftx
	red_y1=ryrot2+self.det_red_shifty
	blue_x1=bxrot2+self.det_blue_shiftx
	blue_y1=byrot2+self.det_blue_shifty
	
	rxrot2 = redo_x1*np.cos(self.det_red_rot)-redo_y1*np.sin(self.det_red_rot)
        ryrot2=  redo_x1*np.sin(self.det_red_rot)+redo_y1*np.cos(self.det_red_rot)
        #print 'last',red_x[1:13]
	bxrot2 = blueo_x1*np.cos(self.det_blue_rot)-blueo_y1*np.sin(self.det_blue_rot)
        byrot2=  blueo_x1*np.sin(self.det_blue_rot)+blueo_y1*np.cos(self.det_blue_rot)
	
	redo_x1=rxrot2+self.det_red_shiftx
	redo_y1=ryrot2+self.det_red_shifty
	blueo_x1=bxrot2+self.det_blue_shiftx
	blueo_y1=byrot2+self.det_blue_shifty
   
	#err=(sum(np.sqrt((redo_x-redo_x1)**2+(redo_y-redo_y1)**2))+ sum(np.sqrt((blueo_x-blueo_x1)**2+(blueo_y-blueo_y1)**2)))/(len(blueo_x)+len(redo_x))
	#print "flex :", err


	#print self.act_coll_dx
	self.panell.flexscale=int(self.tab2.selectFlexScale.GetValue())
	self.panell.drawredshift()
	self.panell2.flexscale=int(self.tab2.selectFlexScale.GetValue())
	self.panell2.drawblueshift()
	

    def onConfigSelect(self, event):
	global redo_x,redo_y,redo_coll_dx_x,redo_coll_dy_x,redo_coll_dz_x,redo_coll_dx_y,redo_coll_dy_y,redo_coll_dz_y,redo_coll_tx_x,redo_coll_ty_x,redo_coll_tz_x,redo_coll_tx_y,redo_coll_ty_y,redo_coll_tz_y,redo_x1,redo_y1,redo_x2,redo_y2
	global blueo_x,blueo_y,blueo_coll_dx_x,blueo_coll_dy_x,blueo_coll_dz_x,blueo_coll_dx_y,blueo_coll_dy_y,blueo_coll_dz_y,blueo_coll_tx_x,blueo_coll_ty_x,blueo_coll_tz_x,blueo_coll_tx_y,blueo_coll_ty_y,blueo_coll_tz_y,blueo_x1,blueo_y1,blueo_x2,blueo_y2
        global redo_x,redo_y,redo_mirr_dx_x,redo_mirr_dy_x,redo_mirr_dz_x,redo_mirr_dx_y,redo_mirr_dy_y,redo_mirr_dz_y,redo_mirr_tx_x,redo_mirr_ty_x,redo_mirr_tz_x,redo_mirr_tx_y,redo_mirr_ty_y,redo_mirr_tz_y
	global blueo_x,blueo_y,blueo_dchr_dx_x,blueo_dchr_dy_x,blueo_dchr_dz_x,blueo_dchr_dx_y,blueo_dchr_dy_y,blueo_dchr_dz_y,blueo_dchr_tx_x,blueo_dchr_ty_x,blueo_dchr_tz_x,blueo_dchr_tx_y,blueo_dchr_ty_y,blueo_dchr_tz_y
        global redo_x,redo_y,redo_gra_dx_x,redo_gra_dy_x,redo_gra_dz_x,redo_gra_dx_y,redo_gra_dy_y,redo_gra_dz_y,redo_gra_tx_x,redo_gra_ty_x,redo_gra_tz_x,redo_gra_tx_y,redo_gra_ty_y,redo_gra_tz_y
	global blueo_x,blueo_y,blueo_gra_dx_x,blueo_gra_dy_x,blueo_gra_dz_x,blueo_gra_dx_y,blueo_gra_dy_y,blueo_gra_dz_y,blueo_gra_tx_x,blueo_gra_ty_x,blueo_gra_tz_x,blueo_gra_tx_y,blueo_gra_ty_y,blueo_gra_tz_y
	
	global redo1_x,redo1_y,redo1_coll_dx_x,redo1_coll_dy_x,redo1_coll_dz_x,redo1_coll_dx_y,redo1_coll_dy_y,redo1_coll_dz_y,redo1_coll_tx_x,redo1_coll_ty_x,redo1_coll_tz_x,redo1_coll_tx_y,redo1_coll_ty_y,redo1_coll_tz_y,redo1_x1,redo1_y1,redo1_x2,redo1_y2
	global blueo1_x,blueo1_y,blueo1_coll_dx_x,blueo1_coll_dy_x,blueo1_coll_dz_x,blueo1_coll_dx_y,blueo1_coll_dy_y,blueo1_coll_dz_y,blueo1_coll_tx_x,blueo1_coll_ty_x,blueo1_coll_tz_x,blueo1_coll_tx_y,blueo1_coll_ty_y,blueo1_coll_tz_y,blueo1_x1,blueo1_y1,blueo1_x2,blueo1_y2
        global redo1_x,redo1_y,redo1_mirr_dx_x,redo1_mirr_dy_x,redo1_mirr_dz_x,redo1_mirr_dx_y,redo1_mirr_dy_y,redo1_mirr_dz_y,redo1_mirr_tx_x,redo1_mirr_ty_x,redo1_mirr_tz_x,redo1_mirr_tx_y,redo1_mirr_ty_y,redo1_mirr_tz_y
	global blueo1_x,blueo1_y,blueo1_dchr_dx_x,blueo1_dchr_dy_x,blueo1_dchr_dz_x,blueo1_dchr_dx_y,blueo1_dchr_dy_y,blueo1_dchr_dz_y,blueo1_dchr_tx_x,blueo1_dchr_ty_x,blueo1_dchr_tz_x,blueo1_dchr_tx_y,blueo1_dchr_ty_y,blueo1_dchr_tz_y
        global redo1_x,redo1_y,redo1_gra_dx_x,redo1_gra_dy_x,redo1_gra_dz_x,redo1_gra_dx_y,redo1_gra_dy_y,redo1_gra_dz_y,redo1_gra_tx_x,redo1_gra_ty_x,redo1_gra_tz_x,redo1_gra_tx_y,redo1_gra_ty_y,redo1_gra_tz_y
	global blueo1_x,blueo1_y,blueo1_gra_dx_x,blueo1_gra_dy_x,blueo1_gra_dz_x,blueo1_gra_dx_y,blueo1_gra_dy_y,blueo1_gra_dz_y,blueo1_gra_tx_x,blueo1_gra_ty_x,blueo1_gra_tz_x,blueo1_gra_tx_y,blueo1_gra_ty_y,blueo1_gra_tz_y
	global config	
	global red_x,red_y,red_coll_dx_x,red_coll_dy_x,red_coll_dz_x,red_coll_dx_y,red_coll_dy_y,red_coll_dz_y,red_coll_tx_x,red_coll_ty_x,red_coll_tz_x,red_coll_tx_y,red_coll_ty_y,red_coll_tz_y,red_x1,red_y1,red_x2,red_y2
	global blue_x,blue_y,blue_coll_dx_x,blue_coll_dy_x,blue_coll_dz_x,blue_coll_dx_y,blue_coll_dy_y,blue_coll_dz_y,blue_coll_tx_x,blue_coll_ty_x,blue_coll_tz_x,blue_coll_tx_y,blue_coll_ty_y,blue_coll_tz_y,blue_x1,blue_y1,blue_x2,blue_y2
        global red_x,red_y,red_mirr_dx_x,red_mirr_dy_x,red_mirr_dz_x,red_mirr_dx_y,red_mirr_dy_y,red_mirr_dz_y,red_mirr_tx_x,red_mirr_ty_x,red_mirr_tz_x,red_mirr_tx_y,red_mirr_ty_y,red_mirr_tz_y
	global blue_x,blue_y,blue_dchr_dx_x,blue_dchr_dy_x,blue_dchr_dz_x,blue_dchr_dx_y,blue_dchr_dy_y,blue_dchr_dz_y,blue_dchr_tx_x,blue_dchr_ty_x,blue_dchr_tz_x,blue_dchr_tx_y,blue_dchr_ty_y,blue_dchr_tz_y
        global red_x,red_y,red_gra_dx_x,red_gra_dy_x,red_gra_dz_x,red_gra_dx_y,red_gra_dy_y,red_gra_dz_y,red_gra_tx_x,red_gra_ty_x,red_gra_tz_x,red_gra_tx_y,red_gra_ty_y,red_gra_tz_y
	global blue_x,blue_y,blue_gra_dx_x,blue_gra_dy_x,blue_gra_dz_x,blue_gra_dx_y,blue_gra_dy_y,blue_gra_dz_y,blue_gra_tx_x,blue_gra_ty_x,blue_gra_tz_x,blue_gra_tx_y,blue_gra_ty_y,blue_gra_tz_y
        config =self.tab2.configcmbx.GetValue()		
	if config=='14out9prx':
	    print "Current Config : 14out9 paraxial camera" 	
	    red_coll_file='./14out9/red_coll_14.9.txt'
	    red_mirr_file='./14out9/red_mirr_14.9.txt'
	    red_gra_file='./14out9/red_gra_14.9.txt'
	    blue_coll_file='./14out9/blue_coll_14.9.txt'
	    blue_dchr_file='./14out9/blue_dchr_14.9.txt'
	    blue_gra_file='./14out9/blue_gra_14.9.txt'	
	if config=='15out6prx':
	    print "Current Config : 15out6 paraxial camera" 		
	    red_coll_file='./15out6/red_coll_15.6.txt'
	    red_mirr_file='./15out6/red_mirr_15.6.txt'
	    red_gra_file='./15out6/red_gra_15.6.txt'
	    blue_coll_file='./15out6/blue_coll_15.6.txt'
	    blue_dchr_file='./15out6/blue_dchr_15.6.txt'
	    blue_gra_file='./15out6/blue_gra_15.6.txt'	
	if config=='MOBIE':
	    print "Current Config : MOBIE" 		
	    red_coll_file='./mobie/red_coll.txt'
	    red_mirr_file='./mobie/red_mirr.txt'
	    red_gra_file='./mobie/red_gra.txt'
	    blue_coll_file='./mobie/blue_coll.txt'
	    blue_dchr_file='./mobie/blue_dchr.txt'
	    blue_gra_file='./mobie/blue_gra.txt'	
	if config=='15out6cam':
	    print "Current Config : 15out6 A7/BBA1 Camera" 		
	    red_coll_file='./15out6cam/red_coll_15out6cam.txt'
	    red_mirr_file='./15out6cam/red_mirr_15out6cam.txt'
	    red_gra_file='./15out6cam/red_gra_15out6cam.txt'
	    blue_coll_file='./15out6cam/blue_coll_15out6cam.txt'
	    blue_dchr_file='./15out6cam/blue_dchr_15out6cam.txt'
	    blue_gra_file='./15out6cam/blue_gra_15out6cam.txt'	
	    
	red_datalist= pd.read_csv(red_coll_file,sep='\s+')
	red_ofield=[]
	red_infield=[]
	
	
	rx=np.array(red_datalist.X)
	ry=np.array(red_datalist.Y)
	rele=red_datalist.Element
	rcomp=red_datalist.Peturbation
	rdelx=np.array(red_datalist.DelX)
	rdely=np.array(red_datalist.DelY)
	
	
	red_x=[]
	red_y=[]
	red_coll_dx_x=[]
	red_coll_dx_y=[]
	red_coll_dy_x=[]
	red_coll_dy_y=[]
	red_coll_dz_x=[]
	red_coll_dz_y=[]
	red_coll_tx_x=[]
	red_coll_tx_y=[]
	red_coll_ty_x=[]
	red_coll_ty_y=[]
	red_coll_tz_x=[]
	red_coll_tz_y=[]
	
	
	for i in range(len(rx)):
		if (i%6 == 0):
			red_x.append(rx[i])
		if (i%6 == 0):
			red_y.append(ry[i])
		if (i%6 == 0):
			red_coll_dx_x.append(rdelx[i])
		if (i%6 == 0):
			red_coll_dx_y.append(rdely[i])
		if (i%6 == 1):
			red_coll_dy_x.append(rdelx[i])
		if (i%6 == 1):
			red_coll_dy_y.append(rdely[i])
		if (i%6 == 2):
			red_coll_dz_x.append(rdelx[i])
		if (i%6 == 2):
			red_coll_dz_y.append(rdely[i])
		if (i%6 == 3):
			red_coll_tx_x.append(rdelx[i])
		if (i%6 == 3):
			red_coll_tx_y.append(rdely[i])
		if (i%6 == 4):
			red_coll_ty_x.append(rdelx[i])
		if (i%6 == 4):
			red_coll_ty_y.append(rdely[i])
		if (i%6 == 5):
			red_coll_tz_x.append(rdelx[i])
		if (i%6 == 5):
			red_coll_tz_y.append(rdely[i])
	
	
	red_x=numpy.array(red_x)
	red_y=numpy.array(red_y)
	red_coll_dx_x=numpy.array(red_coll_dx_x)*1e-6
	red_coll_dx_y=numpy.array(red_coll_dx_y)*1e-6
	red_coll_dy_x=numpy.array(red_coll_dy_x)*1e-6
	red_coll_dy_y=numpy.array(red_coll_dy_y)*1e-6
	red_coll_dz_x=numpy.array(red_coll_dz_x)*1e-6
	red_coll_dz_y=numpy.array(red_coll_dz_y)*1e-6
	red_coll_tx_x=numpy.array(red_coll_tx_x)*1e-6
	red_coll_tx_y=numpy.array(red_coll_tx_y)*1e-6
	red_coll_ty_x=numpy.array(red_coll_ty_x)*1e-6
	red_coll_ty_y=numpy.array(red_coll_ty_y)*1e-6
	red_coll_tz_x=numpy.array(red_coll_tz_x)*1e-6
	red_coll_tz_y=numpy.array(red_coll_tz_y)*1e-6
	
	
	red_datalist= pd.read_csv(red_mirr_file,sep='\s+')
	red_ofield=[]
	red_infield=[]
	
	
	rmele=red_datalist.Element
	rmcomp=red_datalist.Peturbation
	rmdelx=np.array(red_datalist.DelX)
	rmdely=np.array(red_datalist.DelY)
	
	
	
	red_mirr_dx_x=[]
	red_mirr_dx_y=[]
	red_mirr_dy_x=[]
	red_mirr_dy_y=[]
	red_mirr_dz_x=[]
	red_mirr_dz_y=[]
	red_mirr_tx_x=[]
	red_mirr_tx_y=[]
	red_mirr_ty_x=[]
	red_mirr_ty_y=[]
	red_mirr_tz_x=[]
	red_mirr_tz_y=[]
	
	
	for i in range(len(rmdelx)):
		if (i%6 == 0):
			red_mirr_dx_x.append(rmdelx[i])
		if (i%6 == 0):
			red_mirr_dx_y.append(rmdely[i])
		if (i%6 == 1):
			red_mirr_dy_x.append(rmdelx[i])
		if (i%6 == 1):
			red_mirr_dy_y.append(rmdely[i])
		if (i%6 == 2):
			red_mirr_dz_x.append(rmdelx[i])
		if (i%6 == 2):
			red_mirr_dz_y.append(rmdely[i])
		if (i%6 == 3):
			red_mirr_tx_x.append(rmdelx[i])
		if (i%6 == 3):
			red_mirr_tx_y.append(rmdely[i])
		if (i%6 == 4):
			red_mirr_ty_x.append(rmdelx[i])
		if (i%6 == 4):
			red_mirr_ty_y.append(rmdely[i])
		if (i%6 == 5):
			red_mirr_tz_x.append(rmdelx[i])
		if (i%6 == 5):
			red_mirr_tz_y.append(rmdely[i])
	
	
	red_mirr_dx_x=numpy.array(red_mirr_dx_x)*1e-6
	red_mirr_dx_y=numpy.array(red_mirr_dx_y)*1e-6
	red_mirr_dy_x=numpy.array(red_mirr_dy_x)*1e-6
	red_mirr_dy_y=numpy.array(red_mirr_dy_y)*1e-6
	red_mirr_dz_x=numpy.array(red_mirr_dz_x)*1e-6
	red_mirr_dz_y=numpy.array(red_mirr_dz_y)*1e-6
	red_mirr_tx_x=numpy.array(red_mirr_tx_x)*1e-6
	red_mirr_tx_y=numpy.array(red_mirr_tx_y)*1e-6
	red_mirr_ty_x=numpy.array(red_mirr_ty_x)*1e-6
	red_mirr_ty_y=numpy.array(red_mirr_ty_y)*1e-6
	red_mirr_tz_x=numpy.array(red_mirr_tz_x)*1e-6
	red_mirr_tz_y=numpy.array(red_mirr_tz_y)*1e-6
	
	
	
	red_datalist= pd.read_csv(red_gra_file,sep='\s+')
	rgele=red_datalist.Element
	rgcomp=red_datalist.Peturbation
	rgdelx=np.array(red_datalist.DelX)
	rgdely=np.array(red_datalist.DelY)
	
	red_gra_dx_x=[]
	red_gra_dx_y=[]
	red_gra_dy_x=[]
	red_gra_dy_y=[]
	red_gra_dz_x=[]
	red_gra_dz_y=[]
	red_gra_tx_x=[]
	red_gra_tx_y=[]
	red_gra_ty_x=[]
	red_gra_ty_y=[]
	red_gra_tz_x=[]
	red_gra_tz_y=[]
	
	
	for i in range(len(rgdelx)):
		if (i%6 == 0):
			red_gra_dx_x.append(rgdelx[i])
		if (i%6 == 0):
			red_gra_dx_y.append(rgdely[i])
		if (i%6 == 1):
			red_gra_dy_x.append(rgdelx[i])
		if (i%6 == 1):
			red_gra_dy_y.append(rgdely[i])
		if (i%6 == 2):
			red_gra_dz_x.append(rgdelx[i])
		if (i%6 == 2):
			red_gra_dz_y.append(rgdely[i])
		if (i%6 == 3):
			red_gra_tx_x.append(rgdelx[i])
		if (i%6 == 3):
			red_gra_tx_y.append(rgdely[i])
		if (i%6 == 4):
			red_gra_ty_x.append(rgdelx[i])
		if (i%6 == 4):
			red_gra_ty_y.append(rgdely[i])
		if (i%6 == 5):
			red_gra_tz_x.append(rgdelx[i])
		if (i%6 == 5):
			red_gra_tz_y.append(rgdely[i])
	
	red_gra_dx_x=numpy.array(red_gra_dx_x)*1e-6
	red_gra_dx_y=numpy.array(red_gra_dx_y)*1e-6
	red_gra_dy_x=numpy.array(red_gra_dy_x)*1e-6
	red_gra_dy_y=numpy.array(red_gra_dy_y)*1e-6
	red_gra_dz_x=numpy.array(red_gra_dz_x)*1e-6
	red_gra_dz_y=numpy.array(red_gra_dz_y)*1e-6
	red_gra_tx_x=numpy.array(red_gra_tx_x)*1e-6
	red_gra_tx_y=numpy.array(red_gra_tx_y)*1e-6
	red_gra_ty_x=numpy.array(red_gra_ty_x)*1e-6
	red_gra_ty_y=numpy.array(red_gra_ty_y)*1e-6
	red_gra_tz_x=numpy.array(red_gra_tz_x)*1e-6
	red_gra_tz_y=numpy.array(red_gra_tz_y)*1e-6
	
	
	
	
	
	#input bluedata
	
	blue_datalist= pd.read_csv(blue_coll_file,sep='\s+')
	blue_ofield=[]
	blue_infield=[]
	bx=np.array(blue_datalist.X)
	by=np.array(blue_datalist.Y)
	bele=blue_datalist.Element
	bcomp=blue_datalist.Peturbation
	bdelx=np.array(blue_datalist.DelX)
	bdely=np.array(blue_datalist.DelY)
	
	
	blue_x=[]
	blue_y=[]
	blue_coll_dx_x=[]
	blue_coll_dx_y=[]
	blue_coll_dy_x=[]
	blue_coll_dy_y=[]
	blue_coll_dz_x=[]
	blue_coll_dz_y=[]
	blue_coll_tx_x=[]
	blue_coll_tx_y=[]
	blue_coll_ty_x=[]
	blue_coll_ty_y=[]
	blue_coll_tz_x=[]
	blue_coll_tz_y=[]
	
	for i in range(len(bx)):
		if (i%6 == 0):
			blue_x.append(bx[i])
		if (i%6 == 0):
			blue_y.append(by[i])
		if (i%6 == 0):
			blue_coll_dx_x.append(bdelx[i])
		if (i%6 == 0):
			blue_coll_dx_y.append(bdely[i])
		if (i%6 == 1):
			blue_coll_dy_x.append(bdelx[i])
		if (i%6 == 1):
			blue_coll_dy_y.append(bdely[i])
		if (i%6 == 2):
			blue_coll_dz_x.append(bdelx[i])
		if (i%6 == 2):
			blue_coll_dz_y.append(bdely[i])
		if (i%6 == 3):
			blue_coll_tx_x.append(bdelx[i])
		if (i%6 == 3):
			blue_coll_tx_y.append(bdely[i])
		if (i%6 == 4):
			blue_coll_ty_x.append(bdelx[i])
		if (i%6 == 4):
			blue_coll_ty_y.append(bdely[i])
		if (i%6 == 5):
			blue_coll_tz_x.append(bdelx[i])
		if (i%6 == 5):
			blue_coll_tz_y.append(bdely[i])
	
	blue_x=numpy.array(blue_x)
	blue_y=numpy.array(blue_y)
	blue_coll_dx_x=numpy.array(blue_coll_dx_x)*1e-6
	blue_coll_dx_y=numpy.array(blue_coll_dx_y)*1e-6
	blue_coll_dy_x=numpy.array(blue_coll_dy_x)*1e-6
	blue_coll_dy_y=numpy.array(blue_coll_dy_y)*1e-6
	blue_coll_dz_x=numpy.array(blue_coll_dz_x)*1e-6
	blue_coll_dz_y=numpy.array(blue_coll_dz_y)*1e-6
	blue_coll_tx_x=numpy.array(blue_coll_tx_x)*1e-6
	blue_coll_tx_y=numpy.array(blue_coll_tx_y)*1e-6
	blue_coll_ty_x=numpy.array(blue_coll_ty_x)*1e-6
	blue_coll_ty_y=numpy.array(blue_coll_ty_y)*1e-6
	blue_coll_tz_x=numpy.array(blue_coll_tz_x)*1e-6
	blue_coll_tz_y=numpy.array(blue_coll_tz_y)*1e-6
	
	
	
	
	
	blue_datalist= pd.read_csv(blue_dchr_file,sep='\s+')
	bdele=blue_datalist.Element
	bdcomp=blue_datalist.Peturbation
	bddelx=np.array(blue_datalist.DelX)
	bddely=np.array(blue_datalist.DelY)
	
	blue_dchr_dx_x=[]
	blue_dchr_dx_y=[]
	blue_dchr_dy_x=[]
	blue_dchr_dy_y=[]
	blue_dchr_dz_x=[]
	blue_dchr_dz_y=[]
	blue_dchr_tx_x=[]
	blue_dchr_tx_y=[]
	blue_dchr_ty_x=[]
	blue_dchr_ty_y=[]
	blue_dchr_tz_x=[]
	blue_dchr_tz_y=[]
	
	
	for i in range(len(bddelx)):
		if (i%6 == 0):
			blue_dchr_dx_x.append(bddelx[i])
		if (i%6 == 0):
			blue_dchr_dx_y.append(bddely[i])
		if (i%6 == 1):
			blue_dchr_dy_x.append(bddelx[i])
		if (i%6 == 1):
			blue_dchr_dy_y.append(bddely[i])
		if (i%6 == 2):
			blue_dchr_dz_x.append(bddelx[i])
		if (i%6 == 2):
			blue_dchr_dz_y.append(bddely[i])
		if (i%6 == 3):
			blue_dchr_tx_x.append(bddelx[i])
		if (i%6 == 3):
			blue_dchr_tx_y.append(bddely[i])
		if (i%6 == 4):
			blue_dchr_ty_x.append(bddelx[i])
		if (i%6 == 4):
			blue_dchr_ty_y.append(bddely[i])
		if (i%6 == 5):
			blue_dchr_tz_x.append(bddelx[i])
		if (i%6 == 5):
			blue_dchr_tz_y.append(bddely[i])
	
	blue_dchr_dx_x=numpy.array(blue_dchr_dx_x)*1e-6
	blue_dchr_dx_y=numpy.array(blue_dchr_dx_y)*1e-6
	blue_dchr_dy_x=numpy.array(blue_dchr_dy_x)*1e-6
	blue_dchr_dy_y=numpy.array(blue_dchr_dy_y)*1e-6
	blue_dchr_dz_x=numpy.array(blue_dchr_dz_x)*1e-6
	blue_dchr_dz_y=numpy.array(blue_dchr_dz_y)*1e-6
	blue_dchr_tx_x=numpy.array(blue_dchr_tx_x)*1e-6
	blue_dchr_tx_y=numpy.array(blue_dchr_tx_y)*1e-6
	blue_dchr_ty_x=numpy.array(blue_dchr_ty_x)*1e-6
	blue_dchr_ty_y=numpy.array(blue_dchr_ty_y)*1e-6
	blue_dchr_tz_x=numpy.array(blue_dchr_tz_x)*1e-6
	blue_dchr_tz_y=numpy.array(blue_dchr_tz_y)*1e-6
	
	
	
	
	blue_datalist= pd.read_csv(blue_gra_file,sep='\s+')
	bgele=blue_datalist.Element
	bgcomp=blue_datalist.Peturbation
	bgdelx=np.array(blue_datalist.DelX)
	bgdely=np.array(blue_datalist.DelY)
	
	blue_gra_dx_x=[]
	blue_gra_dx_y=[]
	blue_gra_dy_x=[]
	blue_gra_dy_y=[]
	blue_gra_dz_x=[]
	blue_gra_dz_y=[]
	blue_gra_tx_x=[]
	blue_gra_tx_y=[]
	blue_gra_ty_x=[]
	blue_gra_ty_y=[]
	blue_gra_tz_x=[]
	blue_gra_tz_y=[]
	
	
	for i in range(len(bgdelx)):
		if (i%6 == 0):
			blue_gra_dx_x.append(bgdelx[i])
		if (i%6 == 0):
			blue_gra_dx_y.append(bgdely[i])
		if (i%6 == 1):
			blue_gra_dy_x.append(bgdelx[i])
		if (i%6 == 1):
			blue_gra_dy_y.append(bgdely[i])
		if (i%6 == 2):
			blue_gra_dz_x.append(bgdelx[i])
		if (i%6 == 2):
			blue_gra_dz_y.append(bgdely[i])
		if (i%6 == 3):
			blue_gra_tx_x.append(bgdelx[i])
		if (i%6 == 3):
			blue_gra_tx_y.append(bgdely[i])
		if (i%6 == 4):
			blue_gra_ty_x.append(bgdelx[i])
		if (i%6 == 4):
			blue_gra_ty_y.append(bgdely[i])
		if (i%6 == 5):
			blue_gra_tz_x.append(bgdelx[i])
		if (i%6 == 5):
			blue_gra_tz_y.append(bgdely[i])
	
	blue_gra_dx_x=numpy.array(blue_gra_dx_x)*1e-6
	blue_gra_dx_y=numpy.array(blue_gra_dx_y)*1e-6
	blue_gra_dy_x=numpy.array(blue_gra_dy_x)*1e-6
	blue_gra_dy_y=numpy.array(blue_gra_dy_y)*1e-6
	blue_gra_dz_x=numpy.array(blue_gra_dz_x)*1e-6
	blue_gra_dz_y=numpy.array(blue_gra_dz_y)*1e-6
	blue_gra_tx_x=numpy.array(blue_gra_tx_x)*1e-6
	blue_gra_tx_y=numpy.array(blue_gra_tx_y)*1e-6
	blue_gra_ty_x=numpy.array(blue_gra_ty_x)*1e-6
	blue_gra_ty_y=numpy.array(blue_gra_ty_y)*1e-6
	blue_gra_tz_x=numpy.array(blue_gra_tz_x)*1e-6
	blue_gra_tz_y=numpy.array(blue_gra_tz_y)*1e-6
	
	
	#Input: alldata finished
	
	
	#SelectField
	
	
	for i in range(len(red_x)):
		if (abs(red_x[i]) < 0):   # (abs(red_x[i]) < 70):
			red_infield.append(i)
		else:
			red_ofield.append(i)
	red_allfield=range(len(red_x))
	
	redo_x=red_x[red_ofield]
	redo_y=red_y[red_ofield]
	redo_coll_dx_x=red_coll_dx_x[red_ofield]
	redo_coll_dx_y= red_coll_dx_y[red_ofield]
	redo_coll_dy_x= red_coll_dy_x[red_ofield]
	redo_coll_dy_y= red_coll_dy_y[red_ofield]
	redo_coll_dz_x= red_coll_dz_x[red_ofield]
	redo_coll_dz_y= red_coll_dz_y[red_ofield]
	redo_coll_tx_x= red_coll_tx_x[red_ofield]
	redo_coll_tx_y= red_coll_tx_y[red_ofield]
	redo_coll_ty_x= red_coll_ty_x[red_ofield]
	redo_coll_ty_y= red_coll_ty_y[red_ofield]
	redo_coll_tz_x= red_coll_tz_x[red_ofield]
	redo_coll_tz_y= red_coll_tz_y[red_ofield]
	
	redo_mirr_dx_x=red_mirr_dx_x[red_ofield]
	redo_mirr_dx_y= red_mirr_dx_y[red_ofield]
	redo_mirr_dy_x= red_mirr_dy_x[red_ofield]
	redo_mirr_dy_y= red_mirr_dy_y[red_ofield]
	redo_mirr_dz_x= red_mirr_dz_x[red_ofield]
	redo_mirr_dz_y= red_mirr_dz_y[red_ofield]
	redo_mirr_tx_x= red_mirr_tx_x[red_ofield]
	redo_mirr_tx_y= red_mirr_tx_y[red_ofield]
	redo_mirr_ty_x= red_mirr_ty_x[red_ofield]
	redo_mirr_ty_y= red_mirr_ty_y[red_ofield]
	redo_mirr_tz_x= red_mirr_tz_x[red_ofield]
	redo_mirr_tz_y= red_mirr_tz_y[red_ofield]
	
	redo_gra_dx_x=red_gra_dx_x[red_ofield]
	redo_gra_dx_y= red_gra_dx_y[red_ofield]
	redo_gra_dy_x= red_gra_dy_x[red_ofield]
	redo_gra_dy_y= red_gra_dy_y[red_ofield]
	redo_gra_dz_x= red_gra_dz_x[red_ofield]
	redo_gra_dz_y= red_gra_dz_y[red_ofield]
	redo_gra_tx_x= red_gra_tx_x[red_ofield]
	redo_gra_tx_y= red_gra_tx_y[red_ofield]
	redo_gra_ty_x= red_gra_ty_x[red_ofield]
	redo_gra_ty_y= red_gra_ty_y[red_ofield]
	redo_gra_tz_x= red_gra_tz_x[red_ofield]
	redo_gra_tz_y= red_gra_tz_y[red_ofield]
	
	
	
	for i in range(len(blue_x)):
		if ((abs(blue_x[i]) < 0)): # ((blue_x[i]) < 75)):
			blue_infield.append(i)
		else:
			blue_ofield.append(i)
	blue_allfield=range(len(blue_x))
	
	blueo_x=blue_x[blue_ofield]
	blueo_y=blue_y[blue_ofield]
	blueo_coll_dx_x=blue_coll_dx_x[blue_ofield]
	blueo_coll_dx_y= blue_coll_dx_y[blue_ofield]
	blueo_coll_dy_x= blue_coll_dy_x[blue_ofield]
	blueo_coll_dy_y= blue_coll_dy_y[blue_ofield]
	blueo_coll_dz_x= blue_coll_dz_x[blue_ofield]
	blueo_coll_dz_y= blue_coll_dz_y[blue_ofield]
	blueo_coll_tx_x= blue_coll_tx_x[blue_ofield]
	blueo_coll_tx_y= blue_coll_tx_y[blue_ofield]
	blueo_coll_ty_x= blue_coll_ty_x[blue_ofield]
	blueo_coll_ty_y= blue_coll_ty_y[blue_ofield]
	blueo_coll_tz_x= blue_coll_tz_x[blue_ofield]
	blueo_coll_tz_y= blue_coll_tz_y[blue_ofield]
	
	blueo_dchr_dx_x=blue_dchr_dx_x[blue_ofield]
	blueo_dchr_dx_y= blue_dchr_dx_y[blue_ofield]
	blueo_dchr_dy_x= blue_dchr_dy_x[blue_ofield]
	blueo_dchr_dy_y= blue_dchr_dy_y[blue_ofield]
	blueo_dchr_dz_x= blue_dchr_dz_x[blue_ofield]
	blueo_dchr_dz_y= blue_dchr_dz_y[blue_ofield]
	blueo_dchr_tx_x= blue_dchr_tx_x[blue_ofield]
	blueo_dchr_tx_y= blue_dchr_tx_y[blue_ofield]
	blueo_dchr_ty_x= blue_dchr_ty_x[blue_ofield]
	blueo_dchr_ty_y= blue_dchr_ty_y[blue_ofield]
	blueo_dchr_tz_x= blue_dchr_tz_x[blue_ofield]
	blueo_dchr_tz_y= blue_dchr_tz_y[blue_ofield]
	
	blueo_gra_dx_x=blue_gra_dx_x[blue_ofield]
	blueo_gra_dx_y= blue_gra_dx_y[blue_ofield]
	blueo_gra_dy_x= blue_gra_dy_x[blue_ofield]
	blueo_gra_dy_y= blue_gra_dy_y[blue_ofield]
	blueo_gra_dz_x= blue_gra_dz_x[blue_ofield]
	blueo_gra_dz_y= blue_gra_dz_y[blue_ofield]
	blueo_gra_tx_x= blue_gra_tx_x[blue_ofield]
	blueo_gra_tx_y= blue_gra_tx_y[blue_ofield]
	blueo_gra_ty_x= blue_gra_ty_x[blue_ofield]
	blueo_gra_ty_y= blue_gra_ty_y[blue_ofield]
	blueo_gra_tz_x= blue_gra_tz_x[blue_ofield]
	blueo_gra_tz_y= blue_gra_tz_y[blue_ofield]
       
	blue_ofield=[]
	blue_infield=[]
	red_ofield=[]
	red_infield=[]
	
	for i in range(len(red_x)):
		if (abs(red_x[i]) < 70):   # (abs(red_x[i]) < 70):
			red_infield.append(i)
		else:
			red_ofield.append(i)
	red_allfield=range(len(red_x))
	
	redo1_x=red_x[red_ofield]
	redo1_y=red_y[red_ofield]
	redo1_coll_dx_x=red_coll_dx_x[red_ofield]
	redo1_coll_dx_y= red_coll_dx_y[red_ofield]
	redo1_coll_dy_x= red_coll_dy_x[red_ofield]
	redo1_coll_dy_y= red_coll_dy_y[red_ofield]
	redo1_coll_dz_x= red_coll_dz_x[red_ofield]
	redo1_coll_dz_y= red_coll_dz_y[red_ofield]
	redo1_coll_tx_x= red_coll_tx_x[red_ofield]
	redo1_coll_tx_y= red_coll_tx_y[red_ofield]
	redo1_coll_ty_x= red_coll_ty_x[red_ofield]
	redo1_coll_ty_y= red_coll_ty_y[red_ofield]
	redo1_coll_tz_x= red_coll_tz_x[red_ofield]
	redo1_coll_tz_y= red_coll_tz_y[red_ofield]
	
	redo1_mirr_dx_x=red_mirr_dx_x[red_ofield]
	redo1_mirr_dx_y= red_mirr_dx_y[red_ofield]
	redo1_mirr_dy_x= red_mirr_dy_x[red_ofield]
	redo1_mirr_dy_y= red_mirr_dy_y[red_ofield]
	redo1_mirr_dz_x= red_mirr_dz_x[red_ofield]
	redo1_mirr_dz_y= red_mirr_dz_y[red_ofield]
	redo1_mirr_tx_x= red_mirr_tx_x[red_ofield]
	redo1_mirr_tx_y= red_mirr_tx_y[red_ofield]
	redo1_mirr_ty_x= red_mirr_ty_x[red_ofield]
	redo1_mirr_ty_y= red_mirr_ty_y[red_ofield]
	redo1_mirr_tz_x= red_mirr_tz_x[red_ofield]
	redo1_mirr_tz_y= red_mirr_tz_y[red_ofield]
	
	redo1_gra_dx_x=red_gra_dx_x[red_ofield]
	redo1_gra_dx_y= red_gra_dx_y[red_ofield]
	redo1_gra_dy_x= red_gra_dy_x[red_ofield]
	redo1_gra_dy_y= red_gra_dy_y[red_ofield]
	redo1_gra_dz_x= red_gra_dz_x[red_ofield]
	redo1_gra_dz_y= red_gra_dz_y[red_ofield]
	redo1_gra_tx_x= red_gra_tx_x[red_ofield]
	redo1_gra_tx_y= red_gra_tx_y[red_ofield]
	redo1_gra_ty_x= red_gra_ty_x[red_ofield]
	redo1_gra_ty_y= red_gra_ty_y[red_ofield]
	redo1_gra_tz_x= red_gra_tz_x[red_ofield]
	redo1_gra_tz_y= red_gra_tz_y[red_ofield]
	
	
	
	for i in range(len(blue_x)):
		if ((abs(blue_x[i]) < 75)): # ((blue_x[i]) < 75)):
			blue_infield.append(i)
		else:
			blue_ofield.append(i)
	blue_allfield=range(len(blue_x))
	
	blueo1_x=blue_x[blue_ofield]
	blueo1_y=blue_y[blue_ofield]
	blueo1_coll_dx_x=blue_coll_dx_x[blue_ofield]
	blueo1_coll_dx_y= blue_coll_dx_y[blue_ofield]
	blueo1_coll_dy_x= blue_coll_dy_x[blue_ofield]
	blueo1_coll_dy_y= blue_coll_dy_y[blue_ofield]
	blueo1_coll_dz_x= blue_coll_dz_x[blue_ofield]
	blueo1_coll_dz_y= blue_coll_dz_y[blue_ofield]
	blueo1_coll_tx_x= blue_coll_tx_x[blue_ofield]
	blueo1_coll_tx_y= blue_coll_tx_y[blue_ofield]
	blueo1_coll_ty_x= blue_coll_ty_x[blue_ofield]
	blueo1_coll_ty_y= blue_coll_ty_y[blue_ofield]
	blueo1_coll_tz_x= blue_coll_tz_x[blue_ofield]
	blueo1_coll_tz_y= blue_coll_tz_y[blue_ofield]
	
	blueo1_dchr_dx_x=blue_dchr_dx_x[blue_ofield]
	blueo1_dchr_dx_y= blue_dchr_dx_y[blue_ofield]
	blueo1_dchr_dy_x= blue_dchr_dy_x[blue_ofield]
	blueo1_dchr_dy_y= blue_dchr_dy_y[blue_ofield]
	blueo1_dchr_dz_x= blue_dchr_dz_x[blue_ofield]
	blueo1_dchr_dz_y= blue_dchr_dz_y[blue_ofield]
	blueo1_dchr_tx_x= blue_dchr_tx_x[blue_ofield]
	blueo1_dchr_tx_y= blue_dchr_tx_y[blue_ofield]
	blueo1_dchr_ty_x= blue_dchr_ty_x[blue_ofield]
	blueo1_dchr_ty_y= blue_dchr_ty_y[blue_ofield]
	blueo1_dchr_tz_x= blue_dchr_tz_x[blue_ofield]
	blueo1_dchr_tz_y= blue_dchr_tz_y[blue_ofield]
	
	blueo1_gra_dx_x=blue_gra_dx_x[blue_ofield]
	blueo1_gra_dx_y= blue_gra_dx_y[blue_ofield]
	blueo1_gra_dy_x= blue_gra_dy_x[blue_ofield]
	blueo1_gra_dy_y= blue_gra_dy_y[blue_ofield]
	blueo1_gra_dz_x= blue_gra_dz_x[blue_ofield]
	blueo1_gra_dz_y= blue_gra_dz_y[blue_ofield]
	blueo1_gra_tx_x= blue_gra_tx_x[blue_ofield]
	blueo1_gra_tx_y= blue_gra_tx_y[blue_ofield]
	blueo1_gra_ty_x= blue_gra_ty_x[blue_ofield]
	blueo1_gra_ty_y= blue_gra_ty_y[blue_ofield]
	blueo1_gra_tz_x= blue_gra_tz_x[blue_ofield]
	blueo1_gra_tz_y= blue_gra_tz_y[blue_ofield]

	
	
	
	
    def onOpButton(self, event):
	global active,coll_comp,mirr_comp,dchr_comp,rgra_comp,bgra_comp
	global redo_x,redo_y,redo_coll_dx_x,redo_coll_dy_x,redo_coll_dz_x,redo_coll_dx_y,redo_coll_dy_y,redo_coll_dz_y,redo_coll_tx_x,redo_coll_ty_x,redo_coll_tz_x,redo_coll_tx_y,redo_coll_ty_y,redo_coll_tz_y,redo_x1,redo_y1,redo_x2,redo_y2
	global blueo_x,blueo_y,blueo_coll_dx_x,blueo_coll_dy_x,blueo_coll_dz_x,blueo_coll_dx_y,blueo_coll_dy_y,blueo_coll_dz_y,blueo_coll_tx_x,blueo_coll_ty_x,blueo_coll_tz_x,blueo_coll_tx_y,blueo_coll_ty_y,blueo_coll_tz_y,blueo_x1,blueo_y1,blueo_x2,blueo_y2
        global redo_x,redo_y,redo_mirr_dx_x,redo_mirr_dy_x,redo_mirr_dz_x,redo_mirr_dx_y,redo_mirr_dy_y,redo_mirr_dz_y,redo_mirr_tx_x,redo_mirr_ty_x,redo_mirr_tz_x,redo_mirr_tx_y,redo_mirr_ty_y,redo_mirr_tz_y
	global blueo_x,blueo_y,blueo_dchr_dx_x,blueo_dchr_dy_x,blueo_dchr_dz_x,blueo_dchr_dx_y,blueo_dchr_dy_y,blueo_dchr_dz_y,blueo_dchr_tx_x,blueo_dchr_ty_x,blueo_dchr_tz_x,blueo_dchr_tx_y,blueo_dchr_ty_y,blueo_dchr_tz_y
        global redo_x,redo_y,redo_gra_dx_x,redo_gra_dy_x,redo_gra_dz_x,redo_gra_dx_y,redo_gra_dy_y,redo_gra_dz_y,redo_gra_tx_x,redo_gra_ty_x,redo_gra_tz_x,redo_gra_tx_y,redo_gra_ty_y,redo_gra_tz_y
	global blueo_x,blueo_y,blueo_gra_dx_x,blueo_gra_dy_x,blueo_gra_dz_x,blueo_gra_dx_y,blueo_gra_dy_y,blueo_gra_dz_y,blueo_gra_tx_x,blueo_gra_ty_x,blueo_gra_tz_x,blueo_gra_tx_y,blueo_gra_ty_y,blueo_gra_tz_y
	
	global redo1_x,redo1_y,redo1_coll_dx_x,redo1_coll_dy_x,redo1_coll_dz_x,redo1_coll_dx_y,redo1_coll_dy_y,redo1_coll_dz_y,redo1_coll_tx_x,redo1_coll_ty_x,redo1_coll_tz_x,redo1_coll_tx_y,redo1_coll_ty_y,redo1_coll_tz_y,redo1_x1,redo1_y1,redo1_x2,redo1_y2
	global blueo1_x,blueo1_y,blueo1_coll_dx_x,blueo1_coll_dy_x,blueo1_coll_dz_x,blueo1_coll_dx_y,blueo1_coll_dy_y,blueo1_coll_dz_y,blueo1_coll_tx_x,blueo1_coll_ty_x,blueo1_coll_tz_x,blueo1_coll_tx_y,blueo1_coll_ty_y,blueo1_coll_tz_y,blueo1_x1,blueo1_y1,blueo1_x2,blueo1_y2
        global redo1_x,redo1_y,redo1_mirr_dx_x,redo1_mirr_dy_x,redo1_mirr_dz_x,redo1_mirr_dx_y,redo1_mirr_dy_y,redo1_mirr_dz_y,redo1_mirr_tx_x,redo1_mirr_ty_x,redo1_mirr_tz_x,redo1_mirr_tx_y,redo1_mirr_ty_y,redo1_mirr_tz_y
	global blueo1_x,blueo1_y,blueo1_dchr_dx_x,blueo1_dchr_dy_x,blueo1_dchr_dz_x,blueo1_dchr_dx_y,blueo1_dchr_dy_y,blueo1_dchr_dz_y,blueo1_dchr_tx_x,blueo1_dchr_ty_x,blueo1_dchr_tz_x,blueo1_dchr_tx_y,blueo1_dchr_ty_y,blueo1_dchr_tz_y
        global redo1_x,redo1_y,redo1_gra_dx_x,redo1_gra_dy_x,redo1_gra_dz_x,redo1_gra_dx_y,redo1_gra_dy_y,redo1_gra_dz_y,redo1_gra_tx_x,redo1_gra_ty_x,redo1_gra_tz_x,redo1_gra_tx_y,redo1_gra_ty_y,redo1_gra_tz_y
	global blueo1_x,blueo1_y,blueo1_gra_dx_x,blueo1_gra_dy_x,blueo1_gra_dz_x,blueo1_gra_dx_y,blueo1_gra_dy_y,blueo1_gra_dz_y,blueo1_gra_tx_x,blueo1_gra_ty_x,blueo1_gra_tz_x,blueo1_gra_tx_y,blueo1_gra_ty_y,blueo1_gra_tz_y
		
	global red_x,red_y,red_coll_dx_x,red_coll_dy_x,red_coll_dz_x,red_coll_dx_y,red_coll_dy_y,red_coll_dz_y,red_coll_tx_x,red_coll_ty_x,red_coll_tz_x,red_coll_tx_y,red_coll_ty_y,red_coll_tz_y,red_x1,red_y1,red_x2,red_y2
	global blue_x,blue_y,blue_coll_dx_x,blue_coll_dy_x,blue_coll_dz_x,blue_coll_dx_y,blue_coll_dy_y,blue_coll_dz_y,blue_coll_tx_x,blue_coll_ty_x,blue_coll_tz_x,blue_coll_tx_y,blue_coll_ty_y,blue_coll_tz_y,blue_x1,blue_y1,blue_x2,blue_y2
        global red_x,red_y,red_mirr_dx_x,red_mirr_dy_x,red_mirr_dz_x,red_mirr_dx_y,red_mirr_dy_y,red_mirr_dz_y,red_mirr_tx_x,red_mirr_ty_x,red_mirr_tz_x,red_mirr_tx_y,red_mirr_ty_y,red_mirr_tz_y
	global blue_x,blue_y,blue_dchr_dx_x,blue_dchr_dy_x,blue_dchr_dz_x,blue_dchr_dx_y,blue_dchr_dy_y,blue_dchr_dz_y,blue_dchr_tx_x,blue_dchr_ty_x,blue_dchr_tz_x,blue_dchr_tx_y,blue_dchr_ty_y,blue_dchr_tz_y
        global red_x,red_y,red_gra_dx_x,red_gra_dy_x,red_gra_dz_x,red_gra_dx_y,red_gra_dy_y,red_gra_dz_y,red_gra_tx_x,red_gra_ty_x,red_gra_tz_x,red_gra_tx_y,red_gra_ty_y,red_gra_tz_y
	global blue_x,blue_y,blue_gra_dx_x,blue_gra_dy_x,blue_gra_dz_x,blue_gra_dx_y,blue_gra_dy_y,blue_gra_dz_y,blue_gra_tx_x,blue_gra_ty_x,blue_gra_tz_x,blue_gra_tx_y,blue_gra_ty_y,blue_gra_tz_y
	
	fullfield=self.tab1.cb.GetValue()
	if fullfield==True:
	    redo_x= copy.copy(redo1_x)
	    redo_y= copy.copy(redo1_y)
	    redo_coll_dx_x= copy.copy(redo1_coll_dx_x)
	    redo_coll_dx_y= copy.copy(redo1_coll_dx_y)
	    redo_coll_dy_x= copy.copy(redo1_coll_dy_x)
	    redo_coll_dy_y= copy.copy(redo1_coll_dy_y)
	    redo_coll_dz_x= copy.copy(redo1_coll_dz_x)
	    redo_coll_dz_y= copy.copy(redo1_coll_dz_y)
	    redo_coll_tx_x= copy.copy(redo1_coll_tx_x)
	    redo_coll_tx_y= copy.copy(redo1_coll_tx_y)
	    redo_coll_ty_x= copy.copy(redo1_coll_ty_x)
	    redo_coll_ty_y= copy.copy(redo1_coll_ty_y)
	    redo_coll_tz_x= copy.copy(redo1_coll_tz_x)
	    redo_coll_tz_y= copy.copy(redo1_coll_tz_y)
	    
	    redo_mirr_dx_x= copy.copy(redo1_mirr_dx_x)
	    redo_mirr_dx_y= copy.copy(redo1_mirr_dx_y)
	    redo_mirr_dy_x= copy.copy(redo1_mirr_dy_x)
	    redo_mirr_dy_y= copy.copy(redo1_mirr_dy_y)
	    redo_mirr_dz_x= copy.copy(redo1_mirr_dz_x)
	    redo_mirr_dz_y= copy.copy(redo1_mirr_dz_y)
	    redo_mirr_tx_x= copy.copy(redo1_mirr_tx_x)
	    redo_mirr_tx_y= copy.copy(redo1_mirr_tx_y)
	    redo_mirr_ty_x= copy.copy(redo1_mirr_ty_x)
	    redo_mirr_ty_y= copy.copy(redo1_mirr_ty_y)
	    redo_mirr_tz_x= copy.copy(redo1_mirr_tz_x)
	    redo_mirr_tz_y= copy.copy(redo1_mirr_tz_y)
	    
	    redo_gra_dx_x= copy.copy(redo1_gra_dx_x)
	    redo_gra_dx_y= copy.copy(redo1_gra_dx_y)
	    redo_gra_dy_x= copy.copy(redo1_gra_dy_x)
	    redo_gra_dy_y= copy.copy(redo1_gra_dy_y)
	    redo_gra_dz_x= copy.copy(redo1_gra_dz_x)
	    redo_gra_dz_y= copy.copy(redo1_gra_dz_y)
	    redo_gra_tx_x= copy.copy(redo1_gra_tx_x)
	    redo_gra_tx_y= copy.copy(redo1_gra_tx_y)
	    redo_gra_ty_x= copy.copy(redo1_gra_ty_x)
	    redo_gra_ty_y= copy.copy(redo1_gra_ty_y)
	    redo_gra_tz_x= copy.copy(redo1_gra_tz_x)
	    redo_gra_tz_y= copy.copy(redo1_gra_tz_y)
	    
	    blueo_x= copy.copy(blueo1_x)
	    blueo_y= copy.copy(blueo1_y)
	    blueo_coll_dx_x= copy.copy(blueo1_coll_dx_x)
	    blueo_coll_dx_y= copy.copy(blueo1_coll_dx_y)
	    blueo_coll_dy_x= copy.copy(blueo1_coll_dy_x)
	    blueo_coll_dy_y= copy.copy(blueo1_coll_dy_y)
	    blueo_coll_dz_x= copy.copy(blueo1_coll_dz_x)
	    blueo_coll_dz_y= copy.copy(blueo1_coll_dz_y)
	    blueo_coll_tx_x= copy.copy(blueo1_coll_tx_x)
	    blueo_coll_tx_y= copy.copy(blueo1_coll_tx_y)
	    blueo_coll_ty_x= copy.copy(blueo1_coll_ty_x)
	    blueo_coll_ty_y= copy.copy(blueo1_coll_ty_y)
	    blueo_coll_tz_x= copy.copy(blueo1_coll_tz_x)
	    blueo_coll_tz_y= copy.copy(blueo1_coll_tz_y)
	    
	    blueo_dchr_dx_x= copy.copy(blueo1_dchr_dx_x)
	    blueo_dchr_dx_y= copy.copy(blueo1_dchr_dx_y)
	    blueo_dchr_dy_x= copy.copy(blueo1_dchr_dy_x)
	    blueo_dchr_dy_y= copy.copy(blueo1_dchr_dy_y)
	    blueo_dchr_dz_x= copy.copy(blueo1_dchr_dz_x)
	    blueo_dchr_dz_y= copy.copy(blueo1_dchr_dz_y)
	    blueo_dchr_tx_x= copy.copy(blueo1_dchr_tx_x)
	    blueo_dchr_tx_y= copy.copy(blueo1_dchr_tx_y)
	    blueo_dchr_ty_x= copy.copy(blueo1_dchr_ty_x)
	    blueo_dchr_ty_y= copy.copy(blueo1_dchr_ty_y)
	    blueo_dchr_tz_x= copy.copy(blueo1_dchr_tz_x)
	    blueo_dchr_tz_y= copy.copy(blueo1_dchr_tz_y)
	    
	    blueo_gra_dx_x= copy.copy(blueo1_gra_dx_x)
	    blueo_gra_dx_y= copy.copy(blueo1_gra_dx_y)
	    blueo_gra_dy_x= copy.copy(blueo1_gra_dy_x)
	    blueo_gra_dy_y= copy.copy(blueo1_gra_dy_y)
	    blueo_gra_dz_x= copy.copy(blueo1_gra_dz_x)
	    blueo_gra_dz_y= copy.copy(blueo1_gra_dz_y)
	    blueo_gra_tx_x= copy.copy(blueo1_gra_tx_x)
	    blueo_gra_tx_y= copy.copy(blueo1_gra_tx_y)
	    blueo_gra_ty_x= copy.copy(blueo1_gra_ty_x)
	    blueo_gra_ty_y= copy.copy(blueo1_gra_ty_y)
	    blueo_gra_tz_x= copy.copy(blueo1_gra_tz_x)
	    blueo_gra_tz_y= copy.copy(blueo1_gra_tz_y)
	
	if fullfield==False:
	    redo_x= copy.copy(red_x)
	    redo_y= copy.copy(red_y)
	    redo_coll_dx_x= copy.copy(red_coll_dx_x)
	    redo_coll_dx_y= copy.copy(red_coll_dx_y)
	    redo_coll_dy_x= copy.copy(red_coll_dy_x)
	    redo_coll_dy_y= copy.copy(red_coll_dy_y)
	    redo_coll_dz_x= copy.copy(red_coll_dz_x)
	    redo_coll_dz_y= copy.copy(red_coll_dz_y)
	    redo_coll_tx_x= copy.copy(red_coll_tx_x)
	    redo_coll_tx_y= copy.copy(red_coll_tx_y)
	    redo_coll_ty_x= copy.copy(red_coll_ty_x)
	    redo_coll_ty_y= copy.copy(red_coll_ty_y)
	    redo_coll_tz_x= copy.copy(red_coll_tz_x)
	    redo_coll_tz_y= copy.copy(red_coll_tz_y)
	    
	    redo_mirr_dx_x= copy.copy(red_mirr_dx_x)
	    redo_mirr_dx_y= copy.copy(red_mirr_dx_y)
	    redo_mirr_dy_x= copy.copy(red_mirr_dy_x)
	    redo_mirr_dy_y= copy.copy(red_mirr_dy_y)
	    redo_mirr_dz_x= copy.copy(red_mirr_dz_x)
	    redo_mirr_dz_y= copy.copy(red_mirr_dz_y)
	    redo_mirr_tx_x= copy.copy(red_mirr_tx_x)
	    redo_mirr_tx_y= copy.copy(red_mirr_tx_y)
	    redo_mirr_ty_x= copy.copy(red_mirr_ty_x)
	    redo_mirr_ty_y= copy.copy(red_mirr_ty_y)
	    redo_mirr_tz_x= copy.copy(red_mirr_tz_x)
	    redo_mirr_tz_y= copy.copy(red_mirr_tz_y)
	    
	    redo_gra_dx_x= copy.copy(red_gra_dx_x)
	    redo_gra_dx_y= copy.copy(red_gra_dx_y)
	    redo_gra_dy_x= copy.copy(red_gra_dy_x)
	    redo_gra_dy_y= copy.copy(red_gra_dy_y)
	    redo_gra_dz_x= copy.copy(red_gra_dz_x)
	    redo_gra_dz_y= copy.copy(red_gra_dz_y)
	    redo_gra_tx_x= copy.copy(red_gra_tx_x)
	    redo_gra_tx_y= copy.copy(red_gra_tx_y)
	    redo_gra_ty_x= copy.copy(red_gra_ty_x)
	    redo_gra_ty_y= copy.copy(red_gra_ty_y)
	    redo_gra_tz_x= copy.copy(red_gra_tz_x)
	    redo_gra_tz_y= copy.copy(red_gra_tz_y)
	    
	    blueo_x= copy.copy(blue_x)
	    blueo_y= copy.copy(blue_y)
	    blueo_coll_dx_x= copy.copy(blue_coll_dx_x)
	    blueo_coll_dx_y= copy.copy(blue_coll_dx_y)
	    blueo_coll_dy_x= copy.copy(blue_coll_dy_x)
	    blueo_coll_dy_y= copy.copy(blue_coll_dy_y)
	    blueo_coll_dz_x= copy.copy(blue_coll_dz_x)
	    blueo_coll_dz_y= copy.copy(blue_coll_dz_y)
	    blueo_coll_tx_x= copy.copy(blue_coll_tx_x)
	    blueo_coll_tx_y= copy.copy(blue_coll_tx_y)
	    blueo_coll_ty_x= copy.copy(blue_coll_ty_x)
	    blueo_coll_ty_y= copy.copy(blue_coll_ty_y)
	    blueo_coll_tz_x= copy.copy(blue_coll_tz_x)
	    blueo_coll_tz_y= copy.copy(blue_coll_tz_y)
	    
	    blueo_dchr_dx_x= copy.copy(blue_dchr_dx_x)
	    blueo_dchr_dx_y= copy.copy(blue_dchr_dx_y)
	    blueo_dchr_dy_x= copy.copy(blue_dchr_dy_x)
	    blueo_dchr_dy_y= copy.copy(blue_dchr_dy_y)
	    blueo_dchr_dz_x= copy.copy(blue_dchr_dz_x)
	    blueo_dchr_dz_y= copy.copy(blue_dchr_dz_y)
	    blueo_dchr_tx_x= copy.copy(blue_dchr_tx_x)
	    blueo_dchr_tx_y= copy.copy(blue_dchr_tx_y)
	    blueo_dchr_ty_x= copy.copy(blue_dchr_ty_x)
	    blueo_dchr_ty_y= copy.copy(blue_dchr_ty_y)
	    blueo_dchr_tz_x= copy.copy(blue_dchr_tz_x)
	    blueo_dchr_tz_y= copy.copy(blue_dchr_tz_y)
	    
	    blueo_gra_dx_x= copy.copy(blue_gra_dx_x)
	    blueo_gra_dx_y= copy.copy(blue_gra_dx_y)
	    blueo_gra_dy_x= copy.copy(blue_gra_dy_x)
	    blueo_gra_dy_y= copy.copy(blue_gra_dy_y)
	    blueo_gra_dz_x= copy.copy(blue_gra_dz_x)
	    blueo_gra_dz_y= copy.copy(blue_gra_dz_y)
	    blueo_gra_tx_x= copy.copy(blue_gra_tx_x)
	    blueo_gra_tx_y= copy.copy(blue_gra_tx_y)
	    blueo_gra_ty_x= copy.copy(blue_gra_ty_x)
	    blueo_gra_ty_y= copy.copy(blue_gra_ty_y)
	    blueo_gra_tz_x= copy.copy(blue_gra_tz_x)
	    blueo_gra_tz_y= copy.copy(blue_gra_tz_y)
		
	
	active=self.tab1.list_param.GetChecked()
	coll_comp=self.tab1.list_coll.GetChecked()
	mirr_comp=self.tab1.list_mirr.GetChecked()
	dchr_comp=self.tab1.list_dchr.GetChecked()
	rgra_comp=self.tab1.list_rgra.GetChecked()
	bgra_comp=self.tab1.list_bgra.GetChecked()				
	print 'Starting Optimization'
	
        self.pp=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
	self.l=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
	self.u=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
	#print self.pp
	if 0 in active :
		if 0 in coll_comp:
			self.pp[0]=0.01
		if 1 in coll_comp:
			self.pp[1]=0.01			
		if 2 in coll_comp:
			self.pp[2]=0.01
		if 3 in coll_comp:
			self.pp[3]=0.01
		if 4 in coll_comp:
			self.pp[4]=0.01
		if 5 in coll_comp:
			self.pp[5]=0.01

	if 1 in active :
		if 0 in mirr_comp:
			self.pp[6]=0.01
		if 1 in mirr_comp:
			self.pp[7]=0.01			
		if 2 in mirr_comp:
			self.pp[8]=0.01
		if 3 in mirr_comp:
			self.pp[9]=0.01
		if 4 in mirr_comp:
			self.pp[10]=0.01
		if 5 in mirr_comp:
			self.pp[11]=0.01
	if 3 in active :
	    self.pp[12]=0.01
	    self.pp[13]=0.01
	    self.pp[14]=0.01
	    self.pp[15]=0.01
	    self.pp[16]=0.01
	    self.pp[17]=0.01
	    self.l[12]=-2
	    self.l[13]=-2
	    self.l[14]=-2
	    self.l[15]=-2
	    self.l[16]=-2
	    self.l[17]=-2
	    self.u[12]=2
	    self.u[13]=2
	    self.u[14]=2
	    self.u[15]=2
	    self.u[16]=2
	    self.u[17]=2
	if 2 in active :
	    self.pp[18]=0.01
	    self.pp[19]=0.01
	    self.pp[20]=0.01
	    self.pp[21]=0.01
	    self.pp[22]=0.01
	    self.pp[23]=0.01
	    self.l[18]=-2
	    self.l[19]=-2
	    self.l[20]=-2
	    self.l[21]=-2
	    self.l[22]=-2
	    self.l[23]=-2
	    self.u[18]=2
	    self.u[19]=2
	    self.u[20]=2
	    self.u[21]=2
	    self.u[22]=2
	    self.u[23]=2
	if 4 in active :
	    self.pp[24]=0.01
	    self.pp[25]=0.01
	    self.pp[26]=0.01
	    self.pp[27]=0.01
	    self.pp[28]=0.01
	    self.pp[29]=0.01
	    self.l[24]=-2
	    self.l[25]=-2
	    self.l[26]=-2
	    self.l[27]=-2
	    self.l[28]=-2
	    self.l[29]=-2
	    self.u[24]=2
	    self.u[25]=2
	    self.u[26]=2
	    self.u[27]=2
	    self.u[28]=2
	    self.u[29]=2
	if 5 in active :
	    self.pp[30]=.02
	    self.pp[31]=.02
	    self.u[30]=9
	    self.u[31]=9
	    self.l[30]=-9
	    self.l[31]=-9
	if 6 in active :
	    self.pp[32]=.0001
	    self.pp[33]=.0001
	    self.pp[34]=.0001
	    self.pp[35]=.0001
	    
	self.l=[-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000,-3000]
	self.u=[3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000,3000]
	bounds = [(low, high) for low, high in zip(self.l, self.u)]
	

	#print(len(self.pp))		
	minimizer_kwargs = dict(method="L-BFGS-B", bounds=bounds)	
	minimizer_kwargs = dict(method="L-BFGS-B")	
	#print self.l
	#print self.u
	#print self.pp
	self.mytakestep = MyTakeStep()
	it=int(self.tab1.selectIter.GetValue())
	#res = optimize.anneal(f, self.pp,schedule='boltzmann',full_output=True,lower=self.l,upper=self.u ,maxiter=50, dwell=100,learn_rate=.8, disp=True)
	res = optimize.basinhopping(f, self.pp,niter=it,stepsize=.5,take_step=self.mytakestep,callback=print_fun)

	#print res
	#res=res[0]
	print 'final',res
	res=res.x
	self.act_coll_dx=res[0] 
	self.act_coll_dy=res[1] 
	self.act_coll_dz=res[2] 
	self.act_coll_tx=res[3] 
	self.act_coll_ty=res[4] 
	self.act_coll_tz=res[5] 
	
	self.act_mirr_dx=res[6] 
	self.act_mirr_dy=res[7] 
	self.act_mirr_dz=res[8] 
	self.act_mirr_tx=res[9] 
	self.act_mirr_ty=res[10] 
	self.act_mirr_tz=res[11] 	

	self.act_rgra_dx=res[12] 
	self.act_rgra_dy=res[13] 
	self.act_rgra_dz=res[14] 
	self.act_rgra_tx=res[15] 
	self.act_rgra_ty=res[16] 
	self.act_rgra_tz=res[17] 
	
	self.act_dchr_dx=res[18] 
	self.act_dchr_dy=res[19] 
	self.act_dchr_dz=res[20] 
	self.act_dchr_tx=res[21] 
	self.act_dchr_ty=res[22] 
	self.act_dchr_tz=res[23] 
	
	self.act_bgra_dx=res[24] 
	self.act_bgra_dy=res[25] 
	self.act_bgra_dz=res[26] 
	self.act_bgra_tx=res[27] 
	self.act_bgra_ty=res[28] 
	self.act_bgra_tz=res[29] 
	
	self.det_red_rot=res[30]
	self.det_blue_rot=res[31]
 

	self.det_red_shiftx = res[32]
	self.det_red_shifty = res[33]
	self.det_blue_shiftx = res[34]
	self.det_blue_shifty = res[35]			
 

	red_x2= red_x1+self.act_coll_dx*red_coll_dx_x +self.act_coll_dy*red_coll_dy_x	+ self.act_coll_dz*red_coll_dz_x + self.act_coll_tx*red_coll_tx_x +self.act_coll_ty*red_coll_ty_x+ self.act_coll_tz*red_coll_tz_x + \
		 self.act_mirr_dx*red_mirr_dx_x +self.act_mirr_dy*red_mirr_dy_x	+ self.act_mirr_dz*red_mirr_dz_x + self.act_mirr_tx*red_mirr_tx_x +self.act_mirr_ty*red_mirr_ty_x+ self.act_mirr_tz*red_mirr_tz_x + \
		 self.act_rgra_dx*red_gra_dx_x +self.act_rgra_dy*red_gra_dy_x	+ self.act_rgra_dz*red_gra_dz_x + self.act_rgra_tx*red_gra_tx_x +self.act_rgra_ty*red_gra_ty_x+ self.act_rgra_tz*red_gra_tz_x
	red_y2= red_y1+self.act_coll_dx*red_coll_dx_y +self.act_coll_dy*red_coll_dy_y	+ self.act_coll_dz*red_coll_dz_y + self.act_coll_tx*red_coll_tx_y +self.act_coll_ty*red_coll_ty_y+ self.act_coll_tz*red_coll_tz_y + \
		 self.act_mirr_dx*red_mirr_dx_y +self.act_mirr_dy*red_mirr_dy_y	+ self.act_mirr_dz*red_mirr_dz_y + self.act_mirr_tx*red_mirr_tx_y +self.act_mirr_ty*red_mirr_ty_y+ self.act_mirr_tz*red_mirr_tz_y + \
		 self.act_rgra_dx*red_gra_dx_y +self.act_rgra_dy*red_gra_dy_y	+ self.act_rgra_dz*red_gra_dz_y + self.act_rgra_tx*red_gra_tx_y +self.act_rgra_ty*red_gra_ty_y+ self.act_rgra_tz*red_gra_tz_y

	blue_x2= blue_x1+self.act_coll_dx*blue_coll_dx_x +self.act_coll_dy*blue_coll_dy_x	+ self.act_coll_dz*blue_coll_dz_x + self.act_coll_tx*blue_coll_tx_x +self.act_coll_ty*blue_coll_ty_x+ self.act_coll_tz*blue_coll_tz_x + \
		 self.act_dchr_dx*blue_dchr_dx_x +self.act_dchr_dy*blue_dchr_dy_x	+ self.act_dchr_dz*blue_dchr_dz_x + self.act_dchr_tx*blue_dchr_tx_x +self.act_dchr_ty*blue_dchr_ty_x+ self.act_dchr_tz*blue_dchr_tz_x + \
		 self.act_bgra_dx*blue_gra_dx_x +self.act_bgra_dy*blue_gra_dy_x	+ self.act_bgra_dz*blue_gra_dz_x + self.act_bgra_tx*blue_gra_tx_x +self.act_bgra_ty*blue_gra_ty_x+ self.act_bgra_tz*blue_gra_tz_x
	blue_y2= blue_y1+self.act_coll_dx*blue_coll_dx_y +self.act_coll_dy*blue_coll_dy_y	+ self.act_coll_dz*blue_coll_dz_y + self.act_coll_tx*blue_coll_tx_y +self.act_coll_ty*blue_coll_ty_y+ self.act_coll_tz*blue_coll_tz_y + \
		  self.act_dchr_dx*blue_dchr_dx_y +self.act_dchr_dy*blue_dchr_dy_y+ self.act_dchr_dz*blue_dchr_dz_y + self.act_dchr_tx*blue_dchr_tx_y +self.act_dchr_ty*blue_dchr_ty_y+ self.act_dchr_tz*blue_dchr_tz_y + \
		  self.act_bgra_dx*blue_gra_dx_y +self.act_bgra_dy*blue_gra_dy_y+ self.act_bgra_dz*blue_gra_dz_y + self.act_bgra_tx*blue_gra_tx_y +self.act_bgra_ty*blue_gra_ty_y+ self.act_bgra_tz*blue_gra_tz_y
        
	rxrot2 = red_x2*np.cos((self.det_red_rot))-red_y2*np.sin((self.det_red_rot))
        ryrot2=  red_x2*np.sin((self.det_red_rot))+red_y2*np.cos((self.det_red_rot))

	bxrot2 = blue_x2*np.cos(self.det_blue_rot)-blue_y2*np.sin(self.det_blue_rot)
        byrot2=  blue_x2*np.sin(self.det_blue_rot)+blue_y2*np.cos(self.det_blue_rot)

	red_x2=rxrot2+self.det_red_shiftx
	red_y2=ryrot2+self.det_red_shifty
	blue_x2=bxrot2+self.det_blue_shiftx
	blue_y2=byrot2+self.det_blue_shifty

	
	self.tab3.myGrid.SetCellValue(0,0,str("%0.4f" %self.act_coll_dx))
	self.tab3.myGrid.SetCellValue(0,1,str("%0.4f" %self.act_coll_dy))
	self.tab3.myGrid.SetCellValue(0,2,str("%0.4f" %self.act_coll_dz))
	self.tab3.myGrid.SetCellValue(0,3,str("%0.4f" %self.act_coll_tx))
	self.tab3.myGrid.SetCellValue(0,4,str("%0.4f" %self.act_coll_ty))
	self.tab3.myGrid.SetCellValue(0,5,str("%0.4f" %self.act_coll_tz))

	self.tab3.myGrid.SetCellValue(1,0,str("%0.4f" %self.act_mirr_dx))
	self.tab3.myGrid.SetCellValue(1,1,str("%0.4f" %self.act_mirr_dy))
	self.tab3.myGrid.SetCellValue(1,2,str("%0.4f" %self.act_mirr_dz))
	self.tab3.myGrid.SetCellValue(1,3,str("%0.4f" %self.act_mirr_tx))
	self.tab3.myGrid.SetCellValue(1,4,str("%0.4f" %self.act_mirr_ty))
	self.tab3.myGrid.SetCellValue(1,5,str("%0.4f" %self.act_mirr_tz))
	
	self.tab3.myGrid.SetCellValue(2,0,str("%0.4f" %self.act_dchr_dx))
	self.tab3.myGrid.SetCellValue(2,1,str("%0.4f" %self.act_dchr_dy))
	self.tab3.myGrid.SetCellValue(2,2,str("%0.4f" %self.act_dchr_dz))
	self.tab3.myGrid.SetCellValue(2,3,str("%0.4f" %self.act_dchr_tx))
	self.tab3.myGrid.SetCellValue(2,4,str("%0.4f" %self.act_dchr_ty))
	self.tab3.myGrid.SetCellValue(2,5,str("%0.4f" %self.act_dchr_tz))
	
	self.tab3.myGrid.SetCellValue(3,0,str("%0.4f" %self.act_rgra_dx))
	self.tab3.myGrid.SetCellValue(3,1,str("%0.4f" %self.act_rgra_dy))
	self.tab3.myGrid.SetCellValue(3,2,str("%0.4f" %self.act_rgra_dz))
	self.tab3.myGrid.SetCellValue(3,3,str("%0.4f" %self.act_rgra_tx))
	self.tab3.myGrid.SetCellValue(3,4,str("%0.4f" %self.act_rgra_ty))
	self.tab3.myGrid.SetCellValue(3,5,str("%0.4f" %self.act_rgra_tz))
	
	self.tab3.myGrid.SetCellValue(4,0,str("%0.4f" %self.act_bgra_dx))
	self.tab3.myGrid.SetCellValue(4,1,str("%0.4f" %self.act_bgra_dy))
	self.tab3.myGrid.SetCellValue(4,2,str("%0.4f" %self.act_bgra_dz))
	self.tab3.myGrid.SetCellValue(4,3,str("%0.4f" %self.act_bgra_tx))
	self.tab3.myGrid.SetCellValue(4,4,str("%0.4f" %self.act_bgra_ty))
	self.tab3.myGrid.SetCellValue(4,5,str("%0.4f" %self.act_bgra_tz))
	
	self.tab3.myGrid.SetCellValue(8,0,str("%0.4f" %(np.degrees(self.det_red_rot))))
	self.tab3.myGrid.SetCellValue(8,1,str("%0.4f" %(np.degrees(self.det_blue_rot))))
	
	self.tab3.myGrid.SetCellValue(6,0,str("%0.4f" %(self.det_red_shiftx*1e3)))
	self.tab3.myGrid.SetCellValue(6,1,str("%0.4f" %(self.det_red_shifty*1e3)))
	self.tab3.myGrid.SetCellValue(6,2,str("%0.4f" %(self.det_blue_shiftx*1e3)))
	self.tab3.myGrid.SetCellValue(6,3,str("%0.4f" %(self.det_blue_shifty*1e3)))

	self.panell.resscale=int(self.tab1.selectResScale.GetValue())	
	self.panell2.resscale=int(self.tab1.selectResScale.GetValue())		
	self.panell.drawredres()
	self.panell2.drawblueres()
	print 'final'
	print 'red rms residual: %0.4f microns' %(1000*(sum(np.sqrt((red_x-red_x2)**2+(red_y-red_y2)**2))/len(red_y2)))
	print 'blue rms residual: %0.4f microns' %(1000*(sum(np.sqrt((blue_x-blue_x2)**2+(blue_y-blue_y2)**2))/len(blue_y2)))
	print 'red min residual : %0.4f microns' %(1000*(min((np.sqrt((red_x-red_x2)**2+(red_y-red_y2)**2)))))
	print 'red max residual : %0.4f microns' %(1000*max((np.sqrt((red_x-red_x2)**2+(red_y-red_y2)**2))))
        print 'blue min residual : %0.4f microns' %(1000*min((np.sqrt((blue_x-blue_x2)**2+(blue_y-blue_y2)**2))))
	print 'blue max residual : %0.4f microns' %(1000*max((np.sqrt((blue_x-blue_x2)**2+(blue_y-blue_y2)**2))))
        self.redrmsres=1000*(sum(np.sqrt((red_x-red_x2)**2+(red_y-red_y2)**2))/len(red_y2))
	self.bluermsres=1000*(sum(np.sqrt((blue_x-blue_x2)**2+(blue_y-blue_y2)**2))/len(blue_y2))
	print 'finished'
	#self.panell2.drawblue()
	
    def onFCSButton(self, event):
	self.panell.drawred()
	self.panell2.drawblue()
    def onCompButton(self, event):
	self.panell.drawredcomp()
	self.panell2.drawbluecomp()
    def onResButton(self, event):
	self.panell.resscale=int(self.tab1.selectResScale.GetValue())	
	self.panell2.resscale=int(self.tab1.selectResScale.GetValue())		
	self.panell.drawredres()
	self.panell2.drawblueres()
    def onFlexSave(self, event):
	flog = open(self.tab2.selectFile.GetValue(), 'wb')
	flog.write("Element dx dy dz tx ty tz ")
	flog.write("\nCollimator " + self.tab2.myGrid.GetCellValue(0,0)+" " +self.tab2.myGrid.GetCellValue(0,1)+" " +self.tab2.myGrid.GetCellValue(0,2)+" " +self.tab2.myGrid.GetCellValue(0,3)+" " +self.tab2.myGrid.GetCellValue(0,4)+" " +self.tab2.myGrid.GetCellValue(0,5)) 
	flog.write("\nMirror " + self.tab2.myGrid.GetCellValue(1,0)+" " +self.tab2.myGrid.GetCellValue(1,1)+" " +self.tab2.myGrid.GetCellValue(1,2)+" " +self.tab2.myGrid.GetCellValue(1,3)+" " +self.tab2.myGrid.GetCellValue(1,4)+" " +self.tab2.myGrid.GetCellValue(1,5)) 
	flog.write("\nDichroic " + self.tab2.myGrid.GetCellValue(2,0)+" " +self.tab2.myGrid.GetCellValue(2,1)+" " +self.tab2.myGrid.GetCellValue(2,2)+" " +self.tab2.myGrid.GetCellValue(2,3)+" " +self.tab2.myGrid.GetCellValue(2,4)+" " +self.tab2.myGrid.GetCellValue(2,5)) 
	flog.write("\nRed_Grating " + self.tab2.myGrid.GetCellValue(3,0)+" " +self.tab2.myGrid.GetCellValue(3,1)+" " +self.tab2.myGrid.GetCellValue(3,2)+" " +self.tab2.myGrid.GetCellValue(3,3)+" " +self.tab2.myGrid.GetCellValue(3,4)+" " +self.tab2.myGrid.GetCellValue(3,5)) 
	flog.write("\nBlue_Grating " + self.tab2.myGrid.GetCellValue(4,0)+" " +self.tab2.myGrid.GetCellValue(4,1)+" " +self.tab2.myGrid.GetCellValue(4,2)+" " +self.tab2.myGrid.GetCellValue(4,3)+" " +self.tab2.myGrid.GetCellValue(4,4)+" " +self.tab2.myGrid.GetCellValue(4,5)) 
	flog.write("\nRed_Detector_Shift " + self.tab2.myGrid.GetCellValue(6,0)+" "+ self.tab2.myGrid.GetCellValue(6,1)) 
	flog.write("\nBlue_Detector_Shift " + self.tab2.myGrid.GetCellValue(6,2)+" "+ self.tab2.myGrid.GetCellValue(6,3)) 
	flog.write("\nRed_Detector_Rotation_(deg) " + self.tab2.myGrid.GetCellValue(8,0))
	flog.write("\nBlue_Detector_Rotation_(deg) " + self.tab2.myGrid.GetCellValue(8,1))
	flog.write("\n")
	flog.close()	

    def onFlexLoad(self, event):
	flog = open(self.tab2.selectFile.GetValue(), 'r')
	val= pd.read_csv(self.tab2.selectFile.GetValue(),sep='\s+')

	self.tab2.myGrid.SetCellValue(0,0,str(val.values[0][1]))
	self.tab2.myGrid.SetCellValue(0,1,str(val.values[0][2]))
	self.tab2.myGrid.SetCellValue(0,2,str(val.values[0][3]))
	self.tab2.myGrid.SetCellValue(0,3,str(val.values[0][4]))
	self.tab2.myGrid.SetCellValue(0,4,str(val.values[0][5]))
	self.tab2.myGrid.SetCellValue(0,5,str(val.values[0][6]))
	
	self.tab2.myGrid.SetCellValue(1,0,str(val.values[1][1]))
	self.tab2.myGrid.SetCellValue(1,1,str(val.values[1][2]))
	self.tab2.myGrid.SetCellValue(1,2,str(val.values[1][3]))
	self.tab2.myGrid.SetCellValue(1,3,str(val.values[1][4]))
	self.tab2.myGrid.SetCellValue(1,4,str(val.values[1][5]))
	self.tab2.myGrid.SetCellValue(1,5,str(val.values[1][6]))
	
	self.tab2.myGrid.SetCellValue(2,0,str(val.values[2][1]))
	self.tab2.myGrid.SetCellValue(2,1,str(val.values[2][2]))
	self.tab2.myGrid.SetCellValue(2,2,str(val.values[2][3]))
	self.tab2.myGrid.SetCellValue(2,3,str(val.values[2][4]))
	self.tab2.myGrid.SetCellValue(2,4,str(val.values[2][5]))
	self.tab2.myGrid.SetCellValue(2,5,str(val.values[2][6]))
	
	self.tab2.myGrid.SetCellValue(3,0,str(val.values[3][1]))
	self.tab2.myGrid.SetCellValue(3,1,str(val.values[3][2]))
	self.tab2.myGrid.SetCellValue(3,2,str(val.values[3][3]))
	self.tab2.myGrid.SetCellValue(3,3,str(val.values[3][4]))
	self.tab2.myGrid.SetCellValue(3,4,str(val.values[3][5]))
	self.tab2.myGrid.SetCellValue(3,5,str(val.values[3][6]))
	
	self.tab2.myGrid.SetCellValue(4,0,str(val.values[4][1]))
	self.tab2.myGrid.SetCellValue(4,1,str(val.values[4][2]))
	self.tab2.myGrid.SetCellValue(4,2,str(val.values[4][3]))
	self.tab2.myGrid.SetCellValue(4,3,str(val.values[4][4]))
	self.tab2.myGrid.SetCellValue(4,4,str(val.values[4][5]))
	self.tab2.myGrid.SetCellValue(4,5,str(val.values[4][6]))

	flog.close()	
	
	
    def onCompSave(self, event):
	flog = open(self.tab3.selectFile.GetValue(), 'wb')
	flog.write("Compensation_Data   \n")
	flog.write("Element dx dy dz tx ty tz ")
	flog.write("\nCollimator " + self.tab3.myGrid.GetCellValue(0,0)+" " +self.tab3.myGrid.GetCellValue(0,1)+" " +self.tab3.myGrid.GetCellValue(0,2)+" " +self.tab3.myGrid.GetCellValue(0,3)+" " +self.tab3.myGrid.GetCellValue(0,4)+" " +self.tab3.myGrid.GetCellValue(0,5)) 
	flog.write("\nMirror " + self.tab3.myGrid.GetCellValue(1,0)+" " +self.tab3.myGrid.GetCellValue(1,1)+" " +self.tab3.myGrid.GetCellValue(1,2)+" " +self.tab3.myGrid.GetCellValue(1,3)+" " +self.tab3.myGrid.GetCellValue(1,4)+" " +self.tab3.myGrid.GetCellValue(1,5)) 
	flog.write("\nDichroic " + self.tab3.myGrid.GetCellValue(2,0)+" " +self.tab3.myGrid.GetCellValue(2,1)+" " +self.tab3.myGrid.GetCellValue(2,2)+" " +self.tab3.myGrid.GetCellValue(2,3)+" " +self.tab3.myGrid.GetCellValue(2,4)+" " +self.tab3.myGrid.GetCellValue(2,5)) 
	flog.write("\nRed_Grating " + self.tab3.myGrid.GetCellValue(3,0)+" " +self.tab3.myGrid.GetCellValue(3,1)+" " +self.tab3.myGrid.GetCellValue(3,2)+" " +self.tab3.myGrid.GetCellValue(3,3)+" " +self.tab3.myGrid.GetCellValue(3,4)+" " +self.tab3.myGrid.GetCellValue(3,5)) 
	flog.write("\nBlue_Grating " + self.tab3.myGrid.GetCellValue(4,0)+" " +self.tab3.myGrid.GetCellValue(4,1)+" " +self.tab3.myGrid.GetCellValue(4,2)+" " +self.tab3.myGrid.GetCellValue(4,3)+" " +self.tab3.myGrid.GetCellValue(4,4)+" " +self.tab3.myGrid.GetCellValue(4,5)) 
	flog.write("\nRed_Detector_Shift " + self.tab3.myGrid.GetCellValue(6,0)+" "+ self.tab3.myGrid.GetCellValue(6,1)) 
	flog.write("\nBlue_Detector_Shift " + self.tab3.myGrid.GetCellValue(6,2)+" "+ self.tab3.myGrid.GetCellValue(6,3)) 
	flog.write("\nRed_Detector_Rotation_(deg) " + self.tab3.myGrid.GetCellValue(8,0))
	flog.write("\nBlue_Detector_Rotation_(deg) " + self.tab3.myGrid.GetCellValue(8,1))
	flog.write("\nRed_RMS_Residual %f"  %(self.redrmsres))
	flog.write("\nBlue_RMS_Residual %f"  %(self.bluermsres))
	flog.write("\n")
	flog.close()
    def onFullButton(self, event):
	self.panell.drawredfull()
	self.panell2.drawbluefull()
    def onCmbxSelection(self, event):
	if (self.tab2.cmbx.GetValue()=='Reset'):
	        print '-> Resetting All Values'
        	for i in range(5):
	            for j in range(6):
	               self.tab2.myGrid.SetCellValue(i,j,'0')

	        for i in range(6,9,2):
		    for j in range(2):
	               self.tab2.myGrid.SetCellValue(i,j,'0')
	        self.tab2.myGrid.SetCellValue(6,2,'0')   
	        self.tab2.myGrid.SetCellValue(6,3,'0')     
	

	if (self.tab2.cmbx.GetValue()=='Random'):
	        print '-> Resetting All Values to Random'
		self.disrnge=int(self.tab2.selectDisRange.GetValue())
		self.rotrnge=int(self.tab2.selectRotRange.GetValue())
        	for i in range(5):
	            for j in range(3):
	               self.tab2.myGrid.SetCellValue(i,j,str("%0.4f" %np.random.uniform(-1*self.disrnge, self.disrnge)))
        	for i in range(5):
	            for j in range(3,6,1):
	               self.tab2.myGrid.SetCellValue(i,j,str("%0.4f" %np.random.uniform(-1*self.rotrnge, self.rotrnge)))
	        
	        for j in range(4):
	               self.tab2.myGrid.SetCellValue(6,j,str("%0.4f" %np.random.uniform(-0, 0))) 
	        self.tab2.myGrid.SetCellValue(8,0,str("%0.4f" %np.random.uniform(-0, 0)))   
	        self.tab2.myGrid.SetCellValue(8,1,str("%0.4f" %np.random.uniform(-0, 0)))     
	if (self.tab2.cmbx.GetValue()=='Collimator'):
	        print '-> Resetting Collimator Values to 500'
	    	for i in range(5):
	          for j in range(6):
	           self.tab2.myGrid.SetCellValue(i,j,'0')
		for j in range(6):
	           self.tab2.myGrid.SetCellValue(0,j,'500')
		  
	if (self.tab2.cmbx.GetValue()=='Grating'):
	        print '-> Resetting grating Values to 500'
	    	for i in range(5):
	          for j in range(6):
	           self.tab2.myGrid.SetCellValue(i,j,'0')
		for j in range(6):
	           self.tab2.myGrid.SetCellValue(3,j,'500')
		   self.tab2.myGrid.SetCellValue(4,j,'500')	
	if (self.tab2.cmbx.GetValue()=='All'):
	        print '-> Resetting All Values to 500'
	    	for i in range(5):
	          for j in range(6):
	           self.tab2.myGrid.SetCellValue(i,j,'500')
		for j in range(6):
	           self.tab2.myGrid.SetCellValue(2,j,'500')
		   
		   	   
	#self.panell.drawredfull()
	#self.panell2.drawbluefull()
       
    #----------------------------------------------------------------------
    
class MyFrame(wx.Frame):
    """"""
 
    #----------------------------------------------------------------------
    def __init__(self):
        """Constructor"""
        wx.Frame.__init__(self, parent=None, title="WFOS flexure Comp Tool", size=(1500, 1400))
	
	#p = wx.Panel(self)
        self.panel = MyPanel(self)
        self.Show()
 
if __name__ == "__main__":
    app = wx.App(False)
    frame = MyFrame()
    print"WFOS Flexure Compensation Tool"
    print"_______________________________"
    print"Ver 1.6.0 + 01/05/2017"
    print"Indian Institute of Astrophysics"   
    print 'arunsuryaoffice@gmail.com' 
    if config=='14out9':
	print "Current Config : 14out9 paraxial camera" 	
	red_coll_file='./14out9/red_coll_14.9.txt'
	red_mirr_file='./14out9/red_mirr_14.9.txt'
	red_gra_file='./14out9/red_gra_14.9.txt'
	blue_coll_file='./14out9/blue_coll_14.9.txt'
	blue_dchr_file='./14out9/blue_dchr_14.9.txt'
	blue_gra_file='./14out9/blue_gra_14.9.txt'	
    if config=='15out6':
	print "Current Config : 15out6 paraxial camera" 		
	red_coll_file='./15out6/red_coll_15.6.txt'
	red_mirr_file='./15out6/red_mirr_15.6.txt'
	red_gra_file='./15out6/red_gra_15.6.txt'
	blue_coll_file='./15out6/blue_coll_15.6.txt'
	blue_dchr_file='./15out6/blue_dchr_15.6.txt'
	blue_gra_file='./15out6/blue_gra_15.6.txt'	
    if config=='mobie':
	print "Current Config : MOBIE" 		
	red_coll_file='./mobie/red_coll.txt'
	red_mirr_file='./mobie/red_mirr.txt'
	red_gra_file='./mobie/red_gra.txt'
	blue_coll_file='./mobie/blue_coll.txt'
	blue_dchr_file='./mobie/blue_dchr.txt'
	blue_gra_file='./mobie/blue_gra.txt'	
    if config=='15out6cam':
	print "Current Config : 15out6 A7/BBA1 Camera" 		
	red_coll_file='./15out6cam/red_coll_15out6cam.txt'
	red_mirr_file='./15out6cam/red_mirr_15out6cam.txt'
	red_gra_file='./15out6cam/red_gra_15out6cam.txt'
	blue_coll_file='./15out6cam/blue_coll_15out6cam.txt'
	blue_dchr_file='./15out6cam/blue_dchr_15out6cam.txt'
	blue_gra_file='./15out6cam/blue_gra_15out6cam.txt'	
	
    red_datalist= pd.read_csv(red_coll_file,sep='\s+')
    red_ofield=[]
    red_infield=[]
    
    
    rx=np.array(red_datalist.X)
    ry=np.array(red_datalist.Y)
    rele=red_datalist.Element
    rcomp=red_datalist.Peturbation
    rdelx=np.array(red_datalist.DelX)
    rdely=np.array(red_datalist.DelY)
    
    
    red_x=[]
    red_y=[]
    red_coll_dx_x=[]
    red_coll_dx_y=[]
    red_coll_dy_x=[]
    red_coll_dy_y=[]
    red_coll_dz_x=[]
    red_coll_dz_y=[]
    red_coll_tx_x=[]
    red_coll_tx_y=[]
    red_coll_ty_x=[]
    red_coll_ty_y=[]
    red_coll_tz_x=[]
    red_coll_tz_y=[]
    
    
    for i in range(len(rx)):
	    if (i%6 == 0):
		    red_x.append(rx[i])
	    if (i%6 == 0):
		    red_y.append(ry[i])
	    if (i%6 == 0):
		    red_coll_dx_x.append(rdelx[i])
	    if (i%6 == 0):
		    red_coll_dx_y.append(rdely[i])
	    if (i%6 == 1):
		    red_coll_dy_x.append(rdelx[i])
	    if (i%6 == 1):
		    red_coll_dy_y.append(rdely[i])
	    if (i%6 == 2):
		    red_coll_dz_x.append(rdelx[i])
	    if (i%6 == 2):
		    red_coll_dz_y.append(rdely[i])
	    if (i%6 == 3):
		    red_coll_tx_x.append(rdelx[i])
	    if (i%6 == 3):
		    red_coll_tx_y.append(rdely[i])
	    if (i%6 == 4):
		    red_coll_ty_x.append(rdelx[i])
	    if (i%6 == 4):
		    red_coll_ty_y.append(rdely[i])
	    if (i%6 == 5):
		    red_coll_tz_x.append(rdelx[i])
	    if (i%6 == 5):
		    red_coll_tz_y.append(rdely[i])
    
    
    red_x=numpy.array(red_x)
    red_y=numpy.array(red_y)
    red_coll_dx_x=numpy.array(red_coll_dx_x)*1e-6
    red_coll_dx_y=numpy.array(red_coll_dx_y)*1e-6
    red_coll_dy_x=numpy.array(red_coll_dy_x)*1e-6
    red_coll_dy_y=numpy.array(red_coll_dy_y)*1e-6
    red_coll_dz_x=numpy.array(red_coll_dz_x)*1e-6
    red_coll_dz_y=numpy.array(red_coll_dz_y)*1e-6
    red_coll_tx_x=numpy.array(red_coll_tx_x)*1e-6
    red_coll_tx_y=numpy.array(red_coll_tx_y)*1e-6
    red_coll_ty_x=numpy.array(red_coll_ty_x)*1e-6
    red_coll_ty_y=numpy.array(red_coll_ty_y)*1e-6
    red_coll_tz_x=numpy.array(red_coll_tz_x)*1e-6
    red_coll_tz_y=numpy.array(red_coll_tz_y)*1e-6
    
    
    red_datalist= pd.read_csv(red_mirr_file,sep='\s+')
    red_ofield=[]
    red_infield=[]
    
    
    rmele=red_datalist.Element
    rmcomp=red_datalist.Peturbation
    rmdelx=np.array(red_datalist.DelX)
    rmdely=np.array(red_datalist.DelY)
    
    
    
    red_mirr_dx_x=[]
    red_mirr_dx_y=[]
    red_mirr_dy_x=[]
    red_mirr_dy_y=[]
    red_mirr_dz_x=[]
    red_mirr_dz_y=[]
    red_mirr_tx_x=[]
    red_mirr_tx_y=[]
    red_mirr_ty_x=[]
    red_mirr_ty_y=[]
    red_mirr_tz_x=[]
    red_mirr_tz_y=[]
    
    
    for i in range(len(rmdelx)):
	    if (i%6 == 0):
		    red_mirr_dx_x.append(rmdelx[i])
	    if (i%6 == 0):
		    red_mirr_dx_y.append(rmdely[i])
	    if (i%6 == 1):
		    red_mirr_dy_x.append(rmdelx[i])
	    if (i%6 == 1):
		    red_mirr_dy_y.append(rmdely[i])
	    if (i%6 == 2):
		    red_mirr_dz_x.append(rmdelx[i])
	    if (i%6 == 2):
		    red_mirr_dz_y.append(rmdely[i])
	    if (i%6 == 3):
		    red_mirr_tx_x.append(rmdelx[i])
	    if (i%6 == 3):
		    red_mirr_tx_y.append(rmdely[i])
	    if (i%6 == 4):
		    red_mirr_ty_x.append(rmdelx[i])
	    if (i%6 == 4):
		    red_mirr_ty_y.append(rmdely[i])
	    if (i%6 == 5):
		    red_mirr_tz_x.append(rmdelx[i])
	    if (i%6 == 5):
		    red_mirr_tz_y.append(rmdely[i])
    
    
    red_mirr_dx_x=numpy.array(red_mirr_dx_x)*1e-6
    red_mirr_dx_y=numpy.array(red_mirr_dx_y)*1e-6
    red_mirr_dy_x=numpy.array(red_mirr_dy_x)*1e-6
    red_mirr_dy_y=numpy.array(red_mirr_dy_y)*1e-6
    red_mirr_dz_x=numpy.array(red_mirr_dz_x)*1e-6
    red_mirr_dz_y=numpy.array(red_mirr_dz_y)*1e-6
    red_mirr_tx_x=numpy.array(red_mirr_tx_x)*1e-6
    red_mirr_tx_y=numpy.array(red_mirr_tx_y)*1e-6
    red_mirr_ty_x=numpy.array(red_mirr_ty_x)*1e-6
    red_mirr_ty_y=numpy.array(red_mirr_ty_y)*1e-6
    red_mirr_tz_x=numpy.array(red_mirr_tz_x)*1e-6
    red_mirr_tz_y=numpy.array(red_mirr_tz_y)*1e-6
    
    
    
    red_datalist= pd.read_csv(red_gra_file,sep='\s+')
    rgele=red_datalist.Element
    rgcomp=red_datalist.Peturbation
    rgdelx=np.array(red_datalist.DelX)
    rgdely=np.array(red_datalist.DelY)
    
    red_gra_dx_x=[]
    red_gra_dx_y=[]
    red_gra_dy_x=[]
    red_gra_dy_y=[]
    red_gra_dz_x=[]
    red_gra_dz_y=[]
    red_gra_tx_x=[]
    red_gra_tx_y=[]
    red_gra_ty_x=[]
    red_gra_ty_y=[]
    red_gra_tz_x=[]
    red_gra_tz_y=[]
    
    
    for i in range(len(rgdelx)):
	    if (i%6 == 0):
		    red_gra_dx_x.append(rgdelx[i])
	    if (i%6 == 0):
		    red_gra_dx_y.append(rgdely[i])
	    if (i%6 == 1):
		    red_gra_dy_x.append(rgdelx[i])
	    if (i%6 == 1):
		    red_gra_dy_y.append(rgdely[i])
	    if (i%6 == 2):
		    red_gra_dz_x.append(rgdelx[i])
	    if (i%6 == 2):
		    red_gra_dz_y.append(rgdely[i])
	    if (i%6 == 3):
		    red_gra_tx_x.append(rgdelx[i])
	    if (i%6 == 3):
		    red_gra_tx_y.append(rgdely[i])
	    if (i%6 == 4):
		    red_gra_ty_x.append(rgdelx[i])
	    if (i%6 == 4):
		    red_gra_ty_y.append(rgdely[i])
	    if (i%6 == 5):
		    red_gra_tz_x.append(rgdelx[i])
	    if (i%6 == 5):
		    red_gra_tz_y.append(rgdely[i])
    
    red_gra_dx_x=numpy.array(red_gra_dx_x)*1e-6
    red_gra_dx_y=numpy.array(red_gra_dx_y)*1e-6
    red_gra_dy_x=numpy.array(red_gra_dy_x)*1e-6
    red_gra_dy_y=numpy.array(red_gra_dy_y)*1e-6
    red_gra_dz_x=numpy.array(red_gra_dz_x)*1e-6
    red_gra_dz_y=numpy.array(red_gra_dz_y)*1e-6
    red_gra_tx_x=numpy.array(red_gra_tx_x)*1e-6
    red_gra_tx_y=numpy.array(red_gra_tx_y)*1e-6
    red_gra_ty_x=numpy.array(red_gra_ty_x)*1e-6
    red_gra_ty_y=numpy.array(red_gra_ty_y)*1e-6
    red_gra_tz_x=numpy.array(red_gra_tz_x)*1e-6
    red_gra_tz_y=numpy.array(red_gra_tz_y)*1e-6
    
    
    
    
    
    #input bluedata
    
    blue_datalist= pd.read_csv(blue_coll_file,sep='\s+')
    blue_ofield=[]
    blue_infield=[]
    bx=np.array(blue_datalist.X)
    by=np.array(blue_datalist.Y)
    bele=blue_datalist.Element
    bcomp=blue_datalist.Peturbation
    bdelx=np.array(blue_datalist.DelX)
    bdely=np.array(blue_datalist.DelY)
    
    
    blue_x=[]
    blue_y=[]
    blue_coll_dx_x=[]
    blue_coll_dx_y=[]
    blue_coll_dy_x=[]
    blue_coll_dy_y=[]
    blue_coll_dz_x=[]
    blue_coll_dz_y=[]
    blue_coll_tx_x=[]
    blue_coll_tx_y=[]
    blue_coll_ty_x=[]
    blue_coll_ty_y=[]
    blue_coll_tz_x=[]
    blue_coll_tz_y=[]
    
    for i in range(len(bx)):
	    if (i%6 == 0):
		    blue_x.append(bx[i])
	    if (i%6 == 0):
		    blue_y.append(by[i])
	    if (i%6 == 0):
		    blue_coll_dx_x.append(bdelx[i])
	    if (i%6 == 0):
		    blue_coll_dx_y.append(bdely[i])
	    if (i%6 == 1):
		    blue_coll_dy_x.append(bdelx[i])
	    if (i%6 == 1):
		    blue_coll_dy_y.append(bdely[i])
	    if (i%6 == 2):
		    blue_coll_dz_x.append(bdelx[i])
	    if (i%6 == 2):
		    blue_coll_dz_y.append(bdely[i])
	    if (i%6 == 3):
		    blue_coll_tx_x.append(bdelx[i])
	    if (i%6 == 3):
		    blue_coll_tx_y.append(bdely[i])
	    if (i%6 == 4):
		    blue_coll_ty_x.append(bdelx[i])
	    if (i%6 == 4):
		    blue_coll_ty_y.append(bdely[i])
	    if (i%6 == 5):
		    blue_coll_tz_x.append(bdelx[i])
	    if (i%6 == 5):
		    blue_coll_tz_y.append(bdely[i])
    
    blue_x=numpy.array(blue_x)
    blue_y=numpy.array(blue_y)
    blue_coll_dx_x=numpy.array(blue_coll_dx_x)*1e-6
    blue_coll_dx_y=numpy.array(blue_coll_dx_y)*1e-6
    blue_coll_dy_x=numpy.array(blue_coll_dy_x)*1e-6
    blue_coll_dy_y=numpy.array(blue_coll_dy_y)*1e-6
    blue_coll_dz_x=numpy.array(blue_coll_dz_x)*1e-6
    blue_coll_dz_y=numpy.array(blue_coll_dz_y)*1e-6
    blue_coll_tx_x=numpy.array(blue_coll_tx_x)*1e-6
    blue_coll_tx_y=numpy.array(blue_coll_tx_y)*1e-6
    blue_coll_ty_x=numpy.array(blue_coll_ty_x)*1e-6
    blue_coll_ty_y=numpy.array(blue_coll_ty_y)*1e-6
    blue_coll_tz_x=numpy.array(blue_coll_tz_x)*1e-6
    blue_coll_tz_y=numpy.array(blue_coll_tz_y)*1e-6
    
    
    
    
    
    blue_datalist= pd.read_csv(blue_dchr_file,sep='\s+')
    bdele=blue_datalist.Element
    bdcomp=blue_datalist.Peturbation
    bddelx=np.array(blue_datalist.DelX)
    bddely=np.array(blue_datalist.DelY)
    
    blue_dchr_dx_x=[]
    blue_dchr_dx_y=[]
    blue_dchr_dy_x=[]
    blue_dchr_dy_y=[]
    blue_dchr_dz_x=[]
    blue_dchr_dz_y=[]
    blue_dchr_tx_x=[]
    blue_dchr_tx_y=[]
    blue_dchr_ty_x=[]
    blue_dchr_ty_y=[]
    blue_dchr_tz_x=[]
    blue_dchr_tz_y=[]
    
    
    for i in range(len(bddelx)):
	    if (i%6 == 0):
		    blue_dchr_dx_x.append(bddelx[i])
	    if (i%6 == 0):
		    blue_dchr_dx_y.append(bddely[i])
	    if (i%6 == 1):
		    blue_dchr_dy_x.append(bddelx[i])
	    if (i%6 == 1):
		    blue_dchr_dy_y.append(bddely[i])
	    if (i%6 == 2):
		    blue_dchr_dz_x.append(bddelx[i])
	    if (i%6 == 2):
		    blue_dchr_dz_y.append(bddely[i])
	    if (i%6 == 3):
		    blue_dchr_tx_x.append(bddelx[i])
	    if (i%6 == 3):
		    blue_dchr_tx_y.append(bddely[i])
	    if (i%6 == 4):
		    blue_dchr_ty_x.append(bddelx[i])
	    if (i%6 == 4):
		    blue_dchr_ty_y.append(bddely[i])
	    if (i%6 == 5):
		    blue_dchr_tz_x.append(bddelx[i])
	    if (i%6 == 5):
		    blue_dchr_tz_y.append(bddely[i])
    
    blue_dchr_dx_x=numpy.array(blue_dchr_dx_x)*1e-6
    blue_dchr_dx_y=numpy.array(blue_dchr_dx_y)*1e-6
    blue_dchr_dy_x=numpy.array(blue_dchr_dy_x)*1e-6
    blue_dchr_dy_y=numpy.array(blue_dchr_dy_y)*1e-6
    blue_dchr_dz_x=numpy.array(blue_dchr_dz_x)*1e-6
    blue_dchr_dz_y=numpy.array(blue_dchr_dz_y)*1e-6
    blue_dchr_tx_x=numpy.array(blue_dchr_tx_x)*1e-6
    blue_dchr_tx_y=numpy.array(blue_dchr_tx_y)*1e-6
    blue_dchr_ty_x=numpy.array(blue_dchr_ty_x)*1e-6
    blue_dchr_ty_y=numpy.array(blue_dchr_ty_y)*1e-6
    blue_dchr_tz_x=numpy.array(blue_dchr_tz_x)*1e-6
    blue_dchr_tz_y=numpy.array(blue_dchr_tz_y)*1e-6
    
    
    
    
    blue_datalist= pd.read_csv(blue_gra_file,sep='\s+')
    bgele=blue_datalist.Element
    bgcomp=blue_datalist.Peturbation
    bgdelx=np.array(blue_datalist.DelX)
    bgdely=np.array(blue_datalist.DelY)
    
    blue_gra_dx_x=[]
    blue_gra_dx_y=[]
    blue_gra_dy_x=[]
    blue_gra_dy_y=[]
    blue_gra_dz_x=[]
    blue_gra_dz_y=[]
    blue_gra_tx_x=[]
    blue_gra_tx_y=[]
    blue_gra_ty_x=[]
    blue_gra_ty_y=[]
    blue_gra_tz_x=[]
    blue_gra_tz_y=[]
    
    
    for i in range(len(bgdelx)):
	    if (i%6 == 0):
		    blue_gra_dx_x.append(bgdelx[i])
	    if (i%6 == 0):
		    blue_gra_dx_y.append(bgdely[i])
	    if (i%6 == 1):
		    blue_gra_dy_x.append(bgdelx[i])
	    if (i%6 == 1):
		    blue_gra_dy_y.append(bgdely[i])
	    if (i%6 == 2):
		    blue_gra_dz_x.append(bgdelx[i])
	    if (i%6 == 2):
		    blue_gra_dz_y.append(bgdely[i])
	    if (i%6 == 3):
		    blue_gra_tx_x.append(bgdelx[i])
	    if (i%6 == 3):
		    blue_gra_tx_y.append(bgdely[i])
	    if (i%6 == 4):
		    blue_gra_ty_x.append(bgdelx[i])
	    if (i%6 == 4):
		    blue_gra_ty_y.append(bgdely[i])
	    if (i%6 == 5):
		    blue_gra_tz_x.append(bgdelx[i])
	    if (i%6 == 5):
		    blue_gra_tz_y.append(bgdely[i])
    
    blue_gra_dx_x=numpy.array(blue_gra_dx_x)*1e-6
    blue_gra_dx_y=numpy.array(blue_gra_dx_y)*1e-6
    blue_gra_dy_x=numpy.array(blue_gra_dy_x)*1e-6
    blue_gra_dy_y=numpy.array(blue_gra_dy_y)*1e-6
    blue_gra_dz_x=numpy.array(blue_gra_dz_x)*1e-6
    blue_gra_dz_y=numpy.array(blue_gra_dz_y)*1e-6
    blue_gra_tx_x=numpy.array(blue_gra_tx_x)*1e-6
    blue_gra_tx_y=numpy.array(blue_gra_tx_y)*1e-6
    blue_gra_ty_x=numpy.array(blue_gra_ty_x)*1e-6
    blue_gra_ty_y=numpy.array(blue_gra_ty_y)*1e-6
    blue_gra_tz_x=numpy.array(blue_gra_tz_x)*1e-6
    blue_gra_tz_y=numpy.array(blue_gra_tz_y)*1e-6
    
    
    #Input: alldata finished
    
    
    #SelectField
    
    
    for i in range(len(red_x)):
	    if (abs(red_x[i]) < 0):   # (abs(red_x[i]) < 70):
		    red_infield.append(i)
	    else:
		    red_ofield.append(i)
    red_allfield=range(len(red_x))
    
    redo_x=red_x[red_ofield]
    redo_y=red_y[red_ofield]
    redo_coll_dx_x=red_coll_dx_x[red_ofield]
    redo_coll_dx_y= red_coll_dx_y[red_ofield]
    redo_coll_dy_x= red_coll_dy_x[red_ofield]
    redo_coll_dy_y= red_coll_dy_y[red_ofield]
    redo_coll_dz_x= red_coll_dz_x[red_ofield]
    redo_coll_dz_y= red_coll_dz_y[red_ofield]
    redo_coll_tx_x= red_coll_tx_x[red_ofield]
    redo_coll_tx_y= red_coll_tx_y[red_ofield]
    redo_coll_ty_x= red_coll_ty_x[red_ofield]
    redo_coll_ty_y= red_coll_ty_y[red_ofield]
    redo_coll_tz_x= red_coll_tz_x[red_ofield]
    redo_coll_tz_y= red_coll_tz_y[red_ofield]
    
    redo_mirr_dx_x=red_mirr_dx_x[red_ofield]
    redo_mirr_dx_y= red_mirr_dx_y[red_ofield]
    redo_mirr_dy_x= red_mirr_dy_x[red_ofield]
    redo_mirr_dy_y= red_mirr_dy_y[red_ofield]
    redo_mirr_dz_x= red_mirr_dz_x[red_ofield]
    redo_mirr_dz_y= red_mirr_dz_y[red_ofield]
    redo_mirr_tx_x= red_mirr_tx_x[red_ofield]
    redo_mirr_tx_y= red_mirr_tx_y[red_ofield]
    redo_mirr_ty_x= red_mirr_ty_x[red_ofield]
    redo_mirr_ty_y= red_mirr_ty_y[red_ofield]
    redo_mirr_tz_x= red_mirr_tz_x[red_ofield]
    redo_mirr_tz_y= red_mirr_tz_y[red_ofield]
    
    redo_gra_dx_x=red_gra_dx_x[red_ofield]
    redo_gra_dx_y= red_gra_dx_y[red_ofield]
    redo_gra_dy_x= red_gra_dy_x[red_ofield]
    redo_gra_dy_y= red_gra_dy_y[red_ofield]
    redo_gra_dz_x= red_gra_dz_x[red_ofield]
    redo_gra_dz_y= red_gra_dz_y[red_ofield]
    redo_gra_tx_x= red_gra_tx_x[red_ofield]
    redo_gra_tx_y= red_gra_tx_y[red_ofield]
    redo_gra_ty_x= red_gra_ty_x[red_ofield]
    redo_gra_ty_y= red_gra_ty_y[red_ofield]
    redo_gra_tz_x= red_gra_tz_x[red_ofield]
    redo_gra_tz_y= red_gra_tz_y[red_ofield]
    
    
    
    for i in range(len(blue_x)):
	    if ((abs(blue_x[i]) < 0)): # ((blue_x[i]) < 75)):
		    blue_infield.append(i)
	    else:
		    blue_ofield.append(i)
    blue_allfield=range(len(blue_x))
    
    blueo_x=blue_x[blue_ofield]
    blueo_y=blue_y[blue_ofield]
    blueo_coll_dx_x=blue_coll_dx_x[blue_ofield]
    blueo_coll_dx_y= blue_coll_dx_y[blue_ofield]
    blueo_coll_dy_x= blue_coll_dy_x[blue_ofield]
    blueo_coll_dy_y= blue_coll_dy_y[blue_ofield]
    blueo_coll_dz_x= blue_coll_dz_x[blue_ofield]
    blueo_coll_dz_y= blue_coll_dz_y[blue_ofield]
    blueo_coll_tx_x= blue_coll_tx_x[blue_ofield]
    blueo_coll_tx_y= blue_coll_tx_y[blue_ofield]
    blueo_coll_ty_x= blue_coll_ty_x[blue_ofield]
    blueo_coll_ty_y= blue_coll_ty_y[blue_ofield]
    blueo_coll_tz_x= blue_coll_tz_x[blue_ofield]
    blueo_coll_tz_y= blue_coll_tz_y[blue_ofield]
    
    blueo_dchr_dx_x=blue_dchr_dx_x[blue_ofield]
    blueo_dchr_dx_y= blue_dchr_dx_y[blue_ofield]
    blueo_dchr_dy_x= blue_dchr_dy_x[blue_ofield]
    blueo_dchr_dy_y= blue_dchr_dy_y[blue_ofield]
    blueo_dchr_dz_x= blue_dchr_dz_x[blue_ofield]
    blueo_dchr_dz_y= blue_dchr_dz_y[blue_ofield]
    blueo_dchr_tx_x= blue_dchr_tx_x[blue_ofield]
    blueo_dchr_tx_y= blue_dchr_tx_y[blue_ofield]
    blueo_dchr_ty_x= blue_dchr_ty_x[blue_ofield]
    blueo_dchr_ty_y= blue_dchr_ty_y[blue_ofield]
    blueo_dchr_tz_x= blue_dchr_tz_x[blue_ofield]
    blueo_dchr_tz_y= blue_dchr_tz_y[blue_ofield]
    
    blueo_gra_dx_x=blue_gra_dx_x[blue_ofield]
    blueo_gra_dx_y= blue_gra_dx_y[blue_ofield]
    blueo_gra_dy_x= blue_gra_dy_x[blue_ofield]
    blueo_gra_dy_y= blue_gra_dy_y[blue_ofield]
    blueo_gra_dz_x= blue_gra_dz_x[blue_ofield]
    blueo_gra_dz_y= blue_gra_dz_y[blue_ofield]
    blueo_gra_tx_x= blue_gra_tx_x[blue_ofield]
    blueo_gra_tx_y= blue_gra_tx_y[blue_ofield]
    blueo_gra_ty_x= blue_gra_ty_x[blue_ofield]
    blueo_gra_ty_y= blue_gra_ty_y[blue_ofield]
    blueo_gra_tz_x= blue_gra_tz_x[blue_ofield]
    blueo_gra_tz_y= blue_gra_tz_y[blue_ofield]
   
    blue_ofield=[]
    blue_infield=[]
    red_ofield=[]
    red_infield=[]
    
    for i in range(len(red_x)):
	    if (abs(red_x[i]) < 70):   # (abs(red_x[i]) < 70):
		    red_infield.append(i)
	    else:
		    red_ofield.append(i)
    red_allfield=range(len(red_x))
    
    redo1_x=red_x[red_ofield]
    redo1_y=red_y[red_ofield]
    redo1_coll_dx_x=red_coll_dx_x[red_ofield]
    redo1_coll_dx_y= red_coll_dx_y[red_ofield]
    redo1_coll_dy_x= red_coll_dy_x[red_ofield]
    redo1_coll_dy_y= red_coll_dy_y[red_ofield]
    redo1_coll_dz_x= red_coll_dz_x[red_ofield]
    redo1_coll_dz_y= red_coll_dz_y[red_ofield]
    redo1_coll_tx_x= red_coll_tx_x[red_ofield]
    redo1_coll_tx_y= red_coll_tx_y[red_ofield]
    redo1_coll_ty_x= red_coll_ty_x[red_ofield]
    redo1_coll_ty_y= red_coll_ty_y[red_ofield]
    redo1_coll_tz_x= red_coll_tz_x[red_ofield]
    redo1_coll_tz_y= red_coll_tz_y[red_ofield]
    
    redo1_mirr_dx_x=red_mirr_dx_x[red_ofield]
    redo1_mirr_dx_y= red_mirr_dx_y[red_ofield]
    redo1_mirr_dy_x= red_mirr_dy_x[red_ofield]
    redo1_mirr_dy_y= red_mirr_dy_y[red_ofield]
    redo1_mirr_dz_x= red_mirr_dz_x[red_ofield]
    redo1_mirr_dz_y= red_mirr_dz_y[red_ofield]
    redo1_mirr_tx_x= red_mirr_tx_x[red_ofield]
    redo1_mirr_tx_y= red_mirr_tx_y[red_ofield]
    redo1_mirr_ty_x= red_mirr_ty_x[red_ofield]
    redo1_mirr_ty_y= red_mirr_ty_y[red_ofield]
    redo1_mirr_tz_x= red_mirr_tz_x[red_ofield]
    redo1_mirr_tz_y= red_mirr_tz_y[red_ofield]
    
    redo1_gra_dx_x=red_gra_dx_x[red_ofield]
    redo1_gra_dx_y= red_gra_dx_y[red_ofield]
    redo1_gra_dy_x= red_gra_dy_x[red_ofield]
    redo1_gra_dy_y= red_gra_dy_y[red_ofield]
    redo1_gra_dz_x= red_gra_dz_x[red_ofield]
    redo1_gra_dz_y= red_gra_dz_y[red_ofield]
    redo1_gra_tx_x= red_gra_tx_x[red_ofield]
    redo1_gra_tx_y= red_gra_tx_y[red_ofield]
    redo1_gra_ty_x= red_gra_ty_x[red_ofield]
    redo1_gra_ty_y= red_gra_ty_y[red_ofield]
    redo1_gra_tz_x= red_gra_tz_x[red_ofield]
    redo1_gra_tz_y= red_gra_tz_y[red_ofield]
    
    
    
    for i in range(len(blue_x)):
	    if ((abs(blue_x[i]) < 75)): # ((blue_x[i]) < 75)):
		    blue_infield.append(i)
	    else:
		    blue_ofield.append(i)
    blue_allfield=range(len(blue_x))
    
    blueo1_x=blue_x[blue_ofield]
    blueo1_y=blue_y[blue_ofield]
    blueo1_coll_dx_x=blue_coll_dx_x[blue_ofield]
    blueo1_coll_dx_y= blue_coll_dx_y[blue_ofield]
    blueo1_coll_dy_x= blue_coll_dy_x[blue_ofield]
    blueo1_coll_dy_y= blue_coll_dy_y[blue_ofield]
    blueo1_coll_dz_x= blue_coll_dz_x[blue_ofield]
    blueo1_coll_dz_y= blue_coll_dz_y[blue_ofield]
    blueo1_coll_tx_x= blue_coll_tx_x[blue_ofield]
    blueo1_coll_tx_y= blue_coll_tx_y[blue_ofield]
    blueo1_coll_ty_x= blue_coll_ty_x[blue_ofield]
    blueo1_coll_ty_y= blue_coll_ty_y[blue_ofield]
    blueo1_coll_tz_x= blue_coll_tz_x[blue_ofield]
    blueo1_coll_tz_y= blue_coll_tz_y[blue_ofield]
    
    blueo1_dchr_dx_x=blue_dchr_dx_x[blue_ofield]
    blueo1_dchr_dx_y= blue_dchr_dx_y[blue_ofield]
    blueo1_dchr_dy_x= blue_dchr_dy_x[blue_ofield]
    blueo1_dchr_dy_y= blue_dchr_dy_y[blue_ofield]
    blueo1_dchr_dz_x= blue_dchr_dz_x[blue_ofield]
    blueo1_dchr_dz_y= blue_dchr_dz_y[blue_ofield]
    blueo1_dchr_tx_x= blue_dchr_tx_x[blue_ofield]
    blueo1_dchr_tx_y= blue_dchr_tx_y[blue_ofield]
    blueo1_dchr_ty_x= blue_dchr_ty_x[blue_ofield]
    blueo1_dchr_ty_y= blue_dchr_ty_y[blue_ofield]
    blueo1_dchr_tz_x= blue_dchr_tz_x[blue_ofield]
    blueo1_dchr_tz_y= blue_dchr_tz_y[blue_ofield]
    
    blueo1_gra_dx_x=blue_gra_dx_x[blue_ofield]
    blueo1_gra_dx_y= blue_gra_dx_y[blue_ofield]
    blueo1_gra_dy_x= blue_gra_dy_x[blue_ofield]
    blueo1_gra_dy_y= blue_gra_dy_y[blue_ofield]
    blueo1_gra_dz_x= blue_gra_dz_x[blue_ofield]
    blueo1_gra_dz_y= blue_gra_dz_y[blue_ofield]
    blueo1_gra_tx_x= blue_gra_tx_x[blue_ofield]
    blueo1_gra_tx_y= blue_gra_tx_y[blue_ofield]
    blueo1_gra_ty_x= blue_gra_ty_x[blue_ofield]
    blueo1_gra_ty_y= blue_gra_ty_y[blue_ofield]
    blueo1_gra_tz_x= blue_gra_tz_x[blue_ofield]
    blueo1_gra_tz_y= blue_gra_tz_y[blue_ofield]


    app.MainLoop() 






