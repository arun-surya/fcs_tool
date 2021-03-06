Flexure Compensation Simulation Tool
====================================
Introduction
------------
The FCS tool takes sensitivity maps of Optical elements in the design and allows to simulate manual flexure and corresponding compensation from active elements.

The Input sensitivity files are given in the zip file.

blue_coll.txt : Collimator Sensitivity in Blue Channel 
blue_gra.txt  : Grating Sensitivity in Blue Channel 
red_gra.txt   : Grating Sensitivity in Red Channel 
blue_dchr.txt : Dichroic Sensitivity in Blue Channel 
red_coll.txt  : Collimator Sensitivity in Red Camera 
red_mirr.txt  : Fold Mirror Sensitivity in Red Channel 


The tool allows to enter the dx,dy,dz - displacement in x,y,z  and tx,ty,tz - tilt in x,y,z axis.

dx,dy.dz units are in .001 mm
tx,ty,tz units are in .00001 degree
Range of values are (-2000 <-> 2000 units)

You can also add Detector Translation (x,y) (units in .001 mm)
Detector Rotation (in Degree) 


Installation
------------

The tool is written in python and works with following libraries

wxpython (3.0.2.0)
scipy  (0.14.0)
pandas
numpy  (1.11.2)
threading
matplotlib  (1.5.3) 

It should work with new versions as well, but try for the versions given in case things don't work.

Compensation after flexure is done using scipy.optimize module using basin hopping. This routine is present in new versions of scipy as well. In the compensation tab you can select active elements to compensate for the flexure and view the optimized compensation from the algorithm. the results are displayed in the plot window.

Update 19/01/2017
-----------------
1. Range for Random flex values can be specified in the Range field next to Flex Button
2. The compensation values for the active elements show up in a new tab called Compensation Values

I have made a few more changes so there is less hard coded elements and it is more customizable

3. There is a Flex.View.Scale field in Flexure Tab. Default value is 1000. This is the amount by which spot movements are exaggerated  in plot window after each flex. You can check it by changing the value and flexing it again. The visual scale of spot movements change. The pixel box for reference changes also.

4. There is a Residual.View.Scale in Compensation tab. This has the same effect as the previous field. It is the scale by which residuals are enhanced in the plot. after each compensation, you can try changing it and pressing the button-> Show Residuals. Default value is 10000.

5. The number of iterations for the algorithm can be specified in the field in Compensation Tab. 

FCS Tool 1.2.1 Update 06/02/2017
--------------------------------
1. A new save feature is introduced which allows you to save the values in the Flex and Compensation Tabs. To use this , after flex operation, enter a filename and press the save button in the Flex Tab. The same procedure is used to save the Compensation values. In the compensation file the rms residuals are also logged in addition to the compensation flexures.

2. The feature to select Edge Field or Full Field for optimization is provided in the software now. You have to select the checkbox in Compensation Tab (' Use edge field for optimization'). After checking you have to Flex and Compensate to obtain the results with new setting. Default optimization uses full field of View. 

FCS Tool 1.3 Update 10/03/2017
------------------------------
1. The new 15out6 and 14out9 configurations have been used for the simulations. 
2. The configuration can be changed in the code at line 36 changing '14.9' to '15.6' and vice versa.


FCS Tool 1.4 Update 12/04/2017
------------------------------
1. Four configurations can be selected. Currently this has to be done before the fcstool is run . One need to open the fcs_tool.py script with an editor and uncomment the required line in  line no 4,5,6,7. Set 'config =' accordingly to select the configuration. Given options are '14out9','15out6','mobie','15out6cam'.

FCS Tool 1.5 Update 12/04/2017
------------------------------
1. Four configurations can be selected from the Select Config dropbox. Currently these configurations are [15out6prx,15out6cam,14out9prx,MOBIE]
2. The flexure files can be saved and loaded back from text files. These fields are just below the flexure box. 
3. Seperate Randomiser ranges can be set for the dx,dy,dz elements and tx,ty,tz . There are two new fields for this in the flexure panel.  

 
FCS Tool 1.6 Update 01/05/2017
------------------------------
1. Component controls for collimator and mirror in compensation. The new feature adds ability to select active components (dx,dy,dz,tx,ty,tz) for   collimator and fold mirror for the compensation routine. This feature is currently a test feature and will be available for all the active components from next version.

2. New analysis routines. The earlier'Compensation Values' Tab is now 'Analyse' Tab. One can view the compensation values in this tab as earlier. There is an extra Analyse Button in the bottom with check list options of 'Flexure', 'Compensation' and 'Resultent'. The Analyse Routine allows to see the  mean motion component of each element in the movement in the detector. According to the checkbox controls it shows the flexure components and compensation components together or seperately. The resultant mean motion is also shown if the 'Resultant' box is checked.
