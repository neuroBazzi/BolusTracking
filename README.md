# BolusTracking
Matlab GUI to load multi-images TIF stacks, draw ROIs and fit the signal time course from the ROIs with a gamma variate function. 

Introduction:

This code is to extract information from timeseries of 2D images from experiments in which a fluorescent dye bolus passes through the vasculature. Example data is provided.

Instructions:

1) In your matlab command windows, type BolusTrack;
2) When the GUI opens, use the "load data" button to navigate your system and locate the file (use the BolusDataTest file as an example);
3) the code will load the images, and display the Maximum Intensity Projection (MIP) of the stack (after a 3x3 gaussian blurring) in the bottom part of the GUI;
4) Below the load button, input the frame rate and the "number of vessels" (i.e. ROIs) you want to analyze, then press SelectROIs;
5) This will enable to select the ROIs directly from the image below. Click on one poin to begin drawing a polygon, once you're done, close the polygon by double-clicking on the first point;
6)  Once you're done drawing all the ROIs indicated in the "number of vessels" field, the average signal within the first ROI will appear in the centre-top of the figure. At this point you need to initialize the parameters for the gamma variate fitting. You need to input 6 parameters, extracted from 4 features of the signal intensity time-course. The 6 parameters are initialized as NaN, and you can change them easily by selecting the corresponding point from the figure. To do so, move the mouse cursor over the peak of the signal and click on it. Matlab will display the X and Y cohordinates of the point. To acquire those two parameters, click on the editable field of the "Amplitude" and press "enter". Then do the same in the "time to peak" field. Then move the cursor at the signal take off, click on the trace, then move to the "fit start" field and press enter. To the same for the end of the fit. The FWHM and baseline drift are estimated automatically by the previous parameters,so just click with the cursor on the editbale field of FWHM and press enter and then in the baseline shift editable field and press enter. once all the parameters have non-zeros (and non-NaN) vlaues, click on the Fit Gamma function button. f the initial parameters have been entered correctly, the fitted trace will be overlaid in the plot in red. You can improve the fit by adjusing the initial parameters, mainly the fit start and end. If you change those, temeber to refresh the FWHM and baseline shift fields (by clicking on the field with the mouse and pressing enter) as those do not update automatically. Once you're satisfied with the fit, press Save value/next trace, this will store all the fitting parameters, clear the figure, re-set the initial parameters to NaN and plot the next trace.
7) Repeat the steps in point 6 for all the ROIs;
b) If some of the signal time-courses are too noisy to be fitted, just press on the "Save values/next trace" button, the program will handle the NaN and drop them before exporting;
8) Once you're done fitting the signal timecourse from all the ROIs, the GUI will display an empty figure. At this point, press the Export Results button. This will pop up a window that will allow you to select the destination and name of the output file. The program will output two files: the first file contains the initial fitting parameters, estimated fitted parameters and time to peak (with 95% C.I.) in Excel format, the seocnd file is a matlab file (.m) which will contain all the other fitting paramers (it will be named BolusTrack and will be saved in the same folder selected by the user). 

