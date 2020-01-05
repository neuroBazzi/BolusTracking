# BolusTracking
Matlab GUI to load multi-images TIF stacks, draw ROIs and fit the signal time course from the ROIs with a gamma variate function. 

Introduction:

This Matlab program allows to select multiple ROIs from .tiff image stacks, calculate the time course of the average signal intensity within the ROIs and fit them with a gamma-variate function. It then estimate the Onset Time of the fluorescence in the field of view and estimates the time-to-peak of the fluorescence from the gamma function. It allows to specify the experimental condition and to store the ROIs and reload them onto a different image stack (for longitudinal experiments). It exports the results in .csv (compatible with R and Excel) file and in a .mat file (a structure variable containing also all the raw signal time course, the masks, and full fitting results).


Instructions:


1)	Launch the GUI by typing BolusTrack in Matlab command window;
2)	Use the Load Data button to navigate the computer folders, find the file and select it; once the file is loaded, the Maximum Intensity Projection of the gaussian blurred stack will appear in the bottom of the screen;
3)	Input the Frame rate of the acquisition;
4)	There are two options:
a.	Start a new analysis: input the number of ROIs in the editable field called “Number of Vessels” and then click on select ROIs. This will enable the interaction with the image. Star drawing the ROIs on the picture, to close an ROI return to the first point and single click on it. Once an ROI is close you can drag it by positioning the mouse in the center of the ROI, clicking it and move it around. The program will allow to draw as many ROIs as inputted in the “Number of Vessels” field. 
b.	Load the ROIs from a previous analysis: ignore the “Number of Vessels” field and slick on the Import ROIs button, locate the file where the previous ROIs have been saved and select it. This will load the ROIs onto the image. The ROIs are editable and can be shrank, starched, sheared and transposed in the figure so as to correct for small movement of the image during the experiment. Once you’re done, click on “show the ROIs tc”; 
5)	Once all the ROIs are in place (either inputted ex novo or loaded from a previous experiment), the time course of the average intensity signal will appear in the central part of the panel. Use the mouse to select the point from the plot to use as initial parameters. You will need to input 4 parameters, the last 2 will be derived from the previous;
6)	All initial fitting parameters are initialized as “NaN”, click with the mouse on the point of the trace that you want to register, then click in the editbale field of the parameter and press enter. The value from the plot will be automatically acquired from the plot and displayed in the field. For example for the Amplitude and time to peak, click with the mouse on the signal peak and then in the editable field of amplitude, then press enter. Then click with the mouse in the Time to peak field and then press enter again. Do the same for Fit Start and Fit End (there are the boundaries within which the code will look for the optimal place where to start and end the fitting). For the FWHM and Baseline shift field there is no need to click on the image, just select the editable field and press enter. 
7)	Once done entering the initial parameters, click on the Fit Gamma function, this will run the fitting and overlay the fitting over the raw data. Here you have two options:
a.	You like the fitting: click on Save Value/Next trace;
b.	You dislike the fitting: you can try to change the initial parameters and then fit again until you like the fitting (and then click on Save Value/Next Trace) or just give up on this ROI, at this point click on the button NaN Initial Values and then on Save Value/Next Trace. This will result having this ROIs exported as NaN in the final report;
8)	If the time course is too noisy to be fitted (or you inputted the wrong ROI), just skip it by using the Save Value/Next Trace button;
9)	Once you’re done with all the vessels, the plotting quadrant will return white, type in the experimental condition (e.g. “air” or “co2”) and the subject number (e.g. “001”) and then click on “Export Results”. Here you can select the output folder and output file name. The output consists of a .csv file for R, a .mat file with all the useful workspace variables and a mask file (the name has something like “maskObk”) which is the one you will have to use if you want to reload the masks over a new image file. 
