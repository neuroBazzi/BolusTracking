function BolusTrack

% GUI tool to analyze changes in signal intensity 
% (e.g. fluorescence) in a time-series. Through the panel, the user can
% load the data (.tiff files), select the ROIs, input the Frame Rate, and
% estimate the initial fitting parameters using the cursor. The initial
% parameters will then be used to fit a gamma funciton through the data
% points. 

% Author: Paolo Bazzigaluppi; Last Update: January 2019.

%% Set-up main figure and buttons
f = figure('numbertitle', 'off', ...
    'name', 'Bolus Tracking Analysis', ...
    'color','w',...
    'menubar','none', ...
    'toolbar','none', ...
    'resize', 'on', ...
    'tag','main', ...
    'renderer','painters', ...
    'position',[400 200 800 600]);
     
% Import and preview control
Ip = uicontrol;
Ip.Position =  [10 550 80 25];
Ip.String = 'Load Data';
Ip.Callback = @LoadButtonPushed;
    
% Set the Frame Rate
FrText = uicontrol ('Style','text','String','Frame Rate');
FrText.Position = [10 510 80 20];
persistent Fr
Fr = uicontrol('Style','edit');  
Fr.Position = [10 485 80 25];
Fr.Callback = @FrCall;
    
        
% Vessel Number: text, editable field and send it to base WS
VnumTxt = uicontrol ('Style','text','String','Number of Vessels');
VnumTxt.Position = [10 440 80 25];
Vnum = uicontrol('Style','edit');  
Vnum.Position = [10 415 80 25];
Vnum.Callback = @VnumCall;   

% Select ROIs control
Rs = uicontrol;
Rs.Position = [10 390 80 25];
Rs.String = 'Select ROIs';
Rs.Callback = @ROIButtonPushed;
 
% Import ROIs control
Rs = uicontrol;
Rs.Position = [10 350 80 25];
Rs.String = 'Import ROIs';
Rs.Callback = @ROIImportButtonPushed;

% Import ROIs control
Rs = uicontrol;
Rs.Position = [10 325 80 25];
Rs.String = 'Show ROIs tc';
Rs.Callback = @ROItcShowButtonPushed;

% "save/next" button control (callback function at the end of the
% script)
e = uicontrol;
e.Position = [615 510 120 25];
e.String = 'Save values/next trace';
e.Callback = @SaveValueButtonPushed;

% "gamma fit" button control (callback function at the end of the
% script)
g = uicontrol;
g.String = 'Fit Gamma function';
g.Callback = @GammaFit;
g.Position = [615 550 120 25];

% Set the Experimental condition
ExpConCallText = uicontrol ('Style','text','String','Experiment:');
ExpConCallText.Position = [615 480 80 20];
ExpCon = uicontrol('Style','edit');  
ExpCon.Position = [690 480 40 25];
ExpCon.Callback = @ExpConCall;

% Set the Subject number
SubjIDCallText = uicontrol ('Style','text','String','Subject ID:');
SubjIDCallText.Position = [615 450 80 20];
SubjID = uicontrol('Style','edit');  
SubjID.Position = [690 450 40 25];
SubjID.Callback = @SubjIDCall;

% "Export" button control (assigns NaN to all the fields of both output variables)
k = uicontrol;
k.String = 'Export results';
k.Callback = @ExportButtonPushed;
k.Position = [615 415 120 25];

% Clear Initial Values
k = uicontrol;
k.String = 'NaN initial Values';
k.Callback = @ClearInitValsButtonPushed;
k.Position = [630 310 100 25];

% Clear all
k = uicontrol;
k.String = 'Clear All';
k.Callback = @ClearAllButtonPushed;
k.Position = [490 310 100 25];

% parameter 1 controls
text1 = uicontrol ('Style','text','String','Amplitude');
text1.Position = [490 550 60 20];
param1 = uicontrol('Style','edit');
param1.Position = [490 530 60 20];
param1.Callback = @param1In;    

% parameter 2 controls
text2 = uicontrol ('Style','text','String','Time to Peak');
text2.Position = [490 490 60 30];
param2 = uicontrol('Style','edit');
param2.Position = [490 470 60 20];
param2.Callback = @param2In;


% parameter 3 controls
text3 = uicontrol ('Style','text','String','Fit Start');
text3.Position = [490 430 60 20];
param3 = uicontrol('Style','edit');
param3.Position = [490 410 60 20];
param3.Callback = @param3In;

% parameter 4 controls
text4 = uicontrol ('Style','text','String','Fit End');
text4.Position = [490 370 60 20];
param4 = uicontrol('Style','edit');
param4.Position = [490 350 60 20];
param4.Callback = @param4In;

% parameter 5 controls
text5 = uicontrol ('Style','text','String','FWHM');
text5.Position = [600 370 60 20];
param5 = uicontrol('Style','edit');
param5.Position = [600 350 60 20];
param5.Callback = @param5In;

% parameter 6 controls
text6 = uicontrol ('Style','text','String','Baseline shift');
text6.Position = [670 370 60 30];
param6 = uicontrol('Style','edit');
param6.Position = [670 350 60 20];
param6.Callback = @param6In;


% preallocate plot space for time-series visualization

ax = axes(f);
ax.Units = 'pixels';
ax.Position = [160 350 250 220];
ax.XLabel.String = 'Time(s)';
ax.YLabel.String = 'Signal intensity (a.u.)';
assignin('base','ax',ax);

%% declare global variables
clear global idx
clear global fitOut
global idx
idx = 1;
global UpF
UpF = 20;
global fitOut

%% local functions 

    function ClearInitValsButtonPushed(~,~)
        param1.String = NaN;
        param2.String = NaN;
        param3.String = NaN;
        param4.String = NaN;
        param5.String = NaN;
        param6.String = NaN;
    end

    function ClearAllButtonPushed(~,~)
        
        ax = evalin ('base','ax');
        ax1 = evalin('base','ax1');
        delete(ax)
        delete(ax1)
        clear global idx
        clear global fitOut 
        evalin('base', 'clearvars *')
    end

    function ROIButtonPushed(~,~)
  
    % import necessary data from the base WS, move the t dimension in 3rd position
    VesNum = str2double(evalin('base','VesNum'));
    fileIn = evalin('base','fileIn');
    Fr = str2double(evalin('base','Fr'));
    ImageIn = permute(fileIn,[2,3,1]); % required for the dot product below
    
    
    
    % define ROIs and overlays them on the image, saves the timeseries
    maskNum = NaN(VesNum,size(fileIn,2),size(fileIn,3));
    mask = NaN(VesNum,size(fileIn,2),size(fileIn,3));
    persistent maskObj
    
    for ii = 1 : VesNum
        
        fitOut(ii).VNum = ii;
        
        % get the poligon and apexes of the ROI (use drawpolygon to create
        % and object that can be saved and edited)
        %[foo1,foo2,foo3] = roipoly;
        foo = drawpolygon;
        foo1 = poly2mask(foo.Position(:,1),foo.Position(:,2),size(fileIn,2),size(fileIn,3));
        maskObj(ii).poli = foo;
        
        % send the poligon and the coordinates of the apexes to the base WS
        % with dummy names
        assignin('base','foo1',foo1);
        assignin('base','foo2',foo.Position(:,1));
        assignin('base','foo3',foo.Position(:,2));

        % re-import them as indexed variables
        maskNum(ii,:,:) = evalin('base','foo1');
        xx(ii).c = evalin('base','foo2');
        yy(ii).c = evalin('base','foo3');
        
        % plot the ROI contour in red on the structural image using apexes 
        plot(xx(ii).c, yy(ii).c, 'r.-', 'MarkerSize', 15);
        mask(ii,:,:) = squeeze(maskNum(ii,:,:));
        maskNum(ii,:,:) = mask(ii,:,:) * ii;
        
        % creates timeseries of the signal within each ROI by averaging the 
        % point product of the data with the logical mask in every frame
        fitOut(ii).yRaw = squeeze(mean(mean(ImageIn.*squeeze(mask(ii,:,:))))); % raw timeseries
        fitOut(ii).tlRaw = linspace(0,length(fitOut(ii).yRaw)/Fr,length(fitOut(ii).yRaw)); % timeline for the raw timeseries
        
        % Cubic spline interpolation to upsample the timeseries 
        fitOut(ii).tlUs = linspace(0,length(fitOut(ii).yRaw)/Fr,length(fitOut(ii).yRaw)*UpF); % upsampled timeline
        foo1 = fitOut(ii).tlRaw';
        foo2 = fitOut(ii).yRaw;
        fitOut(ii).y = spline(foo1,foo2,fitOut(ii).tlUs); % upsampled data
        
        
    end
    
    % create a figure and plot the first trace
    plot(fitOut(1).tlUs,fitOut(1).y,'Parent',ax);
    ax.XLabel.String = 'Time(s)';
    ax.YLabel.String = 'Signal intensity (a.u.)';
    dcm_obj = datacursormode(gcf); 
    set(dcm_obj,'Enable','on'); 
    
    % send loaded data, fitOut, and ROIs to the base workspace
    assignin('base','mask',mask)
    assignin('base','maskNum',maskNum)
    assignin('base','xx',xx)
    assignin('base','maskObj',maskObj)
    
    % Set initial fitting parameters to NaN in the GUI
    param1.String = NaN;
    param2.String = NaN;
    param3.String = NaN;
    param4.String = NaN;
    param5.String = NaN;
    param6.String = NaN;  
    
end

    function ROIImportButtonPushed(~,~)
    
    % select and load in a temp file the mask from a previous experiment
    [fileNm,pathNm] = uigetfile;
    temp =  load([pathNm fileNm]);
    
    % throw the loaded masks on the image, making them editable to adjust
    % the position, size and vertexes of the ROI
    for tt = 1 : length(temp.maskObj)
        NewROI(tt) = images.roi.Polygon(gca,'Position',[temp.maskObj(tt).poli.Position]);
    end
    
    % send the new ROIs to the base WS
    assignin('base','NewROI',NewROI)
end
 
    function ROItcShowButtonPushed(~,~)
    
    % import new masks and raw data to create the ROI timecourses
    NewROI = evalin('base','NewROI');
    VesNum = length(NewROI);
    fileIn = evalin('base','fileIn');
    ImageIn = permute(fileIn,[2,3,1]); % required for the dot product below
    Fr = str2double(evalin('base','Fr'));
    
    
    % loop through the ROIs, and calculatu ROI timecourse
    for tt = 1 : VesNum
                
        % create a mask out of the new ROI vertexes
        mask(tt,:,:) = poly2mask(NewROI(tt).Position(:,1),NewROI(tt).Position(:,2), ...
            size(fileIn,2),size(fileIn,3));
        maskNum(tt,:,:) = mask(tt,:,:) * tt;
        xx(tt).c = NewROI(tt).Position(:,1);
        yy(tt).c = NewROI(tt).Position(:,2);
           
        % creates timeseries of the signal within each ROI by averaging the 
        % point product of the data with the logical mask in every frame
        fitOut(tt).yRaw = squeeze(mean(mean(ImageIn.*squeeze(mask(tt,:,:))))); % raw timeseries
        fitOut(tt).tlRaw = linspace(0,length(fitOut(tt).yRaw)/Fr,length(fitOut(tt).yRaw)); % timeline for the raw timeseries
        
        % Cubic spline interpolation to upsample the timeseries 
        fitOut(tt).tlUs = linspace(0,length(fitOut(tt).yRaw)/Fr,length(fitOut(tt).yRaw)*UpF); % upsampled timeline
        foo1 = fitOut(tt).tlRaw';
        foo2 = fitOut(tt).yRaw;
        fitOut(tt).y = spline(foo1,foo2,fitOut(tt).tlUs); % upsampled data
                
    end
    
    % create a figure and plot the first trace
    plot(fitOut(1).tlUs,fitOut(1).y,'Parent',ax);
    ax.XLabel.String = 'Time(s)';
    ax.YLabel.String = 'Signal intensity (a.u.)';
    dcm_obj = datacursormode(gcf); 
    set(dcm_obj,'Enable','on'); 
    
    % send loaded data, fitOut, and ROIs to the base workspace
    assignin('base','mask',mask)
    assignin('base','maskNum',maskNum)
    assignin('base','maskObj',NewROI)
    
    % Set initial fitting parameters to NaN in the GUI
    param1.String = NaN;
    param2.String = NaN;
    param3.String = NaN;
    param4.String = NaN;
    param5.String = NaN;
    param6.String = NaN; 
     
    
end   
      
    function SaveValueButtonPushed(~,~)
        
    % reset the values in the editable fields to zero
    param1.String = NaN;
    param2.String = NaN;
    param3.String = NaN;
    param4.String = NaN;
    param5.String = NaN;
    param6.String = NaN;  
    
    % if the counter is not the last (trace), keeps working
    if idx+1 <= length(fitOut)
        
        % spin the counter
        idx = idx + 1;
        
        % plot next trace
        plot(fitOut(idx).tlUs,fitOut(idx).y, 'Parent',ax);
        ax.XLabel.String = 'Time(s)';
        ax.YLabel.String = 'Signal intensity (a.u.)';
        dcm_obj = datacursormode(gcf); 
        set(dcm_obj,'Enable','on'); 
                        
    elseif idx+1 > length(fitOut)
        
        plot(NaN,'Parent',ax);
        assignin('base','idx',idx);
        ax.XLabel.String = 'Time(s)';
        ax.YLabel.String = 'Signal intensity (a.u.)';
        
        OnTzero = min([fitOut.OnT]); 
        for oo = 1 : length(fitOut);fitOut(oo).OnTSc = fitOut(oo).OnT - OnTzero;end
        %arrayfun(@(x) x.OnT - OnTzero, fitOut, 'UniformOutput', false);
        
        assignin('base','fitOut',fitOut);

    end

end

    function GammaFit(~,~)
    
    % import idx and output variables
    Fr = str2double(evalin('base','Fr'));
   
    % get the begin and end of the fit form param 3 and 4
    st = round(Fr*UpF*str2double(param3.String));
    en = round(Fr*UpF*str2double(param4.String));
     
    % crop the trace to fit only the bolus
    tr = fitOut(idx).y(st:en);
    
    % options
    options = statset('RobustWgtFun','cauchy');
    
    % create time line
    xxx = linspace(0,length(tr)/(Fr*UpF),length(tr));
    
    % import the initial parameters of the fit from the GUI
    paramsInit = [str2double(param1.String)...
        (str2double(param2.String)- (st/(Fr*UpF)))...
        (str2double(param5.String))...
        str2double(param6.String)];
    
    % save the initial params
    fitOut(idx).InitAmp = str2double(param1.String);
    fitOut(idx).InitT2p = str2double(param2.String)- (st/(Fr*UpF));
    fitOut(idx).InitFWHM = str2double(param5.String);
    fitOut(idx).InitM = str2double(param6.String);
    
    % outputs the parameters of non-linear fitting
    [fitOut(idx).beta,fitOut(idx).res,fitOut(idx).J,fitOut(idx).covb,fitOut(idx).mse] = ...
        nlinfit(xxx,tr,@gammaFun,paramsInit,options);
    
    % get the CI
    fitOut(idx).CI = nlparci(fitOut(idx).beta,fitOut(idx).res,'covar',fitOut(idx).covb);
    
    % fit the data and estimate AUC
    fitOut(idx).fitTr = feval(@gammaFun,fitOut(idx).beta,xxx);
    fitOut(idx).fitTrN = ((fitOut(idx).fitTr - min(fitOut(idx).fitTr)) / ...
        (max(fitOut(idx).fitTr) - min(fitOut(idx).fitTr)));
    fitOut(idx).AUC = trapz(fitOut(idx).fitTr);
    fitOut(idx).AUCn = trapz(fitOut(idx).fitTrN);
    
    % estimate onset time (OnT) and transit time (TT)
    fitOut(idx).OnT = (find((diff(find(fitOut(idx).fitTrN' < 0.1))) == 1, 1, 'last') + 1) / (Fr*UpF); % threshold: two consecutive frames smaller that 10% of Peak amplitude
    fitOut(idx).TTm = fitOut(idx).beta(2) - 1 - fitOut(idx).OnT ;
    fitOut(idx).TTlb = fitOut(idx).CI(2,1) - 1 - fitOut(idx).OnT;
    fitOut(idx).TThb = fitOut(idx).CI(2,2) - 1 - fitOut(idx).OnT;
    
    
    % build a dummy NaN variable, then overwrite only the values where
    % fitted (with the fitted values),for display purposes
    fitOut(idx).tr2pl = NaN(size(fitOut(idx).y));
    fitOut(idx).tr2pl(1,st:en) = fitOut(idx).fitTr;
    
    % export the segment of raw data that has been fitted
    fitOut(idx).rawDs = fitOut(idx).y(st:en);
    
    % plot the fit over the existing raw-data plot
    hold on
    plot(fitOut(idx).tlUs,fitOut(idx).tr2pl,'Parent',ax);
    ax.XLabel.String = 'Time(s)';
    ax.YLabel.String = 'Signal intensity (a.u.)';
    hold off
    
    % display some joyful message in the command window
    disp('fitting ok')

end

    function param1In(~,~)
    dcm_obj = datacursormode(gcf);
    c_info = getCursorInfo(dcm_obj);
    param1.String = c_info.Position(2);
    assignin('base','param1',param1.String);
end

    function param2In(~,~)
    dcm_obj = datacursormode(gcf);
    c_info = getCursorInfo(dcm_obj);
    param2.String = c_info.Position(1);
    assignin('base','param2',param2.String);
end

    function param3In(~,~)
    dcm_obj = datacursormode(gcf);
    c_info = getCursorInfo(dcm_obj);
    param3.String = c_info.Position(1);
    assignin('base','param3',param3.String); 
    assignin('base','FitStAmp',c_info.Position(2));
end

    function param4In(~,~)
    dcm_obj = datacursormode(gcf);
    c_info = getCursorInfo(dcm_obj);
    param4.String = c_info.Position(1);
    assignin('base','param4',param4.String);
    assignin('base','FitEndAmp',c_info.Position(2));
end

    function param5In(~,~)
    X1 = str2double(evalin('base','param3'));
    X2 = str2double(evalin('base','param4'));
    pk = str2double(evalin('base','param2'));
    HalfWidthU = pk-(X1+((pk - X1)/2));
    HalfWidthD = X2 -pk;
    fwhm =  (pk + HalfWidthD) - (X1 + HalfWidthU);
    param5.String = fwhm;
    assignin('base','param5',param5.String);
end

    function param6In(~,~)
    FitStAmp = (evalin('base','FitStAmp'));
    FitEndAmp = (evalin('base','FitEndAmp'));
    BslnShift = FitEndAmp - FitStAmp;
    param6.String = BslnShift;
    assignin('base','param6',param6.String);
end

    function ExportButtonPushed (~,~)
    
    % import data to export
    fitOut = evalin('base','fitOut');
    ExpCon = evalin('base','ExpCon');
    subj_num = str2double(evalin('base','subj_num'));
    mask = evalin('base','mask');
    
    % locate the vessels that do not have a fit (i.e. empty mse field) and
    % throws in NaN.
    emptyIndex = find(arrayfun(@(fitOut) isempty(fitOut.mse),fitOut));
    fName = fieldnames(fitOut);
    
    idxNaN = find(~ismember(fName,'beta'));
    idxNaN = idxNaN(idxNaN>5);
    
    for nn = 1 : length(emptyIndex)
        
        fitOut(emptyIndex(nn)).(fName{find(ismember(fName,'beta'))}) = [NaN,NaN,NaN,NaN]; %#ok<*FNDSB>
        
        for ii = 1 : length(idxNaN)
        
            fitOut(emptyIndex(nn)).(fName{idxNaN(ii)}) = NaN;

        end
    end

    % set the export location
    [TheFile,ThePath] = uiputfile('*.csv');
          
    % create a csv file and export values of the fitted vessels
    if  exist([ThePath TheFile], 'file') == 0
        
        ff = fopen([ThePath TheFile], 'w');
        
        fprintf(ff,'%s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s',...
            'subj_num,', 'ves_num,','exp,','InitAmp,','InitiT2p,','InitiFWHM,','InitM,','F_Amp,',...
            'F_T2p,','F_FWHM,','F_M,','AUC,','AUCn,','TTlb,','TTm,','TThb,','OnTSc,','ROI size,');
        
        fprintf(ff,'\n');
        
        for jj = 1 : length(fitOut)
            
            fprintf(ff,'%d,%d,%s,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%d,\n', ...
                subj_num, jj, ExpCon, fitOut(jj).InitAmp, fitOut(jj).InitT2p, fitOut(jj).InitFWHM, fitOut(jj).InitM, ...
                fitOut(jj).beta(1), fitOut(jj).beta(2), fitOut(jj).beta(3), fitOut(jj).beta(4), ...
                fitOut(jj).AUC, fitOut(jj).AUCn, fitOut(jj).TTlb, fitOut(jj).TTm, fitOut(jj).TThb, ...  
                fitOut(jj).OnTSc, length(nonzeros(mask(jj,:,:))));
            
        end
        
        fclose(ff);
    end
    
    % save all useful variables as mat files
    saveName = [TheFile(1:end-4) '.mat'];
    saveNameMaskObj = [TheFile(1:end-4) '_MaskObj.mat'];
    mask = evalin('base','mask');
    maskNum = evalin('base','maskNum');
    maskObj = evalin('base','maskObj');
    save([ThePath saveName], 'fitOut','mask','maskNum');
    save([ThePath saveNameMaskObj], 'maskObj');
     
    
end            

    function FrCall(~,~)
            assignin('base','Fr',Fr.String); % send it to base WS
        end

    function LoadButtonPushed(~,~)
    
    % opens GUI to locate file and load metadata
    [fileInName,pathIn] = uigetfile({'*.*'},[],'/media/bazzi/paolo/');
    tiff_info = imfinfo([pathIn,fileInName]); % return tiff structure, one element per image
    
    % preallocate and load tiff stack; 
    fileIn = NaN(length(tiff_info),tiff_info(1).Height, tiff_info(1).Width);
    for jj = 1 : length(tiff_info); fileIn(jj,:,:) = imread([pathIn,fileInName], jj); end   
    
    % blurr with 3x3 kernel
    krnSize= [3 3];
    for tt = 1 : size(fileIn,1); fileIn(tt,:,:) = medfilt2(squeeze(fileIn(tt,:,:)),krnSize); end
   
    % visualize the Maximum Intensity Projection (MIP)
    ax1 = subplot(2,2,3:4);
    imagesc(squeeze(max(max(fileIn,2))),'Parent',ax1)
    assignin('base','ax1',ax1);
    %imagesc(squeeze(mean(fileIn(1:10,:,:))));
    hold on
    xlabel('pixel')
    ylabel('pixel')
    colormap gray
    
    % send the input data to the base workspace
    assignin('base','fileIn',fileIn);
    
end

    function VnumCall(~,~)
        assignin('base','VesNum',Vnum.String); % send it to base WS
    end

    function ExpConCall(~,~)
        assignin('base','ExpCon',ExpCon.String);
    end

    function SubjIDCall(~,~)
        assignin('base','subj_num',SubjID.String);
    end
    
% end of the function
end 