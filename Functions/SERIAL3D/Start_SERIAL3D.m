function app = Start_SERIAL3D(app)
     % Step 0 initiate
            % get the current analysis object
            pos = app.CurrentAnalysisNum;
            obj = app.SavedAnalysis{pos,1};
            
            % load in the options for the analysis
            obj.MPTPara.f_o_s = app.SearchFieldEditField_SERIAL.Value;              % Size of search field: max(|u|,|v|)
            obj.MPTPara.n_neighborsMax = app.MaxNeighborsEditField_SERIAL.Value;     % Max # of neighboring particles
            obj.MPTPara.n_neighborsMin = app.MinNeighborsEditField_SERIAL.Value;      % Min # of neighboring particles
            obj.MPTPara.locSolver_type = app.localsolvermethodDropDown_SERIAL.Value;
            obj.MPTPara.gbSolver_type = app.GlobalstepsolverDropDown_SERIAL.Value;
            obj.MPTPara.smoothness = app.regularizationcoefEditField_SERIAL.Value;       % Coefficient of regularization
            obj.MPTPara.outlrThres = app.outlierthreshEditField_SERIAL.Value;          % Threshold for removing outliers in TPT
            obj.MPTPara.maxIterNum = app.MaxADMMiterEditField_SERIAL.Value;         % Max ADMM iteration number
            obj.MPTPara.iterStopThres = app.ADMMthreshEditField_SERIAL.Value;    % ADMM iteration stopping threshold
            obj.MPTPara.strain_n_neighbors = app.MaxstrainNeighborsEditField_SERIAL.Value; % # of neighboring particles used in strain gauge
            obj.MPTPara.strain_f_o_s = app.Maxstrain_fos_SERIAL.Value;       % Size of virtual strain gauge
            obj.MPTPara.usePrevResults = app.UsepreviousresultstodetectnewparticlesCheckBox_SERIAL.Value;      % use previous results or not: 0-no; 1-yes;
            obj.MPTPara.xstep = app.xstepEditField_SERIAL.Value;
            obj.MPTPara.tstep = app.tstepEditField_SERIAL.Value;
            obj.MPTPara.mode = app.TrackingModeDropDown_SERIAL.Value;
            obj.MPTPara.parType = app.ParticleTypeDropDown_SERIAL.Value;
            obj.MPTPara.detectionMethod = app.DetectionmethodDropDown_SERIAL.Value;
            
            % create the BeadPara structure
            obj.BeadPara.thres = app.beadthresholdEditField_SERIAL.Value ;           % Threshold for detecting particles
            obj.BeadPara.beadSize = app.avgbeadsizeEditField_SERIAL.Value;          % Estimated radius of a single particle (um)
            obj.BeadPara.minSize = app.minbeadsizeEditField_SERIAL.Value;           % Minimum radius of a single particle
            obj.BeadPara.maxSize = app.maxbeadsizeEditField_SERIAL.Value;          % Maximum radius of a single particle (px)
            obj.BeadPara.winSize = [app.maxwinsizeEditField_1_SERIAL.Value, app.maxwinsizeEditField_1_SERIAL.Value, app.maxwinsizeEditField_1_SERIAL.Value];      % By default
            obj.BeadPara.PSF_type = app.PSFtypeDropDown_SERIAL.Value;        % PSF function  
            obj.BeadPara.distMissing = app.distMissingEditField_SERIAL.Value;       % Distance threshold to check whether particle has a match or not
            obj.BeadPara.color = app.BeadColorDropDown_SERIAL.Value;       % Bead color
            
     % Step 1 load in images
             obj.Imges_raw = app.AllImages; %store the images
             obj.NumberOfImages  = length(obj.Imges_raw); %store the number of images
             obj.ImgMaskFiles = app.AllROIMasks; %store the masks
             obj.MPTPara.ImgSize = size(obj.Imges_raw{1});
             
             % Step 2 Apply ROI's
             obj = findgridxyz(obj);
             obj = ApplyROI(obj);
             
             % Step 3 Deconvolute using a PSF
             switch obj.BeadPara.PSF_type
                 case 'None'
                     obj.BeadPara.PSF = [];
                 case 'poisson'                                                                                                             %% ask JIN 
                     obj.BeadPara.PSF = [];
%                      obj.BeadPara.PSF = fspecial('disk', obj.BeadPara.beadSize-1 );
             end
             
             % Do stuff to initial image
             obj = Deconvolute(obj,1);
             obj.parCoordA = DetectParticles(obj,1);
             
     % loop through all images and do stuff
    for ImgSeqNum = 2 : obj.NumberOfImages
         
        % Step 3 Deconvolute using a PSF
             obj = Deconvolute(obj,ImgSeqNum);

%          % Step 4 If beads are black then adjust the image                 % not implemented yet
%                 if obj.BeadPara.color == 'black'
%                     obj = BlackParticles(obj,ImgSeqNum);
%                 end
                
         % Step 5 & 6 detect the particles
                [parCoord, obj] = DetectParticles(obj,ImgSeqNum);
                obj.parCoordB = parCoord;
         % Step 8 track the particle locations
                obj = runSERIAL3D(obj,ImgSeqNum);
                
                %  plot trajectories
                obj = plotTrajectories(obj,ImgSeqNum);
                
         % figure out something to plot for each sequence num maybe the
         % particle locations?
         coneplot = obj.im_cone{ImgSeqNum};
         beadtrack = obj.im_beadtrack{ImgSeqNum};
         
         func_PlaceFigure(app.A_disp_SERIAL,coneplot);
         func_PlaceFigure(app.A_TP_SERIAL,beadtrack);
    end

    % Step 9 implement tracking type (cum, inc, etc)
          obj = CumulativeTrack(obj);

    % Step 10 Compute images
          obj = ComputeTrajectory(obj);
          
          if license('test','PDE_Toolbox')
            % Step 11 plot everything
            obj = plotDispAndStrainMaps(obj);
          else
              disp('sorry we cant display the strains at this time')
          end
    
    % Step 11 Save the object 
     app.SavedAnalysis{pos,1} = obj;

end