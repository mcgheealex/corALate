classdef SERIAL3D
    % ---------------------------------------------
    % Augmented Lagrangian Digital Image Correlation (SERIAL-TRACK)
    %
    % Author of algorithms: Jin Yang, PhD @Caltech
    % Author of GUI and adapted UI: Alexander McGhee, PhD @UF
    % Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
    % Date: 2015.04,06,07; 2016.03,04; 2020.11
    % ---------------------------------------------
    % The following functions are modified/adapted from the original
    % ReadImage -> moved to the Functions/ALDIC folder
    % IntegerSearch -> InitSearchSize
    
    % %%%%%%%%%%%%%%%%%% SerialTrack (2D cumulative mode) %%%%%%%%%%%%%%%%%
    % Main file of code "SerialTrack"
    % ***********************************************
    %
    % -----------------------------------------------
    % References
    % [1] M Patel, SE Leggett, AK Landauer, IY Wong, C Franck. Rapid,
    %     topology-based particle tracking for high-resolution measurements of
    %     large complex 3D motion fields. Scientific Reports. 8:5581 (2018).
    % [2] J Yang, L Hazlett, AK Landauer, C Franck. Augmented Lagrangian
    %     Digital Volume Correlation (ALDVC). Experimental Mechanics (2020).
    % [3] T Janke, R Schwarze, K Bauer. Part2Track: A MATLAB package for double
    %     frame and time resolved Particle Tracking Velocimetry. 11, 100413, SoftwareX (2020).
    % [4] J Heyman. TracTrac: a fast multi-object tracking algorithm for motion
    %     estimation. Computers & Geosciences, vol 128, 11-18 (2019).
    % [5] https://www.mathworks.com/matlabcentral/fileexchange/77347-gridded-interpolation-and-gradients-of-3d-scattered-data
    % [6] https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend
    % -----------------------------------------------
    % Author: Jin Yang
    % Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
    % Date: 2020.12.
    % Author of GUI and adapted UI: Alexander McGhee, PhD
    % Contact and support: amcghee2@wisc.edu
    % Date: 05/2022
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        Type = 'SERIAL3D';
        ProgramFolder = 'SerialTrack3D';
        Name = '';
        SaveFolderPath = '';
        
        % parameters
        MPTPara  = [];
        BeadPara = [];
        
        % images used for analysis
        Imges_raw = [];
        Imges_used = [];
        NumberOfImages = [];
        ImgMaskFiles = [];
        
        % coordinates of the beads
        parCoordB = [];
        parCoordA = [];
        
        % results
        parCoord_prev =[];
        uvw_B2A_prev = [];  % cumulative displacement
        track_A2B_prev =[];
        resultDisp = [];
        resultDefGrad= [];
        parCoordATraj = [];
        
        % cum track
        track_ratio = [];
        DefType = [];
        defList = [];
        track_A2B = [];
        
        % save pictures of the figure for display in the GUI
        im_cone = [];
        im_disp = [];
        im_strain =[];
        im_beadtrack = [];
        im_traj = [];
        
        ResultDisp={};
        ResultDefGrad = {};
        ResultFEMeshEachFrame_coordinatesFEM = {};
        ResultFEMeshEachFrame_elementsFEM = {};
        
        % result variables for plotting
        elementsFEM = [];
        coordinatesFEM = [];
        U = [];
        F = [];
        
    end
    
    methods
        % these are the main functions to use the ALDIC software
        function obj = SERIAL3D()
            %SERIAL3D Constructor
            % add the SERIAL3D functions to the path and remove others
            AddAndRemovePathFolders(obj,1)
            
            % create the MPTPara structure
            MPTPara.ImgRefMask = [];
            MPTPara.gridxyzROIRange.gridy = [];
            MPTPara.gridxyzROIRange.gridx = [];
            MPTPara.gridxyzROIRange.gridz = [];
            MPTPara.LoadImgMethod = 0;
            MPTPara.ImgSize = [];
            MPTPara.f_o_s = 60;              % Size of search field: max(|u|,|v|,|w|)
            MPTPara.n_neighborsMax = 25;     % Max # of neighboring particles
            MPTPara.n_neighborsMin = 1;      % Min # of neighboring particles
            MPTPara.locSolver = 1;           % Local solver: 1-topology-based feature; 2-histogram-based feature first and then topology-based feature;
            MPTPara.gbSolver = 2;            % Global step solver: 1-moving least square fitting; 2-global regularization; 3-ADMM iterations
            MPTPara.smoothness = 1e-2;       % Coefficient of regularization
            MPTPara.outlrThres = 5;          % Threshold for removing outliers in TPT
            MPTPara.maxIterNum = 20;         % Max ADMM iteration number
            MPTPara.iterStopThres = 1e-2;    % ADMM iteration stopping threshold
            MPTPara.strain_n_neighbors = 20; % # of neighboring particles used in strain gauge
            MPTPara.strain_f_o_s = 60;       % Size of virtual strain gauge
            MPTPara.usePrevResults = 0;      % Whether use previous results or not: 0-no; 1-yes;
            MPTPara.xstep = 1;
            MPTPara.tstep = 1;
            MPTPara.DIM = 2;
            MPTPara.mode = 'cum';
            MPTPara.parType = 'hard';
            MPTPara.detectionMethod = 'LoG';
            
            % create the BeadPara structure
            BeadPara.thres = 0.5 ;           % Threshold for detecting particles
            BeadPara.beadSize = 3;          % Estimated radius of a single particle (um)
            BeadPara.minSize = 2;           % Minimum radius of a single particle
            BeadPara.maxSize = 20;          % Maximum radius of a single particle (px)
            BeadPara.winSize = [5, 5, 5];      % By default
            BeadPara.dccd = [1,1,1];          % By default
            BeadPara.abc = [1,1,1];           % By default
            BeadPara.forloop = 1;           % By default
            BeadPara.randNoise = 1e-7;      % By default
            BeadPara.PSF = [];              % PSF function; Example: PSF = fspecial('disk', BeadPara.beadSize-1 ); % Disk blur
            BeadPara.distMissing = 5;       % Distance threshold to check whether particle has a match or not
            BeadPara.color = 'white';       % Bead color
            
            obj.MPTPara = MPTPara;
            obj.BeadPara = BeadPara;
            
        end
        
        function obj = runSERIAL3D(obj,ImgSeqNum)
            
            % SerialTrack particle tracking
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %  Particle Linking
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ====== Track A & B neighbor-topology match ======
            %
            %  Directly run this section if coordinates of detected particles are known:
            %  Coordinates of particles in reference image:   parCoordA
            %  Coordinates of particles in deformed image:    parCoordB
            %
            %    \  |  /                  \  |  /
            %     \ | /                    \ | /
            %   --- A ---       v.s.     --- B ---
            %      / \                      / \
            %     /   \                    /   \
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            [matches_A2B,u_B2A_curr_refB,obj.track_A2B] = f_track_serial_match3D( obj.parCoordA, obj.parCoordB, ...
                'f_o_s',obj.MPTPara.f_o_s, 'n_neighborsMax', obj.MPTPara.n_neighborsMax, 'n_neighborsMin', obj.MPTPara.n_neighborsMin, ...
                'gbSolver', obj.MPTPara.gbSolver, 'smoothness', obj.MPTPara.smoothness, ...
                'outlrThres', obj.MPTPara.outlrThres, 'maxIterNum', obj.MPTPara.maxIterNum, ...
                'iterStopThres', obj.MPTPara.iterStopThres, 'usePrevResults', obj.MPTPara.usePrevResults, ...
                'strain_n_neighbors',obj.MPTPara.strain_n_neighbors, 'strain_f_o_s',obj.MPTPara.strain_f_o_s, ...
                'gridxyzROIRange',obj.MPTPara.gridxyzROIRange, 'parCoordB_prev',obj.parCoord_prev(2:end), ...
                'uvw_B2A_prev',obj.uvw_B2A_prev, 'ImgSeqNum',ImgSeqNum, ...
                'BeadParaDistMissing',obj.BeadPara.distMissing);
            
            %%%%% Compute track_B2A %%%%%
            matches_A2B = matches_A2B(matches_A2B(:,2)>0,:); % Build matches_A2B_temp
            track_B2A = zeros(size(obj.parCoordB,1), 1);
            for tempi = 1:size(matches_A2B)
                track_B2A(matches_A2B(tempi,2)) = matches_A2B(tempi,1);
            end
            
            
            %% %%%%% Plotting %%%%%%
            % Compute displacement from tracked particles on deformed frame
            disp_A2B_parCoordB = -u_B2A_curr_refB;
            
            h = figure('visible','off');
            
            plotCone3(obj.parCoordB(:,1),obj.parCoordB(:,2),obj.parCoordB(:,3),disp_A2B_parCoordB(:,1),disp_A2B_parCoordB(:,2),disp_A2B_parCoordB(:,3));
            set(gca,'fontsize',18); box on; axis equal; axis tight; view(3); set(gca,'YDir','reverse');
            title('Tracked displacements','fontweight','normal');
            xlabel(''); ylabel(''); cb = colorbar; set(cb,'fontsize',18);
            
            frame = getframe(h);
            obj.im_cone{ImgSeqNum} = frame2im(frame);
            
            
            %% %%%%% Compute F = def grad = grad_u = (grad_x X - I) %%%%%
            
            % strain_f_o_s: size of virtual strain gauge
            % strain_n_neighbors: # of neighboring particles used in strain gauge
            
            %%%%% Strain based on Moving Least Square Fitting in deformed configuration %%%%%
            [XYZ_refB,U_B2A_refB,F_B2A_refB] = funCompDefGrad3(disp_A2B_parCoordB, obj.parCoordB, obj.MPTPara.strain_f_o_s, obj.MPTPara.strain_n_neighbors);
            
            %%%%% Strain based on Moving Least Square Fitting in reference configuration %%%%%
            [XYZ_refA,U_A2B_refA,F_A2B_refA] = funCompDefGrad3(disp_A2B_parCoordB, obj.parCoordB-disp_A2B_parCoordB, obj.MPTPara.f_o_s, obj.MPTPara.strain_n_neighbors);
            
            %% %%%%% Store results %%%%%
            uvw_B2A_refB = u_B2A_curr_refB;
            
            resultDisp_tmp = struct('parCoordA',obj.parCoordA,'parCoordB',obj.parCoordB,'track_A2B', ...
                obj.track_A2B,'disp_A2B_parCoordB',disp_A2B_parCoordB);
            
            resultDefGrad_tmp = struct('XY_refA',XYZ_refA,'U_A2B_refA',U_A2B_refA,'F_A2B_refA',F_A2B_refA, ...
                'XY_refB',XYZ_refB,'U_B2A_refB',U_B2A_refB,'F_B2A_refB',F_B2A_refB);
            
            % store results
            obj.parCoord_prev{ImgSeqNum} = obj.parCoordB;
            obj.uvw_B2A_prev{ImgSeqNum-1} = uvw_B2A_refB;  % cumulative displacement
            obj.track_A2B_prev{ImgSeqNum-1} = obj.track_A2B;
            obj.resultDisp{ImgSeqNum-1} = resultDisp_tmp;
            obj.resultDefGrad{ImgSeqNum-1} = resultDefGrad_tmp;
            
        end
        
        function [parCoord,obj] = DetectParticles(obj,ImgSeqNum)
            
            currImg = obj.Imges_used{ImgSeqNum};
            method = obj.MPTPara.detectionMethod;
            
            %%%%% Several methods to detect particles %%%%%
            switch method
                %%%%% Method 1: TPT code %%%%%
                case 'TPT'
                    x{1}{ImgSeqNum} = locateParticles(double(currImg)/max(double(currImg(:))),obj.BeadPara); % Detect particles
                    x{1}{ImgSeqNum} = radialcenter3dvec(double(currImg),x{1}{ImgSeqNum},obj.BeadPara); % Localize particles
                    
                    %%%%% Method 2: LoG operator (modified TracTrac code) %%%%%
                case 'LoG'
                    x{1}{ImgSeqNum} = f_detect_particles3(double(currImg)/max(double(currImg(:))),obj.BeadPara);
                    
                    %%%%% Method 3: sift code %%%%%
                case  'sift' % % Not fully implemented yet
                    disp('not implemented')
                    
            end
            
            % Add MPTPara.gridxyzROIRange left-bottom corner point coordinates
            x{1}{ImgSeqNum} = x{1}{ImgSeqNum} + [obj.MPTPara.gridxyzROIRange.gridx(1)-1, ...
                obj.MPTPara.gridxyzROIRange.gridy(1)-1, ...
                obj.MPTPara.gridxyzROIRange.gridz(1)-1];
            parCoord = x{1}{ImgSeqNum};
            
            % Remove bad coordinates that are out of image ROI
            for tempi=1:3, parCoord( parCoord(:,tempi)>size(currImg,tempi), : ) = []; end
            for tempi=1:3, parCoord( parCoord(:,tempi)<1, : ) = []; end
            
            % make a figure and plot the background image
            fig1=figure('visible','off');
            plot3(parCoord(:,1),parCoord(:,2),parCoord(:,3),'ro');
            view(3); box on; axis equal; axis tight; set(gca,'fontsize',18);
            title('Detected particles in image','fontweight','normal');
            frame = getframe(fig1);
            obj.im_beadtrack{ImgSeqNum} = frame2im(frame);
            
        end
        
        function obj = CumulativeTrack(obj)
            %%%%% Cumulative tracking ratio %%%%%
            obj.track_ratio = zeros(obj.NumberOfImages-1,1);
            obj.DefType = 'exp';
            obj.defList = [2:1:obj.NumberOfImages]';
            
            for ImgSeqNum = 2 : obj.NumberOfImages
                obj.track_A2B = obj.track_A2B_prev{ImgSeqNum-1};
                obj.track_ratio(ImgSeqNum-1) = length(obj.track_A2B(obj.track_A2B>0))/size(obj.parCoord_prev{ImgSeqNum},1);
            end
            
        end
        
        function obj = BlackParticles(obj)
            disp('not implemented yet')
            % this function takes in an image with white background and
            % black beads and saves that image into obj.Imges_used for
            % future use
            %             currImg = obj.Imges_used{ImgSeqNum};
            %             ImgGauss = imgaussfilt(imgaussfilt(currImg,1),1); % figure, imshow(uint16(ImgGauss));
            %             ImgGauss(ImgGauss > obj.BeadPara.thres*max(double(currImg(:)))) = 0;
            %             bw = imbinarize(uint16(ImgGauss),'adaptive','ForegroundPolarity','dark','Sensitivity',0.8); % figure, imshow(bws2);
            %             bws2 = bwareaopen(bw,round(pi*obj.BeadPara.minSize^2)); % remove all object containing fewer than BeadPara.minSize
            %             removeobjradius = obj.BeadPara.minSize; % fill a gap in the pen's cap
            %             se = strel('disk',removeobjradius);
            %             bws2 = imclose(bws2,se);
            %             obj.Imges_used{ImgSeqNum} = double(bws2); % figure, imshow(uint8(currImg2));
        end
        
        function obj = ApplyROI(obj)
            obj.MPTPara.ImgRefMask = obj.ImgMaskFiles{1};
            for ImgSeqNum = 1 : obj.NumberOfImages
                % apply a mask to the image to ensure no particles in these
                % areas are detected
                Img =obj.Imges_raw{ImgSeqNum}.*obj.ImgMaskFiles{ImgSeqNum};
                
                obj.Imges_used{ImgSeqNum} = Img(obj.MPTPara.gridxyzROIRange.gridx(1):obj.MPTPara.gridxyzROIRange.gridx(2), ...
                    obj.MPTPara.gridxyzROIRange.gridy(1):obj.MPTPara.gridxyzROIRange.gridy(2), ...
                    obj.MPTPara.gridxyzROIRange.gridz(1):obj.MPTPara.gridxyzROIRange.gridz(2));
            end
            
        end
        
        function obj = Deconvolute(obj,ImgSeqNum)
            % modify the image to deconvolute the beads
            
            currImg = obj.Imges_used{ImgSeqNum};
            if ~isempty(obj.BeadPara.PSF) % If PSF is non-empty, SerialTrack performs deconvolution %
                currImg = deconvlucy(currImg,obj.BeadPara.PSF,6);
            end
            % save the image
            obj.Imges_used{ImgSeqNum} = currImg;
            
        end
        
        function obj = ComputeTrajectory(obj)
            
            %%%%% Initialization %%%%%
            resultDispCurr = obj.resultDisp{1};
            obj.parCoordA = resultDispCurr.parCoordA;
            obj.parCoordATraj = cell(size(obj.parCoordA,1),1);
            
            %%%%% Compute and collect all trajectory segments %%%%%
            for parInd = 1:size(obj.parCoordA,1)
                
                for ImgSeqNum = 2:(size(obj.resultDisp,1)+1)
                    
                    resultDispCurr = obj.resultDisp{ImgSeqNum-1};
                    obj.parCoordB = resultDispCurr.parCoordB;
                    obj.track_A2B = resultDispCurr.track_A2B;
                    
                    if obj.track_A2B(parInd) > 0
                        obj.parCoordATraj{parInd}(ImgSeqNum-1,1:3) = obj.parCoordB(obj.track_A2B(parInd),1:3);
                    else
                        obj.parCoordATraj{parInd}(ImgSeqNum-1,1:3) = [nan,nan,nan];
                    end
                end
                
            end
            
        end
        
        function AddAndRemovePathFolders(obj,stepNum)
            
            % folders to add to path
            Fa = obj.ProgramFolder;
            
            % get the path info for this class script.
            [ThisPath,~,~] = fileparts(mfilename('fullpath'));
            CodePath = ThisPath(1:end-3);
            
            % remove all folders from the programs folder
            warning('off','all')
            rmpath(genpath([CodePath,'corALate Programs\']))
            warning('on','all')
            
            if isequal(stepNum, 1)
                % add ALDIC folders to path
                addpath(genpath([CodePath,'corALate Programs\',Fa]));
            end
        end
        
        %%%%%  PLOT Functions  %%%%%%%
        
        function obj = plotTrajectories(obj,ImgSeqNum)
            xstep = obj.MPTPara.xstep;
            tstep = obj.MPTPara.tstep;
            %%%%% Plot tracked trajectories %%%%%
            
            h = figure('visible','off');% create a hidden figure
            for parInd = 1:size(obj.parCoordA,1)
                try
                    wayPoints = obj.parCoordATraj{parInd};
                    if (size(obj.resultDisp,1)+1)<4
                        hold on;
                        line(wayPoints(isnan(wayPoints(:,1))<1,1),wayPoints(isnan(wayPoints(:,1))<1,2),wayPoints(isnan(wayPoints(:,1))<1,3));
                        view(3); % straight lines
                    else
                        hold on;
                        fnplt(cscvn(wayPoints(isnan(wayPoints(:,1))<1,:)'),'',1);
                    end
                catch
                end
            end
            
            set(gca,'fontsize',18); view(3); box on; axis equal; axis tight;
            title('Tracked particle trajectory','fontweight','normal');
            xlabel('x'); ylabel('y'); zlabel('z');
            axis(xstep*[obj.MPTPara.gridxyzROIRange.gridx(1), obj.MPTPara.gridxyzROIRange.gridx(2), ...
                obj.MPTPara.gridxyzROIRange.gridy(1), obj.MPTPara.gridxyzROIRange.gridy(2), ...
                obj.MPTPara.gridxyzROIRange.gridz(1), obj.MPTPara.gridxyzROIRange.gridz(2)]);
            frame = getframe(h);
            obj.im_traj{ImgSeqNum} = frame2im(frame);
            
        end
        
        function obj = plotDispAndStrainMaps(obj)
            
            % make this a for loop for all images
            xstep = obj.MPTPara.xstep;
                tstep = obj.MPTPara.tstep;
            for ImgSeqNum = 2:obj.NumberOfImages
                
                %%%%% Previously tracked displacement field %%%%%
                resultDispCurr = obj.resultDisp{ImgSeqNum-1};
                resultDefGradCurr = obj.resultDefGrad{ImgSeqNum-1};
                disp_A2B_parCoordB = resultDispCurr.disp_A2B_parCoordB;
                obj.parCoordB = resultDispCurr.parCoordB;
                
                %%%%% Shift rigid body translations %%%%%
                disp_A2B_parCoordB(:,1) = disp_A2B_parCoordB(:,1) - median(disp_A2B_parCoordB(:,1));
                disp_A2B_parCoordB(:,2) = disp_A2B_parCoordB(:,2) - median(disp_A2B_parCoordB(:,2));
                disp_A2B_parCoordB(:,3) = disp_A2B_parCoordB(:,3) - median(disp_A2B_parCoordB(:,3));
                
                %%%%% Interpolate scatterred data to gridded data %%%%%
                sxyz = min([round(0.5*obj.MPTPara.f_o_s),20])*[1,1,1]; % Step size for griddata
                smoothness = 1e-3; % Smoothness for regularization; "smoothness=0" means no regularization
                
                [x_Grid_refB,y_Grid_refB,z_Grid_refB,u_Grid_refB]=funScatter2Grid3D(obj.parCoordB(:,1),obj.parCoordB(:,2),obj.parCoordB(:,3),disp_A2B_parCoordB(:,1),sxyz,smoothness);
                [~,~,~,v_Grid_refB]=funScatter2Grid3D(obj.parCoordB(:,1),obj.parCoordB(:,2),obj.parCoordB(:,3),disp_A2B_parCoordB(:,2),sxyz,smoothness);
                [~,~,~,w_Grid_refB]=funScatter2Grid3D(obj.parCoordB(:,1),obj.parCoordB(:,2),obj.parCoordB(:,3),disp_A2B_parCoordB(:,3),sxyz,smoothness);
                
                obj.MPTPara.ImgRefMask = obj.ImgMaskFiles{1};
                
                % Build a displacement vector
                uvw_Grid_refB_Vector=[u_Grid_refB(:),v_Grid_refB(:),w_Grid_refB(:)]';
                uvw_Grid_refB_Vector=uvw_Grid_refB_Vector(:);
                
                % Calculate deformation gradient
                D_Grid = funDerivativeOp3(size(x_Grid_refB,1),size(x_Grid_refB,2),size(x_Grid_refB,3),sxyz); % Central finite difference operator
                F_Grid_refB_Vector=D_Grid*uvw_Grid_refB_Vector; % {F}={D}{U}
                
                %%%%% Generate an FE-mesh %%%%%
                [coordinatesFEM_refB,elementsFEM_refB] = funMeshSetUp3(x_Grid_refB*xstep,y_Grid_refB*xstep,z_Grid_refB*xstep);
                
                obj.elementsFEM{ImgSeqNum} = elementsFEM_refB;
                obj.coordinatesFEM{ImgSeqNum} = coordinatesFEM_refB*xstep;
                obj.U{ImgSeqNum} = uvw_Grid_refB_Vector*xstep;
                obj.F{ImgSeqNum} = F_Grid_refB_Vector;
                
                %%%%% Cone plot grid data: displecement %%%%%
                h = figure('visible','off');% create a hidden figure
                plotCone3(x_Grid_refB*xstep,y_Grid_refB*xstep,z_Grid_refB*xstep,u_Grid_refB*xstep,v_Grid_refB*xstep,w_Grid_refB*xstep);
                set(gca,'fontsize',18); view(3); box on; axis equal; axis tight;
                title('Tracked incremental displacement','fontweight','normal');
                axis(xstep*[obj.MPTPara.gridxyzROIRange.gridx(1), obj.MPTPara.gridxyzROIRange.gridx(2), ...
                    obj.MPTPara.gridxyzROIRange.gridy(1), obj.MPTPara.gridxyzROIRange.gridy(2), ...
                    obj.MPTPara.gridxyzROIRange.gridz(1), obj.MPTPara.gridxyzROIRange.gridz(2)]);
                frame = getframe(h);
                obj.im_cone{ImgSeqNum} = frame2im(frame);
                
                %%%%% Cone plot grid data: displacement %%%%%
                [obj.im_disp.u{ImgSeqNum},obj.im_disp.v{ImgSeqNum},obj.im_disp.w{ImgSeqNum}] = obj.Plotdisp_show_obj(obj,obj.U{ImgSeqNum}, obj.coordinatesFEM{ImgSeqNum}, obj.elementsFEM{ImgSeqNum},[],'NoEdgeColor');
                
                %%%%% Cone plot grid data: infinitesimal strain %%%%%
                [obj.im_strain.exx{ImgSeqNum}, obj.im_strain.eyy{ImgSeqNum}, obj.im_strain.ezz{ImgSeqNum},obj.im_strain.exy{ImgSeqNum}, obj.im_strain.exz{ImgSeqNum}, obj.im_strain.eyz{ImgSeqNum}] = Plotstrain_show_obj(obj,obj.F{ImgSeqNum}, obj.coordinatesFEM{ImgSeqNum} , obj.elementsFEM{ImgSeqNum},[],'NoEdgeColor',xstep,tstep);
            end
            obj.ResultDisp = obj.U;
            obj.ResultDefGrad = obj.F;
            obj.ResultFEMeshEachFrame_coordinatesFEM = obj.coordinatesFEM;
            obj.ResultFEMeshEachFrame_elementsFEM = obj.elementsFEM;
            close all
        end
        
        function obj = findgridxyz(obj)
            
            for ImgSeqNum = 1 : obj.NumberOfImages
                ROI = obj.ImgMaskFiles{ImgSeqNum};
                bounds = regionprops(ROI,'BoundingBox');
                if length(bounds)== 1
                    b = bounds.BoundingBox;
                    gridx = round([b(1), b(1)+b(4)-1]);
                    gridy = round([b(2), b(2)+b(5)-1]);
                    gridz = round([b(3), b(3)+b(6)-1]);
                end
                
                imgszx = obj.MPTPara.ImgSize(1);
                imgszy = obj.MPTPara.ImgSize(2);
                imgszz = obj.MPTPara.ImgSize(3);
                
                %%%%% Modify gridxy to make sure it is within the image domain %%%%%
                if gridx(1)<1, gridx(1) = 1; end
                if gridy(1)<1, gridy(1) = 1; end
                if gridz(1)<1, gridz(1) = 1; end
                if gridx(2)>imgszx, gridx(2) = imgszx; end
                if gridy(2)>imgszy, gridy(2) = imgszy; end
                if gridz(2)>imgszz, gridz(2) = imgszz; end
                
                % Store TPTpara
                obj.MPTPara.gridxyzROIRange.gridy = gridy;
                obj.MPTPara.gridxyzROIRange.gridx = gridx;
                obj.MPTPara.gridxyzROIRange.gridz = gridz;
            end
            
        end
        
        function  [im_exx,im_eyy,im_ezz,im_exy,im_exz,im_eyz] = Plotstrain_show_obj( obj, F,coordinatesFEM,varargin)
            %FUNCTION Plotstrain_show(F,coordinatesFEM,elementsFEM)
            % To plot DVC solved 3D strain components
            % ----------------------------------------------
            %
            %   INPUT: F                 Deformation gradient tensor:
            %                            F = [F11_node1, F21_node1, F31_node1, ...
            %                                 F12_node1, F22_node1, F32_node1, ...
            %                                 F13_node1, F23_node1, F33_node1, ...
            %                                 ... , ...
            %                                 F11_nodeN, F21_nodeN, F31_nodeN, ...
            %                                 F12_nodeN, F22_nodeN, F32_nodeN, ...
            %                                 F13_nodeN, F23_nodeN, F33_nodeN]';
            %          coordinatesFEM    FE mesh coordinates
            %          DICpara           chosen DVC parameters
            %          EdgeColorOrNot    show edge color or not
            
            %   OUTPUT: Plots of exx, eyy, and exy strain fields.
            %
            % ----------------------------------------------
            % Author: Jin Yang.
            % Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
            % Last date modified: 2020.12
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Generate model
            dvcVOI = createpde(1);
            
            % Apply mesh
            DT = delaunayTriangulation(coordinatesFEM(:,1), coordinatesFEM(:,2), coordinatesFEM(:,3));
            geometryFromMesh(dvcVOI,DT.Points',DT.ConnectivityList');
            % ------ FEMesh has structure ------
            % FEMesh with properties:
            %
            %              Nodes: [3x10003 double]
            %           Elements: [10x5774 double]
            %     MaxElementSize: 9.7980
            %     MinElementSize: 4.8990
            %      MeshGradation: 1.5000
            %     GeometricOrder: 'quadratic'
            % ----------------------------------
            
            % plot dvcZOI domain
            % ---------------------------------
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ====== 1) strain exx ======
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h = figure('visible','off');% create a hidden figure
            pdeplot3D(dvcVOI,'ColorMapData',F(1:9:end),'FaceAlpha',0.5);
            
            title('Strain $e_{11}$','fontweight','normal','Interpreter','latex');
            set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
            xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            a = gca; a.TickLabelInterpreter = 'latex';
            b = colorbar; b.TickLabelInterpreter = 'latex';
            set(gcf,'color','w'); colormap jet; colorbar; box on
            
            frame = getframe(h);
            im_exx = frame2im(frame);
            
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ====== 2) strain exy ======
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h = figure('visible','off');% create a hidden figure
            pdeplot3D(dvcVOI,'ColorMapData',F(5:9:end),'FaceAlpha',0.5);
            
            title('Strain $e_{22}$','fontweight','normal','Interpreter','latex');
            set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
            xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            a = gca; a.TickLabelInterpreter = 'latex';
            b = colorbar; b.TickLabelInterpreter = 'latex';
            set(gcf,'color','w'); colormap jet; colorbar; box on
            
            frame = getframe(h);
            im_eyy = frame2im(frame);
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ====== 3) strain ezz ======
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h = figure('visible','off');% create a hidden figure
            pdeplot3D(dvcVOI,'ColorMapData',F(9:9:end),'FaceAlpha',0.5);
            
            title('Strain $e_{33}$','fontweight','normal','Interpreter','latex');
            set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
            xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            a = gca; a.TickLabelInterpreter = 'latex';
            b = colorbar; b.TickLabelInterpreter = 'latex';
            set(gcf,'color','w'); colormap jet; colorbar; box on
            
            frame = getframe(h);
            im_ezz = frame2im(frame);
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ====== 3) strain exy ======
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h = figure('visible','off');% create a hidden figure
            pdeplot3D(dvcVOI,'ColorMapData',0.5*(F(2:9:end)+F(4:9:end)),'FaceAlpha',0.5);
            
            title('Strain $e_{12}$','fontweight','normal','Interpreter','latex');
            set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
            xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            a = gca; a.TickLabelInterpreter = 'latex';
            b = colorbar; b.TickLabelInterpreter = 'latex';
            set(gcf,'color','w'); colormap jet; colorbar; box on
            
            frame = getframe(h);
            im_exy = frame2im(frame);
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ====== 3) strain exz ======
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h = figure('visible','off');% create a hidden figure
            pdeplot3D(dvcVOI,'ColorMapData',0.5*(F(3:9:end)+F(7:9:end)),'FaceAlpha',0.5);
            
            title('Strain $e_{13}$','fontweight','normal','Interpreter','latex');
            set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
            xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            a = gca; a.TickLabelInterpreter = 'latex';
            b = colorbar; b.TickLabelInterpreter = 'latex';
            set(gcf,'color','w'); colormap jet; colorbar; box on
            
            frame = getframe(h);
            im_exz = frame2im(frame);
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ====== 3) strain eyz ======
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h = figure('visible','off');% create a hidden figure
            pdeplot3D(dvcVOI,'ColorMapData',0.5*(F(6:9:end)+F(8:9:end)),'FaceAlpha',0.5);
            
            title('Strain $e_{23}$','fontweight','normal','Interpreter','latex');
            set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
            xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            a = gca; a.TickLabelInterpreter = 'latex';
            b = colorbar; b.TickLabelInterpreter = 'latex';
            set(gcf,'color','w'); colormap jet; colorbar; box on
            
            frame = getframe(h);
            im_eyz = frame2im(frame);
            
        end
        
        function  [im_u, im_v, im_w] = Plotdisp_show_obj(obj, U,coordinatesFEM,varargin)
            %PLOTDISP_SHOW: to plot DIC solved displacement components
            %   Plotdisp_show(U,coordinatesFEM,elementsFEM)
            % ----------------------------------------------
            %
            %   INPUT: U                 Displacement vector:
            %                            U = [Ux_node1, Uy_node1, Ux_node2, Uy_node2, ... , Ux_nodeN, Uy_nodeN]';
            %          coordinatesFEM    FE mesh coordinates
            %          elementsFEM       FE mesh elements
            %          DICpara           chosen DIC parameters
            %          EdgeColorOrNot    show edge color or not
            %
            %   OUTPUT: Plots of x-displacement field and y-displacement field.
            %
            %   TODO: users could change caxis range based on their own choices.
            %
            % ----------------------------------------------
            % Author: Jin Yang.
            % Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
            % Last date modified: 2020.12
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %% Generate model
            dvcZOI = createpde(1);
            
            % Apply mesh
            DT = delaunayTriangulation(coordinatesFEM(:,1), coordinatesFEM(:,2), coordinatesFEM(:,3));
            geometryFromMesh(dvcZOI,DT.Points',DT.ConnectivityList');
            % ------ FEMesh has structure ------
            % FEMesh with properties:
            %
            %              Nodes: [3x10003 double]
            %           Elements: [10x5774 double]
            %     MaxElementSize: 9.7980
            %     MinElementSize: 4.8990
            %      MeshGradation: 1.5000
            %     GeometricOrder: 'quadratic'
            % ----------------------------------
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ====== 1) dispx u ======
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h = figure('visible','off');% create a hidden figure
            pdeplot3D(dvcZOI,'ColorMapData',U(1:3:end),'FaceAlpha',0.5);
            
            title('$x$-displacement $u$','FontWeight','Normal','Interpreter','latex');
            set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
            xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            
            a = gca; a.TickLabelInterpreter = 'latex';
            b = colorbar; b.TickLabelInterpreter = 'latex';
            
            set(gcf,'color','w'); colormap jet; colorbar; box on
            
            frame = getframe(h);
            im_u = frame2im(frame);
            
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ====== 2) dispx v ======
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h = figure('visible','off');% create a hidden figure
            pdeplot3D(dvcZOI,'ColorMapData',U(2:3:end),'FaceAlpha',0.5);
            
            title('$y$-displacement $v$','FontWeight','Normal','Interpreter','latex');
            set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
            xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            
            a = gca; a.TickLabelInterpreter = 'latex';
            b = colorbar; b.TickLabelInterpreter = 'latex';
            
            set(gcf,'color','w'); colormap jet; colorbar; box on
            
            frame = getframe(h);
            im_v = frame2im(frame);
            
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ====== 3) dispx z ======
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h = figure('visible','off');% create a hidden figure
            pdeplot3D(dvcZOI,'ColorMapData',U(3:3:end),'FaceAlpha',0.5);
            
            title('$z$-displacement $w$','FontWeight','Normal','Interpreter','latex');
            set(gca,'fontsize',18); axis on;  axis equal;  % view([90 -90])
            xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            
            a = gca; a.TickLabelInterpreter = 'latex';
            b = colorbar; b.TickLabelInterpreter = 'latex';
            
            set(gcf,'color','w'); colormap jet; colorbar; box on
            
            frame = getframe(h);
            im_w = frame2im(frame);
            
            
        end
        
    end
    
end