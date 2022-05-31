classdef SERIAL2D
    
    % %%%%%%%%%%%%%%%%%% SerialTrack (2D cumulative mode) %%%%%%%%%%%%%%%%%
    % Main file of code "SerialTrack"
    % ***********************************************
    % Dimension:            2D
    % Particle rigidity:    hard
    % Tracking mode:        cumulative
    % -----------------------------------------------
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
        Type = 'SERIAL2D';
        ProgramFolder = 'SerialTrack2D';
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
        uv_B2A_prev = [];  % cumulative displacement
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
        function obj = SERIAL2D()
            %SERIAL2D Constructor
            % add the SERIAL2D functions to the path and remove others
            AddAndRemovePathFolders(obj,1)
            
            % create the MPTPara structure
            MPTPara.ImgRefMask = [];
            MPTPara.gridxyROIRange.gridy = [];
            MPTPara.gridxyROIRange.gridx = [];
            MPTPara.LoadImgMethod = 0;
            MPTPara.ImgSize = [];
            MPTPara.f_o_s = Inf;              % Size of search field: max(|u|,|v|)
            MPTPara.n_neighborsMax = 25;     % Max # of neighboring particles
            MPTPara.n_neighborsMin = 1;      % Min # of neighboring particles
            MPTPara.locSolver = 1;           % Local solver: 1-topology-based feature; 2-histogram-based feature first and then topology-based feature;
            MPTPara.gbSolver = 3;            % Global step solver: 1-moving least square fitting; 2-global regularization; 3-ADMM iterations
            MPTPara.smoothness = 1e-2;       % Coefficient of regularization
            MPTPara.outlrThres = 2;          % Threshold for removing outliers in TPT
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
            BeadPara.winSize = [5, 5];      % By default
            BeadPara.dccd = [1,1];          % By default
            BeadPara.abc = [1,1];           % By default
            BeadPara.forloop = 1;           % By default
            BeadPara.randNoise = 1e-7;      % By default
            BeadPara.PSF = [];              % PSF function; Example: PSF = fspecial('disk', BeadPara.beadSize-1 ); % Disk blur
            BeadPara.distMissing = 2;       % Distance threshold to check whether particle has a match or not
            BeadPara.color = 'white';       % Bead color
            
            obj.MPTPara = MPTPara;
            obj.BeadPara = BeadPara;
            
        end
        
        function obj = runSERIAL2D(obj,ImgSeqNum)
            
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
            
            [matches_A2B,u_B2A_curr_refB,obj.track_A2B] = f_track_serial_match2D( obj.parCoordA, obj.parCoordB, ...
                'f_o_s',obj.MPTPara.f_o_s, 'n_neighborsMax', obj.MPTPara.n_neighborsMax, 'n_neighborsMin', obj.MPTPara.n_neighborsMin, ...
                'locSolver',obj.MPTPara.locSolver,'gbSolver', obj.MPTPara.gbSolver, 'smoothness', obj.MPTPara.smoothness, ...
                'outlrThres', obj.MPTPara.outlrThres, 'maxIterNum', obj.MPTPara.maxIterNum, ...
                'iterStopThres', obj.MPTPara.iterStopThres, 'usePrevResults', obj.MPTPara.usePrevResults, ...
                'strain_n_neighbors',obj.MPTPara.strain_n_neighbors, 'strain_f_o_s',obj.MPTPara.strain_f_o_s, ...
                'gridxyROIRange',obj.MPTPara.gridxyROIRange, 'parCoordB_prev',obj.parCoord_prev(2:end), ...
                'uv_B2A_prev',obj.uv_B2A_prev, 'ImgSeqNum',ImgSeqNum, ...
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
             
            plotCone2(obj.parCoordB(:,1),obj.parCoordB(:,2),disp_A2B_parCoordB(:,1),disp_A2B_parCoordB(:,2));
            set(gca,'fontsize',18); box on; axis equal; axis tight; view(2); set(gca,'YDir','reverse');
            title('Tracked displacements','fontweight','normal');
            xlabel(''); ylabel(''); cb = colorbar; set(cb,'fontsize',18);
           
            frame = getframe(h); 
            obj.im_cone{ImgSeqNum} = frame2im(frame);
            
            
            %% %%%%% Compute F = def grad = grad_u = (grad_x X - I) %%%%%
            
            % strain_f_o_s: size of virtual strain gauge
            % strain_n_neighbors: # of neighboring particles used in strain gauge
            
            %%%%% Strain based on Moving Least Square Fitting in deformed configuration %%%%%
            [XY_refB,U_B2A_refB,F_B2A_refB] = funCompDefGrad2(disp_A2B_parCoordB, obj.parCoordB, obj.MPTPara.strain_f_o_s, obj.MPTPara.strain_n_neighbors);
            
            %%%%% Strain based on Moving Least Square Fitting in reference configuration %%%%%
            [XY_refA,U_A2B_refA,F_A2B_refA] = funCompDefGrad2(disp_A2B_parCoordB, obj.parCoordB-disp_A2B_parCoordB, obj.MPTPara.f_o_s, obj.MPTPara.strain_n_neighbors);
            
            
            
            %% %%%%% Store results %%%%%
            uv_B2A_refB = u_B2A_curr_refB;
            
            resultDisp_tmp = struct('parCoordA',obj.parCoordA,'parCoordB',obj.parCoordB,'track_A2B', ...
                obj.track_A2B,'disp_A2B_parCoordB',disp_A2B_parCoordB);
            
            resultDefGrad_tmp = struct('XY_refA',XY_refA,'U_A2B_refA',U_A2B_refA,'F_A2B_refA',F_A2B_refA, ...
                'XY_refB',XY_refB,'U_B2A_refB',U_B2A_refB,'F_B2A_refB',F_B2A_refB);

            % store results
            obj.parCoord_prev{ImgSeqNum} = obj.parCoordB;
            obj.uv_B2A_prev{ImgSeqNum-1} = uv_B2A_refB;  % cumulative displacement
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
                    x{1}{ImgSeqNum} = locateBeads(double(currImg)/max(double(currImg(:))),obj.BeadPara); % Detect particles
                    x{1}{ImgSeqNum} = radial2center(double(currImg)/max(double(currImg(:))),x{1}{ImgSeqNum},obj.BeadPara); % Localize particles
                    
                %%%%% Method 2: LoG operator (modified TracTrac code) %%%%%
                case 'LoG'
                    x{1}{ImgSeqNum} = f_detect_particles(double(currImg)/max(double(currImg(:))),obj.BeadPara);
                    
                %%%%% Method 3: sift code %%%%%
                case  'sift' % % Not fully implemented yet
                    [~,descriptors,locs] = sift(currImg);
                    x{1}{ImgSeqNum} = locs(:,1:2) +0*descriptors;
                    
            end

            x{1}{ImgSeqNum} = x{1}{ImgSeqNum} + [obj.MPTPara.gridxyROIRange.gridx(1)-1, obj.MPTPara.gridxyROIRange.gridy(1)-1];
            parCoord = x{1}{ImgSeqNum};
            
            %%%%% Remove parCoord outside the image area %%%%%
            for tempi=1:2
                parCoord( parCoord(:,tempi)>size(currImg,tempi), : ) = [];
            end
            
            for tempi=1:2
                parCoord(parCoord(:,tempi)<1, : ) = [];
            end
            
            % make a figure and plot the background image
            fig1=figure('visible','off');
            ax1=axes; ax1.Visible = 'off';
            try h1=imshow(imrotate(currImg,90),[],'InitialMagnification','fit');
            catch h1 =surf(imrotate(currImg,90),'EdgeColor','none','LineStyle','none');
            end
            % fix up the image properties
            axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');

            % make another figure axes and plot the color picture on top
            hold on;
            ax2=axes; ax2.Visible = 'off'; 
            h2=plot(parCoord(:,1)-5,parCoord(:,2),'ro');
            
            % fix up the image properties
            set(gca,'fontsize',18); box on; axis equal; axis tight; view(2); set(gca,'YDir','reverse');
            title('Beads Found','fontweight','normal');
            xlabel(''); ylabel('');
            
            % set the alpha of the colorimage with an alphamap or constant value
            alpha(h2,0.5);
            
            % now combine the two axes together
            linkaxes([ax1,ax2]);  % Link axes together
            ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
            colormap(ax1,'gray'); % Give each one its own colormap
            set([ax1,ax2],'Position',[.17 .11 .685 .815]);
            ax1.Visible = 'off';
            
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
        
        function obj = BlackParticles(obj,ImgSeqNum)
            
            % this function takes in an image with white background and
            % black beads and saves that image into obj.Imges_used for
            % future use
            currImg = obj.Imges_used{ImgSeqNum};
            ImgGauss = imgaussfilt(imgaussfilt(currImg,1),1); % figure, imshow(uint16(ImgGauss));
            ImgGauss(ImgGauss > obj.BeadPara.thres*max(double(currImg(:)))) = 0;
            bw = imbinarize(uint16(ImgGauss),'adaptive','ForegroundPolarity','dark','Sensitivity',0.8); % figure, imshow(bws2);
            bws2 = bwareaopen(bw,round(pi*obj.BeadPara.minSize^2)); % remove all object containing fewer than BeadPara.minSize
            removeobjradius = obj.BeadPara.minSize; % fill a gap in the pen's cap
            se = strel('disk',removeobjradius);
            bws2 = imclose(bws2,se);
            obj.Imges_used{ImgSeqNum} = double(bws2); % figure, imshow(uint8(currImg2));
        end
        
        function obj = ApplyROI(obj)
            obj.MPTPara.ImgRefMask = obj.ImgMaskFiles{1};
            for ImgSeqNum = 1 : obj.NumberOfImages
                % apply a mask to the image to ensure no particles in these
                % areas are detected
                Img =obj.Imges_raw{ImgSeqNum}.*obj.ImgMaskFiles{ImgSeqNum};
                
                %%%%% Crop the ROI region %%%%%
                obj.Imges_used{ImgSeqNum} = Img(obj.MPTPara.gridxyROIRange.gridx(1):obj.MPTPara.gridxyROIRange.gridx(2), ...
                    obj.MPTPara.gridxyROIRange.gridy(1):obj.MPTPara.gridxyROIRange.gridy(2));
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
                        obj.parCoordATraj{parInd}(ImgSeqNum-1,1:2) = obj.parCoordB(obj.track_A2B(parInd),1:2);
                    else
                        obj.parCoordATraj{parInd}(ImgSeqNum-1,1:2) = [nan,nan];
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
            
            %%%%% Plot tracked trajectories %%%%%
            
            h = figure('visible','off');% create a hidden figure 
            xstep = obj.MPTPara.xstep;
            
            for parInd = 1:size(obj.parCoordA,1)
                try
                    wayPoints = obj.parCoordATraj{parInd};
                    if length(wayPoints(isnan(wayPoints(:,1))==1,1)) < 3
                        if (size(obj.resultDisp,1)+1)<4
                            hold on; 
                            line(xstep*wayPoints(isnan(wayPoints(:,1))<1,1),xstep*wayPoints(isnan(wayPoints(:,1))<1,2),'linewidth',1.2); 
                            view(3); % straight lines
                        else
                            hold on; 
                            fnplt(cscvn(xstep*wayPoints(isnan(wayPoints(:,1))<1,:)'),'',1.);
                        end
                    end
                    % if sum(1-isnan(wayPoints(:,1)))>1 % Don't show if there is only one point on the trajectory
                    %     hold on; plot(xstep*wayPoints(:,1),xstep*wayPoints(:,2),'.','markersize',5);
                    % end
                catch
                end
            end
                        
            set(gca,'fontsize',18); view(2); box on; axis equal; axis tight;
            title('Tracked particle trajectory','fontweight','normal');
            set(gca,'YDir','reverse'); xlabel('x'); ylabel('y');
            axis(xstep*[obj.MPTPara.gridxyROIRange.gridx(1), obj.MPTPara.gridxyROIRange.gridx(2), ...
                obj.MPTPara.gridxyROIRange.gridy(1), obj.MPTPara.gridxyROIRange.gridy(2) ]);
            
            frame = getframe(h);
            obj.im_traj{ImgSeqNum} = frame2im(frame);
            
        end
        
        function obj = plotDispAndStrainMaps(obj)
            
            % make this a for loop for all images
            
            for ImgSeqNum = 2:obj.NumberOfImages
                xstep = obj.MPTPara.xstep;
                tstep = obj.MPTPara.tstep;
                
                %%%%% Previously tracked displacement field %%%%%
                resultDispCurr = obj.resultDisp{ImgSeqNum-1};
                resultDefGradCurr = obj.resultDefGrad{ImgSeqNum-1};
                disp_A2B_parCoordB = resultDispCurr.disp_A2B_parCoordB;
                obj.parCoordB = resultDispCurr.parCoordB;
                
                %%%%% Interpolate scatterred data to gridded data %%%%%
                sxy = min([round(0.5*obj.MPTPara.f_o_s),20])*[1,1]; % Step size for griddata
                smoothness = 1e-3; % Smoothness for regularization; "smoothness=0" means no regularization
                
                [x_Grid_refB,y_Grid_refB,u_Grid_refB]=funScatter2Grid2D(obj.parCoordB(:,1),obj.parCoordB(:,2),disp_A2B_parCoordB(:,1),sxy,smoothness);
                [~,~,v_Grid_refB]=funScatter2Grid2D(obj.parCoordB(:,1),obj.parCoordB(:,2),disp_A2B_parCoordB(:,2),sxy,smoothness);
                obj.MPTPara.ImgRefMask = obj.ImgMaskFiles{1};
                % Apply ROI image mask
                %             [u_Grid_refB, v_Grid_refB] = funRmROIOutside(x_Grid_refB,y_Grid_refB,obj.MPTPara.ImgRefMask,u_Grid_refB,v_Grid_refB);
                % Build a displacement vector
                uv_Grid_refB_Vector=[u_Grid_refB(:),v_Grid_refB(:)]';
                uv_Grid_refB_Vector=uv_Grid_refB_Vector(:);
                % Calculate deformation gradient
                D_Grid = funDerivativeOp(size(x_Grid_refB,1),size(x_Grid_refB,2),mean(sxy)); % Central finite difference operator
                F_Grid_refB_Vector=D_Grid*uv_Grid_refB_Vector; % {F}={D}{U}
                
                %%%%% Generate an FE-mesh %%%%%
                [coordinatesFEM_refB,elementsFEM_refB] = funMeshSetUp(x_Grid_refB*xstep,y_Grid_refB*xstep);
                
                
                %%%%% Cone plot grid data: displecement %%%%%
                h = figure('visible','off');% create a hidden figure
                plotCone2(x_Grid_refB*xstep,y_Grid_refB*xstep,u_Grid_refB*xstep,v_Grid_refB*xstep );
                set(gca,'fontsize',18); view(2); box on; axis equal; axis tight; set(gca,'YDir','reverse');
                title('Tracked incremental displacement','fontweight','normal');
                axis(xstep*[obj.MPTPara.gridxyROIRange.gridx(1), obj.MPTPara.gridxyROIRange.gridx(2), ...
                    obj.MPTPara.gridxyROIRange.gridy(1), obj.MPTPara.gridxyROIRange.gridy(2) ]);
                
                frame = getframe(h);
                obj.im_cone{ImgSeqNum} = frame2im(frame);
                
                obj.elementsFEM{ImgSeqNum} = elementsFEM_refB;
                obj.coordinatesFEM{ImgSeqNum} = coordinatesFEM_refB*xstep;
                obj.U{ImgSeqNum} = uv_Grid_refB_Vector*xstep;
                obj.F{ImgSeqNum} = F_Grid_refB_Vector;
                
                %%%%% Cone plot grid data: displacement %%%%%
                [obj.im_disp.u{ImgSeqNum},obj.im_disp.v{ImgSeqNum}] = obj.Plotdisp_show_obj(obj,obj.U{ImgSeqNum}, obj.coordinatesFEM{ImgSeqNum}, obj.elementsFEM{ImgSeqNum},[],'NoEdgeColor');
                
                %%%%% Cone plot grid data: infinitesimal strain %%%%%
                [obj.im_strain.exx{ImgSeqNum}, obj.im_strain.exy{ImgSeqNum}, obj.im_strain.eyy{ImgSeqNum}] = Plotstrain_show_obj(obj,obj.F{ImgSeqNum} , obj.coordinatesFEM{ImgSeqNum}, obj.elementsFEM{ImgSeqNum},[],'NoEdgeColor',xstep,tstep);
            end
            obj.ResultDisp = obj.U;
            obj.ResultDefGrad = obj.F;
            obj.ResultFEMeshEachFrame_coordinatesFEM = obj.coordinatesFEM;
            obj.ResultFEMeshEachFrame_elementsFEM = obj.elementsFEM;
        end
        
        function obj = findgridxy(obj)
            
            for ImgSeqNum = 1 : obj.NumberOfImages
                ROI = obj.ImgMaskFiles{ImgSeqNum};
                bounds = regionprops(ROI,'BoundingBox');
                if length(bounds)== 1
                    b = bounds.BoundingBox;
                    gridx = round([b(1), b(1)+b(3)-1]);
                    gridy = round([b(2), b(2)+b(4)-1]);
                end
                
                imgszx = obj.MPTPara.ImgSize(1);
                imgszy = obj.MPTPara.ImgSize(2);
                
                %%%%% Modify gridxy to make sure it is within the image domain %%%%%
                if gridx(1)<1, gridx(1) = 1; end
                if gridy(1)<1, gridy(1) = 1; end
                if gridx(2)>imgszx, gridx(2) = imgszx; end
                if gridy(2)>imgszy, gridy(2) = imgszy; end
                
                % Store TPTpara
                obj.MPTPara.gridxyROIRange.gridy = gridy;
                obj.MPTPara.gridxyROIRange.gridx = gridx;
            end
            
        end
                
        function  [im_exx,im_exy,im_eyy] = Plotstrain_show_obj( obj, F,coordinatesFEM,elementsFEM,varargin)
            %FUNCTION Plotstrain_show(F,coordinatesFEM,elementsFEM)
            % To plot DIC solved strain components
            % ----------------------------------------------
            %
            %   INPUT: F                 Deformation gradient tensor:
            %                            F = [F11_node1, F21_node1, F12_node1, F22_node1, ... , F11_nodeN, F21_nodeN, F12_nodeN, F22_nodeN]';
            %          coordinatesFEM    FE mesh coordinates
            %          elementsFEM       FE mesh elements
            %          DICpara           chosen DIC parameters
            %          EdgeColorOrNot    show edge color or not
            
            %   OUTPUT: Plots of exx, eyy, and exy strain fields.
            %
            % ----------------------------------------------
            % Author: Jin Yang.
            % Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
            % Last date modified: 2020.12
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
            %% Initialization
            warning off;  F = full(F);
            
            %%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            DICpara=[]; EdgeColorOrNot=[];
            
            try
                DICpara=vargin{1};
                try
                    EdgeColorOrNot=vargin{2};
                catch
                end
            catch
            end
            
            %%%%% convert pixel unit to the physical world unit %%%%%
            try um2px = DICpara.um2px;
            catch um2px = 1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            try EdgeColorOrNot = EdgeColorOrNot;
            catch EdgeColorOrNot = 'EdgeColor';
            end
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ====== 1) strain exx ======
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h = figure('visible','off');% create a hidden figure
            show([],elementsFEM(:,1:4),coordinatesFEM,F(1:4:end),EdgeColorOrNot);
            title('Strain $e_{11}$','fontweight','normal','Interpreter','latex'); set(gca,'fontsize',18);
            view(2); axis tight; axis equal;  box on;  set(gcf,'color','w');
            colorbar; colormap jet;
            
            if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
            end
            
            a = gca; a.TickLabelInterpreter = 'latex';
            b = colorbar; b.TickLabelInterpreter = 'latex';
            
            frame = getframe(h); 
            im_exx = frame2im(frame);
            
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ====== 2) strain exy ======
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             h = figure('visible','off');% create a hidden figure
             show([],elementsFEM(:,1:4),coordinatesFEM,0.5*(F(2:4:end)+F(3:4:end)),EdgeColorOrNot);
            title('Strain $e_{12}$','fontweight','normal','Interpreter','latex'); set(gca,'fontsize',18);
            view(2); axis tight; axis equal;  box on;  set(gcf,'color','w');
            colorbar; colormap jet;
            
            if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
            end
            
            a = gca; a.TickLabelInterpreter = 'latex';
            b = colorbar; b.TickLabelInterpreter = 'latex';
            
            frame = getframe(h); 
            im_exy = frame2im(frame);
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ====== 3) strain eyy ======
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h = figure('visible','off');% create a hidden figure
            show([],elementsFEM(:,1:4),coordinatesFEM,F(4:4:end),EdgeColorOrNot);
            title('Strain $e_{22}$','fontweight','normal','Interpreter','latex'); set(gca,'fontsize',18);
            view(2); axis tight; axis equal;  box on;  set(gcf,'color','w');
            colorbar; colormap jet;
            
            if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
            end
            
            a = gca; a.TickLabelInterpreter = 'latex';
            b = colorbar; b.TickLabelInterpreter = 'latex';
            
            frame = getframe(h); 
            im_eyy = frame2im(frame);
            
        end
        
        function  [im_u, im_v] = Plotdisp_show_obj(obj, U,coordinatesFEM,elementsFEM,varargin)
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
            
            
            %% Initialization
            warning off;  U = full(U);
            % run('./img_Vito/Black_rainbow.m');
            
            %%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            DICpara=[]; EdgeColorOrNot=[];
            
            try
                DICpara=vargin{1};
                try
                    EdgeColorOrNot=vargin{2};
                catch
                end
            catch
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%% convert pixel unit to the physical world unit %%%%%
            try um2px = DICpara.um2px;
            catch um2px = 1;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            try EdgeColorOrNot = EdgeColorOrNot;
            catch EdgeColorOrNot = 'EdgeColor';
            end
            
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ====== 1) dispx u ======
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h = figure('visible','off');% create a hidden figure
            show([],elementsFEM(:,1:4),coordinatesFEM,U(1:2:end),EdgeColorOrNot);
            
            title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');
            view(2); set(gca,'fontsize',18); axis tight; axis equal; colorbar; % view([90 -90])
            if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
            end
            set(gcf,'color','w'); colormap turbo; box on;
            a = gca; a.TickLabelInterpreter = 'latex';
            b = colorbar; b.TickLabelInterpreter = 'latex';
            
            frame = getframe(h); 
            im_u = frame2im(frame);
            
            %%%%%% TODO: manually modify colormap and caxis %%%%%%
            % colormap(jet);  % D Sample
            % caxis([-0.1,0.1]) % foam
            % caxis([-0.004,0.004]); % Sample 12
            % caxis([0,0.1]);
            % caxis([-0.56,0.56]); % colorturbo = turbo(64); colorturbo(:,2)=0; colormap(turbo)
            %  colormap(black_rainbow_plus );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ====== 2) dispx v ======
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h = figure('visible','off');% create a hidden figure
            show([],elementsFEM(:,1:4),coordinatesFEM,U(2:2:end),EdgeColorOrNot);
            
            title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
            view(2); set(gca,'fontsize',18); axis tight; axis equal; colorbar;
            if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
            end
            set(gcf,'color','w'); colormap turbo; box on;
            a = gca; a.TickLabelInterpreter = 'latex';
            b = colorbar; b.TickLabelInterpreter = 'latex';
            
            frame = getframe(h); 
            im_v = frame2im(frame);
            
            %%%%%% TODO: manually modify colormap and caxis %%%%%%
            % colormap(jet);  % D Sample
            % caxis([-0.1,0.1]) % foam
            % caxis([-0.004,0.004]); % Sample 12
            % caxis([0,0.1]);
            % caxis([-0.12,0.12]); colormap(black_rainbow_plus);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
        end
        
    end
    
end