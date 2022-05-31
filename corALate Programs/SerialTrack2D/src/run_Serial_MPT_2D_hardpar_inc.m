% %%%%%%%%%%%%%%%%%% SerialTrack (2D incremental mode) %%%%%%%%%%%%%%%%%
% Main file of code "SerialTrack"
% ***********************************************
% Dimension:            2D
% Particle rigidity:    hard 
% Tracking mode:        incremental
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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
%%%%% Load 2D images %%%%%
[file_name,Img,MPTPara_temp] = funReadImage2; close all;
disp('%%%%%% Load images: Done! %%%%%%'); fprintf('\n');

%%%%% Update MPTPara %%%%%
MPTPara.gridxyROIRange = MPTPara_temp.gridxyROIRange;
MPTPara.LoadImgMethod = MPTPara_temp.LoadImgMethod;
MPTPara.ImgSize = MPTPara_temp.ImgSize;

%%%%% Load image mask file(s) %%%%%
try
    if MaskFileLoadingMode == 1
        %%%%% Load one image mask file for all frames %%%%%
        [file_mask_name,ImgMaskFiles] = funReadImageMask2; close all;
    elseif MaskFileLoadingMode == 2
        %%%%% Load multiple image mask files %%%%%
        [file_mask_name,ImgMaskFiles] = funReadImageMask2; close all;
    elseif MaskFileLoadingMode == 3
        try load(im_roi_mask_file_path); catch; end
        try MPTPara.ImgRefMask = im_roi'; % Load stored image roi if existed
        catch, MPTPara.ImgRefMask = ones(size(Img{1})); % Set up default image mask file
        end
    else
        MPTPara.ImgRefMask = ones(size(Img{1})); % Set up default image mask file
    end
    disp('%%%%%% Load image mask file: Done! %%%%%%'); fprintf('\n');
catch
    MPTPara.ImgRefMask = ones(size(Img{1})); % Set up default image mask file
end
  
  
   
%% ====== Detect particles ======
%%%%% Particle detection parameters %%%%%
%%%%%%%%%% Atom dance %%%%%%%%%%%%%
%%%%% Bead Parameter %%%%%
% BeadPara.thres = 0.35 ;           % Threshold for detecting particles
% BeadPara.beadSize = 0;          % Estimated radius of a single particle
% BeadPara.minSize = 2;           % Minimum radius of a single particle
% BeadPara.maxSize = 40;          % Maximum radius of a single particle
% BeadPara.winSize = [5, 5];      % By default
% BeadPara.dccd = [1,1];          % By default
% BeadPara.abc = [1,1];           % By default
% BeadPara.forloop = 1;           % By default
% BeadPara.randNoise = 1e-7;      % By default
% BeadPara.PSF = [];              % PSF function; Example: PSF = fspecial('disk', BeadPara.beadSize-1 ); % Disk blur
% BeadPara.distMissing = 2;       % Distance threshold to check whether particle has a match or not 
% BeadPara.color = 'black';       % Bead color 

%%%%%%%%%% Pipe %%%%%%%%%%%%%
%%%%% Bead Parameter %%%%%
% BeadPara.thres = 0.5;           % Threshold for detecting particles
% BeadPara.beadSize = 3;          % Estimated radius of a single particle
% BeadPara.minSize = 2;           % Minimum radius of a single particle
% BeadPara.maxSize = 20;          % Maximum radius of a single particle
% BeadPara.winSize = [5, 5];      % By default
% BeadPara.dccd = [1,1];          % By default
% BeadPara.abc = [1,1];           % By default
% BeadPara.forloop = 1;           % By default
% BeadPara.randNoise = 1e-7;      % By default
% BeadPara.PSF = [];              % PSF function; Example: PSF = fspecial('disk', BeadPara.beadSize-1 ); % Disk blur
% BeadPara.distMissing = 2;       % Distance threshold to check whether particle has a match or not 
% BeadPara.color = 'white';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ImgSeqNum = 1; % First reference image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Apply the uploaded mask file %%%%%
try
    if (MaskFileLoadingMode==1) || (MaskFileLoadingMode==3)
        currImg = Img{ImgSeqNum}.*MPTPara.ImgRefMask;
    elseif MaskFileLoadingMode==2
        currImg = Img{ImgSeqNum}.*ImgMaskFiles{ImgSeqNum};
    else
        currImg = Img{ImgSeqNum};  
    end
catch
    currImg = Img{ImgSeqNum};  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currImg = currImg(MPTPara.gridxyROIRange.gridx(1):MPTPara.gridxyROIRange.gridx(2), ...
                  MPTPara.gridxyROIRange.gridy(1):MPTPara.gridxyROIRange.gridy(2));

%%%%% If PSF is non-empty, SerialTrack performs deconvolution %%%%%
if ~isempty(BeadPara.PSF)
    currImg = deconvlucy(currImg,BeadPara.PSF,6);
    disp(['----- Deconvolution frame #',num2str(ImgSeqNum),' ------']);
end
%%%% figure, imshow(currImg);

 
%%%%% Pre-process bead image if bead color is black %%%%%
if strcmp(BeadPara.color,'black')
    ImgGauss = imgaussfilt(imgaussfilt(currImg,1),1); % figure, imshow(uint16(ImgGauss));
    ImgGauss(ImgGauss > BeadPara.thres*max(double(currImg(:)))) = 0;
    bw = imbinarize(uint16(ImgGauss),'adaptive','ForegroundPolarity','dark','Sensitivity',0.8); % figure, imshow(bws2);
    bws2 = bwareaopen(bw,round(pi*BeadPara.minSize^2)); % remove all object containing fewer than BeadPara.minSize
    removeobjradius = BeadPara.minSize; % fill a gap in the pen's cap
    se = strel('disk',removeobjradius);
    bws2 = imclose(bws2,se);
    currImg2 = double(bws2); % figure, imshow(uint8(currImg2));
else
    currImg2 = currImg;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Several methods to detect particles %%%%%
%%%%% Method 1: TPT code %%%%%
% x{1}{ImgSeqNum} = locateBeads(double(currImg2)/max(double(currImg2(:))),BeadPara); % Detect particles
% x{1}{ImgSeqNum} = radial2center(double(currImg2)/max(double(currImg2(:))),x{1}{ImgSeqNum},BeadPara); % Localize particles
% ----------------------------
%%%%% Method 2: LoG operator (modified TracTrac code) %%%%%
x{1}{ImgSeqNum} = f_detect_particles(double(currImg2)/max(double(currImg2(:))),BeadPara);
% ----------------------------
%%%%% Method 3: sift code %%%%%
% [~,descriptors,locs] = sift(currImg2);
% x{1}{ImgSeqNum} = locs(:,1:2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% Store particle positions as "parCoordA" %%%%%
x{1}{ImgSeqNum} = x{1}{ImgSeqNum} + [MPTPara.gridxyROIRange.gridx(1)-1, MPTPara.gridxyROIRange.gridy(1)-1];
parCoordA = x{1}{ImgSeqNum};

%%%%% Remove bad parCoord outside the image area %%%%%
for tempi=1:2, parCoordA( parCoordA(:,tempi)>size(Img{ImgSeqNum},tempi), : ) = []; end
for tempi=1:2, parCoordA( parCoordA(:,tempi)<1, : ) = []; end
 
%%%%% Plot %%%%%
figure, imshow(imread(file_name{1,1}))
hold on; plot( parCoordA(:,1), parCoordA(:,2), 'r.');
view(2); box on; axis equal; axis tight; set(gca,'fontsize',18);
title('Detected particles in ref image','fontweight','normal');
 
%%%%% Report detected beads # %%%%%
disp(['Detected particle # in ref image: ',num2str(size(parCoordA,1))]);
disp('%%%%%% Detect particles: Done! %%%%%%'); fprintf('\n');
 


%% %%%%% Initialization %%%%%
%%%%%  MPT Parameter %%%%%
% MPTPara.f_o_s = 30;              % Size of search field: max(|u|,|v|)
% MPTPara.n_neighborsMax = 25;     % Max # of neighboring particles
% MPTPara.n_neighborsMin = 1;      % Min # of neighboring particles
% MPTPara.locSolver = 1;           % Local solver: 1-topology-based feature; 2-histogram-based feature first and then topology-based feature;
% MPTPara.gbSolver = 3;            % Global step solver: 1-moving least square fitting; 2-global regularization; 3-ADMM iterations
% MPTPara.smoothness = 1e-2;       % Coefficient of regularization
% MPTPara.outlrThres = 2;          % Threshold for removing outliers in TPT
% MPTPara.maxIterNum = 20;         % Max ADMM iteration number
% MPTPara.iterStopThres = 1e-2;    % ADMM iteration stopping threshold
% MPTPara.strain_n_neighbors = 20; % # of neighboring particles used in strain gauge
% MPTPara.strain_f_o_s = 60;       % Size of virtual strain gauge
% MPTPara.usePrevResults = 0;      % Whether use previous results or not: 0-no; 1-yes;

%%%%%% To store results %%%%%
parCoord_prev = cell(length(Img)-1,1);      parCoord_prev{1} = parCoordA;
track_A2B_prev = cell(length(Img)-1,1);     track_B2A_prev = cell(length(Img)-1,1);
uv_B2A_prev = cell(length(Img)-1,1);

 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ImgSeqNum = 2 : length(Img)  % "ImgSeqNum" is the frame index
    
    disp(['====== Frame #',num2str(ImgSeqNum),' ======']);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Apply the uploaded mask file %%%%%
    clear defImg;
    try 
        if (MaskFileLoadingMode==1) || (MaskFileLoadingMode==3)
            defImg = Img{ImgSeqNum}.*MPTPara.ImgRefMask;
        elseif MaskFileLoadingMode==2
            defImg = Img{1,ImgSeqNum}.*ImgMaskFiles{1,ImgSeqNum};
        else
            defImg = Img{ImgSeqNum};
        end
    catch
       defImg = Img{ImgSeqNum};
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% SerialTrack particle tracking %%%%%
    [parCoordB_temp,uv_B2A_temp,~,~,track_A2B_temp,track_B2A_temp] = fun_SerialTrack_2D_HardPar( ...
        ImgSeqNum,defImg,BeadPara,MPTPara,parCoord_prev{ImgSeqNum-1},parCoord_prev(2:end),uv_B2A_prev);
     
    %%%%% Store results %%%%%
    parCoord_prev{ImgSeqNum} = parCoordB_temp;
    uv_B2A_prev{ImgSeqNum-1} = uv_B2A_temp; % incremental displacement
    track_A2B_prev{ImgSeqNum-1} = track_A2B_temp;
    track_B2A_prev{ImgSeqNum-1} = track_B2A_temp;
     
end
  
  

%% %%%%% Incremental tracking ratio %%%%%
disp('%%%%% Calculate incremental tracking ratio %%%%%'); fprintf('\n');
track_ratio = zeros(length(Img)-1,1);
DefType = 'exp'; defList = [2:1:length(Img)]';
  
for ImgSeqNum = 2 : length(Img)
    track_A2B = track_A2B_prev{ImgSeqNum-1}; 
    track_ratio(ImgSeqNum-1) = length(track_A2B(track_A2B>0))/size(parCoord_prev{ImgSeqNum},1);      
end
 
fig=figure; ax=axes; hold on; plot(defList,track_ratio,'r^-.','linewidth',1);
adjust_fig(fig,ax,'','',''); box on; title('');
xlabel('Frame #'); ylabel('Incremental tracking ratio');
axis([2,length(file_name),0,1]);


%%%%% Save results %%%%%
disp('%%%%%% ALTPT hard particle tracking: Done! %%%%%%'); fprintf('\n');
results_file_name = 'results_hardpar.mat';
mkdir results
save(['./results/' results_file_name],'parCoord_prev','uv_B2A_prev','track_A2B_prev','track_B2A_prev');
 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%% Visualize tracked incremental displacement of each frame %%%%%
disp('%%%%% Plot tracked incremental deformations %%%%%'); fprintf('\n');

%%%%% Experimental parameters %%%%%
try xstep = MPTPara.xstep; catch, xstep = 1; end % unit: um/px
try tstep = MPTPara.tstep; catch, tstep = 1; end % unit: us  
% ImgSeqNum  % Frame #
 
%%%%% Plot tracked incremental displacement field %%%%%
%%%%% Make a video %%%%%
v = VideoWriter('video_2D_inc.avi'); v.FrameRate = 5; open(v); figure,

for ImgSeqNum = 2:length(Img)
    
    % Displacement from tracked particles on deformed frame
    disp_A2B_parCoordB = -uv_B2A_prev{ImgSeqNum-1};
    parCoordB = parCoord_prev{ImgSeqNum};

    %%%%% Plot displacements %%%%%
    clf, plotCone2(parCoordB(:,1)*xstep,parCoordB(:,2)*xstep,disp_A2B_parCoordB(:,1)*xstep/tstep ,disp_A2B_parCoordB(:,2)*xstep/tstep );
    set(gca,'fontsize',18); view(2); box on; axis equal; axis tight; set(gca,'YDir','reverse');
    title(['Tracked velocity (#',num2str(ImgSeqNum),')'],'fontweight','normal');
    xlabel('x'); ylabel('y');
    axis(xstep*[MPTPara.gridxyROIRange.gridx(1), MPTPara.gridxyROIRange.gridx(2), ...
          MPTPara.gridxyROIRange.gridy(1), MPTPara.gridxyROIRange.gridy(2) ]);
     
    frame = getframe(gcf);
    writeVideo(v,frame);
    
end
close(v);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%% Compute trajectory %%%%%
disp('%%%%% Merge tracked trajectory segments %%%%%'); fprintf('\n');

%%%%% Initialization %%%%%
try xstep = MPTPara.xstep; catch xstep = 1; end % um2px ratio
try tstep = MPTPara.tstep; catch tstep = 1; end % time gap between consecutive frames
trajInd = 0; % index of trajectory segments

%%%%% These variables are defined in the main execute file(s) %%%%%
% distThres = 3; % distance threshold to connect split trajectory segments
% extrapMethod = 'pchip';  % extrapolation scheme to connect split trajectory segments
%                          % suggestion: 'nearest' for Brownian motion                          
% minTrajSegLength = 4;    % the minimum length of trajectory segment that will be extrapolate 
% maxGapTrajSeqLength = 0; % the max frame# gap between connected trajectory segments


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Compute and collect all trajectory segments %%%%%
for tempk = 1 : size(parCoord_prev,1)  % Find trajectories passing through particles in frame #tempk
    
    parCoordCurr = parCoord_prev{tempk}; % Find trajectories passing parCoordCurr
    clear parCoordTrajCurr; parCoordTrajCurr = cell(size(parCoordCurr,1),1); % Initialize a cell structure to store #tempk trajectories
    
    % Add particles in frame #tempk to "parCoordTrajCurr"
    for tempj = 1:length(parCoordCurr)
        parCoordTrajCurr{tempj}(tempk,1:2) = parCoordCurr(tempj,1:2);
    end
    
    for ImgSeqNum = (1+tempk) : size(parCoord_prev,1) % For each later tracked incremental deformation
        
        parCoordB = parCoord_prev{ImgSeqNum}; % Current particle coordinates in frame #ImgSeqNum
        parCoordB_Ind = []; parCoordCurr_Ind = []; % Initialize index
        track_B2A_Curr = track_B2A_prev{ImgSeqNum-1}; % Current tracking in frame #ImgSeqNum
        for tempi = 1:length(track_B2A_Curr)
            try
                % Back propogation
                if ImgSeqNum > 2
                    for tempj = (ImgSeqNum-2) : -1 : tempk
                        track_B2A_Last = track_B2A_prev{tempj};
                        track_B2A_Curr(tempi) = track_B2A_Last(track_B2A_Curr(tempi));
                    end
                end
                % Check whether the trajectory has already been plotted previously
                if tempk>1
                    track_B2A_SecLast = track_B2A_prev{tempk-1};
                    track_B2A_SecLast_temp = track_B2A_SecLast(track_B2A_Curr(tempi)); % if this is not 0, means it's already been plotted
                else % Trajectories from first two frames (tempk=1) will be added
                    track_B2A_SecLast_temp = 0;
                end
                % Assign index values
                if (track_B2A_Curr(tempi)>0) && (track_B2A_SecLast_temp==0)
                    parCoordB_Ind = [parCoordB_Ind; tempi];
                    parCoordCurr_Ind = [parCoordCurr_Ind; track_B2A_Curr(tempi)];
                end
            catch
            end
        end
        
        for tempj = 1:length(parCoordCurr_Ind) % Add found trajectories to cell structure "parCoordTraj"
            parCoordTrajCurr{parCoordCurr_Ind(tempj)}(ImgSeqNum,1:2) = parCoordB(parCoordB_Ind(tempj),1:2);
        end
    end
    
    for parInd = 1:size(parCoordTrajCurr,1)
        wayPoints = parCoordTrajCurr{parInd};
        if ~isempty(wayPoints)
            wayPoints(wayPoints(:,1)==0,:) = wayPoints(wayPoints(:,1)==0,:)*nan;
            wayPoints = [wayPoints; nan(size(parCoord_prev,1)-size(wayPoints,1),2)];
            trajInd = trajInd + 1;
            parCoordTraj{trajInd} = wayPoints;
        end
    end
      
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Merge trajectory segments %%%%%
% Find the starting point and length of each trajectory segment
parCoordTrajPara = []; parCoordTraj = parCoordTraj(:);
for tempi = 1 : size(parCoordTraj,1)
    parCoordTrajCurr = parCoordTraj{tempi};
    parCoordTrajPara(tempi,1) = find(isnan(parCoordTrajCurr(:,1))==0, 1, 'first');
    parCoordTrajPara(tempi,2) = sum(1-isnan(parCoordTrajCurr(:,1)));
end

%%%%% Try to merge trajectory segments by extrapolating the particle position %%%%%
hbar = waitbar(0,'wait');
for tempMergeTime = 1:4 % Try to merge four times
    
    for tempm = 0:maxGapTrajSeqLength % tempm is the # of missing particles between trajectory segments
        
        for tempk = (size(parCoord_prev,1)-1) : -1 : minTrajSegLength % For trajectory segments with length [(size(parCoord_prev,1)-1) : -1 : minTrajSegLength]
            
            [row,col] = find( parCoordTrajPara(:,2)==tempk ); % Find trajectory segments whose "length==tempk"
            [row1,~] = find( parCoordTrajPara(:,2)<size(parCoord_prev,1)+1-tempk & parCoordTrajPara(:,2)>0 ); % Find trajectory setments with length requirement
             
            for tempi = 1:length(row) % For each trajectory segment whose "length==tempk"
                 
                tempWatibar = tempi/length(row)/4/(maxGapTrajSeqLength+1)/((size(parCoord_prev,1)-1)-minTrajSegLength+1) + ...
                    ((size(parCoord_prev,1)-1)-tempk+1)/4/(maxGapTrajSeqLength+1)/((size(parCoord_prev,1)-1)-minTrajSegLength+1) + ...
                    (tempm)/4/(maxGapTrajSeqLength+1) + (tempMergeTime-1)/4;
                
                waitbar(tempWatibar);
                 
                parCoordTrajMat = cell2mat( parCoordTraj );
                parCoordTrajCurr = parCoordTraj{row(tempi)}; % For each trajectory segment whose "length==tempk"
                parCoordTrajCurrPred_x = fillmissing(parCoordTrajCurr(:,1),extrapMethod); % fill missing using 'pchip' method
                parCoordTrajCurrPred_y = fillmissing(parCoordTrajCurr(:,2),extrapMethod); % fill missing using 'pchip' method
                
                % figure, plot(parCoordTrajCurrPred_x,parCoordTrajCurrPred_y,'-o');
                
                
                %%%%% Find all probable trajectory segments in the positive direction %%%%%
                if (sum(parCoordTrajPara(row(tempi),1:2))+tempm<(size(parCoord_prev,1)+1)) && (sum(parCoordTrajPara(row(tempi),1:2))+tempm>0)
                    
                    [row2,col2] = find( parCoordTrajPara(:,1) == sum(parCoordTrajPara(row(tempi),1:2))+tempm ); % starting point requirement
                    row3 = intersect(row1,row2);
                    tempxx = size(parCoord_prev,1)*(row3-1) +(sum(parCoordTrajPara(row(tempi),1:2)))+tempm; tempxx=tempxx'; tempxx=tempxx(:); % Index in parCoordTrajMat
                    tempyy = parCoordTrajMat(tempxx, 1:2);
                    
                    % hold on; plot(tempyy(:,1),tempyy(:,2),'.'); pause;
                    
                    tempzz = sqrt((tempyy(:,1)-parCoordTrajCurrPred_x( sum(parCoordTrajPara(row(tempi),1:2))+tempm )).^2 + ...
                        (tempyy(:,2)-parCoordTrajCurrPred_y( sum(parCoordTrajPara(row(tempi),1:2))+tempm )).^2);
                    
                    [tempzzmin,tempzzminind] = min(tempzz);
                    if tempzzmin < distThres % Find the continuing trajectory segment %JY!!!! threshold distance 3
                        
                        % Merge trajectory segment
                        parCoordTraj{row(tempi)}( parCoordTrajPara(row3(tempzzminind),1) : parCoordTrajPara(row3(tempzzminind),1)+parCoordTrajPara(row3(tempzzminind),2)-1, 1:2 ) = ...
                            parCoordTraj{row3(tempzzminind)}( parCoordTrajPara(row3(tempzzminind),1) : parCoordTrajPara(row3(tempzzminind),1)+parCoordTrajPara(row3(tempzzminind),2)-1, 1:2 );
                        % Remove repeated trajectory segment
                        parCoordTraj{row3(tempzzminind)} = parCoordTraj{row3(tempzzminind)}*nan;
                        
                        % Update varaible "parCoordTrajPara" for trajectory segment length
                        parCoordTrajPara(row(tempi),2) = parCoordTrajPara(row(tempi),2) + parCoordTrajPara(row3(tempzzminind),2) + tempm;
                        parCoordTrajPara(row3(tempzzminind),1:2) = [0,0];
                        
                        %Fillmissing parCoordTraj{row(tempi)}
                        temp = parCoordTraj{row(tempi)}(parCoordTrajPara(row(tempi),1) : sum(parCoordTrajPara(row(tempi),1:2))-1, 1:2);
                        temp_x = fillmissing(temp(:,1),extrapMethod); % fill missing
                        temp_y = fillmissing(temp(:,2),extrapMethod); % fill missing
                        parCoordTraj{row(tempi)}(parCoordTrajPara(row(tempi),1) : sum(parCoordTrajPara(row(tempi),1:2))-1, 1:2) = [temp_x, temp_y];
                        
                    end
                    
                end
                
                
                %%%%% Find all probable trajectory segments in the negative direction %%%%%
                if (sum(parCoordTrajPara(row(tempi),1:2))+tempm<(size(parCoord_prev,1)+1)) && ((parCoordTrajPara(row(tempi),1)-1)-tempm>0)
                    
                    [row2,col2] = find( sum(parCoordTrajPara,2) == parCoordTrajPara(row(tempi),1)-tempm ); % ending point requirement
                    row3 = intersect(row1,row2);
                    tempxx = size(parCoord_prev,1)*(row3-1) + (parCoordTrajPara(row(tempi),1)-1)-tempm; tempxx=tempxx'; tempxx=tempxx(:); % Index in parCoordTrajMat
                    tempyy = parCoordTrajMat(tempxx, 1:2);
                    
                    % hold on; plot(tempyy(:,1),tempyy(:,2),'.');
                    tempzz = sqrt( ( tempyy(:,1)-parCoordTrajCurrPred_x( (parCoordTrajPara(row(tempi),1)-1)-tempm ) ).^2 + ...
                        ( tempyy(:,2)-parCoordTrajCurrPred_y( (parCoordTrajPara(row(tempi),1)-1)-tempm ) ).^2 );
                    
                    [tempzzmin,tempzzminind] = min(tempzz);
                    if tempzzmin < distThres % Find the continuing trajectory segment %JY!!!! threshold distance 3
                        % Merge trajectory segment
                        parCoordTraj{row(tempi)}( parCoordTrajPara(row3(tempzzminind),1) : parCoordTrajPara(row3(tempzzminind),1)+parCoordTrajPara(row3(tempzzminind),2)-1, 1:2 ) = ...
                            parCoordTraj{row3(tempzzminind)}( parCoordTrajPara(row3(tempzzminind),1) : parCoordTrajPara(row3(tempzzminind),1)+parCoordTrajPara(row3(tempzzminind),2)-1, 1:2 );
                        % Remove repeated trajectory segment
                        parCoordTraj{row3(tempzzminind)} = parCoordTraj{row3(tempzzminind)}*nan;
                        % Update varaible "parCoordTrajPara" for both trajectory segment starting point and its length
                        parCoordTrajPara(row(tempi),2) = parCoordTrajPara(row(tempi),2) + parCoordTrajPara(row3(tempzzminind),2) + tempm;
                        parCoordTrajPara(row(tempi),1) = parCoordTrajPara(row3(tempzzminind),1);
                        parCoordTrajPara(row3(tempzzminind),1:2) = [0,0];
                        
                        %Fillmissing parCoordTraj{row(tempi)}
                        temp = parCoordTraj{row(tempi)}(parCoordTrajPara(row(tempi),1) : sum(parCoordTrajPara(row(tempi),1:2))-1, 1:2);
                        temp_x = fillmissing(temp(:,1),extrapMethod); % fill missing
                        temp_y = fillmissing(temp(:,2),extrapMethod); % fill missing
                        parCoordTraj{row(tempi)}(parCoordTrajPara(row(tempi),1) : sum(parCoordTrajPara(row(tempi),1:2))-1, 1:2) = [temp_x, temp_y];
                        
                    end 
                    
                end
                
            end % End of "for each trajectory segment whose "length==tempk" "
            
            
            
        end % End of "for trajectory segments with length [(size(parCoord_prev,1)-1) : -1 : minTrajSegLength]"
        
    end % End of tempm
    
end % End of "for tempMergeTime = 1:5 % Try to merge four times"
close(hbar);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Plot tracked trajectories %%%%%
disp('%%%%% Plot tracked trajectories %%%%%'); fprintf('\n');

figure,
% CN = [36 63 146; 20 76 156; 28 100 175; 9 115 186; 23 128 196; 33 148 210;
%     25 158 218; 19 172 226; 25 186 229; 69 196 221; 118 205 214; 157 216 215;
%     169 220 217; 193 229 224; 216 237 233; 234 246 245]/255;
% CN = CN(end-size(parCoord_prev,1):end, 1:3);


for tempi = 1:size(parCoordTrajPara,1)
    
    wayPoints = parCoordTraj{tempi};
    
    if (length(wayPoints(isnan(wayPoints(:,1))<1,1))+1)<4
        hold on; line(xstep*wayPoints(isnan(wayPoints(:,1))<1,1),xstep*wayPoints(isnan(wayPoints(:,1))<1,2),'linewidth',1.2); view(2); % straight lines
    else
        hold on; fnplt(cscvn(xstep*wayPoints(isnan(wayPoints(:,1))<1,:)'),'',1.2);
    end
     
    %%%%% Codes to plot trajectories with frame-dependent color %%%%% 
    % if sum(1-isnan(wayPoints(:,1)))>1  % Don't show if there is only one point on the trajectory
    %     hold on; plot(xstep*wayPoints(:,1),xstep*wayPoints(:,2),'r.','markersize',5);
    % end
    % 
    % for tempj = 1:size(parCoord_prev,1)-1
    %     hold on; line(xstep*[wayPoints(tempj,1),wayPoints(tempj+1,1)], ...
    %         xstep*[wayPoints(tempj,2),wayPoints(tempj+1,2)],'linewidth',1.2, 'color', CN(tempj,:) );
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
end
 
set(gca,'fontsize',18); view(2); box on; axis equal; axis tight;  
title('Tracked particle trajectory','fontweight','normal');
set(gca,'YDir','reverse'); xlabel('x'); ylabel('y');
axis(xstep*[MPTPara.gridxyROIRange.gridx(1), MPTPara.gridxyROIRange.gridx(2), ...
      MPTPara.gridxyROIRange.gridy(1), MPTPara.gridxyROIRange.gridy(2) ]);
   

  
%%
%%%%% Compute cumulative tracking ratio from emerged trajectories %%%%%
disp('%%%%% Plot tracked cumulative displacements %%%%%'); fprintf('\n');
parCoordTrajMat = cell2mat( parCoordTraj );

[row1,col1] = find(isnan(parCoordTrajMat(1:length(Img):end,1))==0);
trackParCum_ind = row1;
trackParCum_track_ratio = [];

for ImgSeqNum = 2:length(Img)
    [row2,col2] = find(isnan(parCoordTrajMat(ImgSeqNum:length(Img):end,1))==0);
    trackParCum_ind = intersect(row2,trackParCum_ind);
    trackParCum_track_ratio(ImgSeqNum-1) = length(trackParCum_ind) /size(parCoord_prev{1},1);
end

defList = [2:1:length(Img)];
fig=figure; ax=axes; hold on; plot(defList,trackParCum_track_ratio,'bs--','linewidth',1);
adjust_fig(fig,ax,'','',''); box on; title('');
xlabel('Image #'); ylabel('Cumulative tracking ratio');
axis([2,length(Img),0,1]);


%%%%% Plot tracked cumulative displacement field %%%%%
%%%%% Make a video %%%%%
v = VideoWriter('video_2D_inc_cum.avi'); v.FrameRate = 5; open(v); figure,
for ImgSeqNum = 2:length(Img)
    
    parCoordA = parCoordTrajMat(1:length(Img):end,1:2);
    parCoordB = parCoordTrajMat(ImgSeqNum:length(Img):end,1:2);
    parCoordACum = parCoordA(trackParCum_ind,:);
    parCoordBCum = parCoordB(trackParCum_ind,:);
    disp_A2BCum = parCoordBCum - parCoordACum;
    
    % ----- Cone plot grid data: displecement -----
    clf; plotCone2(xstep*parCoordBCum(:,1),xstep*parCoordBCum(:,2),xstep*disp_A2BCum(:,1),xstep*disp_A2BCum(:,2));
    set(gca,'fontsize',18); view(2); box on; axis equal; axis tight; set(gca,'YDir','reverse');
    title(['Tracked cumulative disp (#',num2str(ImgSeqNum),')'],'fontweight','normal');
    xlabel('x'); ylabel('y');
    axis(xstep*[MPTPara.gridxyROIRange.gridx(1), MPTPara.gridxyROIRange.gridx(2), ...
          MPTPara.gridxyROIRange.gridy(1), MPTPara.gridxyROIRange.gridy(2) ]);
     
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Press "Ctrl + C" and modify codes below to plot interpolated displacements and strains on a uniform grid mesh');
pause; 


ImgSeqNum = 2; % Frame #


%%%%% Previously tracked displacement field %%%%%
parCoordA = parCoordTrajMat(1:length(Img):end,1:2);
parCoordB = parCoordTrajMat(ImgSeqNum:length(Img):end,1:2);
parCoordACum = parCoordA(trackParCum_ind,1:2);
parCoordBCum = parCoordB(trackParCum_ind,1:2);
disp_A2BCum = parCoordBCum - parCoordACum;
  
%%%%% Interpolate scatterred data to gridded data %%%%%
sxy = min([round(0.5*MPTPara.f_o_s),20])*[1,1]; % Step size for griddata
smoothness = 1e-3; % Smoothness for regularization; "smoothness=0" means no regularization

[x_Grid_refB,y_Grid_refB,u_Grid_refB]=funScatter2Grid2D(parCoordBCum(:,1),parCoordBCum(:,2),disp_A2BCum(:,1),sxy,smoothness);
[~,~,v_Grid_refB]=funScatter2Grid2D(parCoordBCum(:,1),parCoordBCum(:,2),disp_A2BCum(:,2),sxy,smoothness);

% Apply ROI image mask
[u_Grid_refB, v_Grid_refB] = funRmROIOutside(x_Grid_refB,y_Grid_refB,MPTPara.ImgRefMask,u_Grid_refB,v_Grid_refB);
% Build a displacement vector
uv_Grid_refB_Vector=[u_Grid_refB(:),v_Grid_refB(:)]'; uv_Grid_refB_Vector=uv_Grid_refB_Vector(:);
% Calculate deformation gradient
D_Grid = funDerivativeOp(size(x_Grid_refB,1),size(x_Grid_refB,2),mean(sxy)); % Central finite difference operator
F_Grid_refB_Vector=D_Grid*uv_Grid_refB_Vector; % {F}={D}{U}


%%%%% Cone plot grid data: displecement %%%%%
figure, plotCone2(xstep*x_Grid_refB,xstep*y_Grid_refB,u_Grid_refB*xstep ,v_Grid_refB*xstep );
set(gca,'fontsize',18); view(2); box on; axis equal; axis tight; set(gca,'YDir','reverse');
title('Tracked cumulative displacement','fontweight','normal');
axis(xstep*[MPTPara.gridxyROIRange.gridx(1), MPTPara.gridxyROIRange.gridx(2), ...
      MPTPara.gridxyROIRange.gridy(1), MPTPara.gridxyROIRange.gridy(2) ]);

  
%%%%% Generate a FE-mesh %%%%%
[coordinatesFEM_refB,elementsFEM_refB] = funMeshSetUp(x_Grid_refB*xstep,y_Grid_refB*xstep);

%%%%% Cone plot grid data: displacement %%%%%
Plotdisp_show(uv_Grid_refB_Vector*xstep , coordinatesFEM_refB, elementsFEM_refB,[],'NoEdgeColor');
 
%%%%% Cone plot grid data: infinitesimal strain %%%%%
Plotstrain_show(F_Grid_refB_Vector, coordinatesFEM_refB, elementsFEM_refB,[],'NoEdgeColor',xstep,tstep);
 

 










