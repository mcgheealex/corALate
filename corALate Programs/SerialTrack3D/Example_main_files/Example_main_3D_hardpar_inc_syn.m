%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SerialTrack execute main file
% ===================================================
% Dimension:            3D
% Particle rigidity:    hard 
% Tracking mode:        incremental
% Syn or Exp:           syn
% Deformation mode:     rigid body motions including translations & rotations
%
% ===================================================
% Author: Jin Yang, Ph.D.
% Email: jyang526@wisc.edu -or-  aldicdvc@gmail.com 
% Date: 02/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
close all; clear all; clc; clearvars -global
disp('************************************************');
disp('*** Welcome to SerialTrack Particle Tracking ***');
disp('************************************************');
addpath( './function/','./src/','./Scatter2Grid3D' ); 

 
%% User defined parameters %%%%%

%%%%% Problem dimension and units %%%%%
MPTPara.DIM = 3;    % problem dimension
MPTPara.xstep = 1;  % unit: um/px
MPTPara.tstep = 1;  % unit: us

%%%%% Code mode %%%%%
MPTPara.mode = 'inc'; % {'inc': incremental mode; 
                      %  'cum': cumulative mode}

%%%%% Particle rigidity %%%%%
MPTPara.parType = 'hard'; % {'hard': hard particle; 
                          %  'soft': soft particle}

disp('************************************************');
disp(['Dimention: ',num2str(MPTPara.DIM)]);
disp(['Tracking mode: ',MPTPara.mode]);
disp(['Particle type: ',MPTPara.parType]);
disp('************************************************'); fprintf('\n');

%%%%% SerialTrack path %%%%%
SerialTrackPath = 'D:\MATLAB\SerialTrack3D'; % TODO: modify the path

%%%%% Volumetric image path %%%%%
% fileNameAll = 'VM_3D*.mat';
% fileFolder = '';

%%%%% Synthetic cases %%%%%
DefType = 'stretch';        % {'translation','stretch','simpleshear','rotation'}
SeedingDensityType = 2;     % {1,2,3,4} % particle seeding density
fileNameAll = 'vol_*.mat';  % file name(s)
% ----- file folder name -----
if strcmp(DefType,'translation')==1
    fileFolder = ['./imgFolder/img_syn_hardpar/img_trans_hardpar_sd',num2str(SeedingDensityType)];
elseif strcmp(DefType,'rotation')==1
    fileFolder = ['./imgFolder/img_syn_hardpar/img_rot_hardpar_sd',num2str(SeedingDensityType)];
elseif strcmp(DefType,'stretch')==1
    fileFolder = ['./imgFolder/img_syn_hardpar/img_str_hardpar_sd',num2str(SeedingDensityType)];
elseif strcmp(DefType,'simpleshear')==1
    fileFolder = ['./imgFolder/img_syn_hardpar/img_simshear_hardpar_sd',num2str(SeedingDensityType)];
else
end
 
%%%%% Bead detection method %%%%%
BeadPara.detectionMethod = 1; % {1-deconvolution code; 2-regionprops}

%%%%% Image binary mask file %%%%%
im_roi_mask_file_path = ''; % TODO: leave it as empty if there is no mask file


%%%%% Particle detection parameters %%%%%
%%%%% Bead Parameter %%%%%
BeadPara.thres = 0.5;           % Threshold for detecting particles
BeadPara.beadSize = 3;          % Estimated radius of a single particle
BeadPara.minSize = 4;           % Minimum radius of a single particle
BeadPara.maxSize = 100;         % Maximum radius of a single particle
BeadPara.winSize = [5,5,5];     % By default
BeadPara.dccd = [1,1,1];        % By default
BeadPara.abc = [1,1,1];         % By default
BeadPara.forloop = 1;           % By default
BeadPara.randNoise = 1e-7;      % By default
BeadPara.PSF = [];              % PSF function; Example: PSF = fspecial('disk', BeadPara.beadSize-1 ); % Disk blur
BeadPara.distMissing = 5;       % Distance threshold to check whether particle has a match or not 
BeadPara.color = 'white';       % By default


%% SerialTrack particle tracking

%%%%% Multiple particle tracking (MPT) Parameter %%%%%
MPTPara.f_o_s = 60;              % Size of search field: max(|u|,|v|,|w|)
MPTPara.n_neighborsMax = 25;     % Max # of neighboring particles
MPTPara.n_neighborsMin = 1;      % Min # of neighboring particles
MPTPara.gbSolver = 2;            % Global step solver: 1-moving least square fitting; 2-global regularization; 3-ADMM iterations
MPTPara.smoothness = 1e-1;       % Coefficient of regularization
MPTPara.outlrThres = 5;          % Threshold for removing outliers in MPT
MPTPara.maxIterNum = 20;         % Max ADMM iteration number
MPTPara.iterStopThres = 1e-3;    % ADMM iteration stopping threshold
MPTPara.strain_n_neighbors = 20; % # of neighboring particles used in strain gauge
MPTPara.strain_f_o_s = 60;       % Size of virtual strain gauge
MPTPara.usePrevResults = 0;      % Whether use previous results or not: 0-no; 1-yes;  


%%%% Postprocessing: merge trajectory segments %%%%%
distThres = 1;              % distance threshold to connect split trajectory segments
extrapMethod = 'pchip';     % extrapolation scheme to connect split trajectory segments
                            % suggestion: 'nearest' for Brownian motion                          
minTrajSegLength = 10;      % the minimum length of trajectory segment that will be extrapolate 
maxGapTrajSeqLength = 0;    % the max frame# gap between connected trajectory segments

 
%%%%% Execute SerialTrack particle tracking %%%%%
if strcmp(MPTPara.mode,'inc')==1
    run_Serial_MPT_3D_hardpar_inc;
elseif strcmp(MPTPara.mode,'cum')==1
    run_Serial_MPT_3D_hardpar_cum;    
end
 

