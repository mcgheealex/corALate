%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SerialTrack execute main file
% ===================================================
% Dimension:            2D
% Particle rigidity:    soft 
% Tracking mode:        cumulative
% Syn or Exp:           exp
% Deformation mode:     foam uniaxial compression
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
addpath( './function/','./src/');  

 
%% user defined parameters %%%%%

%%%%% Problem dimension and units %%%%%
MPTPara.DIM = 2;    % problem dimension
MPTPara.xstep = 1;  % unit: um/px
MPTPara.tstep = 1;  % unit: us

%%%%% Code mode %%%%%
MPTPara.mode = 'cum'; % {'inc': incremental mode; 
                      %  'cum': cumulative mode; 
                      %  'dbf': double frame}

MPTPara.parType = 'soft'; % {'hard': hard particle; 
                          %  'soft': soft particle}

disp('************************************************');
disp(['Dimention: ',num2str(MPTPara.DIM)]);
disp(['Tracking mode: ',MPTPara.mode]);
disp(['Particle type: ',MPTPara.parType]);
disp('************************************************'); fprintf('\n');

%%%%% Image binary mask file %%%%%
MaskFileLoadingMode = 0; % {0: No mask file
                         %  1: loading only one mask file for all frames; 
                         %  2: loading one mask file for each single frame;
                         %  3: loading matlab mat file for all frames}

if MaskFileLoadingMode == 3
    im_roi_mask_file_path = '.\img_par2track_lung\im_roi.mat'; % TODO: Path for the loaded MATLAB mat file as the mask file
else
    im_roi_mask_file_path = '';  
end

%%%%% Particle detection parameters %%%%%
%%%%%%%%%% Elastomeric open-cell foam %%%%%%%%%%%%%
%%%%% Bead Parameter %%%%%
BeadPara.thres = 0.3;           % Threshold for detecting particles
BeadPara.beadSize = 0;         % Estimated radius of a single particle
BeadPara.minSize = 3;           % Minimum radius of a single particle
BeadPara.maxSize = 200;         % Maximum radius of a single particle
BeadPara.winSize = [5, 5];      % By default
BeadPara.dccd = [1,1];          % By default
BeadPara.abc = [1,1];           % By default
BeadPara.forloop = 1;           % By default
BeadPara.randNoise = 1e-7;      % By default
BeadPara.PSF = [];              % PSF function; Example: PSF = fspecial('disk', BeadPara.beadSize-1 ); % Disk blur
BeadPara.distMissing = 10;      % Distance threshold to check whether particle has a match or not
BeadPara.color = 'black';       % Bead color: 'white' -or- 'black' 


%% SerialTrack particle tracking

%%%%% Multiple particle tracking (MPT) Parameter %%%%%
MPTPara.f_o_s = 30;              % Size of search field: max(|u|,|v|)
MPTPara.n_neighborsMax = 25;     % Max # of neighboring particles
MPTPara.n_neighborsMin = 1;      % Min # of neighboring particles
MPTPara.locSolver = 1;           % Local solver: 1-topology-based feature; 2-histogram-based feature first and then topology-based feature;
MPTPara.gbSolver = 3;            % Global step solver: 1-moving least square fitting; 2-global regularization; 3-ADMM iterations
MPTPara.smoothness = 1e-1;       % Coefficient of regularization
MPTPara.outlrThres = 2;          % Threshold for removing outliers in TPT
MPTPara.maxIterNum = 20;         % Max ADMM iteration number
MPTPara.iterStopThres = 1e-3;    % ADMM iteration stopping threshold
MPTPara.strain_n_neighbors = 20; % # of neighboring particles used in strain gauge
MPTPara.strain_f_o_s = 60;       % Size of virtual strain gauge
MPTPara.usePrevResults = 1;      % Whether use previous results or not: 0-no; 1-yes;  


%%%% Postprocessing: merge trajectory segments %%%%%
% distThres = 1; % distance threshold to connect split trajectory segments
extrapMethod = 'pchip';  % extrapolation scheme to connect split trajectory segments
                         % suggestion: 'nearest' for Brownian motion                          
minTrajSegLength = 20;   % the minimum length of trajectory segment that will be extrapolate 
% maxGapTrajSeqLength = 0; % the max frame# gap between connected trajectory segments



%%%%% Execute SerialTrack particle tracking %%%%%
if strcmp(MPTPara.mode,'inc')==1
    disp('Not implemented yet.');
    % run_Serial_MPT_2D_softpar_inc;
elseif strcmp(MPTPara.mode,'cum')==1
    run_Serial_MPT_2D_softpar_cum;    
elseif strcmp(MPTPara.mode,'dbf')==1
    disp('Not implemented yet.');
    % run_Serial_MPT_2D_hardpar_dbf;
end

 




