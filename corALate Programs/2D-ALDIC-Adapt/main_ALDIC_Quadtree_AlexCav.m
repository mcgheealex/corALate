% ---------------------------------------------
% Augmented Lagrangian Digital Image Correlation (ALDIC_Quadtree)
% using an adaptive quadtree mesh, which was automatically generated
% based on the DIC raw images
% 
% Author: Jin Yang, PhD @Caltech
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Date: 2015.04,06,07; 2016.03,04; 2020.11
% ---------------------------------------------

%% Section 1: Clear MATLAB environment & mex set up Spline interpolation  
close all; clear; clc; clearvars -global
fprintf('------------ Section 1 Start ------------ \n')
setenv('MW_MINGW64_LOC','C:\TDM-GCC-64');
try mex -O ba_interp2.cpp; catch; end  % mex set up ba_interp2.cpp script
% [Comment]: If this line reports error but it works before, 
% Change line 15 to: "try mex -O ba_interp2.cpp; catch; end"
addpath("./func",'./src','./plotFiles','./func_quadtree','./func_quadtree/refinement','./plotFiles/export_fig-d966721/'); 
addpath('./Images_Quadtree_demo/Images_Sample12'); % TODO: addpath("./YOUR IMAGE FOLDER"); 
addpath('./func/rbfinterp');
fprintf('------------ Section 1 Done ------------ \n \n')


%%%%%%%%% TODO: %%%%%%%%%%%%

% ====== First Alex_cav_dic dataset ======
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_alex_cav\R_data.mat');
% disp_ur_eq = 5.06; % unit: px

% ====== 150um_Top_2Mfps_1 =========
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_70121_DICpaper_150um_Top_2Mfps_1\R_data.mat')
% Rnew([2:6]) = [];
% disp_ur_eq = 2.62; % unit: px

% ====== 100um_Top_2Mfps_1 =========
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_70121_DICpaper_100um_Top_2Mfps_1\R_data.mat')
% Rnew([2:6]) = [];
% disp_ur_eq = 3.35 % !!! previous wrong data 5.10; % unit: px


% ===== GE10: 19_53_04 =======
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_19_53_04\R_data.mat');
% disp_ur_eq = 2.19;  % unit: px

% ===== GE10: 20_02_48 ======
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_20_02_48\R_data.mat');
% disp_ur_eq = 0; % unit: px

% ===== GE4shot1 ======
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_92121_4percent_shot1\R_data.mat');
% Rnew([2:7]) = []; 
% disp_ur_eq = 2.3156; % unit: px

% ===== GE6shot3: 12_58_29 ======
load('R_data.mat');
Rnew([2:7]) = []; 
disp_ur_eq = 3.0439; % unit: px

% ===== GE10shot1 ======
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_92121_10percent_shot1\R_data.mat');
% Rnew([2:7]) = []; 
% disp_ur_eq = 2.254; % unit: px

% ===== GE14shot1 ======
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_92121_14percent_shot1\R_data.mat');
% Rnew([2:7]) = []; 
% disp_ur_eq = 3.1402; % unit: px


%% Section 2: Load DIC parameters and set up DIC parameters 
fprintf('------------ Section 2 Start ------------ \n')
% ====== Read images ====== 
[file_name,Img,DICpara] = ReadImageQuadtree; % Load DIC raw images

%%%%%%%% !!!mask START %%%%%%%
% [DICpara] = ReadImageMask(DICpara); % Load and define an image mask (which is a binary image)
[mask_file_name,ImgMask] = ReadImageMasks; % Load DIC image mask files
%%%%%%%% !!!mask END %%%%%%%

% %%%%%% Uncomment lines below to change the DIC computing region (ROI) manually %%%%%%
% DICpara.gridxROIRange = [gridxROIRange1,gridxROIRange2]; DICpara.gridyROIRange = [Val1, Val2];
% E.g., gridxROIRange = [224,918]; gridyROIRange = [787,1162];

% ====== Normalize images: fNormalized = (f-f_avg)/(f_std) ======
[ImgNormalized,DICpara.gridxyROIRange] = funNormalizeImg(Img,DICpara.gridxyROIRange); 
 
% ====== Initialize variable storage ======
ResultDisp = cell(length(ImgNormalized)-1,1);    ResultDefGrad = cell(length(ImgNormalized)-1,1);
ResultStrain = cell(length(ImgNormalized)-1,1);  ResultStress = cell(length(ImgNormalized)-1,1);
ResultFEMeshEachFrame = cell(length(ImgNormalized)-1,1); % To store FE-mesh for each frame: needs future improvment to make it more efficient.
ResultFEMesh = cell(ceil((length(ImgNormalized)-1)/DICpara.ImgSeqIncUnit),1); % For incremental DIC mode
DICpara.SizeOfFFTSearchRegion = [8,8];
fprintf('------------ Section 2 Done ------------ \n \n')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To solve each frame in an image sequence
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ImgSeqNum =  2 :  length(ImgNormalized) 

    close all; 
    
    if ImgSeqNum < 7
        DICpara.NewFFTSearch = 1;
    else
        DICpara.NewFFTSearch = 0;
    end

    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);
    
    % ====== Compute image gradients ======
    %%%%%%%% !!!mask START %%%%%%%
    % fNormalized = ImgNormalized{length(ImgNormalized)}; % Load the first referece image 
    % DICpara.ImgRefMask = ImgMask{length(ImgMask)};
    
    fNormalizedMask = double( ImgMask{1} ) ;
    fNormalized = ImgNormalized{1} .* fNormalizedMask; % Load the first referece image 
    % fNormalized = fNormalized .* double(fNormalizedMask);
    %%%%%%%% !!!mask END %%%%%%%
    Df = funImgGradient(fNormalized,fNormalized,fNormalizedMask); % Finite difference to compute image grayscale gradients;
    %%%%%%%% !!!mask END %%%%%%%
    
    %%%%%%%% !!!mask START %%%%%%%
    % gNormalized = ImgNormalized{length(ImgNormalized)+1-ImgSeqNum} ; % Load current deformed image frame 
    % gNormalizedMask = double(ImgMask{length(ImgNormalized)+1-ImgSeqNum});
    gNormalized = ImgNormalized{ ImgSeqNum} ; % Load current deformed image frame 
    gNormalizedMask = double(ImgMask{ ImgSeqNum});
    gNormalized = gNormalized .* gNormalizedMask ;
    %%%%%%%% !!!mask END %%%%%%%
    
    DICpara.ImgRefMask = fNormalizedMask;
    
    figure, subplot(2,2,1); imshow(fNormalized'); title('fNormalized')
    subplot(2,2,2); imshow(gNormalized'); title('gNormalized')
    subplot(2,2,3); imshow(fNormalizedMask'); title('f mask')
    subplot(2,2,4); imshow(gNormalizedMask'); title('g mask')
    
    
    %% Section 3: Compute an initial guess of the unknown displacement field
    fprintf('\n'); fprintf('------------ Section 3 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to find or update an initial guess of the unknown displacements.
    % The key idea is to either to use a new FFT-based cross correlation peak fitting,
    % or use the results from the last frame as the new initial guess for the next frame;
    % Particularly in the incremental mode DIC, the reference image can also be updated, e.g.,
    % " fNormalized = ImgNormalized{ImgSeqNum-mod(ImgSeqNum-1,ImgSeqIncUnit)}; "
    %
    % DICpara.NewFFTSearch = 0; % If you want to apply the FFT-based cross correlation to 
    % compute the initial guess for each frame, please make sure that "DICpara.NewFFTSearch = 0". 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ImgSeqNum == 2 || DICpara.NewFFTSearch == 1 % Apply FFT-based cross correlation to compute the initial guess 
        
        % JY!!!dddddd
        % DICpara.winsize = 20 ; 
        
        % ====== FFT-based cross correlation ======
        %JY!!! [DICpara,x0temp,y0temp,u,v,cc]= IntegerSearch(fNormalized,gNormalized,file_name,DICpara);
        [DICpara,x0temp_g,y0temp_g,u_g,v_g,cc]= IntegerSearch(gNormalized,fNormalized,file_name,DICpara);
        u_f = -u_g; v_f = -v_g; x0temp_f = x0temp_g+u_g; y0temp_f = y0temp_g+v_g;
        
        %%%%% Interpolate to f %%%%%
        xnodes = max([1+0.5*DICpara.winsize,DICpara.gridxyROIRange.gridx(1)])  ...
            : DICpara.winstepsize : min([size(fNormalized,1)-0.5*DICpara.winsize-1,DICpara.gridxyROIRange.gridx(2)]);
        ynodes = max([1+0.5*DICpara.winsize,DICpara.gridxyROIRange.gridy(1)])  ...
            : DICpara.winstepsize : min([size(fNormalized,2)-0.5*DICpara.winsize-1,DICpara.gridxyROIRange.gridy(2)]);
         
        % ceil(min(x0temp_f(:))) : DICpara.winstepsize : floor(max(x0temp_f(:)));
        % ynodes = % ceil(min(y0temp_f(:))) : DICpara.winstepsize : floor(max(y0temp_f(:)));
        [x0temp,y0temp] = ndgrid(xnodes,ynodes);
        u_f_NotNanInd = find(~isnan(u_f(:)));
        
        %u = griddata(x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd),u_f(u_f_NotNanInd),x0temp,y0temp,'linear');
        %v = griddata(x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd),v_f(u_f_NotNanInd),x0temp,y0temp,'linear');
%         F_u = scatteredInterpolant( x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd),u_f(u_f_NotNanInd),'natural','nearest' );
%         F_v = scatteredInterpolant( x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd),v_f(u_f_NotNanInd),'natural','nearest' );
%         u = F_u( x0temp,y0temp ); v = F_v( x0temp,y0temp );
%         
        
%     load('E:\Jin\kb305_JinYang\2D_ALDIC_v4.1\img_alex_cav\New folder (2)\Med Bubble Top\R_data.mat')
% 
%     bubble_y_img = (CircleFitPar(ImgSeqNum,1))*1; %DICpara.um2px;
%     bubble_x_img = CircleFitPar(ImgSeqNum,2)*1; % JY!!!
% 
%     delta_r = mean(Rnew(end-30:end));
%     theta_list = linspace(-pi,pi,201);
%     x_at_delta_r = bubble_x_img + cos(theta_list) * delta_r;
%     y_at_delta_r = bubble_y_img - sin(theta_list) * delta_r;
% 
%     disp_r_at_delta_r = Rnew(ImgSeqNum)-delta_r;
%     disp_x_at_delta_r = cos(theta_list) * disp_r_at_delta_r;
%     disp_y_at_delta_r = - sin(theta_list) * disp_r_at_delta_r;

    
    op1 = rbfcreate( [x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd)]',[u_f(u_f_NotNanInd)]','RBFFunction', 'thinplate'); rbfcheck(op1);
    u = rbfinterp([x0temp(:),y0temp(:)]', op1 );
    
%     op1_001 = rbfcreate( [x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd);  x_at_delta_r(:),y_at_delta_r(:)]', ...
%         [u_f(u_f_NotNanInd); disp_x_at_delta_r(:)]','RBFFunction', 'thinplate'); rbfcheck(op1);
%     u = rbfinterp([x0temp(:),y0temp(:)]', op1_001 );
%     
%     op1_002 = rbfcreate( [x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd);  x_at_delta_r(:),y_at_delta_r(:)]', ...
%         [v_f(u_f_NotNanInd); disp_y_at_delta_r(:)]','RBFFunction', 'thinplate'); rbfcheck(op1);
%     v = rbfinterp([x0temp(:),y0temp(:)]', op1_002 );


        op2 = rbfcreate( [x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd)]',[v_f(u_f_NotNanInd)]','RBFFunction', 'thinplate'); rbfcheck(op2);
        v = rbfinterp([x0temp(:),y0temp(:)]', op2 );

% op1 = rbfcreate_mod_img_mask( [x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd)]',[u_f(u_f_NotNanInd)]',fNormalizedMask,'RBFFunction', 'thinplate'); rbfcheck(op1);
% u = rbfinterp([x0temp(:),y0temp(:)]', op1 );
% op2 = rbfcreate_mod_img_mask( [x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd)]',[v_f(u_f_NotNanInd)]',fNormalizedMask,'RBFFunction', 'thinplate'); rbfcheck(op2);
% v = rbfinterp([x0temp(:),y0temp(:)]', op2 );

 
% figure, scatter3([x0temp(:) ], [y0temp(:)],  [u(:)],  '.' );
%  hold on, scatter3([x0temp_f(u_f_NotNanInd) ], [y0temp_f(u_f_NotNanInd)],  [u_f(u_f_NotNanInd)],  'k.' );
% hold on,  scatter3([x0temp(:) ], [y0temp(:)],  [u_000(:)],  'r.' );
% hold on; scatter3(x_at_delta_r(:),y_at_delta_r(:),disp_y_at_delta_r(:),'r.');


        u = regularizeNd([x0temp(:),y0temp(:)],u(:),{xnodes',ynodes'},1e-3);
        v = regularizeNd([x0temp(:),y0temp(:)],v(:),{xnodes',ynodes'},1e-3);

        % ====== DIC uniform FE-mesh set up ======
        [DICmesh] = MeshSetUp(x0temp,y0temp,DICpara); % clear x0temp y0temp;
        % ====== Initial Value ======
        U0 = Init(u,v,cc.max,DICmesh.x0,DICmesh.y0,0); % [Temp code:] PlotuvInit;
        
 
        
        for tempi = 1:size(u,1)
            for tempj = 1:size(u,2)
                try
                    if ~fNormalizedMask(x0temp(tempi,tempj),y0temp(tempi,tempj)) % || ...
                          %   ~gNormalizedMask( floor(x0temp(tempi,tempj)+u(tempi,tempj)), floor(y0temp(tempi,tempj)+v(tempi,tempj)) )
                         
                        U0(2*(tempj+(tempi-1)*(size(u,2)))) = nan;
                        U0(2*(tempj+(tempi-1)*(size(u,2)))-1) = nan;

                    end
                catch
                end
            end
        end
        
        
        % ====== Deal with incremental mode ======
        fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
        if DICpara.ImgSeqIncUnit == 1, fNormalizedNewIndex = fNormalizedNewIndex-1; end
        ResultFEMesh{1+floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info
            struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ...
            'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );
         
%         % ====== Generate a quadtree mesh considering sample's complex geometry ======
        DICmesh.elementMinSize = 2; % min element size in the refined quadtree mesh
%         %%%%%%%% !!!mask START %%%%%%%
%         %%%% Update mask file of fNormalized %%%%%
%         xnodes_Img = 1:1:size(gNormalizedMask,1); % min(DICmesh.coordinatesFEM(:,1)):1:max(DICmesh.coordinatesFEM(:,1));
%         ynodes_Img = 1:1:size(gNormalizedMask,2); % min(DICmesh.coordinatesFEM(:,2)):1:max(DICmesh.coordinatesFEM(:,2));
%         [XX_Img,YY_Img] = ndgrid(xnodes_Img, ynodes_Img);
%         U0NotNanInd = find(~isnan(U0(1:2:end)));
%         F_u = scatteredInterpolant( DICmesh.coordinatesFEM(U0NotNanInd,1),DICmesh.coordinatesFEM(U0NotNanInd,2),U0(2*U0NotNanInd-1),'linear','linear' );
%         F_v = scatteredInterpolant( DICmesh.coordinatesFEM(U0NotNanInd,1),DICmesh.coordinatesFEM(U0NotNanInd,2),U0(2*U0NotNanInd),'linear','linear' );
%         u220 = F_u( XX_Img, YY_Img ); v220 = F_v( XX_Img, YY_Img );
%         
% %         u220 = regularizeNd([DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2)],U0(1:2:end),{xnodes_Img',ynodes_Img'},1e-3);
% %         v220 = regularizeNd([DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2)],U0(2:2:end),{xnodes_Img',ynodes_Img'},1e-3);
%             
%         %u220 = griddata(DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),U0(1:2:end),XX,YY,'linear');
%         %v220 = griddata(DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),U0(2:2:end),XX,YY,'linear');
%         u22 = u220 + XX_Img; v22 = v220 + YY_Img;
%         temp = ba_interp2(gNormalizedMask, v22, u22, 'cubic'); temp(isnan(temp(:))) = 0;
%         DICpara.ImgRefMask = logical(temp); Df.ImgRefMask = DICpara.ImgRefMask;
%         % fNormalized = ImgNormalized{length(ImgNormalized)}; % Load the first referece image
%         fNormalized = ImgNormalized{1}; % Load the first referece image
%         % fNormalized = fNormalized .* double(DICpara.ImgRefMask);
        
        GenerateQuadtreeMesh; % Generate a quadtree mesh
        
        %%%%% Update search region %%%%%
        DICpara.SizeOfFFTSearchRegion = [ ceil( max( [max(3+abs(U0(1:2:end))), 3] ) ), ...
                                          ceil( max( [max(3+abs(U0(2:2:end))), 3] ) ) ];
        %%%%%%%% !!!mask END %%%%%%%
        
        %%%%%%%% !!!mask START %%%%%%%
%         DICmesh.markCoordHoleEdge = [];
%         DICmesh.dirichlet = DICmesh.markCoordHoleEdge;
        %%%%%%%% !!!mask END %%%%%%%
         
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) == 0 % To update ref image in incremental mode
        fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
        if DICpara.ImgSeqIncUnit == 1,  fNormalizedNewIndex = fNormalizedNewIndex-1; end
        fNormalized = ImgNormalized{fNormalizedNewIndex}; % Update reference
        [DICpara,DICmesh] = ReadImageRefUpdate(file_name,ImgSeqNum,ResultDisp{ImgSeqNum-2}.U,DICpara,DICmesh); % Update reference image if needed;
        U0 = zeros(2*size(DICmesh.coordinatesFEM,1),1); % [Temporary code: " PlotuvInit; "] 
        ResultFEMesh{1+floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info
            struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ...
            'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else % Use the solved results from the last frame as the new initial guess
        if ImgSeqNum < 7 %|| ImgSeqNum==12 % Import previous U for ImgSeqNum [2,6] 
            U0 = ResultDisp{ImgSeqNum-2}.U;
             
        else % When ImgSeqNum > 6: POD predicts next disp U0 from previous results of (ImgSeqNum+[-5:1:-1])
            nTime = 5; 
            
            xnodes = DICpara.winstepsize+DICpara.gridxyROIRange.gridx(1) : DICpara.winstepsize : -DICpara.winstepsize+DICpara.gridxyROIRange.gridx(2);
            ynodes = DICpara.winstepsize+DICpara.gridxyROIRange.gridy(1) : DICpara.winstepsize : -DICpara.winstepsize+DICpara.gridxyROIRange.gridy(2);
            [x0temp,y0temp] = ndgrid(xnodes,ynodes);
            np = size(x0temp,1)*size(x0temp,2);
            T_data_u = zeros(nTime,np); T_data_v = zeros(nTime,np);
            for tempi = 1:nTime
                
                tempx = ResultFEMeshEachFrame{ImgSeqNum-(2+nTime)+tempi, 1}.coordinatesFEM(:,1);
                tempy = ResultFEMeshEachFrame{ImgSeqNum-(2+nTime)+tempi, 1}.coordinatesFEM(:,2);
                tempu = ResultDisp{ImgSeqNum-(2+nTime)+tempi, 1}.U(1:2:end);
                tempv = ResultDisp{ImgSeqNum-(2+nTime)+tempi, 1}.U(2:2:end);
                
%                 F_u = scatteredInterpolant( tempx,tempy,tempu,'natural','nearest' );
%                 F_v = scatteredInterpolant( tempx,tempy,tempv,'natural','nearest' );
%                 unodes = F_u( x0temp,y0temp );
%                 vnodes = F_v( x0temp,y0temp );


                op1 = rbfcreate( [tempx(:),tempy(:)]',[tempu(:)]','RBFFunction', 'thinplate'); rbfcheck(op1);
                unodes = rbfinterp([x0temp(:),y0temp(:)]', op1 );
                op2 = rbfcreate( [tempx(:),tempy(:)]',[tempv(:)]','RBFFunction', 'thinplate'); rbfcheck(op2);
                vnodes = rbfinterp([x0temp(:),y0temp(:)]', op2 );       
 
 
                T_data_u(tempi,:) = unodes(:);
                T_data_v(tempi,:) = vnodes(:);
                
            end
            nB = 3;
            t_train = [ImgSeqNum-1-nTime:ImgSeqNum-2]'; 
            t_pre = [ImgSeqNum-1]';
            [u_pred,~,~,~] = funPOR_GPR(T_data_u,t_train,t_pre,nB);
            [v_pred,~,~,~] = funPOR_GPR(T_data_v,t_train,t_pre,nB);
            trained_u = reshape(u_pred(1,:),size(x0temp)); 
            trained_v = reshape(v_pred(1,:),size(x0temp));
            
            % ====== DIC uniform FE-mesh set up ======
            [DICmesh] = MeshSetUp(x0temp,y0temp,DICpara); % clear x0temp y0temp;
            % ====== Initial Value ======
            U0 = Init(trained_u,trained_v,[],DICmesh.x0,DICmesh.y0,0); % [Temp code:] PlotuvInit;
        
            %%%%%% U0 = [trained_u(:),trained_v(:)]'; U0 = U0(:);
            
            for tempi = 1:size(trained_u,1)
                for tempj = 1:size(trained_u,2)
                    try
                        if ~fNormalizedMask(x0temp(tempi,tempj),y0temp(tempi,tempj)) % || ...
                             %  ~gNormalizedMask( floor(x0temp(tempi,tempj)+u(tempi,tempj)), floor(y0temp(tempi,tempj)+v(tempi,tempj)) )
                            
                            U0(2*(tempj+(tempi-1)*(size(u,2)))) = nan;
                            U0(2*(tempj+(tempi-1)*(size(u,2)))-1) = nan;
                            
                        end
                    catch
                    end
                end
            end
            
            
             
%             np = length(ResultDisp{ImgSeqNum-2}.U)/2; % "nTime" value 5 is an empirical value, can be changed.
%             T_data_u = zeros(nTime,np); T_data_v = zeros(nTime,np); 
%             for tempi = 1:nTime
%                 T_data_u(tempi,:) = ResultDisp{ImgSeqNum-(2+nTime)+tempi, 1}.U(1:2:np*2)';
%                 T_data_v(tempi,:) = ResultDisp{ImgSeqNum-(2+nTime)+tempi, 1}.U(2:2:np*2)';
%             end
%             nB = 3; t_train = [ImgSeqNum-1-nTime:ImgSeqNum-2]'; t_pre = [ImgSeqNum-1]';
%             [u_pred,~,~,~] = funPOR_GPR(T_data_u,t_train,t_pre,nB);
%             [v_pred,~,~,~] = funPOR_GPR(T_data_v,t_train,t_pre,nB);
%             tempu = u_pred(1,:); tempv = v_pred(1,:);
%             U0 = [tempu(:),tempv(:)]'; U0 = U0(:);
             
        end
        
        
        % ====== Generate a quadtree mesh considering sample's complex geometry ======
        DICmesh.elementMinSize = 2; % min element size in the refined quadtree mesh
        %%%%%%%% !!!mask START %%%%%%%
        %%%% Update mask file of fNormalized %%%%%
        xnodes_Img = 1:1:size(gNormalizedMask,1); % min(DICmesh.coordinatesFEM(:,1)):1:max(DICmesh.coordinatesFEM(:,1));
        ynodes_Img = 1:1:size(gNormalizedMask,2); % min(DICmesh.coordinatesFEM(:,2)):1:max(DICmesh.coordinatesFEM(:,2));
        [XX_Img,YY_Img] = ndgrid(xnodes_Img, ynodes_Img);
        U0NotNanInd = find(~isnan(U0(1:2:end)));
        %F_u = scatteredInterpolant( DICmesh.coordinatesFEM(U0NotNanInd,1),DICmesh.coordinatesFEM(U0NotNanInd,2),U0(2*U0NotNanInd-1),'linear','linear' );
        %F_v = scatteredInterpolant( DICmesh.coordinatesFEM(U0NotNanInd,1),DICmesh.coordinatesFEM(U0NotNanInd,2),U0(2*U0NotNanInd),'linear','linear' );
        
        u220 = regularizeNd([DICmesh.coordinatesFEM(U0NotNanInd,1),DICmesh.coordinatesFEM(U0NotNanInd,2)], ...
            U0(2*U0NotNanInd-1),{xnodes_Img',ynodes_Img'},1e-1);
        v220 = regularizeNd([DICmesh.coordinatesFEM(U0NotNanInd,1),DICmesh.coordinatesFEM(U0NotNanInd,2)], ...
            U0(2*U0NotNanInd ),{xnodes_Img',ynodes_Img'},1e-1);
        
        for tempij = 1:size(u220,1)*size(u220,2)
            try
                if ~fNormalizedMask(XX_Img(tempij),YY_Img(tempij)) %|| ...
                      % ~gNormalizedMask( floor(XX_Img(tempij)+u220(tempij)), floor(YY_Img(tempij)+v220(tempij)) )
                    u220(tempij) = nan;
                    v220(tempij) = nan;
                    
                end
            catch
            end
            
        end
        
        
%         %u220 = griddata(DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),U0(1:2:end),XX,YY,'linear');
%         %v220 = griddata(DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),U0(2:2:end),XX,YY,'linear');
%         u22 = u220 + XX_Img; v22 = v220 + YY_Img;
%         temp = ba_interp2(gNormalizedMask, v22, u22, 'cubic'); temp(isnan(temp(:))) = 0;
%         
%         DICpara.ImgRefMask = logical(temp); Df.ImgRefMask = DICpara.ImgRefMask;
%         % fNormalized = ImgNormalized{length(ImgNormalized)}; % Load the first referece image
%         fNormalized = ImgNormalized{1}; % Load the first referece image
%         % fNormalized = fNormalized .* double(DICpara.ImgRefMask);
        
        GenerateQuadtreeMesh; % Generate a quadtree mesh
        
        
    end
    
    % ====== Compute f(X)-g(x+u) ======
    % PlotImgDiff(x0,y0,u,v,fNormalized,gNormalized);
    ResultFEMeshEachFrame{ImgSeqNum-1} = struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM,'markCoordHoleEdge',DICmesh.markCoordHoleEdge );
    fprintf('------------ Section 3 Done ------------ \n \n')

 
    %% Section 4: ALDIC Subproblem 1 -or- Local ICGN Subset DIC
    fprintf('------------ Section 4 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve the first local step in ALDIC: Subproblem 1
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % JY!!!
    % DICpara.winsize = 40; 
     
    % ====== ALStep 1 Subproblem1: Local Subset DIC ======
    mu=0; beta=0; tol=1e-2; ALSolveStep=1; ALSub1Time=zeros(6,1); ALSub2Time=zeros(6,1); 
    ConvItPerEle=zeros(size(DICmesh.coordinatesFEM,1),6); 
    ALSub1BadPtNum=zeros(6,1);
    
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ Start Local DIC IC-GN iteration ------
    [USubpb1,FSubpb1,HtempPar,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp,markCoordHoleStrain] = ...
        LocalICGNQuadtree(U0,DICmesh.coordinatesFEM,Df,fNormalized,gNormalized,DICpara,'GaussNewton',tol);
    
    ALSub1Time(ALSolveStep) = ALSub1Timetemp; 
    ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; 
    ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp; 
    toc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------  Manually find some bad points from Local Subset ICGN step ------
    % Comment these lines below if you don't have local bad points
    % %%%%% Comment START %%%%%%
    % USubpb1(2*DICmesh.markCoordHoleEdge-1:2*DICmesh.markCoordHoleEdge) = nan;
    % FSubpb1(4*DICmesh.markCoordHoleEdge-3:4*DICmesh.markCoordHoleEdge) = nan;
%     [USubpb1,FSubpb1] = funRemoveOutliersQuadtree(DICmesh,DICpara,USubpb1,FSubpb1);
%     disp('--- Remove bad points done ---')
    % %%%%% Comment END %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %%%%%%%% !!!mask START %%%%%%%
%         DICmesh.markCoordHoleEdge = [];
%         DICmesh.dirichlet = DICmesh.markCoordHoleEdge;
        %%%%%%%% !!!mask END %%%%%%%
        % figure, scatter3(DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),USubpb1(2:2:end))
        
        
%         load('E:\Jin\kb305_JinYang\2D_ALDIC_v4.1\img_alex_cav\New folder (2)\Med Bubble Top\R_data.mat')
%         
%         bubble_y = (251-CircleFitPar(ImgSeqNum+1,1))*DICpara.um2px;
%         bubble_x = CircleFitPar(ImgSeqNum+1,2)*DICpara.um2px;
%         
%         
%         r = sqrt( (coordinatesFEMWorldDef(:,1)-bubble_x).^2 + (coordinatesFEMWorldDef(:,2)-bubble_y).^2  );
%         theta = atan2( -(coordinatesFEMWorldDef(:,2)-bubble_y), coordinatesFEMWorldDef(:,1)-bubble_x);
%         
%         disp_r = cos(theta).*disp_u - sin(theta).*disp_v;
%         disp_t = sin(theta).*disp_u + cos(theta).*disp_v;
       
coordinatesFEM = DICmesh.coordinatesFEM;
U = USubpb1; F = FSubpb1; 
nanindex = find(isnan(U(1:2:end))==1); 
notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindex);
nanindexF = find(isnan(F(1:4:end))==1);
notnanindexF = setdiff([1:1:size(coordinatesFEM,1)],nanindexF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),U(2*notnanindex-1),'natural','linear');
% U1 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
% Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),U(2*notnanindex),'natural','linear');
% U2 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
% Ftemp = scatteredInterpolant(coordinatesFEM(notnanindexF,1),coordinatesFEM(notnanindexF,2),F(4*notnanindexF-3),'linear','linear');
% F11 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
% Ftemp = scatteredInterpolant(coordinatesFEM(notnanindexF,1),coordinatesFEM(notnanindexF,2),F(4*notnanindexF-2),'linear','linear');
% F21 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
% Ftemp = scatteredInterpolant(coordinatesFEM(notnanindexF,1),coordinatesFEM(notnanindexF,2),F(4*notnanindexF-1),'linear','linear');
% F12 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
% Ftemp = scatteredInterpolant(coordinatesFEM(notnanindexF,1),coordinatesFEM(notnanindexF,2),F(4*notnanindexF-0),'linear','linear');
% F22 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
% 
% U_scattInterp = [U1(:),U2(:)]'; U_scattInterp = U_scattInterp(:);
% F_scattInterp = [F11(:),F21(:),F12(:),F22(:)]'; F_scattInterp = F_scattInterp(:);
         



 
bubble_y_img = (251-mean(CircleFitPar(end-20:end,1)))*1;
bubble_x_img = mean(CircleFitPar(end-20:end,2))*1;
 
delta_r = mean(Rnew(end-30:end)) - disp_ur_eq; % Req = 5.06 px is computed from comparing the first and the last frames
theta_list = linspace(-pi,pi,201);
x_at_delta_r = bubble_x_img + cos(theta_list) * delta_r;
y_at_delta_r = bubble_y_img - sin(theta_list) * delta_r;

disp_r_at_delta_r = Rnew(ImgSeqNum) - delta_r;
disp_x_at_delta_r = cos(theta_list) * disp_r_at_delta_r;
disp_y_at_delta_r = - sin(theta_list) * disp_r_at_delta_r;
%  
% % op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
% %     U(2*notnanindex-1)','RBFFunction', 'multiquadric', 'RBFConstant', 2); rbfcheck(op);
% % fi1 = rbfinterp([coordinatesFEM(:,1:2)]', op );
% % op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
% %     U(2*notnanindex)','RBFFunction', 'multiquadric', 'RBFConstant', 2); rbfcheck(op);
% % fi2 = rbfinterp([coordinatesFEM(:,1:2)]', op );         
% % U_rbf_multiquadratic = [fi1(:),fi2(:)]';    U_rbf_multiquadratic = U_rbf_multiquadratic(:);
% 
% % figure, scatter3([coordinatesFEM(notnanindex,1 ); [x_at_delta_r(:) ]], ...
% %   [coordinatesFEM(notnanindex,2 ); [y_at_delta_r(:) ]], ...
% %   [U(2*notnanindex-1); disp_x_at_delta_r(:)],  '.' );
% 
op1 =rbfcreate( [coordinatesFEM(notnanindex,1:2); [x_at_delta_r(:), y_at_delta_r(:)]]', ...
    [U(2*notnanindex-1);  disp_x_at_delta_r(:)]','RBFFunction', 'thinplate'); 
rbfcheck(op1);
fi1 = rbfinterp([coordinatesFEM(:,1:2)]', op1 );
op2 =rbfcreate( [coordinatesFEM(notnanindex,1:2); [x_at_delta_r(:), y_at_delta_r(:)]]', ...
    [U(2*notnanindex); disp_y_at_delta_r(:)]','RBFFunction', 'thinplate'); rbfcheck(op2);
fi2 = rbfinterp([coordinatesFEM(:,1:2)]', op2 );         
U_rbf_thinplate = [fi1(:),fi2(:)]';    U_rbf_thinplate = U_rbf_thinplate(:);

% figure, scatter3( coordinatesFEM( : ,1 ) ,  coordinatesFEM(: ,2 )  ,  U_rbf_thinplate(1:2:end)  ,  '.' );

Convtemp = [ConvItPerEletemp, ConvItPerEletemp]'; Convtemp = Convtemp(:);
 Plotdisp_show(Convtemp,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');

% Plotdisp_show(U0,DICmesh.coordinatesFEM,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
% USubpb1ICGN = USubpb1; Plotdisp_show(USubpb1ICGN,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor'); pause;
% Plotdisp_show(U_scattInterp,DICmesh.coordinatesFEM,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
% Plotdisp_show(U_rbf_multiquadratic,DICmesh.coordinatesFEM,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
%  Plotdisp_show(U_rbf_thinplate,DICmesh.coordinatesFEM,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
 

% op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
%     F(4*notnanindex-3)','RBFFunction', 'multiquadric', 'RBFConstant', 2); rbfcheck(op);
% fi11 = rbfinterp([coordinatesFEM(:,1:2)]', op );
% op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
%     F(4*notnanindex-2)','RBFFunction', 'multiquadric', 'RBFConstant', 2); rbfcheck(op);
% fi21 = rbfinterp([coordinatesFEM(:,1:2)]', op );   
% op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
%     F(4*notnanindex-1)','RBFFunction', 'multiquadric', 'RBFConstant', 2); rbfcheck(op);
% fi12 = rbfinterp([coordinatesFEM(:,1:2)]', op );
% op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
%     F(4*notnanindex-0)','RBFFunction', 'multiquadric', 'RBFConstant', 2); rbfcheck(op);
% fi22 = rbfinterp([coordinatesFEM(:,1:2)]', op );
% 
% F_rbf_multiquadratic = [fi11(:),fi21(:),fi12(:),fi22(:)]';  
% F_rbf_multiquadratic = F_rbf_multiquadratic(:);

op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
    F(4*notnanindex-3)','RBFFunction', 'thinplate'); rbfcheck(op);
fi11 = rbfinterp([coordinatesFEM(:,1:2)]', op );
op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
    F(4*notnanindex-2)','RBFFunction', 'thinplate'); rbfcheck(op);
fi21 = rbfinterp([coordinatesFEM(:,1:2)]', op );    
op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
    F(4*notnanindex-1)','RBFFunction', 'thinplate'); rbfcheck(op);
fi12 = rbfinterp([coordinatesFEM(:,1:2)]', op ); 
op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
    F(4*notnanindex-0)','RBFFunction', 'thinplate'); rbfcheck(op);
fi22 = rbfinterp([coordinatesFEM(:,1:2)]', op ); 
 

F_rbf_thinplate = [fi11(:),fi21(:),fi12(:),fi22(:)]';    F_rbf_thinplate = F_rbf_thinplate(:);


% Plotstrain_show(FSubpb1,DICmesh.coordinatesFEM,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
% Plotstrain_show(F_scattInterp,DICmesh.coordinatesFEM,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
% Plotstrain_show(F_rbf_multiquadratic,DICmesh.coordinatesFEM,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
% Plotstrain_show(F_rbf_thinplate,DICmesh.coordinatesFEM,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');


    USubpb1 = U_rbf_thinplate;
    FSubpb1 = F_rbf_thinplate;

 

    % ------ Plot ------
    USubpb1World = USubpb1; USubpb1World(2:2:end) = -USubpb1(2:2:end); FSubpb1World = FSubpb1; 
%     close all;
Plotdisp_show(USubpb1World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
%     Plotstrain_show(FSubpb1World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
    save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
    fprintf('------------ Section 4 Done ------------ \n \n')

    %%%%% 
% Convtemp = [ConvItPerEletemp, ConvItPerEletemp]'; Convtemp = Convtemp(:);
%  
% PlotMesh_show(Convtemp,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
%     
% pause;

   %%
%     %%%%% START Mask %%%%%
%     % ====== Generate a quadtree mesh considering sample's complex geometry ======
%     DICmesh.elementMinSize = 1; % min element size in the refined quadtree mesh
%     %%%%%%%% !!!mask START %%%%%%%
%     %%%% Update mask file of fNormalized %%%%%
%     xnodes_Img = 1:1:size(gNormalizedMask,1); % min(DICmesh.coordinatesFEM(:,1)):1:max(DICmesh.coordinatesFEM(:,1));
%     ynodes_Img = 1:1:size(gNormalizedMask,2); % min(DICmesh.coordinatesFEM(:,2)):1:max(DICmesh.coordinatesFEM(:,2));
%     [XX_Img,YY_Img] = ndgrid(xnodes_Img, ynodes_Img);
%     USubpb1NotNanInd = find(~isnan(USubpb1(1:2:end)));
%     F_u = scatteredInterpolant( DICmesh.coordinatesFEM(USubpb1NotNanInd,1),DICmesh.coordinatesFEM(USubpb1NotNanInd,2),USubpb1(2*USubpb1NotNanInd-1),'linear','linear' );
%     F_v = scatteredInterpolant( DICmesh.coordinatesFEM(USubpb1NotNanInd,1),DICmesh.coordinatesFEM(USubpb1NotNanInd,2),USubpb1(2*USubpb1NotNanInd),'linear','linear' );
%     u220 = F_u( XX_Img, YY_Img ); v220 = F_v( XX_Img, YY_Img );
%     
%     % u220 = regularizeNd([DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2)],U0(1:2:end),{xnodes_Img',ynodes_Img'},1e-3);
%     % v220 = regularizeNd([DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2)],U0(2:2:end),{xnodes_Img',ynodes_Img'},1e-3);
%     
%     %u220 = griddata(DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),U0(1:2:end),XX,YY,'linear');
%     %v220 = griddata(DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),U0(2:2:end),XX,YY,'linear');
%     u22 = u220 + XX_Img; v22 = v220 + YY_Img;
%     temp = ba_interp2(gNormalizedMask, v22, u22, 'cubic'); temp(isnan(temp(:))) = 0;
%     DICpara.ImgRefMask = logical(temp); Df.ImgRefMask = DICpara.ImgRefMask;
    %%%%% Mask END %%%%%
    
        
   
    %% Section 5: Subproblem 2 -- solve the global compatible displacement field
    fprintf('------------ Section 5 Start ------------ \n'); tic;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve the global step in ALDIC Subproblem 2
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ======= ALStep 1 Subproblem 2: Global constraint =======
    % ------ Smooth displacements for a better F ------
    DICpara.DispFilterSize=0; DICpara.DispFilterStd=0; DICpara.StrainFilterSize=0; DICpara.StrainFilterStd=0; LevelNo=1;
    DICpara.DispSmoothness = 0; DICpara.StrainSmoothness = 1e-4;
    if DICpara.DispSmoothness>1e-6, USubpb1 = funSmoothDispQuadtree(USubpb1,DICmesh,DICpara); end
    if DICpara.StrainSmoothness>1e-6, FSubpb1 = funSmoothStrainQuadtree(FSubpb1,DICmesh,DICpara); end
    
	% ====== Define penalty parameter ======
    mu = 1e-3; udual = 0*FSubpb1; vdual = 0*USubpb1; 
    betaList = [1e-3,1e-2,1e-1]*mean(DICpara.winstepsize).^2.*mu; % Tune beta in the betaList 
    Err1 = zeros(length(betaList),1); Err2 = Err1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****']);
    DICpara.GaussPtOrder = 2; alpha = 0;  % No regularization added
    % ====== Solver using finite element method ======
    if ImgSeqNum == 2
        for tempk = 1:length(betaList)
            beta = betaList(tempk); display(['Try #',num2str(tempk),' beta = ',num2str(beta)]);
            GaussPtOrder=3; alpha=0; [USubpb2] = Subpb2Quadtree(DICmesh,DICpara.GaussPtOrder,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,mean(DICpara.winstepsize),0);
            [FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2,DICpara.GaussPtOrder,0);

            Err1(tempk) = norm(USubpb1-USubpb2,2);
            Err2(tempk) = norm(FSubpb1-FSubpb2,2);
        end
        
        Err1Norm = (Err1-mean(Err1))/std(Err1); % figure, plot(Err1Norm);
        Err2Norm = (Err2-mean(Err2))/std(Err2); % figure, plot(Err2Norm);
        ErrSum = Err1Norm+Err2Norm; % figure, plot(ErrSum); title('Tune the best \beta in the subproblem 2'); 
        [~,indexOfbeta] = min(ErrSum); 
     
        try % Tune the best beta by a quadratic polynomial 0fitting
            [fitobj] = fit(log10(betaList(indexOfbeta-1:1:indexOfbeta+1))',ErrSum(indexOfbeta-1:1:indexOfbeta+1),'poly2');
            p = coeffvalues(fitobj); beta = 10^(-p(2)/2/p(1));
        catch, beta = betaList(indexOfbeta);
        end
        display(['Best beta = ',num2str(beta)]);
    else 
        try beta = DICpara.beta;
        catch, beta = 1e-3*mean(DICpara.winstepsize).^2.*mu;
        end
    end
      
    % Using the optimal beta to solve the ALDIC Subproblem 2 again
    if abs(beta-betaList(end))>abs(eps)
        [USubpb2] = Subpb2Quadtree(DICmesh,DICpara.GaussPtOrder,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,mean(DICpara.winstepsize),0);
        [FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2,DICpara.GaussPtOrder,0);
        ALSub2Time(ALSolveStep) = toc; toc
    end
    
    % ------- Smooth strain field --------
    if DICpara.DispSmoothness>1e-6, USubpb2 = funSmoothDispQuadtree(USubpb2,DICmesh,DICpara); end
    % ------- Don't smooth strain fields near the boundary --------
    for tempk=0:3, FSubpb2(4*DICmesh.markCoordHoleEdge-tempk) = FSubpb1(4*DICmesh.markCoordHoleEdge-tempk); end
    if DICpara.StrainSmoothness>1e-6, FSubpb2 = funSmoothStrainQuadtree(0.1*FSubpb2+0.9*FSubpb1,DICmesh,DICpara); end
    for tempk=0:1, USubpb2(2*markCoordHoleStrain-tempk) = USubpb1(2*markCoordHoleStrain-tempk); end
    for tempk=0:3, FSubpb2(4*markCoordHoleStrain-tempk) = FSubpb1(4*markCoordHoleStrain-tempk); end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ------- Save data ------
    save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2');

    % ------ Plot ------
    USubpb2World = USubpb2; USubpb2World(2:2:end) = -USubpb2(2:2:end); FSubpb2World = FSubpb2;  
%     close all; Plotdisp_show(USubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
%     Plotstrain_show(FSubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');

    % ======= Update dual variables =======
    udual = FSubpb2 - FSubpb1; vdual = USubpb2 - USubpb1;
	save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');
    fprintf('------------ Section 5 Done ------------ \n \n')
 

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Section 6: ADMM iterations
    fprintf('------------ Section 6 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is the ADMM iteration, where both Subproblems 1 & 2 are solved iteratively.
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ==================== ADMM AL Loop ==========================
    ALSolveStep = 1; tol2 = 1e-2; UpdateY = 1e4;  
    HPar = cell(21,1); for tempj = 1:21, HPar{tempj} = HtempPar(:,tempj); end

    while (ALSolveStep < 3)
        ALSolveStep = ALSolveStep + 1;  % Update using the last step
        
        
%         Ftemp1 = FSubpb2(1:2:end); Ftemp2 = FSubpb2(2:2:end);
%         [DFtemp1,~,~] = funGlobalNodalStrainQuadtree(DICmesh,Ftemp1,DICpara.GaussPtOrder,0);
%         [DFtemp2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,Ftemp2,DICpara.GaussPtOrder,0);
%         
%         winsize_x_ub1 = abs(2*FSubpb2(1:4:end)./DFtemp1(1:4:end));
%         winsize_x_ub2 = abs(2*FSubpb2(3:4:end)./DFtemp1(3:4:end));
%         winsize_y_ub1 = abs(2*FSubpb2(1:4:end)./DFtemp1(2:4:end));
%         winsize_y_ub2 = abs(2*FSubpb2(3:4:end)./DFtemp1(4:4:end));
%         
%         winsize_x_ub3 = abs(2*FSubpb2(2:4:end)./DFtemp2(1:4:end));
%         winsize_x_ub4 = abs(2*FSubpb2(4:4:end)./DFtemp2(3:4:end));
%         winsize_y_ub3 = abs(2*FSubpb2(2:4:end)./DFtemp2(2:4:end));
%         winsize_y_ub4 = abs(2*FSubpb2(4:4:end)./DFtemp2(4:4:end));
%         
%         winsize_x_ub = round(min([winsize_x_ub1,winsize_x_ub2,winsize_x_ub3,winsize_x_ub4,DICpara.winsize*ones(length(winsize_x_ub1),1)],[],2));
%         winsize_x_List = max([winsize_x_ub, 10*ones(length(winsize_x_ub1),1)],[],2);
%         winsize_y_ub = round(min([winsize_y_ub1,winsize_y_ub2,winsize_y_ub3,winsize_y_ub4,DICpara.winsize*ones(length(winsize_y_ub1),1)],[],2));
%         winsize_y_List = max([winsize_y_ub, 10*ones(length(winsize_y_ub1),1)],[],2);
%         winsize_List = 2*ceil([winsize_x_List,winsize_y_List]/2);
         winsize_List = DICpara.winsize*ones(size(DICmesh.coordinatesFEM,1),2);
        DICpara.winsize_List = winsize_List;
        

        %%%%%%%%%%%%%%%%%%%%%%% Subproblem 1 %%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****']);
        tic; [USubpb1,~,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp] = Subpb1Quadtree(...
                                            USubpb2,FSubpb2,udual,vdual,DICmesh.coordinatesFEM,...
                                            Df,fNormalized,gNormalized,mu,beta,HPar,ALSolveStep,DICpara,'GaussNewton',tol);
        FSubpb1 = FSubpb2; toc 
        % for tempk=0:1, USubpb1(2*markCoordHoleStrain-tempk) = USubpb2(2*markCoordHoleStrain-tempk); end
        ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ------  Manually find some bad points from Local Subset ICGN step ------
        % disp('--- Start to manually remove bad points --- \n')
        % disp('    Comment codes here if you do not have bad local points \n')
        % %%%%% Comment START %%%%%
        %  [USubpb1,FSubpb1] = funRemoveOutliersQuadtree(DICmesh,DICpara,USubpb1,FSubpb1);
        %  disp('--- Remove bad points done ---')
        % %%%%% Comment END %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ============== Subproblem 2 ==============
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****'])
        tic; [USubpb2] = Subpb2Quadtree(DICmesh,DICpara.GaussPtOrder,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,mean(DICpara.winstepsize),0);
		[FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2,DICpara.GaussPtOrder,0);
        ALSub2Time(ALSolveStep) = toc; toc
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ------- Smooth strain field --------
        if DICpara.DispSmoothness>1e-6, USubpb2 = funSmoothDispQuadtree(USubpb2,DICmesh,DICpara); end
        % ------- Don't change strain fields near the boundary --------
        for tempk=0:3, FSubpb2(4*DICmesh.markCoordHoleEdge-tempk) = FSubpb1(4*DICmesh.markCoordHoleEdge-tempk); end
        if DICpara.StrainSmoothness>1e-6, FSubpb2 = funSmoothStrainQuadtree(0.1*FSubpb2+0.9*FSubpb1,DICmesh,DICpara); end
        for tempk=0:1, USubpb2(2*markCoordHoleStrain-tempk) = USubpb1(2*markCoordHoleStrain-tempk); end
        for tempk=0:3, FSubpb2(4*markCoordHoleStrain-tempk) = FSubpb1(4*markCoordHoleStrain-tempk); end
         
		save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2');
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute norm of UpdateY
        USubpb2_Old = load(['Subpb2_step',num2str(ALSolveStep-1)],'USubpb2');
        USubpb2_New = load(['Subpb2_step',num2str(ALSolveStep)],'USubpb2');
        USubpb1_Old = load(['Subpb1_step',num2str(ALSolveStep-1)],'USubpb1');
        USubpb1_New = load(['Subpb1_step',num2str(ALSolveStep)],'USubpb1');
        if (mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) ~= 0 && (ImgSeqNum>2)) || (ImgSeqNum < DICpara.ImgSeqIncUnit)
            UpdateY = norm((USubpb2_Old.USubpb2 - USubpb2_New.USubpb2), 2)/sqrt(size(USubpb2_Old.USubpb2,1));
            try
                UpdateY2 = norm((USubpb1_Old.USubpb1 - USubpb1_New.USubpb1), 2)/sqrt(size(USubpb1_Old.USubpb1,1));
            catch
            end
        end
        try
            disp(['Update local step  = ',num2str(UpdateY2)]);
            disp(['Update global step = ',num2str(UpdateY)]);
        catch
        end
        fprintf('*********************************** \n \n');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update dual variables------------------------------
        udual = FSubpb2 - FSubpb1; vdual = USubpb2 - USubpb1; 
		 
        save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');
        try
        if UpdateY < tol2 || UpdateY2 < tol2
            break
        end
        catch
        end
         
    end
    fprintf('------------ Section 6 Done ------------ \n \n')
 
    % Save data
    ResultDisp{ImgSeqNum-1}.U = full(USubpb2);
    ResultDisp{ImgSeqNum-1}.ALSub1BadPtNum = ALSub1BadPtNum;
    ResultDefGrad{ImgSeqNum-1}.F = full(FSubpb2);  


%     % Save data
%     ResultDisp{ImgSeqNum-1}.U = full(USubpb1);
%     ResultDisp{ImgSeqNum-1}.ALSub1BadPtNum = ALSub1BadPtNum;
%     ResultDefGrad{ImgSeqNum-1}.F = full(FSubpb1); 
    

end    


%% ------ Plot ------
USubpb2World = USubpb2; USubpb2World(2:2:end) = -USubpb2(2:2:end); FSubpb2World = FSubpb2; 
close all; Plotdisp_show(USubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
Plotstrain_show(FSubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');

% winsize_xy = [winsize_x_List, winsize_y_List]'; winsize_xy = winsize_xy(:);
% Plotdisp_show(winsize_xy,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
%     
%         
        
% ------ Save results ------
% Find img name and save all the results 
[~,imgname,imgext] = fileparts(file_name{1,end});
results_name = ['results_',imgname,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultFEMesh','ResultFEMeshEachFrame','ALSub1Time','ALSub2Time','ALSolveStep');


%% Section 7: Check convergence
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to check convergence of ADMM
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('------------ Section 7 Start ------------ \n')
% ====== Check convergence ======
fprintf('***** Check convergence ***** \n');
ALSolveStep1 = min(6,ALSolveStep);
disp('==== uhat^(k) - u^(k) ====');
for ALSolveStep = 1:ALSolveStep1
    USubpb2 = load(['Subpb2_step',num2str(ALSolveStep )],'USubpb2');
    USubpb1 = load(['Subpb1_step',num2str(ALSolveStep )],'USubpb1');
    UpdateY = norm((USubpb2.USubpb2 - USubpb1.USubpb1), 2)/sqrt(length(USubpb2.USubpb2));
    disp(num2str(UpdateY));
end
disp('==== Fhat^(k) - F^(k) ====');
for ALSolveStep = 1:ALSolveStep1
    FSubpb1 = load(['Subpb1_step',num2str(ALSolveStep )],'FSubpb1');
    FSubpb2 = load(['Subpb2_step',num2str(ALSolveStep )],'FSubpb2');
    UpdateF = norm((FSubpb1.FSubpb1 - FSubpb2.FSubpb2), 2)/sqrt(length(FSubpb1.FSubpb1));
    disp(num2str(UpdateF));
end
disp('==== uhat^(k) - uhat^(k-1) ====');
for ALSolveStep = 2:ALSolveStep1
    USubpb2_Old = load(['Subpb2_step',num2str(ALSolveStep-1)],'USubpb2');
    USubpb2_New = load(['Subpb2_step',num2str(ALSolveStep)],'USubpb2');
    UpdateY = norm((USubpb2_Old.USubpb2 - USubpb2_New.USubpb2), 2)/sqrt(length(USubpb2.USubpb2));
    disp(num2str(UpdateY));
end
disp('==== udual^(k) - udual^(k-1) ====');
for ALSolveStep = 2:ALSolveStep1
    uvdual_Old = load(['uvdual_step',num2str(ALSolveStep-1)],'udual');
    uvdual_New = load(['uvdual_step',num2str(ALSolveStep)],'udual');
    UpdateW = norm((uvdual_Old.udual - uvdual_New.udual), 2)/sqrt(length(uvdual_Old.udual));
    disp(num2str(UpdateW));
end
disp('==== vdual^(k) - vdual^(k-1) ====');
for ALSolveStep = 2:ALSolveStep1
    uvdual_Old = load(['uvdual_step',num2str(ALSolveStep-1)],'vdual');
    uvdual_New = load(['uvdual_step',num2str(ALSolveStep)],'vdual');
    Updatev = norm((uvdual_Old.vdual - uvdual_New.vdual), 2)/sqrt(length(uvdual_Old.vdual));
    disp(num2str(Updatev));
end
fprintf('------------ Section 7 Done ------------ \n \n')

% ------ Delete temp files ------
%%%%% Comment START %%%%%
% Uncomment these lines to delete temporary files
%   for tempi = 1:ALSolveStep
%       file_name_Subpb1 = ['Subpb1_step',num2str(tempi),'.mat'];
%       file_name_Subpb2 = ['Subpb2_step',num2str(tempi),'.mat'];
%       file_name_dual = ['uvdual_step',num2str(tempi),'.mat'];
%       delete(file_name_Subpb1); delete(file_name_Subpb2); delete(file_name_dual);
%   end
%%%%% Comment END %%%%%

% ------ clear temp variables ------
clear a ALSub1BadPtNum ALSub1Timetemp atemp b btemp cc ConvItPerEletemp hbar Hbar 
clear coordinatesFEMQuadtree elementsFEMQuadtree 


%% Section 8: Compute strains
fprintf('------------ Section 8 Start ------------ \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to compute strain fields and plot disp and strain results
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Convert units from pixels to the physical world ------
DICpara.um2px = funParaInput('ConvertUnit');
% ------ Smooth displacements ------
DICpara.DoYouWantToSmoothOnceMore = 1; % No need to smooth disp fields
DICpara.smoothness = funParaInput('RegularizationSmoothness'); % Regularization to smooth strain fields           
% ------ Choose strain computation method ------
DICpara.MethodToComputeStrain = 2; %funParaInput('StrainMethodOp'); 
if DICpara.MethodToComputeStrain == 2 % Compute strain method II: Use Plane Fitting method
    prompt = 'What is your half window size (unit: px): ';
    Rad = input(prompt);      
end
% ------ Choose strain type (infinitesimal, Eulerian, Green-Lagrangian) ------
DICpara.StrainType = funParaInput('StrainType');
% ------ Choose image to plot (first only, second and next images) ------
if length(ImgNormalized)==2, DICpara.Image2PlotResults = funParaInput('Image2PlotResults');
else DICpara.Image2PlotResults = 1; % Plot over current, deformed image by default
end
% ------ Save fig format ------
DICpara.MethodToSaveFig = funParaInput('SaveFigFormat');
% ------ Choose overlay image transparency ------
DICpara.OrigDICImgTransparency = 1; 
if DICpara.MethodToSaveFig == 1  
    DICpara.OrigDICImgTransparency = funParaInput('OrigDICImgTransparency');         
end

%%
 

%% ====== Start main part ======
% This section is to calculate strain fields based on the transformed
% cumulative displacements: [F] = [D][U]
% [F] = [..., F11_nodei, F21_nodei, F12_nodei, F22_nodei, ...]';
% [u] = [..., U1_nodei, U2_nodei, ...]';
% [D]: finite difference/finite element operator to compute first derivatives

% v = VideoWriter('video_cav_DIC.mp4');
% v.FrameRate = 10;
% open(v);

ResultCavDIC = cell(length(ImgNormalized)-1,1);
for ImgSeqNum =  [  2: length(ImgNormalized) ] 
    
    close all; disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);
     
    %%%%% Load deformed image %%%%%%%
    gNormalizedMask = double( ImgMask{ImgSeqNum} ); % Load the mask file of current deformed frame
    gNormalized = ImgNormalized{ImgSeqNum} .* gNormalizedMask ; % Load current deformed frame 
    Dg = funImgGradient(gNormalized,gNormalized,gNormalizedMask); % Finite difference to compute image grayscale gradients;
    
    fNormalizedMask = double( ImgMask{1} ); % 
    DICpara.ImgRefMask = fNormalizedMask;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ====== Old version ======
    % fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
    % if DICpara.ImgSeqIncUnit > 1
    %     FEMeshIndLast = floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit);
    % elseif DICpara.ImgSeqIncUnit == 1
    %     FEMeshIndLast = floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)-1;
    % end
    % FEMeshInd = FEMeshIndLast + 1;
    % 
    % if FEMeshInd == 1
    %     USubpb2 = ResultDisp{ImgSeqNum-1}.U; %+ ResultDisp{10}.U + ResultDisp{20}.U;
    %     coordinatesFEM = ResultFEMesh{1}.coordinatesFEM; 
    %     elementsFEM = ResultFEMesh{1}.elementsFEM;
    %     if (ImgSeqNum-1 == 1) || (DICpara.ImgSeqIncROIUpdateOrNot==1), UFEMesh = 0*USubpb2; end
    % else
    %     USubpb2 = ResultDisp{ImgSeqNum-1}.U;
    %     if mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) == 0
    %         coordinatesFEM = ResultFEMesh{FEMeshInd}.coordinatesFEM;
    %         elementsFEM = ResultFEMesh{FEMeshInd}.elementsFEM;
    %         coordinatesFEMLast = ResultFEMesh{FEMeshIndLast}.coordinatesFEM;
    %         UFEMeshLast = ResultDisp{ImgSeqNum-2}.U + UFEMesh;
    %         xq = coordinatesFEM(:,1); yq = coordinatesFEM(:,2);
    %         UFEMesh = 0*USubpb2;
    %         UFEMesh(1:2:end) = griddata(coordinatesFEMLast(:,1),coordinatesFEMLast(:,2),UFEMeshLast(1:2:end),xq,yq,'v4');
    %         UFEMesh(2:2:end) = griddata(coordinatesFEMLast(:,1),coordinatesFEMLast(:,2),UFEMeshLast(2:2:end),xq,yq,'v4');
    %     end
    %     USubpb2 = USubpb2 + UFEMesh;
    % end
    %
    % FSubpb2 = ResultDefGrad{ImgSeqNum-1}.F;
    % coordinatesFEM = ResultFEMeshEachFrame{ImgSeqNum-1}.coordinatesFEM;
    % elementsFEM = ResultFEMeshEachFrame{ImgSeqNum-1}.elementsFEM;
    
    USubpb2 = ResultDisp{ImgSeqNum-1}.U;
    coordinatesFEM = ResultFEMeshEachFrame{ImgSeqNum-1}.coordinatesFEM;
    elementsFEM = ResultFEMeshEachFrame{ImgSeqNum-1}.elementsFEM;
 
    try markCoordHoleEdge = ResultFEMeshEachFrame{ImgSeqNum-1}.markCoordHoleEdge; catch; end
    DICmesh.coordinatesFEM = coordinatesFEM;
    DICmesh.elementsFEM = elementsFEM;
    coordinatesFEMWorld = DICpara.um2px*[coordinatesFEM(:,1),size(ImgNormalized{1},2)+1-coordinatesFEM(:,2)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ------ Plotting and Compute Strain-------
    if size(USubpb2,1) == 1
        ULocal = USubpb2_New.USubpb2; FLocal = FSubpb2.FSubpb2; 
    else
        ULocal = USubpb2; FLocal = FSubpb2;
    end
    UWorld = DICpara.um2px*ULocal; UWorld(2:2:end) = -UWorld(2:2:end);   % close all; Plotuv(UWorld,x0,y0World);
     
    % ------ Smooth displacements ------
    %prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
    %DoYouWantToSmoothOnceMore = input(prompt); 
    SmoothTimes = 0;
    try
        while DICpara.DoYouWantToSmoothOnceMore == 0 && SmoothTimes < 3
            ULocal = funSmoothDispQuadtree(ULocal,DICmesh,DICpara);
            %close all; Plotuv(ULocal,x0,y0); %DICpara.DoYouWantToSmoothOnceMore = input(prompt);
            SmoothTimes = SmoothTimes + 1;
        end
    catch
    end
    
    % ----- Compute strain field ------
    ComputeStrainQuadtree;  
    
    % ------ Smooth strain fields ------
    DICpara.DoYouWantToSmoothOnceMore = 0; % Smooth strain fields if necessary
    %prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
    %DoYouWantToSmoothOnceMore = input(prompt); 
    SmoothTimes = 0;
    try
        while DICpara.DoYouWantToSmoothOnceMore == 0 && SmoothTimes < 3
            FStraintemp = funSmoothStrainQuadtree(FStraintemp,DICmesh,DICpara);
            %close all; Plotuv(ULocal,x0,y0); %DICpara.DoYouWantToSmoothOnceMore = input(prompt);
            SmoothTimes = SmoothTimes + 1;
        end
    catch
    end
    FStrainWorld = FStraintemp; FStrainWorld(2:4:end) = -FStrainWorld(2:4:end); FStrainWorld(3:4:end) = -FStrainWorld(3:4:end); 
     
    
    % ------ Plot disp and strain ------
    if DICpara.OrigDICImgTransparency == 1
        Plotdisp_show(UWorld,coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'NoEdgeColor');
        [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min,strain_maxshear,strain_vonMises] = ...
                   Plotstrain0Quadtree(FStrainWorld,coordinatesFEMWorld,elementsFEM(:,1:4),DICpara);
    
    else % Plot over raw DIC images
        if DICpara.Image2PlotResults == 0 % Plot over the first image; "file_name{1,1}" corresponds to the first image	
            PlotdispQuadtree(UWorld,coordinatesFEMWorld,elementsFEM(:,1:4),file_name{1,1},DICpara);
            [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
                strain_maxshear,strain_vonMises] = PlotstrainQuadtree(UWorld,FStrainWorld, ...
                coordinatesFEMWorld,elementsFEM(:,1:4),file_name{1,1},DICpara);

        else % Plot over second or next deformed images
            
            %%%%%% Old codes: without applying mask files %%%%%%
%             PlotdispQuadtree(UWorld,coordinatesFEMWorld,elementsFEM(:,1:4),...
%                file_name{1,ImgSeqNum},DICpara);
%             
%             [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
%                strain_maxshear,strain_vonMises] = PlotstrainQuadtree(UWorld,FStraintemp, ...
%                coordinatesFEMWorld,elementsFEM(:,1:4),file_name{1,ImgSeqNum},DICpara);

            %%%%%% New codes: applying mask files %%%%%%
%             PlotdispQuadtreeMasks(UWorld,coordinatesFEMWorld,elementsFEM(:,1:4),...
%                 file_name{1, ImgSeqNum}, ...
%                 ImgMask{ ImgSeqNum },DICpara);
%             
%             [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
%                strain_maxshear,strain_vonMises] = PlotstrainQuadtreeMasks(UWorld,FStraintemp, ...
%                coordinatesFEMWorld,elementsFEM(:,1:4),file_name{1, ImgSeqNum}, ...
%                ImgMask{ ImgSeqNum },DICpara);
           

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO 
            PlotdispQuadtreePolarMasks(UWorld,coordinatesFEMWorld,elementsFEM(:,1:4),...
                file_name{1, ImgSeqNum}, ...
                ImgMask{ ImgSeqNum },DICpara,ImgSeqNum);
            
            [strain_err,strain_ert,strain_ett,strain_principal_max,strain_principal_min, ...
               strain_maxshear,strain_vonMises] = PlotstrainQuadtreePolarMasks(UWorld,FStrainWorld, ...
               coordinatesFEMWorld,elementsFEM(:,1:4),file_name{1, ImgSeqNum}, ...
               ImgMask{ ImgSeqNum },DICpara, [] ,ImgSeqNum);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO           
           

 %%%%%%%%% cav-dic part %%%%%%%%%% TODO 
            disp_u = UWorld(1:2:end); disp_v = UWorld(2:2:end);
           coordinatesFEMWorldDef = [coordinatesFEMWorld(:,1)+ disp_u, ...
                                     coordinatesFEMWorld(:,2)+ disp_v];
 

         bubble_y = (251-mean(CircleFitPar(end-20:end-1,1)))*DICpara.um2px;
         bubble_x = mean(CircleFitPar(end-20:end-1,2))*DICpara.um2px;
         
         r = sqrt( (coordinatesFEMWorldDef(:,1)-bubble_x).^2 + (coordinatesFEMWorldDef(:,2)-bubble_y).^2  );
         theta = atan2(  (coordinatesFEMWorldDef(:,2)-bubble_y), coordinatesFEMWorldDef(:,1)-bubble_x);
         
         disp_r = cos(theta).*disp_u + sin(theta).*disp_v; % JY!!! correction on 08/15/2021
         disp_t = - sin(theta).*disp_u + cos(theta).*disp_v;
         
      
         strain_logEtt = log(1+strain_ett);
         strain_logErr = log(1+strain_err);
         
         Jacobian = (1+strain_err).*(1+strain_ett).^2; 
         
         % frame = getframe(gcf);
         % writeVideo(v,frame);
    
        %%%%%%%%% cav-dic part %%%%%%%%%%

        end
    end
  
    % ----- Save strain results ------
%     ResultStrain{ImgSeqNum-1} = struct('strainxCoord',coordinatesFEMWorld(:,1),'strainyCoord',coordinatesFEMWorld(:,2), ...
%             'dispu',UWorld(1:2:end),'dispv',UWorld(2:2:end), ...
%             'dudx',FStraintemp(1:4:end),'dvdx',FStraintemp(2:4:end),'dudy',FStraintemp(3:4:end),'dvdy',FStraintemp(4:4:end), ...
%             'strain_exx',strain_exx,'strain_exy',strain_exy,'strain_eyy',strain_eyy, ...
%             'strain_principal_max',strain_principal_max,'strain_principal_min',strain_principal_min, ...
%             'strain_maxshear',strain_maxshear,'strain_vonMises',strain_vonMises);
     
    ResultCavDIC{ImgSeqNum-1} = struct('bubble_center_x',bubble_x,'bubble_center_y',bubble_y, ...
            'r',r,'theta',theta,'disp_r',disp_r,'disp_t',disp_t, ...
            'R',Rnew(ImgSeqNum)*DICpara.um2px,'DICwinsizePhy',DICpara.winsize*DICpara.um2px,  ...
            'dudx',dudx,'dvdx',dvdx,'dudy',dudy,'dvdy',dvdy, ...
            'strain_err',strain_err,'strain_ert',strain_ert,'strain_ett',strain_ett, ...
            'strain_logErr',strain_logErr,'strain_logEtt',strain_logEtt, ...
            'Jacobian',Jacobian');


% % ------ Save figures for tracked displacement and strain fields ------
[~,imgname,imgext] = fileparts(file_name{1,ImgSeqNum}); % Find img name
  SaveFigFilesDispAndStrainQuadtree;
   % pause;
    
    
end
% ------ END of for-loop {ImgSeqNum = 2:length(ImgNormalized)} ------
fprintf('------------ Section 8 Done ------------ \n \n')


% ------ Save data again including solved strain fields ------
% results_name = ['results_',imgname,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
% save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultFEMesh','ResultFEMeshEachFrame',...
%                    'ALSub1Time','ALSub2Time','ALSolveStep','ResultStrain');

save('results_cav_dic.mat', 'ResultCavDIC');
% close(v);


%% Section 9: Compute stress
fprintf('------------ Section 9 Start ------------ \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to compute stress fields and plot stress fields
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Choose material model ------ 
DICpara.MaterialModel = funParaInput('MaterialModel');
% ------ Define parameters in material models ------
if (DICpara.MaterialModel == 1) || (DICpara.MaterialModel == 2) % Linear elasticity
    fprintf('Define Linear elasticity parameters \n')
    fprintf("Young's modulus (unit: Pa). \n "); prompt = 'Input here (e.g., 69e9): '; 
    DICpara.MaterialModelPara.YoungsModulus = input(prompt); 
    fprintf("Poisson's ratio \n"); prompt = 'Input here (e.g., 0.3): '; 
    DICpara.MaterialModelPara.PoissonsRatio = input(prompt);
    fprintf('------------------------------------- \n');
end

% ------ Start main part ------
for ImgSeqNum = 2 : length(ImgNormalized)
    
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]); close all;
    
    coordinatesFEM = ResultFEMeshEachFrame{ImgSeqNum-1}.coordinatesFEM;
    elementsFEM = ResultFEMeshEachFrame{ImgSeqNum-1}.elementsFEM;
    coordinatesFEMWorldDef = DICpara.um2px*[coordinatesFEM(:,1),size(ImgNormalized{1},2)+1-coordinatesFEM(:,2)] + ...
                             DICpara.Image2PlotResults*[ResultStrain{ImgSeqNum-1}.dispu, ResultStrain{ImgSeqNum-1}.dispv];
     
    % ------ Plot stress ------
    if DICpara.OrigDICImgTransparency == 1
        [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
                stress_principal_min_xyplane, stress_maxshear_xyplane, ...
                stress_maxshear_xyz3d, stress_vonMises]  =  Plotstress0Quadtree( ...
            DICpara,ResultStrain{ImgSeqNum-1},coordinatesFEMWorldDef,elementsFEM(:,1:4)); 
        
    else % Plot over raw DIC images
        if DICpara.Image2PlotResults == 0 % Plot over the first image; "file_name{1,1}" corresponds to the first image	
            [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
                stress_principal_min_xyplane, stress_maxshear_xyplane, ...
                stress_maxshear_xyz3d, stress_vonMises] = PlotstressQuadtree( ...
                DICpara,ResultStrain{ImgSeqNum-1},coordinatesFEMWorldDef,elementsFEM(:,1:4),file_name{1,1});
             
        else % Plot over second or next deformed images
           [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
                stress_principal_min_xyplane, stress_maxshear_xyplane, ...
                stress_maxshear_xyz3d, stress_vonMises] = PlotstressQuadtree( ...
                DICpara,ResultStrain{ImgSeqNum-1},coordinatesFEMWorldDef,elementsFEM(:,1:4),file_name{1,ImgSeqNum});
 
        end
    end
    
    
    % ------ Save figures for computed stress fields ------
    SaveFigFilesStress; 
    
    % ----- Save strain results ------
    ResultStress{ImgSeqNum-1} = struct('stressxCoord',ResultStrain{ImgSeqNum-1}.strainxCoord,'stressyCoord',ResultStrain{ImgSeqNum-1}.strainyCoord, ...
            'stress_sxx',stress_sxx,'stress_sxy',stress_sxy,'stress_syy',stress_syy, ...
            'stress_principal_max_xyplane',stress_principal_max_xyplane, 'stress_principal_min_xyplane',stress_principal_min_xyplane, ...
            'stress_maxshear_xyplane',stress_maxshear_xyplane,'stress_maxshear_xyz3d',stress_maxshear_xyz3d, ...
            'stress_vonMises',stress_vonMises);
        
end
% ------ END of for-loop {ImgSeqNum = 2:length(ImgNormalized)} ------
fprintf('------------ Section 9 Done ------------ \n \n')

% ------ Save data again including solved stress fields ------
results_name = ['results_',imgname,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultFEMesh','ResultFEMeshEachFrame', ...
                   'ALSub1Time','ALSub2Time','ALSolveStep','ResultStrain','ResultStress');


%% Section 10: Plot the generated quadtree mesh 
v = VideoWriter('video_mesh.mp4');
v.FrameRate = 5;
open(v);
figure,
for ImgSeqNum = 2 : (1+size(ResultDisp,1))
    
    clf; patch('Faces', DICmesh.elementsFEM(:,1:4), 'Vertices', DICmesh.coordinatesFEMWorld + ...
        [ResultDisp{ImgSeqNum-1}.U(1:2:end), -ResultDisp{ImgSeqNum-1}.U(2:2:end)], 'Facecolor','none','linewidth',1)
    xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    tt = title(['Frame #',num2str(ImgSeqNum)],'fontweight','normal');
    set(tt,'Interpreter','latex','fontsize',10);
    axis equal; axis tight; set(gca,'fontsize',18); set(gcf,'color','w'); box on;
    a = gca; a.TickLabelInterpreter = 'latex';
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    
end
close(v);

