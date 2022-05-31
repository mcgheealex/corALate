% ---------------------------------------------
% Augmented Lagrangian Digital Image Correlation (ALDIC_Quadtree_v2)
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
addpath("./func",'./plotFiles','./func_quadtree','./plotFiles/export_fig-d966721/', ...
    './func_quadtree/refinement','./func_quadtree/rbfinterp'); 
addpath('./Images_Quadtree_demo/Images_Sample12'); % TODO: addpath("./YOUR IMAGE FOLDER"); 
addpath('./rbfinterp');
fprintf('------------ Section 1 Done ------------ \n \n')


%% Section 2: Load DIC parameters and set up DIC parameters 
fprintf('------------ Section 2 Start ------------ \n')
% ====== Read images ====== 
close all; [file_name,Img,DICpara] = ReadImageQuadtree; % Load DIC raw images

% ====== Read image mask files ======
[mask_file_name,ImgMask] = ReadImageMasks; % Load DIC image mask files

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
fprintf('------------ Section 2 Done ------------ \n \n')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To solve each frame in an image sequence
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ImgSeqNum = 64 % 2 : length(ImgNormalized) 
 

    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);
    
    % ====== Load image mask files & Compute image gradients ======
    %%%%% Load reference image %%%%%
    fNormalizedMask = double( ImgMask{1} );
    fNormalized = ImgNormalized{1} .* fNormalizedMask; % Load the first referece image 
    Df = funImgGradient(fNormalized,fNormalized,fNormalizedMask); % Finite difference to compute image grayscale gradients;
    %%%%% Load deformed image %%%%%%%
    gNormalized = ImgNormalized{ ImgSeqNum } ; % Load current deformed image frame 
    try
        gNormalizedMask = double( ImgMask{ImgSeqNum} );
        gNormalized = gNormalized .* gNormalizedMask ;
    catch
        gNormalizedMask = ones(size(gNormalized));
    end
    DICpara.ImgRefMask = fNormalizedMask;
    
    figure, subplot(2,2,1); imshow(fNormalized'); title('Ref img f')
    subplot(2,2,2); imshow(gNormalized'); title('Def img g')
    subplot(2,2,3); imshow(fNormalizedMask'); title('f mask')
    subplot(2,2,4); imshow(gNormalizedMask'); title('g mask')
    
    
    %% Section 3: Compute an initial guess of the unknown displacement field
    fprintf('\n'); fprintf('------------ Section 3 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to find or update an initial guess of the unknown displacements.
    % The key idea is to either to use a new FFT-based cross correlation peak fitting,
    % or use the results from previous frames to compute a new initial guess for the next frame;
    % Particularly in the incremental mode DIC, the reference image can also be updated, e.g.,
    % " fNormalized = ImgNormalized{ImgSeqNum-mod(ImgSeqNum-1,ImgSeqIncUnit)}; "
    %
    % DICpara.NewFFTSearch = 0; % If you want to apply the FFT-based cross correlation to 
    % compute the initial guess for each frame, please make sure that "DICpara.NewFFTSearch = 0". 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% One practical strategy is to let first 7 frames do the FFT-based 
    %%%%% cross correlation and then let data driven method to estimate new 
    %%%%% initial guesses for other frames 
%     DICpara.NewFFTSearch = 1;
%     if ImgSeqNum < 7
%         DICpara.NewFFTSearch = 1; % Use FFT-based cross correlation to compute the initial guess
%     else
%         DICpara.NewFFTSearch = 0; % Apply data driven method to estimate initial guesses for later frames
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ImgSeqNum == 2 || DICpara.NewFFTSearch == 1 % Apply FFT-based cross correlation to compute the initial guess 
         
        % ====== FFT-based cross correlation ======
        % Switch the order of Img f and Img g
        [DICpara,x0temp_g,y0temp_g,u_g,v_g,cc]= IntegerSearch(gNormalized,fNormalized,file_name,DICpara);
        u_f = -u_g; v_f = -v_g; x0temp_f = x0temp_g+u_g; y0temp_f = y0temp_g+v_g;
        
        %%%%% Interpolate to f %%%%%
        xnodes = max([4+0.5*DICpara.winsize,DICpara.gridxyROIRange.gridx(1)])  ...
            : DICpara.winstepsize : min([size(fNormalized,1)-0.5*DICpara.winsize-3,DICpara.gridxyROIRange.gridx(2)]);
        ynodes = max([4+0.5*DICpara.winsize,DICpara.gridxyROIRange.gridy(1)])  ...
            : DICpara.winstepsize : min([size(fNormalized,2)-0.5*DICpara.winsize-3,DICpara.gridxyROIRange.gridy(2)]);
         
        [x0temp,y0temp] = ndgrid(xnodes,ynodes);
        u_f_NotNanInd = find(~isnan(u_f(:)));
         
        op1 = rbfcreate( [x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd)]',[u_f(u_f_NotNanInd)]','RBFFunction', 'thinplate'); rbfcheck(op1);
        u = rbfinterp( [x0temp(:),y0temp(:)]', op1 );
        op2 = rbfcreate( [x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd)]',[v_f(u_f_NotNanInd)]','RBFFunction', 'thinplate'); rbfcheck(op2);
        v = rbfinterp([x0temp(:),y0temp(:)]', op2 );
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % figure, plot3(x0temp(:),y0temp(:),u(:),'.');
        % figure, plot3(x0temp(:),y0temp(:),v(:),'.');
        %%%%% Do some regularization to further decrease the noise %%%%%
        u = regularizeNd([x0temp(:),y0temp(:)],u(:),{xnodes',ynodes'},1e-3);
        v = regularizeNd([x0temp(:),y0temp(:)],v(:),{xnodes',ynodes'},1e-3);

        % ====== DIC uniform FE-mesh set up ======
        [DICmesh] = MeshSetUp(x0temp,y0temp,DICpara); % clear x0temp y0temp;
        % ====== Initial Value ======
        U0 = Init(u,v,cc.max,DICmesh.x0,DICmesh.y0,0);  
         
        for tempi = 1:size(u,1)
            for tempj = 1:size(u,2)
                try
                    if ~fNormalizedMask(x0temp(tempi,tempj),y0temp(tempi,tempj)) || ...
                           ~gNormalizedMask( floor(x0temp(tempi,tempj)+u(tempi,tempj)), floor(y0temp(tempi,tempj)+v(tempi,tempj)) )
                         
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
         
        % ====== Generate a quadtree mesh considering sample's complex geometry ======
        DICmesh.elementMinSize = 2; % min element size in the refined quadtree mesh
 
        GenerateQuadtreeMesh; % Generate a quadtree mesh
        
        %%%%% Update search region %%%%%
        % DICpara.SizeOfFFTSearchRegion = [ ceil( max( [max(3+abs(U0(1:2:end))), 3] ) ), ...
        %                                   ceil( max( [max(3+abs(U0(2:2:end))), 3] ) ) ];
          
  
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
        if ImgSeqNum < 7 % Import previous U for ImgSeqNum [2,6] 
            U0 = ResultDisp{ImgSeqNum-2}.U;
             
        else % When ImgSeqNum > 6: POD predicts next disp U0 from previous results of (ImgSeqNum+[-5:1:-1])
            nTime = 5; np = length(ResultDisp{ImgSeqNum-2}.U)/2; % "nTime" value 5 is an empirical value, can be changed.
            T_data_u = zeros(nTime,np); T_data_v = zeros(nTime,np); 
            for tempi = 1:nTime
                T_data_u(tempi,:) = ResultDisp{ImgSeqNum-(2+nTime)+tempi, 1}.U(1:2:np*2)';
                T_data_v(tempi,:) = ResultDisp{ImgSeqNum-(2+nTime)+tempi, 1}.U(2:2:np*2)';
            end
            nB = 3; t_train = [ImgSeqNum-1-nTime:ImgSeqNum-2]'; t_pre = [ImgSeqNum-1]';
            [u_pred,~,~,~] = funPOR_GPR(T_data_u,t_train,t_pre,nB);
            [v_pred,~,~,~] = funPOR_GPR(T_data_v,t_train,t_pre,nB);
            tempu = u_pred(1,:); tempv = v_pred(1,:);
            U0 = [tempu(:),tempv(:)]'; U0 = U0(:);
             
        end
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
     
    % ====== ALStep 1 Subproblem1: Local Subset DIC ======
    mu=0; beta=0; tol=1e-2; ALSolveStep=1; ALSub1Time=zeros(6,1); ALSub2Time=zeros(6,1); 
    ConvItPerEle=zeros(size(DICmesh.coordinatesFEM,1),6); ALSub1BadPtNum=zeros(6,1);
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ Start Local DIC IC-GN iteration ------
    [USubpb1,FSubpb1,HtempPar,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp,markCoordHoleStrain] = ...
        LocalICGNQuadtree(U0,DICmesh.coordinatesFEM,Df,fNormalized,gNormalized,DICpara,'GaussNewton',tol);
    ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp; toc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------  Manually find some bad points from Local Subset ICGN step ------
    % Comment these lines below if you don't have local bad points
    % %%%%% Comment START %%%%%%
    %    [USubpb1,FSubpb1] = funRemoveOutliersQuadtree(DICmesh,DICpara,USubpb1,FSubpb1);
    %    disp('--- Remove bad points done ---')
    % %%%%% Comment END %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    % ------ Plot ------
    USubpb1World = USubpb1; USubpb1World(2:2:end) = -USubpb1(2:2:end); 
    FSubpb1World = FSubpb1; FSubpb1World(2:4:end) = -FSubpb1(2:4:end); FSubpb1World(3:4:end) = -FSubpb1(3:4:end); 
    close all; Plotdisp_show(USubpb1World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
    Plotstrain_show(FSubpb1World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
    save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
    fprintf('------------ Section 4 Done ------------ \n \n')

    %%%%% Plot ICGN iteration step %%%%%
    % Convtemp = [ConvItPerEletemp, ConvItPerEletemp]'; Convtemp = Convtemp(:);
    % PlotMesh_show(Convtemp,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
   
     
   
    %% Section 5: Subproblem 2 -- solve the global compatible displacement field
    fprintf('------------ Section 5 Start ------------ \n'); tic;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve the global step in ALDIC Subproblem 2
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ======= ALStep 1 Subproblem 2: Global constraint =======
    % ------ Smooth displacements for a better F ------
    DICpara.DispFilterSize=0; DICpara.DispFilterStd=0; DICpara.StrainFilterSize=0; DICpara.StrainFilterStd=0; LevelNo=1;
    DICpara.DispSmoothness = 0; DICpara.StrainSmoothness = 1e-4;
    if DICpara.DispSmoothness>1e-6, USubpb1 = funSmoothDispQuadtreeRBF(USubpb1,DICmesh,DICpara); end
    if DICpara.StrainSmoothness>1e-6, FSubpb1 = funSmoothStrainQuadtreeRBF(FSubpb1,DICmesh,DICpara); end
    
	% ====== Define penalty parameter ======
    mu = 1e-3; udual = 0*FSubpb1; vdual = 0*USubpb1; 
    betaList = [1e-3,1e-2,1e-1]*mean(DICpara.winstepsize).^2.*mu; % Tune beta in the betaList 
    Err1 = zeros(length(betaList),1); Err2 = Err1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****']);
    DICpara.GaussPtOrder = 2; alpha = 0;  % No regularization added
    % ====== Solver using finite element method ======
    % if ImgSeqNum == 2
        for tempk = 1:length(betaList)
            beta = betaList(tempk); display(['Try #',num2str(tempk),' beta = ',num2str(beta)]);
            GaussPtOrder=3; alpha=0; [USubpb2] = Subpb2Quadtree(DICmesh,DICpara.GaussPtOrder,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,mean(DICpara.winstepsize),0);
            % [FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2,DICpara.GaussPtOrder,0);
           [FSubpb2] = funGlobalNodalStrainRBF(DICmesh,DICpara,USubpb2);
             
%         [~,~,FSubpb2] = funCompDefGrad2([USubpb2(1:2:end),USubpb2(2:2:end)], ...
%                  DICmesh.coordinatesFEM, DICpara.winstepsize, 10, DICpara.ImgRefMask);
             
            
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
    %     else 
    %         try beta = DICpara.beta;
    %         catch, beta = 1e-3*mean(DICpara.winstepsize).^2.*mu;
    %         end
    %     end
      
    % Using the optimal beta to solve the ALDIC Subproblem 2 again
    if abs(beta-betaList(end))>abs(eps)
        [USubpb2] = Subpb2Quadtree(DICmesh,DICpara.GaussPtOrder,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,mean(DICpara.winstepsize),0);
        % [FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2,DICpara.GaussPtOrder,0);
        [FSubpb2] = funGlobalNodalStrainRBF(DICmesh,DICpara,USubpb2);
        ALSub2Time(ALSolveStep) = toc; toc
    end
    
    % ------- Smooth strain field --------
    if DICpara.DispSmoothness>1e-6, USubpb2 = funSmoothDispQuadtreeRBF(USubpb2,DICmesh,DICpara); end
    % ------- Don't smooth strain fields near the boundary --------
    for tempk=0:3, FSubpb2(4*DICmesh.markCoordHoleEdge-tempk) = FSubpb1(4*DICmesh.markCoordHoleEdge-tempk); end
    if DICpara.StrainSmoothness>1e-6, FSubpb2 = funSmoothStrainQuadtreeRBF(0.1*FSubpb2+0.9*FSubpb1,DICmesh,DICpara); end
    for tempk=0:1, USubpb2(2*markCoordHoleStrain-tempk) = USubpb1(2*markCoordHoleStrain-tempk); end
    for tempk=0:3, FSubpb2(4*markCoordHoleStrain-tempk) = FSubpb1(4*markCoordHoleStrain-tempk); end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ------- Save data ------
    save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2');

    % ------ Plot ------
    USubpb2World = USubpb2; USubpb2World(2:2:end) = -USubpb2(2:2:end);
    FSubpb2World = FSubpb2; FSubpb2World(2:4:end) = -FSubpb2World(2:4:end); FSubpb2World(3:4:end) = -FSubpb2World(3:4:end); 
    close all; Plotdisp_show(USubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
    Plotstrain_show(FSubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');

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
        
        %%%%% If wen want to adjust the DIC subset size %%%%%
        % Ftemp1 = FSubpb2(1:2:end); Ftemp2 = FSubpb2(2:2:end);
        % [DFtemp1,~,~] = funGlobalNodalStrainQuadtree(DICmesh,Ftemp1,DICpara.GaussPtOrder,0);
        % [DFtemp2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,Ftemp2,DICpara.GaussPtOrder,0);
        % 
        % winsize_x_ub1 = abs(2*FSubpb2(1:4:end)./DFtemp1(1:4:end));
        % winsize_x_ub2 = abs(2*FSubpb2(3:4:end)./DFtemp1(3:4:end));
        % winsize_y_ub1 = abs(2*FSubpb2(1:4:end)./DFtemp1(2:4:end));
        % winsize_y_ub2 = abs(2*FSubpb2(3:4:end)./DFtemp1(4:4:end));
        % 
        % winsize_x_ub3 = abs(2*FSubpb2(2:4:end)./DFtemp2(1:4:end));
        % winsize_x_ub4 = abs(2*FSubpb2(4:4:end)./DFtemp2(3:4:end));
        % winsize_y_ub3 = abs(2*FSubpb2(2:4:end)./DFtemp2(2:4:end));
        % winsize_y_ub4 = abs(2*FSubpb2(4:4:end)./DFtemp2(4:4:end));
        % 
        % winsize_x_ub = round(min([winsize_x_ub1,winsize_x_ub2,winsize_x_ub3,winsize_x_ub4,DICpara.winsize*ones(length(winsize_x_ub1),1)],[],2));
        % winsize_x_List = max([winsize_x_ub, 10*ones(length(winsize_x_ub1),1)],[],2);
        % winsize_y_ub = round(min([winsize_y_ub1,winsize_y_ub2,winsize_y_ub3,winsize_y_ub4,DICpara.winsize*ones(length(winsize_y_ub1),1)],[],2));
        % winsize_y_List = max([winsize_y_ub, 10*ones(length(winsize_y_ub1),1)],[],2);
        % winsize_List = 2*ceil([winsize_x_List,winsize_y_List]/2);
        
        %%%%% Otherwise, we use the fixed DIC subset size %%%%%
        winsize_List = DICpara.winsize*ones(size(DICmesh.coordinatesFEM,1),2);
        DICpara.winsize_List = winsize_List;
        

        %%%%%%%%%%%%%%%%%%%%%%% Subproblem 1 %%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****']);
        tic; [USubpb1,~,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp] = Subpb1Quadtree(...
                                            USubpb2,FSubpb2,udual,vdual,DICmesh.coordinatesFEM,...
                                            Df,fNormalized,gNormalized,mu,beta,HPar,ALSolveStep,DICpara,'GaussNewton',tol);
        FSubpb1 = FSubpb2; toc 
        for tempk=0:1, USubpb1(2*markCoordHoleStrain-tempk) = USubpb2(2*markCoordHoleStrain-tempk); end
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
		% [FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2,DICpara.GaussPtOrder,0);
        [FSubpb2] = funGlobalNodalStrainRBF(DICmesh,DICpara,USubpb2);
        
%         [~,~,FSubpb2] = funCompDefGrad2([USubpb2(1:2:end),USubpb2(2:2:end)], ...
%                  DICmesh.coordinatesFEM, DICpara.winstepsize, 10, DICpara.ImgRefMask);
             
        ALSub2Time(ALSolveStep) = toc; toc
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ------- Smooth strain field --------
        if DICpara.DispSmoothness>1e-6, USubpb2 = funSmoothDispQuadtreeRBF(USubpb2,DICmesh,DICpara); end
        % ------- Don't change strain fields near the boundary --------
        for tempk=0:3, FSubpb2(4*DICmesh.markCoordHoleEdge-tempk) = FSubpb1(4*DICmesh.markCoordHoleEdge-tempk); end
        if DICpara.StrainSmoothness>1e-6, FSubpb2 = funSmoothStrainQuadtreeRBF(0.1*FSubpb2+0.9*FSubpb1,DICmesh,DICpara); end
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
USubpb2World = USubpb2; USubpb2World(2:2:end) = -USubpb2(2:2:end); 
FSubpb2World = FSubpb2; FSubpb2World(2:4:end) = -FSubpb2World(2:4:end); FSubpb2World(3:4:end) = -FSubpb2World(3:4:end); 
close all; Plotdisp_show(USubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'NoEdgeColor');
Plotstrain_show(FSubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'NoEdgeColor');

% winsize_xy = [winsize_x_List, winsize_y_List]'; winsize_xy = winsize_xy(:);
% Plotdisp_show(winsize_xy,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
        
% ------ Save results ------
% Find img name and save all the results 
[~,imgname,imgext] = fileparts(file_name{1,end});
results_name = ['results_',imgname,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultFEMesh','ResultFEMeshEachFrame','ALSub1Time','ALSub2Time','ALSolveStep');


figure(1); colormap(black_rainbow_plus); caxis([-0.5,0.5]);
savefig(['fig_',num2str(LevelNo),'_dispx.fig']);
figure(2); colormap(black_rainbow_plus); caxis([-0.1,0.1]);
savefig(['fig_',num2str(LevelNo),'_dispy.fig']);
figure(3); colormap(black_rainbow_plus); caxis([-0.004,0.004]);
savefig(['fig_',num2str(LevelNo),'_strain_e11.fig']);
figure(4); colormap(black_rainbow_plus); caxis([-0.002 ,0.002 ]);
savefig(['fig_',num2str(LevelNo),'_strain_e12.fig']);
figure(5); colormap(black_rainbow_plus); caxis([-0.001 ,0.001 ]);
savefig(['fig_',num2str(LevelNo),'_strain_e22.fig']);











%% ============ Level 1 --> Level 2 ==============
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ULevel1 = USubpb2; FLevel1 = FSubpb2; udualLevel1 = udual;  vdualLevel1 = vdual;
elementsFEM = elementsFEMQuadtree;
coordinatesFEM = coordinatesFEMQuadtree;
irregularEdge = irregular;

betavar = beta; muvar = mu; alphavar = 0;

CrackOrNot = 0; CrackPath1 = [0,0]; CrackPath2 = [0,0]; CrackTipOrNot = 0; CrackTip = [0,0]; 
EnrHAndTipEleIndex = []; EnrTipEleIndex = []; dirichlet = []; neumann= []; imgPyramidUnit = 1;

winstepsize = DICpara.winstepsize;

eleGeneration = ones(size(elementsFEM,1),1);
for eleInd = 1:size(elementsFEM,1)
    eleLength = abs(coordinatesFEM(elementsFEM(eleInd,1),1) - coordinatesFEM(elementsFEM(eleInd,3),1));
    eleGeneration(eleInd) = 1 + log2(winstepsize/eleLength);
end

%% ============== Estimate: Level1 ==============
% ------ Compute rho to find where to refine ------
rhoLevel1Vector = zeros(size(elementsFEM,1),1); rhoLevel1Vectortemp1 = rhoLevel1Vector; rhoLevel1Vectortemp2 = rhoLevel1Vector;
ULevel1temp = [ULevel1;0;0]; hbar = parfor_progressbar(size(elementsFEM,1),'Please wait for a posterior error estimate!'); tic

for j = 1:size(elementsFEM,1)
    % eleNeighborIndexAndEdge = findEleNeighbors(elementsFEMIter, coordinatesFEMIter, 23 ,eleGenerationIter, winstepsize);
    point1x = coordinatesFEM(elementsFEM(j,1),1); point3x = coordinatesFEM(elementsFEM(j,3),1);lengthOfElement = (point3x-point1x)/winstepsize;
    % ---- Compute error estimator ----
    rhoLevel1Vectortemp1(j) = (lengthOfElement)^2 * aPostErrEstInt_ALSq(coordinatesFEM,elementsFEM,j,betavar,muvar,ULevel1,FLevel1,udualLevel1,vdualLevel1,imgPyramidUnit,CrackOrNot);
    rhoLevel1Vectortemp2(j) = (lengthOfElement) * (betavar+alphavar)^2 * aPostErrEstJump_ALSq(coordinatesFEM,elementsFEM,j,eleGeneration,winstepsize,ULevel1, ...
                               EnrHAndTipEleIndex,EnrTipEleIndex,CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot,DICmesh.markCoordHoleEdge);
    rhoLevel1Vector(j) = rhoLevel1Vectortemp1(j) + rhoLevel1Vectortemp2(j);
    hbar.iterate(1);
end
close(hbar); EstimateTime = toc;
Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(1,size(elementsFEM,1)); Sqctemp1 = Sqc; Sqctemp2 = Sqc;
for tempj = 1:size(elementsFEM,1)
    Sqx(1:4,tempj) = coordinatesFEM(elementsFEM(tempj,1:4),1);Sqy(1:4,tempj) = coordinatesFEM(elementsFEM(tempj,1:4),2); 
    Sqc(tempj) = rhoLevel1Vector(tempj);Sqctemp1(tempj) = rhoLevel1Vectortemp1(tempj);Sqctemp2(tempj) = rhoLevel1Vectortemp2(tempj);
end
close all;
figure; patch(Sqx,Sqy,Sqc,'edgecolor','none' ); c=colorbar; title('Ele total err estimator','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal; % set(gca,'XTick',[]);
figure; patch(Sqx,Sqy,Sqctemp1,'edgecolor','none');c=colorbar; title('Ele interior err estimator','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal; % set(gca,'XTick',[]);
figure; patch(Sqx,Sqy,Sqctemp2,'edgecolor','none');c=colorbar; title('Ele jump err estimator','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal; % set(gca,'XTick',[]);
        
figure(1); set(gcf,'color','w'); savefig(['fig_L',num2str(LevelNo),'_err_all.fig']);
figure(2); set(gcf,'color','w'); savefig(['fig_L',num2str(LevelNo),'_err_i.fig']);
figure(3); set(gcf,'color','w'); savefig(['fig_L',num2str(LevelNo),'_err_j.fig']);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DorflerSatisfied = 1; 
while DorflerSatisfied == 1
    fprintf('--- Dorfler strategy to mark elements to refine? --- \n')
    fprintf('Input a real number between 0 and 1 here: ');
    prompt = 'Input here: '; thetaDorfler = input(prompt);
     
    %% ================ Mark: Level1 ==================
    % ======= Dorfler's Strategy to Mark elements to refine =======
    tic; [rhoLevel1VectorSorted,rhoLevel1VectorSortIndex] = sort(rhoLevel1Vector,'descend');
    rhoLevel1SumM = 0; rhoLevel1SumAll = sum(rhoLevel1Vector(:)); rhoLevel1SumMIndex = 0;
    while rhoLevel1SumM < 0.999999*thetaDorfler*rhoLevel1SumAll
        rhoLevel1SumMIndex = rhoLevel1SumMIndex + 1;
        rhoLevel1SumM = rhoLevel1SumM + rhoLevel1VectorSorted(rhoLevel1SumMIndex);
    end
    refinementPatch = rhoLevel1VectorSortIndex(1:rhoLevel1SumMIndex);
    if CrackOrNot > 0
        refinementPatch = union(refinementPatch,[EleCrackTop(end)+1:1:EleCrackBottom(1)-1]');
    end
    Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(1,size(elementsFEM,1));
    for j = 1:size(elementsFEM,1); Sqx(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),1); Sqy(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),2); end
    for j = 1:size(refinementPatch,1); Sqc(refinementPatch(j)) = 1; end
    MarkTime = toc;
    figure; patch(Sqx,Sqy,Sqc,'EdgeColor','none'); caxis([0,1]); colorbar;
    title('Marked ele (yellow)','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal; 
    set(gcf,'color','w');

    savefig(['fig_L',num2str(LevelNo+1),'_mark.fig']); % JY!!!
 
    
    %% ============== Refine: Level1 ==============
    if thetaDorfler>0 && size(coordinatesFEM,1)>1e3  
    
        tic; [coordinatesFEMLevel2,elementsFEMLevel2,irregularEdgeLevel2,dirichletLevel2,neumannLevel2] = ...
            QrefineRwBCNormal(coordinatesFEM,elementsFEM(:,1:4),irregularEdge,dirichlet,[],refinementPatch);
        elementsFEMLevel2 = [elementsFEMLevel2,zeros(size(elementsFEMLevel2,1),4)];
        
        % Assign element generation #
        eleGenerationLevel2 = ones(size(elementsFEMLevel2,1),1);
        for tempj=1:size(elementsFEMLevel2,1)
            eleGenerationLevel2(tempj) = 1-log2(norm(coordinatesFEMLevel2(elementsFEMLevel2(tempj,1),:) - coordinatesFEMLevel2(elementsFEMLevel2(tempj,2),:))/winstepsize);
        end
         
        % Reorder element nodes
        for tempj=1:length(eleGenerationLevel2)
            [~,row2] = max(  sum(coordinatesFEMLevel2(elementsFEMLevel2(tempj,1:4)',:)')  );
            if row2==1
                elementsFEMLevel2(tempj,1:4) = elementsFEMLevel2(tempj,[3,4,1,2]);
            elseif row2==2
                elementsFEMLevel2(tempj,1:4) = elementsFEMLevel2(tempj,[4,1,2,3]);
            elseif row2==4
                elementsFEMLevel2(tempj,1:4) = elementsFEMLevel2(tempj,[2,3,4,1]);
            end
        end
        
        % For hanging nodes in FE-mesh
        for tempj=1:size(irregularEdgeLevel2,1)
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,1:2), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,8)=irregularEdgeLevel2(tempj,3); end
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,2:3), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,5)=irregularEdgeLevel2(tempj,3); end
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,3:4), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,6)=irregularEdgeLevel2(tempj,3); end
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,[4,1]), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,7)=irregularEdgeLevel2(tempj,3); end
            
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,[2,1]), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,8)=irregularEdgeLevel2(tempj,3); end
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,[3,2]), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,5)=irregularEdgeLevel2(tempj,3); end
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,[4,3]), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,6)=irregularEdgeLevel2(tempj,3); end
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,[1,4]), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,7)=irregularEdgeLevel2(tempj,3); end
        end
        RefineTime = toc;
        
        % Initialize variables F0_Level2 and U0_Level2
        F_dispu = scatteredInterpolant( coordinatesFEM(:,1),coordinatesFEM(:,2),ULevel1(1:2:end) );
        F_dispv = scatteredInterpolant( coordinatesFEM(:,1),coordinatesFEM(:,2),ULevel1(2:2:end) );
        F_F11 = scatteredInterpolant( coordinatesFEM(:,1),coordinatesFEM(:,2),FLevel1(1:4:end) );
        F_F21 = scatteredInterpolant( coordinatesFEM(:,1),coordinatesFEM(:,2),FLevel1(2:4:end) );
        F_F12 = scatteredInterpolant( coordinatesFEM(:,1),coordinatesFEM(:,2),FLevel1(3:4:end) );
        F_F22 = scatteredInterpolant( coordinatesFEM(:,1),coordinatesFEM(:,2),FLevel1(4:4:end) );
        
        U0Level2 = 0*coordinatesFEMLevel2(:); F0Level2=0*[U0Level2(:);U0Level2(:)];
        temp = F_dispu(coordinatesFEMLevel2(:,1),coordinatesFEMLevel2(:,2)); U0Level2(1:2:end)=temp(:);
        temp = F_dispv(coordinatesFEMLevel2(:,1),coordinatesFEMLevel2(:,2)); U0Level2(2:2:end)=temp(:);
        temp = F_F11(coordinatesFEMLevel2(:,1),coordinatesFEMLevel2(:,2)); F0Level2(1:4:end)=temp(:);
        temp = F_F21(coordinatesFEMLevel2(:,1),coordinatesFEMLevel2(:,2)); F0Level2(2:4:end)=temp(:);
        temp = F_F12(coordinatesFEMLevel2(:,1),coordinatesFEMLevel2(:,2)); F0Level2(3:4:end)=temp(:);
        temp = F_F22(coordinatesFEMLevel2(:,1),coordinatesFEMLevel2(:,2)); F0Level2(4:4:end)=temp(:);
     
    else  
             
        tic; refinedEleIDListLevel1 = 0;
        coordinatesFEMLevel2 = coordinatesFEM; elementsFEMLevel2 = elementsFEM; U0Level2 = ULevel1; F0Level2 = FLevel1;
        dirichletLevel2 = unique(dirichlet(:)); neumannLevel2 = neumann; eleGenerationLevel2 = eleGeneration;
        h=waitbar(0,'Wait for mesh refinement.');
        %tempNeuMat = [];
        for j =  1:size(refinementPatch,1)
            checkAlreadyRefinedOrNot = ismember(refinementPatch(j),refinedEleIDListLevel1);
            if checkAlreadyRefinedOrNot == 0
                FinishRefineRecursive = 1;
                [coordinatesFEMLevel2,elementsFEMLevel2,U0Level2,F0Level2,refinedEleIDListNew,dirichletNew,neumannNew,eleGenerationLevel2] = ...
                    refineRecursiveSq(coordinatesFEMLevel2,elementsFEMLevel2,eleGenerationLevel2,U0Level2,F0Level2,refinementPatch(j),winstepsize,FinishRefineRecursive,...
                    CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot);
                
                refinedEleIDListLevel1 = unique([refinedEleIDListLevel1;refinedEleIDListNew]);
                dirichletLevel2 = unique([dirichletLevel2;dirichletNew']);
                % this step is wrong. newmann size is L*4:
                if size(neumannNew,1)>1
                    %tempNeuMat = [tempNeuMat; refinementPatch(j)];
                    [rowtemp1,~] = find( (neumannLevel2(:,1) == neumannNew(2,1))  );
                    [rowtemp2,~] = find( (neumannLevel2(:,2) == neumannNew(2,2))  );
                    [rowtemp3,~] = find( (neumannLevel2(:,2) == neumannNew(2,1))  );
                    [rowtemp4,~] = find( (neumannLevel2(:,1) == neumannNew(2,2))  );
                    rowtemp = union(intersect(rowtemp1,rowtemp2), intersect(rowtemp3,rowtemp4))  ;
                    if neumannLevel2(rowtemp,1) ~= neumannNew(2,3)
                        LengthOfNeumanntemp = size(neumannLevel2,1);
                        neumannLevel2(LengthOfNeumanntemp+1,1:4) = [neumannLevel2(rowtemp,1),neumannNew(2,3),neumannLevel2(rowtemp,3:4)];
                        neumannLevel2(rowtemp,1) = [neumannNew(2,3)];
                    end
                end
                % neumannLevel2 = unique([neumannLevel2;neumannNew']);
            end
            waitbar(j/size(refinementPatch,1));
        end
        RefineTime = toc; close(h); dirichletLevel2 = dirichletLevel2(2:end);
        if CrackOrNot > 0
            [CoordCrackTopNextIter,~,CoordCrackBottomNextIter,~] = funCrackTopOrBottom(coordinatesFEMLevel2,elementsFEMLevel2,CrackPath1,CrackPath2,CrackTip);
            F0Level2 = funSmoothStrainCrack(F0Level2,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,CoordCrackTopNextIter,CoordCrackBottomNextIter,LevelNo+1,0,0);
        end
        
        % Update irregularEdgeLevel2
        irregularEdgeLevel2 = [0,0,0];
        [row,~]=find(elementsFEMLevel2(:,5)>0);
        for tempi=1:length(row)
            irregularEdgeLevel2 = [irregularEdgeLevel2; elementsFEMLevel2(row(tempi),2:3),  elementsFEMLevel2(row(tempi),5) ];
        end
        [row,~]=find(elementsFEMLevel2(:,6)>0);
        for tempi=1:length(row)
            irregularEdgeLevel2 = [irregularEdgeLevel2; elementsFEMLevel2(row(tempi),3:4),  elementsFEMLevel2(row(tempi),6) ];
        end
        [row,~]=find(elementsFEMLevel2(:,7)>0);
        for tempi=1:length(row)
            irregularEdgeLevel2 = [irregularEdgeLevel2; elementsFEMLevel2(row(tempi),[4,1]),  elementsFEMLevel2(row(tempi),7) ];
        end
        [row,~]=find(elementsFEMLevel2(:,8)>0);
        for tempi=1:length(row)
            irregularEdgeLevel2 = [irregularEdgeLevel2; elementsFEMLevel2(row(tempi),1:2),  elementsFEMLevel2(row(tempi),8) ];
        end
        
        dirichletLevel2 = [];
         
    end % END of mesh refinement
    
    

    % ====== Show refined mesh ======
    Sqx = zeros(4,size(elementsFEMLevel2,1)); Sqy = zeros(4,size(elementsFEMLevel2,1)); Sqc = zeros(1,size(elementsFEMLevel2,1));
    for j = 1:size(elementsFEMLevel2,1)
        Sqx(1:4,j) = coordinatesFEMLevel2(elementsFEMLevel2(j,1:4),1);
        Sqy(1:4,j) = coordinatesFEMLevel2(elementsFEMLevel2(j,1:4),2); Sqc(j) =  1;
    end
    figure; patch(Sqx,Sqy,Sqc,'FaceColor','none');  axis equal; axis tight; title(['Level ',num2str(LevelNo+1),' Mesh'],'fontweight','normal'); 
    set(gca,'fontsize',16); set(gcf,'color','w');

    savefig(['fig_L',num2str(LevelNo+1),'_refine.fig']); % JY!!!
     
    fprintf('--- Are you satisfied with current Dorfler refinement number (0-yes; 1-no)? ---  \n')
    prompt = 'Input here: ';
    DorflerSatisfied = input(prompt);
    if DorflerSatisfied == 0
        fprintf('--- Mesh refinement done! ---  \n')
    end

    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LevelNo = LevelNo + 1; 
coordinatesFEMIter = coordinatesFEMLevel2;  
coordinatesFEMWorldIter = coordinatesFEMIter; coordinatesFEMWorldIter(:,2) = DICpara.ImgSize(2)+1-coordinatesFEMWorldIter(:,2);
elementsFEMIter = elementsFEMLevel2;
irregularEdgeIter = irregularEdgeLevel2;  

DICmeshIter.coordinatesFEM = coordinatesFEMLevel2;
DICmeshIter.elementsFEM = elementsFEMLevel2;
DICmeshIter.dirichlet = []; DICmeshIter.neumann = [];
DICmeshIter.coordinatesFEMWorld = coordinatesFEMWorldIter;
DICmeshIter.elementMinSize = DICmesh.elementMinSize;
DICmeshIter.irregular = irregularEdgeLevel2;

U0Iter = U0Level2; ConvItPerEleIter = zeros(size(coordinatesFEMIter,1),6);



    %% Section 4: ALDIC Subproblem 1 -or- Local ICGN Subset DIC
    fprintf('------------ Section 4 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve the first local step in ALDIC: Subproblem 1
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ====== ALStep 1 Subproblem1: Local Subset DIC ======
    mu=0; beta=0; tol=1e-2; ALSolveStep=1; ALSub1Time=zeros(6,1); ALSub2Time=zeros(6,1);
    ConvItPerEleIter=zeros(size(coordinatesFEMIter,1),6); ALSub1BadPtNum=zeros(6,1);
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ Start Local DIC IC-GN iteration ------
    [USubpb1Iter,FSubpb1Iter,HtempPar,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp,markCoordHoleStrain] = ...
        LocalICGNQuadtree(U0Iter,coordinatesFEMIter,Df,fNormalized,gNormalized,DICpara,'GaussNewton',tol);
    ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEleIter(:,ALSolveStep) = ConvItPerEletemp; ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp; toc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------  Manually find some bad points from Local Subset ICGN step ------
    % Comment these lines below if you don't have local bad points
    % %%%%% Comment START %%%%%%
    %    [USubpb1,FSubpb1] = funRemoveOutliersQuadtree(DICmesh,DICpara,USubpb1,FSubpb1);
    %    disp('--- Remove bad points done ---')
    % %%%%% Comment END %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ------ Plot ------
    USubpb1WorldIter = USubpb1Iter; USubpb1WorldIter(2:2:end) = -USubpb1Iter(2:2:end); 

    close all; Plotdisp_show(USubpb1WorldIter, coordinatesFEMIter, elementsFEMIter(:,1:4),DICpara,'EdgeColor');
    Plotstrain_show(FSubpb1Iter, coordinatesFEMIter, elementsFEMIter(:,1:4),DICpara,'EdgeColor');
    save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
    fprintf('------------ Section 4 Done ------------ \n \n')
 
   
    %% Section 5: Subproblem 2 -- solve the global compatible displacement field
    fprintf('------------ Section 5 Start ------------ \n'); tic;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve the global step in ALDIC Subproblem 2
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ======= ALStep 1 Subproblem 2: Global constraint =======
    % ------ Smooth displacements for a better F ------
    DICpara.DispFilterSize=0; DICpara.DispFilterStd=0; DICpara.StrainFilterSize=0; DICpara.StrainFilterStd=0;  
    DICpara.DispSmoothness = 0; DICpara.StrainSmoothness = 1e-4;
    if DICpara.DispSmoothness>1e-6, USubpb1Iter = funSmoothDispQuadtreeRBF(USubpb1Iter,DICmeshIter,DICpara); end %%%%% should be Iter
    if DICpara.StrainSmoothness>1e-6, FSubpb1Iter = funSmoothStrainQuadtreeRBF(FSubpb1Iter,DICmeshIter,DICpara); end  %%%%% should be Iter
    
	% ====== Define penalty parameter ======
    muIter = 1e-3; udualIter = 0*FSubpb1Iter; vdualIter = 0*USubpb1Iter; 
    betaListIter = [1e-3,1e-2,1e-1]*mean(DICpara.winstepsize).^2.*muIter; % Tune beta in the betaList 
    Err1 = zeros(length(betaList),1); Err2 = Err1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****']);
    DICpara.GaussPtOrder = 2; alpha = 0;  % No regularization added
    % ====== Solver using finite element method ======
    % if ImgSeqNum == 2
        for tempk = 1:length(betaListIter)
            betaIter = betaListIter(tempk); display(['Try #',num2str(tempk),' beta = ',num2str(betaIter)]);
            alphaIter = 100*betaIter;
            GaussPtOrder=3; alpha=0; [USubpb2Iter] = Subpb2Quadtree(DICmeshIter,DICpara.GaussPtOrder, ...
                betaIter,muIter,USubpb1Iter,FSubpb1Iter,udualIter,vdualIter,alphaIter,mean(DICpara.winstepsize),0);
            
            %  [FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2,DICpara.GaussPtOrder,0);
        [FSubpb2Iter] = funGlobalNodalStrainRBF(DICmeshIter,DICpara,USubpb2Iter);
        %[~,~,FSubpb2Iter] = funCompDefGrad2([USubpb2Iter(1:2:end),USubpb2Iter(2:2:end)], ...
        %         DICmeshIter.coordinatesFEM, DICpara.winstepsize, 20, DICpara.ImgRefMask);
            
            Err1(tempk) = norm(USubpb1Iter-USubpb2Iter,2);
            Err2(tempk) = norm(FSubpb1Iter-FSubpb2Iter,2);
        end
        
        Err1Norm = (Err1-mean(Err1))/std(Err1); % figure, plot(Err1Norm);
        Err2Norm = (Err2-mean(Err2))/std(Err2); % figure, plot(Err2Norm);
        ErrSum = Err1Norm+Err2Norm; % figure, plot(ErrSum); title('Tune the best \beta in the subproblem 2'); 
        [~,indexOfbeta] = min(ErrSum); 
     
        try % Tune the best beta by a quadratic polynomial 0fitting
            [fitobj] = fit(log10(betaList(indexOfbeta-1:1:indexOfbeta+1))',ErrSum(indexOfbeta-1:1:indexOfbeta+1),'poly2');
            p = coeffvalues(fitobj); betaIter = 10^(-p(2)/2/p(1));
        catch, betaIter = betaListIter(indexOfbeta);
        end
        display(['Best beta = ',num2str(betaIter)]);
    % else 
    %    try betaIter = DICpara.beta;
    %    catch, betaIter = 1e-3*mean(DICpara.winstepsize).^2.*muIter;
    %    end
    %end
      
    % Using the optimal beta to solve the ALDIC Subproblem 2 again
    if abs(betaIter-betaListIter(end))>abs(eps)
        alphaIter = 100*betaIter;
        [USubpb2Iter] = Subpb2Quadtree(DICmeshIter,DICpara.GaussPtOrder,betaIter,muIter,USubpb1Iter,FSubpb1Iter,udualIter,vdualIter,alphaIter,mean(DICpara.winstepsize),0);
        % [FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2,DICpara.GaussPtOrder,0);
        [FSubpb2Iter] = funGlobalNodalStrainRBF(DICmeshIter,DICpara,USubpb2Iter);
        
       % [~,~,FSubpb2Iter] = funCompDefGrad2([USubpb2Iter(1:2:end),USubpb2Iter(2:2:end)], ...
       %          DICmeshIter.coordinatesFEM, DICpara.winstepsize, 20, DICpara.ImgRefMask);
             
        ALSub2Time(ALSolveStep) = toc; toc
    end
    
    % ------- Smooth strain field --------
    if DICpara.DispSmoothness>1e-6, USubpb2Iter = funSmoothDispQuadtreeRBF(USubpb2Iter,DICmeshIter,DICpara); end
    % ------- Don't smooth strain fields near the boundary --------
    for tempk=0:3, FSubpb2Iter(4*DICmesh.markCoordHoleEdge-tempk) = FSubpb1Iter(4*DICmesh.markCoordHoleEdge-tempk); end
    if DICpara.StrainSmoothness>1e-6, FSubpb2Iter = funSmoothStrainQuadtreeRBF(0.1*FSubpb2Iter+0.9*FSubpb1Iter,DICmeshIter,DICpara); end
    for tempk=0:1, USubpb2Iter(2*markCoordHoleStrain-tempk) = USubpb1Iter(2*markCoordHoleStrain-tempk); end
    for tempk=0:3, FSubpb2Iter(4*markCoordHoleStrain-tempk) = FSubpb1Iter(4*markCoordHoleStrain-tempk); end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ------- Save data ------
    save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2');

    % ------ Plot ------
    USubpb2WorldIter = USubpb2Iter; USubpb2WorldIter(2:2:end) = -USubpb2Iter(2:2:end);
    FSubpb2WorldIter = FSubpb2Iter; FSubpb2WorldIter(2:4:end) = -FSubpb2WorldIter(2:4:end); FSubpb2WorldIter(3:4:end) = -FSubpb2WorldIter(3:4:end); 
    close all; Plotdisp_show(USubpb2WorldIter,DICmeshIter.coordinatesFEMWorld,DICmeshIter.elementsFEM(:,1:4),DICpara,'EdgeColor');
    Plotstrain_show(FSubpb2WorldIter,DICmeshIter.coordinatesFEMWorld,DICmeshIter.elementsFEM(:,1:4),DICpara,'EdgeColor');

    % ======= Update dual variables =======
    udualIter = FSubpb2Iter - FSubpb1Iter; vdualIter = USubpb2Iter - USubpb1Iter;
	save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');
    fprintf('------------ Section 5 Done ------------ \n \n')
    
    
    
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
         
        %%%%% Otherwise, we use the fixed DIC subset size %%%%%
        winsize_List = DICpara.winsize*ones(size(DICmeshIter.coordinatesFEM,1),2);
        DICpara.winsize_List = winsize_List;
        

        %%%%%%%%%%%%%%%%%%%%%%% Subproblem 1 %%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****']);
        tic; [USubpb1Iter,~,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp] = Subpb1Quadtree(...
                                            USubpb2Iter,FSubpb2Iter,udualIter,vdualIter,DICmeshIter.coordinatesFEM,...
                                            Df,fNormalized,gNormalized,muIter,betaIter,HPar,ALSolveStep,DICpara,'GaussNewton',tol);
        FSubpb1Iter = FSubpb2Iter; toc 
        for tempk=0:1, USubpb1Iter(2*markCoordHoleStrain-tempk) = USubpb2Iter(2*markCoordHoleStrain-tempk); end
        ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEleIter(:,ALSolveStep) = ConvItPerEletemp; ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp;
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
        tic; [USubpb2Iter] = Subpb2Quadtree(DICmeshIter,DICpara.GaussPtOrder,betaIter,muIter, ...
            USubpb1Iter,FSubpb1Iter,udualIter,vdualIter,alphaIter,mean(DICpara.winstepsize),0);
		% [FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2,DICpara.GaussPtOrder,0);
        [FSubpb2Iter] = funGlobalNodalStrainRBF(DICmeshIter,DICpara,USubpb2Iter);
       % [~,~,FSubpb2Iter] = funCompDefGrad2([USubpb2Iter(1:2:end),USubpb2Iter(2:2:end)], ...
       %          DICmeshIter.coordinatesFEM, DICpara.winstepsize, 10, DICpara.ImgRefMask);
        ALSub2Time(ALSolveStep) = toc; toc
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ------- Smooth strain field --------
        if DICpara.DispSmoothness>1e-6, USubpb2Iter = funSmoothDispQuadtreeRBF(USubpb2Iter,DICmeshIter,DICpara); end
        % ------- Don't change strain fields near the boundary --------
        for tempk=0:3, FSubpb2Iter(4*DICmesh.markCoordHoleEdge-tempk) = FSubpb1Iter(4*DICmesh.markCoordHoleEdge-tempk); end
        if DICpara.StrainSmoothness>1e-6, FSubpb2Iter = funSmoothStrainQuadtreeRBF(0.1*FSubpb2Iter+0.9*FSubpb1Iter,DICmeshIter,DICpara); end
        for tempk=0:1, USubpb2Iter(2*markCoordHoleStrain-tempk) = USubpb1Iter(2*markCoordHoleStrain-tempk); end
        for tempk=0:3, FSubpb2Iter(4*markCoordHoleStrain-tempk) = FSubpb1Iter(4*markCoordHoleStrain-tempk); end
         
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
        udualIter = FSubpb2Iter - FSubpb1Iter; vdualIter = USubpb2Iter - USubpb1Iter; 
		 
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
    % ResultDisp{ImgSeqNum-1}.U = full(USubpb2);
    % ResultDisp{ImgSeqNum-1}.ALSub1BadPtNum = ALSub1BadPtNum;
    % ResultDefGrad{ImgSeqNum-1}.F = full(FSubpb2);  
    % close all; Plotdisp_show(USubpb2Iter,DICmeshIter.coordinatesFEM,DICmeshIter.elementsFEM(:,1:4),DICpara,'EdgeColor');
    % Plotstrain_show(FSubpb2Iter,DICmeshIter.coordinatesFEM,DICmeshIter.elementsFEM(:,1:4),DICpara,'EdgeColor');

    % ------ Plot ------
    USubpb2WorldIter = USubpb2Iter; USubpb2WorldIter(2:2:end) = -USubpb2Iter(2:2:end);
    FSubpb2WorldIter = FSubpb2Iter; FSubpb2WorldIter(2:4:end) = -FSubpb2WorldIter(2:4:end); FSubpb2WorldIter(3:4:end) = -FSubpb2WorldIter(3:4:end);
    close all; Plotdisp_show(USubpb2WorldIter,DICmeshIter.coordinatesFEMWorld,DICmeshIter.elementsFEM(:,1:4),DICpara,'NoEdgeColor');
    Plotstrain_show(FSubpb2WorldIter,DICmeshIter.coordinatesFEMWorld,DICmeshIter.elementsFEM(:,1:4),DICpara,'NoEdgeColor');
    
    

figure(1); colormap(black_rainbow_plus); caxis([-0.5,0.5]);
savefig(['fig_',num2str(LevelNo),'_dispx.fig']);
figure(2); colormap(black_rainbow_plus); caxis([-0.1,0.1]);
savefig(['fig_',num2str(LevelNo),'_dispy.fig']);
figure(3); colormap(black_rainbow_plus); caxis([-0.004,0.004]);
savefig(['fig_',num2str(LevelNo),'_strain_e11.fig']);
figure(4); colormap(black_rainbow_plus); caxis([-0.002 ,0.002 ]);
savefig(['fig_',num2str(LevelNo),'_strain_e12.fig']);
figure(5); colormap(black_rainbow_plus); caxis([-0.001 ,0.001 ]);
savefig(['fig_',num2str(LevelNo),'_strain_e22.fig']);
    
    
    
    
    
    
    
    
    
    
    
    
    
    %% ============ Level 2 --> Level 3 ==============
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ULevel1 = USubpb2Iter; FLevel1 = FSubpb2Iter; udualLevel1 = udualIter;  vdualLevel1 = vdualIter;
elementsFEM = DICmeshIter.elementsFEM;
coordinatesFEM = DICmeshIter.coordinatesFEM;
irregularEdge = DICmeshIter.irregular;

betavar = betaIter; muvar = muIter; alphavar = alphaIter;

CrackOrNot = 0; CrackPath1 = [0,0]; CrackPath2 = [0,0]; CrackTipOrNot = 0; CrackTip = [0,0]; 
EnrHAndTipEleIndex = []; EnrTipEleIndex = []; dirichlet = []; neumann= []; imgPyramidUnit = 1;

winstepsize = DICpara.winstepsize;

eleGeneration = ones(size(elementsFEM,1),1);
for eleInd = 1:size(elementsFEM,1)
    eleLength = abs(coordinatesFEM(elementsFEM(eleInd,1),1) - coordinatesFEM(elementsFEM(eleInd,3),1));
    eleGeneration(eleInd) = 1 + log2(winstepsize/eleLength);
end


%% ============== Estimate: Level 1 ==============
% ------ Compute rho to find where to refine ------
rhoLevel1Vector = zeros(size(elementsFEM,1),1); rhoLevel1Vectortemp1 = rhoLevel1Vector; rhoLevel1Vectortemp2 = rhoLevel1Vector;
ULevel1temp = [ULevel1;0;0]; hbar = parfor_progressbar(size(elementsFEM,1),'Please wait for a posterior error estimate!'); tic

for j = 1:size(elementsFEM,1)
    % eleNeighborIndexAndEdge = findEleNeighbors(elementsFEMIter, coordinatesFEMIter, 23 ,eleGenerationIter, winstepsize);
    point1x = coordinatesFEM(elementsFEM(j,1),1); point3x = coordinatesFEM(elementsFEM(j,3),1);lengthOfElement = (point3x-point1x)/winstepsize;
    % ---- Compute error estimator ----
    rhoLevel1Vectortemp1(j) = (lengthOfElement)^2 * aPostErrEstInt_ALSq(coordinatesFEM,elementsFEM,j,betavar,muvar,ULevel1,FLevel1,udualLevel1,vdualLevel1,imgPyramidUnit,CrackOrNot);
    rhoLevel1Vectortemp2(j) = (lengthOfElement) * (betavar+alphavar)^2 * aPostErrEstJump_ALSq(coordinatesFEM,elementsFEM,j,eleGeneration,winstepsize,ULevel1, ...
                               EnrHAndTipEleIndex,EnrTipEleIndex,CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot,DICmesh.markCoordHoleEdge);
    rhoLevel1Vector(j) = rhoLevel1Vectortemp1(j) + rhoLevel1Vectortemp2(j);
    hbar.iterate(1);
end
close(hbar); EstimateTime = toc;
Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(1,size(elementsFEM,1)); Sqctemp1 = Sqc; Sqctemp2 = Sqc;
for tempj = 1:size(elementsFEM,1)
    Sqx(1:4,tempj) = coordinatesFEM(elementsFEM(tempj,1:4),1);Sqy(1:4,tempj) = coordinatesFEM(elementsFEM(tempj,1:4),2); 
    Sqc(tempj) = rhoLevel1Vector(tempj);Sqctemp1(tempj) = rhoLevel1Vectortemp1(tempj);Sqctemp2(tempj) = rhoLevel1Vectortemp2(tempj);
end
close all; 
figure; patch(Sqx,Sqy,Sqc,'edgecolor','none' ); c=colorbar; title('Ele total err estimator','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal; % set(gca,'XTick',[]);
figure; patch(Sqx,Sqy,Sqctemp1,'edgecolor','none');c=colorbar; title('Ele interior err estimator','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal; % set(gca,'XTick',[]);
figure; patch(Sqx,Sqy,Sqctemp2,'edgecolor','none');c=colorbar; title('Ele jump err estimator','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal; % set(gca,'XTick',[]);
        
figure(1); set(gcf,'color','w'); savefig(['fig_L',num2str(LevelNo),'_err_all.fig']);
figure(2); set(gcf,'color','w'); savefig(['fig_L',num2str(LevelNo),'_err_i.fig']);
figure(3); set(gcf,'color','w'); savefig(['fig_L',num2str(LevelNo),'_err_j.fig']);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DorflerSatisfied = 1; 
while DorflerSatisfied == 1
    fprintf('--- Dorfler strategy to mark elements to refine? --- \n')
    fprintf('Input a real number between 0 and 1 here: ');
    prompt = 'Input here: '; thetaDorfler = input(prompt);
     
    %% ================ Mark: Level1 ==================
    % ======= Dorfler's Strategy to Mark elements to refine =======
    tic; [rhoLevel1VectorSorted,rhoLevel1VectorSortIndex] = sort(rhoLevel1Vector,'descend');
    rhoLevel1SumM = 0; rhoLevel1SumAll = sum(rhoLevel1Vector(:)); rhoLevel1SumMIndex = 0;
    while rhoLevel1SumM < 0.999999*thetaDorfler*rhoLevel1SumAll
        rhoLevel1SumMIndex = rhoLevel1SumMIndex + 1;
        rhoLevel1SumM = rhoLevel1SumM + rhoLevel1VectorSorted(rhoLevel1SumMIndex);
    end
    refinementPatch = rhoLevel1VectorSortIndex(1:rhoLevel1SumMIndex);
    if CrackOrNot > 0
        refinementPatch = union(refinementPatch,[EleCrackTop(end)+1:1:EleCrackBottom(1)-1]');
    end
    Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(1,size(elementsFEM,1));
    for j = 1:size(elementsFEM,1); Sqx(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),1); Sqy(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),2); end
    for j = 1:size(refinementPatch,1); Sqc(refinementPatch(j)) = 1; end
    MarkTime = toc;
    figure; patch(Sqx,Sqy,Sqc,'EdgeColor','none'); caxis([0,1]); colorbar;
    title('Marked ele (yellow)','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal;  
    set(gcf,'color','w');

    savefig(['fig_L',num2str(LevelNo+1),'_mark.fig']); % JY!!!
 
    
    %% ============== Refine: Level1 ==============
    if thetaDorfler>0 && size(coordinatesFEM,1)>1e3  
    
        tic; [coordinatesFEMLevel2,elementsFEMLevel2,irregularEdgeLevel2,dirichletLevel2,neumannLevel2] = ...
            QrefineRwBCNormal(coordinatesFEM,elementsFEM(:,1:4),irregularEdge,dirichlet,[],refinementPatch);
        elementsFEMLevel2 = [elementsFEMLevel2,zeros(size(elementsFEMLevel2,1),4)];
        
        % Assign element generation #
        eleGenerationLevel2 = ones(size(elementsFEMLevel2,1),1);
        for tempj=1:size(elementsFEMLevel2,1)
            eleGenerationLevel2(tempj) = 1-log2(norm(coordinatesFEMLevel2(elementsFEMLevel2(tempj,1),:) - coordinatesFEMLevel2(elementsFEMLevel2(tempj,2),:))/winstepsize);
        end
         
        % Reorder element nodes
        for tempj=1:length(eleGenerationLevel2)
            [~,row2] = max(  sum(coordinatesFEMLevel2(elementsFEMLevel2(tempj,1:4)',:)')  );
            if row2==1
                elementsFEMLevel2(tempj,1:4) = elementsFEMLevel2(tempj,[3,4,1,2]);
            elseif row2==2
                elementsFEMLevel2(tempj,1:4) = elementsFEMLevel2(tempj,[4,1,2,3]);
            elseif row2==4
                elementsFEMLevel2(tempj,1:4) = elementsFEMLevel2(tempj,[2,3,4,1]);
            end
        end
        
        % For hanging nodes in FE-mesh
        for tempj=1:size(irregularEdgeLevel2,1)
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,1:2), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,8)=irregularEdgeLevel2(tempj,3); end
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,2:3), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,5)=irregularEdgeLevel2(tempj,3); end
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,3:4), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,6)=irregularEdgeLevel2(tempj,3); end
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,[4,1]), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,7)=irregularEdgeLevel2(tempj,3); end
            
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,[2,1]), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,8)=irregularEdgeLevel2(tempj,3); end
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,[3,2]), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,5)=irregularEdgeLevel2(tempj,3); end
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,[4,3]), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,6)=irregularEdgeLevel2(tempj,3); end
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,[1,4]), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,7)=irregularEdgeLevel2(tempj,3); end
        end
        RefineTime = toc;
        
        % Initialize variables F0_Level2 and U0_Level2
        F_dispu = scatteredInterpolant( coordinatesFEM(:,1),coordinatesFEM(:,2),ULevel1(1:2:end) );
        F_dispv = scatteredInterpolant( coordinatesFEM(:,1),coordinatesFEM(:,2),ULevel1(2:2:end) );
        F_F11 = scatteredInterpolant( coordinatesFEM(:,1),coordinatesFEM(:,2),FLevel1(1:4:end) );
        F_F21 = scatteredInterpolant( coordinatesFEM(:,1),coordinatesFEM(:,2),FLevel1(2:4:end) );
        F_F12 = scatteredInterpolant( coordinatesFEM(:,1),coordinatesFEM(:,2),FLevel1(3:4:end) );
        F_F22 = scatteredInterpolant( coordinatesFEM(:,1),coordinatesFEM(:,2),FLevel1(4:4:end) );
        
        U0Level2 = 0*coordinatesFEMLevel2(:); F0Level2=0*[U0Level2(:);U0Level2(:)];
        temp = F_dispu(coordinatesFEMLevel2(:,1),coordinatesFEMLevel2(:,2)); U0Level2(1:2:end)=temp(:);
        temp = F_dispv(coordinatesFEMLevel2(:,1),coordinatesFEMLevel2(:,2)); U0Level2(2:2:end)=temp(:);
        temp = F_F11(coordinatesFEMLevel2(:,1),coordinatesFEMLevel2(:,2)); F0Level2(1:4:end)=temp(:);
        temp = F_F21(coordinatesFEMLevel2(:,1),coordinatesFEMLevel2(:,2)); F0Level2(2:4:end)=temp(:);
        temp = F_F12(coordinatesFEMLevel2(:,1),coordinatesFEMLevel2(:,2)); F0Level2(3:4:end)=temp(:);
        temp = F_F22(coordinatesFEMLevel2(:,1),coordinatesFEMLevel2(:,2)); F0Level2(4:4:end)=temp(:);
     
    else  
             
        tic; refinedEleIDListLevel1 = 0;
        coordinatesFEMLevel2 = coordinatesFEM; elementsFEMLevel2 = elementsFEM; U0Level2 = ULevel1; F0Level2 = FLevel1;
        dirichletLevel2 = unique(dirichlet(:)); neumannLevel2 = neumann; eleGenerationLevel2 = eleGeneration;
        h=waitbar(0,'Wait for mesh refinement.');
        %tempNeuMat = [];
        for j =  1:size(refinementPatch,1)
            checkAlreadyRefinedOrNot = ismember(refinementPatch(j),refinedEleIDListLevel1);
            if checkAlreadyRefinedOrNot == 0
                FinishRefineRecursive = 1;
                [coordinatesFEMLevel2,elementsFEMLevel2,U0Level2,F0Level2,refinedEleIDListNew,dirichletNew,neumannNew,eleGenerationLevel2] = ...
                    refineRecursiveSq(coordinatesFEMLevel2,elementsFEMLevel2,eleGenerationLevel2,U0Level2,F0Level2,refinementPatch(j),winstepsize,FinishRefineRecursive,...
                    CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot);
                
                refinedEleIDListLevel1 = unique([refinedEleIDListLevel1;refinedEleIDListNew]);
                dirichletLevel2 = unique([dirichletLevel2;dirichletNew']);
                % this step is wrong. newmann size is L*4:
                if size(neumannNew,1)>1
                    %tempNeuMat = [tempNeuMat; refinementPatch(j)];
                    [rowtemp1,~] = find( (neumannLevel2(:,1) == neumannNew(2,1))  );
                    [rowtemp2,~] = find( (neumannLevel2(:,2) == neumannNew(2,2))  );
                    [rowtemp3,~] = find( (neumannLevel2(:,2) == neumannNew(2,1))  );
                    [rowtemp4,~] = find( (neumannLevel2(:,1) == neumannNew(2,2))  );
                    rowtemp = union(intersect(rowtemp1,rowtemp2), intersect(rowtemp3,rowtemp4))  ;
                    if neumannLevel2(rowtemp,1) ~= neumannNew(2,3)
                        LengthOfNeumanntemp = size(neumannLevel2,1);
                        neumannLevel2(LengthOfNeumanntemp+1,1:4) = [neumannLevel2(rowtemp,1),neumannNew(2,3),neumannLevel2(rowtemp,3:4)];
                        neumannLevel2(rowtemp,1) = [neumannNew(2,3)];
                    end
                end
                % neumannLevel2 = unique([neumannLevel2;neumannNew']);
            end
            waitbar(j/size(refinementPatch,1));
        end
        RefineTime = toc; close(h); dirichletLevel2 = dirichletLevel2(2:end);
        if CrackOrNot > 0
            [CoordCrackTopNextIter,~,CoordCrackBottomNextIter,~] = funCrackTopOrBottom(coordinatesFEMLevel2,elementsFEMLevel2,CrackPath1,CrackPath2,CrackTip);
            F0Level2 = funSmoothStrainCrack(F0Level2,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,CoordCrackTopNextIter,CoordCrackBottomNextIter,LevelNo+1,0,0);
        end
        
        % Update irregularEdgeLevel2
        irregularEdgeLevel2 = [0,0,0];
        [row,~]=find(elementsFEMLevel2(:,5)>0);
        for tempi=1:length(row)
            irregularEdgeLevel2 = [irregularEdgeLevel2; elementsFEMLevel2(row(tempi),2:3),  elementsFEMLevel2(row(tempi),5) ];
        end
        [row,~]=find(elementsFEMLevel2(:,6)>0);
        for tempi=1:length(row)
            irregularEdgeLevel2 = [irregularEdgeLevel2; elementsFEMLevel2(row(tempi),3:4),  elementsFEMLevel2(row(tempi),6) ];
        end
        [row,~]=find(elementsFEMLevel2(:,7)>0);
        for tempi=1:length(row)
            irregularEdgeLevel2 = [irregularEdgeLevel2; elementsFEMLevel2(row(tempi),[4,1]),  elementsFEMLevel2(row(tempi),7) ];
        end
        [row,~]=find(elementsFEMLevel2(:,8)>0);
        for tempi=1:length(row)
            irregularEdgeLevel2 = [irregularEdgeLevel2; elementsFEMLevel2(row(tempi),1:2),  elementsFEMLevel2(row(tempi),8) ];
        end
        
        dirichletLevel2 = [];
         
    end % END of mesh refinement
    
    

    % ====== Show refined mesh ======
    Sqx = zeros(4,size(elementsFEMLevel2,1)); Sqy = zeros(4,size(elementsFEMLevel2,1)); Sqc = zeros(1,size(elementsFEMLevel2,1));
    for j = 1:size(elementsFEMLevel2,1)
        Sqx(1:4,j) = coordinatesFEMLevel2(elementsFEMLevel2(j,1:4),1);
        Sqy(1:4,j) = coordinatesFEMLevel2(elementsFEMLevel2(j,1:4),2); Sqc(j) =  1;
    end
    figure; patch(Sqx,Sqy,Sqc,'FaceColor','none');  axis equal;  title(['Level ',num2str(LevelNo+1),' Mesh'],'fontweight','normal'); 
    set(gca,'fontsize',16); set(gcf,'color','w'); axis on; axis tight;

    savefig(['fig_L',num2str(LevelNo+1),'_refine.fig']); % JY!!!
     
    fprintf('--- Are you satisfied with current Dorfler refinement number (0-yes; 1-no)? ---  \n')
    prompt = 'Input here: ';
    DorflerSatisfied = input(prompt);
    if DorflerSatisfied == 0
        fprintf('--- Mesh refinement done! ---  \n')
    end

    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LevelNo = LevelNo + 1; 
coordinatesFEMIter = coordinatesFEMLevel2;  
coordinatesFEMWorldIter = coordinatesFEMIter; coordinatesFEMWorldIter(:,2) = DICpara.ImgSize(2)+1-coordinatesFEMWorldIter(:,2);
elementsFEMIter = elementsFEMLevel2;
irregularEdgeIter = irregularEdgeLevel2;  

DICmeshIter.coordinatesFEM = coordinatesFEMLevel2;
DICmeshIter.elementsFEM = elementsFEMLevel2;
DICmeshIter.dirichlet = []; DICmeshIter.neumann = [];
DICmeshIter.coordinatesFEMWorld = coordinatesFEMWorldIter;
DICmeshIter.elementMinSize = DICmesh.elementMinSize;
DICmeshIter.irregular = irregularEdgeLevel2;

U0Iter = U0Level2; ConvItPerEleIter = zeros(size(coordinatesFEMIter,1),6);



    %% Section 4: ALDIC Subproblem 1 -or- Local ICGN Subset DIC
    fprintf('------------ Section 4 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve the first local step in ALDIC: Subproblem 1
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ====== ALStep 1 Subproblem1: Local Subset DIC ======
    mu=0; beta=0; tol=1e-2; ALSolveStep=1; ALSub1Time=zeros(6,1); ALSub2Time=zeros(6,1);
    ConvItPerEleIter=zeros(size(coordinatesFEMIter,1),6); ALSub1BadPtNum=zeros(6,1);
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ Start Local DIC IC-GN iteration ------
    [USubpb1Iter,FSubpb1Iter,HtempPar,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp,markCoordHoleStrain] = ...
        LocalICGNQuadtree(U0Iter,coordinatesFEMIter,Df,fNormalized,gNormalized,DICpara,'GaussNewton',tol);
    ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEleIter(:,ALSolveStep) = ConvItPerEletemp; ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp; toc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------  Manually find some bad points from Local Subset ICGN step ------
    % Comment these lines below if you don't have local bad points
    % %%%%% Comment START %%%%%%
    %    [USubpb1,FSubpb1] = funRemoveOutliersQuadtree(DICmesh,DICpara,USubpb1,FSubpb1);
    %    disp('--- Remove bad points done ---')
    % %%%%% Comment END %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ------ Plot ------
    USubpb1WorldIter = USubpb1Iter; USubpb1WorldIter(2:2:end) = -USubpb1Iter(2:2:end); 
    FSubpb1WorldIter = FSubpb1Iter; FSubpb1WorldIter(2:4:end) = -FSubpb1WorldIter(2:4:end); FSubpb1WorldIter(3:4:end) = -FSubpb1WorldIter(3:4:end); 
    
    close all; Plotdisp_show(USubpb1WorldIter, DICmeshIter.coordinatesFEMWorld, DICmeshIter.elementsFEM(:,1:4),DICpara,'NoEdgeColor');
    Plotstrain_show(FSubpb1WorldIter, DICmeshIter.coordinatesFEMWorld, DICmeshIter.elementsFEM(:,1:4),DICpara,'NoEdgeColor');
    save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
    fprintf('------------ Section 4 Done ------------ \n \n')
 
   
    %% Section 5: Subproblem 2 -- solve the global compatible displacement field
    fprintf('------------ Section 5 Start ------------ \n'); tic;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve the global step in ALDIC Subproblem 2
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ======= ALStep 1 Subproblem 2: Global constraint =======
    % ------ Smooth displacements for a better F ------
    DICpara.DispFilterSize=0; DICpara.DispFilterStd=0; DICpara.StrainFilterSize=0; DICpara.StrainFilterStd=0;  
    DICpara.DispSmoothness = 0; DICpara.StrainSmoothness = 1e-4;
    if DICpara.DispSmoothness>1e-6, USubpb1Iter = funSmoothDispQuadtreeRBF(USubpb1Iter,DICmeshIter,DICpara); end %%%%% should be Iter
    if DICpara.StrainSmoothness>1e-6, FSubpb1Iter = funSmoothStrainQuadtreeRBF(FSubpb1Iter,DICmeshIter,DICpara); end  %%%%% should be Iter
    
	% ====== Define penalty parameter ======
    muIter = 1e-3; udualIter = 0*FSubpb1Iter; vdualIter = 0*USubpb1Iter; 
    betaListIter = [1e-3,1e-2,1e-1]*mean(DICpara.winstepsize).^2.*muIter; % Tune beta in the betaList 
    Err1 = zeros(length(betaList),1); Err2 = Err1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****']);
    DICpara.GaussPtOrder = 2; alpha = 0;  % No regularization added
    % ====== Solver using finite element method ======
    % if ImgSeqNum == 2
        for tempk = 1:length(betaListIter)
            betaIter = betaListIter(tempk); display(['Try #',num2str(tempk),' beta = ',num2str(betaIter)]);
            alphaIter = 100*betaIter;
            GaussPtOrder=3; alpha=0; [USubpb2Iter] = Subpb2Quadtree(DICmeshIter,DICpara.GaussPtOrder, ...
                betaIter,muIter,USubpb1Iter,FSubpb1Iter,udualIter,vdualIter,alphaIter,mean(DICpara.winstepsize),0);
            % [FSubpb2Iter,~,~] = funGlobalNodalStrainQuadtree(DICmeshIter,USubpb2Iter,DICpara.GaussPtOrder,0);
           %  [FSubpb2Iter] = funGlobalNodalStrainRBF(DICmeshIter,DICpara,USubpb2Iter);
            [~,~,FSubpb2Iter] = funCompDefGrad2([USubpb2Iter(1:2:end),USubpb2Iter(2:2:end)], ...
                DICmeshIter.coordinatesFEM, DICpara.winstepsize, 20, DICpara.ImgRefMask);
            
            Err1(tempk) = norm(USubpb1Iter-USubpb2Iter,2);
            Err2(tempk) = norm(FSubpb1Iter-FSubpb2Iter,2);
        end
        
        Err1Norm = (Err1-mean(Err1))/std(Err1); % figure, plot(Err1Norm);
        Err2Norm = (Err2-mean(Err2))/std(Err2); % figure, plot(Err2Norm);
        ErrSum = Err1Norm+Err2Norm; % figure, plot(ErrSum); title('Tune the best \beta in the subproblem 2'); 
        [~,indexOfbeta] = min(ErrSum); 
     
        try % Tune the best beta by a quadratic polynomial 0fitting
            [fitobj] = fit(log10(betaList(indexOfbeta-1:1:indexOfbeta+1))',ErrSum(indexOfbeta-1:1:indexOfbeta+1),'poly2');
            p = coeffvalues(fitobj); betaIter = 10^(-p(2)/2/p(1));
        catch, betaIter = betaListIter(indexOfbeta);
        end
        display(['Best beta = ',num2str(betaIter)]);
    % else 
    %    try betaIter = DICpara.beta;
    %    catch, betaIter = 1e-3*mean(DICpara.winstepsize).^2.*muIter;
    %    end
    %end
      
    % Using the optimal beta to solve the ALDIC Subproblem 2 again
    alphaIter = 100*betaIter;
    if abs(betaIter-betaListIter(end))>abs(eps)
        [USubpb2Iter] = Subpb2Quadtree(DICmeshIter,DICpara.GaussPtOrder,betaIter,muIter,USubpb1Iter,FSubpb1Iter,udualIter,vdualIter,alphaIter,mean(DICpara.winstepsize),0);
        % [FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2,DICpara.GaussPtOrder,0);
        % [FSubpb2Iter] = funGlobalNodalStrainRBF(DICmeshIter,DICpara,USubpb2Iter);
        [~,~,FSubpb2Iter] = funCompDefGrad2([USubpb2Iter(1:2:end),USubpb2Iter(2:2:end)], ...
                 DICmeshIter.coordinatesFEM, DICpara.winstepsize, 20, DICpara.ImgRefMask);
        
        ALSub2Time(ALSolveStep) = toc; toc
    end
    
    % ------- Smooth strain field --------
    if DICpara.DispSmoothness>1e-6, USubpb2Iter = funSmoothDispQuadtreeRBF(USubpb2Iter,DICmeshIter,DICpara); end
    % ------- Don't smooth strain fields near the boundary --------
    for tempk=0:3, FSubpb2Iter(4*DICmesh.markCoordHoleEdge-tempk) = FSubpb1Iter(4*DICmesh.markCoordHoleEdge-tempk); end
    if DICpara.StrainSmoothness>1e-6, FSubpb2Iter = funSmoothStrainQuadtreeRBF(0.1*FSubpb2Iter+0.9*FSubpb1Iter,DICmeshIter,DICpara); end
    for tempk=0:1, USubpb2Iter(2*markCoordHoleStrain-tempk) = USubpb1Iter(2*markCoordHoleStrain-tempk); end
    for tempk=0:3, FSubpb2Iter(4*markCoordHoleStrain-tempk) = FSubpb1Iter(4*markCoordHoleStrain-tempk); end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ------- Save data ------
    save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2');

    % ------ Plot ------
    USubpb2WorldIter = USubpb2Iter; USubpb2WorldIter(2:2:end) = -USubpb2Iter(2:2:end);
    FSubpb2WorldIter = FSubpb2Iter; FSubpb2WorldIter(2:4:end) = -FSubpb2WorldIter(2:4:end); FSubpb2WorldIter(3:4:end) = -FSubpb2WorldIter(3:4:end); 
    close all; Plotdisp_show(USubpb2WorldIter, DICmeshIter.coordinatesFEMWorld, DICmeshIter.elementsFEM(:,1:4),DICpara,'EdgeColor');
    Plotstrain_show(FSubpb2WorldIter, DICmeshIter.coordinatesFEMWorld, DICmeshIter.elementsFEM(:,1:4),DICpara,'EdgeColor');
    
    % ======= Update dual variables =======
    udualIter = FSubpb2Iter - FSubpb1Iter; vdualIter = USubpb2Iter - USubpb1Iter;
	save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');
    fprintf('------------ Section 5 Done ------------ \n \n')
    
    
    
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
         
        %%%%% Otherwise, we use the fixed DIC subset size %%%%%
        winsize_List = DICpara.winsize*ones(size(DICmeshIter.coordinatesFEM,1),2);
        DICpara.winsize_List = winsize_List;
        

        %%%%%%%%%%%%%%%%%%%%%%% Subproblem 1 %%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****']);
        tic; [USubpb1Iter,~,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp] = Subpb1Quadtree(...
                                            USubpb2Iter,FSubpb2Iter,udualIter,vdualIter,DICmeshIter.coordinatesFEM,...
                                            Df,fNormalized,gNormalized,muIter,betaIter,HPar,ALSolveStep,DICpara,'GaussNewton',tol);
        FSubpb1Iter = FSubpb2Iter; toc 
        for tempk=0:1, USubpb1Iter(2*markCoordHoleStrain-tempk) = USubpb2Iter(2*markCoordHoleStrain-tempk); end
        ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEleIter(:,ALSolveStep) = ConvItPerEletemp; ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp;
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
        tic; [USubpb2Iter] = Subpb2Quadtree(DICmeshIter,DICpara.GaussPtOrder,betaIter,muIter, ...
            USubpb1Iter,FSubpb1Iter,udualIter,vdualIter,alphaIter,mean(DICpara.winstepsize),0);
		% [FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2,DICpara.GaussPtOrder,0);
        % [FSubpb2Iter] = funGlobalNodalStrainRBF(DICmeshIter,DICpara,USubpb2Iter);
        [~,~,FSubpb2Iter] = funCompDefGrad2([USubpb2Iter(1:2:end),USubpb2Iter(2:2:end)], ...
                 DICmeshIter.coordinatesFEM, DICpara.winstepsize,20, DICpara.ImgRefMask);
        ALSub2Time(ALSolveStep) = toc; toc
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ------- Smooth strain field --------
        if DICpara.DispSmoothness>1e-6, USubpb2Iter = funSmoothDispQuadtreeRBF(USubpb2Iter,DICmeshIter,DICpara); end
        % ------- Don't change strain fields near the boundary --------
        for tempk=0:3, FSubpb2Iter(4*DICmesh.markCoordHoleEdge-tempk) = FSubpb1Iter(4*DICmesh.markCoordHoleEdge-tempk); end
        if DICpara.StrainSmoothness>1e-6, FSubpb2Iter = funSmoothStrainQuadtreeRBF(0.1*FSubpb2Iter+0.9*FSubpb1Iter,DICmeshIter,DICpara); end
        for tempk=0:1, USubpb2Iter(2*markCoordHoleStrain-tempk) = USubpb1Iter(2*markCoordHoleStrain-tempk); end
        for tempk=0:3, FSubpb2Iter(4*markCoordHoleStrain-tempk) = FSubpb1Iter(4*markCoordHoleStrain-tempk); end
         
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
        udualIter = FSubpb2Iter - FSubpb1Iter; vdualIter = USubpb2Iter - USubpb1Iter; 
		 
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
    % ResultDisp{ImgSeqNum-1}.U = full(USubpb2);
    % ResultDisp{ImgSeqNum-1}.ALSub1BadPtNum = ALSub1BadPtNum;
    % ResultDefGrad{ImgSeqNum-1}.F = full(FSubpb2);  
    % close all; Plotdisp_show(USubpb2Iter,DICmeshIter.coordinatesFEM,DICmeshIter.elementsFEM(:,1:4),DICpara,'NoEdgeColor');
    % Plotstrain_show(FSubpb2Iter,DICmeshIter.coordinatesFEM,DICmeshIter.elementsFEM(:,1:4),DICpara,'NoEdgeColor');
    
    % ------ Plot ------
    USubpb2WorldIter = USubpb2Iter; USubpb2WorldIter(2:2:end) = -USubpb2Iter(2:2:end);
    FSubpb2WorldIter = FSubpb2Iter; FSubpb2WorldIter(2:4:end) = -FSubpb2WorldIter(2:4:end); FSubpb2WorldIter(3:4:end) = -FSubpb2WorldIter(3:4:end);
    close all; Plotdisp_show(USubpb2WorldIter,DICmeshIter.coordinatesFEMWorld,DICmeshIter.elementsFEM(:,1:4),DICpara,'NoEdgeColor');
    Plotstrain_show(FSubpb2WorldIter,DICmeshIter.coordinatesFEMWorld,DICmeshIter.elementsFEM(:,1:4),DICpara,'NoEdgeColor');
    
    figure(1); colormap(black_rainbow_plus); caxis([-0.5,0.5]);
savefig(['fig_L',num2str(LevelNo),'_dispx.fig']);
figure(2); colormap(black_rainbow_plus); caxis([-0.1,0.1]);
savefig(['fig_L',num2str(LevelNo),'_dispy.fig']);
figure(3); colormap(black_rainbow_plus); caxis([-0.004,0.004]);
savefig(['fig_L',num2str(LevelNo),'_strain_e11.fig']);
figure(4); colormap(black_rainbow_plus); caxis([-0.002 ,0.002 ]);
savefig(['fig_L',num2str(LevelNo),'_strain_e12.fig']);
figure(5); colormap(black_rainbow_plus); caxis([-0.001 ,0.001 ]);
savefig(['fig_L',num2str(LevelNo),'_strain_e22.fig']);
    





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
DICpara.DoYouWantToSmoothOnceMore = 0; % No need to smooth disp fields
DICpara.smoothness = funParaInput('RegularizationSmoothness'); % Regularization to smooth strain fields           
% ------ Choose strain computation method ------
DICpara.MethodToComputeStrain = 2; % funParaInput('StrainMethodOp'); 
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

%% ------ Start main part ------
for ImgSeqNum = [2 : length(ImgNormalized)] 
    
    close all; 
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
    if DICpara.ImgSeqIncUnit > 1
        FEMeshIndLast = floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit);
    elseif DICpara.ImgSeqIncUnit == 1
        FEMeshIndLast = floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)-1;
    end
    FEMeshInd = FEMeshIndLast + 1;
    
    if FEMeshInd == 1
        USubpb2 = ResultDisp{ImgSeqNum-1}.U; %+ ResultDisp{10}.U + ResultDisp{20}.U;
        coordinatesFEM = ResultFEMesh{1}.coordinatesFEM; 
        elementsFEM = ResultFEMesh{1}.elementsFEM;
        if (ImgSeqNum-1 == 1) || (DICpara.ImgSeqIncROIUpdateOrNot==1), UFEMesh = 0*USubpb2; end
    else
        USubpb2 = ResultDisp{ImgSeqNum-1}.U;
        if mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) == 0
            coordinatesFEM = ResultFEMesh{FEMeshInd}.coordinatesFEM;
            elementsFEM = ResultFEMesh{FEMeshInd}.elementsFEM;
            coordinatesFEMLast = ResultFEMesh{FEMeshIndLast}.coordinatesFEM;
            UFEMeshLast = ResultDisp{ImgSeqNum-2}.U + UFEMesh;
            xq = coordinatesFEM(:,1); yq = coordinatesFEM(:,2);
            UFEMesh = 0*USubpb2;
            UFEMesh(1:2:end) = griddata(coordinatesFEMLast(:,1),coordinatesFEMLast(:,2),UFEMeshLast(1:2:end),xq,yq,'v4');
            UFEMesh(2:2:end) = griddata(coordinatesFEMLast(:,1),coordinatesFEMLast(:,2),UFEMeshLast(2:2:end),xq,yq,'v4');
        end
        USubpb2 = USubpb2 + UFEMesh;
    end

    FSubpb2 = ResultDefGrad{ImgSeqNum-1}.F;
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
    UWorld = DICpara.um2px*ULocal; UWorld(2:2:end) = -UWorld(2:2:end); FWorld = FLocal; 
     
    % ------ Smooth displacements ------
    %prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
    %DoYouWantToSmoothOnceMore = input(prompt); 
    SmoothTimes = 0;
    try
        while DICpara.DoYouWantToSmoothOnceMore == 0 && SmoothTimes < 3
            ULocal = funSmoothDispQuadtreeRBF(ULocal,DICmesh,DICpara);
            %close all; Plotuv(ULocal,x0,y0); %DICpara.DoYouWantToSmoothOnceMore = input(prompt);
            SmoothTimes = SmoothTimes + 1;
        end
    catch
    end
    
    % ----- Compute strain field ------
    %%%%% Load deformed image %%%%%%%
    fNormalized = ImgNormalized{ 1 } ; % Load current deformed image frame 
    fNormalizedMask = double( ImgMask{1} );
    fNormalized = fNormalized .* fNormalizedMask ;
    Dg = funImgGradient(fNormalized,fNormalized,fNormalizedMask); % Finite difference to compute image grayscale gradients;
    
    ComputeStrainQuadtree;  
     
    % ------ Plot disp and strain ------
    if DICpara.OrigDICImgTransparency == 1
        Plotdisp_showMasks(UWorld,coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),ImgMask{1},DICpara,'NoEdgeColor');
        [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min,strain_maxshear,strain_vonMises] = ...
                   Plotstrain0Quadtree(FStrainWorld,coordinatesFEMWorld,elementsFEM(:,1:4),DICpara);
    
    else % Plot over raw DIC images
        if DICpara.Image2PlotResults == 0 % Plot over the first image; "file_name{1,1}" corresponds to the first image	
             
            PlotdispQuadtreeMasks(UWorld,coordinatesFEMWorld,elementsFEM(:,1:4),...
                file_name{1,1}, ImgMask{1},DICpara);
            
            [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
               strain_maxshear,strain_vonMises] = PlotstrainQuadtreeMasks(UWorld,FStrainWorld, ...
               coordinatesFEMWorld,elementsFEM(:,1:4),file_name{1,1}, ImgMask{1},DICpara);
           
           
        else % Plot over second or next deformed images
             
             PlotdispQuadtreeMasks(UWorld,coordinatesFEMWorld,elementsFEM(:,1:4),...
                file_name{1, ImgSeqNum}, [], DICpara);
              
            [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
               strain_maxshear,strain_vonMises] = PlotstrainQuadtreeMasks(UWorld,FStrainWorld, ...
               coordinatesFEMWorld,elementsFEM(:,1:4),file_name{1, ImgSeqNum}, ...
                [],DICpara);
             
        end
    end
 
    % ----- Save strain results ------
    ResultStrain{ImgSeqNum-1} = struct('strainxCoord',coordinatesFEMWorld(:,1),'strainyCoord',coordinatesFEMWorld(:,2), ...
            'dispu',UWorld(1:2:end),'dispv',UWorld(2:2:end), ...
            'dudx',FStraintemp(1:4:end),'dvdx',FStraintemp(2:4:end),'dudy',FStraintemp(3:4:end),'dvdy',FStraintemp(4:4:end), ...
            'strain_exx',strain_exx,'strain_exy',strain_exy,'strain_eyy',strain_eyy, ...
            'strain_principal_max',strain_principal_max,'strain_principal_min',strain_principal_min, ...
            'strain_maxshear',strain_maxshear,'strain_vonMises',strain_vonMises);
      
    % ------ Save figures for tracked displacement and strain fields ------
    SaveFigFilesDispAndStrainQuadtree;
    % pause;
    
    
end
% ------ END of for-loop {ImgSeqNum = 2:length(ImgNormalized)} ------
fprintf('------------ Section 8 Done ------------ \n \n')


% ------ Save data again including solved strain fields ------
results_name = ['results_',imgname,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultFEMesh','ResultFEMeshEachFrame', ...
                   'ALSub1Time','ALSub2Time','ALSolveStep','ResultStrain');

 


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
v = VideoWriter('video_mesh.avi');
v.FrameRate = 10;
open(v);
figure,
for ImgSeqNum = 81 % : (1+size(ResultDisp,1))
    
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


%%



%% ------ Start main part for velocity field ------
for ImgSeqNum = [3 : length(ImgNormalized)-1] 
    
    close all; 
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
    if DICpara.ImgSeqIncUnit > 1
        FEMeshIndLast = floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit);
    elseif DICpara.ImgSeqIncUnit == 1
        FEMeshIndLast = floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)-1;
    end
    FEMeshInd = FEMeshIndLast + 1;
    
    if FEMeshInd == 1
        USubpb2 = ResultDisp{ImgSeqNum-1}.U; %+ ResultDisp{10}.U + ResultDisp{20}.U;
        coordinatesFEM = ResultFEMesh{1}.coordinatesFEM; 
        elementsFEM = ResultFEMesh{1}.elementsFEM;
        if (ImgSeqNum-1 == 1) || (DICpara.ImgSeqIncROIUpdateOrNot==1), UFEMesh = 0*USubpb2; end
    else
        USubpb2 = ResultDisp{ImgSeqNum-1}.U;
        if mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) == 0
            coordinatesFEM = ResultFEMesh{FEMeshInd}.coordinatesFEM;
            elementsFEM = ResultFEMesh{FEMeshInd}.elementsFEM;
            coordinatesFEMLast = ResultFEMesh{FEMeshIndLast}.coordinatesFEM;
            UFEMeshLast = ResultDisp{ImgSeqNum-2}.U + UFEMesh;
            xq = coordinatesFEM(:,1); yq = coordinatesFEM(:,2);
            UFEMesh = 0*USubpb2;
            UFEMesh(1:2:end) = griddata(coordinatesFEMLast(:,1),coordinatesFEMLast(:,2),UFEMeshLast(1:2:end),xq,yq,'v4');
            UFEMesh(2:2:end) = griddata(coordinatesFEMLast(:,1),coordinatesFEMLast(:,2),UFEMeshLast(2:2:end),xq,yq,'v4');
        end
        USubpb2 = USubpb2 + UFEMesh;
    end

    USubpb2Prev = ResultDisp{ImgSeqNum-2}.U;
    USubpb2Next = ResultDisp{ImgSeqNum}.U;
    VSubpb2 = (USubpb2Next-USubpb2Prev)/2;
    %coordinatesFEM = ResultFEMeshEachFrame{ImgSeqNum-1}.coordinatesFEM;
    %elementsFEM = ResultFEMeshEachFrame{ImgSeqNum-1}.elementsFEM;
    try markCoordHoleEdge = ResultFEMeshEachFrame{ImgSeqNum-1}.markCoordHoleEdge; catch; end
    DICmesh.coordinatesFEM = coordinatesFEM;
    DICmesh.elementsFEM = elementsFEM;
    coordinatesFEMWorld = DICpara.um2px*[coordinatesFEM(:,1),size(ImgNormalized{1},2)+1-coordinatesFEM(:,2)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    VLocal = VSubpb2;
    VWorld = DICpara.um2px*VLocal; VWorld(2:2:end) = -VWorld(2:2:end);   % close all; Plotuv(UWorld,x0,y0World);
     
    % ------ Smooth displacements ------
    %prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
    %DoYouWantToSmoothOnceMore = input(prompt); 
    SmoothTimes = 0;
    try
        while DICpara.DoYouWantToSmoothOnceMore == 0 && SmoothTimes < 3
            VLocal = funSmoothDispQuadtreeRBF(VLocal,DICmesh,DICpara);
            %close all; Plotuv(ULocal,x0,y0); %DICpara.DoYouWantToSmoothOnceMore = input(prompt);
            SmoothTimes = SmoothTimes + 1;
        end
    catch
    end
     
    % ------ Plot disp and strain ------
    if DICpara.OrigDICImgTransparency == 1
        Plotdisp_showMasks(UWorld,coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),ImgMask{1},DICpara,'NoEdgeColor');
        
    else % Plot over raw DIC images
        if DICpara.Image2PlotResults == 0 % Plot over the first image; "file_name{1,1}" corresponds to the first image	
             
            PlotdispQuadtreeMasks(UWorld,coordinatesFEMWorld,elementsFEM(:,1:4),...
                file_name{1,1}, ImgMask{1},DICpara);
            
           
        else % Plot over second or next deformed images
             
             PlotdispQuadtreeMasks(VWorld,coordinatesFEMWorld,elementsFEM(:,1:4),...
                file_name{1, ImgSeqNum}, ImgMask{ 1 }, DICpara);
              
        end
    end
 
    
    % ------ Save figures for tracked displacement and strain fields ------
    % SaveFigFilesDispAndStrainQuadtree;
    % pause;
    
    
end
% ------ END of for-loop {ImgSeqNum = 2:length(ImgNormalized)} ------
fprintf('------------ Section 8 Done ------------ \n \n')


% ------ Save data again including solved strain fields ------
% results_name = ['results_',imgname,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
% save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultFEMesh',...
%                    'ALSub1Time','ALSub2Time','ALSolveStep','ResultStrain');
