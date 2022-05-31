classdef ALDIC
    % ---------------------------------------------
    % Augmented Lagrangian Digital Image Correlation (AL-DIC)
    % working with an adaptive quadtree mesh
    %
    % Author of algorithms: Jin Yang, PhD @Caltech
    % Author of GUI and adapted UI: Alexander McGhee, PhD @UF
    % Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
    % Date: 2015.04,06,07; 2016.03,04; 2020.11
    % ---------------------------------------------
    % The following functions are modified/adapted from the original
    % ReadImage -> moved to the Functions/ALDIC folder
    % IntegerSearch -> InitSearchSize
    
    properties
        Type = 'ALDIC';
        ProgramFolder = '2D_ALDIC-master';
        Name = '';
        SaveFolderPath = '';
        quadtree = false;
        
        % Raw Image Data
        Images = {};
        RefImages = [];
        % Raw ROI mask Data
        ROIMasks = {};
        
        % DIC parameters
        DICparawinsize = 50;
        DICparawinstepsize = [5,5];
        DICparagridx = [];
        DICparagridy = [];
        DICparaLoadImgMethod = 0;
        DICparaImgSeqIncUnit = 1;
        DICparaImgSeqIncROIUpdateOrNot = 0;
        DICparaImgSize = [];
        DICparaImgRefMask = [];
        DICparaSubpb2FDOrFEM = 0;
        DICparaNewFFTSearch = 0;
        DICparaClusterNo = [];
        DICparamu = 1e-3;
        DICparabeta = 0;
        DICparatol = 1e-2;
        DICparaDispFilterSize=0;
        DICparaDispFilterStd=0;
        DICparaStrainFilterSize=0;
        DICparaStrainFilterStd=0;
        DICparaInitFFTSearchMethod = 0;
        DICparaSizeOfFFTSearchRegion = 0;
        DICparaum2px = 1;
        DICparaDoYouWantToSmoothOnceMore = 0;
        DICparaMethodToComputeStrain = 0;
        DICparaStrainType = 0;
        DICparaImage2PlotResults = 0;
        DICparaMethodToSaveFig = 0;
        DICparaOrigDICImgTransparency = 0;
        % DIC mesh
        DICmesh_coordinatesFEM = [];
        DICmesh_elementsFEM = [];
        DICmesh_dirichlet = [];
        DICmesh_neumann = [];
        DICmesh_x0 = [];
        DICmesh_y0 = [];
        DICmesh_M = [];
        DICmesh_N = [];
        DICmesh_y0World = [];
        DICmesh_coordinatesFEMWorld = [];
        DICmesh_elementMinSize = 2;
        
        % DIC algo vars
        Df = [];
        U0 = [];
        FSubpb1 = {};
        USubpb1 = {};
        FSubpb2 = {};
        USubpb2 = {};
        udual = {};
        vdual = {};
        HtempPar = [];
        ALSub1Timetemp = [];
        ConvItPerEletemp = [];
        LocalICGNBadPtNumtemp = [];
        temp3 = [];
        temp4 = [];
        ALSolveStep = [];
        pts = [1,1];
        x0temp_f = [];
        y0temp_f = [];
        u_f = [];
        v_f = [];
        ccThreshold = [];
        ccmax = [];
        ccA = [];
        ccqfactors = [];
        u = [];
        v = [];
        Ux = [];
        Vy = [];
        GaussPtOrder = 2;
        alpha = 0;
        
        % convergence results
        Convergence_V = [];
        Convergence_W = [];
        Convergence_U = [];
        Convergence_F = [];
        Convergence_Y = [];
        
        % Result matricies
        ResultDisp = {};
        ResultDefGrad = {};
        ResultFEMeshEachFrame_coordinatesFEM = {};
        ResultFEMeshEachFrame_elementsFEM = {};
        ResultFEMesh_coordinatesFEM = {};
        ResultFEMesh_elementsFEM = {};
        ResultFEMesh_winsize = {};
        ResultFEMesh_winstepsize = {};
        ResultFEMesh_gridx = {};
        ResultFEMesh_gridy = {};
        ResultStrainDirect = {};
        ResultStrainFD = {};
        ResultStrainPlane = {};
        ResultStrainFE = {};
        ResultStrainWorld = {};
        ResultStressWorld = {};
        ResultStrainMap = {};
        ResultDispMapu = {};
        ResultDispMapv = {};
        ResultStrainRad = [];
        ResultStrainT3 = [];
    end
    
    methods
        % these are the main functions to use the ALDIC software
        function obj = ALDIC()
            %ALDIC Constructor
            
            % add the ALDIC functions to the path and remove others
            AddAndRemovePathFolders(obj,1)
            
            % initialize the minGW c++ engine
            obj.initMGW()
            
        end
        
        function obj = runALDIC(obj,ImgRef,ImgDef,ImgSeqNum)
            
            %step5.2 Compute an initial guess
            obj = InitialGuess(obj,ImgRef,ImgDef,ImgSeqNum);
            
            %step5.3 solve Subproblem 1
            obj = Subproblem1(obj,ImgRef,ImgDef);
            
            %step5.4 solve Subproblem 2
            obj = Subproblem2(obj);
            
            %step5.5 apply ADMM iterations
            obj = ADMMiteration(obj,ImgRef,ImgDef,ImgSeqNum);
            
            
        end
        
        function obj = endALDIC(obj)
            %clean up all the variables to only include the results and
            %other important variables.
            
        end
        
        function obj = save(obj,filename)
            save(filename,'obj')
        end
        
        % The initial guess steps FFT based cross correlation
        function obj = InitialGuess(obj,ImgRef,ImgDef,ImgSeqNum)
            % This function computes the FFT-based cross correlation on first 7 frames
            % and then uses the data driven method to estimate the initial guesses for other frames
            if ImgSeqNum == 2
                obj.DICparaNewFFTSearch = 1;
                obj.DICparaInitFFTSearchMethod = [];
            elseif ImgSeqNum < 7
                obj.DICparaNewFFTSearch = 1; % Use FFT-based cross correlation to compute the initial guess
            else
                obj.DICparaNewFFTSearch = 0; % Apply data driven method to estimate initial guesses for later frames
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % zip all the DIC options into a struct
            DICpara = zipDICPara(obj);
            % zip all the DIC mesh parameters into a struct
            DICmesh = zipDICmesh(obj);
            
            if ImgSeqNum == 2 || obj.DICparaNewFFTSearch == 1 % Apply FFT-based cross correlation to compute the initial guess
                
                % ====== FFT-based cross correlation ======
                % Switch the order of Img f and Img g
                [obj,~,obj.x0temp_f,obj.y0temp_f,obj.u_f,obj.v_f,cc] = InitSearchSize(obj,ImgRef,ImgDef);
                obj = unzipcc(obj,cc);
                
                % ====== FEM mesh set up ======
                xnodes = max([1+0.5*obj.DICparawinsize+ obj.DICparaSizeOfFFTSearchRegion(1), obj.DICparagridx(1) ])  ...
                    : obj.DICparawinstepsize : min([size(ImgRef,1)-0.5*obj.DICparawinsize-1- obj.DICparaSizeOfFFTSearchRegion(1),obj.DICparagridx(2) ]);
                ynodes = max([1+0.5*obj.DICparawinsize+ obj.DICparaSizeOfFFTSearchRegion(2),obj.DICparagridy(1) ])  ...
                    : obj.DICparawinstepsize : min([size(ImgRef,2)-0.5*obj.DICparawinsize-1- obj.DICparaSizeOfFFTSearchRegion(2),obj.DICparagridy(2) ]);
                
                [x0temp,y0temp] = ndgrid(xnodes,ynodes);
                u_f_NotNanInd = find(~isnan(obj.u_f(:)));
                
                op1 = rbfcreate( [obj.x0temp_f(u_f_NotNanInd),obj.y0temp_f(u_f_NotNanInd)]',(obj.u_f(u_f_NotNanInd))','RBFFunction', 'thinplate'); %rbfcheck(op1);
                obj.u = rbfinterp( [x0temp(:),y0temp(:)]', op1 );
                op2 = rbfcreate( [obj.x0temp_f(u_f_NotNanInd),obj.y0temp_f(u_f_NotNanInd)]',(obj.v_f(u_f_NotNanInd))','RBFFunction', 'thinplate'); %rbfcheck(op2);
                obj.v = rbfinterp([x0temp(:),y0temp(:)]', op2 );
                x0temp = x0temp'; y0temp = y0temp'; obj.u=obj.u'; obj.v=obj.v';
                
                %%%%% Do some regularization to further decrease the noise %%%%%
                % u = regularizeNd([x0temp(:),y0temp(:)],u(:),{xnodes',ynodes'},1e-3);
                % v = regularizeNd([x0temp(:),y0temp(:)],v(:),{xnodes',ynodes'},1e-3);
                
                [DICMESH] = MeshSetUp(x0temp,y0temp,DICpara);
                obj = unzipDICmesh(obj,DICMESH);
                
                % ====== Initial Value ======
                obj.U0 = Init(obj.u,obj.v,obj.ccmax,obj.DICmesh_x0,obj.DICmesh_y0,0); % PlotuvInit; [x0temp,y0temp,u,v,cc]= IntegerSearchMg(ImgRef,ImgDef,file_name,DICpara);
                % ====== Deal with incremental mode ======
                fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,obj.DICparaImgSeqIncUnit)-1;
                
                if obj.DICparaImgSeqIncUnit == 1, fNormalizedNewIndex = fNormalizedNewIndex-1; end
                
                obj.ResultFEMesh_coordinatesFEM{1+floor(fNormalizedNewIndex/obj.DICparaImgSeqIncUnit)} = obj.DICmesh_coordinatesFEM;
                obj.ResultFEMesh_elementsFEM{1+floor(fNormalizedNewIndex/obj.DICparaImgSeqIncUnit)} = obj.DICmesh_elementsFEM;
                obj.ResultFEMesh_winsize{1+floor(fNormalizedNewIndex/obj.DICparaImgSeqIncUnit)} = obj.DICparawinsize;
                obj.ResultFEMesh_winstepsize{1+floor(fNormalizedNewIndex/obj.DICparaImgSeqIncUnit)} = obj.DICparawinstepsize;
                obj.ResultFEMesh_gridx{1+floor(fNormalizedNewIndex/obj.DICparaImgSeqIncUnit)} = obj.DICparagridx;
                obj.ResultFEMesh_gridy{1+floor(fNormalizedNewIndex/obj.DICparaImgSeqIncUnit)} = obj.DICparagridy;
                
                if obj.quadtree
                    % ====== Generate a quadtree mesh considering sample's complex geometry   83
                    GenerateQuadtreeMesh; % Generate a quadtree mesh
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else % Use the solved results from the last frame as the new initial guess
                if ImgSeqNum < 7 % Import previous U for ImgSeqNum [2,6]
                    obj.U0 = obj.ResultDisp{ImgSeqNum-2};
                    
                else % When ImgSeqNum > 6: POD predicts next disp U0 from previous results of (ImgSeqNum+[-5:1:-1])
                    nTime = 5; np = length(obj.ResultDisp{ImgSeqNum-2})/2; % "nTime" value 5 is an empirical value, can be changed.
                    T_data_u = zeros(nTime,np); T_data_v = zeros(nTime,np);
                    for tempi = 1:nTime
                        T_data_u(tempi,:) = obj.ResultDisp{ImgSeqNum-(2+nTime)+tempi, 1}(1:2:np*2)';
                        T_data_v(tempi,:) = obj.ResultDisp{ImgSeqNum-(2+nTime)+tempi, 1}(2:2:np*2)';
                        
                    end
                    nB = 3; t_train = (ImgSeqNum-1-nTime:ImgSeqNum-2)';
                    t_pre = (ImgSeqNum-1)';
                    [u_pred,~,~,~] = funPOR_GPR(T_data_u,t_train,t_pre,nB);
                    [v_pred,~,~,~] = funPOR_GPR(T_data_v,t_train,t_pre,nB);
                    tempu = u_pred(1,:); tempv = v_pred(1,:);
                    obj.U0 = [tempu(:),tempv(:)]';
                    obj.U0 = obj.U0(:);
                    
                    % %%%%% After running the new ImgSeqNum, you can uncomment these
                    % %%%%% lines to compare how the initial guess has been improved.
                    % Plotdisp_show(U0-ResultDisp{ImgSeqNum-1}.U,DICmesh_coordinatesFEMWorld,DICmesh_elementsFEM(:,1:4));
                    % Plotdisp_show(ResultDisp{ImgSeqNum-2}.U-ResultDisp{ImgSeqNum-1}.U,DICmesh_coordinatesFEMWorld,DICmesh_elementsFEM(:,1:4));
                end
            end
            
            % ====== Compute f(X)-g(x+u) ======
            obj.ResultFEMeshEachFrame_coordinatesFEM{ImgSeqNum-1} = obj.DICmesh_coordinatesFEM;
            obj.ResultFEMeshEachFrame_elementsFEM{ImgSeqNum-1} = obj.DICmesh_elementsFEM;
            
        end
        
        function [obj,DICpara,x0,y0,uIS,vIS,cc] = InitSearchSize(obj,ImgRef,ImgDef)
            
            %function [DICpara,x0,y0,u,v,cc]= IntegerSearch(ImgRef,ImgDef,file_name,DICpara)
            %FUNCTION [DICpara,x0,y0,u,v,cc]= IntegerSearch(ImgRef,ImgDef,file_name,DICpara)
            % Objective: To compute an inititial guess of the unknown displacement
            % field by maximizing the FFT-based cross correlation
            % ----------------------------------------------
            %   INPUT: ImgRef       Reference image
            %          ImgDef       Deformed image
            %          file_name    Loaded DIC raw images file name
            %          DICpara      Current DIC parameters
            %
            %   OUTPUT: DICpara     Updated DIC parameters
            %           x0,y0       DIC subset x- and y- positions
            %           u,v         x- and y- displacements
            %           cc          Cross correlation information
            %
            % ----------------------------------------------
            % Author: Jin Yang.
            % Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
            % Last time updated: 02/2020.
            % ==============================================
            
            %% Initialization
            % zip all the DIC options into a struct
            DICpara = zipDICPara(obj);
            % zip all the DIC mesh parameters into a struct
            DICmesh = zipDICmesh(obj);
            
            gridxROIRange = DICpara.gridxyROIRange.gridx;
            gridyROIRange = DICpara.gridxyROIRange.gridy;
            winsize = DICpara.winsize;
            winstepsize = DICpara.winstepsize;
            
            %% To compute the inititial guess from maximizing the FFT-based cross correlation
            if obj.DICparaInitFFTSearchMethod > 0
                InitialGuessSatisfied = 1;
                while InitialGuessSatisfied == 1
                    
                    if length(obj.DICparaSizeOfFFTSearchRegion) == 1, obj.DICparaSizeOfFFTSearchRegion = obj.DICparaSizeOfFFTSearchRegion*[1,1]; end
                    
                    if (obj.DICparaInitFFTSearchMethod == 1) % whole field for initial guess,
                        [x0,y0,uIS,vIS,cc] = funIntegerSearch(ImgRef,ImgDef,obj.DICparaSizeOfFFTSearchRegion,gridxROIRange,gridyROIRange,winsize,winstepsize,0,winstepsize);
                        
                    else % (InitFFTSearchMethod == 1), several local seeds for initial guess
                        
                        %pts = [row,col]; % this is the ginput selection for the local seed points
                        
                        [x0,y0,uIS,vIS,cc] = funIntegerSearch(ImgRef,ImgDef,obj.DICparaSizeOfFFTSearchRegion,gridxROIRange,gridyROIRange,winsize,winstepsize,1,obj.pts);
                        
                    end
                    % store the cc struct into normal variables under obj
                    obj = unzipcc(obj,cc);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Apply ImgRefMask to make u,v nans if there is a hole
                    try
                        x0y0Ind = sub2ind(DICpara.ImgSize, x0(:), y0(:));
                        temp1 = double(DICpara.ImgRefMask(x0y0Ind));
                        temp1(~logical(temp1))=nan;
                        HolePtIndMat=reshape(temp1,size(x0));
                        uIS = uIS.*HolePtIndMat;
                        vIS = vIS.*HolePtIndMat;
                    catch
                        
                    end
                    
                end
                
                % ======== Find some bad inital guess points ========
                obj.ccThreshold = 1.25; % bad cross-correlation threshold (mean - ccThreshold*stdev for q-factor distribution)
                qDICOrNot = 0.5;
                Thr0 = 100;
                cc = zipcc(obj);
                [uIS,vIS,~] = funRemoveOutliers(uIS,vIS,cc,qDICOrNot,Thr0);
                
                %%
            else % Multigrid search
                
                [x0,y0,uIS,vIS,cc] = funIntegerSearchMg(ImgRef,ImgDef,gridxROIRange,gridyROIRange,winsize,winstepsize,winstepsize);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Apply ImgRefMask to make u,v nans if there is a hole
                try
                    x0y0Ind = sub2ind(DICpara.ImgSize, x0(:), y0(:));
                    temp1 = double(DICpara.ImgRefMask(x0y0Ind));
                    temp1(~logical(temp1))=nan;
                    HolePtIndMat=reshape(temp1,size(x0));
                    uIS = uIS.*HolePtIndMat;
                    vIS = vIS.*HolePtIndMat;
                catch
                end
                
            end
            
        end
        
        % the main DIC subproblem steps
        function obj = Subproblem1(obj,ImgRef,ImgDef)
            % Section 4: Subproblem 1 -or- Local ICGN Subset DIC
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This section is to solve the first local step in ALDIC: Subproblem 1
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ====== ALStep 1 Subproblem1: Local Subset DIC ======
            
            % zip all the DIC options into a struct
            DICpara = zipDICPara(obj);
            % zip all the DIC mesh parameters into a struct
            DICmesh = zipDICmesh(obj);
            
            [obj.USubpb1{1},obj.FSubpb1{1},obj.HtempPar,obj.ALSub1Timetemp,obj.ConvItPerEletemp,obj.LocalICGNBadPtNumtemp] = ...
                LocalICGN(obj.U0,obj.DICmesh_coordinatesFEM,obj.Df,ImgRef,ImgDef,DICpara,'GaussNewton',obj.DICparatol);
            
        end
        
        function obj = Subproblem2(obj)
            % Section 5: Subproblem 2 -- solve the global compatible displacement field
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This section is to solve the global step in ALDIC Subproblem 2
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % init all parameters
            obj.DICparaDispFilterSize=0;
            obj.DICparaDispFilterStd=0;
            obj.DICparaStrainFilterSize=0;
            obj.DICparaStrainFilterStd=0;
            
            % zip all the DIC options into a struct
            DICpara = zipDICPara(obj);
            % zip all the DIC mesh parameters into a struct
            DICmesh = zipDICmesh(obj);
            
            
            % ======= ALStep 1 Subproblem 2: Global constraint =======
            % ------ Smooth displacements for a better F ------
            obj.FSubpb1{1} = funSmoothStrain(obj.FSubpb1{1},DICmesh,DICpara);
            
            % ====== Define penalty parameter ======
            obj.udual{1} = 0*obj.FSubpb1{1};
            obj.vdual{1} = 0*obj.USubpb1{1};
            betaList = [1e-3,sqrt(1e-5),1e-2,sqrt(1e-3),1e-1,sqrt(1e-1)]*mean(obj.DICparawinstepsize).^2.*obj.DICparamu; % Tune beta in the betaList
            Err1 = zeros(length(betaList),1);
            Err2 = Err1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ====== Check to use FD or FE methods to solve Subpb2 step ======
            if obj.DICparaSubpb2FDOrFEM == 1 % Use FD method
                % ====== Build sparse finite difference operator ======
                Rad = 1;
                M = size(obj.DICmesh_x0,1);
                N = size(obj.DICmesh_x0,2);
                D = funDerivativeOp((M-2*Rad),(N-2*Rad),obj.DICparawinstepsize); % D = sparse(4*(M-2*Rad)*(N-2*Rad), 2*(M-2*Rad)*(N-2*Rad));
                D2 = funDerivativeOp(M,N,obj.DICparawinstepsize);
                
                % ===== Solver using finite difference approximation ======
                a = obj.FSubpb1{1}-obj.udual{1};
                b = obj.USubpb1{1}-obj.vdual{1};
                [obj.temp3,obj.temp4] = funFDNeumannBCInd(size(obj.DICmesh_coordinatesFEM,1),M,N,Rad); % Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
                atemp = a(obj.temp3);
                btemp = b(obj.temp4);
                
                for tempk = 1:length(betaList)
                    obj.DICparabeta = betaList(tempk);
                    tempAMatrixSub2 = (obj.DICparabeta*(D')*D) + obj.DICparamu*speye(2*(M-2*Rad)*(N-2*Rad));
                    USubpb2temp = (tempAMatrixSub2) \ (obj.DICparabeta*D'*atemp + obj.DICparamu*btemp ) ;
                    obj.USubpb2{1} = obj.USubpb1{1};
                    obj.USubpb2{1}(obj.temp4) = USubpb2temp;
                    obj.FSubpb2{1} = D2*obj.USubpb2{1};
                    
                    Err1(tempk) = norm(obj.USubpb1{1}-obj.USubpb2{1},2);
                    Err2(tempk) = norm(obj.FSubpb1{1}-obj.FSubpb2{1},2);
                    
                end
                ErrSum = Err1+Err2*mean(obj.DICparawinstepsize)^2;
                [~,indexOfbeta] = min(ErrSum);
                
                try % Tune the best beta by a quadratic polynomial fitting
                    [fitobj] = fit(log10(betaList(indexOfbeta-1:1:indexOfbeta+1))',ErrSum(indexOfbeta-1:1:indexOfbeta+1),'poly2');
                    p = coeffvalues(fitobj);
                    obj.DICparabeta = 10^(-p(2)/2/p(1));
                catch, obj.DICparabeta = betaList(indexOfbeta);
                end
                
                % Using the optimal beta to solve the ALDIC Subproblem 2 again
                tempAMatrixSub2 = (obj.DICparabeta*(D')*D) + obj.DICparamu*speye(2*(M-2*Rad)*(N-2*Rad));
                USubpb2temp = (tempAMatrixSub2) \ (obj.DICparabeta*D'*atemp + obj.DICparamu*btemp) ;
                obj.USubpb2{1} = obj.USubpb1{1};
                obj.USubpb2{1}(obj.temp4) = USubpb2temp;
                
                %%%%%%%%%%%%%% End of using finite difference approximation %%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else %Subpb2FDOrFEM: Using FE method
                M = size(obj.DICmesh_x0,1);
                N = size(obj.DICmesh_x0,2);
                obj.GaussPtOrder = 2;
                obj.alpha = 0;
                % ====== Solver using finite element method ======
                for tempk = 1:length(betaList)
                    obj.DICparabeta = betaList(tempk);
                    obj.GaussPtOrder = 2;
                    obj.alpha = 0;
                    
                    [obj.USubpb2{1}] = Subpb2(DICmesh,obj.DICparabeta,obj.DICparamu,obj.USubpb1{1},obj.FSubpb1{1},obj.udual{1},obj.vdual{1},obj.alpha,obj.GaussPtOrder);
                    [obj.FSubpb2{1}] = funGlobal_NodalStrainAvg(obj.DICmesh_coordinatesFEM,obj.DICmesh_elementsFEM,obj.USubpb2{1},obj.GaussPtOrder);
                    
                    Err1(tempk) = norm(obj.USubpb1{1}-obj.USubpb2{1},2);
                    Err2(tempk) = norm(obj.FSubpb1{1}-obj.FSubpb2{1},2);
                    
                end
                Err1Norm = (Err1-mean(Err1))/std(Err1);
                Err2Norm = (Err2-mean(Err2))/std(Err2);
                ErrSum = Err1Norm+Err2Norm;
                [~,indexOfbeta] = min(ErrSum);
                
                try % Tune the best beta by a quadratic polynomial fitting
                    [fitobj] = fit(log10(betaList(indexOfbeta-1:1:indexOfbeta+1))',ErrSum(indexOfbeta-1:1:indexOfbeta+1),'poly2');
                    p = coeffvalues(fitobj);
                    obj.DICparabeta = 10^(-p(2)/2/p(1));
                catch, obj.DICparabeta = betaList(indexOfbeta);
                end
                
                % Using the optimal beta to solve the ALDIC Subproblem 2 again
                [obj.USubpb2{1}] = Subpb2(DICmesh,obj.DICparabeta,obj.DICparamu,obj.USubpb1{1},obj.FSubpb1{1},obj.udual{1},obj.vdual{1},obj.alpha,obj.GaussPtOrder);
                obj.USubpb2{1} = full(obj.USubpb2{1});
                
            end
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ------- Before computing strain, we smooth the displacement field -------
            %obj.USubpb2{1} = funSmoothDisp(obj.USubpb2{1},coordinatesFEM,elementsFEM,x0,y0,winstepsize,DispFilterSize,DispFilterStd);
            
            % ------- Compute strain field --------
            if obj.DICparaSubpb2FDOrFEM == 1 %FD
                obj.FSubpb2{1} = D2*obj.USubpb2{1}; % D2 = funDerivativeOp(M,N,winstepsize);
            else %FEM
                [obj.FSubpb2{1}] = funGlobal_NodalStrainAvg(obj.DICmesh_coordinatesFEM,obj.DICmesh_elementsFEM,obj.USubpb2{1},obj.GaussPtOrder);
            end
            
            % ------- Smooth strain field --------
            obj.FSubpb2{1} = funSmoothStrain(obj.FSubpb2{1},DICmesh,DICpara);
            
            % ======= Update dual variables =======
            if obj.DICparaSubpb2FDOrFEM == 1 %FD
                udualtemp1 = (obj.FSubpb2{1} - obj.FSubpb1{1});
                udualtemp2 = udualtemp1(obj.temp3);
                vdualtemp1 = (obj.USubpb2{1} - obj.USubpb1{1});
                vdualtemp2 = vdualtemp1(obj.temp4);
                obj.udual{1} = zeros(4*M*N,1);
                obj.vdual{1} = zeros(2*M*N,1);
                obj.udual{1}(obj.temp3) = udualtemp2;
                obj.vdual{1}(obj.temp4) = vdualtemp2;
                
            else  % FEM or other methods
                obj.udual{1} = obj.FSubpb2{1} - obj.FSubpb1{1};
                obj.vdual{1} = obj.USubpb2{1} - obj.USubpb1{1};
                
            end
            
            
        end
        
        function obj = ADMMiteration(obj,ImgRef,ImgDef,ImgSeqNum)
            
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Section 6: ADMM iterations
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This section is to run ADMM iteration: Subproblem 1 & 2
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ==================== ADMM AL Loop ==========================
            
            % zip all the DIC options into a struct
            DICpara = zipDICPara(obj);
            % zip all the DIC mesh parameters into a struct
            DICmesh = zipDICmesh(obj);
            
            obj.ALSolveStep = 1;
            tol2 = 1e-4;
            UpdateY = 1e4;
            %             CrackOrNot = 0;
            %             CrackPath1 = [0,0];
            %             CrackPath2 = [0,0];
            %             CrackTip = [0,0];
            HPar = cell(21,1);
            for tempj = 1:21, HPar{tempj} = obj.HtempPar(:,tempj); end
            
            % preform ADMM
            while (obj.ALSolveStep < 5)
                % update variables
                TempUSubpb2 = obj.USubpb2{obj.ALSolveStep};
                TempFSubpb2 = obj.FSubpb2{obj.ALSolveStep};
                Tempudual = obj.udual{obj.ALSolveStep};
                Tempvdual = obj.vdual{obj.ALSolveStep};
                
                obj.ALSolveStep = obj.ALSolveStep + 1;  % Update using the last step
                
                %%%%%%%%%%%%%%%%%%%%%%% Subproblem 1 %%%%%%%%%%%%%%%%%%%%%%%%%
                
                [TempUSubpb1,~,~,~] = Subpb1( ...
                    TempUSubpb2, TempFSubpb2, Tempudual,Tempvdual, obj.DICmesh_coordinatesFEM,...
                    obj.Df,ImgRef,ImgDef,obj.DICparamu,obj.DICparabeta,...
                    HPar,obj.ALSolveStep, DICpara,'GaussNewton',obj.DICparatol);
                
                TempFSubpb1 = TempFSubpb2;
                
                %Save the results of U and F sub prob 1 to array
                obj.USubpb1{obj.ALSolveStep} = TempUSubpb1;
                obj.FSubpb1{obj.ALSolveStep} = TempFSubpb1;
                
                TempUSubpb1 = funSmoothDisp(TempUSubpb1,DICmesh,DICpara);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % ============== Subproblem 2 ==============
                
                if obj.DICparaSubpb2FDOrFEM == 1 %FD
                    % use finite difference approximation
                    a = TempFSubpb1-Tempudual;
                    b = TempUSubpb1-Tempvdual;
                    atemp = a(obj.temp3);
                    btemp = b(obj.temp4);
                    USubpb2temp = (tempAMatrixSub2) \ (obj.DICparabeta*D'*atemp + obj.DICparamu*btemp ) ;
                    TempUSubpb2 = TempUSubpb1;
                    TempUSubpb2(obj.temp4) = USubpb2temp;
                    % End of using finite difference approximation
                else % FEM
                    % use finite element method
                    [TempUSubpb2] = Subpb2(DICmesh,obj.DICparabeta,obj.DICparamu,TempUSubpb1,TempFSubpb1,Tempudual,Tempvdual,obj.alpha,obj.GaussPtOrder);
                    TempUSubpb2 = full(TempUSubpb2);
                end
                
                % Before computing strain, we smooth the displacement field
                TempUSubpb2 = funSmoothDisp(TempUSubpb2,DICmesh,DICpara);
                
                % Compute strain field
                if obj.DICparaSubpb2FDOrFEM == 1 %FD
                    TempFSubpb2 = D2*TempUSubpb2; % D2 = funDerivativeOp(M,N,winstepsize);
                else %FEM
                    obj.GaussPtOrder = 2;
                    [TempFSubpb2] = funGlobal_NodalStrainAvg(obj.DICmesh_coordinatesFEM,obj.DICmesh_elementsFEM,TempUSubpb2,obj.GaussPtOrder);
                end
                % ------- Smooth strain field --------
                TempFSubpb2 = funSmoothStrain(TempFSubpb2,DICmesh,DICpara);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %Save the results of U and F sub prob 1 to array
                obj.USubpb2{obj.ALSolveStep} = TempUSubpb2;
                obj.FSubpb2{obj.ALSolveStep} = TempFSubpb2;
                
                
                % Compute norm of UpdateY
                USubpb2_Old = obj.USubpb2{obj.ALSolveStep-1};
                USubpb2_New = obj.USubpb2{obj.ALSolveStep};
                USubpb1_Old = obj.USubpb1{obj.ALSolveStep-1};
                USubpb1_New = obj.USubpb1{obj.ALSolveStep};
                
                if (mod(ImgSeqNum-2,obj.DICparaImgSeqIncUnit) ~= 0 && (ImgSeqNum>2)) || (ImgSeqNum < obj.DICparaImgSeqIncUnit)
                    UpdateY = norm((USubpb2_Old.USubpb2 - USubpb2_New.USubpb2), 2)/sqrt(size(USubpb2_Old.USubpb2,1));
                    try
                        UpdateY2 = norm((USubpb1_Old.USubpb1 - USubpb1_New.USubpb1), 2)/sqrt(size(USubpb1_Old.USubpb1,1));
                    catch
                    end
                end
                
                % Update dual variables
                if obj.DICparaSubpb2FDOrFEM == 1 %FD
                    udualtemp1 =  (TempFSubpb2 - TempFSubpb1);
                    udualtemp2 = udualtemp1(obj.temp3);
                    vdualtemp1 =  (TempUSubpb2 - TempUSubpb1);
                    vdualtemp2 = vdualtemp1(obj.temp4);
                    Tempudual(obj.temp3) = Tempudual(obj.temp3)+udualtemp2;
                    Tempvdual(obj.temp4) = Tempvdual(obj.temp4)+vdualtemp2;
                else %FEM
                    Tempudual = TempFSubpb2 - TempFSubpb1;
                    Tempvdual = TempUSubpb2 - TempUSubpb1;
                end
                
                obj.udual{obj.ALSolveStep} = Tempudual;
                obj.vdual{obj.ALSolveStep} = Tempvdual;
                
                try
                    if UpdateY < tol2 || UpdateY2 < tol2
                        break
                    end
                catch
                end
                
            end
            
            % Save data
            obj.ResultDisp{ImgSeqNum-1} = full(TempUSubpb2);
            obj.ResultDefGrad{ImgSeqNum-1} = full(TempFSubpb2); % tempFoamAL;
            % return the u and v displacement fields
            U = DICpara.um2px*obj.ResultDisp{ImgSeqNum-1};
            obj.Ux = U(1:2:end); obj.Vy = U(2:2:end);
            u0 = reshape(obj.Ux,obj.DICmesh_M,obj.DICmesh_N); v0 = reshape(obj.Vy,obj.DICmesh_M,obj.DICmesh_N);
            obj.Ux = u0; obj.Vy = v0;
        end
        
        function obj = checkConvergence(obj)
            %% Section 7: Check convergence
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This section is to check convergence of ADMM
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ====== Check convergence ======
            ALSolveStep1 = min(6,obj.ALSolveStep);
            
            %'==== uhat^(k) - u^(k) ====');
            for ALSs = 1:ALSolveStep1
                TempUSubpb2 = obj.USubpb2{ALSs};
                TempUSubpb1 = obj.USubpb1{ALSs};
                obj.Convergence_Y = norm((TempUSubpb2 - TempUSubpb1), 2)/sqrt(length(TempUSubpb2));
            end
            
            %'==== Fhat^(k) - F^(k) ====');
            for ALSs = 1:ALSolveStep1
                TempFSubpb1 = obj.FSubpb1{ALSs};
                TempFSubpb2 = obj.FSubpb2{ALSs};
                obj.Convergence_F = norm((TempFSubpb1 - TempFSubpb2), 2)/sqrt(length(TempFSubpb1));
                
            end
            
            %'==== uhat^(k) - uhat^(k-1) ====');
            for ALSs = 2:ALSolveStep1
                USubpb2_Old = obj.USubpb2{ALSs-1};
                USubpb2_New = obj.USubpb2{ALSs};
                obj.Convergence_U = norm((USubpb2_Old - USubpb2_New), 2)/sqrt(length(TempUSubpb2));
            end
            
            %'==== udual^(k) - udual^(k-1) ====');
            for ALSs = 2:ALSolveStep1
                uvdual_Old = obj.udual{ALSs-1};
                uvdual_New = obj.udual{ALSs};
                obj.Convergence_W = norm((uvdual_Old - uvdual_New), 2)/sqrt(length(uvdual_Old));
                
            end
            %'==== vdual^(k) - vdual^(k-1) ====');
            for ALSs = 2:ALSolveStep1
                uvdual_Old = obj.vdual{ALSs-1};
                uvdual_New = obj.vdual{ALSs};
                obj.Convergence_V = norm((uvdual_Old - uvdual_New), 2)/sqrt(length(uvdual_Old));
                
            end
            
            
        end
        
        function obj = computeStrains(obj,StrainOptions,ImgSeqNum)
            %% Compute strains
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This section is to compute strain fields and plot disp and strain results
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % zip all the DIC options into a struct
            DICpara = zipDICPara(obj);

            %%%%%%% GUI version of the code %%%%%%%
            % this will allow the user to input options (StrainOptions) to
            % return a strain map for a single image
            
            ULocal = obj.ResultDisp{ImgSeqNum}; % get the displacement result information
            FLocal = obj.ResultDefGrad{ImgSeqNum};
            coordinatesFEM = obj.ResultFEMeshEachFrame_coordinatesFEM{ImgSeqNum};% get the mesh used
            elementsFEM = obj.ResultFEMeshEachFrame_elementsFEM{ImgSeqNum};% get the mesh elements
            winstepsize = DICpara.winstepsize;
            xList = min(coordinatesFEM(:,1)):winstepsize:max(coordinatesFEM(:,1));
            M = length(xList);
            yList = min(coordinatesFEM(:,2)):winstepsize:max(coordinatesFEM(:,2));
            N = length(yList);
                 
            switch StrainOptions.Method
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 0 % ALDIC directly solved deformation gradients
                    FStrain = FLocal;
                    
                    Rad = 0;
                    
                    temp = 1:1:size(coordinatesFEM,1);
                    temp = temp';
                    temp = reshape(temp,M,N);
                    temp2 = temp(Rad+1:M-Rad, Rad+1:N-Rad);
                    temp2 = reshape(temp2, (M-2*Rad)*(N-2*Rad),1);
                    
                    T3 = zeros(4*(M-2*Rad)*(N-2*Rad),1);
                    for i = 1:(M-2*Rad)*(N-2*Rad)
                        T3(4*i-3:4*i) = [4*temp2(i)-3; 4*temp2(i)-2; 4*temp2(i)-1; 4*temp2(i)];
                    end
                    obj.ResultStrainDirect{ImgSeqNum} = FStrain;
                    obj.ResultStrainRad{ImgSeqNum,StrainOptions.Method+1} = Rad;
                    obj.ResultStrainT3{ImgSeqNum,StrainOptions.Method+1} = T3;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 1 % Central finite difference
                    
                    % Compute strain method I: Use Finite difference operator
                    D = funDerivativeOp(M,N,winstepsize); % D = sparse(4*M*N, 2*M*N);
                    FStrain = D*reshape(ULocal,length(ULocal),1);
                    
                    Rad = 1;
                    
                    % Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
                    temp = 1:1:size(coordinatesFEM,1); temp = temp';
                    temp = reshape(temp,M,N); temp2 = temp(Rad+1:M-Rad, Rad+1:N-Rad);
                    temp2 = reshape(temp2, (M-2*Rad)*(N-2*Rad),1);
                    
                    T3 = zeros(4*(M-2*Rad)*(N-2*Rad),1);
                    for i = 1:(M-2*Rad)*(N-2*Rad)
                        T3(4*i-3:4*i) = [4*temp2(i)-3; 4*temp2(i)-2; 4*temp2(i)-1; 4*temp2(i)];
                    end
                    
                    obj.ResultStrainFD{ImgSeqNum} = FStrain;
                    obj.ResultStrainRad{ImgSeqNum,StrainOptions.Method+1} = Rad;
                    obj.ResultStrainT3{ImgSeqNum,StrainOptions.Method+1} = T3;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 2 % Plane fitting
                    D = funDerivativeOp(M,N,winstepsize); % D = sparse(4*M*N, 2*M*N);
                    FStrain = D*reshape(ULocal,length(ULocal),1);
                    
                    % Compute strain method II: Use Plane Fitting method
                    Rad = StrainOptions.PWS;
                    
                    [Uytemp, Uxtemp, ~] = PlaneFit2(reshape(ULocal(1:2:end),M,N), winstepsize, winstepsize,Rad);
                    [Vytemp, Vxtemp, ~] = PlaneFit2(reshape(ULocal(2:2:end),M,N), winstepsize, winstepsize,Rad);
                    
                    FStraintemp = zeros(4*(M-2*Rad)*(N-2*Rad),1);
                    FStraintemp(1:4:end) = reshape(Uxtemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
                    FStraintemp(2:4:end) = reshape(Vxtemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
                    FStraintemp(3:4:end) = reshape(Uytemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
                    FStraintemp(4:4:end) = reshape(Vytemp((Rad+1):M-Rad,(Rad+1):N-Rad), (M-2*Rad)*(N-2*Rad),1);
                    
                    % Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
                    temp = 1:1:size(coordinatesFEM,1); temp = temp';
                    temp = reshape(temp,M,N); temp2 = temp(Rad+1:M-Rad, Rad+1:N-Rad);
                    temp2 = reshape(temp2, (M-2*Rad)*(N-2*Rad),1);
                    
                    T3 = zeros(4*(M-2*Rad)*(N-2*Rad),1);
                    for i = 1:(M-2*Rad)*(N-2*Rad)
                        T3(4*i-3:4*i) = [4*temp2(i)-3; 4*temp2(i)-2; 4*temp2(i)-1; 4*temp2(i)];
                    end
                    
                    FStrain(T3) = FStraintemp;
                    obj.ResultStrainPlane{ImgSeqNum} = FStrain;
                    obj.ResultStrainRad{ImgSeqNum,StrainOptions.Method+1} = Rad;
                    obj.ResultStrainT3{ImgSeqNum,StrainOptions.Method+1} = T3;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case 3 % Finite element method
                    
                    GPO = 2;
                    [FStrain] = funGlobal_NodalStrainAvg(coordinatesFEM,elementsFEM,ULocal,GPO);
                    Rad = 1;
                    
                    % Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
                    temp = 1:1:size(coordinatesFEM,1);
                    temp = temp';
                    temp = reshape(temp,M,N);
                    temp2 = temp(Rad+1:M-Rad, Rad+1:N-Rad);
                    temp2 = reshape(temp2, (M-2*Rad)*(N-2*Rad),1);
                    
                    T3 = zeros(4*(M-2*Rad)*(N-2*Rad),1);
                    for i = 1:(M-2*Rad)*(N-2*Rad)
                        T3(4*i-3:4*i) = [4*temp2(i)-3; 4*temp2(i)-2; 4*temp2(i)-1; 4*temp2(i)];
                    end
                    obj.ResultStrainFE{ImgSeqNum} = FStrain;
                    obj.ResultStrainRad{ImgSeqNum,StrainOptions.Method+1} = Rad;
                    obj.ResultStrainT3{ImgSeqNum,StrainOptions.Method+1} = T3;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                otherwise
                    disp('Wrong Input to compute strain field!')
                    
            end
            
        end
        
        function obj = convertStrainType(obj,StrainOptions,ImgSeqNum)
            %% Update infinitesimal strain to other finite strains
            
            % check to see if the strain has been calculated
            % if not first calculate those strains
            switch StrainOptions.Method
                case 0
                    FStrain = obj.ResultStrainDirect{ImgSeqNum};
                case 1
                    FStrain = obj.ResultStrainFD{ImgSeqNum};
                case 2
                    FStrain = obj.ResultStrainPlane{ImgSeqNum};
                case 3
                    FStrain = obj.ResultStrainFE{ImgSeqNum};
            end
            
            FStrainFinite = FStrain;
            for tempi = 1:4:length(FStrain)
                
                % Obtain each component of def grad tensor
                dudx = FStrain(tempi);
                dvdx = FStrain(tempi+1);
                dudy = FStrain(tempi+2);
                dvdy = FStrain(tempi+3);
                
                switch StrainOptions.StrainType
                    case 0 % Infinitesimal stran
                        % Do nothing
                    case 1 % Eluerian strain
                        FStrainFinite(tempi) = 1/(1-dudx)-1;
                        FStrainFinite(tempi+3) = 1/(1-dvdy)-1;
                        FStrainFinite(tempi+2) = dudy/(1-dvdy);
                        FStrainFinite(tempi+1) = dvdx/(1-dudx);
                    case 2 % Green-Lagrangian strain: E=(C-I)/2
                        FStrainFinite(tempi) = 0.5*(dudx*2-dudx^2-dvdx^2);
                        FStrainFinite(tempi+3) = 0.5*(dvdy*2-dudy^2-dvdy^2);
                        FStrainFinite(tempi+2) = 0.5*(dudy+dvdx-dudx*dudy-dvdx*dvdy);
                        FStrainFinite(tempi+1) = 0.5*(dvdx+dudy-dudy*dudx-dvdy*dvdx);
                    otherwise
                        disp('Wrong strain type!');
                end
            end
            
            T3 =  obj.ResultStrainT3{ImgSeqNum,StrainOptions.Method+1};
            FStraintemp = FStrainFinite(T3);
            FStrainWorld = FStraintemp;
            FStrainWorld(2:4:end) = -FStrainWorld(2:4:end);
            FStrainWorld(3:4:end) = -FStrainWorld(3:4:end);
            
            obj.ResultStrainWorld{ImgSeqNum,StrainOptions.Method+1,StrainOptions.StrainType+1} = FStrainWorld;
        end
        
        function [StrainMap, obj] = getStrainMap(obj,StrainOptions,ImgSeqNum)
            
            % zip all the DIC options into a struct
            DICpara = zipDICPara(obj);
            
            % temp 
            DICpara.um2px = 1;
            
            ULocal = obj.ResultDisp{ImgSeqNum}; % get the displacement result information
            UFEMesh = 0*ULocal;
            U = DICpara.um2px*ULocal;
            F = obj.ResultStrainWorld{ImgSeqNum,StrainOptions.Method+1,StrainOptions.StrainType+1};
            Rad = obj.ResultStrainRad{ImgSeqNum,StrainOptions.Method+1};
            
            coordinatesFEM = obj.ResultFEMeshEachFrame_coordinatesFEM{ImgSeqNum};% get the mesh used
            xList = min(coordinatesFEM(:,1)):DICpara.winstepsize:max(coordinatesFEM(:,1));
            M = length(xList);
            yList = min(coordinatesFEM(:,2)):DICpara.winstepsize:max(coordinatesFEM(:,2));
            N = length(yList);
            [x0,y0] = ndgrid(xList,yList);
            x0 = x0-reshape(UFEMesh(1:2:end),size(x0,1),size(x0,2));
            y0 = y0-reshape(UFEMesh(2:2:end),size(y0,1),size(y0,2));
            x0 = DICpara.um2px*x0;
            y0 = DICpara.um2px*y0; % Ignore this: (size(ImgNormalized{1},2)+1-y0);
            sizeOfImg = DICpara.ImgSize;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Compute strain components
            
            x = x0(1+Rad:end-Rad,1+Rad:end-Rad);
            y = y0(1+Rad:end-Rad,1+Rad:end-Rad);
            
            M = size(x,1); N = size(x,2);
            u_x = F(1:4:end); v_x = F(2:4:end);
            u_y = F(3:4:end); v_y = F(4:4:end);
            
            u_x = reshape(u_x,M,N); v_x = reshape(v_x,M,N);
            u_y = reshape(u_y,M,N); v_y = reshape(v_y,M,N);
            
            UX = U(1:2:end); VX = U(2:2:end);
            u0 = reshape(UX,M+2*Rad,N+2*Rad); v0 = reshape(VX,M+2*Rad,N+2*Rad);
            UX = u0(1+Rad:end-Rad,1+Rad:end-Rad); VX = v0(1+Rad:end-Rad,1+Rad:end-Rad);
            
            % imagesc([x(1,1) x(end,1)], [y(1,1) y(1,end)], flipud(g)); hold on;
            if M < 9, x2 = x(:,1)'; else x2 = linspace(x(1,1),x(end,1),4*(length(x(:,1))-1)+1); x2=x2(:)'; end
            if N < 9, y2 = y(1,:); else y2 = linspace(y(1,1),y(1,end),4*(length(y(1,:))-1)+1); y2=y2(:)'; end
            
            
            %% Compute displacement components to manipulate the reference image
            disp_u = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(UX,M*N,1),x2,y2);
            disp_v = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(VX,M*N,1),x2,y2);
            
            %% Compute strain components
            dudx = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(u_x,M*N,1),x2,y2);
            dvdx = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(v_x,M*N,1),x2,y2);
            dudy = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(u_y,M*N,1),x2,y2);
            dvdy = gridfit(reshape(x,M*N,1),reshape(y,M*N,1),reshape(v_y,M*N,1),x2,y2);
            
            strain_exx = dudx;
            strain_exy = 0.5*(dvdx + dudy);
            strain_eyy = dvdy;
            
            strain_maxshear = sqrt((0.5*(strain_exx-strain_eyy)).^2 + strain_exy.^2);
            % Principal strain
            strain_principal_max = 0.5*(strain_exx+strain_eyy) + strain_maxshear;
            strain_principal_min = 0.5*(strain_exx+strain_eyy) - strain_maxshear;
            % equivalent von Mises strain
            strain_vonMises = sqrt(strain_principal_max.^2 + strain_principal_min.^2 - ...
                strain_principal_max.*strain_principal_min + 3*strain_maxshear.^2);
            
            % Please don't delete this line, to deal with the image and physical world coordinates
            [x2,y2]=ndgrid(x2,y2); x2=x2'; y2=y2';
            
            xplot = x2;
            yplot = DICpara.um2px*(sizeOfImg(2)+1)-y2;
            StrainMap = {xplot,yplot,strain_exx,strain_exy,strain_eyy,strain_maxshear,strain_principal_max,strain_principal_min,strain_vonMises};
            obj.ResultStrainMap{ImgSeqNum,1} = xplot;
            obj.ResultStrainMap{ImgSeqNum,2} = yplot;
            obj.ResultStrainMap{ImgSeqNum,3} = strain_exx;
            obj.ResultStrainMap{ImgSeqNum,4} = strain_exy;
            obj.ResultStrainMap{ImgSeqNum,5} = strain_eyy;
            obj.ResultStrainMap{ImgSeqNum,6} = strain_maxshear;
            obj.ResultStrainMap{ImgSeqNum,7} = strain_principal_max;
            obj.ResultStrainMap{ImgSeqNum,8} = strain_principal_min;
            obj.ResultStrainMap{ImgSeqNum,9} = strain_vonMises;
            obj.ResultDispMapu{ImgSeqNum} = {disp_u};
            obj.ResultDispMapv{ImgSeqNum} = {disp_v};
        end
        
        
        % needs revison to incorperate GUI
        function obj = manualBadPointRemoval(obj)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ------  Manually find some bad points from Local Subset ICGN step ------
            % Comment these lines below if you don't need to manually remove local bad points
            % %%%%% Comment START %%%%%
            
            %             close all;
            %             USubpb1World = USubpb1;
            %             USubpb1World(2:2:end) = -USubpb1(2:2:end);
            %             Plotuv(USubpb1World,DICmesh.x0,DICmesh.y0World);
            %
            %             disp('--- Start to manually remove bad points ---')
            %             u = reshape(USubpb1(1:2:end),size(DICmesh.x0,1),size(DICmesh.x0,2));
            %             v = reshape(USubpb1(2:2:end),size(DICmesh.x0,1),size(DICmesh.x0,2));
            %             [u,v,~,Local_BadptRow,Local_BadptCol,RemoveOutliersList] = funRemoveOutliers(u',v',[],0.5,100);
            %             u=u';   v=v';
            %
            %             USubpb1(1:2:end) = reshape(u,size(DICmesh.coordinatesFEM,1),1);
            %             USubpb1(2:2:end) = reshape(v,size(DICmesh.coordinatesFEM,1),1);
            %
            %             f11 = reshape(FSubpb1(1:4:end),size(DICmesh.x0,1),size(DICmesh.x0,2));
            %             f21 = reshape(FSubpb1(2:4:end),size(DICmesh.x0,1),size(DICmesh.x0,2));
            %             f12 = reshape(FSubpb1(3:4:end),size(DICmesh.x0,1),size(DICmesh.x0,2));
            %             f22 = reshape(FSubpb1(4:4:end),size(DICmesh.x0,1),size(DICmesh.x0,2));
            %
            %             f11=f11';    f11(RemoveOutliersList) = NaN;    f11 = inpaint_nans(f11,4);    f11=f11';
            %             f21=f21';    f21(RemoveOutliersList) = NaN;    f21 = inpaint_nans(f21,4);    f21=f21';
            %             f12=f12';    f12(RemoveOutliersList) = NaN;    f12 = inpaint_nans(f12,4);    f12=f12';
            %             f22=f22';    f22(RemoveOutliersList) = NaN;    f22 = inpaint_nans(f22,4);    f22=f22';
            %
            %             FSubpb1(1:4:end) = reshape(f11,size(DICmesh.coordinatesFEM,1),1);
            %             FSubpb1(2:4:end) = reshape(f21,size(DICmesh.coordinatesFEM,1),1);
            %             FSubpb1(3:4:end) = reshape(f12,size(DICmesh.coordinatesFEM,1),1);
            %             FSubpb1(4:4:end) = reshape(f22,size(DICmesh.coordinatesFEM,1),1);
            %
            %             disp('--- Remove bad points done ---')
            % %%%%% Comment END %%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end
        
        
        % initialization of variables
        function obj = InitResultVars(obj,ImgCount)
            % initializes the results to account for the number of frames
            % to be analyzed
            obj.ResultDisp = cell(ImgCount-1,1);
            obj.ResultDefGrad = cell(ImgCount-1,1);
            obj.ResultStrainWorld = cell(ImgCount-1,4,3); %ImgSeqNum,StrainMethod,StrainType}
            obj.ResultStressWorld = cell(ImgCount-1,1);
            obj.ResultFEMeshEachFrame_coordinatesFEM = cell(ImgCount-1,1);
            obj.ResultFEMeshEachFrame_elementsFEM = cell(ImgCount-1,1);
            obj.ResultFEMesh_coordinatesFEM = cell(ceil((ImgCount-1)/obj.DICparaImgSeqIncUnit),1); % For incremental DIC mode
            obj.ResultFEMesh_elementsFEM = obj.ResultFEMesh_coordinatesFEM;
            obj.ResultFEMesh_winsize = obj.ResultFEMesh_coordinatesFEM;
            obj.ResultFEMesh_winstepsize = obj.ResultFEMesh_coordinatesFEM;
            obj.ResultFEMesh_gridx = obj.ResultFEMesh_coordinatesFEM;
            obj.ResultFEMesh_gridy = obj.ResultFEMesh_coordinatesFEM;
            obj.ResultStrainDirect = cell(ImgCount-1,1);
            obj.ResultStrainFD = cell(ImgCount-1,1);
            obj.ResultStrainPlane = cell(ImgCount-1,1);
            obj.ResultStrainFE = cell(ImgCount-1,1);
            obj.ResultStrainMap = cell(ImgCount-1,9);
            obj.ResultDispMapu = cell(ImgCount-1,1);
            obj.ResultDispMapv = cell(ImgCount-1,1);
        end
        
        
        % GUI interfaceing functions
        function obj = setDICoptions(obj,variable,value)
            % this function will set the DIC options into the objects
            % parameters
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            switch variable
                case 'winsize'
                    obj.DICparawinsize = value;
                case 'winstepsize'
                    obj.DICparawinstepsize = value;
                case 'Subpb2FDOrFEM' % Solve ALDIC Subproblem 2 via a '1) finite difference' or '2) finite element' solver?
                    obj.DICparaSubpb2FDOrFEM = value;
                case 'ClusterNo' % How many parallel pools to open
                    obj.DICparaClusterNo = value;
                case 'NewFFTSearch' % re-run FFT for an initial guess or not 0= last frame, 1=redo initial guess
                    obj.DICparaNewFFTSearch = value;
                case 'DICIncOrNot' % Decide DIC as 0=accumulative or 1 = incremental mode?
                    obj.DICparaDICIncOrNot = value;
                case 'ImgSeqIncUnit'
                    obj.DICparaImgSeqIncUnit = value;
                case 'ImgSeqIncROIUpdateOrNot'
                    obj.DICparaImgSeqIncROIUpdateOrNot = value;
                case 'mu'
                    obj.DICparamu = value;
                case 'beta'
                    obj.DICparabeta = value;
                case 'tol'
                    obj.DICparatol = value;
                case 'DICparagridx'
                    obj.DICparagridx = value;
                case 'DICparagridy'
                    obj.DICparagridy = value;
                case 'InitFFTSearchMethod'
                    obj.DICparaInitFFTSearchMethod = value;
                case 'SizeOfFFTSearchRegion'
                    obj.DICparaSizeOfFFTSearchRegion = value;
                    
            end
            
        end
        
        function value = readDICoptions(obj,variable)
            % this function will read the DIC options from the objects
            % parameters
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            switch variable
                case 'winsize'
                    value = obj.DICparawinsize ;
                case 'winstepsize'
                    value = obj.DICparawinstepsize;
                case 'Subpb2FDOrFEM' % Solve ALDIC Subproblem 2 via a '1) finite difference' or '2) finite element' solver?
                    value = obj.DICparaSubpb2FDOrFEM;
                case 'ClusterNo' % How many parallel pools to open
                    value = obj.DICparaClusterNo;
                case 'NewFFTSearch' % re-run FFT for an initial guess or not 0= last frame, 1=redo initial guess
                    value = obj.DICparaNewFFTSearch;
                case 'DICIncOrNot' % Decide DIC as 0=accumulative or 1 = incremental mode?
                    value = obj.DICparaDICIncOrNot;
                case 'ImgSeqIncUnit'
                    value = obj.DICparaImgSeqIncUnit;
                case 'ImgSeqIncROIUpdateOrNot'
                    value = obj.DICparaImgSeqIncROIUpdateOrNot;
                case 'mu'
                    value = obj.DICparamu;
                case 'beta'
                    value = obj.DICparabeta;
                case 'tol'
                    value = obj.DICparatol;
                case 'DICparagridx'
                    value = obj.DICparagridx;
                case 'DICparagridy'
                    value = obj.DICparagridy;
            end
            
        end
        
        
        % Struct zipping and unzipping
        function DICpara = zipDICPara(obj)
            
            % DIC parameters
            DICpara.winsize = obj.DICparawinsize;
            DICpara.winstepsize = obj.DICparawinstepsize;
            DICpara.gridxyROIRange.gridx = obj.DICparagridx;
            DICpara.gridxyROIRange.gridy = obj.DICparagridy;
            DICpara.LoadImgMethod = obj.DICparaLoadImgMethod;
            DICpara.ImgSeqIncUnit = obj.DICparaImgSeqIncUnit;
            DICpara.ImgSeqIncROIUpdateOrNot = obj.DICparaImgSeqIncROIUpdateOrNot;
            DICpara.Subpb2FDOrFEM = obj.DICparaSubpb2FDOrFEM;
            DICpara.NewFFTSearch = obj.DICparaNewFFTSearch;
            DICpara.ClusterNo = obj.DICparaClusterNo;
            DICpara.ImgSize = obj.DICparaImgSize;
            DICpara.mu = obj.DICparamu;
            DICpara.beta = obj.DICparabeta;
            DICpara.tol = obj.DICparatol;
            DICpara.DispFilterSize = obj.DICparaDispFilterSize;
            DICpara.DispFilterStd = obj.DICparaDispFilterStd;
            DICpara.StrainFilterSize = obj.DICparaStrainFilterSize;
            DICpara.StrainFilterStd = obj.DICparaStrainFilterStd;
            DICpara.InitFFTSearchMethod = obj.DICparaInitFFTSearchMethod;
            DICpara.SizeOfFFTSearchRegion = obj.DICparaSizeOfFFTSearchRegion;
            DICpara.ImgRefMask = obj.DICparaImgRefMask;
            DICpara.um2px = obj.DICparaum2px;
            DICpara.DoYouWantToSmoothOnceMore = obj.DICparaDoYouWantToSmoothOnceMore;
            DICpara.MethodToComputeStrain = obj.DICparaMethodToComputeStrain;
            DICpara.StrainType = obj.DICparaStrainType;
            DICpara.Image2PlotResults = obj.DICparaImage2PlotResults;
            DICpara.MethodToSaveFig = obj.DICparaMethodToSaveFig;
            DICpara.OrigDICImgTransparency = obj.DICparaOrigDICImgTransparency;
            
        end
        
        function obj = unzipDICPara(obj,DICpara)
            
            % DIC parameters
            obj.DICparawinsize = DICpara.winsize;
            obj.DICparawinstepsize = DICpara.winstepsize;
            obj.DICparagridx = DICpara.gridxyROIRange.gridx;
            obj.DICparagridy = DICpara.gridxyROIRange.gridy;
            obj.DICparaLoadImgMethod = DICpara.LoadImgMethod;
            obj.DICparaImgSeqIncUnit = DICpara.ImgSeqIncUnit;
            obj.DICparaImgSeqIncROIUpdateOrNot = DICpara.ImgSeqIncROIUpdateOrNot;
            obj.DICparaSubpb2FDOrFEM = DICpara.Subpb2FDOrFEM;
            obj.DICparaNewFFTSearch = DICpara.NewFFTSearch;
            obj.DICparaClusterNo = DICpara.ClusterNo;
            obj.DICparaImgSize = DICpara.ImgSize;
            obj.DICparamu = DICpara.mu;
            obj.DICparabeta = DICpara.beta;
            obj.DICparatol = DICpara.tol;
            obj.DICparaDispFilterSize = DICpara.DispFilterSize;
            obj.DICparaDispFilterStd = DICpara.DispFilterStd;
            obj.DICparaStrainFilterSize = DICpara.StrainFilterSize;
            obj.DICparaStrainFilterStd = DICpara.StrainFilterStd;
            obj.DICparaInitFFTSearchMethod = DICpara.InitFFTSearchMethod;
            obj.DICparaSizeOfFFTSearchRegion = DICpara.SizeOfFFTSearchRegion;
            obj.DICparaImgRefMask = DICpara.ImgRefMask;
            obj.DICparaum2px = DICpara.um2px;
            obj.DICparaDoYouWantToSmoothOnceMore = DICpara.DoYouWantToSmoothOnceMore;
            obj.DICparaMethodToComputeStrain = DICpara.MethodToComputeStrain;
            obj.DICparaStrainType = DICpara.StrainType;
            obj.DICparaImage2PlotResults = DICpara.Image2PlotResults;
            obj.DICparaMethodToSaveFig = DICpara.MethodToSaveFig;
            obj.DICparaOrigDICImgTransparency = DICpara.OrigDICImgTransparency;
        end
        
        function DICmesh = zipDICmesh(obj)
            
            try DICmesh.coordinatesFEM = obj.DICmesh_coordinatesFEM; catch;  end
            try DICmesh.elementsFEM = obj.DICmesh_elementsFEM; catch;  end
            try DICmesh.dirichlet = obj.DICmesh_dirichlet; catch;  end
            try DICmesh.neumann = obj.DICmesh_neumann; catch;  end
            try DICmesh.x0 = obj.DICmesh_x0;  catch;  end
            try DICmesh.y0 = obj.DICmesh_y0; catch;  end
            try DICmesh.M = obj.DICmesh_M;  catch;  end
            try DICmesh.N = obj.DICmesh_N;  catch;  end
            try DICmesh.y0World = obj.DICmesh_y0World; catch;  end
            try DICmesh.coordinatesFEMWorld = obj.DICmesh_coordinatesFEMWorld; catch;  end
            try DICmesh.elementMinSize = obj.DICmesh_elementMinSize;  catch;  end
        end
        
        function obj = unzipDICmesh(obj,DICmesh)
            
            try obj.DICmesh_coordinatesFEM = DICmesh.coordinatesFEM; catch;  end
            try obj.DICmesh_elementsFEM = DICmesh.elementsFEM; catch;  end
            try obj.DICmesh_dirichlet = DICmesh.dirichlet; catch;  end
            try obj.DICmesh_neumann = DICmesh.neumann; catch;  end
            try obj.DICmesh_x0 = DICmesh.x0; catch;  end
            try obj.DICmesh_y0 = DICmesh.y0; catch;  end
            try obj.DICmesh_M = DICmesh.M; catch;  end
            try obj.DICmesh_N = DICmesh.N; catch; end
            try obj.DICmesh_y0World = DICmesh.y0World; catch; end
            try obj.DICmesh_coordinatesFEMWorld = DICmesh.coordinatesFEMWorld; catch; end
            try obj.DICmesh_elementMinSize = DICmesh.elementMinSize; catch; end
        end
        
        function cc = zipcc(obj)
            try cc.Threshold = obj.ccThreshold; catch; end
            try cc.max = obj.ccmax; catch; end
            try cc.A = obj.ccA; catch; end
            try cc.qfactors = obj.ccqfactors; catch; end
            
        end
        
        function obj = unzipcc(obj,cc)
            try obj.ccThreshold = cc.Threshold; catch; end
            try obj.ccmax = cc.max; catch; end
            try obj.ccA = cc.A; catch; end
            try obj.ccqfactors = cc.qfactors; catch;  end
            
        end
        
        
        
        % helper functions
        function AddAndRemovePathFolders(obj,stepNum)
            %% ADD FUNCTIONS TO PATH IF THEY ARE PART OF THIS ALGO
            
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
        
    end
    
    methods(Static)
        
        function initMGW()
            
            % set the value of an operating system environment MINGW
            setenv('MW_MINGW64_LOC','C:\TDM-Gcc-64');
            try mex -O ba_interp2.cpp; catch; end % mex set up ba_interp2.cpp script
            
        end
        
        function [gridx,gridy] = convertROItogridpoints(ROI)
            
            bounds = regionprops(ROI,'BoundingBox');
            
            if length(bounds)== 1
                b = bounds.BoundingBox;
                gridx = round([b(1), b(1)+b(3)-1]);
                gridy = round([b(2), b(2)+b(4)-1]);
            end
        end
        
    end
    
end