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
    % ReadImage -> InitDICoptions 
    % IntegerSearch -> InitSearchSize
    
    properties
        ProgramFolder = '2D_ALDIC-master';
        Name = '';
        SaveFolderPath = '';
        
        % Raw Image Data
        Images = {};
        RefImages = [];
        % Raw ROI mask Data
        ROIMasks = {};
        
        % DIC parameters
        DICparawinsize = 50;
        DICparawinstepsize = 5;
        DICparagridx = [];
        DICparagridy = [];
        DICparaLoadImgMethod = 0;
        DICparaImgSeqIncUnit = 1;
        DICparaImgSeqIncROIUpdateOrNot = 0;
        DICparaSubpb2FDOrFEM = 0;
        DICparaNewFFTSearch = 0;
        DICparaClusterNo = 1;
        DICparaImgSize = [];
        DICparamu = 1e-3;
        DICparabeta = 0;
        DICparatol = 1e-2;
        DICparaDispFilterSize=0; 
        DICparaDispFilterStd=0; 
        DICparaStrainFilterSize=0; 
        DICparaStrainFilterStd=0;
        DICparaInitFFTSearchMethod = 0;
        DICparaSizeOfFFTSearchRegion = 0;
        
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
        ResultStrainWorld = {};
        ResultStressWorld = {};

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
        
        function obj = saveALDIC(obj,filename)
            %save the resulting data into a .mat file for future use
            
        end
        
        
        % The initial guess steps FFT based cross correlation
        % These functions still need revision for use with GUI
        %  user should be able to iteratively update the search window size
        %  user should be able to view the results of the initial guess
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
            
            % get the DIC parameters
            DICpara = zipDICPara(obj);
            
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
            % get the DIC parameters
            DICpara = zipDICPara(obj);
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
            
            %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                 % Plotting initial guess
            %                 % --------------------------------------
            %                 close all;
            %                 figure;
            %                 surf(u);
            %                 colorbar;
            %                 title('Displacement u','fontweight','normal')
            %                 set(gca,'fontSize',18);
            %                 title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');
            %                 axis tight; %axis equal; % set(gca,'XTick',[] );
            %                 xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            %                 set(gcf,'color','w');
            %                 a = gca;
            %                 a.TickLabelInterpreter = 'latex';
            %                 b = colorbar;
            %                 b.TickLabelInterpreter = 'latex';
            %                 box on; colormap jet;
            %
            %                 figure;
            %                 surf(v);
            %                 colorbar;
            %                 title('Displacement v','fontweight','normal')
            %                 set(gca,'fontSize',18);
            %                 title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
            %                 axis tight; %axis equal; % set(gca,'XTick',[] );
            %                 xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            %                 set(gcf,'color','w');
            %                 a = gca;
            %                 a.TickLabelInterpreter = 'latex';
            %                 b = colorbar;
            %                 b.TickLabelInterpreter = 'latex';
            %                 box on; colormap jet;
            %
            
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
            obj.ResultStrainWorld = cell(ImgCount-1,1);  
            obj.ResultStressWorld = cell(ImgCount-1,1);
            obj.ResultFEMeshEachFrame_coordinatesFEM = cell(ImgCount-1,1);
            obj.ResultFEMeshEachFrame_elementsFEM = cell(ImgCount-1,1);
            obj.ResultFEMesh_coordinatesFEM = cell(ceil((ImgCount-1)/obj.DICparaImgSeqIncUnit),1); % For incremental DIC mode
            obj.ResultFEMesh_elementsFEM = obj.ResultFEMesh_coordinatesFEM;
            obj.ResultFEMesh_winsize = obj.ResultFEMesh_coordinatesFEM;
            obj.ResultFEMesh_winstepsize = obj.ResultFEMesh_coordinatesFEM;
            obj.ResultFEMesh_gridx = obj.ResultFEMesh_coordinatesFEM;
            obj.ResultFEMesh_gridy = obj.ResultFEMesh_coordinatesFEM;
            
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
            DICpara.DICparaInitFFTSearchMethod = obj.DICparaInitFFTSearchMethod;
            DICpara.DICparaSizeOfFFTSearchRegion = obj.DICparaSizeOfFFTSearchRegion;

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
            obj.DICparaInitFFTSearchMethod = DICpara.DICparaInitFFTSearchMethod;
            obj.DICparaSizeOfFFTSearchRegion = DICpara.DICparaSizeOfFFTSearchRegion;

        end
        
        function DICmesh = zipDICmesh(obj)
            
            DICmesh.coordinatesFEM = obj.DICmesh_coordinatesFEM;
            DICmesh.elementsFEM = obj.DICmesh_elementsFEM;
            DICmesh.dirichlet = obj.DICmesh_dirichlet;
            DICmesh.neumann = obj.DICmesh_neumann;
            DICmesh.x0 = obj.DICmesh_x0; 
            DICmesh.y0 = obj.DICmesh_y0;
            DICmesh.M = obj.DICmesh_M; 
            DICmesh.N = obj.DICmesh_N; 
            DICmesh.y0World = obj.DICmesh_y0World;
            DICmesh.coordinatesFEMWorld = obj.DICmesh_coordinatesFEMWorld;
        
        end
        
        function obj = unzipDICmesh(obj,DICmesh)
            
            obj.DICmesh_coordinatesFEM = DICmesh.coordinatesFEM;
            obj.DICmesh_elementsFEM = DICmesh.elementsFEM;
            obj.DICmesh_dirichlet = DICmesh.dirichlet;
            obj.DICmesh_neumann = DICmesh.neumann;
            obj.DICmesh_x0 = DICmesh.x0; 
            obj.DICmesh_y0 = DICmesh.y0;
            obj.DICmesh_M = DICmesh.M; 
            obj.DICmesh_N = DICmesh.N; 
            obj.DICmesh_y0World = DICmesh.y0World;
            obj.DICmesh_coordinatesFEMWorld = DICmesh.coordinatesFEMWorld;
        
        end
        
        function cc = zipcc(obj)
            try cc.Threshold = obj.ccThreshold; catch; cc.ccThreshold = []; end
            try cc.max = obj.ccmax; catch; cc.ccmax = []; end
            try cc.A = obj.ccA; catch; cc.ccA = []; end
            try cc.qfactors = obj.ccqfactors; catch; cc.ccqfactors = []; end
            
        end
        
        function obj = unzipcc(obj,cc)
            try obj.ccThreshold = cc.Threshold; catch; obj.ccThreshold = []; end
            try obj.ccmax = cc.max; catch; obj.ccmax = []; end
            try obj.ccA = cc.A; catch; obj.ccA = []; end
            try obj.ccqfactors = cc.qfactors; catch; obj.ccqfactors = []; end
            
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