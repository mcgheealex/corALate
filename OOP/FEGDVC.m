classdef FEGDVC
    %FEGDVC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Type = 'FEGDVC';
        ProgramFolder = 'FE_Global_DVC-main';
        Name = '';
        SaveFolderPath = '';
        
        % Structs
        DVCparaPack = []; % struct of many DIC parameters
        DVCmeshPack = [];
        ccPack = []; % struct of quality factors
        DfPack = [];  % struct of gradient
        
        % Computation variables
        U0 = [];
        
        
        % Result variables
        ResultDisp = [];
        ResultDefGrad = [];
        ResultStrain = [];
        ResultFEMesh = [];
        ResultAlpha = [];
        ResultNormOfW = [];
        ResultTimeICGN = [];
    end
    
    methods
        function obj = FEGDVC()
            % add the FEGDVC functions to the path and remove others
            AddAndRemovePathFolders(obj,1)
            
            initMinGW();
            obj = initStructs(obj);
            obj = initResultVars(obj);
        end
        
        function obj = initialGuess(obj,ImgRef,ImgDef,ImgSeqNum)
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This section is to find or update initial guess for ALDVC
            % The key idea is to either to use FFT peak fitting, or use last frame
            % results for the next new frame;
            % Particularly in incremental mode, the reference image can also be updated.
            % fNormalized = ImgNormalized{ImgSeqNum-mod(ImgSeqNum-1,ImgSeqIncUnit)};
            
            % this is a way to handle code updates easily
            DVCparas = obj.DVCparaPack;
            cc = obj.ccPack;
            DVCmesh = obj.DVCmeshPack;
            Df = obj.DfPack;
            
            obj.DVCpara.NewFFTSearch = 0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % check to see if this is the first image to guess or if we
            % need a new search
            if ImgSeqNum==2 || DVCparas.NewFFTSearch==1
                % ====== Integer Search ======
                [xyz0,uvw0,cc,SizeOfFFTSearchRegion] = IntegerSearch3Mg({ImgRef,ImgDef},DVCparas); 
                DVCparas.SizeOfFFTSearchRegion = SizeOfFFTSearchRegion;
                
                % ======== Find some bad inital guess points ========
                cc.ccThreshold=1.25; % bad cross-correlation threshold (mean - ccThreshold*stdev for q-factor distribution)
                DVCparas.qDICOrNot=0; 
                DVCparas.Thr0=0;
                
                % remove bad points
                [uvw,cc]=RemoveOutliers3(uvw0,cc,DVCparas.qDICOrNot,DVCparas.Thr0);%Last term is threshold value 
                
                % ====== FEM mesh set up ======
                [DVCmesh] = MeshSetUp3(xyz0,DVCparas);
                
                % ====== Assign initial values ======
                obj.U0 = Init3(uvw,DVCmesh.xyz0); 

            else
                obj.U0 = obj.ResultDisp{ImgSeqNum-2}.U;
            end
            
            % Compute image grayscale value gradients
            Df = funImgGradient3(ImgRef,'stencil7'); 
            Df.imgSize = size(ImgRef); 
            
            % FOR GUI
            % Plotdisp_show3(U0,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM);

            % repack object
            obj.DVCparaPack = DVCparas;
            obj.ccPack = cc;
            obj.DVCmeshPack = DVCmesh;
            obj.DfPack = Df;
        end
        
        function obj = globalDVC(obj,ImgRef,ImgDef,ImgSeqNum)

            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Finite element based global DVC iterations
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            alphaList = DVCpara.alpha; % Set regularization coefficient, alpha, as 100 as an example
            
            % ====== Tune regularization coefficient ======
            % If you don't know the best alpha (coefficient), please run the following
            % codes to tune the best value of the coefficient of the regularizer |grad u|^2):
            %
            % %%%%% Uncomment the following line to tune the best value of alpha %%%%%%
            alphaList = [ 1e-2,1e-1,1e0,1e1,1e2 ]*mean(DVCpara.winstepsize);
            % Err1List = zeros(length(alphaList),1); Err2List=Err1List;
            % UList = cell(length(alphaList),1); FList=UList;
            
            % ------------------------------------------------
            for alphaInd = 1:length(alphaList)
                
                tic; alpha = alphaList(alphaInd);
                % Solve displacement U with each alpha
                [U, normOfW, timeICGN] = funGlobalICGN3(DVCmesh,Df,Img{1},Img{2},U0,alpha,DVCpara.tol,DVCpara.maxIter);
                
                % Compute F deformation gradient with solved U
                DVCpara.GaussPtOrder = 2; [F,~,~] = funGlobal_NodalStrainAvg3(DVCmesh,U,DVCpara.GaussPtOrder);
                
                % UList{alphaInd} = U; FList{alphaInd} = F;
                Err1List(alphaInd) = norm(U-U0,2);
                Err2List(alphaInd) = norm(F,2);
                
            end
            
            % ====== Tune the coefficient of |grad u| regularizer ======
            ErrSumList = Err1List + 1*mean(DVCpara.winstepsize)*Err2List;  % 10 is an empirical number
            [~,indexOfalpha] = min(ErrSumList);
            % figure, plot(alphaList,ErrSumList,'o-'); set(gca,'xscale','log');
            % xlabel('\alpha'); ylabel('Err'); set(gca,'fontsize',14);
            
            try
                [fitobj] = fit(log10(alphaList(indexOfalpha-1:1:indexOfalpha+1))',ErrSumList(indexOfalpha-1:1:indexOfalpha+1),'poly2');
                p = coeffvalues(fitobj); alpha_best = 10^(-p(2)/2/p(1));
            catch
                alpha_best = alphaList(indexOfalpha);
            end
            DVCpara.alpha = alpha_best;
            
            
            % ====== Re-run global DVC iterations with tuned alpha_best ======
            if abs(alpha_best - alpha) > abs(eps)
                [U, normOfW, timeICGN] = funGlobalICGN3(DVCmesh,Df,Img{1},Img{2},U0,alpha,DVCpara.tol,DVCpara.maxIter);
                DVCpara.GaussPtOrder = 2; [F,~,~] = funGlobal_NodalStrainAvg3(DVCmesh,U,DVCpara.GaussPtOrder);
            end
            
            
            % ====== Plot solved results ======
            Plotdisp_show3(U,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM);
            Plotstrain_show3(F,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM);
            
            
            %% ====== Save data ======
            ResultDisp{ImgSeqNum-1}.U = full(U);
            ResultDefGrad{ImgSeqNum-1}.F = full(F);
            ResultAlpha{ImgSeqNum-1}.alpha = alpha_best;
            ResultNormOfW{ImgSeqNum-1}.normOfW = full(normOfW);
            ResultTimeICGN{ImgSeqNum-1}.timeICGN = full(timeICGN);
            
            fprintf('------------ Section 4 Done ------------ \n \n')
            
            
        end
        
        function obj = imagesetup(obj)          
            [ImgNormalized,obj.DVCparaPack.gridRange] = funNormalizeImg3(Img,obj.DVCpara.gridRange,'Normalize');
        end
        
        function obj = initResultVars(obj)

            obj.ResultDisp = cell(length(ImgNormalized)-1,1);
            obj.ResultDefGrad = cell(length(ImgNormalized)-1,1);
            obj.ResultStrain = cell(length(ImgNormalized)-1,1);
            obj.ResultFEMesh = cell(ceil((length(ImgNormalized)-1)/obj.DVCpara.ImgSeqIncUnit),1); % For incremental DIC mode
            obj.ResultAlpha = cell(length(ImgNormalized)-1,1);
            obj.ResultNormOfW = cell(length(ImgNormalized)-1,1);
            obj.ResultTimeICGN = cell(length(ImgNormalized)-1,1);
        end
        
        function obj = initStructs(obj)
            obj.DVCparaPack.interpmethod = 'cubic';   % Grayscale interpolation scheme: choose from {'linear','cubic','spline','default'}
            obj.DVCparaPack.displayIterOrNot = 0;     % Display Section 4 Subpb1 IC-GN convergence info
            obj.DVCparaPack.maxIter = 100;            % Maximum IC-GN iterations in IC-GN iterations
            obj.DVCparaPack.tol = 1e-1;               % Iteration stopping threshold
            obj.DVCparaPack.alpha = 1e2;
            obj.DVCparaPack.winsize = [];
            obj.DVCparaPack.winstepsize = [];
            obj.DVCparaPack.gridRange = [];
            obj.DVCparaPack.ClusterNo = 1;
            obj.DVCparaPack.ImgSize = [0 0 ];
            obj.DVCparaPack.ImgSeqIncUnit = [];
            obj.DVCparaPack.ImgSeqIncROIUpdateOrNot = [];
            obj.DVCparaPack.InitFFTMethod = [];
            obj.DVCparaPack.NewFFTSearch = [];
            obj.DVCparaPack.DIM = 3;
            obj.DVCparaPack.NewFFTSearch = 0;
            obj.DVCparaPack.qDICOrNot=0; 
            obj.DVCparaPack.Thr0=0; 
            obj.DVCparaPack.SizeOfFFTSearchRegion=[];
            
            obj.ccPack.A = {};
            obj.ccPack.max = 0;
            obj.ccPack.qfactors = [];
            obj.ccPack.ccThreshold = 1.25;
            
            obj.DVCmeshPack.coordinatesFEM = [];
            obj.DVCmeshPack.elementsFEM = [];
            obj.DVCmeshPack.dirichlet = [];
            obj.DVCmeshPack.neumann = [];
            obj.DVCmeshPack.xyz0 = [];
            
            obj.DfPack.DfAxis = []; 
            obj.DfPack.imgSize = [];
            obj.DfPack.DfDx = []; 
            obj.DfPack.DfDy = []; 
            obj.DfPack.DfDz = [];
        end
        

        

        
    end
    
    methods(Static)
        
        function initMinGW()
            %% Section 1
            setenv('MW_MINGW64_LOC','C:\TDM-GCC-64')
            mex -O ba_interp3.cpp; warning('off'); % dbstop if error
        end
        
    end
end

