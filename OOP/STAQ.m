classdef STAQ
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
        Type = 'ALDIC_STAQ';
        ProgramFolder = '2D-ALDIC-Adapt';
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
        DICparaImgDefMask = [];
        DICparaSubpb2FDOrFEM = 0;
        DICparaNewFFTSearch = 0;
        DICparaClusterNo = 0;
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
        DICparawinsizeMin = 4;
        DICparawinsize_List = [];
        DICparawinsize_opt = 'constant';
        DICparaDispSmoothness = 0;
        DICparaStrainSmoothness = 1e-4;
            
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
        DICmesh_markCoordHoleEdge = [];
        DICmesh_markCoordHoleStrain = [];
        DICmesh_irregular = [];
        
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
        mu = 1e-3;
        beta = [];
        
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
        ResultFEMesh_markCoordHoleEdge = {};
        ResultFEMesh_elementMinSize = {};
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
        function obj = STAQ()
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
            if ImgSeqNum < 7
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
                [obj,~, x0temp_g, y0temp_g, u_g, v_g, cc] = InitSearchSize(obj,ImgRef,ImgDef);
                obj.u_f = -u_g;
                obj.v_f = -v_g;
                obj.x0temp_f = x0temp_g+u_g;
                obj.y0temp_f = y0temp_g+v_g;
                obj = unzipcc(obj,cc);
                
                %%%%% Interpolate to f %%%%%
                xnodes = max([4+0.5*DICpara.winsize,DICpara.gridxyROIRange.gridx(1)])  ...
                    : DICpara.winstepsize : min([size(ImgRef,1)-0.5*DICpara.winsize-3,DICpara.gridxyROIRange.gridx(2)]);
                
                ynodes = max([4+0.5*DICpara.winsize,DICpara.gridxyROIRange.gridy(1)])  ...
                    : DICpara.winstepsize : min([size(ImgRef,2)-0.5*DICpara.winsize-3,DICpara.gridxyROIRange.gridy(2)]);
                
                [x0temp,y0temp] = ndgrid(xnodes,ynodes);
                u_f_NotNanInd = find(~isnan(obj.u_f(:)));
                
                op1 = rbfcreate( [obj.x0temp_f(u_f_NotNanInd),obj.y0temp_f(u_f_NotNanInd)]',(obj.u_f(u_f_NotNanInd))','RBFFunction', 'thinplate'); rbfcheck(op1);
                obj.u = rbfinterp( [x0temp(:),y0temp(:)]', op1 );
                op2 = rbfcreate( [obj.x0temp_f(u_f_NotNanInd),obj.y0temp_f(u_f_NotNanInd)]',(obj.v_f(u_f_NotNanInd))','RBFFunction', 'thinplate'); rbfcheck(op2);
                obj.v = rbfinterp([x0temp(:),y0temp(:)]', op2 );
%                 x0temp = x0temp'; y0temp = y0temp'; obj.u=obj.u'; obj.v=obj.v';
                
                %%%%% Do some regularization to further decrease the noise %%%%%
                obj.u = regularizeNd([x0temp(:),y0temp(:)],obj.u(:),{xnodes',ynodes'},1e-3);
                obj.v = regularizeNd([x0temp(:),y0temp(:)],obj.v(:),{xnodes',ynodes'},1e-3);
                
                % ====== I create the DIC mesh ====== I
                [DICMESH] = MeshSetUp(x0temp,y0temp,DICpara);
                obj = unzipDICmesh(obj,DICMESH);
                % zip all the DIC mesh parameters into a struct
                DICmesh = zipDICmesh(obj);
                
                
                % ====== Initial Value ======
                obj.U0 = Init(obj.u,obj.v,cc.max,DICmesh.x0,DICmesh.y0,0);
                
                % init U0 with nans
                for tempi = 1:size(obj.u,1)
                    
                    for tempj = 1:size(obj.u,2)
                        try
                            if ~obj.DICparaImgRefMask(x0temp(tempi,tempj),y0temp(tempi,tempj))  || ...
                                   ~obj.DICparaImgDefMask( floor(x0temp(tempi,tempj)+obj.u(tempi,tempj)), floor(y0temp(tempi,tempj)+obj.v(tempi,tempj)) )
                                
                                obj.U0(2*(tempj+(tempi-1)*(size(obj.u,2)))) = nan;
                                obj.U0(2*(tempj+(tempi-1)*(size(obj.u,2)))-1) = nan;
                                
                            end
                        catch
                        end
                    end
                end
                
                % ====== Deal with incremental mode ======
                ImgRefNewIndex = ImgSeqNum-mod(ImgSeqNum-2,obj.DICparaImgSeqIncUnit)-1;
                
                if obj.DICparaImgSeqIncUnit == 1, ImgRefNewIndex = ImgRefNewIndex-1; end
                
                obj.ResultFEMesh_coordinatesFEM{1+floor(ImgRefNewIndex/obj.DICparaImgSeqIncUnit)} = obj.DICmesh_coordinatesFEM;
                obj.ResultFEMesh_elementsFEM{1+floor(ImgRefNewIndex/obj.DICparaImgSeqIncUnit)} = obj.DICmesh_elementsFEM;
                obj.ResultFEMesh_winsize{1+floor(ImgRefNewIndex/obj.DICparaImgSeqIncUnit)} = obj.DICparawinsize;
                obj.ResultFEMesh_winstepsize{1+floor(ImgRefNewIndex/obj.DICparaImgSeqIncUnit)} = obj.DICparawinstepsize;
                obj.ResultFEMesh_gridx{1+floor(ImgRefNewIndex/obj.DICparaImgSeqIncUnit)} = obj.DICparagridx;
                obj.ResultFEMesh_gridy{1+floor(ImgRefNewIndex/obj.DICparaImgSeqIncUnit)} = obj.DICparagridy;
                
                % ====== Generate a quadtree mesh considering sample's complex geometry   83
                DICmesh.elementMinSize = DICpara.winsizeMin; % min element size in the refined quadtree mesh
                
                % Generate a quadtree mesh
                obj = GenerateQuadtreeMesh(obj,ImgSeqNum); 
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CODE IS BROKEN HERE
                
                % zip all the DIC options into a struct
                DICpara = zipDICPara(obj);
                % zip all the DIC mesh parameters into a struct
                DICmesh = zipDICmesh(obj);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %   Incremental mode not established in this code yet
                %
                %                         elseif mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) == 0 % To update ref image in incremental mode
                %                             ImgRefNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
                %                             if DICpara.ImgSeqIncUnit == 1,  ImgRefNewIndex = ImgRefNewIndex-1; end
                %                             ImgRef = ImgNormalized{ImgRefNewIndex}; % Update reference
                %                             [DICpara,DICmesh] = ReadImageRefUpdate(file_name,ImgSeqNum,ResultDisp{ImgSeqNum-2}.U,DICpara,DICmesh); % Update reference image if needed;
                %                             U0 = zeros(2*size(DICmesh.coordinatesFEM,1),1); % [Temporary code: " PlotuvInit; "]
                %                             ResultFEMesh{1+floor(ImgRefNewIndex/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info
                %                                 struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ...
                %                                 'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else % Use the solved results from the last frame as the new initial guess
                if ImgSeqNum < 7 % Import previous U for ImgSeqNum [2,6]
                    obj.U0 = obj.ResultDisp{ImgSeqNum-2};
                    
                else % When ImgSeqNum > 6: POD predicts next disp U0 from previous results of (ImgSeqNum+[-5:1:-1])
                    nTime = 5;
                    np = length(obj.ResultDisp{ImgSeqNum-2})/2; % "nTime" value 5 is an empirical value, can be changed.
                    T_data_u = zeros(nTime,np);
                    T_data_v = zeros(nTime,np);
                    
                    for tempi = 1:nTime
                        T_data_u(tempi,:) = obj.ResultDisp{ImgSeqNum-(2+nTime)+tempi, 1}(1:2:np*2)';
                        T_data_v(tempi,:) = obj.ResultDisp{ImgSeqNum-(2+nTime)+tempi, 1}(2:2:np*2)';
                    end
                    
                    nB = 3;
                    t_train = [ImgSeqNum-1-nTime:ImgSeqNum-2]';
                    t_pre = [ImgSeqNum-1]';
                    [u_pred,~,~,~] = funPOR_GPR(T_data_u,t_train,t_pre,nB);
                    [v_pred,~,~,~] = funPOR_GPR(T_data_v,t_train,t_pre,nB);
                    tempu = u_pred(1,:); tempv = v_pred(1,:);
                    obj.U0 = [tempu(:),tempv(:)]';
                    obj.U0 = obj.U0(:);
                    
                end
            end
            
            % ====== Compute f(X)-g(x+u) ======
            obj.ResultFEMeshEachFrame_coordinatesFEM{ImgSeqNum-1} = obj.DICmesh_coordinatesFEM;
            obj.ResultFEMeshEachFrame_elementsFEM{ImgSeqNum-1} = obj.DICmesh_elementsFEM;
            
            % save current state of DIC mesh and para
            obj = unzipDICmesh(obj,DICmesh);
            obj = unzipDICPara(obj,DICpara);
            
        end
           
        function obj = GenerateQuadtreeMesh(obj,ImgSeqNum)
            % Generate a quadtree mesh considering sample's complex geometry
            %
            % ----------------------------------------------
            % References
            % [1] J Yang, K Bhattacharya. Fast adaptive mesh augmented Lagrangian Digital Image
            % Correlation. Under review.
            % [2] S Funken, A Schmidt. Adaptive mesh refinement in 2D: an efficient
            % implementation in MATLAB. Comp. Meth. Appl. Math. 20:459-479, 2020.
            % [3] Rbfinterp. Matlab File Exchange open source.
            % https://www.mathworks.com/matlabcentral/fileexchange/10056-scattered-data-interpolation-and-approximation-using-radial-base-functions
            % ----------------------------------------------
            % Author: Jin Yang
            % Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
            % Last time updated: 02/2020.
            % ==============================================
            
            % zip all the DIC options into a struct
            DICpara = zipDICPara(obj);
            % zip all the DIC mesh parameters into a struct
            DICmesh = zipDICmesh(obj);
            
            disp('--- Generate a quadtree mesh ---')
            %% ====== Remove finite elements where there is a hole ======
            coordinatesFEMQuadtree = DICmesh.coordinatesFEM;
            elementsFEMQuadtree = DICmesh.elementsFEM(:,1:4);
            irregular = zeros(0,3);
            
            while 1 % Generate a Quadtree mesh
                [~,mark4] = funMarkEdge(coordinatesFEMQuadtree,elementsFEMQuadtree,DICpara.ImgRefMask,DICmesh.elementMinSize); % Don't delete "*2"
                mark4 = find(mark4);
                [coordinatesFEMQuadtree,elementsFEMQuadtree,irregular] = QrefineR(coordinatesFEMQuadtree,elementsFEMQuadtree,irregular,mark4);
                if isempty(mark4)
                    break
                end
            end
            
            %%%%% Re-order node index in elements %%%%%
            for tempj = 1:size(elementsFEMQuadtree,1)
                coordxAll = coordinatesFEMQuadtree(elementsFEMQuadtree(tempj,1:4),1);
                coordyAll = coordinatesFEMQuadtree(elementsFEMQuadtree(tempj,1:4),2);
                coordxAll_sorted = sort(coordxAll(:));
                coordyAll_sorted = sort(coordyAll(:));
                
                temp1 = (coordxAll-coordxAll_sorted(1)).^2 + 2*((coordyAll-coordyAll_sorted(1))).^2 ;
                [~,obj.temp3] = sort(temp1);
                
                elementsFEMQuadtree(tempj,1:4) = elementsFEMQuadtree(tempj,obj.temp3([1,2,4,3]));
                
                % coordinatesFEMQuadtree(elementsFEMQuadtree(tempj,1:4),1:2)
            end
            
            
            %% %%%%% Plot refined mesh %%%%%
            % figure; patch('Faces', elementsFEMQuadtree(:,1:4), 'Vertices', coordinatesFEMQuadtree, 'Facecolor','none','linewidth',1)
            % axis equal; axis tight; set(gca,'fontsize',18); set(gcf,'color','w'); box on;
            % xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
            % title('Quadtree mesh','Interpreter','latex');
            % a = gca; a.TickLabelInterpreter = 'latex';
            
            % Update the quadtree mesh to deal with hanging nodes
            for tempj = 1:size(irregular,1)
                [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,1:2), 'rows' );
                if Lia>0, elementsFEMQuadtree(Locb,8)=irregular(tempj,3); end
                [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,2:3), 'rows' );
                if Lia>0, elementsFEMQuadtree(Locb,5)=irregular(tempj,3); end
                [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,3:4), 'rows' );
                if Lia>0, elementsFEMQuadtree(Locb,6)=irregular(tempj,3); end
                [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,[4,1]), 'rows' );
                if Lia>0, elementsFEMQuadtree(Locb,7)=irregular(tempj,3); end
                
                [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,[2,1]), 'rows' );
                if Lia>0, elementsFEMQuadtree(Locb,8)=irregular(tempj,3); end
                [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,[3,2]), 'rows' );
                if Lia>0, elementsFEMQuadtree(Locb,5)=irregular(tempj,3); end
                [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,[4,3]), 'rows' );
                if Lia>0, elementsFEMQuadtree(Locb,6)=irregular(tempj,3); end
                [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,[1,4]), 'rows' );
                if Lia>0, elementsFEMQuadtree(Locb,7)=irregular(tempj,3); end
            end
            
            % Remove elements within the center hole
            [~,markOutside4] = funMarkInside(coordinatesFEMQuadtree,elementsFEMQuadtree,DICpara.ImgRefMask);
            elementsFEMQuadtree = elementsFEMQuadtree(markOutside4,:);
            
            
            %% Find nodes near the edges
            
            % %%%%% Old codes: Erode ref image mask, and to find elements near holes' edges,
            % nhood = [1 1 1; 1 1 1; 1 1 1];
            % ImgRefMaskErode = DICpara.ImgRefMask;
            % for tempi = 1: floor(0.5*max([20,mean(DICpara.winsize),mean(DICpara.winstepsize)]))-1
            %     ImgRefMaskErode = imerode(ImgRefMaskErode, nhood);
            % end
            % % figure, imshow(ImgRefMaskErode');
            % [markEleHoleEdge4,markEleFarOutside4] = funMarkInside(coordinatesFEMQuadtree,elementsFEMQuadtree,ImgRefMaskErode);
            
            % %%%%% New codes: Find elements which are refined %%%%%%
            elementsFEMQuadtreeSize = sqrt( ( coordinatesFEMQuadtree(elementsFEMQuadtree(:,1),1) - coordinatesFEMQuadtree(elementsFEMQuadtree(:,3),1) ).^2 + ...
                ( coordinatesFEMQuadtree(elementsFEMQuadtree(:,1),2) - coordinatesFEMQuadtree(elementsFEMQuadtree(:,3),2) ).^2 );
            [markEleRefine4,~] = find(elementsFEMQuadtreeSize <  0.99*sqrt(2)*max([DICpara.winstepsize,0*DICpara.winsize]));
            
            % %%%%% New codes: Find elements near the boudary %%%%%%
            xMin = min(DICmesh.coordinatesFEM(:,1)); xMax = max(DICmesh.coordinatesFEM(:,1));
            yMin = min(DICmesh.coordinatesFEM(:,2)); yMax = max(DICmesh.coordinatesFEM(:,2));
            [row1,~] = find(coordinatesFEMQuadtree(elementsFEMQuadtree(:,1),1) < xMin+1.01*DICpara.winstepsize);
            [row2,~] = find(coordinatesFEMQuadtree(elementsFEMQuadtree(:,3),1) > xMax-1.01*DICpara.winstepsize);
            [row3,~] = find(coordinatesFEMQuadtree(elementsFEMQuadtree(:,1),2) < yMin+1.01*DICpara.winstepsize);
            [row4,~] = find(coordinatesFEMQuadtree(elementsFEMQuadtree(:,3),2) > yMax-1.01*DICpara.winstepsize);
            
            markEleHoleEdge4 =  union(row4,union(row3,union(row2,union(row1,markEleRefine4))));
            markCoordHoleEdge = unique(elementsFEMQuadtree(markEleHoleEdge4,:));
            try
                if markCoordHoleEdge(1)==0, markCoordHoleEdge = markCoordHoleEdge(2:end); end
            catch
            end
            
            %%%%%% New codes: Find elements near marked elements %%%%%%
            for tempi = 1:2 % 2+(round( 32 / mean(DICpara.winstepsize) )^2)
                
                markEleHoleEdgeNeigh4 = zeros(size(elementsFEMQuadtree,1),1);
                for eleInd = 1:size(elementsFEMQuadtree,1)
                    markEleHoleEdgeNeigh4(eleInd) = length(intersect(elementsFEMQuadtree(eleInd,:),markCoordHoleEdge));
                end
                [markEleHoleEdgeNeigh4,~] = find(markEleHoleEdgeNeigh4>0);
                %%%%%%%%%
                markCoordHoleEdge = unique(elementsFEMQuadtree(markEleHoleEdgeNeigh4,:)) ;
                try
                    if markCoordHoleEdge(1) == 0, markCoordHoleEdge = markCoordHoleEdge(2:end); end
                catch
                end
                
            end
            
            
            % %%%%% Store data structure %%%%%
            DICmesh.markCoordHoleEdge = markCoordHoleEdge;
            DICmesh.dirichlet = DICmesh.markCoordHoleEdge;
            
                                                                                                                                                        %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                                                                                                        %             % %%%%% Plot %%%%%
                                                                                                                                                        %             figure;
                                                                                                                                                        %             patch('Faces', elementsFEMQuadtree(:,1:4), 'Vertices', coordinatesFEMQuadtree, 'Facecolor','white','linewidth',1);
                                                                                                                                                        %             patch('Faces', elementsFEMQuadtree(markEleHoleEdge4,1:4), 'Vertices', coordinatesFEMQuadtree, 'Facecolor','yellow','linewidth',1);
                                                                                                                                                        %             hold on; patch('Faces', elementsFEMQuadtree(markEleHoleEdgeNeigh4,1:4), 'Vertices', coordinatesFEMQuadtree, 'Facecolor','yellow','linewidth',1);
                                                                                                                                                        %             axis equal; axis tight; set(gca,'fontsize',18); set(gcf,'color','w'); box on;
                                                                                                                                                        %             xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
                                                                                                                                                        %             title('Quadtree mesh','Interpreter','latex');
                                                                                                                                                        %             a = gca; a.TickLabelInterpreter = 'latex';
                                                                                                                                                        %             
                                                                                                                                                        %             lgd = legend('Quadtree mesh elements','Elements near the edge','interpreter','latex','location','northeastoutside');
                                                                                                                                                        %             set(lgd,'fontsize',13);
                                                                                                                                                        %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %% Initialize variable U for the generated quadtree mesh
            U0NotNanInd = find(~isnan(obj.U0(1:2:end)));
            U0Quadtree = 0*coordinatesFEMQuadtree(:);
            
            dilatedI = ( imgaussfilt(double(obj.Df.ImgRefMask),1) );
            dilatedI = logical( dilatedI > 0.01);
            cc = bwconncomp(dilatedI,8);
            indPxAll = sub2ind( obj.Df.imgSize, coordinatesFEMQuadtree(:,1), coordinatesFEMQuadtree(:,2) );
            indPxNotNanAll = sub2ind( obj.Df.imgSize, DICmesh.coordinatesFEM(U0NotNanInd,1), DICmesh.coordinatesFEM(U0NotNanInd,2) );
            stats = regionprops(cc,'Area','PixelList');
            
            for tempi = 1:length(stats)
                
                try
                    
                    %%%%% Find those nodes %%%%%
                    indPxtempi = sub2ind( obj.Df.imgSize, stats(tempi).PixelList(:,2), stats(tempi).PixelList(:,1) );
                    Lia = ismember(indPxAll,indPxtempi); [LiaList,~] = find(Lia==1);
                    Lib = ismember(indPxNotNanAll,indPxtempi); [LibList,~] = find(Lib==1);
                    
                    %%%%% RBF (Radial basis function) works better than "scatteredInterpolant" %%%%%
                    % ------ Disp u ------
                    op1 = rbfcreate( [DICmesh.coordinatesFEM(U0NotNanInd(LibList),1:2)]',[obj.U0(2*U0NotNanInd(LibList)-1)]','RBFFunction', 'thinplate'); rbfcheck(op1);
                    fi1 = rbfinterp( [coordinatesFEMQuadtree(LiaList,1:2)]', op1);
                    U0Quadtree(2*LiaList-1) = fi1(:);
                    
                    % ------ Disp v ------
                    op1 = rbfcreate( [DICmesh.coordinatesFEM(U0NotNanInd(LibList),1:2)]',[obj.U0(2*U0NotNanInd(LibList))]','RBFFunction', 'thinplate'); rbfcheck(op1);
                    fi1 = rbfinterp( [coordinatesFEMQuadtree(LiaList,1:2)]', op1);
                    U0Quadtree(2*LiaList) = fi1(:);
                    
                catch
                end
            end
            
            obj.U0 = U0Quadtree;
            
            % F_dispu = scatteredInterpolant( DICmesh.coordinatesFEM(U0NotNanInd,1),DICmesh.coordinatesFEM(U0NotNanInd,2),U0(2*U0NotNanInd-1),'linear','linear' );
            % F_dispv = scatteredInterpolant( DICmesh.coordinatesFEM(U0NotNanInd,1),DICmesh.coordinatesFEM(U0NotNanInd,2),U0(2*U0NotNanInd),'linear','linear' );
            %
            % U0 = 0*coordinatesFEMQuadtree(:);
            % temp = F_dispu(coordinatesFEMQuadtree(:,1),coordinatesFEMQuadtree(:,2)); U0(1:2:end)=temp(:);
            % temp = F_dispv(coordinatesFEMQuadtree(:,1),coordinatesFEMQuadtree(:,2)); U0(2:2:end)=temp(:);
            
%             Plotdisp_show( -obj.U0,coordinatesFEMQuadtree,elementsFEMQuadtree(:,1:4),DICpara,'NoEdgeColor');
            
            if size(elementsFEMQuadtree,2)<8 % make sure the column# of elementsFEMQuadtree is 8
                elementsFEMQuadtree(:,8) = 0*elementsFEMQuadtree(:,1);
            end
            
            DICmesh.coordinatesFEM = coordinatesFEMQuadtree;
            DICmesh.elementsFEM = elementsFEMQuadtree;
            DICmesh.irregular = irregular;
            DICmesh.coordinatesFEMWorld = [DICmesh.coordinatesFEM(:,1),size(DICpara.ImgRefMask,2)+1-DICmesh.coordinatesFEM(:,2)];
            
            % ====== Deal with incremental mode ======
            ImgRefNewIndex = ImgSeqNum-mod(ImgSeqNum-2,obj.DICparaImgSeqIncUnit)-1;
            
            if obj.DICparaImgSeqIncUnit == 1, ImgRefNewIndex = ImgRefNewIndex-1; end
            
            obj.ResultFEMeshEachFrame_coordinatesFEM{1+floor(ImgRefNewIndex/obj.DICparaImgSeqIncUnit)} = DICmesh.coordinatesFEM;
            obj.ResultFEMeshEachFrame_elementsFEM{1+floor(ImgRefNewIndex/obj.DICparaImgSeqIncUnit)} = DICmesh.elementsFEM;
            obj.ResultFEMesh_coordinatesFEM{1+floor(ImgRefNewIndex/obj.DICparaImgSeqIncUnit)} = {};
            obj.ResultFEMesh_elementsFEM{1+floor(ImgRefNewIndex/obj.DICparaImgSeqIncUnit)} = {};
            obj.ResultFEMesh_winsize{1+floor(ImgRefNewIndex/obj.DICparaImgSeqIncUnit)} = DICpara.winsize;
            obj.ResultFEMesh_winstepsize{1+floor(ImgRefNewIndex/obj.DICparaImgSeqIncUnit)} =DICpara.winstepsize;
            obj.ResultFEMesh_gridx{1+floor(ImgRefNewIndex/obj.DICparaImgSeqIncUnit)} = DICpara.gridxyROIRange.gridy;
            obj.ResultFEMesh_gridy{1+floor(ImgRefNewIndex/obj.DICparaImgSeqIncUnit)} = DICpara.gridxyROIRange.gridy;
            obj.ResultFEMesh_markCoordHoleEdge{1+floor(ImgRefNewIndex/obj.DICparaImgSeqIncUnit)} = DICmesh.markCoordHoleEdge;
            obj.ResultFEMesh_elementMinSize{1+floor(ImgRefNewIndex/obj.DICparaImgSeqIncUnit)} = DICmesh.elementMinSize;
            
            % ===== Remove bad points =====
            % JY!!!
%             [obj.U0,~] = funRemoveOutliersQuadtree(DICmesh,DICpara,obj.U0,[obj.U0;obj.U0]);
            
            % save current state of DIC mesh and para
            obj = unzipDICmesh(obj,DICmesh);
            obj = unzipDICPara(obj,DICpara);

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
            
             % ------ Start Local DIC IC-GN iteration ------
            [obj.USubpb1{1},obj.FSubpb1{1},obj.HtempPar,~,~,~,DICmesh.markCoordHoleStrain] = ...
                LocalICGNQuadtree(obj.U0,DICmesh.coordinatesFEM,obj.Df,ImgRef,ImgDef,DICpara,'GaussNewton',DICpara.tol);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ------  Manually find some bad points from Local Subset ICGN step ------
            % Comment these lines below if you don't have local bad points
            % %%%%% Comment START %%%%%%
            % USubpb1(2*DICmesh.markCoordHoleEdge-1:2*DICmesh.markCoordHoleEdge) = nan;
            % FSubpb1(4*DICmesh.markCoordHoleEdge-3:4*DICmesh.markCoordHoleEdge) = nan;
            % [USubpb1,FSubpb1] = funRemoveOutliersQuadtree(DICmesh,DICpara,USubpb1,FSubpb1);
            % disp('--- Remove bad points done ---')
            % %%%%% Comment END %%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % ====== Thin-plate interpolate bad points =====
            coordinatesFEM = DICmesh.coordinatesFEM;
            U = obj.USubpb1{1};
            F = obj.FSubpb1{1};
            nanindexU = find(isnan(U(1:2:end))==1); notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindexU);
            nanindexF = find(isnan(F(1:4:end))==1); notnanindexF = setdiff([1:1:size(coordinatesFEM,1)],nanindexF);
            
            %%%%% Ux %%%%%%
            op1 = rbfcreate( [coordinatesFEM(notnanindex,1:2)]',[U(2*notnanindex-1)]','RBFFunction', 'thinplate');
            rbfcheck_maxdiff = rbfcheck(op1); % check rbf interpolation
            if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
            fi1 = rbfinterp([coordinatesFEM(:,1:2)]', op1 );
            
            %%%%% Uy %%%%%%
            op2 = rbfcreate( [coordinatesFEM(notnanindex,1:2)]',[U(2*notnanindex)]','RBFFunction', 'thinplate');
            rbfcheck_maxdiff = rbfcheck(op2); % check rbf interpolation
            if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
            fi2 = rbfinterp([coordinatesFEM(:,1:2)]', op2 );
            
            %%%%% Assemble [Ux, Uy] %%%%%
            U_rbf_thinplate = [fi1(:),fi2(:)]';  U_rbf_thinplate = U_rbf_thinplate(:);
            
            %%%%% F11 %%%%%
            op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
                F(4*notnanindex-3)','RBFFunction', 'thinplate');
            rbfcheck_maxdiff = rbfcheck(op);
            if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
            fi11 = rbfinterp([coordinatesFEM(:,1:2)]', op );
            
            %%%%% F21 %%%%%
            op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
                F(4*notnanindex-2)','RBFFunction', 'thinplate');
            rbfcheck_maxdiff = rbfcheck(op);
            if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
            fi21 = rbfinterp([coordinatesFEM(:,1:2)]', op );
            
            %%%%% F12 %%%%%
            op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
                F(4*notnanindex-1)','RBFFunction', 'thinplate');
            rbfcheck_maxdiff = rbfcheck(op);
            if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
            fi12 = rbfinterp([coordinatesFEM(:,1:2)]', op );
            
            %%%%% F22 %%%%%
            op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
                F(4*notnanindex-0)','RBFFunction', 'thinplate');
            rbfcheck_maxdiff = rbfcheck(op);
            if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
            fi22 = rbfinterp([coordinatesFEM(:,1:2)]', op );
            
            %%%%% Assemble [F11,F21,F12,F22] %%%%%
            F_rbf_thinplate = [fi11(:),fi21(:),fi12(:),fi22(:)]';  F_rbf_thinplate = F_rbf_thinplate(:);
            
            % ------ Plot ------
            obj.USubpb1{1} = U_rbf_thinplate;
            obj.FSubpb1{1} = F_rbf_thinplate;
            
            
            obj = unzipDICmesh(obj,DICmesh);
            obj = unzipDICPara(obj,DICpara);
            
        end
        
        function obj = Subproblem2(obj,ImgSeqNum)
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
            
            DICpara.DispFilterSize=0;
            DICpara.DispFilterStd=0;
            DICpara.StrainFilterSize=0;
            DICpara.StrainFilterStd=0;
            DICpara.DispSmoothness = 0;
            DICpara.StrainSmoothness = 1e-4;
            
            if DICpara.DispSmoothness>1e-6
                obj.USubpb1{1} = funSmoothDispQuadtree(obj.USubpb1{1},DICmesh,DICpara);
            end
            
            if DICpara.StrainSmoothness>1e-6
                obj.FSubpb1{1}  = funSmoothStrainQuadtree(obj.FSubpb1{1} ,DICmesh,DICpara);
            end
            
            % ====== Define penalty parameter ======
            obj.udual{1} = 0*obj.FSubpb1{1};
            obj.vdual{1} = 0*obj.USubpb1{1};
            betaList = [1e-3,1e-2,1e-1]*mean(DICpara.winstepsize).^2.*obj.mu; % Tune beta in the betaList   (taken from main ALDIC Quadtree inc.m)
            Err1 = zeros(length(betaList),1);
            Err2 = Err1;
            
            % ====== Solver using finite element method ======
            if ImgSeqNum == 2
                for tempk = 1:length(betaList)
                    obj.beta = betaList(tempk);
                    [obj.USubpb2{1}] = Subpb2Quadtree(DICmesh,obj.GaussPtOrder,obj.beta,obj.mu,obj.USubpb1{1},obj.FSubpb1{1},obj.udual{1},obj.vdual{1},obj.alpha,mean(DICpara.winstepsize),0);
                    [obj.FSubpb2{1},~,~] = funGlobalNodalStrainQuadtree(DICmesh,obj.USubpb2{1},obj.GaussPtOrder,0);
                    
                    Err1(tempk) = norm(obj.USubpb1{1}-obj.USubpb2{1},2);
                    Err2(tempk) = norm(obj.FSubpb1{1}-obj.FSubpb2{1},2);
                end
                
                Err1Norm = (Err1-mean(Err1))/std(Err1); % figure, plot(Err1Norm);
                Err2Norm = (Err2-mean(Err2))/std(Err2); % figure, plot(Err2Norm);
                ErrSum = Err1Norm+Err2Norm; % figure, plot(ErrSum); title('Tune the best \beta in the subproblem 2');
                [~,indexOfbeta] = min(ErrSum);
                
                try % Tune the best beta by a quadratic polynomial 0fitting
                    [fitobj] = fit(log10(betaList(indexOfbeta-1:1:indexOfbeta+1))',ErrSum(indexOfbeta-1:1:indexOfbeta+1),'poly2');
                    p = coeffvalues(fitobj);
                    obj.beta = 10^(-p(2)/2/p(1));
                catch, obj.beta = betaList(indexOfbeta);
                end
                
            else
                
                try obj.beta = DICpara.beta;
                catch, obj.beta = 1e-3*mean(DICpara.winstepsize).^2.*obj.mu;
                end
                
            end
            
            % Using the optimal beta to solve the ALDIC Subproblem 2 again
            if abs(obj.beta-betaList(end))>abs(eps)
                [obj.USubpb2{1}] = Subpb2Quadtree(DICmesh,obj.GaussPtOrder,obj.beta,obj.mu,obj.USubpb1{1},obj.FSubpb1{1},obj.udual{1},obj.vdual{1},obj.alpha,mean(DICpara.winstepsize),0);
                [obj.FSubpb2{1},~,~] = funGlobalNodalStrainQuadtree(DICmesh,obj.USubpb2{1},obj.GaussPtOrder,0);
            end
            
            % ------- Smooth diap & strain field --------
            if DICpara.DispSmoothness>1e-6, obj.USubpb2{1} = funSmoothDispQuadtree(obj.USubpb2{1},DICmesh,DICpara); end
            % ------- Don't smooth strain fields near the boundary --------
            for tempk=0:3, obj.FSubpb2{1}(4*DICmesh.markCoordHoleEdge-tempk) = obj.FSubpb1{1}(4*DICmesh.markCoordHoleEdge-tempk); end
            if DICpara.StrainSmoothness>1e-6, obj.FSubpb2{1} = funSmoothStrainQuadtree(0.1*obj.FSubpb2{1}+0.9*obj.FSubpb1{1},DICmesh,DICpara); end
            for tempk=0:1, obj.USubpb2{1}(2*DICmesh.markCoordHoleStrain-tempk) = obj.USubpb1{1}(2*DICmesh.markCoordHoleStrain-tempk); end
            for tempk=0:3, obj.FSubpb2{1}(4*DICmesh.markCoordHoleStrain-tempk) = obj.FSubpb1{1}(4*DICmesh.markCoordHoleStrain-tempk); end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            obj.udual{1} = obj.FSubpb2{1} - obj.FSubpb1{1};
            obj.vdual{1} = obj.USubpb2{1} - obj.USubpb1{1};
            
            obj = unzipDICmesh(obj,DICmesh);
            obj = unzipDICPara(obj,DICpara);
            
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
            
            % ==================== ADMM AL Loop ==========================
            obj.ALSolveStep = 1;
            tol2 = 1e-2;
            UpdateY = 1e4;
            HPar = cell(21,1);
            for tempj = 1:21, HPar{tempj} = obj.HtempPar(:,tempj); end
            
            % preform ADMM
            while (obj.ALSolveStep < 3)
                
                % update variables
                TempUSubpb2 = obj.USubpb2{obj.ALSolveStep};
                TempFSubpb2 = obj.FSubpb2{obj.ALSolveStep};
                Tempudual = obj.udual{obj.ALSolveStep};
                Tempvdual = obj.vdual{obj.ALSolveStep};
                
                obj.ALSolveStep = obj.ALSolveStep + 1;  % Update using the last step
                
                % constant window size? or adjusted subset size?
                if DICpara.winsize_opt == 'constant'
                    winsize_List = DICpara.winsize*ones(size(DICmesh.coordinatesFEM,1),2);
                    DICpara.winsize_List = winsize_List;
                elseif DICpara.winsize_opt == 'variable'
                    %%%%% If wen want to adjust the DIC subset size %%%%%
                    Ftemp1 = obj.FSubpb2{1}(1:2:end); Ftemp2 = obj.FSubpb2{1}(2:2:end);
                    [DFtemp1,~,~] = funGlobalNodalStrainQuadtree(DICmesh,Ftemp1,DICpara.GaussPtOrder,0);
                    [DFtemp2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,Ftemp2,DICpara.GaussPtOrder,0);
                    
                    winsize_x_ub1 = abs(2*obj.FSubpb2{1}(1:4:end)./DFtemp1(1:4:end));
                    winsize_x_ub2 = abs(2*obj.FSubpb2{1}(3:4:end)./DFtemp1(3:4:end));
                    winsize_y_ub1 = abs(2*obj.FSubpb2{1}(1:4:end)./DFtemp1(2:4:end));
                    winsize_y_ub2 = abs(2*obj.FSubpb2{1}(3:4:end)./DFtemp1(4:4:end));
                    
                    winsize_x_ub3 = abs(2*obj.FSubpb2{1}(2:4:end)./DFtemp2(1:4:end));
                    winsize_x_ub4 = abs(2*obj.FSubpb2{1}(4:4:end)./DFtemp2(3:4:end));
                    winsize_y_ub3 = abs(2*obj.FSubpb2{1}(2:4:end)./DFtemp2(2:4:end));
                    winsize_y_ub4 = abs(2*obj.FSubpb2{1}(4:4:end)./DFtemp2(4:4:end));
                    
                    winsize_x_ub = round(min([winsize_x_ub1,winsize_x_ub2,winsize_x_ub3,winsize_x_ub4,DICpara.winsize*ones(length(winsize_x_ub1),1)],[],2));
                    winsize_x_List = max([winsize_x_ub, 10*ones(length(winsize_x_ub1),1)],[],2);
                    winsize_y_ub = round(min([winsize_y_ub1,winsize_y_ub2,winsize_y_ub3,winsize_y_ub4,DICpara.winsize*ones(length(winsize_y_ub1),1)],[],2));
                    winsize_y_List = max([winsize_y_ub, 10*ones(length(winsize_y_ub1),1)],[],2);
                    winsize_List = 2*ceil([winsize_x_List,winsize_y_List]/2);
                    DICpara.winsize_List = winsize_List;
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%%% Subproblem 1 %%%%%%%%%%%%%%%%%%%%%%%%%
                [TempUSubpb1,~,~,~,~] = Subpb1Quadtree(...
                    TempUSubpb2,TempFSubpb2,Tempudual,Tempvdual,DICmesh.coordinatesFEM,...
                    obj.Df,ImgRef,ImgDef,obj.mu,obj.beta,HPar,obj.ALSolveStep,DICpara,'GaussNewton',DICpara.tol);
                
                TempFSubpb1 = TempFSubpb2;
                
                for tempk=0:1, TempUSubpb1(2*DICmesh.markCoordHoleStrain-tempk) = TempUSubpb2(2*DICmesh.markCoordHoleStrain-tempk); end
                
                
                %Save the results of U and F sub prob 1 to array
                obj.USubpb1{obj.ALSolveStep} = TempUSubpb1;
                obj.FSubpb1{obj.ALSolveStep} = TempFSubpb1;
                
                TempUSubpb1 = funSmoothDisp(TempUSubpb1,DICmesh,DICpara);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % ============== Subproblem 2 ==============
                
                [TempUSubpb2] = Subpb2Quadtree(DICmesh,obj.GaussPtOrder,obj.beta,obj.mu,TempUSubpb1,TempFSubpb1,Tempudual,Tempvdual,obj.alpha,mean(DICpara.winstepsize),0);
                [TempFSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,TempUSubpb2,obj.GaussPtOrder,0);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % ------- Smooth strain field --------
                if DICpara.DispSmoothness>1e-6, TempUSubpb2 = funSmoothDispQuadtree(TempUSubpb2,DICmesh,DICpara); end
                % ------- Don't change strain fields near the boundary --------
                for tempk=0:3, TempFSubpb2(4*DICmesh.markCoordHoleEdge-tempk) = TempFSubpb1(4*DICmesh.markCoordHoleEdge-tempk); end
                if DICpara.StrainSmoothness>1e-6, TempFSubpb2 = funSmoothStrainQuadtree(0.1*TempFSubpb2+0.9*TempFSubpb1,DICmesh,DICpara); end
                for tempk=0:1, TempUSubpb2(2*DICmesh.markCoordHoleStrain-tempk) = TempUSubpb1(2*DICmesh.markCoordHoleStrain-tempk); end
                for tempk=0:3, TempFSubpb2(4*DICmesh.markCoordHoleStrain-tempk) = TempFSubpb1(4*DICmesh.markCoordHoleStrain-tempk); end
                
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
                
                % Update dual variables------------------------------
                Tempudual = TempFSubpb2 - TempFSubpb1;
                Tempvdual = TempUSubpb2 - TempUSubpb1;
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
            
            obj = unzipDICmesh(obj,DICmesh);
            obj = unzipDICPara(obj,DICpara);
            
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
            DICpara.ImgDefMask = obj.DICparaImgDefMask;
            DICpara.um2px = obj.DICparaum2px;
            DICpara.DoYouWantToSmoothOnceMore = obj.DICparaDoYouWantToSmoothOnceMore;
            DICpara.MethodToComputeStrain = obj.DICparaMethodToComputeStrain;
            DICpara.StrainType = obj.DICparaStrainType;
            DICpara.Image2PlotResults = obj.DICparaImage2PlotResults;
            DICpara.MethodToSaveFig = obj.DICparaMethodToSaveFig;
            DICpara.OrigDICImgTransparency = obj.DICparaOrigDICImgTransparency;
            DICpara.winsizeMin = obj.DICparawinsizeMin ;
            DICpara.winsize_List =obj.DICparawinsize_List;
            DICpara.winsize_opt =obj.DICparawinsize_opt;
            DICpara.DispSmoothness = obj.DICparaDispSmoothness;
            DICpara.StrainSmoothness = obj.DICparaStrainSmoothness;
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
            obj.DICparaImgDefMask = DICpara.ImgDefMask;
            obj.DICparaum2px = DICpara.um2px;
            obj.DICparaDoYouWantToSmoothOnceMore = DICpara.DoYouWantToSmoothOnceMore;
            obj.DICparaMethodToComputeStrain = DICpara.MethodToComputeStrain;
            obj.DICparaStrainType = DICpara.StrainType;
            obj.DICparaImage2PlotResults = DICpara.Image2PlotResults;
            obj.DICparaMethodToSaveFig = DICpara.MethodToSaveFig;
            obj.DICparaOrigDICImgTransparency = DICpara.OrigDICImgTransparency;
            obj.DICparawinsizeMin = DICpara.winsizeMin ;
            obj.DICparawinsize_List = DICpara.winsize_List;
            obj.DICparawinsize_opt = DICpara.winsize_opt;
            obj.DICparaDispSmoothness = DICpara.DispSmoothness;
            obj.DICparaStrainSmoothness = DICpara.StrainSmoothness;
        
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
            try DICmesh.markCoordHoleEdge = obj.DICmesh_markCoordHoleEdge;  catch;  end
            try DICmesh.markCoordHoleStrain = obj.DICmesh_markCoordHoleStrain;  catch;  end
            try DICmesh.irregular =obj.DICmesh_irregular;  catch;  end
            
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
            try obj.DICmesh_markCoordHoleEdge = DICmesh.markCoordHoleEdge;  catch;  end
            try obj.DICmesh_markCoordHoleStrain = DICmesh.markCoordHoleStrain;  catch;  end
            try obj.DICmesh_irregular =DICmesh.irregular;  catch;  end
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