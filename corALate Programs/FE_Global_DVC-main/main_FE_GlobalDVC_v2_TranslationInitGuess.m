% ---------------------------------------------------
% Finite-Element-Based Global Digital Volume Correlation (FE-Global-DVC)
% Author: Jin Yang, Postdoc UW-Madison; PhD '19 Caltech;
% Contact: jyang526@wisc.edu;  aldicdvc@gmail.com
% 2019.07, 2020.11
% ---------------------------------------------------

%% Section 0
% ====== Prepare volumetric image data files ======
% Please go to subfolder code: './DVC_images/GenerateVolMatfile.m' to
% transform your volumetric image stacks to a Matlab matfile for the DVC code.


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start each frame in an image sequence
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Img{1} = ImgNormalized{1}; % The first frame is reference 
for ImgSeqNum = 2 % : length(ImgNormalized)
     
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]); 
    
    %% Section 3: Find initial guess
    fprintf('\n'); fprintf('------------ Section 3 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to find or update initial guess for ALDVC
    % The key idea is to either to use FFT peak fitting, or use last frame
    % results for the next new frame;
    % Particularly in incremental mode, the reference image can also be updated.
    % fNormalized = ImgNormalized{ImgSeqNum-mod(ImgSeqNum-1,ImgSeqIncUnit)};
    Img{2} = ImgNormalized{ImgSeqNum}; NewFFTSearchCheck = 0; DVCpara.NewFFTSearch = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ImgSeqNum==2 || DVCpara.NewFFTSearch==1
        
        [xyz0_x,xyz0_y,xyz0_z] = ndgrid( [ max([4,DVCpara.gridRange.gridxRange(1)]) :  DVCpara.winstepsize(1)  :  min([DVCpara.ImgSize(1)-3,DVCpara.gridRange.gridxRange(2)])  ] , ...
            [ max([4,DVCpara.gridRange.gridyRange(1)]) :  DVCpara.winstepsize(2)  :  min([DVCpara.ImgSize(2)-3,DVCpara.gridRange.gridyRange(2)])  ] , ...
            [ max([4,DVCpara.gridRange.gridzRange(1)]) :  DVCpara.winstepsize(3)  :  min([DVCpara.ImgSize(3)-3,DVCpara.gridRange.gridzRange(2)])  ] );
        
        xyz0.x = xyz0_x; xyz0.y = xyz0_y; xyz0.z = xyz0_z;
        
        DVCpara.SizeOfFFTSearchRegion = 0; % Default
        DVCpara.qDICOrNot=0; DVCpara.Thr0=0; % Default
        
        % ====== FEM mesh set up ======
        [DVCmesh] = MeshSetUp3(xyz0,DVCpara);
        
        uvw.u = 3.568*ones(size(xyz0.x)); 
        uvw.v = -2.751*ones(size(xyz0.x)); 
        uvw.w = 0.1286*ones(size(xyz0.x)); 
        
        % ====== Assign initial values ======
        U0 = Init3(uvw,DVCmesh.xyz0); 
        %%%%% Plotdisp_show3(U0,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM); 
          
        
        % ====== Deal with incremental mode ======
        Img1NewIndex = ImgSeqNum-mod(ImgSeqNum-2,DVCpara.ImgSeqIncUnit)-1;
        if DVCpara.ImgSeqIncUnit==1, Img1NewIndex = Img1NewIndex-1; end
        ResultFEMesh{1+floor(Img1NewIndex/DVCpara.ImgSeqIncUnit)} = ... % To save first mesh info
            struct( 'coordinatesFEM',DVCmesh.coordinatesFEM,'elementsFEM',DVCmesh.elementsFEM, ...
            'winsize',DVCpara.winsize,'winstepsize',DVCpara.winstepsize,'gridxyROIRange',DVCpara.gridRange );
       
    elseif mod(ImgSeqNum-2,DVCpara.ImgSeqIncUnit)==0 % TO update ref image in incremental mode
        Img1NewIndex = ImgSeqNum-mod(ImgSeqNum-2,DVCpara.ImgSeqIncUnit)-1;
        if DVCpara.ImgSeqIncUnit==1,  Img1NewIndex = Img1NewIndex-1; end
        Img{1} = ImgNormalized{Img1NewIndex}; % Update reference
        [DVCpara,DVCmesh] = ReadImageRefUpdate(file_name,ImgSeqNum,ResultDisp{ImgSeqNum-2}.U,DVCpara,DVCmesh); % Update reference image if needed;
        U0 = zeros(3*size(DVCmesh.coordinatesFEM,1),1); % PlotuvInit;
        ResultFEMesh{1+floor(Img1NewIndex/DVCpara.ImgSeqIncUnit)} = ... % To save first mesh info
            struct( 'coordinatesFEM',DVCmesh.coordinatesFEM,'elementsFEM',DVCmesh.elementsFEM, ...
            'winsize',DVCpara.winsize,'winstepsize',DVCpara.winstepsize,'gridxyROIRange',DVCpara.gridRange );
    else
        U0 = ResultDisp{ImgSeqNum-2}.U;
    end
    % ====== Spline interpolation images ======
    disp('--- Start to compute image gradients ---');
    Df = funImgGradient3(Img{1},'stencil7'); 
    Df.imgSize = size(Img{1}); % Compute image grayscale value gradients
    disp('--- Computing image gradients done ---');
    
    fprintf('------------ Section 3 Done ------------ \n \n')
   
    
    %% Section 4 Global DVC
    fprintf('------------ Section 4 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Finite element based global DVC iterations
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    alphaList = DVCpara.alpha; % Set regularization coefficient, alpha, as 100 as an example
    
    % ====== Tune regularization coefficient ======
    % If you don't know the best alpha (coefficient), please run the following
    % codes to tune the best value of the coefficient of the regularizer |grad u|^2):
    % 
    % %%%%% Uncomment the following line to tune the best value of alpha %%%%%%
    alphaList = [ 100 ]*mean(DVCpara.winstepsize);
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

% ImgSeqNum = 2; U = ResultDisp{ImgSeqNum-1}.U;
% Plotdisp_show3(U,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM);
% ====== Compute f(X)-g(x+u) ======
% PlotImgDiff(xyz0.x,xyz0.y, ImgNormalized{1},ImgNormalized{2}); % Img grayscale value residual
    


%% Section 5
fprintf('------------ Section 5 Start ------------ \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to compute strain and plot figures
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Smooth displacements ------
DVCpara.DoYouWantToSmoothOnceMore = funParaInput('SmoothDispOrNot');
% ------ Choose strain computation method ------
DVCpara.MethodToComputeStrain = funParaInput('StrainMethodOp'); 
% ------ Choose strain type (infinitesimal, Eulerian, Green-Lagrangian) ------
DVCpara.StrainType = funParaInput('StrainType');
% ------ Plot displacement & strain components individually or all together ------
DVCpara.PlotComponentEachOrAll = funParaInput('PlotComponentEachOrAll');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Start plotting part -----
for ImgSeqNum = 2 %:length(ImgNormalized)
    
    % disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DVCpara.ImgSeqIncUnit)-1;
    if DVCpara.ImgSeqIncUnit > 1
        FEMeshIndLast = floor(fNormalizedNewIndex/DVCpara.ImgSeqIncUnit);
    elseif DVCpara.ImgSeqIncUnit == 1
        FEMeshIndLast = floor(fNormalizedNewIndex/DVCpara.ImgSeqIncUnit)-1;
    end
    FEMeshInd = FEMeshIndLast + 1;
    
    if FEMeshInd == 1
        USubpb2 = ResultDisp{ImgSeqNum-1}.U; %+ ResultDisp{10}.U + ResultDisp{20}.U;
        coordinatesFEM = ResultFEMesh{1}.coordinatesFEM; 
        elementsFEM = ResultFEMesh{1}.elementsFEM;
        if (ImgSeqNum-1 == 1) || (DVCpara.ImgSeqIncROIUpdateOrNot==1), UFEMesh = 0*USubpb2; end
    else
        USubpb2 = ResultDisp{ImgSeqNum-1}.U;
        if mod(ImgSeqNum-2,DVCpara.ImgSeqIncUnit) == 0
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %USubpb2 = ResultDisp{ImgSeqNum-1}.U;
    FSubpb2 = ResultDefGrad{ImgSeqNum-1}.F;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xList = min(coordinatesFEM(:,1)):DVCpara.winstepsize(1):max(coordinatesFEM(:,1)); M = length(xList);
    yList = min(coordinatesFEM(:,2)):DVCpara.winstepsize(2):max(coordinatesFEM(:,2)); N = length(yList);
    zList = min(coordinatesFEM(:,3)):DVCpara.winstepsize(3):max(coordinatesFEM(:,3)); L = length(zList);
    [xGrid,yGrid,zGrid] = ndgrid(xList,yList,zList);
    xGrid = xGrid-reshape(UFEMesh(1:3:end),size(xGrid));
    yGrid = yGrid-reshape(UFEMesh(2:3:end),size(yGrid));
    zGrid = zGrid-reshape(UFEMesh(3:3:end),size(zGrid));
    xyz0.x = xGrid; xyz0.y = yGrid; xyz0.z = zGrid;
 
    if size(USubpb2,1)==1
        ULocal = full(USubpb2_New.USubpb2); FLocal = full(FSubpb2.FSubpb2); 
    else
        ULocal = full(USubpb2); FLocal = full(FSubpb2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ Smooth displacements ------
    %Plotdisp_show3(full(ULocal),coordinatesFEM,elementsFEM);
    % prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
    % DoYouWantToSmoothOnceMore = input(prompt); DispFilterSize=0; DispFilterStd=1;
    SmoothTimes=0;
    try
        while DVCpara.DoYouWantToSmoothOnceMore==0 && SmoothTimes<3
            ULocal = funSmoothDisp3(ULocal,DVCmesh,DVCpara);
            % close all; Plotdisp_show3(full(ULocal),coordinatesFEM,elementsFEM); % Plotuv(ULocal,x0,y0); 
            SmoothTimes = SmoothTimes + 1; %DoYouWantToSmoothOnceMore = input(prompt);
        end
    catch
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ----- Compute strain field ------
    ComputeStrain3;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ----- Plot results -----
    close all;
    
    % ------ Plot disp ------
    Plotdisp03(ULocal,DVCmesh.coordinatesFEM,DVCmesh.elementsFEM,DVCpara.PlotComponentEachOrAll);
     
    % ------ Plot strain ------
    % %%%%% Or use this line %%%%%
    Plotstrain_show3(FLocal,coordinatesFEM,elementsFEM);
    
    % %%%%% Or use this line %%%%%
    % Plotstrain03(full(FStraintemp),xyz0.x(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3)), ...
    %     xyz0.y(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3)), ...
    %     xyz0.z(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3)),size(Img{1}),DVCpara.PlotComponentEachOrAll);
    
    % ------ Store strain data ------
    ResultStrain{ImgSeqNum-1}.Strain = FStraintemp;
    tempx = xyz0.x(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3));
    tempy = xyz0.y(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3));
    tempz = xyz0.z(1+Rad(1):M-Rad(1),1+Rad(2):N-Rad(2),1+Rad(3):L-Rad(3));
    ResultStrain{ImgSeqNum-1}.coordinatesFEMStrain = [tempx(:), tempy(:), tempz(:)]; 
    
    % % ------- Add filter and plot strain field -------
    % Plotstrain_Fij;
    % %caxis auto; load('colormap_RdYlBu.mat'); colormap(cMap)
    % Plotstrain03(FStraintemp,xyz0.x(1+Rad:M-Rad,1+Rad:N-Rad,1+Rad:L-Rad),xyz0.y(1+Rad:M-Rad,1+Rad:N-Rad,1+Rad:L-Rad), ...
    %    xyz0.z(1+Rad:M-Rad,1+Rad:N-Rad,1+Rad:L-Rad),size(Img{1}),DVCpara.PlotComponentEachOrAll);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ Save figures ------
    % Write down your own codes to save figures! E.g.: % print(['fig_dispu'],'-dpdf');
    if strcmp(DVCpara.PlotComponentEachOrAll,'All')==1
        figure(1); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_disp.fig']);
        figure(2); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_strain.fig']);
    elseif strcmp(DVCpara.PlotComponentEachOrAll,'Individual')==1
        figure(1); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_dispx.fig']);
        figure(2); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_dispy.fig']);
        figure(3); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_dispz.fig']);
        figure(4); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_strainexx.fig']);
        figure(5); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_straineyy.fig']);
        figure(6); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_strainezz.fig']);
        figure(7); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_strainexy.fig']);
        figure(8); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_strainexz.fig']);
        figure(9); saveas(gcf,['fig_ImgSeqNum_',num2str(ImgSeqNum),'_straineyz.fig']);
    else
        disp('=== Wrong input in DVCpara.PlotComponentEachOrAll! ===')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
fprintf('------------ Section 5 Done ------------ \n \n')


%% Save data again including the computed strain 
results_name = ['results_FE_globalDVC_st',num2str(DVCpara.winstepsize(1)),'_alpha',num2str(DVCpara.alpha),'.mat'];
save(results_name, 'file_name','DVCpara','DVCmesh','ResultDisp','ResultDefGrad','ResultStrain','ResultFEMesh', ...
     'ResultAlpha','ResultNormOfW','ResultTimeICGN');

  
%% %%%%%%%%%%%%% Extensions for body and slice plottings %%%%%%%%%%%%%%%%
disp('Extensions for body and slice plottings'); pause;  
plotExt_bodyslice; % Feel free to modify this file (./PlotFiles/plotExt_bodyslice.m) on your purpose.


%% %%%%%%%%%%%%% Extensions %%%%%%%%%%%%%%%%
 
% ====== Compute f(X)-g(x+u) ======
% Img grayscale value residual
% tempx = DVCmesh.xyz0.x; tempy = DVCmesh.xyz0.y; tempz = DVCmesh.xyz0.z;

PlotImgDiff3(DVCmesh.coordinatesFEM(:,1), DVCmesh.coordinatesFEM(:,2), ...
    DVCmesh.coordinatesFEM(:,3), -ResultDisp{1,1}.U(1:3:end), ...
    -ResultDisp{1,1}.U(2:3:end), -ResultDisp{1,1}.U(3:3:end), ...
    Img{1}, Img{2});
    

    
    
    
    
    
    
    
    