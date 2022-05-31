function app = StartSTAQ(app)
    
reverse = false;

    %% get the current analysis object
    pos = app.CurrentAnalysisNum;
    obj = app.SavedAnalysis{pos,1};
    
    %% push the user defined settings into the object
    % Type of search
    IGSmethodString = app.RadioGroup_InitialguessMethod_ALDIC.SelectedObject.Text;
    switch IGSmethodString
        case 'Multigrid search'
            method = 0;
        case 'Whole field search'
            method = 1;
        otherwise
            method = 0;
    end
    obj.DICparaInitFFTSearchMethod = method; % always MG search right now
    
    % win search size
    winsize = app.windowsizeEditField_ALDIC.Value;
    obj.DICparaSizeOfFFTSearchRegion = [0, 0];
    obj.DICparawinsize = winsize;
    
    % win step size
    winstep = app.stepsizeEditField_ALDIC.Value;
    obj.DICparawinstepsize = [winstep];
    
    % ALDIC options
    obj.DICparatol = app.toleranceEditField_ALDIC.Value;
    obj.DICparamu = app.muEditField_ALDIC.Value;
    obj.DICparaClusterNo = app.numberofparallelpoolsSpinner_ALDIC.Value;
    
    if app.UsemostrecentframeasguessCheckBox_ALDIC.Value
        Value = 0;
    else
        Value = 1;
    end
    
    obj.DICparaNewFFTSearch = Value;
    
    %% Begin ALDIC code
    
    % Gather the ROI for each image
    ROIarray = app.AllROIMasks;
    % Gather the images to preform the initial guess
    
    %step1 load images between the first selected reference image and the next
    IMarray = app.AllImages;
    N = length(IMarray);
    refIndicies = app.RefSelction;
    
    % save all ROI data and Image Data into the object
    obj.Images = IMarray;
    obj.ROIMasks = ROIarray;
    
    obj.DICparaImgSize = size(IMarray{1});
    
    %step2 normalize images
    [y,x] = obj.convertROItogridpoints(ROIarray{1});
    obj.DICparagridx = x;
    obj.DICparagridy = y;
    DICpara = obj.zipDICPara();
    [ImgNormalized,DICpara.gridxyROIRange] = funNormalizeImg(IMarray,DICpara.gridxyROIRange);
    
    obj = obj.unzipDICPara(DICpara);
    
    %step3 compute image gradient of the reference image and multiply by mask
    ReffNormalizedMask = double( ROIarray{1} ) ; % Load the mask file of previous frame
    ImgRef = ImgNormalized{1} .* ReffNormalizedMask; % Load previous frame
    obj.Df = funImgGradient(ImgRef,ImgRef,ReffNormalizedMask); % Finite difference to compute image grayscale gradients;
    % store the reference image mask
     obj.DICparaImgRefMask = double( ROIarray{1}) ;

    %step4 Initialize variable storage
    obj = obj.InitResultVars(length(ImgNormalized));

    %step5 solve each frame in the image sequence
    for i = 2: N
        
        % load deformaed image and multiply by mask
        obj.DICparaImgDefMask = double(ROIarray{i}); % Load the mask file of current frame
        ImgDef = ImgNormalized{i} .* obj.DICparaImgDefMask ; % Load current deformed image frame
        
%         figure, subplot(2,2,1); imshow(ImgRef'); title('Ref img f')
%         subplot(2,2,2); imshow(ImgDef'); title('Def img g')
%         subplot(2,2,3); imshow(obj.DICparaImgRefMask'); title('f mask')
%         subplot(2,2,4); imshow(obj.DICparaImgDefMask'); title('g mask')
        
        %step5.2 Compute an initial guess
        obj = InitialGuess(obj,ImgRef,ImgDef,i);

        %step5.3 solve Subproblem 1
        obj = Subproblem1(obj,ImgRef,ImgDef);

        %step5.4 solve Subproblem 2
        obj = Subproblem2(obj,i);

        %step5.5 apply ADMM iterations
        obj = ADMMiteration(obj,ImgRef,ImgDef,i);
        
        u = obj.Ux;
        v = obj.Vy;
        
        x = [min(obj.DICmesh_y0,[],'all'),max(obj.DICmesh_y0,[],'all')]; % not a mistake for x = y0
        y = [min(obj.DICmesh_x0,[],'all'),max(obj.DICmesh_x0,[],'all')];
        
        % remake ROI array based on the DICmesh 
        U = zeros(size(ROIarray{i}));
        U(y(1):y(2),x(1):x(2))=1;
        obj.ROIMasks{i} = U;
        ROIarray{i} = U;
        uu = imresize(u,[y(2)-y(1)+1,x(2)-x(1)+1]);
        U(y(1):y(2),x(1):x(2)) = uu;
        
        V = ROIarray{i};
        vv = imresize(v,[y(2)-y(1)+1,x(2)-x(1)+1]);
        V(y(1):y(2),x(1):x(2)) = vv;
        
        scaleBar = { 10, 200, 'px', 1.6, [30,30]} ;
        colorScale = {min(U,[],'all'),max(U,[],'all'),' px '};
        [ImgoutU] = Imageoverlay_withcolorbar_STAQ(IMarray{i}, U', double(ROIarray{i}), 0.5, ' ', scaleBar , colorScale);
        colorScale = {min(V,[],'all'),max(V,[],'all'),' px '};
        [ImgoutV] = Imageoverlay_withcolorbar_STAQ(IMarray{i}, V', double(ROIarray{i}), 0.5, ' ', scaleBar , colorScale);
        
        func_PlaceFigure(app.A_Udisp_ALDIC,ImgoutU);
        func_PlaceFigure(app.A_Vdisp_ALDIC,ImgoutV);
    end
     
    %step6 check convergence
    %obj = obj.checkConvergence(obj);
     app.SavedAnalysis{pos,1} = obj;
     script_UpdateListOfObjects
end