function app = StartALDIC(app)
    
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
    obj.DICparawinstepsize = winstep;
    
    % ALDIC options
    obj.DICparatol = app.toleranceEditField_ALDIC.Value;
    obj.DICparamu = app.muEditField_ALDIC.Value;
    
    
    SP1methodString = app.Subproblem1MethodButtonGroup_ALDIC.SelectedObject.Text;
    switch SP1methodString
        case 'Levenberg Marquardt'
            SP1method = 0;
        case 'Gauss Newton'
            SP1method = 1;
            
        otherwise
            SP1method = 0;
    end
      
    SP2methodString = app.Subproblem1MethodButtonGroup_ALDIC.SelectedObject.Text;
    switch SP2methodString
        case 'Finite Difference'
            SP2method = 0;
        case 'Finite Element Method'
            SP2method = 1;
            
        otherwise
            SP2method = 0;
    end
    obj.DICparaSubpb2FDOrFEM = SP2method;
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
    obj.DICparaImgSize = size(IMarray{1});
    
    %step2 normalize images
    [x,y] = obj.convertROItogridpoints(ROIarray{1});
    obj.DICparagridx = x;
    obj.DICparagridy = y;
    DICpara = obj.zipDICPara();
    [ImgNormalized,DICpara.gridxyROIRange] = funNormalizeImg(IMarray,DICpara.gridxyROIRange);
    obj = obj.unzipDICPara(DICpara);

    %step3 compute image gradient of the reference image
    ImgRef = ImgNormalized{1};
    obj.Df = funImgGradient(ImgRef,ImgRef); 

    %step4 Initialize variable storage 
    obj = obj.InitResultVars(length(ImgNormalized));

    %step5 solve each frame in the image sequence
    
    for i = 2:length(IMarray)
       
        % Def image
        ImgDef = ImgNormalized{i};

        obj = obj.runALDIC(ImgRef,ImgDef,i);
        
        u = obj.u_f;
        v = obj.v_f;
        x0=obj.DICmesh_x0;
        y0=obj.DICmesh_y0;
        
        U = double(ROIarray{i});
        uu = imresize(flipud(imrotate(u,90)),[y(2)-y(1)+1,x(2)-x(1)+1]);
        U(y(1):y(2),x(1):x(2)) = uu;
        
        V = double(ROIarray{i});
        vv = imresize(flipud(imrotate(v,90)),[y(2)-y(1)+1,x(2)-x(1)+1]);
        V(y(1):y(2),x(1):x(2)) = vv;
        
        scaleBar = { 10, 200, 'px', 1.6, [30,30]} ;
        colorScale = {min(U,[],'all'),max(U,[],'all'),' px '};
        [ImgoutU] = Imageoverlay_withcolorbar3(IMarray{i}, U, double(ROIarray{i}), 0.5, ' ', scaleBar , colorScale);
        colorScale = {min(V,[],'all'),max(V,[],'all'),' px '};
        [ImgoutV] = Imageoverlay_withcolorbar3(IMarray{i}, V, double(ROIarray{i}), 0.5, ' ', scaleBar , colorScale);
        
        func_PlaceFigure(app.A_Udisp_ALDIC,ImgoutU);
        func_PlaceFigure(app.A_Vdisp_ALDIC,ImgoutV);
        
    end
    
    %step6 check convergence
    %obj = obj.checkConvergence(obj);
    
end