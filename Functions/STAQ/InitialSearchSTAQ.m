function app = InitialSearchSTAQ(app)
    
    % get the current analysis object
    pos = app.CurrentAnalysisNum;
    obj = app.SavedAnalysis{pos,1};
    
    % push the initial search settings into the object
    % Type of search
    methodString = app.RadioGroup_InitialguessMethod_ALDIC.SelectedObject.Text;
    switch methodString
        case 'Multigrid search'
            method = 0;
        case 'Whole field search'
            method = 1;
            
        otherwise
            method = 0;
    end
    obj.DICparaInitFFTSearchMethod = method; % always MG search right now
    
    % win search size
    winsize = app.searchregionsizeEditField_ALDIC.Value;
    obj.DICparaSizeOfFFTSearchRegion = [0, 0];
    obj.DICparawinsize = winsize;
    
    % win step size
    winstep = app.searchstepsizeEditField_ALDIC.Value;
    obj.DICparawinstepsize = winstep;
    
    % Gather the ROI for each image
    ROIarray = app.AllROIMasks;
    
    % Gather the images to preform the initial guess
    IMarray = app.AllImages;
    obj.DICparaImgSize = size(IMarray{1});
    
    ImgRef = IMarray{1};
    func_PlaceFigure(app.IDG_ReferenceImage_ALDIC,ImgRef);
    
    muv = zeros(1,6);
    
    for i = 2:min([6,length(IMarray)])
        % ROI grid points
        [x,y] = obj.convertROItogridpoints(ROIarray{i});
        obj.DICparagridx = x;
        obj.DICparagridy = y;
        
        % Def image
        ImgDef = IMarray{i};
        
        obj = obj.InitialGuess(ImgRef,ImgDef,i);
        u = obj.u_f;
        v = obj.v_f;
        x0=obj.DICmesh_x0;
        y0=obj.DICmesh_y0;
        
        UV = double(ROIarray{i});
        uv = imresize(flipud(imrotate(sqrt(u.^2 + v.^2),90)),[y(2)-y(1)+1,x(2)-x(1)+1]);
%         uv = imresize(imrotate(sqrt(u.^2 + v.^2),90),[y(2)-y(1)+1,x(2)-x(1)+1]);
        UV(y(1):y(2),x(1):x(2)) = uv;
        
        muv(i) = max(uv,[],'all');
        
        app.maxdisplacementEditField_ALDIC.Value = max(muv);
        
        %% setup image to display
        scaleBar = { 10, 200, 'px', 1.6, [30,30]} ;
        colorScale = {0,muv(2),' px '};
        [Imgout] = Imageoverlay_withcolorbar3(ImgDef, UV, double(ROIarray{i}), 0.5, ' ', scaleBar , colorScale);
        func_PlaceFigure(app.IDG_DeformedImage_ALDIC,Imgout);
    end
    
end