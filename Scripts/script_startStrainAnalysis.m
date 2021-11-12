
for i = 1:size(Data,1)
    ImgSeqNum = Data{i,1};
    obj = app.SavedAnalysis{Data{i,2}}; % get object
    
    % run the strain data for the object and the deformation fields
    % selected
    
    % get the selected Strain computation method
    switch app.StraincomputationmethodDropDown.Value
        case 'Finite Difference'
            StrainOptions.Method = 1;
        case 'Finite Element'
            StrainOptions.Method = 3;    
        case 'Plane Fitting'
            StrainOptions.Method = 2;
            StrainOptions.PWS = app.WindowSizeEditField_StrainPlane.Value;
        case 'Direct output'
            StrainOptions.Method = 0;
        otherwise
            StrainOptions.Method = 0;
    end
    
    
    % get the selected Strain Type 
    switch app.StraintypeDropDown.Value
        case 'Infinitesimal'
            StrainOptions.StrainType = 0;
        case 'Eulerian'
            StrainOptions.StrainType = 1;    
        case 'Green-Lagrangian'
            StrainOptions.StrainType = 2;
        otherwise
            StrainOptions.StrainType = 0;
    end
    
    
    need2updateConversion = false;
    switch StrainOptions.Method
        case 0
            if isempty(obj.ResultStrainDirect{ImgSeqNum})
                obj = computeStrains(obj,StrainOptions,ImgSeqNum);
                need2updateConversion = true;
            end
        case 1
            if isempty(obj.ResultStrainFD{ImgSeqNum})
                obj = computeStrains(obj,StrainOptions,ImgSeqNum);
                need2updateConversion = true;
            end
        case 2
            if isempty(obj.ResultStrainPlane{ImgSeqNum})
                obj = computeStrains(obj,StrainOptions,ImgSeqNum);
                need2updateConversion = true;
            end
        case 3
            if isempty(obj.ResultStrainFE{ImgSeqNum})
                obj = computeStrains(obj,StrainOptions,ImgSeqNum);
                need2updateConversion = true;
            end
    end
  

    need2updateMaps = false;
    if or(isempty(obj.ResultStrainWorld{ImgSeqNum,StrainOptions.Method+1,StrainOptions.StrainType+1}),need2updateConversion)
        obj = convertStrainType(obj,StrainOptions,ImgSeqNum);
        need2updateMaps = true;
    end
    
    if need2updateMaps
        [StrainMap,obj] = getStrainMap(obj,StrainOptions,ImgSeqNum);
    else
        StrainMap = obj.ResultStrainMap(ImgSeqNum,:);
    end

    switch app.CurrentStrainPlotType
        case 'exx'
            Smap = 3;
        case 'exy'
            Smap = 4;
        case 'eyy'
            Smap = 5;
        case 'max shear'
            Smap = 6;
        case 'max principal'
            Smap = 7;
        case 'min principal'
            Smap = 8;
        case 'vonMises'
            Smap = 9;
    end
    
    app.SavedAnalysis{Data{i,2}} = obj; % return the altered object
    
    e = StrainMap{Smap};
    % initialize the output strain
    x = [min(obj.DICmesh_y0,[],'all'),max(obj.DICmesh_y0,[],'all')];  % not a mistake for x = y0
    y = [min(obj.DICmesh_x0,[],'all'),max(obj.DICmesh_x0,[],'all')];
    % remake ROI array based on the DICmesh 
    E = zeros(size(obj.ROIMasks{1}));
    E(y(1):y(2),x(1):x(2))=1;
    % scale and rotate the output to match the ROI
    ee = imresize(flipud(imrotate(e,90)),[y(2)-y(1)+1,x(2)-x(1)+1]);
    % place that data in the ROI
    E(y(1):y(2),x(1):x(2)) = ee;
    
    scaleBar = { 10, 200, 'px', 1.6, [30,30]} ;
    colorScale = {min(E,[],'all'),max(E,[],'all'),' '};
    [ImgoutE] = Imageoverlay_withcolorbar3(obj.Images{i}, E, double(obj.ROIMasks{i}), 0.5, ' ', scaleBar , colorScale);
    func_PlaceFigure(app.StrainMap,ImgoutE);
    
    % tell the user what data is being plotted
    CurStrainPlot = ['Plot of data set ' Data{i,3} ' from Analysis ' Data{i,4}];
    app.StrainCurrentPlotLabel.Text = CurStrainPlot;
    
    pause(0.5)
    
end