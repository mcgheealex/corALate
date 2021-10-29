function [file_name,Img,DICpara] = func_GetSetup(app)
    DICpara.winsize = app.winsize;
    DICpara.winstepsize = app.winstepsize;
    DICpara.gridxyROIRange = app.gridxyROIRange;
    DICpara.LoadImgMethod = app.LoadImgMethod;
    DICpara.ImgSeqIncUnit = app.ImgSeqIncUnit;
    DICpara.ImgSeqIncROIUpdateOrNot = app.ImgSeqIncROIUpdateOrNot;
    DICpara.Subpb2FDOrFEM = app.Subpb2FDOrFEM;
    DICpara.NewFFTSearch = app.NewFFTSearch;
    DICpara.ClusterNo = app.ClusterNo;
    DICpara.ImgSize = app.ImgSize;
    
    file_name = app.ImageFilePath;
    % read in all images and place them into the cell Img
    for i=1:size(file_name,2)
        Img{1,i} = imread(file_name{1,i});
        if ~isequal(size(Img{1,i},3),1)
            Img{1,i} = im2gray(Img{1,i});
        end
        Img{1,i} = double(Img{1,i});
    end
    
end

