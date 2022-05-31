function projection = func_project(Volume,projax)

switch projax
    case 'XY'
        projection = imresize(squeeze(max(Volume,[],3)),[1000,1000]);
    case 'XZ'
        projection = imresize(squeeze(max(Volume,[],2)),[1000,1000]);
    case 'YZ'
        projection = imresize(squeeze(max(Volume,[],1)),[1000,1000]);
end

end
