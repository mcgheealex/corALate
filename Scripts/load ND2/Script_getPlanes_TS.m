try
    % use the getPlane function from BFToolbox to get the Z=plane,
    % C=color, T=frame, S=section of the image
    
    width = imgObj.width  ;
    height = imgObj.height;
    series = imgObj.series;
    S = imgObj.seriesCount;
    Z = imgObj.sizeZ;
    C = imgObj.sizeC;
    T = imgObj.sizeT;
    channelNames = imgObj.channelNames ;
    pxsz = imgObj.pxSize;
    swapZT = imgObj.swapZandT;
    
    if T > 1 % this ND2 file must have multiple time points
        for s = 1:S
            for c=1:C
                for t=1:T
                    % get all the Z images stacked up
                    for z=1:Z
                        IMstack(:,:,z) = double(getPlane(imgObj, Z, C, T,S));
                    end
                    % save the time,solor,series to the cell array of stacks
                    IMstacks{t,c,s} = IMstack;
                end
            end
        end
    else
        error('The ND2 must have multiple timepoints ')
    end
    
catch
end