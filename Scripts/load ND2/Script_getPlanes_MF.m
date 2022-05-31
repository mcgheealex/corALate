try
    % use the getPlane function from BFToolbox to get the Z=plane,
    % C=color, T=frame, S=section of the image
   alltimes = 1;
   
    for tt = 1:length( IMobjs)
        
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
        
        for s = 1:S
            for c=1:C
                for t = 1:T
                    % get all the Z images stacked up
                    for z=1:Z
                        IMstack(:,:,z) = double(getPlane(imgObj, z, c, t,s));
                    end
                    % save the time,solor,series to the cell array of stacks
                    IMstacks{alltimes,c,s} = IMstack;
                    alltimes = alltimes+1;
                end
            end
        end
        
    end
    
catch
end