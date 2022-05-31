try
    % use the getPlane function from BFToolbox to get the Z=plane,
    % C=color, T=frame, S=section of the image
    warning('off','all') % Suppress all the tiff warnings
    for f = 1:length(fnames)
        filename = fnames{f};
        tstack  = Tiff(filename);
        [I,J] = size(tstack.read());
        K = length(imfinfo(filename));
        data = zeros(I,J,K);
        data(:,:,1)  = tstack.read();
        for n = 2:K
            tstack.nextDirectory()
            data(:,:,n) = tstack.read();
        end
        
        IMstacks{f} = data;
        
    end
    warning('on','all')
catch
end