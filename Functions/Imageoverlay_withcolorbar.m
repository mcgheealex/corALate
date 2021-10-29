function [Imgout] = Imageoverlay_withcolorbar(BG, OL, Mask, alpha, PlotTitle, scaleBar , colorScale)
%  IMOVERLAY_FUN This function takes in two or more 2D matricies and overlays the second image 
% over top the first using some alpha map.
% 
% options include :
% 
% alpha sets the opacity of the overlay 
% 
% alpha should be a number between 0 and 1
% 
% place scale bar at the bottom of the screen in 
% 
% scale bar should be a cell vector of { pixels thick, pixels ling, units , 
% px2unit conversion, location from bottom left edge} example { 5, 100, 'um', 1.6, [200,200]} 
% 
% place a colorscale bar to the right of the overlay
% 
% colorscale should be a cell vector of { min value , max value, text description 
% of units }  example {-0.1,0.1,' eyy '}

ws = 100; % how much white space to add to the edge
padding = 10;
Rez = 7;
BS = ws/2;

% normalize the mask to get 0,1 values
    M = Mask./max(Mask,[],'all');
    
    % ensure both images have the same size
    [rowsa, colsa, ~] = size(BG);
    [rowsb, colsb, ~] = size(OL);
    if rowsa ~= rowsb || colsa ~= colsb
       	OL = imresize(OL, [rowsa, colsa], 'nearest');
    end
    
    Clims = [colorScale{1},colorScale{2}];
    
    %make colormap where the first and last values are used for the overlay
    cmap2 = jet(256);
    cmap2(1,:) = 0;
    cmap2(256,:) = 255;
    
    % prepare overlay (the right side has 50 extea pixels of white space
    % for the image
    OL1 = [OL zeros(size(OL,1),ws)];
    OLmask = [ones(size(OL)) zeros(size(OL,1),ws) + 0.25];
    
    % adjust the overlay to include the extra white space and shift any values at the top and bottom range  
    OL2 = uint8(255 * mat2gray(OL1));
    OL2(OL2==0) = 1;
    OL2(OL2==255) = 254;
    rgbOL = ind2rgb(OL2, cmap2);
    
    % make the squares for the colormap
    CM = ones(size(OL1))+255;
    % make the colored squares for the colormap
    inS = 2; % inset length (the thickness of the black boarder)
    clrs = linspace(1,255,Rez);
    
    % block size 
    % padding
    % Rez
    
    totLen = Rez*(BS+padding); % this is the total length of the colorbar
    sy = round(linspace(padding,totLen,Rez));
    sy = [sy' sy'+BS];
    sx = sy*0;
    sx(:,1) = ws - padding;
    sx(:,2) = sx(:,1) - BS;

    for cs = 1:Rez
        CM(sy(cs,1):sy(cs,2),end-sx(cs,1):end-sx(cs,2)) = 0;
        CM(sy(cs,1)+inS:sy(cs,2)-inS,end-sx(cs,1)+inS:end-sx(cs,2)-inS) = clrs(cs);
    end
    
    CM = ind2rgb(uint8(255 * mat2gray(CM)), cmap2);

    % also adjust the mask and background to accomidate the colorbar
    M1 = [M zeros(size(BG,1),ws)];
    M1=double(M1);
    M1(M1==0)=NaN;
    M1=M1*alpha;
    
    % add a scalebar to the bottom of the figure
    rgbOL(200:200+scaleBar{1},400-scaleBar{2}-20:400-20,:)=0;
    M1(200:200+scaleBar{1},400-scaleBar{2}-20:400-20)=1;
    
    BG1 = [BG uint16(ones(size(BG,1),ws))*uint16(max(BG,[],'all'))];
    
    CBarFig = imshow(CM, []);
    hold on
    BGFig = imshow(BG1, [], 'Colormap', gray(256));
    BGFig.AlphaData = OLmask;
    hold on;
    
    OLFig = imshow(rgbOL);
    % Set opacity/transparency to something less than 1 (alpha).
    % 1 is the default and it means the last image is opaque and the image below can't be seen.
    OLFig.AlphaData = M1;
        
    % add the colorbar values to the right of the corresponding color
    [Tx,Ty] = size(M1);
    text(totLen,Ty-ws/2,num2str(min(Clims)))
    text(padding,Ty-ws/2,num2str(max(Clims)))
    
    % place the plot title at the center of the background image
    title(PlotTitle)
    set(get(gca,'title'),'Position',[200 1 1.00011],'HorizontalAlignment',"center")
   
    % add scalebar text
    sbX = 0;%round((rowsa-scaleBar{5}(1)-scaleBar{2}-40)/2);
    sbY = 0;%scaleBar{5}(2);
    text(sbX,sbY,[num2str(round(scaleBar{2}*scaleBar{4})) ' ' scaleBar{3}],'Color','black','HorizontalAlignment',"center",'FontWeight',"bold")
    
    F = getframe;
    Imgout = frame2im(F);
    
end