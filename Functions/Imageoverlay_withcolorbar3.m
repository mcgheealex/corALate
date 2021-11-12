function [Imgout] = Imageoverlay_withcolorbar3(BG, OL, Mask, alpha, PlotTitle, scaleBar , colorScale)

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

    ws = 150; % how much white space to add to the edge
    sidepadding = 10;
    toppadding = 0;
    CBpadding = 50;
    Rez = 7;
    BS = ws/3;
    fontsize = 18;
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
    cmap2(256,:) = 0.94;

    % prepare overlay (the right side has extea pixels of white
    % space (ws) to the right
    OL1 = [OL zeros(size(OL,1),ws)];
    OL1 = [zeros(toppadding,size(OL1,2)); OL1 ];

    OLmask = [ones(size(OL)) zeros(size(OL,1),ws)];
    OLmask = [ zeros(toppadding,size(OLmask,2)); ones(size(OLmask))];
    OLmask(:,end-ws:end)=0.25;
    OLmask(1:toppadding,:)=0.25;
    % adjust the overlay to include the extra white space and shift any values at the top and bottom range
    OL2 = uint8(255 * mat2gray(OL1));
    OL2(OL2==0) = 1;
    OL2(OL2==255) = 254;
    rgbOL = ind2rgb(OL2, cmap2);

    % also adjust the mask and background to accomidate the colorbar
    M1 = [M zeros(size(BG,1),ws)];
    M1 = [zeros(toppadding,size(M1,2)); M1 ];
    M1=double(M1);
    M1(M1==0)=NaN;
    M1=M1*alpha;

    scf = uint16(max(BG,[],'all')); % max scale factor of background
    BG1 = [BG uint16(ones(size(BG,1),ws))*scf];
    BG1 = [uint16(ones(toppadding,size(BG1,2)))*scf; BG1 ];

    % create a scalebar made of discreet colored blocks

    % make the squares for the colormap
    CM = ones(size(OL1))+255;
    % make the colored squares for the colormap
    inS = 2; % inset length (the thickness of the black boarder)
    clrs = linspace(1,255,Rez);

    totLen = Rez*(BS+sidepadding); % this is the total length of the colorbar
    sy = round(linspace(toppadding+CBpadding,totLen,Rez));
    sy = [sy' sy'+BS];
    sx = sy*0;
    sx(:,1) = ws - sidepadding;
    sx(:,2) = sx(:,1) - BS;

    for cs = 1:Rez
        CM(sy(cs,1):sy(cs,2),end-sx(cs,1):end-sx(cs,2)) = 0;
        CM(sy(cs,1)+inS:sy(cs,2)-inS,end-sx(cs,1)+inS:end-sx(cs,2)-inS) = clrs(cs);
    end

    CM = ind2rgb(uint8(255 * mat2gray(CM)), cmap2);

    [Tx,Ty] = size(M1);

    % add a scalebar to the bottom left of the figure

    sLx = Tx - scaleBar{5}(1);% this is the start location of the scale bar
    sLy = scaleBar{5}(2); % this is the start location of the scale bar
    st = scaleBar{1}; % thickness
    sL = scaleBar{2}; % length in pixels

    rgbOL(sLx:sLx+st,sLy:sLy+sL)=0;
    M1(sLx:sLx+st,sLy:sLy+sL)=1;
    
    ogim = size(OL1);
    B1 = figure('visible','off');
    B1.Position = [1,1,ogim(2)+100,ogim(1)+100];
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
    text(Ty-ws/1.5+2*sidepadding, toppadding+CBpadding + (BS)/2,num2str(round(max(Clims),1)),'FontSize',fontsize) %top color
    text(Ty-ws/1.5+2*sidepadding, (toppadding+CBpadding)/2 + totLen/2 + (BS)/2,num2str(round(mean(Clims),1)),'FontSize',fontsize) % middle color
    text(Ty-ws/1.5+2*sidepadding, totLen + (BS)/2,num2str(round(min(Clims),1)),'FontSize',fontsize) %bottom color

    text(Ty-ws/1.5+2*sidepadding, (toppadding+CBpadding)/2 ,colorScale{3},'FontSize',fontsize) % scale units

    % place the plot title at the center of the background image
    text(Ty/2,toppadding/2,PlotTitle,'FontSize',fontsize,'FontWeight',"bold",'HorizontalAlignment','center')

    % add scalebar text
    text(sLy+sL/2,sLx - 4*st,[num2str(round(scaleBar{2}*scaleBar{4})) ' ' scaleBar{3}],'FontSize',fontsize,'Color','black','HorizontalAlignment',"center")

    F = getframe;
    Imgout = frame2im(F);

end