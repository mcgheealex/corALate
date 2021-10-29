function [Imgout,CBout] = Imageoverlay_andCB(BG, OL, Mask, alpha, PlotTitle, scaleBar , colorScale)


    ws = 150; % how much white space to add to the edge
    sidepadding = 10;
    toppadding = 50;
    Rez = 7;
    BS = ws/3;

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

    % prepare overlay (the right side has extea pixels of white
    % space (ws) to the right
    OL1 = zeros(size(OL,1));
    OLmask = ones(size(OL));
    
    % adjust the overlay to include the extra white space and shift any values at the top and bottom range
    OL2 = uint8(255 * mat2gray(OL1));
    OL2(OL2==0) = 1;
    OL2(OL2==255) = 254;
    rgbOL = ind2rgb(OL2, cmap2);

    % also adjust the mask and background to accomidate the colorbar
    M1=double(M);
    M1(M1==0)=NaN;
    M1=M1*alpha;

    BG1 = BG;
    % create a scalebar made of discreet colored blocks

    % make the squares for the colormap
    totLen = Rez*(BS+sidepadding); % this is the total length of the colorbar
    CM = ones(ws,totLen)+255;
    
    % make the colored squares for the colormap
    inS = 2; % inset length (the thickness of the black boarder)
    clrs = linspace(1,255,Rez);

    sy = round(linspace(toppadding,totLen,Rez));
    sy = [sy' sy'+BS];
    sx = sy*0;
    sx(:,1) = ws - sidepadding;
    sx(:,2) = sx(:,1) - BS;

    for cs = 1:Rez
        CM(sy(cs,1):sy(cs,2),end-sx(cs,1):end-sx(cs,2)) = 0;
        CM(sy(cs,1)+inS:sy(cs,2)-inS,end-sx(cs,1)+inS:end-sx(cs,2)-inS) = clrs(cs);
    end

    CM = ind2rgb(uint8(255 * mat2gray(CM)), cmap2);
    [~,Ty] = size(CM);
    figure
    imshow(CM, []);
    % add the colorbar values to the right of the corresponding color
    text(Ty-ws/1.5+2*sidepadding, toppadding + (BS)/2,num2str(round(max(Clims),1)),'FontSize',14) %top color
    text(Ty-ws/1.5+2*sidepadding, toppadding/2 + totLen/2 + (BS)/2,num2str(round(mean(Clims),1)),'FontSize',14) % middle color
    text(Ty-ws/1.5+2*sidepadding, totLen + (BS)/2,num2str(round(min(Clims),1)),'FontSize',14) %bottom color

    text(Ty-ws/1.5+2*sidepadding, toppadding/2 ,colorScale{3},'FontSize',14) % scale units
    
    F = getframe;
    CBout = frame2im(F);

%     [Tx,Ty] = size(M1);
%     % add a scalebar to the bottom left of the figure
% 
%     sLx = Tx - scaleBar{5}(1);% this is the start location of the scale bar
%     sLy = scaleBar{5}(2); % this is the start location of the scale bar
%     st = scaleBar{1}; % thickness
%     sL = scaleBar{2}; % length in pixels
% 
%     rgbOL(sLx:sLx+st,sLy:sLy+sL)=0;
%     M1(sLx:sLx+st,sLy:sLy+sL)=1;

    BGFig = imshow(BG1, [], 'Colormap', gray(256));
    BGFig.AlphaData = OLmask;
    hold on;
    OLFig = imshow(rgbOL);
    % Set opacity/transparency to something less than 1 (alpha).
    % 1 is the default and it means the last image is opaque and the image below can't be seen.
    OLFig.AlphaData = M1;
    
    F = getframe;
    Imgout = frame2im(F);

end