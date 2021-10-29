b2 = fir1(3, [0.5, 0.999], 'high');
h2 = b2' * b2;
g2 = imfilter(Im, h2);
imagesc(g2)

mask = imgaussfilt(g2,8);
mask(mask<65000)=0;
mask(mask~=0)=10000;
mask = imgaussfilt(mask,10);

mask(mask<10)=0;
mask(mask~=0)=1;
mask=mask+1;
mask(mask==2)=0;
imagesc(mask)