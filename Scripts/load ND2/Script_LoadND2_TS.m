try

    try
        imgObj = BioformatsImage(fname);
        reader = imgObj.bfReader;
        s=reader.getSeriesMetadata();
        g=reader.getGlobalMetadata();
        
        javaMethod('merge', 'loci.formats.MetadataTools', g, s, 'Global ');
        metadata = s;
        % get the time of each scan
        if imgObj.sizeT > 1
            % how many digits ar in the numbers
            n = num2str(floor(log10(imgObj.sizeT))+1); % the number of digits in the number (i.e 279 has 3 digits)
            for i=1:imgObj.sizeT
                ND2.Times{i} = metadata.get(['timestamp #' sprintf(['%0' n 'd'],i)]);
            end
        end

        key = {
            'Global CH1ChannelDyeName'
            'Global CH1PMTHighVoltage'
            'Global CH1LaserPower'
            'Global CH1ChannelLaserIndex'
            
            'Global CH2ChannelDyeName'
            'Global CH2PMTHighVoltage'
            'Global CH2LaserPower'
            'Global CH2ChannelLaserIndex'
            
            'Global CH3ChannelDyeName'
            'Global CH3PMTHighVoltage'
            'Global CH3LaserPower'
            'Global CH3ChannelLaserIndex'
            
            'Global CH4ChannelDyeName'
            'Global CH4PMTHighVoltage'
            'Global CH4LaserPower'
            'Global CH4ChannelLaserIndex'
            
            'Global dExposureTime'
            'Global {Pinhole Size(um)} #1'
            'Global m_sMicroscopePhysShortName'
            };
        
        for i=1:length(key)
            k = key{i};
            META{i} = metadata.get(k);
        end
        
        if imgObj.sizeC > 1
            ND2.LTX = 1;
            ND2.LTY = 2;
        else
            ND2.LTX = 1;
            ND2.LTY = 1;
        end
        
        ND2.meta = META;
        
        ND2.CH1_Name.Text = ['CH 1 ' META{1} ];
        ND2.CH1_Gain.Value = str2num(META{2});
        ND2.CH1_Power.Value = str2num(META{3});
        
        ND2.CH2_Name.Text = ['CH 2 ' META{5} ];
        ND2.CH2_Gain.Value = str2num(META{6});
        ND2.CH2_Power.Value = str2num(META{7});
        
        ND2.CH3_Name.Text = ['CH 3 ' META{9} ];
        ND2.CH3_Gain.Value = str2num(META{10});
        ND2.CH3_Power.Value = str2num(META{11});
        
        ND2.CH4_Name.Text = ['CH 4 ' META{13} ];
        ND2.CH4_Gain.Value = str2num(META{14});
        ND2.CH4_Power.Value = str2num(META{15});
        
        ND2.PinholeSizeEditField.Value = [META{18} ' ?m'];
        ND2.MicroscopeName.Text = META{19};
    catch ME
        fprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ME.stack(1).name, ME.stack(1).line, ME.message)
    end
    
    ND2.ResolutionEditField.Value = [num2str(imgObj.width) 'x' num2str(imgObj.height)];
    ND2.meta{20} = [num2str(imgObj.width) 'x' num2str(imgObj.height)];
    ND2.ZstackEditField.Value = num2str(imgObj.sizeZ);
    ND2.meta{21} = num2str(imgObj.sizeZ);
    ND2.LasertypeEditField.Value = num2str(imgObj.sizeC);
    ND2.meta{22} = num2str(imgObj.sizeC);
    ND2.TimeseriesEditField.Value = num2str(imgObj.sizeT);
    ND2.meta{23} = num2str(imgObj.sizeT);
    ND2.PixelsizeEditField.Value = [num2str(imgObj.pxSize) ' ' imgObj.pxUnits] ;
    ND2.meta{24} = [num2str(imgObj.pxSize) ' ' imgObj.pxUnits] ;
    ND2.bitDepthEditField.Value = num2str(imgObj.bitDepth);
    ND2.meta{25} = num2str(imgObj.bitDepth);
    ND2.LaserTypes = imgObj.channelNames;
    ND2.meta{26} = imgObj.channelNames;
    % update all controls
    
    ND2s{1} = ND2;
    IMobjs{1} = imgObj;
    
catch ME
    
end