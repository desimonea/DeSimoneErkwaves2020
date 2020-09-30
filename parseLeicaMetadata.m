function [m] = parseLeicaMetadata(file)
%  Copyright (C) 2020  Alessandro De Simone
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% [m]= parseLeicaMetadata(file)
% it open the 'file' metadata and extract selected information -> stored in
% the struct m -> everythign is pretty self-explanatory
% The voxel is calculated from magnification, zoom and image size


display(['Now parsing:' file]);

% reads xml file and converts it into a matlab struct
f=xml2struct(file);

% checks which "Attachment" element contains HardwareSetting
indexHS=[];
for i=1:numel(f.Data.Image.Attachment)
    if(strcmp(f.Data.Image.Attachment{i}.Attributes.Name,'HardwareSetting'))
        indexHS=i;
        break
    end
end
if(isempty(indexHS))
    error('HardwareSetting not found in metatada file!')
end
att=f.Data.Image.Attachment{indexHS}.ATLConfocalSettingDefinition.Attributes;

% stores various information - everything is put in um
m=[];
m.Magnification=str2num(att.Magnification);
m.sizeX=str2num(f.Data.Image.ImageDescription.Dimensions.DimensionDescription{1}.Attributes.NumberOfElements);
m.sizeY=str2num(f.Data.Image.ImageDescription.Dimensions.DimensionDescription{2}.Attributes.NumberOfElements);

m.ScanSpeed=str2num(att.ScanSpeed);
m.PosX=str2num(att.StagePosX)*10^6; %in um
m.PosY=str2num(att.StagePosY)*10^6; %in um
m.PosZ=str2num(att.Begin)*10^6; %in um
m.Zoom=str2num(att.Zoom);

% calculates voxel size
voxX= 0.6054688 * (1024/m.sizeX)* (0.75/m.Zoom) * (25/m.Magnification);   % in um
voxY= 0.6054688 * (1024/m.sizeY)* (0.75/m.Zoom) * (25/m.Magnification);   % in um
voxZ=str2num(f.Data.Image.Attachment{indexHS}.ATLConfocalSettingDefinition.Quantity.Attributes.Value)*10^6; %in um
m.voxel=[voxX voxY voxZ];


% read number of sequentials
if(isfield(f.Data.Image.Attachment{indexHS},'LDM_Block_Sequential'))
    nSeq=numel(f.Data.Image.Attachment{indexHS}.LDM_Block_Sequential.LDM_Block_Sequential_List.ATLConfocalSettingDefinition);
    if(nSeq>1)
        atl=f.Data.Image.Attachment{indexHS}.LDM_Block_Sequential.LDM_Block_Sequential_List.ATLConfocalSettingDefinition{1};
    else
        atl=f.Data.Image.Attachment{indexHS}.LDM_Block_Sequential.LDM_Block_Sequential_List.ATLConfocalSettingDefinition;

    end
else
    nSeq=1; 
    atl=f.Data.Image.Attachment{indexHS}.ATLConfocalSettingDefinition;
end

m.nSeq = nSeq;

%reads laser power


nlines=numel(atl.AotfList(1).Aotf{2}.LaserLineSetting);

lasers=nan(nlines,2,nSeq);

for j=1:nSeq

     if(isfield(f.Data.Image.Attachment{indexHS},'LDM_Block_Sequential'))
        if(nSeq>1)
            atl=f.Data.Image.Attachment{indexHS}.LDM_Block_Sequential.LDM_Block_Sequential_List.ATLConfocalSettingDefinition{j};
        else
            atl=f.Data.Image.Attachment{indexHS}.LDM_Block_Sequential.LDM_Block_Sequential_List.ATLConfocalSettingDefinition;
        end
     else
        atl=f.Data.Image.Attachment{indexHS}.ATLConfocalSettingDefinition;
    end

      
    lasers(1,1,j)=str2num(atl.AotfList(1).Aotf{1}.LaserLineSetting.Attributes.LaserLine);
    lasers(1,2,j)=str2num(atl.AotfList(1).Aotf{1}.LaserLineSetting.Attributes.IntensityDev);

    for i=1:nlines
        lasers(i+1,1,j)=str2num(atl.AotfList(1).Aotf{2}.LaserLineSetting{i}.Attributes.LaserLine);
        lasers(i+1,2,j)=str2num(atl.AotfList(1).Aotf{2}.LaserLineSetting{i}.Attributes.IntensityDev);
    end
end

m.lasers=lasers;

% read detectors
nDetectors=numel(atl.DetectorList.Detector);

detectors=nan(nDetectors,3);
for i=1:nDetectors
    detectors(i,1)=i;
   
    if(isfield(atl.DetectorList.Detector{i}.Attributes,'Gain'));
        detectors(i,2)=str2num(atl.DetectorList.Detector{i}.Attributes.Gain);
    else
        detectors(i,2)=NaN;
    end
    
    if(isfield(atl.DetectorList.Detector{i}.Attributes,'IsActive'));
        detectors(i,2)=str2num(atl.DetectorList.Detector{i}.Attributes.IsActive);
    else
        detectors(i,3)=NaN;
    end
end

m.detectors=detectors;
m.ActiveDetectors=detectors(detectors(:,3)==1,1);

% time-stamp

fp=xml2struct([file(1:end-4) '_Properties.xml']);

try
timestamp=fp.Data.Image.ImageDescription.StartTime.Text;
m.timestamp=timestamp;

delim=regexp(m.timestamp,' ');
datestring=timestamp(1:delim(2));

dt=datetime(datestring);
ampm=timestamp(delim(2)+1:delim(2)+2);

if(strcmp(ampm,'AM'))
    if(hour(dt)==12)
        dt=datetime(year(dt),month(dt),day(dt),0,minute(dt),second(dt));
    end
else
    if(~(hour(dt)==12))
        dt=datetime(year(dt),month(dt),day(dt),hour(dt)+12,minute(dt),second(dt));
    end
end

% old and to be deleted
% if(size(fp.Data.Image.ImageDescription.Dimensions.DimensionDescription,2)==4)
%     tst=fp.Data.Image.ImageDescription.Dimensions.DimensionDescription{4}.Attributes;
%     
%     timeFrames=tst.NumberofElements
%     
%     if(~strcmp(tst.DimID,'T'))
%         error('Some mess with time in the time-stamp');
%     end
%     tLength=tst.Length;
%     
%     delimiterH=regexp(tLength,'h');
%     if(~isempty(delimiterH))
%         hh=str2num(tLength(1:delimiterH-1));
%         startMin=delimiterH+1;
%     else
%         hh=0;
%         startMin=1;
%     end
%     
%     delimiterM=regexp(tLength,'m');
%     if(~isempty(delimiterM))
%         mm=str2num(tLength(startMin:delimiterM-1));
%         startS=delimiterM+1;
%     else
%         mm=0;
%         startS=1;
%     end
%     
%     delimiterS=regexp(tLength,'s');
%     if(~isempty(delimiterS))
%         ss=str2num(tLength(startS:delimiterS-1));
%     else
%         ss=0;
%     end
%     hTimeStep=datenum(0,0,0,hh,mm,ss)*24;
% else
%     hTimeStep=0; 
% end

m.dt=dt;
m.hTimeStamp=datenum(year(dt),month(dt),day(dt),hour(dt),minute(dt),second(dt))*24;
end

end

