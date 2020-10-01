%% De Simone et al., Control of osteoblast regeneration by a train of Erk activity waves 
%% Sample code for image processing
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
% 
% This code will perform image processing to quantify Erk activity and
% tissue flows from images of regenerating zebrafish scales.
%
% The code requires MATLAB (MathWorks).
% The code runs in MATLAB_R2016b. The code runs on a 
% iMac (27-inch, Late 2013) macOS Sierra
% Processor 3.5 GHz Intel Core i7
% Memory 32 GB 1600 MHz DDR3
%
% INSTALLATION AND USAGE
% The user can run the code by changing 
% MATLAB's present working directory to 'DeSimoneErkwaves2020-master'.
% MS Windows users have to adapt paths to MS Windows sintax.
% 
% Data files are too large to include within this Github repo and must be
% downloaded separately from:
% https://drive.google.com/file/d/117tJhklJNEQJ3ND4D76V_xl22kcwfzZ4/view?usp=sharing
% The downloaded folder must be unzipped and tne 'data' folder must be
% placed in 'DeSimoneErkwaves2020-master'.
%
% MATLAB Toolboxes required are:
% - Image Processing Toolbox (version 9.5)
% - Statistics and Machine Learning Toolbox (version 11.0)
% - Financial Toolbox (version 5.8)
% - Curve Fitting Toolbox (3.5.4)
%
% The code requires external MATLAB functions:
% loadtiff.m, saveastiff.m, weightedcov.m, xml2struct.m, nanconv.m,
% brewermap.m
% Those functions be downloaded from MathWorks. We attach them to the code together with
% their licenses.
%
% The main image processing routine is the script
% CREATE_H2A_ERKKTR_EXAMPLE.M
% The script can run altogether or each section sequentially.
% In some instances, user input would be required, but
% saved parameters are used for demo purposes. Users wishing to 
% use custom parameters can change USESAVEDPARAMS to false
% in the appropriate sections.
%
% Sample outputs for each steps are provided in the data folder.
%
% Scale data are saved in objects files in the "objects" folder (myScale
% data struct). Stacks are saved as tiff files.
%
% Variables that an user could change are briefly
% described. 

% DATA SAMPLE
%
% A sample time-series is provided. Erk activity and tissue flow maps can
% be generated using the code. 
%
% Time-points include two channels:
% ch1 -> Erk sensor (osx:ERKKtr-mCerulean)
% ch2 -> osteoblast nuclear marker (osx:H2A-mCherry)
%
% Expected output for each step is provided. 
% 16 Gb RAM is required
% Data sample runtime: ~3h

%% Stitch scales
% The regenerating scale is usually image in 4-9 position that require
% need to be stitched in 3D. This section stitches the positions in a
% single stack.

pluckTime=datetime('2019-03-02 12:01:00'); %time of scsle plucking

rootnames={'fish1_','fish2_','fish3_','fish4_'}; %fish names 
scalenames={'scale1'}; %scale names

hours=[96:111]; %hpp labels

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/'; %folder where data is stored
paths.inFolder=  [paths.masterFolder 'raw/']; %input folder 
paths.outFolder= [paths.masterFolder 'stitched/']; %output folder
paths.objFolder= [paths.masterFolder 'objects/']; %objects folder

opts=[];
opts.refCh=2; %reference channel to use for stitching 
opts.targetVox=[0.6055,0.6055,0.6055]; %desidered voxel size (resizing may be applied)
opts.cropImage=0; %crop both stacks during stitching (0:no, 1:yes)
opts.verbose=0; %display images while stitching (0:no, 1:yes)
opts.inFolders=1; %files from the same experiment are contained in root folders with the same name (0:no, 1:yes)
opts.allSubs=0; % consider all position for stitching (0:no, 1:yes)
opts.timeStep = 1;
opts.resizeFactors=[0.2 0.2]; %image resizing factor used for stitching
opts.zBest=[]; % z to use for stitching ([position z]; [] for automatic choice)

mkdir(paths.objFolder);

for f=1:numel(rootnames)
    rootname=rootnames{f};
    for s=1:numel(scalenames)
        scalename=scalenames{s};
        for h=1:numel(hours)
          
            folder=[paths.inFolder rootname  scalename '_' num2str(hours(h)) 'hpp/']        
            basename=[rootname scalename '_' num2str(hours(h)) 'hpp'];
            if(~exist(folder,'dir'))
                continue
            end

                if(exist([folder 'MetaData/' basename '_zBest.txt'],'file'))
                    opts.zBest=textread([folder 'MetaData/' basename '_zBest.txt']);
                elseif(exist([folder 'MetaData/zBest.txt'],'file'))
                    opts.zBest=textread([folder 'MetaData/zBest.txt']);
                end
                
                [myScale,stitched] = stitchScale(paths,basename,opts);
                
                try
                    % write correct hpp
                    myScale.pluckTime=pluckTime;
                    myScale.hppTrue=myScale.metadata{1}.hTimeStamp-datenum(pluckTime)*24;
                    wpathmat=[paths.objFolder basename '.mat'];  
                    save(wpathmat,'myScale');          
                end
            
        end
    end
end

%% Add time, scale and fish
% This section reads the time-point filename and stores information on fish
% number, scale number and time labels.

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/'; %folder where data is stored
paths.objFolder= [paths.masterFolder 'objects/']; %objects folder

fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
scale='1'; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' '.mat']);
 
for i=1:numel(st_dir)

   disp('------------');
   name=st_dir(i).name;    
    
   display([st_dir(i).folder filesep st_dir(i).name]);
   load([st_dir(i).folder filesep st_dir(i).name]);
    
   name=myScale.basename;
   
   delimiters=[regexp(name,'_') numel(name)+1];
   % time - hpp
   poshpp=regexp(name,'hpp');
   if(~isempty(poshpp))
        postime=delimiters(find(delimiters<poshpp,1,'last'))+1;
        myScale.hpp =str2num(name(postime:poshpp-1));
   else
        myScale.hpp =[]; 
   end
   
   % time - dpp
   posdpp=regexp(name,'dpp');
   if(~isempty(posdpp))
        postime=delimiters(find(delimiters<posdpp,1,'last'))+1;
        myScale.dpp =str2num(name(postime:posdpp-1));
        myScale.hpp = myScale.dpp*24;
   else
        myScale.dpp =[]; 
   end

   % time - t
   post=regexp(name,'_t');
   if(~isempty(post))
        posEndtime=delimiters(find(delimiters>post,1,'first'))-1;
        myScale.t =str2num(name(post+2:posEndtime));
   else
        myScale.t=[]; 
   end
   
   % Mark and Find
   posMF=regexp(name,'Mark_and_Find_');
   if(~isempty(posMF))
       posEndMF=delimiters(find(delimiters>posMF,4,'first'))-1;
       posEndMF=posEndMF(end);
       myScale.MF=str2num(name(posMF+14:posEndMF));
   else
       myScale.MF=[];
   end

   % fish
   posfish=regexp(name,'fish');
   posnumber=delimiters(find(delimiters>posfish,1,'first'))-1;
   fish =name((posfish+4):posnumber);

   % scale
   posscale=regexp(name,'scale');
   posnumber=delimiters(find(delimiters>posscale,1,'first'))-1;
   scale =name((posscale+5):posnumber);

   myScale.fish = str2num(fish);
   myScale.scale= str2num(scale);

   save([st_dir(i).folder filesep st_dir(i).name],'myScale');
   
end

%% Resize stitched scales
% This section resizes stitched images to fasten calculation of optimal
% rotation angles.

chToResize = [2]; %reference channel to resize 

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/'; %folder where data is stored
paths.inFolder=  [paths.masterFolder 'stitched/']; %input images to resize
paths.outFolder= [paths.masterFolder 'stitched/resized4x/']; %output folder
paths.objFolder= [paths.masterFolder 'objects/']; %objects folder
paths.fieldName=  'resizeStitched'; %field in myScale struct where resizing information is stored

fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
scale='1'; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' '.mat']);

opts=[];
opts.resizeFactor = 0.25; %resizing factor

optsSave=[];
optsSave.compress='lzw';
optsSave.overwrite=true;

mkdir(paths.outFolder);

% resize 
for i=1:numel(st_dir)
    
    display([paths.objFolder filesep st_dir(i).name]);
    load([paths.objFolder filesep st_dir(i).name]);
 
    for ch=1:numel(chToResize)
        s=loadtiff([paths.inFolder filesep st_dir(i).name(1:end-4) '_ch' num2str(chToResize(ch)) '.tif']);
        startSize=size(s);
        s1 = imresize(s,opts.resizeFactor);
       
        xSize = size(s1,1);
        s1=permute(s1,[1 3 2]);
        s2 = imresize(s1,[xSize round(opts.resizeFactor.*startSize(3))]);
        s2=permute(s2,[1 3 2]);
        saveastiff(s2,[paths.outFolder filesep st_dir(i).name(1:end-4) '_ch' num2str(chToResize(ch)) '.tif'],optsSave);
    end
    
    myScale.(paths.fieldName).startSize=startSize;
    myScale.(paths.fieldName).resizeFactor=opts.resizeFactor;
    save([paths.objFolder filesep st_dir(i).name],'myScale');

end


%% Rotate scales - Create rotation paramaters using resized images - USER INPUT REQUIRED (SEE BELOW) 
% The stack is rotated so that scales are as parallel as possibile
% to the xy plane. Angles are automatically calculated and the user chooses
% cropping rectangles to contain only the rotated scales. 
%
% This section requires user input. When data sample is tested, this
% section is skipped as we already provide parameters.
% If users want to chose parameters, please set useSavedParams= false

useSavedParams = true;

if(~useSavedParams)
    
    refCh = 2;
    chToRotate = []; %channels to rotate at this step ([] for none)
    ths=[]; %rotation angles ([th1 th2 th3]; [] to calculate angles automatically)
    rects=[NaN NaN NaN NaN; NaN NaN NaN NaN; NaN NaN NaN NaN; 0 0 0 0];
    % cropping rectangles in different projections ([positionX positionY width height]; [NaN NaN NaN
    % NaN] to ask user input; [0 0 0 0] to not crop)
    zPos=0; % planes to use ([zMin zMax]; 0 to keep them all)

    paths=[];
    paths.masterFolder='data/H2A_ERKKTR_test/'; % folder where data is stored
    paths.inFolder=  [paths.masterFolder 'stitched/resized4x/']; %images to use to calculate rotation
    paths.outFolder= [paths.masterFolder 'rotated/']; % output folder (not used here)
    paths.objFolder= [paths.masterFolder 'objects/']; % objects folder

    paths.fieldName=  'resizeStitched'; %struct field where parameters are read
    paths.rotFieldName = 'rotScale';    %struct field where rotation parameters are stored

    fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
    scale='1'; 
    hpp='';
    st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' '.mat']);

    % calculate rotation
    for i=1:numel(st_dir)

        display([paths.objFolder filesep st_dir(i).name]);

        load([paths.objFolder filesep st_dir(i).name]);

        [myScale] = rotateScale(myScale,refCh,chToRotate,ths,rects,zPos,paths);
        myScale.(paths.rotFieldName).rotAngles = myScale.(paths.fieldName).rotAngles;
        myScale.(paths.rotFieldName).zPos =      round(myScale.(paths.fieldName).zPos./(myScale.(paths.fieldName).resizeFactor));

        resFactor=myScale.(paths.fieldName).resizeFactor;

        %rescale position and size
        rectsHere(:,1:2)= ((myScale.(paths.fieldName).cropRects(:,1:2)-0.51)./resFactor)+0.51;
        rectsHere(:,3:4)= ((myScale.(paths.fieldName).cropRects(:,3:4)+0.02)./resFactor)+4-1-0.02;   
        rectsHere(all(~myScale.(paths.fieldName).cropRects,2),:)=0;

        myScale.(paths.rotFieldName).cropRects = rectsHere;
        save([paths.objFolder filesep st_dir(i).name],'myScale');
    end
else
   disp('Saved parameters are used, please continue to next section'); 
end
%% Rotate scales - Apply rotation using stored parameters

refCh=2;
chToRotate = [1 2]; %channel to be processed

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/'; % folder where data is stored
paths.inFolder=  [paths.masterFolder 'stitched/']; % input folder for images to rotate
paths.outFolder= [paths.masterFolder 'rotated/'];  % output folder for rotated images
paths.objFolder= [paths.masterFolder 'objects/'];  % objects folder
paths.objFolderSaved= [paths.masterFolder 'objectsSaved/']; % provided objects folder
paths.fieldName=  'rotScale'; %struct where rotation parameters are read

useSavedParams = true; % use provided parameters (true: use provided parameters; false: use user parameters)

fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
scale='1'; 
hpp='';

st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' '.mat']);

% apply rotation
for i=1:numel(st_dir)
    
    display([paths.objFolder filesep st_dir(i).name]);
   
    load([paths.objFolderSaved filesep st_dir(i).name]);
    myScaleSaved = myScale;
    load([paths.objFolder filesep st_dir(i).name]);

    if(useSavedParams)
        ths=myScaleSaved.rotScale.rotAngles; 
        rects=myScaleSaved.rotScale.cropRects;
        zPos=myScaleSaved.rotScale.zPos;
    else
        ths = myScale.rotScale.rotAngles; 
        rects=myScale.rotScale.cropRects;
        zPos=myScale.rotScale.zPos;
    end
    
    [myScale] = rotateScale(myScale,refCh,chToRotate,ths,rects,zPos,paths);

end

%% Cleans up  scales - Create ROI file - USER INPUT REQUIRED (SEE BELOW)
% Scale stack include neighboring scales portions. Manual data curation is
% required to eliminate neighbouring scales. In this step, the user selects the
% region where the  selected scale is located in three different views. 
%
% This section requires user input. When data sample is tested, this
% section is skipped as we already provide parameters.
% If users want to chose parameters, please set useSavedParams= false

useSavedParams = true;

refCh=[2]; % reference channels to use for ROI selections
chToClean=[2];  % reference channels to use for manual cleaning step
roi=[]; % roi to use (stack; [] to ask user input)

paths=[];
paths.masterFolder = ['data/H2A_ERKKTR_test/']; % folder where data is stored
paths.inFolder=    [paths.masterFolder 'rotated/']; % input folder for images to cure
paths.outFolder=   [paths.masterFolder 'cleaned/']; % output folder for ROIs
paths.outFolderSaved=   [paths.masterFolder 'cleanedSaved/']; % output folder for ROIs
paths.objFolder=   [paths.masterFolder 'objects/']; %objects folder

if(~useSavedParams)
    
    fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
    scale='1'; 
    hpp='';
    st_dir=dir([paths.inFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' 'ch2.tif']);

    opts=[];
    opts.type='manual'; 
    opts.cleanupPlanes = 0;
    opts.adj = [0 0.2]; % the user can change image display adjustment ( 0 < [min max] < 1)

    for i=1:numel(st_dir)

        disp([paths.objFolder filesep st_dir(i).name(1:end-8) '.mat']);
        load([paths.objFolder filesep st_dir(i).name(1:end-8) '.mat']);
        [myScale] = cleanupScale(myScale,refCh,chToClean,roi,paths,opts);
        close all

    end
else
    
   mkdir(paths.outFolder);
   copyfile([paths.outFolderSaved '*ROI.tif'],[paths.outFolder])
   copyfile([paths.outFolderSaved '*ch' num2str(refCh) '*.tif'],[paths.outFolder])
   disp('Saved parameters are used, please continue to next section'); 

end

%% Create file to cleanup planes with ImageJ/Fiji

refCh=[2]; 
roi=[];

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/';  % folder where data is stored
paths.inFolder=  [paths.masterFolder 'cleaned/']; % input folder for stack to prepare for data curation
paths.outFolder= [paths.masterFolder 'cleaned/']; % output folder
paths.saveSuffix='_ROI2';

fish='1';  % put the number of the fish/scale/hpp to process or '' to choose them all
scale='1'; 
hpp='';
st_dir=dir([paths.inFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*'  'ch' num2str(refCh) '.tif']);

opts=[];
opts.adj=[0 0.2];  % the user can change image display adjustment ( 0 < [min max] < 1)

options=[];
options.compress='lzw';
options.overwrite=true;

for i=1:numel(st_dir)
    
    display([paths.inFolder st_dir(i).name]);
    
    stack=loadtiff([paths.inFolder st_dir(i).name]);
    
    adj = opts.adj*255;
    
    stack=stack+1;
    stack=double(stack);
    stack= 255.*(stack-adj(1))./(adj(2)-adj(1));
    stack=uint8(stack);
    
    saveastiff(stack,[paths.inFolder st_dir(i).name(1:end-8) paths.saveSuffix '.tif'],options);
    
end

%% Manual data curation using Fiji
% Here the user opens every ROI2 image using ImageJ/Fiji and manually
% removes neighboring scales plane by plane. 

% This section requires user input. When data sample is tested, this
% section is skipped as we already provide parameters.
% If users want to chose parameters, please set userSavedParams= false

useSavedParams = true;

if(useSavedParams)
   
   paths=[];
   paths.masterFolder = ['data/H2A_ERKKTR_test/']; % folder where data is stored
   paths.outFolder=   [paths.masterFolder 'cleaned/']; % output folder for ROIs
   paths.outFolderSaved=   [paths.masterFolder 'cleanedSaved/']; % output folder for ROIs

   copyfile([paths.outFolderSaved '*ROI2.tif'],[paths.outFolder])
   disp('Saved parameters are used, please continue to next section'); 

end

%% Merge ROI1 and ROI2 
% This step merges the ROIs selected during the two stages of data
% curation.

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/';  % folder where data is stored
paths.inFolder=  [paths.masterFolder 'cleaned/']; % folder where ROIs are stored
paths.outFolder= [paths.masterFolder 'cleaned/']; % output folder

fish='1';  % put the number of the fish/scale/hpp to process or '' to choose them all
scale='1'; 
hpp='';
st_dir=dir(strcat(paths.inFolder,['fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' ],'ROI.tif'));

saveOptions=[];
saveOptions.overwrite='true';
saveOptions.compress='lzw';

roi1=[];
roi2=[];

for i=1:numel(st_dir)
    
    display([paths.inFolder filesep st_dir(i).name]);
    
    roi1=loadtiff([paths.outFolder filesep st_dir(i).name]);
    roi2=loadtiff([paths.outFolder filesep st_dir(i).name(1:end-4) '2.tif']);
    saveastiff(im2uint8(roi2),[paths.outFolder filesep st_dir(i).name(1:end-4) '3.tif'],saveOptions);

    % binarize  
    roi1 = roi1>0;
    roi2=  roi2>0;
    
    % combine
    roinew=roi1&roi2;
    
    %save
    saveastiff(im2uint8(roinew),[paths.outFolder filesep st_dir(i).name(1:end-4) '3.tif'],saveOptions);
    
end

%% Apply ROI3 to clean images 
% Apply ROI3 to rotated images to get cleaned images

refCh=2; % reference channel (not used since you are providing ROIs)
chToClean=[1 2]; % channels to apply cleaning to
roi=[];

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/'; % folder where data is stored
paths.inFolder  =  [paths.masterFolder 'rotated/']; % input folder stacks to rotate
paths.roiFolder  =  [paths.masterFolder 'cleanedSaved/']; % folder provided ROIs
paths.outFolder =  [paths.masterFolder 'cleaned/']; % output folder
paths.objFolder =  [paths.masterFolder 'objects/']; % objects folder 

fish='1';  % put the number of the fish/scale/hpp to process or '' to choose them all
scale='1'; 
hpp='';
st_dir=dir(strcat(paths.roiFolder,['fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' ],'ROI3.tif'));

opts=[];
opts.type='manual';
 
for i=1:numel(st_dir)
    
    disp([paths.objFolder filesep st_dir(i).name(1:end-9) '.mat']);
    load([paths.objFolder filesep st_dir(i).name(1:end-9) '.mat']);
    
    roi=loadtiff([paths.roiFolder filesep st_dir(i).name]);
    [myScale] = cleanupScale(myScale,refCh,chToClean,roi,paths,opts);
end

%% Resize cleaned scales 
% This section resizes cleaned images to fasten calculation of optimal
% re-rotation angles.

chToResize = [2]; %reference channel for angles calculation

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/'; % folder where data is stored
paths.inFolder= [paths.masterFolder 'cleaned/']; % input folder cleaned images
paths.outFolder= [paths.masterFolder 'cleaned/resized4x/']; % output folder resized images
paths.objFolder= [paths.masterFolder 'objects/']; % objects folder
paths.fieldName=  'resizeCleaned'; % struct field resizing parameters

fish='1';  % put the number of the fish/scale/hpp to process or '' to choose them all
scale='1'; 
hpp=' ';
st_dir=dir([paths.inFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' 'ch' num2str(chToResize) '.tif']);

opts=[];
opts.resizeFactor = 0.25; % resizing factor 

optsSave=[];
optsSave.compress='lzw';
optsSave.overwrite=true;

mkdir(paths.outFolder);

for i=1:numel(st_dir)
    
    display([paths.objFolder filesep st_dir(i).name]);
    load([paths.objFolder filesep st_dir(i).name(1:end-8) '.mat']);
 
    for ch=1:numel(chToResize)
        s=loadtiff([paths.inFolder filesep st_dir(i).name(1:end)]);
        startSize=size(s);
        s1 = imresize(s,opts.resizeFactor);
       
        xSize = size(s1,1);
        s1=permute(s1,[1 3 2]);
        s2 = imresize(s1,[xSize round(opts.resizeFactor.*startSize(3))]);
        s2=permute(s2,[1 3 2]);
        saveastiff(s2,[paths.outFolder filesep st_dir(i).name(1:end)],optsSave);
    end
    
    myScale.(paths.fieldName).startSize=startSize;
    myScale.(paths.fieldName).resizeFactor=opts.resizeFactor;
    save([paths.objFolder filesep st_dir(i).name(1:end-8) '.mat'],'myScale');

end

%% Re-rotate scales - Calculate rotation parameters - USER INPUT REQUIRED (SEE BELOW)
% The stack is re-rotated so that the cleaned scale is as parallel as possibile
% to the xy plane. Angles are automatically calculated and the user chooses
% cropping rectangles to contain only the rotated scales. 
%
% This section requires user input. When data sample is tested, this
% section is skipped as we already provide parameters.
% If users want to chose parameters, please set useSavedParams= false

useSavedParams= true;

if(~useSavedParams)

    chToRotate = []; %channels to rotate at this step ([] for none)
    ths=[]; %rotation angles ([th1 th2 th3]; [] to calculate angles automatically)
    rects=[NaN NaN NaN NaN; NaN NaN NaN NaN; NaN NaN NaN NaN; 0 0 0 0];
    % cropping rectangles in different projections ([positionX positionY width height]; [NaN NaN NaN
    % NaN] to ask user input; [0 0 0 0] to not crop)
    zPos=0; % planes to use ([zMin zMax]; 0 to keep them all)

    paths=[];
    paths.masterFolder='data/H2A_ERKKTR_test/'; % folder where data is stored
    paths.inFolder=  [paths.masterFolder 'cleaned/resized4x/']; %images to use to calculate rotation
    paths.outFolder= [paths.masterFolder 'cleanedRerot/']; % output folder (not used here)
    paths.objFolder= [paths.masterFolder 'objects/']; % objects folder

    paths.fieldName=  'resizeCleaned'; %struct field where parameters are read
    paths.rotFieldName = 'rerotScale';    %struct field where rotation parameters are stored

    fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
    scale='1'; 
    hpp='';
    st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' '.mat']);

    resizeFactor = 0.25; % resizing factor

    % calculate rotations
    for i=1:numel(st_dir)
        display([paths.objFolder filesep st_dir(i).name]);
        load([paths.objFolder filesep st_dir(i).name]);

        myScale.(paths.fieldName).resizeFactor = resizeFactor;
        [myScale] = rotateScale(myScale,refCh,chToRotate,ths,rects,zPos,paths);
        myScale.(paths.rotFieldName).rotAngles = myScale.(paths.fieldName).rotAngles;
        myScale.(paths.rotFieldName).zPos = round(myScale.(paths.fieldName).zPos./(myScale.(paths.fieldName).resizeFactor));

        resFactor=myScale.(paths.fieldName).resizeFactor;

        %rescale position and size
        rectsHere(:,1:2)= ((myScale.(paths.fieldName).cropRects(:,1:2)-0.51)./resFactor)+0.51;
        rectsHere(:,3:4)= ((myScale.(paths.fieldName).cropRects(:,3:4)+0.02)./resFactor)+4-1-0.02;   
        rectsHere(all(~myScale.(paths.fieldName).cropRects,2),:)=0;

        myScale.(paths.rotFieldName).cropRects = rectsHere;
        save([paths.objFolder filesep st_dir(i).name],'myScale');
    end
else
    disp('Saved parameters are used, please continue to next section'); 
end

%% Re-rotate scales - Apply rotation using provided parameters

refCh=2;
chToRotate = [1 2]; %channel to be processed

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/'; % folder where data is stored
paths.inFolder=  [paths.masterFolder 'cleaned/']; % input folder for images to re-rotate
paths.outFolder= [paths.masterFolder 'cleanedRerot/'];  % output folder for re-rotated images
paths.objFolder= [paths.masterFolder 'objects/'];  % objects folder
paths.objFolderSaved= [paths.masterFolder 'objectsSaved/']; % provided objects folder
paths.fieldName=  'rerotScale'; %struct where re-rotation parameters are read

useSavedParams = true; % use provided parameters (true: use provided parameters; false: use user parameters)

fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
scale='1'; 
hpp='';

st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' '.mat']);

% apply re-rotations
for i=1:numel(st_dir)
    
    disp([paths.objFolder filesep st_dir(i).name]);
    
    load([paths.objFolderSaved filesep st_dir(i).name]);
    myScaleSaved = myScale;
    load([paths.objFolder filesep st_dir(i).name]);
   
    if(useSavedParams)
        ths   = myScaleSaved.rerotScale.rotAngles; 
        rects = myScaleSaved.rerotScale.cropRects;
        zPos  = myScaleSaved.rerotScale.zPos;
    else
        ths   = myScale.rerotScale.rotAngles; 
        rects = myScale.rerotScale.cropRects;
        zPos  = myScale.rerotScale.zPos;
    end
    
    [myScale] = rotateScale(myScale,refCh,chToRotate,ths,rects,zPos,paths);
    save([paths.objFolder filesep st_dir(i).name],'myScale');
end

%% Create time-series legend 

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/'; % folder where data is stored
paths.objFolder=    [paths.masterFolder 'objects/']; % objects folder
paths.legendFolder=  [paths.masterFolder 'objects/legend/']; % legend folder
mkdir(paths.legendFolder);

fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
scale='1'; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' '.mat']);

savename = 'fish1_scale1_96hpp_series'; %name time series legend
 
ts=[];MFs=[];hppTrues=[]; hpp=[];
for i=1:numel(st_dir)
    load([st_dir(i).folder filesep st_dir(i).name]);
    
    if(isfield('myScale','t'))
        ts(i)=myScale.t;    
    else
        ts(i)=NaN;
    end
    
     if(isfield('myScale','hpp'))
       hpp(i)=myScale.hpp;    
    else
       hpp(i)=NaN;
    end
    
    if(isfield('myScale','MF'))
        MFs(i)=myScale.MF;    
    else
        MFs(i)=NaN;
    end

    hppTrues(i)=myScale.hppTrue;
  
end

[hppTrues,idxs]=sort(hppTrues);

myScaleSeries=[];
myScaleSeries.ts=ts(idxs);
myScaleSeries.hppTrues=hppTrues;
myScaleSeries.hpp=hpp;
myScaleSeries.MFs=MFs(idxs);
myScaleSeries.objNames={st_dir(idxs).name};

save([paths.legendFolder savename ],'myScaleSeries');

for i=1:numel(st_dir)
    load([st_dir(i).folder filesep st_dir(i).name]);
    myScale.time_idx=idxs(i);
    save([st_dir(i).folder filesep st_dir(i).name],'myScale');
end

%% Register images
% This section registers the images for subsequent nuclei tracking. The
% procedure is repeated for times. 

% --- FIRST REGISTRATION

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/'; % folder where data is stored
paths.inFolder=    [paths.masterFolder    'cleanedRerot/']; % input folder to calculate first registration parameters
paths.toRegFolder= [paths.masterFolder    'cleanedRerot/']; % output folder firt registration
paths.outFolder=   [paths.masterFolder    'registered/']; % input folder stacks to register
paths.objFolder=   [paths.masterFolder    'objects/']; % objects folder
paths.legendFolder=[paths.masterFolder    'objects/legend/']; % legend folder

paths.inSuffix =    '.tif';
paths.outSuffix =   '.tif';
paths.toRegSuffix = '.tif';

opts.refCh=2; % channel to use to calculate registration parameters
opts.chToRegister=[1 2]; % channels to register
opts.verbose =0; % display intermediate steps registration (0: no; 1: yes)
opts.useCrossCorr = 1; % calculate a preliminary registration using cross correlation

[optimizer, metric] = imregconfig('multimodal'); % parameters registration
opts.optimizer = optimizer;
opts.metric = metric;

load([paths.legendFolder 'fish1_scale1_96hpp_series.mat']);
[regInfo]=imregisterStackSeriesXY3D(myScaleSeries,paths,opts); % registration step

% --- SECOND REGISTRATION 

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/';
paths.inFolder=    [paths.masterFolder    'registered/'];
paths.toRegFolder= [paths.masterFolder    'registered/'];
paths.outFolder=   [paths.masterFolder    'registered2/'];
paths.objFolder=   [paths.masterFolder    'objects/'];
paths.legendFolder=[paths.masterFolder    'objects/legend/'];

paths.inSuffix =    '.tif';
paths.outSuffix =   '.tif';
paths.toRegSuffix = '.tif';

opts.refCh=1;
opts.chToRegister=[1 2];
opts.verbose =0;
opts.useCrossCorr = 0;

[optimizer, metric] = imregconfig('monomodal');
opts.optimizer = optimizer;
opts.metric = metric;

[regInfo]=imregisterStackSeriesXY3D(myScaleSeries,paths,opts);

% --- THIRD REGISTRATION 

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/';
paths.inFolder=    [paths.masterFolder    'registered2/'];
paths.toRegFolder= [paths.masterFolder    'registered2/'];
paths.outFolder=   [paths.masterFolder    'registered3/'];
paths.objFolder=   [paths.masterFolder    'objects/'];
paths.legendFolder=[paths.masterFolder    'objects/legend/'];


paths.inSuffix =    '.tif';
paths.outSuffix =   '.tif';
paths.toRegSuffix = '.tif';

opts.refCh=1;
opts.chToRegister=[1 2];
opts.verbose =0;
opts.useCrossCorr = 0;

[optimizer, metric] = imregconfig('monomodal');
opts.optimizer = optimizer;
opts.metric = metric;

[regInfo]=imregisterStackSeriesXY3D(myScaleSeries,paths,opts);

% --- FOURTH REGISTRATION 

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/';
paths.inFolder=    [paths.masterFolder    'registered3/'];
paths.toRegFolder= [paths.masterFolder    'registered3/'];
paths.outFolder=   [paths.masterFolder    'registered4/'];
paths.objFolder=   [paths.masterFolder    'objects/'];
paths.legendFolder=[paths.masterFolder    'objects/legend/'];

paths.inSuffix =    '.tif';
paths.outSuffix =   '.tif';
paths.toRegSuffix = '.tif';

opts.refCh=1;
opts.chToRegister=[1 2];
opts.verbose =0;
opts.useCrossCorr = 0';

[optimizer, metric] = imregconfig('monomodal');
opts.optimizer = optimizer;
opts.metric = metric;

[regInfo]=imregisterStackSeriesXY3D(myScaleSeries,paths,opts);

%% Equalization nuclear signal for segmentation and tracking

chToEq=[2]; % channel to equalize

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/'; % folder where data is stored
paths.inFolder=     [paths.masterFolder 'registered4/']; % input folder stack to equalize
paths.outFolder=    [paths.masterFolder 'equalized/']; % output folder equalized stacks
paths.objFolder=    [paths.masterFolder 'objects/']; % objects folder

fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
scale='1'; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' '.mat']);

 % equalization
for i=1:numel(st_dir)
    display([st_dir(i).folder filesep st_dir(i).name]);
    load([st_dir(i).folder filesep st_dir(i).name]);
    myScale = equalizeScale(myScale,chToEq,paths);
end
%%  Entire scale segmentation in 3D

refCh=2; % reference channels for segmentation
paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/';  % folder where data is stored
paths.inFolder=     [paths.masterFolder 'registered4/']; % input folder stack to segment
paths.outFolder= [paths.masterFolder 'segmented/']; %segmented 3D ROIs
paths.objFolder=    [paths.masterFolder 'objects/']; % objects folder

fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
scale='1'; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' '.mat']);
verbose =0; % display segmentation (0: no; 1: yes)

for i=1:numel(st_dir)
   
    disp('------------');
    display([paths.objFolder filesep st_dir(i).name(1:end-4) '.mat']);
    load([paths.objFolder filesep st_dir(i).name(1:end-4) '.mat']);
    [myScale] = segmentMyScale(myScale,refCh,paths,[],verbose);
end
%% Scale flattening - Equalized stacks
% The scale is flattened by bringing its segmented hyposquamal border at the same
% z-position. 

toflattenCh=[2]; % channel to flatten

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/'; % folder where data is stored
paths.inFolder= [paths.masterFolder 'equalized/']; % input folder equalized stacks
paths.roiFolder= [paths.masterFolder 'segmented/']; % segmented ROI folder
paths.outFolder= [paths.masterFolder 'flatten_eq/']; % output folder flattened stacks
paths.objFolder=    [paths.masterFolder 'objects/']; % objects folder
paths.fitFieldName='findHypo'; % struct field flattening parameters

fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
scale='1'; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' '.mat']);

opts=[];
opts.method='fit'; % with 'fit' ROI is not flattened - launch giving toFlattenCh='ROI' to flatten ROI
 
for i=1:numel(st_dir)
    
    display([st_dir(i).folder filesep st_dir(i).name]);
    load([st_dir(i).folder filesep st_dir(i).name]);  
    [myScale] = flattenLayersScaleSegm(myScale,toflattenCh,paths,opts);
    
end

%% Scale flattening - Not equalized stacks
% The scale is flattened by bringing its segmented hyposquamal border at the same
% z-position. 

toflattenCh=[1 2]; % channel to flatten

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/'; % folder where data is stored
paths.inFolder=     [paths.masterFolder 'registered4/']; % input folder stacks
paths.roiFolder= [paths.masterFolder 'segmented/']; % segmented ROI folder
paths.outFolder= [paths.masterFolder 'flatten/']; % output folder flattened stacks
paths.objFolder=    [paths.masterFolder 'objects/']; % objects folder
paths.fitFieldName='findHypo';

fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
scale='1'; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' '.mat']);

opts=[];
opts.method='fit'; % with 'fit' ROI is not flattened - launch giving toFlattenCh='ROI' to flatten ROI
 
for i=1:numel(st_dir)
    
    display([st_dir(i).folder filesep st_dir(i).name]);
    load([st_dir(i).folder filesep st_dir(i).name]);  
    [myScale] = flattenLayersScaleSegm(myScale,toflattenCh,paths,opts);
    
end

%% Isolation hyposquamal layer - Equalized stacks
% The hyposquamal layer is computationally dissected . 
% To this end, the peak of the total plane intensity
% is calculated. The hyposquamal layer is taken from the start of the stack 
% to 15 planes past the z-peak of total plane intensity.

refCh=2; % reference channel for calculation z-peak
todivideCh=[2]; % channels to divide

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/'; % folder where data is stored
paths.inFolder=     [paths.masterFolder 'flatten_eq/']; % input folder stacks to divide
paths.refFolder=    [paths.masterFolder 'flatten_eq/']; % reference folder stack to calculate z-peak
paths.outFolder=    [paths.masterFolder 'divided_eq/']; % output folder for isolated hyposquamal layers
paths.objFolder=    [paths.masterFolder 'objects/']; % objects folders
paths.fitFieldName='findHypo'; % struct field name for hyposquaml layer isolation parameters

fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
scale='1'; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' '.mat']);
 
divOpts=[];
divOpts.verbose=0;    % display z-intensities
divOpts.methodEpi=15; % z-limit hyposquamal layer - planes past the peak in total z-intensity

for i=1:numel(st_dir)
    
    display([st_dir(i).folder filesep st_dir(i).name]);
    load([paths.objFolder st_dir(i).name]);

    if(myScale.hpp<48)
        divOpts.methodHypo='peak';
        [myScale,zprofile] = divideLayersScaleSegm(myScale,refCh,todivideCh,paths,divOpts);
    elseif(myScale.hpp<85)
        divOpts.methodHypo='peak2';
        [myScale,zprofile] = divideLayersScaleSegm(myScale,refCh,todivideCh,paths,divOpts);
    else
        divOpts.methodHypo='peak';
        [myScale,zprofile] = divideLayersScaleSegm(myScale,refCh,todivideCh,paths,divOpts);
    end

end

%% Isolation hyposquamal layer - Non equalized stacks
% The hyposquamal layer is computationally dissected . 
% To this end, the peak of the total plane intensity
% is calculated. The hyposquamal layer is taken from the start of the stack 
% to 15 planes past the z-peak of total plane intensity.

refCh=2;
todivideCh=[1 2];

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/'; % folder where data is stored
paths.inFolder=     [paths.masterFolder 'flatten/'];  % input folder stacks to divide
paths.refFolder=    [paths.masterFolder 'flatten_eq/']; % reference folder stack to calculate z-peak
paths.outFolder=    [paths.masterFolder 'divided_refEq/']; % output folder for isolated hyposquamal layers
paths.objFolder=    [paths.masterFolder 'objects/']; % objects folders

fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
scale='1'; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' '.mat']);
 
divOpts=[];
divOpts.verbose=0;
divOpts.methodEpi=15; % z-limit hyposquamal layer - planes past the peak in total z-intensity
  
for i=1:numel(st_dir)
    
    display([st_dir(i).folder filesep st_dir(i).name]);
    load([paths.objFolder st_dir(i).name]);

    if(myScale.hpp<48)
        divOpts.methodHypo='peak';
        [myScale,zprofile] = divideLayersScaleSegm(myScale,refCh,todivideCh,paths,divOpts);
    elseif(myScale.hpp<80)
        divOpts.methodHypo='peak2';
        [myScale,zprofile] = divideLayersScaleSegm(myScale,refCh,todivideCh,paths,divOpts);
    else
        divOpts.methodHypo='peak';
        [myScale,zprofile] = divideLayersScaleSegm(myScale,refCh,todivideCh,paths,divOpts);
    end
   
end

%% Prepare files for TGMM - Create config files
% This section prepares stacks and configuration files for TGMM
% segmentation. 

% TGMM has been published in:
% Fast, accurate reconstruction of cell lineages from large-scale fluorescence microscopy data (Amat et al. 2014, Nature Methods)
% and can be found at this link:
% https://www.janelia.org/lab/keller-lab/software/fast-accurate-reconstruction-cell-lineages-large-scale-fluorescence
%
% Please follow TGMM installation and usage manual.
%
% Configuration files paths are customized for MS Windows as
% we ran TGMM runs on that operative system. 
%
% TGMM must run on each config file. 
% TGMM can be run sequentially on each config file using a batch script.
% An example is provided: BATH_TGMM_H2A_ERKKTR_TEST.BAT 
% 
% When data sample is tested, this section is skipped as we already provide
% TGMM output in Matlab form.  Users can go to "START ERK  QUANTIFICATION"
% If users want to use this section, please set useSavedParams = false

useSavedParams = true;

if(~useSavedParams)

    paths=[];
    paths.nameFolder='H2A_ERKKTR_test/'; % folder where data is stored
    paths.masterFolder=[pwd '/data/' paths.nameFolder];

    tpl=           [paths.masterFolder 'TGMM_configFile.txt']; % config file template path
    configfolder=  [paths.masterFolder 'TGMM_hypo_eq_ch2/TGMMconfig/']; % config file output folder

    datasource   = [paths.masterFolder  'divided_eq/']; % input folder images to segment
    datapathmac  = [paths.masterFolder  'TGMM_hypo_eq_ch2/data/']; % output folder images to segment
    respathmac   = [paths.masterFolder  'TGMM_hypo_eq_ch2/results/']; % output folder TGMM output

    % paths on Windows computer
    datapathwin= ['G:\' paths.nameFolder(1:end-1) '\TGMM_hypo_eq_ch2\data\']; % input folder stacks on Windows computer
    respathwin=  ['G:\' paths.nameFolder(1:end-1) '\TGMM_hypo_eq_ch2\results\']; % output folder TGMM output on Windows computer

    fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
    scale='1'; 
    hpp='';
    st_dir=dir([datasource 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' '_ch2_hypo.tif']);

    mkdir(configfolder);
    mkdir(respathmac);
    mkdir(datapathmac);

    for i=1:numel(st_dir)
        prepareTGMM(st_dir(i).name,tpl,configfolder,datasource,datapathmac,respathmac,datapathwin,respathwin)
    end
else
    disp('Saved parameters are used, please continue to next section'); 
end

%% Here TGMM must be used to segment nuclei. 
% This step may be skipped as we provide TGMM output for the sample dataset.
% Users can go to "START ERK  QUANTIFICATION"

%% Organize TGMM output in Matlab structs
% When data sample is tested, this section is skipped as we already provide
% TGMM output in Matlab form.  Users can go to "START ERK  QUANTIFICATION"
% 
% If users want to use this section, please set useSavedParams = false

useSavedParams = true;

if(~useSavedParams)        

    respath='data/H2A_ERKKTR_test/TGMM_hypo_eq_ch2/results/'; % input folder TGMM output
    objpath='data/H2A_ERKKTR_test/TGMM_hypo_eq_ch2/objects/'; % output folder TGMM output in Matlab format
    mkdir(objpath);

    st_dir=dir([respath 'fish*scale*hpp*']);

    for i=1:numel(st_dir)

       gme=dir([respath st_dir(i).name filesep 'G*'])

       if(numel(gme)>0)
         f=dir([respath st_dir(i).name filesep gme(end).name  filesep 'XML_finalResult_lht' filesep 'GME*0000.xml']);

         if(numel(f)>0)

            fxml=([respath st_dir(i).name filesep gme(end).name filesep 'XML_finalResult_lht' filesep f(1).name]);
            obj=readXMLmixtureGaussians(fxml);
            [svList,sizeIm]=readListSupervoxelsFromBinaryFile([fxml(1:end-4) '.svb']);
         end
       end
       save([objpath st_dir(i).name],'obj','svList','sizeIm');   
    end

else
    disp('Saved parameters are used, please continue to next section'); 
end

%% START ERK  QUANTIFICATION

%% Measure Erk activity in individual cells

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/'; % folder where data is stored
paths.inFolder=     [paths.masterFolder 'divided_refEq/']; % input folder non-equalized Erk sensor signal
paths.refFolder=    [paths.masterFolder 'flatten_eq/']; % input folder equalized nuclear signals
paths.objFolder=    [paths.masterFolder 'objects/']; % myScale objects folder
paths.tgmmFolder=   [paths.masterFolder 'TGMM_hypo_eq_ch2/objects/']; % TGMM objects folder
paths.KTRsuffix = '_ch1_hypo'; % suffix Erk sensor signal stacks
paths.H2Asuffix = '_ch2_hypo'; % suffix nuclear signal stacks

opts=[];
opts.infoName='TGMM_hypo_eq_ch2'; % struct where TGMM information must be stored
opts.saveName='ktr'; % name ktr field
opts.verbose = 0; % display steps processing (0: none; 1: some; 2: all)

fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
scale='1'; 
hpp='';
st_dir=dir([paths.tgmmFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' '.mat']);
  
for i=1:numel(st_dir)

   disp('------------');
   display([st_dir(i).folder filesep st_dir(i).name]); 
 
   load([st_dir(i).folder filesep st_dir(i).name]); %load TGMM
   load([paths.objFolder filesep st_dir(i).name(1:end-13) '.mat']); %load myScale
   
   [myScale,obj] = assignKTRReduced(myScale,obj,svList,paths,opts);

   save([paths.tgmmFolder filesep st_dir(i).name],'obj','-append');
   save([paths.objFolder filesep st_dir(i).name(1:end-13) '.mat'],'myScale'); 
end

%% Select landmark to orient scales. USER INPUT REQUIRED (SEE BELOW)
% In this section, the user selects a point on the posterior (right edge) of
% the scale for orientation. Since time-points are registered, user input is
% required only for the first time-point. 

% When data sample is tested, users can use the already provided
% orientation parameters.
% If users want to choose orientation parameters, please set useSavedParams = false

useSavedParams = true;

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/'; % folder where data is stored
paths.inFolder=  [paths.masterFolder  'registered4/']; % input folder for stacks to orient
paths.objFolder= [paths.masterFolder  'objects/'     ]; % objects folder
paths.suffix='_ch1maxproj'; % channel to use for orientation
paths.fieldName= 'frontDirection'; % struct field name for orientation informarmation

fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
scale='1'; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' '.mat']);
  
i = 1;
obj=[];
disp(['------------ ' st_dir(i).folder filesep st_dir(i).name ' ------------']);
load([st_dir(i).folder filesep st_dir(i).name]) %myScale    

if(useSavedParams)
    myScale.(paths.fieldName).landmark=[1583  798];
else
    stack=loadtiff([paths.inFolder st_dir(i).name(1:end-4) paths.suffix '.tif']);
    img=max(stack,[],3);
    imshow(img,[]);
    [xd,yd,~]=impixel;
    myScale.(paths.fieldName).landmark=[xd(1) yd(1)];
end    

landmarkAll=myScale.(paths.fieldName).landmark;
save([st_dir(i).folder filesep st_dir(i).name],'myScale'); %myScale    
 
for i=2:numel(st_dir)
   obj=[];
   disp(['------------ ' st_dir(i).folder filesep st_dir(i).name ' ------------']);
   load([st_dir(i).folder filesep st_dir(i).name]) %myScale    
   myScale.(paths.fieldName).landmark= landmarkAll;
   save([st_dir(i).folder filesep st_dir(i).name],'myScale'); %myScale    
end

%% Generate Erk activity maps

paths=[];
paths.masterFolder='data/H2A_ERKKTR_test/'; % folder where data is stored
paths.tgmmFolder=   [paths.masterFolder 'TGMM_hypo_eq_ch2/objects/']; % tgmm output folder
paths.objFolder=   [paths.masterFolder 'objects/']; % myScale objects folder
paths.stackFolder= [paths.masterFolder 'divided_refEq/']; % input folder Erk sensor images
paths.outFolder=   [paths.masterFolder 'results/quantifications/']; % output folder quantifications
paths.suffix='_ch1_hypo';

fish='1'; % put the number of the fish/scale/hpp to process or '' to choose them all
scale='1'; 
hpp='';
st_dir=dir([paths.objFolder 'fish'  num2str(fish) '*scale' num2str(scale) '*' num2str(hpp) 'hpp*' '.mat']);
  
opts=[];
opts.infoName='TGMM_hypo_eq_ch2'; % struct field name for TGMM information
opts.statName='stat_hypo_eq_ch2'; % struct field name for Erk information
opts.adj=[0 254];
opts.print=true; % print images
opts.visualType='nuclei'; %  'centers' vs 'nuclei';
opts.limValue=[0.7 1.6]; % min and max values log(Erk_ratio) for visualization
opts.imgbox=[2000 2000]; % quantification image size
opts.dirField='frontDirection'; % struct field name with scale orientation information
opts.chooseSource=false; % the user can choose manually the source location
opts.infoName='TGMM_hypo_eq_ch2';
opts.valueToPlot='ktr'; 
opt.thetaLandmark=0; % rotation to apply

hppTrues=nan(1,numel(st_dir));
for i=1:numel(st_dir)
    load([st_dir(i).folder filesep st_dir(i).name]); %myScale
    hppTrues(i)=myScale.hppTrue;
end

[~,idxHpp]=sort(hppTrues);

for i=idxHpp
   obj=[];
   disp('------------');
   display([st_dir(i).folder filesep st_dir(i).name]);
   load([st_dir(i).folder filesep st_dir(i).name]); %myScale
   load([paths.tgmmFolder filesep st_dir(i).name(1:end-4) '_ch2_hypo.mat']); %load TGMM
   
   opts.obj=obj;
   opts.svList=svList;
     
   [myScale] = visualizeERKKTR(myScale,paths,opts);
   save([st_dir(i).folder filesep st_dir(i).name],'myScale');
   close all;
   
end

%% TISSUE FLOW CALCULATION

%% Create assembled stack for Ilastik (z-stack and time-points)
% This step generates stacks to be tracked using Ilastik.
% After generating the stack the user should use ImageJ/Fiji to organize
% the stack in a Hyperstack

paths=[];
paths.nameFolder = 'H2A_ERKKTR_test/'; % folder where data is stored
paths.masterFolder=['data/' paths.nameFolder];
paths.inFolder=    [paths.masterFolder  'divided_eq/']; % input folder stacks for tracking
paths.outFolder=    [paths.masterFolder 'divided_eqIlastik/']; % output folder stack s for ilastik 
paths.objFolder=    [paths.masterFolder 'objects/']; % objects folder

paths.inSuffix = '_ch2_hypo.tif'; % suffix images to be prepared for ilastik
paths.scaleName = 'fish1_scale1_';  
paths.outName = [paths.scaleName 'ch2_hypo_all.tif']; % name output stack

hpps = [96:3:108 111]; % time-points to be assembled

nz =18; % planes per z-stack

sAll = [];
saveOptions=[];
saveOptions.overwrite = true;

for h=1:numel(hpps)
  disp([paths.inFolder paths.scaleName num2str(hpps(h)) 'hpp' paths.inSuffix]);
  s=loadtiff([paths.inFolder paths.scaleName num2str(hpps(h)) 'hpp' paths.inSuffix]);
  snew = zeros([size(s,1) size(s,2) nz],'uint8');
  snew(:,:,end-size(s,3)+1:end)=s;
  sAll = cat(3,sAll,snew);
end

saveastiff(sAll,[paths.outFolder paths.outName],saveOptions);

%% ilastik segmentation and tracking.
% Here the user should use ilastik to segment and track cells, using 
%first the Pixel Classification and then Animal tracking pipelines. 
% ilastik was published in:
% ilastik: interactive machine learning for (bio)image analysis
% Stuart Berg, Dominik Kutra, Thorben Kroeger, Christoph N. Straehle, Bernhard X. Kausler, Carsten Haubold, Martin Schiegg, Janez Ales, Thorsten Beier, Markus Rudy, Kemal Eren, Jaime I Cervantes, Buote Xu, Fynn Beuttenmueller, Adrian Wolny, Chong Zhang, Ullrich Koethe, Fred A. Hamprecht & Anna Kreshuk
% in: Nature Methods, (2019)
% and can be downloaded at https://www.ilastik.org/
% Please refer to ilastik installation guide.

%% Plot tissue flows and growth maps

paths=[];
paths.nameFolder = 'H2A_ERKKTR_test/'; % folder where data is stored
paths.masterFolder=['data/' paths.nameFolder];
paths.nameSeries = 'fish1_scale1_96hpp_series';

opts = [];
opts.printHD = false; % resolution printing
opts.pxsize= 0.606; 
opts.xStart = -1000; opts.xEnd = 1000; %coordinates output image from scale centroid
opts.yStart = -1000; opts.yEnd = 1000;
opts.dx = 50; opts.dy= 50; % step size for velocity averaging
opts.plotRange = [-0.01 0.01]; % velocity display range
opts.dt = 3; % in hours
 
% calculate flows
load([paths.masterFolder '/objects/legend/' paths.nameSeries '.mat'])
opts.frameMin =1;
opts.frameMax = 6;
opts.dFrame = 3;
paths.outFolder =     [paths.masterFolder 'results/flowmap_' num2str(opts.dFrame.*opts.dt) 'h/'] ; % output folder quantifications
track=readtable([paths.masterFolder '/divided_eqIlastik/fish1_scale1_ch2_hypo_all_CSV-Table.csv']);
[divMsmooth_fish11_06mar19,coord06mar19] = flowsScaleReduced(paths,myScaleSeries,track,opts);
