function [myScale,stitched] = stitchScale(paths,basename,opts)
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
%   OPTS fiels
%   opts.stitchMethod -> "3DsinglePlane" or "3Dproj"

% ---------

inFolder = paths.inFolder;
outFolder= paths.outFolder;
objFolder= paths.objFolder;

targetVox=opts.targetVox;
cropImage=opts.cropImage;
verbose=opts.verbose;
refCh=opts.refCh;

% zBest is the z to use for stitching
if(isfield(opts,'zBest'))
    zBest=opts.zBest;
else
    display('zBest not provided');
    zBest=[];
end

optsParse=[];
optsParse.timeStep=opts.timeStep;


if(isfield(opts,'folderName'))
   folderName = opts.folderName; 
else
   folderName = basename;
end


% parse metadata
if(isfield(opts,'inFolders'))
    if(opts.inFolders)
        myScale=parse_scale_folder([inFolder folderName filesep],basename,opts);
    else
        myScale=parse_scale_folder([inFolder filesep],basename,opts);        
    end
else
        myScale=parse_scale_folder([inFolder filesep],basename,opts);        
end

myScale.targetVox=targetVox;
myscale.refCh=refCh;

% assign bitDepth
info=imfinfo([myScale.folder myScale.dir(1).name]);
for k=1:numel(myScale.Position_list)
    myScale.metadata{k}.BitDepth=info.BitDepth;
end

% assign zBest if present
for k=1:numel(myScale.Position_list)
    myScale.metadata{k}.zBest=[];
    if(~isempty(zBest))
        idx=find(zBest(:,1)==myScale.Position_list(k));
        if(~isempty(idx))
            myScale.metadata{k}.zBest=zBest(idx,2);
        end
    end
end


%% we construct the stack

for n=1:numel(myScale.dir)
    
    if(all(isnan(myScale.Position_list)))
        pos=1;
    else
        pos=find(myScale.Position_list==myScale.dir(n).Position);
    end
    ch=myScale.dir(n).ch+1;
    z=myScale.dir(n).z+1;
    dir_paths{pos,ch}{z}=[myScale.folder myScale.dir(n).name];    
end

myScale.dir_paths=dir_paths;


%% stitch stacks
%stitched=stitchStacksCoord(myScale);


if(~isfield(opts,'stitchMethod'))

    stitched=stitchStacksCoordCoor2(myScale,opts);%refCh,targetVox,cropImage,verbose);
else

    if(strcmp(opts.stitchMethod,'3DsinglePlane'))
        stitched=stitchStacksCoordCoor2(myScale,opts);%refCh,targetVox,cropImage,verbose);
    elseif(strcmp(opts.stitchMethod,'3Dproj'))
        stitched=stitchStacksCoordCoor2Proj(myScale,opts);%refCh,targetVox,cropImage,verbose);
    end
end

%% write stitched stack

mkdir(outFolder);

% myScale.stScale=[];
% myScale.stScale.folder=stitched_folder;
% myScale.stScale.paths=[];

for ch=1:numel(stitched)

    wpath=[outFolder basename '_ch' num2str(ch) '.tif']; 
    myScale.stScale.paths{ch}=wpath;
   
%   % with imwrite
%   imwrite(stitched{ch}(:,:,1),wpath,'Compression',compression);
%   for z=2:size(stitched{ch},3)
%         imwrite(stitched{ch}(:,:,z),wpath, 'writemode', 'append','Compression',compression);
%   end       
%   %write projections
%   wpath=[stitched_folder basename '_ch' num2str(ch) 'maxproj.tif'];   
%   imwrite(max(stitched{ch},[],3),wpath,'Compression',compression); 

   %with saveastiff
   options=[];
   options.overwrite='true';
   options.compress='lzw';
   saveastiff(stitched{ch},wpath,options);
   % projections
   wpath=[outFolder basename '_ch' num2str(ch) 'maxproj.tif'];   
   saveastiff(max(stitched{ch},[],3),wpath,options);   
   
end

   %save mat file
   wpathmat=[objFolder basename '.mat'];  
   save(wpathmat,'myScale');


end

