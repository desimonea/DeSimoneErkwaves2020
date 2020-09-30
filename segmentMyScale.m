function [myScale] = segmentMyScale(myScale,refCh,paths,segOpts,verbose)
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
% it gets myScale and runs segmentScale passing a stack to it

   inFolder = paths.inFolder;
   outFolder= paths.outFolder;
   objFolder= paths.objFolder;

   pxArea=myScale.targetVox(1)*myScale.targetVox(2);
   voxVolume=prod(myScale.targetVox);

   basename=myScale.basename;

   ch=refCh;
   spath=[inFolder basename '_ch' num2str(ch) '.tif'];
  
   mkdir(outFolder);
   myScale.segScale=struct;
   myScale.segScale.folder=outFolder;
   
   options=[];
   options.overwrite='true';
   options.compress='lzw';
   
   %using loadtiff
   stack=loadtiff(spath);
          
   % segmentation
   [mask,segOpts] = segmentScale2(stack,segOpts,verbose); 
   
   proj=max(mask,[],3);
   
   area=sum(proj(:));
   volume=sum(mask(:));
  
   myScale.segScale.area=area*pxArea;
   myScale.segScale.volume=volume*voxVolume;
   myScale.segScale.segOpts=segOpts;
   
   %write segmentation  (ROI)
   wpath=[outFolder basename '_ROI.tif']; 
   myScale.segScale.paths_ROI=wpath;
   saveastiff(im2uint8(mask), wpath, options);
   
   %write segmented
   wpath=[outFolder basename '_ch' num2str(refCh) '_seg.tif']; 
   myScale.segScale.paths_segmented=wpath;
   segStack=im2uint8(stack); segStack(~logical(mask))=255;
   saveastiff(segStack, wpath, options);
   clear('segStack');
   
   %write not segmented
   wpath=[outFolder basename '_ch' num2str(refCh) '_notseg.tif']; 
   myScale.segScale.paths_notSegmented=wpath;
   notsegStack=im2uint8(stack); notsegStack(logical(mask))=0;
   saveastiff(notsegStack, wpath, options);
   clear('notsegStack');
   
   %save mat file
   wpathmat=[objFolder basename '.mat'];  
   save(wpathmat,'myScale');

end