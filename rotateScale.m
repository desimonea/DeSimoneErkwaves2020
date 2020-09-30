function [myScale] = rotateScale(myScale,refCh,chToRotate,ths,rects,zPos,paths)
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
%% It rotates each channel in the scale  

   inFolder = paths.inFolder;
   outFolder= paths.outFolder;
   objFolder= paths.objFolder;
   fieldName=paths.fieldName;

   basename=myScale.basename;

   ch=refCh;
   spath=[inFolder basename '_ch' num2str(refCh) '.tif'];
  
   mkdir(outFolder);
   
   myScale.(fieldName).folder=outFolder;
   
   options=[];
   options.overwrite='true';
   options.compress='lzw';
   
   %using loadtiff
   stack=loadtiff(spath);
          
   % find angles
   [outStack,ths,rects,zPos] = rotateStack(stack,ths,rects,zPos); 
   myScale.(fieldName).rotAngles=ths;
   myScale.(fieldName).cropRects=rects;
   myScale.(fieldName).zPos=zPos;
   
   if(isfield(paths,'saveName'))
        saveName=paths.saveName;
   else       
       saveName = basename;
   end
   
   if(isfield(paths,'objSaveName'))
       objSaveName=paths.objSaveName;
   else       
       objSaveName = basename;
   end
   
   
   
   if(any(chToRotate==refCh))
       
        %write channel
        wpath=[outFolder saveName '_ch' num2str(ch) '.tif']; 
        myScale.(fieldName).paths{ch}=wpath;
        
        saveastiff(outStack, wpath, options)
          
        %write projections
        wpath=[outFolder saveName '_ch' num2str(ch) 'maxproj.tif'];   
        saveastiff(max(outStack,[],3),wpath,options); 
   
   end
   
   for ch=chToRotate
      
       if(ch==refCh)
          continue; 
       end
       spath=[inFolder basename '_ch' num2str(ch) '.tif'];      
       infoStack=imfinfo(spath);
       
       %using loadtiff
       stack=loadtiff(spath);
   
       [outStack] = rotateStack(stack,ths,rects,zPos);

       %write channel
       wpath=[outFolder saveName '_ch' num2str(ch) '.tif']; 
       myScale.(fieldName).paths{ch}=wpath;
        
       saveastiff(outStack, wpath, options);
          
       %write projections
       wpath=[outFolder saveName '_ch' num2str(ch) 'maxproj.tif'];   
       saveastiff(max(outStack,[],3),wpath,options); 

   end

   %save mat file
   wpathmat=[objFolder objSaveName '.mat'];  
   save(wpathmat,'myScale');

end

