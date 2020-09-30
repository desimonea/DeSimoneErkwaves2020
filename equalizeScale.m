function [myScale] = equalizeScale(myScale,chToEq,paths)
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
   inFolder=paths.inFolder;
   outFolder=paths.outFolder;
   objFolder=paths.objFolder;

   basename=myScale.basename;
   
   options=[];
   options.overwrite='true';
   options.compress='lzw';
   
   for ch=chToEq
      
       spath=[inFolder basename '_ch' num2str(ch) '.tif'];
      
       %using loadtiff
       stack=loadtiff(spath);
       
       stack=equalizeStack(stack);

       %write channel
       wpath=[outFolder basename '_ch' num2str(ch) '.tif']; 
        
       saveastiff(stack, wpath, options);
          
       %write projections
       wpath=[outFolder basename '_ch' num2str(ch) 'maxproj.tif'];   
       saveastiff(max(stack,[],3),wpath,options); 

   end
   
   %save mat file
   wpathmat=[objFolder basename '.mat'];  
   save(wpathmat,'myScale');

end

