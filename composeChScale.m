function [myScale] = composeChScale(myScale,chToCompose,rescaleFactors,nameComposed,paths)
%     Copyright (C) 2020  Alessandro De Simone
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
    if(strcmp(rescaleFactors,'adjust'))
        rescaleFactors=nan(numel(chToCompose),1);
    end

   inFolder = paths.inFolder;
   outFolder= paths.outFolder;
   objFolder= paths.objFolder;

   basename=myScale.basename;
   %rescaleFactors=rescaleFactors./sum(rescaleFactors); %normalization
   
   for i=1:numel(chToCompose)
        spath=[inFolder basename '_ch' num2str(chToCompose(i)) '.tif'];
        s=loadtiff(spath);
         
        if(any(isnan(rescaleFactors)))
            rescaleFactors(i)=255/double(max(s(:)));
        end
        
        if(i==1)
            snew=s.*rescaleFactors(1);
        else
           snew2=snew+s.*rescaleFactors(i);
           snew=snew2;
        end
   end
  
   options=[];
   options.overwrite='true';
   options.compress='lzw';
   
   % save new stack
   wpath=[outFolder basename '_ch' num2str(nameComposed) '.tif'];
   saveastiff(snew, wpath, options);
          
   %write projections
   wpath=[outFolder basename '_ch' num2str(nameComposed) 'maxproj.tif'];
   saveastiff(max(snew,[],3),wpath,options); 
   
end

