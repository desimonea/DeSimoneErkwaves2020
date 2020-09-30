function [myScale] = cleanupScale(myScale,refCh,chToClean,roi,paths,opts)
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
% [myScale] = cleanupScale(myScale,refCh,chToClean,roi,paths,opts)
% 
% myScale is (obviously) the myScale object
% refCh is the reference channel to use for cleanup
% chToClean chooses what channels to clean
% roi is the roi to use for cleanup (no re-calculation of the roi)

%%

   inFolder = paths.inFolder;
   outFolder= paths.outFolder;
   objFolder= paths.objFolder;

   basename=myScale.basename;
   
   mkdir(outFolder);
   
   ch=refCh;
   
   % parsing opts field
   if(isfield(opts,'type'))
       type=opts.type;
   else
       type='manual';
   end
       
   if(isfield(opts,'segOpts'))
        segOpts=opts.segOpts;
   else
         segOpts=[];
   end
   
   if(isfield(opts,'verbose'))
       verbose=opts.verbose;
   else
       verbose=[]; 
   end
    
   % writing options
   options=[];
   options.overwrite='true';
   options.compress='lzw';
   
   % read old roi
   if(isfield(opts,'readROI'))
       if(~isempty(opts.readROI))
            roipath=[paths.(opts.readROI) basename '_ROI'  '.tif'];
            if(exist(roipath,'file'))
                oldroi= loadtiff(roipath);
                oldroi=logical(oldroi);
            else
                oldroi=[];
            end
       end
   else
        oldroi=[];
   end
      
   
   if(isempty(roi))
       if(ch==0) % sum all channels
           for c=1:myScale.nch

               spath=[inFolder basename '_ch' num2str(c) '.tif'];

                %using loadtiff
                stack_ch=loadtiff(spath);  
                if(c==1)
                    stack=stack_ch;
                else
                    stack=max(stack,stack_ch);
                end
           end
           
           if(~isempty(oldroi))
                stack(~oldroi)=0;
           end
           
           if(strcmp(type,'manual'))
               [stack,roi]=cleanupStack(stack,roi,opts); 
           elseif(strcmp(type,'KTR'))
               [stack,roi]=cleanupScaleSegKTR(stack,segOpts,verbose); 
           end
       else

           spath=[inFolder basename '_ch' num2str(ch) '.tif'];

           myScale.cleanedScale=struct;
           myScale.cleanedScale.folder=outFolder;

           %using loadtiff
           stack=loadtiff(spath);      
           
           if(~isempty(oldroi))
                stack(~oldroi)=0;
           end
           
           
           if(strcmp(type,'manual'))
               [stack,roi]=cleanupStack(stack,oldroi,opts); 
           elseif(strcmp(type,'KTR'))
               [stack,roi]=cleanupScaleSegKTR(stack,segOpts,verbose); 
           end

           if(any(chToClean==refCh))

               %write channel
               wpath=[outFolder basename '_ch' num2str(ch) '.tif']; 
               myScale.cleanedScale.paths{ch}=wpath;

               saveastiff(stack, wpath, options);

               %write projections
               wpath=[outFolder basename '_ch' num2str(ch) 'maxproj.tif'];   
               saveastiff(max(stack,[],3),wpath,options); 
           end
       end
   end

   newroi=roi;
   clear('roi'); clear('oldroi');
   
   %write roi
   wpath=[outFolder basename '_ROI'  '.tif']; 
   myScale.cleanedScale.roipath=wpath;
   saveastiff(im2uint8(newroi), wpath, options);
   
   for ch=chToClean
      
%        if(ch==refCh)
%           continue; 
%        end
       
       spath=[inFolder basename '_ch' num2str(ch) '.tif'];
       infoStack=imfinfo(spath);
       
       %using loadtiff
       stack=loadtiff(spath);
       stack(~newroi)=0;   

       %write channel
       wpath=[outFolder basename '_ch' num2str(ch) '.tif']; 
       myScale.cleanedScale.paths{ch}=wpath;
        
       saveastiff(stack, wpath, options);
          
       %write projections
       wpath=[outFolder basename '_ch' num2str(ch) 'maxproj.tif'];   
       saveastiff(max(stack,[],3),wpath,options); 

   end
   
   %save mat file
   wpathmat=[objFolder basename '.mat'];  
   save(wpathmat,'myScale');

end

