function [myScale] = flattenLayersScaleSegm(myScale,toflattenCh,paths,opts)
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

   options=[];
   options.overwrite='true';
   options.compress='lzw';
   
   inFolder = paths.inFolder;
   outFolder= paths.outFolder;
   objFolder= paths.objFolder;
   
   basename=myScale.basename;
   
   if(strcmp(opts.method,'fit'))
     fitFieldName=paths.fitFieldName;
   end
       
   roiFolder= paths.roiFolder;
   spath=[roiFolder basename '_ROI.tif'];
   
   sroi=loadtiff(spath);
   sroi=logical(sroi);     
   
   roughnessFit=assignOpt(opts,'roughnessFit',0.0001); %0.0001
   zshift=assignOpt(opts,'zShift',10); 
   
   
   if(strcmp(opts.method,'segmentation'))

       
       zstarts=nan([size(sroi,1) size(sroi,2)]);
       zends=zstarts;  
   
       mask=max(sroi,[],3);
   
       for i=1:size(sroi,1)
            for j=1:size(sroi,2)
            
                if(mask(i,j))
   
                    zstart=1;
                    zend=min(find([sroi(i,j,:)==1],1,'last')+zshift,size(sroi,3)); 
            
                    if(isempty(zend))
                        zstart=1;
                        zend=1;
                    end
            
                    zstarts(i,j)=zstart;
                    zends(i,j)=zend;    
                end
            end
        end
       
   elseif(strcmp(opts.method,'fit'))
       
       %% calculate zEnd
         
       %find points where there are nuclei
       nucleiReg=max(sroi,[],3);
       [XX,YY]=meshgrid(1:size(nucleiReg,2),1:size(nucleiReg,1));
       dotList=[YY(nucleiReg(:)) XX(nucleiReg(:))];

       zEndData=nan(size(dotList,1),1);
   
       zIncr=repmat(reshape(1:size(sroi,3),[1 1 size(sroi,3)]),[size(nucleiReg,1) size(nucleiReg,2)]);    
       zMax=max(zIncr.*sroi,[],3);
       zEndData=zMax(nucleiReg(:));
       zEndData=zEndData+zshift;
          
        %fitHypo=fit([dotList(:,1) dotList(:,2)],zEndData(:),'poly55'); 
        %%

         F = scatteredInterpolant([dotList(:,1) dotList(:,2)],zEndData(:),'natural','nearest'); 
%          
%          isScaleY=logical(sum(nucleiReg,2));
%          isScaleX=logical(sum(nucleiReg,1));
%          
%          startY=find(isScaleY,1,'first');
%          endY=find(isScaleY,1,'last');
%          
%          startX=find(isScaleX,1,'first');
%          endX=find(isScaleX,1,'last');
%          
% 
         startX=1; endX=size(nucleiReg,2);
         startY=1; endY=size(nucleiReg,1);

         dstep=1;
         [XXrough,YYrough]=meshgrid(startX:dstep:endX,startY:dstep:endY);
         interpZ = F(YYrough,XXrough);
   
         xrough={startY:dstep:endY,startX:dstep:endX};
         fitHypo =  csaps(xrough,interpZ,roughnessFit); %,[],x);   0.00001
%    
         mask=ones([size(sroi,1) size(sroi,2)],'logical'); %we have a value for all points
         zstarts=ones(size(YY));
         zends= nan(size(YY));  
% 
         x={startY:endY,startX:endX};
         zendsCrop= fnval(fitHypo,x); 
         
         zends(startY:endY,startX:endX)=zendsCrop;
%       
         zends=round(zends);
         mask(zends(:)> size(sroi,3))= 0; %if larger the z-stack we set it to nan
         mask(zends(:)< 1)= 0; %if below if we set it to nan
         mask(isnan(zends))=0;
         %%
%          figure
%          for j=1:10:size(zends,2)
%             disp(j);
%             idxs=dotList(:,2)==j;
%             plot(dotList(idxs,1),zEndData(idxs),'ok'); hold on
%             %plot((1:size(interpZ,1)).*dstep,interpZ(:,ceil(j./dstep)),'ob'); hold on
%             plot(zends(:,j),'r-'); hold off; 
%             ylim([0 150]);
%             waitforbuttonpress;
%             drawnow;
%             hold off;
%          end         

   else
       error('Please provide a valid method: segmentation or fit'); 
   end
  
   for ch=toflattenCh
       
    %%  read stack      
         
        spath=[inFolder basename '_ch' num2str(ch) '.tif'];
        s=loadtiff(spath);
  
   %%   calculate stack
  
        snew=zeros(size(s),'uint8');
        for i=1:size(snew,1)
            for j=1:size(snew,2)      
                if(mask(i,j))
                    zend=zends(i,j);
                    zstart=zstarts(i,j);
                    snew(i,j,1:1+zend-zstart)= s(i,j,zend:-1:zstart);
                    %snew(i,j,1:1+zend-zstart)= s(i,j,zstart:zend);
                end        
            end
        end 
        
       %% write stack
       
       wpath=[outFolder myScale.basename '_ch' num2str(ch) '.tif'];
       saveastiff(snew,wpath,options);
       
       %% write zend -> z' = zend - z +1 
       
       wpath=[outFolder myScale.basename '_zend.tif'];
       saveastiff(uint16(zends),wpath,options);
  
   end
   
   %% flatten and write ROI
           
       snew=zeros(size(sroi),'uint8');

       for i=1:size(snew,1)
        for j=1:size(snew,2)      
            if(mask(i,j))
                        zend=zends(i,j);
                        zstart=zstarts(i,j);
                        snew(i,j,1:1+zend-zstart)= sroi(i,j,zend:-1:zstart);
            end        
         end
        end
        wpath=[outFolder myScale.basename '_ROI.tif'];
        saveastiff(snew,wpath,options);

  
   %% save mat file
   wpathmat=[objFolder basename '.mat'];  
   save(wpathmat,'myScale');
   
end 

function [opt]=assignOpt(opts,nameOpt,defaultValue)

    if(isfield(opts,nameOpt))
        opt=opts.(nameOpt);
    else
        opt=defaultValue;
    end
end