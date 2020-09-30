function [outStack,ths,rects,zPos] = rotateStack(stack,ths,rects,zPos,method)
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
% [outStack,plane,hyperrect] = rotateStack(stack,plane,hyperrect)
% 
% It rotates the stack to align the selected plane to the z-axis and
% optionally crops it
% 
% STACK is the stack, no surprise here
% PLANE is the 3D vector perpendicular to the plane
% HYPERRECT is the 3D rectangle to crop [xmin ymin zmin width height depth]


     if(isempty(ths))
         ths=[NaN NaN];
     end
         
     if(isempty(rects))
        rects=nan(4,4); 
     end
     
     if(strcmp(zPos,'all'))
        zPos=[1 size(stack,3)];
     end     

     % rotation around y
     [cropStack,th,rect]=imrotateStack(stack,0,rects(1,:));
     rects(1,:)=rect;   
     
     permStack1=permute(cropStack,[2 3 1]);
     
     [rotStack1,th,rect]=imrotateStack(permStack1,ths(1),rects(2,:));
     ths(1)=th;
     rects(2,:)=rect;  
     
     permStack2= permute(rotStack1,[3 1 2]);
     
     % rotation around x
     permStack2= permute(permStack2,[1 3 2]); 
     
     [rotStack2,th,rect]=imrotateStack(permStack2,ths(2),rects(3,:));
     ths(2)=th;
     rects(3,:)=rect;  
     
     outStack=permute(rotStack2,[1 3 2]); 
     
     %slice selection
     if(isempty(zPos))
        zPos=planeSelect(outStack,[0 255]);
     end
     
     if(numel(zPos)>1)
        outStack=outStack(:,:,min(zPos):max(zPos));
     elseif(zPos==0)
        disp('Taking the whole z stack');
     else
        disp('Not enough zPos elements, we will not crop anything!'); 
     end
     
     % final cropping
     [outStack,th,rect]=imrotateStack(outStack,0,rects(4,:));
     rects(4,:)=rect; 
     
     %imshow(max(outStack,[],3),[]);
     
end

function [rotStack,th,rect]=imrotateStack(stack,th,rect)

     h1=figure;

     %proj=mean(stack,3); %edited 16feb18
     %proj=mean(single(stack),3);
     proj=max(stack,[],3);
     proj=imadjust(proj,[0 0.2]);
    
     if(isnan(th))
         
%     % by clicking 
%       disp('Select points along the scale line');
%       [x,y,~]=impixel(proj);
%         
%       if(numel(x)>1)      
%             l=fit(x,y,'poly1');
%             th=atan(l.p1);
%       else
%             th=pi/2; 
%       end
%         
      % fitting the image
        [X,Y]=meshgrid(1:size(proj,2),1:size(proj,1));
        w=im2double(proj(:));
        idxs=find(w>0);
        data=[X(idxs) Y(idxs)];
        covv=weightedcov(data,w(idxs)');
        th=atan(1/(covv(2,1)./covv(4)));
        
         if(th>0)
             th= pi/2-th;
         else
             th=-(pi/2+th);
         end
         
         display(th);
         
     end
     
     projrot=imrotate(proj,-rad2deg(th));
     
     if(isnan(rect(1)))
        disp('Select the cropping rectangle');
        imshow(projrot,[]);
        [projrotcrop,rect]=imcrop;
        if(isempty(rect))
            rect=[0 0 0 0];
            projrotcrop=projrot;
        end
     elseif(all(rect==0))
        projrotcrop=projrot;
     else
        [projrotcrop,rect]=imcrop(projrot,rect);
     end
     
     rotStack=zeros([size(projrotcrop) size(stack,3)],'uint8');
     
     for n=1:size(stack,3)
        img=stack(:,:,n); 
        img2=imrotate(img,-rad2deg(th),'bilinear');
        if(all(rect==0))
            rotStack(:,:,n)=img2;
        else
            rotStack(:,:,n)=imcrop(img2,rect);
        end
     end
   
     close(h1)
     
end

