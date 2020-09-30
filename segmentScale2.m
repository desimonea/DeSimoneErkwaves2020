function [outStack,segOpts] = segmentScale2(s,segOpts,verbose)
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
% it gets a stack with nuclei and segments it out
% verbose = 0 (nothing), 1 (steps), 2 (every z stack and stops to show you)
% TO DO: implement segOpts


    if(verbose==2)
        stshow(s,[]);
    end

    if(verbose)
        disp('Thresholding and cleaning');
    end

    
    %T = adaptthresh(mean(s,3), 0.4,'NeighborhoodSize',[251 251]);

    levels=s(:);
    levels=levels(~(levels==0));
    thresh=graythresh(levels);
    
    %%
    
    for z=1:size(s,3)
        %T = adaptthresh(max(s(:,:,z:min(z+10,size(s,3))),[],3), 0.4,'NeighborhoodSize',[251 251]);
                
        if(verbose==2)
            disp(['Processing ' num2str(z) '/' num2str(size(s,3))]);
        end
        
        img=s(:,:,z); 
        %img=imbinarize(im2uint8(img),T);
        img=imbinarize(im2uint8(img),thresh*0.4); 
         s(:,:,z)=img;
    end

    
    if(verbose==2)
        stshow(s,[]);
    end

    for z=1:size(s,3)
        img=s(:,:,z); 
        img=bwmorph(img,'clean');
        s(:,:,z)=img;
    end
    
    s=permute(s,[1 3 2]);
    for z=1:size(s,3)
        img=s(:,:,z); 
        img=bwmorph(img,'clean');
        s(:,:,z)=img;
    end
    s=permute(s,[1 3 2]);
    
    
    s=permute(s,[3 2 1]);
    for z=1:size(s,3)
        img=s(:,:,z); 
        img=bwmorph(img,'clean');
        s(:,:,z)=img;
    end
    s=permute(s,[3 2 1]);
    
    for z=1:size(s,3)
        img=s(:,:,z); 
        img=bwmorph(img,'majority',3);
        seg3(:,:,z)=logical(img); 
    end

    if(verbose==2)
        stshow(seg3,[]);
    end
    

    %%
    
    if(verbose)
        disp('Opening a bit the image');
    end

    % ABw: fare bwareaopne instead of opening
    
    %seg4=imopen(seg3,strel('sphere',3));
    seg4=imopen(seg3,strel('disk',3));
    seg4=imopen(seg4,strel('disk',3));
    seg4=imopen(seg4,strel('disk',3));
    
    seg4=permute(seg4,[1 3 2]);
    seg4=imopen(seg4,strel('disk',3));
    seg4=permute(seg4,[1 3 2]);
    
    clear('seg3');
    
    if(verbose==2)
        stshow(seg4,[]);
    end

   %%
     
    if(verbose)
        disp('Dilating and closing');
    end

    seg5=imdilate(seg4,strel('disk',5));
    clear('seg4');
 
    seg5=imdilate(seg5,strel('disk',5));
    seg5=imclose(seg5,strel('disk',20));

    
    if(verbose==2)
        stshow(seg5,[]);
    end
    
    %%
    if(verbose)
        disp('Filling holes');
    end

    seg5neg=~seg5;
    seg6neg=bwareaopen(seg5neg,10000,8);
    
    seg6=~seg6neg;
    
    if(verbose==2)
        stshow(seg6,[]);
    end
 

   outStack=seg6;
%%


    %%
    if(verbose==2)
        masked=ones(size(seg6))*30;masked(seg6)=s(seg6);
        stshow(masked,[0 30]);

        notmasked=s;notmasked(seg6)=0;
        imshow(max(notmasked,[],3),[0 30]);
    end

end
