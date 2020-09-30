function [stitched,offsets] = stitchStacksCoordCoor2(mys,opts)%refCh,targetVox,cropImage,verbose)
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

vox1=mys.metadata{1}.voxel;
%targetVox=vox1;

bitDepth=mys.metadata{1}.BitDepth;

% read options
refCh=opts.refCh;
targetVox=opts.targetVox;
cropImage=opts.cropImage;
verbose=opts.verbose;

if(isfield(opts,'resizeFactors'))
    resizeFactors=opts.resizeFactors;
else
    resizeFactors=[1 1];
end


if(isfield(opts,'overlaps'))
    overlaps=opts.overlaps;
else
    overlaps=[];
end

% build stack1
k=1;

if(size(mys.dir_paths{k,refCh},2)==1)
        stack1=im2uint8(loadtiff(mys.dir_paths{k,refCh}{1}));
else
    stack1=zeros([mys.metadata{k}.sizeY mys.metadata{k}.sizeX mys.nz(k)],'uint8');
    for n=1:size(mys.dir_paths{k,refCh},2)
        img_temp=loadtiff(mys.dir_paths{k,refCh}{n});
        stack1(:,:,n)=im2uint8(img_temp);
    end
end

% we rescale the stack
stack1=rescaleStack(stack1,vox1,targetVox);

%the first image fills the entire stack
subImages=[1 1 1 size(stack1)];

% we convert positions in pixels - posZ now is the first stack
X1=[-mys.metadata{1}.PosY/targetVox(1) mys.metadata{1}.PosX/targetVox(2) mys.metadata{1}.PosZ/targetVox(3)];

stitched{refCh}=stack1;

for ch=1:size(mys.dir_paths,2)
    if(ch==refCh)
        continue
    else
        
        %build stack1
        if(size(mys.dir_paths{1,ch},2)==1)
             stack_ch=im2uint8(loadtiff(mys.dir_paths{1,ch}{1}));
        else
             stack_ch=zeros([mys.metadata{k}.sizeY mys.metadata{k}.sizeX mys.nz(k)],'uint8');
             for n=1:size(mys.dir_paths{1,refCh},2)
                img_temp=loadtiff(mys.dir_paths{1,ch}{n});
                stack_ch(:,:,n)=im2uint8(img_temp);
             end
        end
        % we rescale the stack
        stack_ch=rescaleStack(stack_ch,vox1,targetVox);
        stitched{ch}=stack_ch;
    end
end

offsets=[0 0 0];

for k=2:size(mys.dir_paths,1)

    display(['-----------']);
    display(['Now stitching Position ' num2str(mys.Position_list(k))]);
    
    stack1=stitched{refCh};
    
    vox2=mys.metadata{k}.voxel;
     
    % build stack2
    if(size(mys.dir_paths{k,refCh},2)==1)
        stack2=im2uint8(loadtiff(mys.dir_paths{k,refCh}{1}));
    else
        stack2=zeros([mys.metadata{k}.sizeY mys.metadata{k}.sizeX mys.nz(k)],'uint8');
        for n=1:size(mys.dir_paths{k,refCh},2)
            img_temp=loadtiff(mys.dir_paths{k,refCh}{n});
            stack2(:,:,n)=im2uint8(img_temp);
        end
    end

    % we rescale the stack
    stack2=rescaleStack(stack2,vox2,targetVox);
   
    if(~mys.metadata{k}.voxel==vox1)
        error('Images do not have the same voxel size!');
    end

    X2=[-mys.metadata{k}.PosY/targetVox(1) mys.metadata{k}.PosX/targetVox(2) mys.metadata{k}.PosZ/targetVox(3)];
    size1=size(stack1);
    size2=size(stack2);
   
    %offsets using coordinates
    %yOffset=round(X2(1)-X1(1));
    %xOffset=round(X2(2)-X1(2));
    %Offset=round(X2(3)-X1(3));
    
    % we transform in discrete coordinates
    X1=floor(X1);X2=floor(X2);
    
    if(~isempty(overlaps))
        if(~isempty(overlaps{k}))
            XStart1=overlaps{k}.XStart1;
            XEnd1=overlaps{k}.XEnd1;
            XStart2=overlaps{k}.XStart2;
            XEnd2=overlaps{k}.XEnd2;
            % TO DO: modify so that you provide XStart2 and XEnd2 and the
            % subImage to use
        else
            % we calulate the overlap between stack2 and a subImage in stack1
            [XStart1,XEnd1,XStart2,XEnd2,usedSub]=findOverlap(X1,X2,subImages,size1,size2,cropImage);
            disp(['Subimage used: Position' num2str(mys.Position_list(usedSub))]);
        end
    else
            if(opts.allSubs==0 & isempty(mys.metadata{k}.zBest))
                % we calulate the overlap between stack2 and a subImage in stack1
                [XStart1,XEnd1,XStart2,XEnd2,usedSub]=findOverlap(X1,X2,subImages,size1,size2,cropImage);
                disp(['Subimage used: Position' num2str(mys.Position_list(usedSub))]);

                if(isempty(XStart1))
                    disp('There is no overlap, so I check the overlap with all subs');
                end
            else
                XStart1=[];XEnd1=[];
                XStart2=[];XEnd2=[];
            end
            
            if(~isempty(mys.metadata{k}.zBest))
                   disp('You have provided the best z, so I check the overlap with all subs');
            end
        
    end
    
    if(~isempty(XStart1))
    
        yStart1=XStart1(1);    xStart1=XStart1(2);    zStart1=XStart1(3);
        yStart2=XStart2(1);    xStart2=XStart2(2);    zStart2=XStart2(3);

        yEnd1=XEnd1(1);    xEnd1=XEnd1(2);    zEnd1=XEnd1(3);
        yEnd2=XEnd2(1);    xEnd2=XEnd2(2);    zEnd2=XEnd2(3);

        %we carve out the templates
        template1=stack1(yStart1:yEnd1,xStart1:xEnd1,zStart1:zEnd1);
        template2=stack2(yStart2:yEnd2,xStart2:xEnd2,zStart2:zEnd2);
        
        %first round of alignment using a coarse-grained image

        template1small=imresize(template1,resizeFactors(1));
        template2small=imresize(template2,resizeFactors(1));
        [relOffset,corrValue,zbest2]=correlateStacks(template1small,template2small,[]);
        relOffset=[relOffset(1)/resizeFactors(1) relOffset(2)/resizeFactors(1) relOffset(3)];
        
        % second round using only zbest1 and zbest2 and using the finer
        % movements
        zbest1=zbest2+relOffset(3);
        corrhere=normxcorr2(template2(:,:,zbest2),template1(:,:,zbest1));
        [ypeak, xpeak] = find(corrhere==max(corrhere(:)),1,'first');
        relOffset(1) = ypeak-size(template2,1);
        relOffset(2) = xpeak-size(template2,2);
        corrValue=[corrhere(ypeak,xpeak)];

    else
        
       corrValue=-Inf;
       
       for sub=1:size(subImages,1) 
          
            subHereRel=subImages(sub,:);
            
            XStart1= subHereRel(1:3);
            XEnd1  = subHereRel(4:6);
            
            XStart2= [1 1 1];
            XEnd2  = size2;
             
            yStart1=XStart1(1);    xStart1=XStart1(2);    zStart1=XStart1(3);
            yStart2=XStart2(1);    xStart2=XStart2(2);    zStart2=XStart2(3);

            yEnd1=XEnd1(1);    xEnd1=XEnd1(2);    zEnd1=XEnd1(3);
            yEnd2=XEnd2(1);    xEnd2=XEnd2(2);    zEnd2=XEnd2(3);

            %we carve out the templates
            template1=stack1(yStart1:yEnd1,xStart1:xEnd1,zStart1:zEnd1);
            template2=stack2(yStart2:yEnd2,xStart2:xEnd2,zStart2:zEnd2);
            
            template1small=imresize(template1,resizeFactors(2));
            template2small=imresize(template2,resizeFactors(2));
            [relOffsetHere,corrValueHere,zbest2]=correlateStacks(template1small,template2small,mys.metadata{k}.zBest);

            relOffsetHere=[relOffsetHere(1)/resizeFactors(2) relOffsetHere(2)/resizeFactors(2) relOffsetHere(3)]; 
            
            zbest1=zbest2+relOffsetHere(3);
            corrhere=normxcorr2(template2(:,:,zbest2),template1(:,:,zbest1));
            [ypeak, xpeak] = find(corrhere==max(corrhere(:)),1,'first');
            relOffsetHere(1) = ypeak-size(template2,1);
            relOffsetHere(2) = xpeak-size(template2,2);
            corrValueHere=[corrhere(ypeak,xpeak)];
            
            sub
            relOffsetHere
            corrValueHere
            
            if(corrValueHere > corrValue)
                chosenSub = sub;
                relOffset = relOffsetHere;
                corrValue =corrValueHere;
            end           
       end 
       
       disp(['Subimage used: Position' num2str(mys.Position_list(chosenSub))]);
       subHereRel=subImages(chosenSub,:);        
       XStart1= subHereRel(1:3); XEnd1  = subHereRel(4:6);     
       XStart2= [1 1 1]; XEnd2  = size2;
             
       yStart1=XStart1(1);    xStart1=XStart1(2);    zStart1=XStart1(3);
       yStart2=XStart2(1);    xStart2=XStart2(2);    zStart2=XStart2(3);

       yEnd1=XEnd1(1);    xEnd1=XEnd1(2);    zEnd1=XEnd1(3);
       yEnd2=XEnd2(1);    xEnd2=XEnd2(2);    zEnd2=XEnd2(3);
       
    end
    clear('template1');clear('template2');
    
    disp('Stitching ...');
    
    % we calculate the offset with respect to the total stack1
    yOffset=relOffset(1)+yStart1-1-(yStart2-1);
    xOffset=relOffset(2)+xStart1-1-(xStart2-1);
    zOffset=relOffset(3)+zStart1-1-(zStart2-1);
    
    [newImage,X1new,subImagesNew]=stitch(stitched{refCh},stack2,[yOffset,xOffset,zOffset],X1,subImages);
    stitched{refCh}=newImage;
    
    % time to stitch images (finally)  
    for ch=1:size(mys.dir_paths,2)
        if(ch==refCh)
            continue
        end
        
        if(size(mys.dir_paths{k,refCh})==1)
            stack2stitch=im2uint8(loadtiff(mys.dir_paths{k,ch}{1}));
        else
            stack2stitch=zeros([mys.metadata{k}.sizeY mys.metadata{k}.sizeX mys.nz(k)],'uint8');
            for n=1:size(mys.dir_paths{k,refCh},2)
                img_temp=loadtiff(mys.dir_paths{k,ch}{n});
                stack2stitch(:,:,n)=im2uint8(img_temp);
            end
        end
        
        stack2stitch=rescaleStack(stack2stitch,vox2,targetVox);
        [newImage,X1new,subImagesNew]=stitch(stitched{ch},stack2stitch,[yOffset,xOffset,zOffset],X1,subImages);
        stitched{ch}=newImage;
    end

    % update coordinates and subImages
    offsets(k,:)=[yOffset,xOffset,zOffset];
    X1=X1new;
    subImages= subImagesNew;
    
    if(verbose==1)
      h1=figure;
      imshow(max(stitched{refCh},[],3),[]);
      h2=figure;
      imshow(squeeze(max(stitched{refCh},[],2)),[]);
      %waitforbuttonpress;
      %close(h1);close(h2);
    end
end
    
if(verbose==2)
      h1=figure;
      imshow(max(stitched{refCh},[],3),[]);
      h2=figure;
      imshow(squeeze(max(stitched{refCh},[],2)),[]);
      %waitforbuttonpress;
      %close(h1);close(h2);
end

end

function [newStack,newCoord,subImages]=stitch(stack1,stack2,offsets,coord,subImages)

    yOffset=offsets(1);xOffset=offsets(2);zOffset=offsets(3);
     
    size1=size(stack1);
    size2=size(stack2);
    
    yStart=min(1,1+yOffset);
    xStart=min(1,1+xOffset);
    zStart=min(1,1+zOffset);
    
    yEnd=max(size1(1),size2(1)+yOffset);
    xEnd=max(size1(2),size2(2)+xOffset); 
    zEnd=max(size1(3),size2(3)+zOffset);
    
    newCoord=coord+[yStart-1,xStart-1,zStart-1];
    sizeNew=[yEnd-yStart+1,xEnd-xStart+1,zEnd-zStart+1];   
    
    % subImages - [yStartSubIm1 xStartSubIm1 zStartSubIm1 yEndSubIm1 xEndSubIm1 zEndSubIm1; ...] 
    % we move subImages of stack1
    subImages(:,1)=subImages(:,1)+1-yStart;subImages(:,4)=subImages(:,4)+1-yStart;
    subImages(:,2)=subImages(:,2)+1-xStart;subImages(:,5)=subImages(:,5)+1-xStart;
    subImages(:,3)=subImages(:,3)+1-zStart;subImages(:,6)=subImages(:,6)+1-zStart;
    
    %we create a logical image to know where the subimages are
    isNew=ones(sizeNew,'logical');
    
    for sub=1:size(subImages,1)
        isNew(subImages(sub,1):subImages(sub,4),subImages(sub,2):subImages(sub,5),subImages(sub,3):subImages(sub,6))=0;
    end

    newStack=zeros(sizeNew,'uint8');
    newStack12=zeros(sizeNew,'uint8');
    
    %we position stack2
    newStack(1+yOffset-yStart+1:size2(1)+yOffset-yStart+1,1+xOffset-xStart+1:size2(2)+xOffset-xStart+1,1+zOffset-zStart+1:size2(3)+zOffset-zStart+1)=stack2;
    % we position stack1 - it will make some black regions
    newStack(1-yStart+1:size1(1)-yStart+1,1-xStart+1:size1(2)-xStart+1,1-zStart+1:size1(3)-zStart+1)=stack1;
   
    % we position stack1 - it has the problem of bleaching  
    newStack12(1-yStart+1:size1(1)-yStart+1,1-xStart+1:size1(2)-xStart+1,1-zStart+1:size1(3)-zStart+1)=stack1;
    %we position stack2
    newStack12(1+yOffset-yStart+1:size2(1)+yOffset-yStart+1,1+xOffset-xStart+1:size2(2)+xOffset-xStart+1,1+zOffset-zStart+1:size2(3)+zOffset-zStart+1)=stack2;
   
    % we use newStack, but substitute the new regions with newStack12
    newStack(isNew)=newStack12(isNew);
    clear('newStack12'); clear('isNew');

    %we attach the subImage of stack2
    newSubImage=[1+yOffset-yStart+1, 1+xOffset-xStart+1,1+zOffset-zStart+1];
    newSubImage=[newSubImage, size2(1)+yOffset-yStart+1,size2(2)+xOffset-xStart+1, size2(3)+zOffset-zStart+1];
    subImages=vertcat(subImages,newSubImage);
    
end

function [offset,corrValue,zbest2]=correlateStacks(stack1,stack2,zbest2)
   
    if(isempty(zbest2))
        
%         % we need to choose which plane to use to compare with stack1
%         % we believe that stack2 is overlapping with stack 1
% 
%         % aligns the projections
%         s1proj=max(stack1,[],3);
%         s2proj=max(stack2,[],3);
% 
%         corrhere=normxcorr2(s2proj,s1proj);
%         [ypeak, xpeak] = find(corrhere==max(corrhere(:)),1,'first');
%         xOffset = xpeak-size(s2proj,2);
%         xOverStart= max(1-xOffset,1); % this is the ref system of stack2
%         xOverEnd  = min(size(s1proj,2)-xOffset,size(s2proj,2));
% 
%         xbest = round((xOverStart+xOverEnd)./2);
% 
%         % this method aligns in z
%         s1_xz=squeeze(stack1(:,xbest+xOffset,:));
%         s2_xz=squeeze(stack2(:,xbest,:));    
% 
%         if(size(s1_xz,2)>size(s2_xz,2))
%             corrhere=normxcorr2(s2_xz,s1_xz);
%             [ypeak, zpeak] = find(corrhere==max(corrhere(:)),1,'first');
%             zOffset = zpeak-size(s2_xz,2);
%         else
%             corrhere=normxcorr2(s1_xz,s2_xz);
%             [ypeak, zpeak] = find(corrhere==max(corrhere(:)),1,'first');
%             zOffset = -zpeak+size(s1_xz,2);
%         end
% 
%         zOverStart= max(1-zOffset,1); % this is the ref system of stack2
%         zOverEnd  = min(size(s1_xz,2)-zOffset,size(s2_xz,2));
%      
%         % choose the brightest plane in the overlap
%         [m,zbest2]=max(squeeze(mean(mean(stack2(:,:,zOverStart:zOverEnd),1),2)));
%         zbest2 =zbest2+zOverStart-1;

          [m,zbest2]=max(squeeze(mean(mean(stack2(:,:,:),1),2)));
          zbest2 = floor(size(stack2,3)./2);

    end
     
    best2=stack2(:,:,zbest2);
    
    if(std(im2double(best2(:)))==0)
       error('The selected portion of stack2 is completely flat!');
    end
    
   for z=1:size(stack1,3)
        if(std(im2double(stack1(:,:,z)))==0)
            display(['We skip z= ' num2str(z) 'because it is flat']);
            continue;
        end
        
        corrhere=normxcorr2(best2,stack1(:,:,z));
        [ypeak, xpeak] = find(corrhere==max(corrhere(:)),1,'first');
        yOffset = ypeak-size(best2,1);
        xOffset = xpeak-size(best2,2);
        corrValueZ(z,:)=[yOffset xOffset corrhere(ypeak,xpeak)];
        
    end
    
    [corrValue,zpeak]=max(corrValueZ(:,3));
    zOffset=zpeak-1-zbest2+1;
    yOffset=corrValueZ(zpeak,1); xOffset=corrValueZ(zpeak,2);
    
    offset=[yOffset,xOffset,zOffset];

end

function [XStart1,XEnd1,XStart2,XEnd2,usedSub]=findOverlap(X1,X2,subImages,size1,size2,cropImage)
% X1, X2 are the absolute coordinates of the two stacks
% subImages contains the relative coordinates of the subImages

    XStart1=[];
    XEnd1=[];

    XStart2=[];
    XEnd2=[];

    usedSub=[];

    for sub=size(subImages,1):-1:1

        subHereRel=subImages(sub,:);
       
        %we calculate absolute coordinates subImage   
        subHereAbs(1:3)=subHereRel(1:3)-1+X1(1:3);
        subHereAbs(4:6)=subHereRel(4:6)-1+X1(1:3);

        % lets calculate overlap with subImages in absolute coordinates
        yStart=max(subHereAbs(1),X2(1));
        xStart=max(subHereAbs(2),X2(2));
        zStart=max(subHereAbs(3),X2(3));
    
        yEnd=min(subHereAbs(4),X2(1)+size2(1)-1);
        xEnd=min(subHereAbs(5),X2(2)+size2(2)-1);
        zEnd=min(subHereAbs(6),X2(3)+size2(3)-1);    
        
        if(((yEnd-yStart+1)>0)&((xEnd-xStart+1)>0)&((zEnd-zStart+1)>0))
            
            if(isempty(XStart1)|(yEnd-yStart+1)*(xEnd-xStart+1)>(XEnd2(1)-XStart2(1)+1)*(XEnd2(2)-XStart2(2)+1))
  
                    %we convert in relative coordinates in stack1 and stack2
                    yStart1=yStart-X1(1)+1;xStart1=xStart-X1(2)+1;zStart1=zStart-X1(3)+1;
                    yStart2=yStart-X2(1)+1;xStart2=xStart-X2(2)+1;zStart2=zStart-X2(3)+1;

                    yEnd1=yEnd-X1(1)+1;xEnd1=xEnd-X1(2)+1;zEnd1=zEnd-X1(3)+1;
                    yEnd2=yEnd-X2(1)+1;xEnd2=xEnd-X2(2)+1;zEnd2=zEnd-X2(3)+1; 
                    
                    XStart1=[yStart1,xStart1,zStart1];
                    XEnd1=[yEnd1,xEnd1,zEnd1];

                    XStart2=[yStart2,xStart2,zStart2];
                    XEnd2=[yEnd2,xEnd2,zEnd2];
                    
                    if(cropImage==0)
                        XStart1= [subHereRel(1:3)];
                        XEnd1  = [subHereRel(4:6)];
                    end
                    
                    usedSub=sub;
            end   
        end
    end
    if(isempty(XStart1))
       display('This position has no overlap with the other position! '); 
    end

end


function [newStack]=rescaleStack(stack,voxSize,targetVoxSize)
   
 rescale_factor=voxSize./targetVoxSize;
  
 if(any(abs((rescale_factor-[1 1 1]))>0.05))
    display('Now resampling image');
    display(['Rescale factor:'  num2str(rescale_factor)])
    display(['Size before:'  num2str(size(stack))])
    tform=affine3d(diag([rescale_factor 1]));
    newStack=imwarp(stack,tform,'nearest','SmoothEdges',true);
    display(['Size after:'  num2str(size(newStack))]);
    newStack=newStack(2:end-1,2:end-1,2:end-1);
 else
    newStack=stack; 
 end

end





