function [regInfo]=imregisterStackSeriesXY3D(myScaleSeries,paths,opts,regInfo)
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
    if(nargin<4)
       txyAbs = [];
       imSizeNew = [];
    else
       txyAbs = regInfo.txyAbs;
       imSizeNew =  regInfo.imSizeNew;
       RNew =  regInfo.RNew;
    end

    refCh=opts.refCh;
    chToRegister=opts.chToRegister;

    nFrames=numel(myScaleSeries.objNames);

    inSuffix=paths.inSuffix;
    outSuffix=paths.outSuffix;
    toRegSuffix=paths.toRegSuffix;

    if(isempty(txyAbs))

            WnewAll=nan(4,2,nFrames);
            txyAbs=nan(3,3,nFrames);
            txyAbs(:,:,end)=diag([1 1 1]);

            load([paths.objFolder myScaleSeries.objNames{end}]);
            s2 = loadtiff([paths.inFolder myScale.basename '_ch' num2str(refCh) inSuffix]);

            s2proj= adapthisteq(max(s2,[],3),'NumTiles',[16 16],'ClipLimit',0.01);
            R1 = imref2d([size(s2,1) size(s2,1)],1,1);
            W1 = [R1.YWorldLimits(1) R1.XWorldLimits(1)];   
            W2 = [R1.YWorldLimits(1) R1.XWorldLimits(2)]; 
            W3 = [R1.YWorldLimits(2) R1.XWorldLimits(1)];  
            W4 = [R1.YWorldLimits(2) R1.XWorldLimits(2)];  
            WnewAll(:,:,end)=[W1; W2; W3; W4];

            verbose=0;

            for t=nFrames-1:-1:1

                disp(['Calculating registration for t = ' num2str(t)]);

                load([paths.objFolder myScaleSeries.objNames{t}]);     myScale1=myScale;

                s1=loadtiff([paths.inFolder myScale1.basename '_ch' num2str(refCh) inSuffix]);

                s1proj= adapthisteq(max(s1,[],3),'NumTiles',[16 16],'ClipLimit',0.01);

                s1proj=double(imadjust(s1proj,[]));
                s2proj=double(imadjust(s2proj,[]));

                [~,txy,Wnew,R1,R2]=registerImg(s1proj,s2proj,opts);
                txyAbs(:,:,t)=txy;
                WnewAll(:,:,t)=Wnew;

                % apply complete transformation to original stack 
                %s1New=zeros([size(s2,1) size(s2,2) size(s1,3)],'uint8');
                s1New = imwarp(s1,R1,affine2d(txyAbs(:,:,t)),'OutputView',R2);
                s1projNew=max(s1New,[],3);

                s2= s1New;
                s2proj= s1projNew;

            end

            Y1WorldNew(1,1)= round(min(min(WnewAll(1:4,1,1:nFrames),[],1),[],3))+0.5;
            Y1WorldNew(1,2)= round(max(max(WnewAll(1:4,1,1:nFrames),[],1),[],3))+0.5;

            X1WorldNew(1,1)= round(min(min(WnewAll(1:4,2,1:nFrames),[],1),[],3))+0.5;
            X1WorldNew(1,2)= round(max(max(WnewAll(1:4,2,1:nFrames),[],1),[],3))+0.5;

            imSizeNew=([Y1WorldNew(2)-Y1WorldNew(1) X1WorldNew(2)-X1WorldNew(1)]);
            Rnew=imref2d(imSizeNew,X1WorldNew,Y1WorldNew);
    end

    options=[];
    options.overwrite=true;
    options.compress='lzw';

    for t=1:nFrames

        disp(['Registering t = ' num2str(t)]);

        load([paths.objFolder myScaleSeries.objNames{t  }]); 

        for ch = chToRegister 

            s1=loadtiff([paths.toRegFolder myScale.basename '_ch' num2str(ch) toRegSuffix]);

            s1New=zeros([imSizeNew size(s1,3)],'uint8');

            R1   = imref2d([size(s1,1) size(s1,2)],1,1);
            s1New = imwarp(s1,R1,affine2d(txyAbs(:,:,t)),'OutputView',Rnew);

            saveastiff(s1New,          [paths.outFolder myScale.basename '_ch' num2str(ch) outSuffix],options);
            saveastiff(max(s1New,[],3),[paths.outFolder myScale.basename  '_ch' num2str(ch) 'maxproj' outSuffix],options);
        end
    end

    regInfo=[];
    regInfo.txyAbs=txyAbs;
    regInfo.Rnew=Rnew;
    regInfo.imSizeNew=imSizeNew;


end

%%
function [s1projNew,txy,Wnew,R1,R2]=registerImg(s1proj,s2proj,opts)
    
    verbose = opts.verbose;
    useCrossCorr = opts.useCrossCorr;
    optimizer = opts.optimizer;
    metric = opts.metric;

    s1projlarge=zeros(max(size(s1proj,1),size(s2proj,1)),max(size(s1proj,2),size(s2proj,2)));

    x1clarge=floor(size(s1projlarge,2)./2); 
    y1clarge=floor(size(s1projlarge,1)./2);

    xStart=x1clarge-floor(size(s1proj,2)./2);
    yStart=y1clarge-floor(size(s1proj,1)./2);

    s1projlarge(yStart+1:yStart+size(s1proj,1),xStart+1:xStart+size(s1proj,2))=s1proj;  %edit 05jan18
    t_large=[1 0 0; 0 1 0; xStart-1 yStart-1 1];

    idxNull1=s1projlarge(:)<4;
    idxNull2=s2proj(:)<4;

    s1projlarge(idxNull1)=floor(rand(sum(idxNull1),1).*0);
    s2proj(idxNull2)=     floor(rand(sum(idxNull2),1).*0);

    if(useCrossCorr)
        % translation with correlation
        ymin1=1; 
        ymax1=min(size(s1projlarge,1),size(s2proj,1)-1);
        xmin1=1; 
        xmax1=min(size(s1projlarge,2),size(s2proj,2)-1);

        s1crop=s1projlarge(ymin1:ymax1,xmin1:xmax1);

        c = normxcorr2(s1crop,s2proj);
        [ypeak, xpeak] = find(c==max(c(:)));
        yoffSet = ypeak-size(s1crop,1)-ymin1;
        xoffSet = xpeak-size(s1crop,2)-xmin1;
    else
        xoffSet=0; yoffSet=0;
    end
    
    R1large   =  imref2d(size(s1projlarge),1,1);
    R2   =  imref2d(size(s2proj),1,1);

    t1_2a=[1 0 0 ;0 1 0; xoffSet yoffSet 1];
    s1projNewa = imwarp(s1projlarge,R1large,affine2d(t1_2a),'OutputView',R1large);
    R1a   = imref2d(size(s1projNewa),1,1);
   
    %[optimizer, metric] = imregconfig('monomodal');
    t1_2=imregtform(s1projNewa,R1a,s2proj,R2,'rigid',optimizer, metric);
    txy=t1_2.T*t1_2a*t_large; % complete transformation

    % calculate new limits
    R1   = imref2d([size(s1proj,1) size(s1proj,2)],1,1);
    W1 = [R1.YWorldLimits(1) R1.XWorldLimits(1)];   
    W2 = [R1.YWorldLimits(1) R1.XWorldLimits(2)]; 
    W3 = [R1.YWorldLimits(2) R1.XWorldLimits(1)];  
    W4 = [R1.YWorldLimits(2) R1.XWorldLimits(2)];
    Wnew = transformPointsForward(affine2d(txy),[W1; W2; W3; W4]);   

    s1projNew = imwarp(s1proj,R1,affine2d(txy),'OutputView',R2);
    
    if(verbose)
        % visualize process
        figure;
        imshowpair(s2proj,s1proj);
        figure;
        imshowpair(s2proj,s1projlarge)
        figure;
        imshowpair(s2proj,s1projNewa);
        figure;
        imshowpair(s2proj,s1projNew);
        waitforbuttonpress;
    end
end


