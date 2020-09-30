function [divMsmooth,coord,vvLast] = flowsScale(paths,myScaleSeries,track,opts)
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
%   plot growth maps

    nameSeries = paths.nameSeries;

    hppTrues = myScaleSeries.hppTrues;

    frames = track.frame+1;
    frameList = unique(frames);
    xc = track.Center_of_the_object_0;
    yc = track.Center_of_the_object_1;
    zc = track.Center_of_the_object_2;
    trackId = track.trackId;

    printHD = opts.printHD;

    frameMin =opts.frameMin;
    frameMax = opts.frameMax;
    dFrame = opts.dFrame;

    mkdir(paths.outFolder);
    pxsize= opts.pxsize;

    xStart = opts.xStart; yStart = opts.yStart; 
    xEnd = opts.xEnd; yEnd = opts.yEnd; 
    dx = opts.dx;
    dy = opts.dy;
    plotRange = opts.plotRange;
    
    quiverStep =2;

    saveOptions= [];
    saveOptions.overwrite = 'true';
    saveOptions.compress = 'lzw';

    vx_2all = [];
    vy_2all = [];

    for frameStart = frameMin:frameMax-dFrame

        disp(['Processing frame: ' num2str(frameStart)]);
        idx0 = find(frames== frameStart & trackId >1);

        v= [];

        dt=hppTrues(frameStart+dFrame)-hppTrues(frameStart);

        xcenter = mean(xc(idx0)); ycenter = mean(yc(idx0));

        for i=1:numel(idx0)
            frameAll = frames(find(trackId == trackId(idx0(i))));
            idx1 = find(trackId == trackId(idx0(i)) & frames== frameStart+dFrame);

            if(isempty(idx1) | ~setdiff(frameList(frameMin:frameMax),frameAll))
                continue
            end
            v = vertcat(v, [xc(idx0(i))-xcenter yc(idx0(i))-ycenter (xc(idx1)- xc(idx0(i))).*pxsize./dt (yc(idx1)-yc(idx0(i))).*pxsize./dt]);
        end  

        vNorm = ((sqrt(v(:,3).^2+v(:,4).^2)));
        idxGood = vNorm<2.4; %4 without normalization
        v = v(idxGood,:);

        xMin=xStart-dx/2; xMax = yEnd+dx/2;
        yMin=xStart-dx/2; yMax=  yEnd+dx/2;
        Xedges = xMin:dx:xMax;
        Yedges = yMin:dy:yMax;

        [N,Xedges,Yedges,binX,binY] = histcounts2(v(:,1),v(:,2),Xedges,Yedges);
        [Xi,Yi]= meshgrid(1:numel(Xedges)-1,1:numel(Yedges)-1);
        [Xireal,Yireal]= meshgrid(Xedges(1:end-1)+dx/2,Yedges(1:end-1)+dy/2);
        [XiHD,YiHD]= meshgrid(xMin:xMax,yMin:yMax);
        vCoarse = nan(numel(Xi),4);

        for k=1:numel(Xi)  
           idxHere = binX==Xi(k) &  binY==Yi(k);
           vCoarse(k,1:4) = [Xireal(k) Yireal(k) nanmean(v(idxHere,3)) nanmean(v(idxHere,4))];
        end

        idxConvHull = convhull(v(:,1),v(:,2));
        mask=inpolygon(Xireal,Yireal,v(idxConvHull,1),v(idxConvHull,2));

        vxM = reshape(vCoarse(:,3),[numel(Xedges)-1 numel(Yedges)-1]);
        vyM = reshape(vCoarse(:,4),[numel(Xedges)-1 numel(Yedges)-1]);

        lmean1=2;
        K1 = ones(lmean1)./lmean1^2;
        vxM = nanconv(vxM,K1,'same','noedge');
        vyM = nanconv(vyM,K1,'same','noedge');

        divM = divergence(Xi,Yi,vxM,vyM)./(pxsize.*dx);
        
        lmean2 = 5;
        K2 = fspecial('gaussian',5,3);
        divMsmooth = nanconv(divM,K2,'same','noedge');
        divMsmooth(~mask)=NaN;

        if(printHD)
            figure;
            idxNotNan = ~isnan(divMsmooth);
            divMHD = griddata(Xireal(idxNotNan),Yireal(idxNotNan),divMsmooth(idxNotNan),XiHD,YiHD,'cubic');
            maskHD=inpolygon(XiHD,YiHD,v(idxConvHull,1),v(idxConvHull,2));
            divMHD(~maskHD) = NaN;
            colormap([0 0 0; parula(100)]);    
            imagesc(divMHD,plotRange); hold on;
            plot(v(idxConvHull,1),v(idxConvHull,2),'w','LineWidth',2);
            axis off
            hold off;
            xlim = [XStart xEnd]/dx;
            ylim = [yStart yEnd]/dy;
            print([paths.outFolder nameSeries '_div_t' num2str(frameStart) ''],'-r100','-dtiff');
        else
            figure;
            colormap([0 0 0; parula(100)]);    
            colormap([0 0 0; brewermap(1000,'RdYlGn')]);
            imagesc(divMsmooth,plotRange); hold on; %
            axis off
            xlim = [xStart xEnd]/dx;
            ylim = [yStart yEnd]/dy;
           
            axis equal;
            axis tight;
            ax=gca;
            ax.Position=[0 0 1 1];
            set(gcf, 'InvertHardCopy', 'off');
            %print([paths.outFolder nameSeries '_div_t' num2str(frameStart) ''],'-r100','-dpng');   
        end

           vxMPlot = vxM; vyMPlot = vyM;
           vxMPlot(~mask)=NaN; vyMPlot(~mask)=NaN;
           scaleFactor = 60;
           hold on;
 
           myquiver(Xi,Yi,vxMPlot.*scaleFactor./dx,vyMPlot.*scaleFactor./dy,quiverStep);

           axis off
           xlim = [xStart./dx xEnd./dx];
           ylim = [yStart./dy yEnd./dy];
           gcfHere = gcf;
           set(gcf,'units','centimeters','Position',[gcfHere.Position(1) gcfHere.Position(2) 44 44]);
           ax=gca;
           axis equal;
           axis tight;
           set(gcf, 'InvertHardCopy', 'off');
           ax.Position=[0 0 1 1];
           print([paths.outFolder nameSeries '_vel_t' num2str(frameStart) ''],'-r100','-dpng');
           vx_2all = cat(3,vx_2all , vxM);
           vy_2all = cat(3,vy_2all , vyM);
           
           coord = cat(3, Xireal, Yireal);
           vvLast = cat(3, vxM, vyM);
           
    end

end

function myquiver(Xi,Yi,vx,vy,quiverStep)

   lmean1=3;
   K1 = ones(lmean1)./lmean1^2;
   myvx =  nanconv(vx,K1,'same','nanout');
   myvy =  nanconv(vy,K1,'same','nanout');

   lmean1=3;
   K1 = ones(lmean1)./lmean1^2;
   myvx =  nanconv(vx,K1,'same','nanout');
   myvy =  nanconv(vy,K1,'same','nanout');

   quiver(Xi(1:quiverStep:end),Yi(1:quiverStep:end),myvx(1:quiverStep:end).*quiverStep,myvy(1:quiverStep:end).*quiverStep,0,'color','b','LineWidth',2);

end


