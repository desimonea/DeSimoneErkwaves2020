function [myScale] = visualizeERKKTR(myScale,paths,opts)
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
      saveOptionsColor=[];
      saveOptionsColor.compress='lzw';
      saveOptionsColor.color=true;
      saveOptionsColor.overwrite=true;
      
      saveOptions=[];
      saveOptions.compress='lzw';
      saveOptions.overwrite=true;

      stackFolder=paths.stackFolder;
      basename=myScale.basename;
      suffix=paths.suffix;
        
      infoName=opts.infoName;
  
      centers=myScale.(infoName).centers;
      value= myScale.(infoName).(opts.valueToPlot); 
      
      idxs=~isnan(value);
      
      %setup colors- logarithmic
      colorLevels=100;
      colorRescale = 1.2;
      colors=parula(2*colorLevels+1);
      
      colors = colors.*repmat((((1:size(colors,1))./numel(colors).*(colorRescale-1)+1).^2)',1,3);
      %colors(:) = colors(:)+0.2;%colorRescale;
      colors(colors(:)>1)=1;
      
      
      if(isfield(opts,'limValue'))
          limLogValue=log(opts.limValue);
      else
          limLogValue=[min(log(value))  max(log(value))];
      end
      
      colorIdx=round(2*colorLevels*(log(value)-limLogValue(1))/(limLogValue(2)-limLogValue(1))+1);
      
      % deal with saturation - logarithmic
      colorIdx(colorIdx>2*colorLevels+1)=2*colorLevels+1;
      colorIdx(colorIdx<1)=1;
      
      %colors=parula(colorLevels);
      %limAct=max(erkAct);
      %colorIdx=floor(erkAct*colorLevels./limAct);  
      
      imgbox=opts.imgbox;
      x0img=imgbox(2)/2; y0img=imgbox(1)/2;
      
      if(strcmp(opts.visualType,'centers'))
          
        % read stack
        spath=[stackFolder basename suffix '.tif'];
        s=loadtiff(spath); 
        s=max(s,[],3);  
                   
        figure;
        x=1:size(s,2); 
        %plot(x,(x-f_ell.X0_in)*tan(angle)+f_ell.Y0_in,'r','LineWidth',3);
        imshow(s,opts.adj); hold on;
        scatter(centers(:,1),centers(:,2),100,colors(colorIdx',:),'.');
      
      elseif(strcmp(opts.visualType,'nuclei'))
          
        % read stack
        spath=[stackFolder basename suffix '.tif'];
        s=loadtiff(spath);   
        s=permute(s,[2 1 3]);
        
        % create color stacks
        scolor_g=zeros(size(s));
        scolor_r=zeros(size(s));
        scolor_b=zeros(size(s));
        obj=opts.obj;
        svList=opts.svList;
        
        % position center and angle
        centers=myScale.TGMM_hypo_eq_ch2.centers;
        x0=round(mean(centers(:,1))); y0=round(mean(centers(:,2)));
   
        if(isfield(opts,'thetaLandmark'))
            theta=opts.thetaLandmark;
        else
            xd=myScale.(opts.dirField).landmark(1)-x0; yd=myScale.(opts.dirField).landmark(2)-y0;  
            theta=atan2(yd,xd);
        end

        for c=1:numel(obj)    
       
            svIdx=obj(c).svIdx;
            cell_indices=[];
            
            for sv=1:numel(svIdx)
                cell_indices=vertcat(cell_indices,svList{svIdx(sv)+1});
            end
                
            colorIdx(c);
            
            if(~isnan(colorIdx(c)))
                scolor_r(cell_indices)=colors(colorIdx(c),1);
                scolor_g(cell_indices)=colors(colorIdx(c),2);
                scolor_b(cell_indices)=colors(colorIdx(c),3);
            else
                scolor_r(cell_indices)=1;
                scolor_g(cell_indices)=1;
                scolor_b(cell_indices)=1;
            end
            
        end
        
        scolor_r_proj=max(scolor_r,[],3);
        scolor_g_proj=max(scolor_g,[],3);
        scolor_b_proj=max(scolor_b,[],3);
        scolor_rgb=cat(3,scolor_r_proj,scolor_g_proj,scolor_b_proj);
        scolor_rgb=permute(scolor_rgb,[2 1 3]);

        sproj=imadjust(max(s,[],3),opts.adj'./255,[0;1]);
        sproj=repmat(sproj,[1 1 3]);
        sproj=permute(sproj,[2 1 3]);      
        
        % create img projection
        sprojBig      =zeros(imgbox(1),imgbox(2),3);
        scolor_rgbBig =zeros(imgbox(1),imgbox(2),3);
        
        
        Rimg = imref2d(size(sproj));
        Rout= imref2d(imgbox);
        
        transl = affine2d([1 0 0; 0 1 0; -x0+x0img -y0+y0img 1]);
          
        % translation x'=x-x0+x0img
%         xStart=max(1-x0+x0img,1);
%         yStart=max(1-y0+y0img,1);
%         
%         xEnd=max(size(sproj,2)-x0+x0img,imgbox(2));
%         yEnd=max(size(sproj,1)-y0+y0img,imgbox(1));
%         
%         sprojBig(yStart:yEnd,xStart:xEnd,1:3)=sproj;
%         scolor_rgbBig(yStart:yEnd,xStart:xEnd,1:3)=scolor_rgb;
        
        sprojBig = imwarp(sproj,Rimg,transl,'OutputView',Rout);
        scolor_rgbBig = imwarp(scolor_rgb,Rimg,transl,'OutputView',Rout);

        % rotation
        %sprojBig=     imrotate(sprojBig,theta*180/pi,'crop');
        %scolor_rgbBig=imrotate(scolor_rgbBig,theta*180/pi,'crop');
        
        sprojBig=uint8(sprojBig); 
        sprojBig=max(sprojBig,[],3);
        sprojBig = adapthisteq(sprojBig,'ClipLimit',0.004);
     
        scolor_rgbBig=uint8(scolor_rgbBig.*255);
               
        sshow=cat(2,sproj,scolor_rgb.*255);
           
        centers=myScale.(opts.infoName).centers;
        x0=mean(centers(:,1)); y0=mean(centers(:,2));
        
        if(isempty(opts.chooseSource) | ~opts.chooseSource)
                myScale.(opts.dirField).source=[-1 0];
                myScale.(opts.dirField).angleSource=0; 
        else(opts.chooseSource)
                imshow(scolor_rgb);
                [xd,yd,~]=impixel;
                myScale.(opts.dirField).source=[xd(1) yd(1)];
                myScale.(opts.dirField).angleSource=atan2(yd(1)-y0,xd(1)-x0);
        end
        
        figure;
        imshow(sprojBig,'Border','tight'); 
        if(opts.print)
            saveastiff(sprojBig,[paths.outFolder myScale.basename '_hypoxquant.tif'],saveOptions);
        end
        
        figure;
        imshow(scolor_rgbBig,'Border','tight');
        hold on;
        plot(x0img+[0 1000.*cos(myScale.(opts.dirField).angleSource)],y0img+[0 1000.*sin(myScale.(opts.dirField).angleSource)] ,'w--','LineWidth',10);
        if(opts.print)
            saveastiff(scolor_rgbBig,[paths.outFolder myScale.basename '_quant.tif'],saveOptionsColor);
        end

      else
         error('I do not know this opts.visualType - centers or nuclei'); 
      end
   
%          figure;
%          c=colorbar('Ticks',[0 0.5 1],'TickLabels',{-1,0,1},'FontSize',24,'Color','white');
%          c.Label.String='Renormalized log cytoplasmic/nuclear signal';
%          c.Label.Color='white';
%          set(gcf,'color','black','InvertHardCopy','off')
   
         myScale.(opts.dirField).thetaLandmark=theta;

end
