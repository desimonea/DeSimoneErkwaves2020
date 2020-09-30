function [myScale,obj] = assignKTR(myScale,obj,svList,paths,opts)
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
   
   inFolder=paths.inFolder;
   basename=myScale.basename;
   
   infoName=opts.infoName;
   
   KTRsuffix=  paths.KTRsuffix;
   H2Asuffix=  paths.H2Asuffix;
   
   sktr= loadtiff([inFolder basename KTRsuffix '.tif']);
   sh2a= loadtiff([inFolder basename H2Asuffix '.tif']); 
   sktr=permute(sktr,[2 1 3]); sh2a=permute(sh2a,[2 1 3]);
   
    sizes=size(sktr);
   
  % run through objects and assigns colors
  verbose = opts.verbose;
  
  colorLevels=100;
  
  if(verbose)
      scolor_g=zeros(size(sktr));
      scolor_r=zeros(size(sktr));
      scolor_b=zeros(size(sktr));
      colors=parula(2*colorLevels);
      limLogAct=2;
  end

  cent=round(vertcat(obj.m));
  
  
  if(isfield(obj,'area'))
    obj =rmfield(obj,'area');
  end
  
  if(isfield(obj,'volume'))
    obj =rmfield(obj,'volume');
  end
  
  if(isfield(obj,'volume'))
    obj =rmfield(obj,'ktr');
  end

  
  
  for c=1:numel(obj)    
     
     svIdx=obj(c).svIdx;
     cell_indices=[];
            
     for sv=1:numel(svIdx)
        cell_indices=vertcat(cell_indices,svList{svIdx(sv)+1});
     end
 
     %% calculate cytoplasmic and nuclear signal ERK KTR

     % cut sub-windows
     smask=zeros(size(sh2a),'logical');
     smask(cell_indices)=1;
     dx = round(10); dy = round(10); dz = round(15);
     cent=round(obj(c).m);     

     ymin=max(cent(1)-dy,1);
     ymax=min(cent(1)+dy,size(smask,1));
     xmin=max(cent(2)-dx,1);
     xmax=min(cent(2)+dx,size(smask,2));
     zmin=max(cent(3)-dz,1);
     zmax=max(min(cent(3)+dz,size(smask,3)),1);

     minimask=  smask(ymin:ymax,xmin:xmax,zmin:zmax);
     mini_sh2a= double(sh2a(ymin:ymax,xmin:xmax,zmin:zmax));
     mini_sktr= double(sktr(ymin:ymax,xmin:xmax,zmin:zmax));

     area = squeeze(sum(sum(max(minimask,[],3),1),2));
     if(isempty(area))
         area=0;
     end

     [~,zbest]=max(squeeze(mean(mean(double(minimask),1),2)));
     minimask =  minimask(:,:,zbest);
     mini_sh2a=  mini_sh2a(:,:,zbest);
     mini_sktr=  mini_sktr(:,:,zbest);

     % calculate cyt and nuc masks
     mm_nuc=imclose(minimask,strel('disk',3));
     mm_cyt=logical(imdilate(minimask,strel('disk',3))-mm_nuc);

     % calculate signals
     volNuc=sum(double(mm_nuc(:)));
     volCyt=sum(double(mm_cyt(:)));
     signalCyt=sum(mini_sktr(mm_cyt(:)))./volCyt;
     signalNuc=sum(mini_sktr(mm_nuc(:)))./volNuc;
     averageH2ANuc=sum(mini_sh2a(mm_nuc(:)))./volNuc;
     averageH2ACyt=sum(mini_sh2a(mm_cyt(:)))./volCyt;
     ratioH2ANucCyt=averageH2ANuc./averageH2ACyt;
     ktr=signalCyt./signalNuc;

    % visualize masks
    if(verbose==2)   

        sizexy=[size(minimask,1) size(minimask,2)];
        red =   cat(3, ones(sizexy),zeros(sizexy),zeros(sizexy));
        green = cat(3, zeros(sizexy),ones(sizexy),zeros(sizexy));     

        figure;
        % cyto mask on nuclei
        imshow(max(mini_sh2a,[],3),[],'InitialMagnification',1000);
        hold on 
        %h = imshow(green); 
        %set(h,'AlphaData',double(max(mm_nuc,[],3))*0.1)
        h2 = imshow(red); 
        set(h2,'AlphaData',double(max(mm_cyt,[],3))*0.5)
        hold off
        waitforbuttonpress

        figure;
        %  mask on nuclei 
        imshow(max(mini_sh2a,[],3),[],'InitialMagnification',1000);
        hold on 
        h = imshow(red); 
        set(h,'AlphaData',bwperim(imdilate(double(max(mm_nuc,[],3)),strel('disk',1)))*0.8);
        %h2 = imshow(red); 
        %set(h2,'AlphaData',double(max(mm_cyt,[],3))*0.1)
        hold off
        waitforbuttonpress

        figure;
        %  nuclear mask on sensor 
        imshow(max(mini_sktr,[],3),[],'InitialMagnification',1000);
        hold on 
        h = imshow(red); 
        set(h,'AlphaData',double(max(mm_nuc,[],3))*0.2)
        %h2 = imshow(red); 
        %set(h2,'AlphaData',double(max(mm_cyt,[],3))*0.1)
        hold off
        waitforbuttonpress

        figure;
        % cyto mask on sensor
        imshow(max(mini_sktr,[],3),[],'InitialMagnification',1000);
        hold on 
        %h = imshow(green); 
        %set(h,'AlphaData',double(max(mm_nuc,[],3))*0.1)
        h2 = imshow(green); 
        set(h2,'AlphaData',double(max(mm_cyt,[],3))*0.2)
        hold off
        waitforbuttonpress

        ktr
    end
     
%%     SAVE STUFF INTO OBJ
     
       obj(c).volume=numel(cell_indices);
       obj(c).ktr=ktr;
       obj(c).signalCyt=signalCyt;
       obj(c).signalNuc=signalNuc;
       obj(c).area=area;
       obj(c).volNuc= volNuc;
       obj(c).volCyt= volCyt;
       obj(c).averageH2ANuc=averageH2ANuc;
       obj(c).averageH2ACyt=averageH2ACyt;
       obj(c).ratioH2ANucCyt=ratioH2ANucCyt;

  end   
   
  if(verbose)
    scolor_r_proj=max(scolor_r,[],3);
    scolor_g_proj=max(scolor_g,[],3);
    scolor_b_proj=max(scolor_b,[],3);
    imshow(cat(3,scolor_r_proj,scolor_g_proj,scolor_b_proj));
  end

    opts.saveNameCyt = [opts.saveName 'Cyt'];
    opts.saveNameNuc = [opts.saveName 'Nuc'];

    myScale.(infoName).centers = vertcat(obj.m);
    myScale.(infoName).volume           =    vertcat(obj.volume);
    myScale.(infoName).areaNucleus      =    vertcat(obj.volNuc);
    myScale.(infoName).areaCyt          =    vertcat(obj.volCyt);

    myScale.(infoName).(opts.saveName)  =    vertcat(obj.ktr);
    myScale.(infoName).(opts.saveNameCyt)  = vertcat(obj.signalCyt);
    myScale.(infoName).(opts.saveNameNuc)  = vertcat(obj.signalNuc);

    myScale.(infoName).area             =    vertcat(obj.area);

    myScale.(infoName).averageH2ANuc           =     vertcat(obj.averageH2ANuc);
    myScale.(infoName).averageH2ACyt           =     vertcat(obj.averageH2ACyt);
    myScale.(infoName).ratioH2ANucCyt           =    vertcat(obj.ratioH2ANucCyt);

   
end
   


