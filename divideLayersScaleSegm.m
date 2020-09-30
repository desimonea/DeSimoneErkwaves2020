function [myScale,zprofile] = flattenLayersScaleSegm(myScale,refCh,todivideCh,paths,divOpts)
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
   verbose=divOpts.verbose;

   basename=myScale.basename;
   
   inFolder = paths.inFolder;
   outFolder= paths.outFolder;
   objFolder= paths.objFolder;
   refFolder= paths.refFolder;
   
   depth_fit=15; %in planes for now
   zshiftHypo=0;
   
   options=[];
   options.overwrite='true';
   options.compress='lzw';
   
   spathRef=[refFolder basename '_ch' num2str(refCh) '.tif'];   
   s=loadtiff(spathRef);
   
   zprofiles=[]; counter=1;
   
   if(strcmp(divOpts.methodHypo,'manual'))
        display('Please choose the last plane of the hypo projection');
        [zHypo] = planeSelect(s);
   elseif(strcmp(divOpts.methodHypo,'peak'))
        zprofile=squeeze(sum(sum(s,1),2));
        
        if(verbose)
            figure;
            plot(zprofile);
        end
        [m,zHypo]=max(zprofile);
        zHypo=min(zHypo)+zshiftHypo;  
   elseif(strcmp(divOpts.methodHypo,'peak2'))
        
        zprofile=squeeze(sum(sum(s,1),2));
        
        [zpeak,zcenter]=max(zprofile);
        
        xx=max(1,zcenter-depth_fit):min(numel(zprofile)+1,zcenter+depth_fit);
        zp=zprofile(xx); 
        
        fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[max(zp)/10,max(zp)/10,5,5,0,0],...
               'Upper',[max(zp),max(zp),numel(zp),numel(zp),25,25],...
               'StartPoint',[max(zp) max(zp) 0 numel(zp) 13 13]);
        ft = fittype('a1*exp(-((x-m1)/s1)^2)+a2*exp(-((x-m2)/s2)^2)','options',fo);

        [curve2,gof2] = fit(xx',zp,ft);
       
        if(verbose)
            curve2
            figure;
            plot(curve2,xx,zp);
        end
        
        zHypo=round(max(min(curve2.m1,curve2.m2),0));
    
   elseif(isnum(divOpts.methodHypo))
        zHypo=method;
   else
       error('I do not know what method to use to choose the stack!');
   end 
        
    img_hypo=max(s(:,:,1:zHypo),[],3);
%          
    if(strcmp(divOpts.methodEpi,'manual'))
        display('Please choose the first plane of the epi projection');
        [zEpi] = planeSelect(s);   
    elseif(isnumeric(divOpts.methodEpi))     
         zEpi=zHypo+divOpts.methodEpi;
    elseif(strcmp(divOpts.methodEpi,'correlative'))
        for z=zHypo:numel(zprofile)
            img_shared=max(s(:,:,zHypo:z),[],3);
            img_epi=max(s(:,:,z:end),[],3);
            corrz1(z-zHypo+1)=corr2(img_shared,img_epi); % to be maximized
            corrz2(z-zHypo+1)=corr2(img_hypo,img_epi); % to be minimized
            corrtot(z-zHypo+1)=corrz2(z-zHypo+1)*3-corrz1(z-zHypo+1);
        end
        [~,zEpi]=min(corrtot);zEpi=zEpi+zHypo-1;
    end

    if(verbose==2)
        img_hypo=max(s(:,:,1:zHypo),[],3);
        img_epi=max(s(:,:,zEpi:end),[],3);
        img_shared=max(s(:,:,zHypo:zEpi),[],3);
        %img_tot=(img_hypo+img_epi)*2; img_tot(:,:,2)=img_shared*2; img_tot(:,:,3)=0;
        %imshow(img_tot);
        figure;imshow(img_hypo);
        figure;imshow(img_epi);
    end  
        
   clear('s');
   
   
   %%
   options=[];
   options.compress='lzw';
   options.overwrite='true';
     
   for ch=todivideCh
       
       spathToDiv=[inFolder basename '_ch' num2str(ch) '.tif'];

       sToDiv=loadtiff(spathToDiv);
       img_hypo= max(sToDiv(:,:,1:zHypo(1)),[],3);
       img_epi=max(sToDiv(:,:,zEpi:end),[],3);
       img_shared=max(sToDiv(:,:,zHypo:zHypo),[],3);
       
       % projections
       wpathHypo=[outFolder myScale.basename '_ch' num2str(ch) '_hypomaxproj.tif'];
       saveastiff(img_hypo,wpathHypo,options);
       wpathEpi=[outFolder myScale.basename '_ch' num2str(ch) '_epimaxproj.tif'];
       saveastiff(img_epi,wpathEpi,options);
       wpathShared=[outFolder myScale.basename '_ch' num2str(ch) '_sharedmaxproj.tif'];
       saveastiff(img_shared,wpathShared,options);

       % stacks
       wpathHypo=[outFolder myScale.basename '_ch' num2str(ch) '_hypo.tif'];
       saveastiff(sToDiv(:,:,1:zHypo(1)),wpathHypo,options);
       wpathEpi=[outFolder myScale.basename '_ch' num2str(ch) '_epi.tif'];
       saveastiff(sToDiv(:,:,zEpi:end),wpathEpi,options);
       wpathShared=[outFolder myScale.basename '_ch' num2str(ch) '_shared.tif'];
       saveastiff(sToDiv(:,:,zHypo:zHypo),wpathShared,options);
       
   end
   
end
   
   
   
   