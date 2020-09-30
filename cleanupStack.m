function [stackOut,roi] = cleanupStack(stackIn,roi,opts)
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

stay=1;
if(isfield(opts,'cleanupPlanes'))
    if(opts.cleanupPlanes==0)
        stay=0;
    end
end

permutation = [1 2 3];

stack=stackIn;
sizeIn=size(stackIn);

largestep=10;
    
n=1;  
 
if(isempty(roi))
   roi=ones(size(stack),'logical');
end

[stack,roi] = selectROIstack(stack,roi,opts);

adj=[0 0.1];
size_old=size(stack);
size_new=[600 size_old(2)*600/size_old(1)];
h = figure('Position', [0 0 size_new(2) size_new(1)]);
fpos = get(h,'Position');
axOffset = (fpos(3:4)-[size_new(2) size_new(1)])/2;
ha = axes('Parent',h,'Units','pixels','Position',[axOffset size_new(2) size_new(1)]);
imshow(imadjust(stack(:,:,n),adj,[0 1]),'Parent',ha);

forbidUndo=0;
while(stay==1)
  
        display('Please select regions to clean (h for help)');
        imshow(imadjust(stack(:,:,n),adj,[0 1]),'Parent',ha);
       
        display(n);
        waitforbuttonpress;
        button = get(h,'CurrentCharacter');
        display(button);
        
        if (button == 'h')
           
            disp(' a -> clean projections ')
            disp(' c -> clean region selected')
            disp(' m -> clean same region as previous')
            disp(' u -> undo')
            disp(' s -> move to z = start')
            disp(' e -> move to z = end')
            disp(' ,/< -> step backward')
            disp(' ,/< -> step forward')
            disp(' -/_ -> fast backward')
            disp(' -/+ -> fast forward')
            disp(' r -> resize ');
            disp(' o -> return to original size ');
            disp(' p -> permute');
            disp(' q -> quit when finished')
              
        elseif (button == 'c') % clean
        
            n_old=n;
            img_old=stack(:,:,n);
            img_roi_old=roi(:,:,n);
            
            mask=roipoly;
        
            img_roi=roi(:,:,n);
            img_roi(mask)=0;
            roi(:,:,n)=img_roi;
        
            img=stack(:,:,n);
            img(mask)=0;
            stack(:,:,n)=img;  
            
            %h=figure;
            forbidUndo=0;
            
        elseif (button == 'm') % use the same mask used before 
            
            n_old=n;

            img_old=stack(:,:,n);
            img_roi_old=roi(:,:,n);
            
            img_roi=roi(:,:,n);
            img_roi(mask)=0;
            roi(:,:,n)=img_roi;
        
            img=stack(:,:,n);
            img(mask)=0;
            stack(:,:,n)=img;    
            
            forbidUndo=0;
             
        elseif (button == 'u') % undo
            if(forbidUndo)
                disp('Undo is forbidden now');
            else  
                roi(:,:,n_old)=img_roi_old;
                stack(:,:,n_old)=img_old;
            end
        elseif (button == ',' | button == '<') % <
             n=max(n-1,1);
             forbidUndo=0;
        elseif (button == '.' | button == '>') % >
             n=min(n+1,size(stack,3));
             forbidUndo=0;
        elseif button == '-' % -
             n=max(n-largestep,1);
             forbidUndo=0;
        elseif (button == '+' | button == '=') % +
             n=min(n+largestep,size(stack,3));
             forbidUndo=0;
        elseif button == 's' % +
             n=1;
             forbidUndo=0;
        elseif button == 'e' % +
             n=size(stack,3);
             forbidUndo=0;
        elseif button == 'p'     
             vie=input('Which view would you like? (1: xy - 2: yz - 3: xz) ','s');
             forbidUndo=1;
             
              %perform counterpermutation
             
             if(vie == '1') 
                 roi=permute(roi,permutation);
                 stack=permute(stack,permutation);
                 permutation=[1 2 3]; 
             elseif(vie == '2')
                 roi=permute(roi,permutation);
                 stack=permute(stack,permutation);
                 permutation=[1 3 2];
                 stack=permute(stack,permutation);
                 roi=permute(roi,permutation);
             elseif(vie =='3')
                 roi=permute(roi,permutation);
                 stack=permute(stack,permutation);
                 permutation=[3 2 1];
                 stack=permute(stack,permutation);
                 roi=permute(roi,permutation);
             else
                disp('I do not know this permutation :(');
             end
              
             n=1;
             mask=[];
  
        elseif button == 'r'
             
            inp=input('Please indicate the three rescale factors in [i j k] : ');

            if(numel(inp)<3)
                disp('Wrong number of inputs (they must be three)');
            else

                sizeNew=round([inp(1)*size(stack,1),inp(2)*size(stack,2),inp(3)*size(stack,3)]);
                
                stack=imresize(stack,[sizeNew(permutation(1)) sizeNew(permutation(2))]);
                roi=  imresize(roi,  [sizeNew(permutation(1)) sizeNew(permutation(2))]);
                  
                roi=permute(roi,[1 3 2]);
                stack=permute(stack,[1 3 2]);
                
                stack=imresize(stack,[sizeNew(permutation(1)) sizeNew(permutation(3))]);
                roi=  imresize(roi,  [sizeNew(permutation(1)) sizeNew(permutation(3))]);
                
                roi=permute(roi,[1 3 2]);
                stack=permute(stack,[1 3 2]);
                  
            end
            
           
             n=1;
             mask=[];
             forbidUndo=1;
            
        elseif button == 'o'
             
            roi=permute(roi,permutation);

            roi=double(roi);
            roi=   imresize(roi,sizeIn(1:2));

            roi   = permute(roi  ,[1 3 2]);

            roi=   imresize(roi,[sizeIn(1) sizeIn(3)]);

            roi   = permute(roi  ,[1 3 2]);
            
            roi=round(roi);
            roi=logical(roi);

            stack=stackIn;
            stack(~roi)=0;
            
            
        elseif button== 'a'    
            
            [stack,roi] = selectROIstack(stack,roi);
            mask=[];
            forbidUndo=1;
             
        elseif button == 'q'
          
            roi=permute(roi,permutation);

            roi=double(roi);
            roi=   imresize(roi,sizeIn(1:2));

            roi   = permute(roi  ,[1 3 2]);

            roi=   imresize(roi,[sizeIn(1) sizeIn(3)]);

            roi   = permute(roi  ,[1 3 2]);
            
            roi=round(roi);
            roi=logical(roi);

            stack=stackIn;
            stack(~roi)=0;
            
            imshow(imadjust(max(stack,[],3),adj,[0 1]),'Parent',ha);
           
            answ=input('Are you satisfied?','s');
            
            if(strcmp(answ,'y'))
                stay=0;                
            end
        end
end

stackOut=stack;

end