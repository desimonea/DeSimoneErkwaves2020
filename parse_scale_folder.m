function [mystack] = parse_scale_folder(folder,basename,opts)
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
stack_dir=dir([folder basename '*ch*' '.tif']);

nfiles=numel(stack_dir);

Position_list=[];
size_list=[];
nch=0;
nz=[];

for n=1:nfiles
   name=stack_dir(n).name;
   delimiters=[regexp(name,'_') regexp(name,'.tif')];
   
   % time - hpp
   %poshpp=regexp(name,'hpp_');
   %postime=delimiters(find(delimiters<poshpp,1,'last'))+1;
   %hpp=name(postime:poshpp-1);
   %stack_dir(n).hpp=hpp;
      
%    % Mark and Find
%    posMF=regexp(name,'Mark_and_Find_')
%    if(~isempty(posMF))
%        posEndMF=delimiters(find(delimiters>posMF,4,'first'))-1;
%        posEndMF=posEndMF(end);
%        MF=str2num(name(posMF+14:posEndMF));
%        stack_dir(n).MF=MF;
%    else
%        stack_dir(n).MF=[];
%    end
% 
%    % time
    posT=regexp(name,'_t');
    if(~isempty(posT))
       posEndT=delimiters(find(delimiters>posT,1,'first'))-1;
       T=str2num(name(posT+2:posEndT));
       stack_dir(n).T=T;
       basenameMetadata=basename(1:posT-1);
    else
       stack_dir(n).T=[];
       basenameMetadata=basename;
       T=[];
    end
   
   % position - Position
   posPosition=regexp(name,'_Position');
   
   if(isempty(posPosition)) %if there is no string Position, then there is only 1 positions - called NaN
    Position=NaN;
   else
    posEndPos=delimiters(find(delimiters>posPosition,1,'first'))-1;
    Position=str2num(name(posPosition+9:posEndPos));
   end
   
   stack_dir(n).Position=Position; %save Position 
   whichPos=find(Position_list==Position | isnan(Position_list)); %we wonder whether the Position is already present and create it
   if(isempty(whichPos))
        Position_list=[Position_list; Position];
        imgsize=size(imread([folder,stack_dir(n).name]));
        size_list=vertcat(size_list,imgsize);
        whichPos=numel(Position_list);
        nz(whichPos)=0;
   end
   
   % channel - ch
   posch=regexp(name,'_ch');
   posEndch=delimiters(find(delimiters>posch,1,'first'))-1;
   ch=str2num(name(posch+3:posEndch));
   stack_dir(n).ch=ch;
   nch=max(ch+1,nch);
   
   % z
   posz=regexp(name,'_z');
   if(~isempty(posz))
        posEndz=delimiters(find(delimiters>posz,1,'first'))-1;
        z=str2num(name(posz+2:posEndz));
        stack_dir(n).z=z;
        nz(whichPos)=max(z+1,nz(whichPos));
   else
        stack_dir(n).z=0;
        nz(whichPos)=1;
   end
end


if(sum(isnan(Position_list))>1)
    display('Position_list is:');
    Position_list
    error('There is more than a NaN in the PositionList!');
end
    
if(isempty(Position_list))
    display('PositionList is empty!');
end    
    

% find metadata file
if(isnan(Position_list(1)))
    metadata_dir=dir(strcat([folder 'Metadata/'],[basenameMetadata],'.xml'));
    metadataPath{1}=[folder 'Metadata/' metadata_dir.name];
else
    for p=1:numel(Position_list)
        metadata_dir=dir(strcat([folder 'Metadata/'],[basenameMetadata '_Position' num2str(Position_list(p), '%03u')],'.xml'))
        metadataPath{p}=[folder 'Metadata/' metadata_dir.name];
        metadata{p}=parseLeicaMetadata(metadataPath{p});
    end
end


%update metadata timestamp with time
if(~isempty(T))
    if(~isfield(opts,'timeStep'))
       error('Please provide timeStep in opts (in s)'); 
    end
    for p=1:numel(metadata)
            metadata{p}.hTimeStep=opts.timeStep/3600; %conversion to h
            % we update the timestamp to take care of the time-shift in
            % successive time-points
            dt=metadata{p}.dt;
            metadata{p}.hTimeStamp=metadata{p}.hTimeStep*T+datenum(year(dt),month(dt),day(dt),hour(dt),minute(dt),second(dt))*24;
            metadata{p}.dt=datetime(year(dt),month(dt),day(dt),hour(dt),minute(dt),second(dt)+metadata{p}.hTimeStep*T*3600);
    end
end


mystack=[];
mystack.basename=basename;
mystack.dir=stack_dir;
mystack.metadata=metadata;
mystack.nz=nz;
mystack.nch=nch;
mystack.Position_list=Position_list;
mystack.size_list=size_list;
mystack.folder=folder;
mystack.metadata=metadata;
mystack.metadataPath=metadataPath;

end

