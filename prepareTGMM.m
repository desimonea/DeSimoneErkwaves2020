function [] = prepareTGMM(fname,config_tpl,configfolder,datasource,datapathmac,respathmac,datapathwin,respathwin)
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
    basename=fname(1:end-4);

    copyfile([datasource filesep fname], [datapathmac filesep fname(1:end-4) '_t00.tif']);
    copyfile([datasource filesep fname], [datapathmac filesep fname(1:end-4) '_t01.tif']);

        
    folderpath=[respathmac basename];
    mkdir(folderpath);
    
    fin  = fopen(config_tpl,'r');
    fout = fopen([configfolder filesep  basename '_config.txt'],'w');
        
    idk=0;
    while ~feof(fin)
    
        idk=idk+1;
        s = fgetl(fin);

        if(idk==10)
            s=['imgFilePattern=' datapathwin basename '_t??'];
        end

        if(idk==13)
            s=['debugPathPrefix=' respathwin basename '\'];
        end

        fprintf(fout,'%s\n',s);
            
     end
  
     fclose(fin);
     fclose(fout);    
        
end
    

