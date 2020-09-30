function [outStack] = equalizeStack(s)
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
    
    eq=zeros(size(s));

    for z=1:size(s,3)

        img=im2double(s(:,:,z));

        if(std(img)==0)
            img(:,:)=0;
        else
            img=adapthisteq(img,'NumTiles',[10 10]);
        end
        eq(:,:,z)=img;

    end

    outStack=im2uint8(eq);

end