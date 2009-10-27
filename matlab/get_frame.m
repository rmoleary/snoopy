%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
%	This file is part of the Snoopy code.

%   Snoopy code is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    Snoopy code is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.

%    You should have received a copy of the GNU General Public License
%    along with Snoopy code.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Example File: how to read VTK files in Matlab
% This script actually calls several subroutine which
% take care of the I/O and display. 

% In this example, the data are loaded in the structure
% V, which contains the fields and the coordinate system.
% Displays are done throuth the plotXY, plotXZ and plotYZ
% routines

%SNOOPY ROOT directory

rep='../';

% frame number
number=0001;

%Read VTK files
V=readVTK([rep,'data/v',num2str(number,'%0.4d'),'.vtk']);

plotYZ(V.vx, V.y, V.z, 1, 1);
title('Vx','fontsize',18);

plotYZ(V.vy, V.y, V.z, 2, 1);
title('Vy','fontsize',18);

plotXY(V.vz,V.x, V.y, 3, 1) 
title('Vz','fontsize',18);
