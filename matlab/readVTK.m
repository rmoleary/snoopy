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

function V = readVTK(vtkfile)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Usage: V = readVTK(vtkfile)
%
%   V:       Structure in which the data are stored
%   vtkfile: The filename
%   notes:   Only reads binary STRUCTURED_POINTS
%
% Erik Vidholm 2006
% Geoffroy Lesur 2009
%       Extended to include several fields in the same file (as Snoopy
%       does)
%       The output is now a structure.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% open file (OBS! big endian format)
fid = fopen(vtkfile,'r','b');

if( fid == -1 )
  disp('File not found');
  return
end

fgetl(fid); % # vtk DataFile Version x.x
fgetl(fid); % comments
fgetl(fid); % BINARY
fgetl(fid); % DATASET STRUCTURED_POINTS

s = fgetl(fid); % DIMENSIONS NX NY NZ
sz = sscanf(s, '%*s%d%d%d').';

s=fgetl(fid); % ORIGIN OX OY OZ
origin = sscanf(s, '%*s%f%f%f').';
s=fgetl(fid); % SPACING SX SY SZ
spacing = sscanf(s, '%*s%f%f%f').';

% Generate coordinates

for i = 1:sz(1)
    x(i)=origin(1)+(i-1)*spacing(1);
end
for i = 1:sz(2)
    y(i)=origin(2)+(i-1)*spacing(2);
end
for i = 1:sz(3)
    z(i)=origin(3)+(i-1)*spacing(3);
end

%Inialize the output structure
V = struct('x',x,'y',y,'z',z);
V.nx = sz(1);
V.ny = sz(2);
V.nz = sz(3);

fgetl(fid); % POINT_DATA NXNYNZ
s = fgetl(fid); % SCALARS/VECTORS name data_type (ex: SCALARS imagedata unsigned_char)
varname = sscanf(s, '%*s%s%*s');

% The first one is a scalar...
fgetl(fid); % the lookup table

% read data
Q = fread(fid,prod(sz),'*single');
V = setfield(V,varname,reshape(Q,sz));

%Let's now fight with the FIELD region for the other components.
s=fgetl(fid); %FIELD Fieldname num_field
num_field = sscanf(s, '%*s%*s%d');

for i = 1:num_field     % Loop on all the remaining fields.
    s = fgetl(fid);     % fieldname dimensionlity num_point data_type
    varname = sscanf(s, '%s',1);
    % read data
    Q = fread(fid,prod(sz),'*single');
    V = setfield(V,varname,reshape(Q,sz));
end

fclose(fid);