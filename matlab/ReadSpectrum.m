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



% Read spectrum written by latest versions of Snoopy (this is a beta
% version).

rep='../';

%rep='/Volumes/Idris_Workdir/snoopy/mri/768x384x192/Pm1/';

mean_start=1;

nspec=21;    % 9 for old spectrum (no transfer). 15 for recent versions, 21 for v5.0 version
spectruml=importdata([rep,'spectrum.dat']);
spectrum.k=spectruml(1,1:(end-1));
spectrum.n=spectruml(2,1:(end-1));
spectrum.vx=spectruml(3:nspec:end,2:end);
spectrum.vy=spectruml(4:nspec:end,2:end);
spectrum.vz=spectruml(5:nspec:end,2:end);
spectrum.bx=spectruml(6:nspec:end,2:end);
spectrum.by=spectruml(7:nspec:end,2:end);
spectrum.bz=spectruml(8:nspec:end,2:end);
spectrum.th=spectruml(9:nspec:end,2:end);
spectrum.vxvy=spectruml(10:nspec:end,2:end);
spectrum.bxby=spectruml(11:nspec:end,2:end);
spectrum.ad_vx=spectruml(12:nspec:end,2:end);
spectrum.ad_vy=spectruml(13:nspec:end,2:end);
spectrum.ad_vz=spectruml(14:nspec:end,2:end);
spectrum.tr_bx=spectruml(15:nspec:end,2:end);
spectrum.tr_by=spectruml(16:nspec:end,2:end);
spectrum.tr_bz=spectruml(17:nspec:end,2:end);
spectrum.tr_vx=spectruml(18:nspec:end,2:end);
spectrum.tr_vy=spectruml(19:nspec:end,2:end);
spectrum.tr_vz=spectruml(20:nspec:end,2:end);
spectrum.hel=spectruml(21:nspec:end,2:end)+spectruml(22:nspec:end,2:end)+spectruml(23:nspec:end,2:end);


spectrum.t=transpose(spectruml(4:nspec:end,1));

clear spectrum1;

% Compute the fluxes
flux.adk=-cumsum(mean(spectrum.ad_vx(mean_start:end,:)+spectrum.ad_vy(mean_start:end,:)+spectrum.ad_vz(mean_start:end,:)));
flux.adm=-cumsum(mean(spectrum.tr_vx(mean_start:end,:)+spectrum.tr_vy(mean_start:end,:)+spectrum.tr_vz(mean_start:end,:)+spectrum.tr_bx(mean_start:end,:)+spectrum.tr_by(mean_start:end,:)+spectrum.tr_bz(mean_start:end,:)));


figure(1)
loglog(spectrum.k,spectrum.n);
title('Number of modes');

figure(2)
loglog(spectrum.k,mean(spectrum.vx(mean_start:end,:)+spectrum.vy(mean_start:end,:)+spectrum.vz(mean_start:end,:)),spectrum.k,mean(spectrum.bx(mean_start:end,:)+spectrum.by(mean_start:end,:)+spectrum.bz(mean_start:end,:)),spectrum.k,1e-1*spectrum.k.^(-5/3));
title('Spectrum');
legend('Kinetic','Magnetic','K43')


figure(3)
semilogx(spectrum.k,mean(spectrum.vxvy(mean_start:end,:)),spectrum.k,mean(spectrum.bxby(mean_start:end,:)));

title('Transport');
legend('Reynolds','Maxwell')


figure(4)
loglog(spectrum.k,spectrum.vx(mean_start:end,:)+spectrum.vy(mean_start:end,:)+spectrum.vz(mean_start:end,:))
title('Kinetic spectrum in time')

figure(5)
loglog(spectrum.k,spectrum.bx(mean_start:end,:)+spectrum.by(mean_start:end,:)+spectrum.bz(mean_start:end,:))
title('Magnetic spectrum in time')

figure(6)
semilogx(spectrum.k,flux.adk,spectrum.k,flux.adm);
title('Fluxes');
legend('Kinetic transfer','Magnetic transfer')

figure(7)
loglog(spectrum.k,spectrum.k.*spectrum.k.*mean(spectrum.vx(mean_start:end,:)+spectrum.vy(mean_start:end,:)+spectrum.vz(mean_start:end,:)),spectrum.k,spectrum.k.*spectrum.k.*mean(spectrum.bx(mean_start:end,:)+spectrum.by(mean_start:end,:)+spectrum.bz(mean_start:end,:)),spectrum.k,1e-1*spectrum.k.^(1/3));
title('dissipation');
legend('Kinetic','Magnetic','K43')

figure(8)
semilogx(spectrum.k,spectrum.k.*mean(spectrum.hel(mean_start:end,:))),
title('Kinetic helicity');

figure(9)
semilogx(spectrum.k,spectrum.hel(mean_start:end,:));
title('Kinetic helicity in time')
