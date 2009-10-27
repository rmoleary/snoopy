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


% Example File: how to read the timevar file in Matlab
% This example reads all the variable in the timevar file
% and display the time history of several quantities...

rep='../';


mean_start=100;     % Where do we start our time-averages?
ncorr=20;

timevar=(transpose(importdata([rep,'timevar'])));
t=timevar(1,:);
vel.e=timevar(2,:);
mag.e=timevar(3,:);
vel.vxmax=timevar(4,:);
vel.vxmin=timevar(5,:);
vel.vymax=timevar(6,:);
vel.vymin=timevar(7,:);
vel.vzmax=timevar(8,:);
vel.vzmin=timevar(9,:);
vel.vxvy=timevar(10,:);
mag.bxmax=timevar(11,:);
mag.bxmin=timevar(12,:);
mag.bymax=timevar(13,:);
mag.bymin=timevar(14,:);
mag.bzmax=timevar(15,:);
mag.bzmin=timevar(16,:);
mag.bxby=timevar(17,:);
vel.thmax=timevar(18,:);
vel.thmin=timevar(19,:);
vel.vort=timevar(20,:);
mag.curr=timevar(21,:);
mag.helicity=timevar(22,:);


budget.alpha=vel.vxvy-mag.bxby;

clear mean_transport;
if length(vel.vxvy) > mean_start
    mean_vel_transport=cumsum(vel.vxvy(mean_start:end))./(1:length(vel.vxvy(mean_start:end)));
    mean_mag_transport=cumsum(mag.bxby(mean_start:end))./(1:length(mag.bxby(mean_start:end)));
    mean_tot_transport=cumsum(budget.alpha(mean_start:end))./(1:length(budget.alpha(mean_start:end)));

    e=mean(vel.e(mean_start:end))
    alpha_v=mean_vel_transport(end)
    alpha_b=mean_mag_transport(end)
    alpha=alpha_v-alpha_b
else
    disp(['Caution, we are not removing initial conditions']);
    e=mean(vel.e)
    alpha_v=mean(vel.vxvy);
    alpha_b=-mean(mag.bxby);
    alpha=alpha_v-alpha_b
end

% Compute the correlation time of alpha
clear corr

alphanorm=budget.alpha-mean(budget.alpha);
for i=1:ncorr
    j=i-1;
    corr(i)=mean(alphanorm(1:end-j).*alphanorm(1+j:end))./mean(alphanorm.^2);
end
corr_time=sum(corr)*(t(2)-t(1))

plot(corr);
% compute transport over 100 time steps (allow statistical computations for
% transport vs Romega)

figure(1)
subplot(2,1,1),  plot(t,log10(vel.e))
title('Velocity field: energy, Reynolds stress (log scale)');
subplot(2,1,2),  plot(t,log10(abs(vel.vxvy)))

figure(2)
subplot(3,1,1), plot(t,vel.vxmax,t,vel.vxmin);
title('Velocity field; extrema vx,vy,vz');
subplot(3,1,2), plot(t,vel.vymax,t,vel.vymin);
subplot(3,1,3), plot(t,vel.vzmax,t,vel.vzmin);

figure(3)
subplot(2,1,1),  plot(t,log10(mag.e))
title('Magnetic field: energy, Maxwell stress (log scale)');
subplot(2,1,2),  plot(t,log10(abs(mag.bxby)))

figure(4)
subplot(3,1,1), plot(t,mag.bxmax,t,mag.bxmin);
title('Magnetic field: extrema bx,by,bz');
subplot(3,1,2), plot(t,mag.bymax,t,mag.bymin);
subplot(3,1,3), plot(t,mag.bzmax,t,mag.bzmin);

figure(5)

plot(t(mean_start:end),-mean_vel_transport,t(mean_start:end),mean_mag_transport,t(mean_start:end),mean_tot_transport);
legend('Reynolds Stress','Maxwell Stress','alpha coefficient');
title('Transport Cumulative average');



figure(6)
subplot(3,1,1),  semilogy(t,vel.e)
title('Energy: Kinetic energy, magnetic energy, total turbulent transport');
subplot(3,1,2),  semilogy(t,mag.e)
subplot(3,1,3),  semilogy(t,abs(budget.alpha))


figure(7)
plot(t,budget.alpha,t,mag.bxby,t,vel.vxvy);
legend('alpha coefficient','Maxwell Stress','Reynolds Stress');
title('Turbulent Transport')

figure(8)
plot(t,budget.alpha,t,mag.bxby,'linewidth',2);
xlabel('t','fontsize',16)
ylabel('\alpha','fontsize',20);
set(gca,'fontsize',16)

figure(9)
semilogy(t,vel.vort,t,mag.curr,'linewidth',2);
legend('Enstrophy','J^2')
xlabel('t','fontsize',16)
set(gca,'fontsize',16)


figure(10)
plot(t,mag.helicity,'linewidth',2);
title('Magnetic helicity')
xlabel('t','fontsize',16)
set(gca,'fontsize',16)


