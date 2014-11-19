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

full_file=importdata([rep,'timevar']);
timevar=transpose(full_file.data);
nvar=size(timevar,1);
var_name=strread(full_file.textdata{1},'%s',nvar);

for i=1:nvar
   assignin('base',var_name{i},timevar(i,:)); 
end


%mag.helicity=timevar(22,:);

% Only when braginskii pressure is present
%brag.p=timevar(23,:);
%brag.stress=timevar(24,:);


%mode1.kx=timevar(23,:);
%mode1.wz=timevar(24,:);
%mode2.kx=timevar(25,:);
%mode2.wz=timevar(26,:);

%vel.shear=timevar(20,:);
%vel.vxmode=timevar(20,:);
%vel.vymode=timevar(21,:);
if exist('dvxvy') 
    vxvy=dvxvy;
end

if exist('vxvy')&&exist('bxby')

    alpha=vxvy-bxby;

    figure(8)
    plot(t,alpha,t,bxby,t,vxvy);
    legend('alpha coefficient','Maxwell Stress','Reynolds Stress');
    title('Turbulent Transport')
    
    clear mean_transport;
    
    if length(vxvy) > mean_start
        mean_vel_transport=cumsum(vxvy(mean_start:end))./(1:length(vxvy(mean_start:end)));
        mean_mag_transport=cumsum(bxby(mean_start:end))./(1:length(bxby(mean_start:end)));
        mean_tot_transport=cumsum(alpha(mean_start:end))./(1:length(alpha(mean_start:end)));

        e=mean(ev(mean_start:end))
        alpha_v=mean_vel_transport(end);
        alpha_b=mean_mag_transport(end);
        alpha=alpha_v-alpha_b
    else
        disp(['Caution, we are not removing initial conditions']);
        e=mean(ev)
        alpha_v=mean(vxvy);
        alpha_b=-mean(bxby);
        alpha=alpha_v-alpha_b
    end


elseif exist('vxvy')
    alpha=vxvy;

    figure(8)
    semilogy(t,alpha);
    title('Turbulent Transport')
    
    clear mean_transport;
    
    if length(vxvy) > mean_start
        mean_vel_transport=cumsum(vxvy(mean_start:end))./(1:length(vxvy(mean_start:end)));
        mean_tot_transport=cumsum(alpha(mean_start:end))./(1:length(alpha(mean_start:end)));

        e=mean(ev(mean_start:end))
        alpha_v=mean_vel_transport(end);
        alpha=alpha_v
    else
        disp(['Caution, we are not removing initial conditions']);
        e=mean(ev)
        alpha_v=mean(vxvy);
        alpha=alpha_v
    end
    
end


if exist('ev')
    figure(1)
    semilogy(t,ev)
    title('Kinetic energy')
end

if exist('em')
    figure(2)
    semilogy(t,em)
    title('Magnetic energy')
end

if exist('vxmax')&&exist('vxmin')&&exist('vymax')&&exist('vymin')&&exist('vzmax')&&exist('vzmin')
    figure(3)
    subplot(3,1,1), plot(t,vxmax,t,vxmin);
    title('Velocity field; extrema vx,vy,vz');
    subplot(3,1,2), plot(t,vymax,t,vymin);
    subplot(3,1,3), plot(t,vzmax,t,vzmin);
end

if exist('bxmax')&&exist('bxmin')&&exist('bymax')&&exist('bymin')&&exist('bzmax')&&exist('bzmin')
    figure(4)
    subplot(3,1,1), plot(t,bxmax,t,bxmin);
    title('Magnetic field; extrema bx,by,bz');
    subplot(3,1,2), plot(t,bymax,t,bymin);
    subplot(3,1,3), plot(t,bzmax,t,bzmin);
end


if exist('j2')&&exist('w2')
    figure(11)
    semilogy(t,w2,t,j2,'linewidth',2);
    legend('Enstrophy','J^2')
    xlabel('t','fontsize',16)
    set(gca,'fontsize',16)
    
elseif exist('w2')
    figure(11)
    semilogy(t,w2,'linewidth',2);
    legend('Enstrophy')
    xlabel('t','fontsize',16)
    set(gca,'fontsize',16)
end

if exist('hm')
    figure(12)
    plot(t,hm,'linewidth',2);
    title('Magnetic helicity')
    xlabel('t','fontsize',16)
    set(gca,'fontsize',16)
end


if exist('dmax')&&exist('dmin')
    figure(13)
    plot(t,dmax,t,dmin);
    title('Density field extrema');
end

if exist('az2')
    figure(14)
    plot(t,az2);
    title('potential vector')
end

if exist('vxaz')
        figure(15)
        plot(t,vxaz);
        title('potential vector source')
end

if exist('hc')
    figure(16)
    plot(t,hc);
    title('Cross helicity')
end
