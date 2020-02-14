function results=dsmceq_two_JAbejo(NumPart,Nstep)
% dsmceq - Dilute gas simulation using DSMC algorithm
% This version illustrates the approach to equilibrium
% Envoke as: dsmceq_two_JAbejo(3000,50) or dsmceq_JAbejo(3000,5)
% Edited code by Jonathan Abejo (214865273), USING MATLAB R2018b

clc, close all; help dsmceq_two_JAbejo;   % Clear memory and print header

%* Initialize constants  (particle mass, diameter, etc.)
boltz = 1.3806e-23;    % Boltzmann's constant (J/K)
mass = 6.63e-26;       % Mass of argon atom (kg)
diam = 3.66e-10;       % Effective diameter of argon atom (m)
T = 273;               % Temperature (K)
density = 1.78;        % Density of argon at STP (kg/m^3)
L = 1e-6;              % System size is one micron
npart = NumPart;
eff_num = density/mass*L^3/npart;
fprintf('Each particle represents %g atoms\n',eff_num);

%* Assign random positions and velocities to particles
rand('state',0);       % Initialize random number generator
x = L*rand(npart,1);   % Assign random positions
y = L*rand(npart,1);
z = L*rand(npart,1);
V = (max(x)-min(x))*(max(y)-min(y))*(max(z)-min(z));
A = V^(1/3);
P = npart*boltz*T/V;
fprintf('The expected pressure is %2.4g Pascals\n',P);
figure(10);
plot3(x,y,z, '*')
v_init = sqrt(3*boltz*T/mass);    % Initial speed
v_mp = sqrt(2*boltz*T/mass);        %v_mp uses most probably velocity (expected)
v = zeros(npart,3);    % Only x-component is non-zero
v2 = zeros(npart,3);
v(:,1) = v_init * (1 - 2*floor(2*rand(npart,1)));
v2(:,1) = v_mp * (1 - 2*floor(2*rand(npart,1)));

%* Plot the initial speed distribution
figure(1);  clf;
vmag = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
vbin = 50:100:1050;    % Bins for histogram
hist(vmag,vbin);  title('Measured Initial speed distribution');
xlabel('Speed (m/s)');  ylabel('Number');

figure(2); clf;
vmag2 = sqrt(v2(:,1).^2 + v2(:,2).^2 + v2(:,3).^2);
vbin2 = 50:100:1050;
hist(vmag2,vbin2);  title('Expected Initial speed distribution');
xlabel('Speed (m/s)');  ylabel('Number');


%* Initialize variables used for evaluating collisions
ncell = 15;                     % Number of cells
tau = 0.2*(L/ncell)/v_init;     % Set timestep tau
tau2 = 0.2*(L/ncell)/v_mp; 
vrmax = 3*v_init*ones(ncell,1); % Estimated max rel. speed
vrmax2 = 3*v_mp*ones(ncell,1); 
selxtra = zeros(ncell,1);       % Used by routine "colider"
selxtra2 = zeros(ncell,1); 
coeff = 0.5*eff_num*pi*diam^2*tau/(L^3/ncell);
coeff2 = 0.5*eff_num*pi*diam^2*tau2/(L^3/ncell);

coltot = 0;                     % Count total collisions
coltot2 = 0;

%* Declare structure for lists used in sorting
sortData = struct('ncell',ncell,'npart',npart,'cell_n',zeros(ncell,1),'index',zeros(ncell,1),'Xref',zeros(npart,1));  
sortData2 = struct('ncell',ncell,'npart',npart,'cell_n',zeros(ncell,1),'index',zeros(ncell,1),'Xref',zeros(npart,1));  

%* Loop for the desired number of time steps
nstep = Nstep;
for istep = 1:nstep
	
  %* Move all the particles ballistically
  x(:) = x(:) + v(:,1)*tau; % Update x position of particle
  x = rem(x+L,L);           % Periodic boundary conditions

  %* Sort the particles into cells
  sortData = sorter(x,L,sortData);
  sortData2 = sorter(x,L,sortData2);
  
  %* Evaluate collisions among the particles
  [v, vrmax, selxtra, col] = ...
                colider(v,vrmax,selxtra,coeff,sortData);
  coltot = coltot + col;
  
  [v2, vrmax2, selxtra2, col2] = ...
                colider(v2,vrmax2,selxtra2,coeff2,sortData2);
  coltot2 = coltot2 + col2;
  
  %* Periodically display the current progress
  if( rem(istep,10) < 1 )
    figure(7); clf;
    vmag = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
    hist(vmag,vbin);
    title(sprintf('Measured: Done %g of %g steps; %g collisions',...
                                      istep,nstep,coltot));
    xlabel('Speed (m/s)');  ylabel('Number');
    drawnow;
    figure(6); clf;
    vmag2 = sqrt(v2(:,1).^2 + v2(:,2).^2 + v2(:,3).^2);
    hist(vmag2,vbin2);
    title(sprintf('Expected: Done %g of %g steps; %g collisions',...
                                      istep,nstep,coltot2));
    xlabel('Speed (m/s)');  ylabel('Number');
    drawnow;

  end
  
  
end

%* Plot the histogram of the final speed distribution

%{
figure(3); clf;
vmag = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
hist(vmag,vbin);  
title(sprintf('Final distrib., Time = %g sec.',nstep*tau));
xlabel('Speed (m/s)'); ylabel('Number');
%}

figure(4); clf;
vmag2 = sqrt(v2(:,1).^2 + v2(:,2).^2 + v2(:,3).^2);
Force = mass*sum(vmag2)/tau;
Press = Force/A;
fprintf('The measured pressure is %2.4g Pascals\n',Press);
histogram(vmag2,vbin2,'DisplayStyle','bar','FaceColor','none','EdgeColor',[0 0.4470 0.7410]); 
hold on
histogram(vmag,vbin,'DisplayStyle','bar','FaceColor','none','EdgeColor',[0.8500 0.3250 0.0980]);  
hold off
title('Expected and Measured Equilibrium Speed Distributions');
legend('Measured', 'Expected')
%title(sprintf('Final distrib., Time = %g sec.',nstep*tau2));
xlabel('Speed (m/s)'); ylabel('Number');
end