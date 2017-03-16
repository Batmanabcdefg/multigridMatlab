mrstModule add coarsegrid;

close all;

%% Multigrid variables
% Presmoothing steps
v1_its = 1;
% Postsmoothing steps
v2_its = 2;

%% Set up model geometry
[nx,ny,nz] = deal( 10,  10, 10);
[Dx,Dy,Dz] = deal(200, 200, 50);
G = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
G = computeGeometry(G);

plotGrid(G); view(3); axis tight

%% Define rock model
rock = makeRock(G, 30*milli*darcy, 0.3);

cr   = 1e-6/barsa;
p_r  = 200*barsa;
pv_r = poreVolume(G, rock);
pv   = @(p) pv_r .* exp( cr * (p - p_r) );

p = linspace(100*barsa,220*barsa,50)';
s = linspace(0,1,50)';
plot(p/barsa, pv_r(1).*exp(cr*(p-p_r)),'LineWidth',2);


%% Define model for two-phase compressible fluid
% Define a water phase
muW    = 1*centi*poise;
cw     = 1e-5/barsa;
rho_rw = 960*kilogram/meter^3;
rhoWS  = 1000*kilogram/meter^3;
rhoW   = @(p) rho_rw .* exp( cw * (p - p_r) );
krW = @(S) S.^2;

% Define a lighter, more viscous oil phase with different relative
% permeability function
muO   = 5*centi*poise;
co      = 1e-4/barsa;
rho_ro = 850*kilogram/meter^3;
rhoOS  = 750*kilogram/meter^3;
krO = @(S) S.^3;

rhoO   = @(p) rho_ro .* exp( co * (p - p_r) );
figure;
plot(p/barsa, [rhoW(p), rhoO(p)],'LineWidth',2);
legend('Water density', 'Oil density')

figure;
plot(p/barsa, [krW(s), krO(s)],'LineWidth',2);
legend('krW', 'krO')
%% Impose vertical equilibrium
gravity reset on, g = norm(gravity);
[z_0, z_max] = deal(0, max(G.cells.centroids(:,3)));
equil  = ode23(@(z,p) g .* rhoO(p), [z_0, z_max], p_r);
p_init = reshape(deval(equil, G.cells.centroids(:,3)), [], 1);  clear equil
sW_init = zeros(G.cells.num, 1);
%% Compute transmissibilities
N  = double(G.faces.neighbors);
intInx = all(N ~= 0, 2);
N  = N(intInx, :);                          % Interior neighbors
hT = computeTrans(G, rock);                 % Half-transmissibilities
cf = G.cells.faces(:,1);
nf = G.faces.num;
T  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]); % Harmonic average
T  = T(intInx);                             % Restricted to interior

%% Define discrete operators
n = size(N,1);
C = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[-1 1], n, G.cells.num);
grad = @(x)C*x; % Discrete gradient
div  = @(x)-C'*x; % Discrete divergence
avg  = @(x) 0.5 * (x(N(:,1)) + x(N(:,2))); % Averaging
upw = @(flag, x) flag.*x(N(:, 1)) + ~flag.*x(N(:, 2)); % Upwinding
spy(C)

gradz  = grad(G.cells.centroids(:,3));
%% Initialize for solution loop
[p_ad, sW_ad] = initVariablesADI(p_init, sW_init);
nc = G.cells.num;
pIx = 1:nc;
sIx = (nc+1):(2*nc);

numSteps = 100;                  % number of time-steps
totTime  = 365*day;             % total simulation time
dt       = totTime / numSteps;  % constant time step
tol      = 1e-5;                % Newton tolerance
maxits   = 15;                  % max number of Newton its

injIndex = 1;
prodIndex = G.cells.num;

inRate = 1*sum(pv(p_init))/totTime;
outRate = 0.5*inRate;

sol = repmat(struct('time',[],'pressure',[], 's', []),[numSteps+1,1]);
sol(1)  = struct('time', 0, 'pressure', double(p_ad), ...
    's', double(sW_ad));

%% Main loop
t = 0; step = 0;
hwb = waitbar(t,'Simulation ..');
while t < totTime
   t = t + dt;
   step = step + 1;
   fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
      step, convertTo(t - dt, day), convertTo(t, day));

  % Multigrid
  [p_ad, sW_ad,nit] = multigridV1(v1_its,v2_its,G, rock,p_ad,sW_ad,tol,maxits,rhoW,rhoO,grad,gradz,pv,krW,krO,muW,muO,T,g,avg,upw,dt,div,injIndex,inRate,rhoWS,prodIndex,pIx,sIx,cr,p_r);

  
  % Newton loop
  % [p_ad,sW_ad,res,nit] = newtonAD(p_ad,sW_ad,tol,maxits,rhoW,rhoO,grad,gradz,pv,krW,krO,muW,muO,T,g,avg,upw,dt,div,injIndex,inRate,rhoWS,prodIndex,pIx,sIx);

   if nit > maxits
      error('Newton solves did not converge')
   else % store solution
      sol(step+1)  = struct('time', t, ...
                            'pressure', double(p_approx), ...
                            's', double(sW_approx));
      waitbar(t/totTime,hwb);
   end
end
close(hwb);

%% Plot pressure evolution

for i = 1:numSteps
    figure(1); clf
    subplot(2, 1, 1)
    plotCellData(G, sol(i).pressure);
    title('Pressure')
    view(30, 40);
    subplot(2, 1, 2)
    plotCellData(G, sol(i).s);
    caxis([0, 1])
    view(30, 40);
    title('Watersaturation')
    drawnow
end

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}