mrstModule add coarsegrid;

close all;

%% Set up model
 % Set up model geometry
[nx,ny,nz] = deal( 5,  5, 1);
[Dx,Dy,Dz] = deal(100, 100, 1);
grid = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
grid = computeGeometry(grid);

plotGrid(grid); view(3); axis tight

% Define rock model
rock = makeRock(grid, 30*milli*darcy, 0.3);


% Initiate complete model
model = initiateModel(grid, rock);

% 
% p = linspace(100*barsa,220*barsa,50)';
% s = linspace(0,1,50)';
% plot(p/barsa, model.rock.pv_r(1).*exp(model.rock.cr*(p-model.rock.p_r)),'LineWidth',2);


%% Plot model for two-phase compressible fluid
% Water phase and a a lighter, more viscous oil phase with different relative
% permeability function

% figure;
% plot(p/barsa, [model.water.rhoW(p), model.oil.rhoO(p)],'LineWidth',2);
% legend('Water density', 'Oil density')
% 
% figure;
% plot(p/barsa, [model.water.krW(s), model.oil.krO(s)],'LineWidth',2);
% legend('krW', 'krO')
% 
% spy(model.operator.C)

%% Impose vertical equilibrium
gravity reset on, g = norm(gravity);
[z_0, z_max] = deal(0, max(model.G.cells.centroids(:,3)));
equil  = ode23(@(z,p) g .* model.oil.rhoO(p), [z_0, z_max], model.rock.p_r);
p_init = reshape(deval(equil, model.G.cells.centroids(:,3)), [], 1);  clear equil
sW_init = zeros(model.G.cells.num, 1);

%% Initialize for solution loop
[p_ad, sW_ad] = initVariablesADI(p_init, sW_init);

numSteps = 100;                  % number of time-steps
totTime  = 365*day;             % total simulation time
dt       = totTime / numSteps;  % constant time step
tol      = 1e-5;                % Newton tolerance
maxits   = 20;                  % max number of Newton its

% Multigrid variables
% Presmoothing steps
v1_iter = 1;
% Postsmoothing steps
v2_iter = 5;


model.well.inRate = 1*sum(model.rock.pv(p_init))/totTime;
model.well.outRate = 0.5*model.well.inRate;

sol = repmat(struct('time',[],'pressure',[], 's', []),[numSteps+1,1]);
sol(1)  = struct('time', 0, 'pressure', double(p_ad), ...
    's', double(sW_ad));
figure(3)
%% Main loop
t = 0; step = 0;
hwb = waitbar(t,'Simulation ..');
while t < totTime
   t = t + dt;
   step = step + 1;
   fprintf('\nTime step %d: Time %.2f -> %.2f days\n', ...
      step, convertTo(t - dt, day), convertTo(t, day));

  % Multigrid
  [p_ad, sW_ad,nit] = multigridCycleV3(v1_iter,v2_iter,model,p_ad,sW_ad,tol,maxits,g,dt);
   
  if mod(step,10) == 0
    figure
    subplot(2, 1, 1); plot(model.G.cells.indexMap,p_ad.val);
    title('Pressure')
    
    subplot(2, 1, 2); plot(model.G.cells.indexMap,sW_ad.val);
    title('Saturation')
    drawnow
   end
   
   if nit > maxits
      error('Newton solves did not converge')
   else % store solution
      sol(step+1)  = struct('time', t, ...
                            'pressure', double(p_ad), ...
                            's', double(sW_ad));
      waitbar(t/totTime,hwb);
   end
end
close(hwb);

%% Plot pressure evolution

for i = 1:numSteps
    figure(1); clf
    subplot(2, 1, 1)
    plotCellData(grid, sol(i).pressure);
    title('Pressure')
    view(30, 40);
    subplot(2, 1, 2)
    plotCellData(grid, sol(i).s);
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