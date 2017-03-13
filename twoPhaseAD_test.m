%%Attempt at combining code from singlePhaseAD.m and the MRST introduction
%%book
%%
gravity reset on, g = norm(gravity);

%% Set up model geometry
%cartDim = [128 128 1];
%domain = [250 250 20];
%grid      = computeGeometry(cartGrid(cartDim,domain));

[nx,ny,nz] = deal( 10,  10, 10);
[Dx,Dy,Dz] = deal(200, 200, 50);
%computeGeometry - struct: grid.{cell{volumes,centroids}, faces.{areas,normals,centroids}}
grid = computeGeometry(cartGrid([nx, ny, nz], [Dx, Dy, Dz]));

%% Define rock and fluid model - incompressible fluids + capillary pressure (last: how - here?)

%makeRock - struct.{perm, poro}
rock   = makeRock(grid, 100*milli*darcy, 0.2);
pv     = poreVolume(grid, rock);

%Initialize fluids - with simplified Corey model (p 306)
%init__Fluid: struct.{properties(),saturation(state),relperm(saturation,state)}
%[mu, rho] = fluid.properties(); %gives mu_w and mu_n etc
%fluid = initCoreyFluid('mu' , [ 1, 10]*centi*poise ,...
%                       'rho' , [1014, 859]*kilogram/meter^3,... 
%                       'n' , [ 3, 2.5] ,...
%                       'sr' , [ 0.2, .15] , ...
%                       'kwm', [ 1, .85]);

fluid = initSimpleFluidPc('mu' , [   1,  10]*centi*poise     , ...
                          'rho', [1014, 859]*kilogram/meter^3, ...
                          'n'  , [   2,   2], ...
                          'pc_scale', 2*barsa);

%% Compute transmissibilities
N  = double(grid.faces.neighbors);
intInx = all(N ~= 0, 2);
N  = N(intInx, :);                          % Interior neighbors
hT = computeTrans(grid, rock);                 % Half-transmissibilities
cf = grid.cells.faces(:,1);
nf = grid.faces.num;
T  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]); % Harmonic average
T  = T(intInx); 

%% Define discrete operators
n = size(N,1);
C = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[-1 1], n, grid.cells.num);
grad = @(x)C*x;
div  = @(x)-C'*x;
avg  = @(x) 0.5 * (x(N(:,1)) + x(N(:,2)));
spy(C) 

%% Define flow equations
gradz  = grad(grid.cells.centroids(:,3));
[mu,rho] = fluid.properties();
v      = @(p)  -(T/mu(1)).*( grad(p) - g*rho(1)*gradz );

%kr - phase relative permeability
kr = @(phase,realperm) realperm(:,phase);
%m - mobility
m = @(phase, state) kr(phase, fluid.realperm(fluid.saturation(state),state))/mu(phase);
%ps - phase saturation
s_phase = @(phase, saturation) saturation(:,phase); 
%phase pressure
pp = @(phase, pressure) pressure(:,phase);

%Well/source present in flow equations? -> Assume this is not the case as
%it is only required at the residual eq.

%saturationEq = @(s_ad,state,dt) (pv/dt)*(s_ad-s_phase(1,fluid.saturation(state))) + ...
%    grad(fluid.relperm(s_phase(1,fluid.saturation(state)),state) * v(pp(1,fluid.pressure(state))));

saturationEq = @(s_ad,s_0,dt) (pv/dt)*(s_ad-s_0) + ...
    grad(fluid.relperm(s_phase(1,fluid.saturation(state)),state) * v(pp(1,fluid.pressure(state))));



pressureEq = @(state) - grad((m(1,state)+ m(2,state)) * T * grad(pp(1,state))) ...
    - grad(m(2,state)* T * grad(fluid.pc(state))) - ...
    grad((m(1,state)*rho(1) + m(2,state)*rho(2)) * T * g*gradz);


%% Wells - secton 5.1.5 p.148,

% Add wells
wells = addWell([],grid, rock, 1, 'Type', 'bhp', ...
    'Val', 100*barsa, 'name', 'I', 'radius', .1, 'Comp_i', [1 0]);
wells = addWell(wells,grid, rock, grid.cells.num, 'Type', 'bhp', ...
    'Val', 0, 'name', 'P', 'radius', .1, 'Comp_i', [0 1]);

% Wells form tuorial_two_phase.m
%rate = 0.5*meter^3/day;
%bhp  = 1*barsa;
%W = verticalWell([], G, rock, 1, 1, 1:nz,          ...
%                  'Type', 'rate', 'Val', rate, ...
%                  'Radius', .1, 'Name', 'I', 'Comp_i', [1 0]);
% W = verticalWell(W, G, rock, nx, ny, 1:nz,     ...
%                  'Type','bhp', 'Val', bhp, ...
%                  'Radius', .1, 'Dir', 'x', 'Name', 'P', 'Comp_i', [0 1]);

%Define well equations
well_cells = wells(1).cells; % connection grid cells
wellsI = wells(1).WI;    % well-indices
dz = wells(1).dZ;    % connection depth relative to bottom-hole

p_conn  = @(bhp)  bhp + g*dz.*rho; %connection pressures
q_conn  = @(p,bhp) wellsI .* (rho / mu) .* (p_conn(bhp) - p(well_cells));

rateEq = @(p,bhp,qS)  qS-sum(q_conn(p, bhp))/rhoS;
ctrlEq = @(bhp)       bhp-100*barsa;

%% Initialize state of the simnulation
%initState: - struct.{pressure, flux, s, facePressure}
state  = initState(grid,wells,100*barsa, [0 1]);

[s_ad,p_ad,bhp_ad,qS_ad] = initVariablesADI(state.s(1),state.pressure,state.pressure(well_cells),0);

%%%%
day = 86400;
numSteps = 52;                  % number of time-steps
totTime  = 365*day;             % total simulation time
dt       = totTime / numSteps;  % constant time step
tol      = 1e-5;                % Newton tolerance
maxits   = 10;                  % max number of Newton its

res_norm = 1;

sol = repmat(struct('time',[],'saturation',[],'pressure',[],'bhp',[],'qS',[]),[numSteps+1,1]);
sol(1)  = struct('time', 0,'saturation',double(s_ad), 'pressure', double(p_ad), ...
   'bhp', double(bhp_ad), 'qS', double(qS_ad));
s_0 = s_ad;
%% Main loop
time = 0; iteration = 0;
hwb = waitbar(time,'Simulation ..');
while time < totTime
    time = time+dt;
    iteration = iteration + 1;
    
    % Newton Transport looop
    res_norm = 1e99;
    init_pres = double(p_ad);
    newtonits = 0;
    while(tol < res_norm) && (maxits > iteration)
        s_eq = saturationEq(s_ad.jac,s_0,dt); 
        %s_eq(well_cells)  = saturationEq(well_cells) - q_conn(saturation_ad, bhp_ad);
        
        %s_eqs = {s_eq, rateEq(s_ad, bhp_ad, qS_ad), ctrlEq(bhp_ad)}; 
        
        %saturation_eq = cat(eq{:});
        J  = s_eq.jac{1};
        res = s_eq.val;
        upd = -(J\res);
        
        % Update variables
        s_ad.val   = s_ad.val   + upd(pIx);
        %bhp_ad.val = bhp_ad.val + upd(bhpIx);
        %qS_ad.val  = qS_ad.val  + upd(qSIx);
      
        res_norm = norm(res);
        newtonits = newtonits + 1;
        
        if iteration > maxits,
            error('Newton transport solves did not converge')
        end
    end
    % Newton pressure loop
    res_norm = 1e99;
    while(tol < res_norm) && (maxits > iteration)
        p_eq = presEq(p_ad,state); 
        pressureEq(well_cells)  = pressureEq(well_cells) - q_conn(p_ad, bhp_ad);
        
        pressure_eqs = {pressureEq, rateEq(p_ad, bhp_ad, qS_ad), ctrlEq(bhp_ad)}; 
        
        pressureEq = cat(eq{:});
        J  = eq.jac{1};
        res = eq.val;
        upd = -(J\res);
        
        % Update variables
        p_ad.val   = p_ad.val   + upd(pIx);
        bhp_ad.val = bhp_ad.val + upd(bhpIx);
        qS_ad.val  = qS_ad.val  + upd(qSIx);

        res_norm = norm(res);
        newtonits = newtonits + 1;
        
    if iteration > maxits,
          error('Newton pressure solves did not converge')
    else % store solution
        sol(step+1)  = struct('time', time,'saturation',double(s_ad),...
            'pressure', double(p_ad), 'bhp', double(bhp_ad), 'qS', double(qS_ad));
        waitbar(time/totTime,hwb);
    end
    end
end
close(hwb);


%% Plot production rate and pressure decay
clf,
[ha,hr,hp] = plotyy(...
   [sol(2:end).time]/day, -[sol(2:end).qS]*day, ...
   [sol(2:end).time]/day, mean([sol(2:end).pressure]/barsa), 'stairs', 'plot');
set(ha,'FontSize',16);
set(hr,'LineWidth', 2);
set(hp,'LineStyle','none','Marker','o','LineWidth', 1);
set(ha(2),'YLim',[100 210],'YTick',100:50:200);
xlabel('time [days]');
ylabel(ha(1), 'rate [m^3/day]');
ylabel(ha(2), 'avg pressure [bar]');

%% Plot pressure evolution
clf;
steps = [2 5 10 20];
for i=1:4
   subplot(2,2,i);
   plotCellData(G, sol(steps(i)).pressure/barsa, show,'EdgeColor',.5*[1 1 1]);
   plotWell(G,W);
   view(-125,20), camproj perspective
   caxis([115 205]);
   axis tight off; zoom(1.4)
   text(200,170,-8,[num2str(round(steps(i)*dt/day)) ' days'],'FontSize',14);
end
h=colorbar('horiz','Position',[.1 .05 .8 .025]);
colormap(jet(55));
