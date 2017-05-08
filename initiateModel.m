function [model] = initiateModel(grid,varargin)
  %% Function description
  %
  % PARAMETERS:
  % G        - initialized grid structure
  % rock     - initialized rock structure
  % # Not currently used: varargin - {1} coarse/true Contains a true statement if model is to be
  %             initialized for a coarsed grid: Do not compute transmissibilities
  %
  % RETURNS:
  % p_ad_coarse - Coarsed version of the pressure eqs. with reinitialized
  %               ADI structure
  % sW_ad_coarse - Coarsed version of the pressure eqs. with reinitialized
  %               ADI structure
  % defect_p_ad_coarse - Coarsed version of the pressure defect
  % defect_s_ad_coarse - Coarsed version of the saturation defect
  % pIx_coarse   - Coarsed version of the pressure index
  % sIx_coarse   - Coarsed version of the saturation index
  %
  % COMMENTS:
  % - This coarsening function is currently not optimized for performance
  % - The function may be bugged in current state
  %
  % SEE ALSO:
  %

  %% Define wells
  injIndex = 1;
  prodIndex = grid.cells.num;

  inRate = 1;
  outRate = 0.5;
  
  well = struct('injIndex',injIndex, 'prodIndex',prodIndex, 'inRate', inRate, 'outRate', outRate);
  
  %% Rock model
  % Call to makeRock(grid, permeability, porosity) to initiate rock model
  % permeability range: {poor: 1-15, moderate: 15-20, good: 50-250, very
  % good: 250-1000
  % porosity range: {fair: 0.25, very low: 0.1}
  rock = makeRock(grid, 30*milli*darcy, 0.25);
  % Compressibility: normally in the range of 10^-6 to 10^-7, assumed to be
  % constant
  cr   = 1e-6/barsa;
  % Reference pressure
  p_r  = 200*barsa;
  pv_r = poreVolume(grid, rock);
  pv   = @(p) pv_r .* exp( cr * (p - p_r) );

  if(strcmp(grid.type{1},'generateCoarseGrid'))
      old_model = varargin{1};
      %pv_r = accumarray(partition,old_model.rock.pv_r);
      
      rock.perm(injIndex) = old_model.rock.perm(old_model.well.injIndex);
      rock.perm(prodIndex) = old_model.rock.perm(old_model.well.prodIndex);
      
      rock.poro(injIndex) = old_model.rock.poro(old_model.well.injIndex);
      rock.poro(prodIndex) = old_model.rock.poro(old_model.well.prodIndex);
          pv_r(injIndex) = old_model.rock.pv_r(old_model.well.injIndex);
      pv_r(prodIndex) = old_model.rock.pv_r(old_model.well.prodIndex);
  end
  %rock struct
  rock = struct('perm',rock.perm,'poro',rock.poro, ...
      'cr', cr, 'p_r',p_r, 'pv_r', pv_r, 'pv',pv);
  
  %% Define model for two-phase compressible fluid
  % Define a water phase
  muW    = 1*centi*poise;
  cw     = 1e-5/barsa;
  rho_rw = 960*kilogram/meter^3;
  rhoWS  = 1000*kilogram/meter^3;
  rhoW   = @(p) rho_rw .* exp( cw * (p - p_r) );
  krW = @(S) S.^2;
  
  water = struct('muW', muW, 'cw', cw, 'rho_rw', rho_rw, 'rhoWS', rhoWS, 'rhoW', rhoW, 'krW', krW);
    
  % Define a lighter, more viscous oil phase with different relative
  % permeability function
  muO   = 5*centi*poise;
  % Compressibility range: {slighly: 10^-5 to 10^-6, compressible: 10^-3 to
  % 10^-4}psi^-1
  co      = 1e-3/barsa; %1e-4
  rho_ro = 1050*kilogram/meter^3; % 850
  rhoOS  = 750*kilogram/meter^3; % 750
  krO = @(S) S.^3;

  rhoO   = @(p) rho_ro .* exp( co * (p - p_r) );

  oil = struct('muO', muO, 'co', co, 'rho_ro', rho_ro, 'rhoOS', rhoOS, 'rhoO', rhoO, 'krO', krO);
  
  
  %% Compute transmissibilities
  N  = double(grid.faces.neighbors);
  intInx = all(N ~= 0, 2);
  N  = N(intInx, :);                          % Interior neighbors
  hT = computeTrans(grid, rock);                 % Half-transmissibilities
  cf = grid.cells.faces(:,1);
  nf = grid.faces.num;
  T  = 1 ./ accumarray(cf, 1 ./ hT, [nf, 1]); % Harmonic average
  T  = T(intInx);                             % Restricted to interior
  
  %% Define discrete operators
  n = size(N,1);
  C = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[-1 1], n, grid.cells.num);
  grad = @(x)C*x; % Discrete gradient
  div  = @(x)-C'*x; % Discrete divergence
  avg  = @(x) 0.5 * (x(N(:,1)) + x(N(:,2))); % Averaging
  upw = @(flag, x) flag.*x(N(:, 1)) + ~flag.*x(N(:, 2)); % Upwinding 

  gradz  = grad(grid.cells.centroids(:,3));

  operator = struct('grad', grad, 'div', div, 'avg', avg, 'upw', upw, 'gradz', gradz, 'C',C);
  
  %% Remaining variables
  %Note: Write a better description
  nc = grid.cells.num;
  pIx = 1:nc;
  sIx = (nc+1):(2*nc);
  p_ad = 0;
  sW_ad = 0;
  gravity reset on, g = norm(gravity);
  
  
  %% Place all model parts and help function i a "modelstruct"
  model = struct('grid',grid,'rock', rock, 'water', water, 'oil',oil, 'T', T, ...
      'operator', operator, 'well', well, 'pIx', pIx,'sIx',sIx,'p_ad',p_ad,'sW_ad',sW_ad,'g',g);

end