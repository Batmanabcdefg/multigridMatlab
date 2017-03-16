function [CG,coarse_rock,coarse_trans, coarse_p_ad, coarse_sW_ad,coarse_defect_p_ad,coarse_defect_s_ad,N, ...
    pv_coarse, pIx_coarse, sIx_coarse, grad_coarse, div_coarse, avg_coarse, upw_coarse, prodIndex_coarse, gradz_coarse] ...
    = coarsening(G, coarse_rock, trans, p_ad, sW_ad,defect,cr,p_r)

  %%Set up of coarse grid
  coarse_dims = ceil(G.cartDims/2);
  partition  = partitionCartGrid(G.cartDims,coarse_dims);
  CG = generateCoarseGrid(G, partition);
  CG = coarsenGeometry(CG);
  
  %Coarsen discrete operator matrix
  N  = double(CG.faces.neighbors);
  intInx = all(N ~= 0, 2);
  N  = N(intInx, :);
  
  n = size(N,1);
  C_coarse = sparse( [(1:n)'; (1:n)'], N, ones(n,1)*[-1 1], n, CG.cells.num);

  % Add fields to the coarse grid to ensure that it passes as a
  % regular grid for our purposes.
  CG.cartDims = coarse_dims;
  
  %% Restrict variables and model
  weighting = accumarray(partition,1);
  
  %Upscale permeability and transmissibility
  %rock.perm = upscalePerm(G,CG, rock,'T',[]); % -> Not working properly
  %rock.poro = accumarray(partition,rock.poro)./weighting;
  
  coarse_rock = makeRock(CG, 30*milli*darcy, 0.3);

  coarse_trans = computeTrans(CG, coarse_rock);  % upscaleTrans(G,CG,trans,'fix_trans',false);
  coarse_trans  = coarse_trans(intInx);                             % Restricted to interior

  %Prolongate variables - averages 
  
  coarse_p_init = accumarray(partition, p_ad.val)./weighting;
  coarse_sW_init = accumarray(partition,sW_ad.val)./weighting;
  % Until a better aproach is found, the ADI varaables is re-initiated 
  [coarse_p_ad, coarse_sW_ad] = initVariablesADI(coarse_p_init, coarse_sW_init);

  %Prolongate residual/defect - sum
  nc = G.cells.num;
  pIx = 1:nc;
  sIx = (nc+1):(2*nc);
  [coarse_defect_p_ad, coarse_defect_s_ad]= initVariablesADI(accumarray(partition,defect(pIx)),accumarray(partition,defect(sIx))); 

  
  %% Coarsening all helpfunctions and variables
   
  pv_r_coarse = poreVolume(CG, coarse_rock);
  pv_coarse   = @(p) pv_r_coarse .* exp( cr * (p - p_r) );
  nc_coarse = CG.cells.num;
  pIx_coarse = 1:nc_coarse;
  sIx_coarse = (nc_coarse+1):(2*nc_coarse);
  grad_coarse = @(x)C_coarse*x; % Discrete gradient
  div_coarse  = @(x)-C_coarse'*x; % Discrete divergence
  avg_coarse  = @(x) 0.5 * (x(N(:,1)) + x(N(:,2))); % Averaging
  upw_coarse = @(flag, x) flag.*x(N(:, 1)) + ~flag.*x(N(:, 2)); % Upwinding
  prodIndex_coarse = CG.cells.num;
  gradz_coarse  = grad_coarse(CG.cells.centroids(:,3));
end