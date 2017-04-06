function [coarse_model,p_ad_coarse, sW_ad_coarse, p_ad_0_coarse, sW_ad_0_coarse, coarse_defect] ...
    = coarseningV2(model, p_ad, sW_ad, p_ad_0, sW_ad_0, varargin)
  %% Function description
  %
  % PARAMETERS:
  % model    - System model structure with grid, rock, phases and operator
  %            substructs
  % p_ad     - ADI struct for the pressure
  % s_ad     - ADI struct for the saturation
  % defect   - The defect of the current approximization. (NOTE: varargin)
  %
  % RETURNS:
  % p_ad_coarse - Coarsed version of the pressure eqs. with reinitialized
  %               ADI structure
  % sW_ad_coarse - Coarsed version of the pressure eqs. with reinitialized
  %               ADI structure
  % defect_p_ad_coarse - Coarsed version of the pressure defect
  % defect_s_ad_coarse - Coarsed version of the saturation defect
  %
  % COMMENTS:
  % - This coarsening function is currently not optimized for performance
  % - The function may be bugged in current state
  %
  % SEE ALSO:
  %

  %% Create new coarse model
  % Set up coarse grid
  coarse_dims = ceil(model.G.cartDims/2);
  
  partition  = partitionCartGrid(model.G.cartDims,coarse_dims);
  
  if(size(partition,1) < model.G.cells.num)
      partition = cat(1,partition,partition(end));
  end
  partition(model.well.prodIndex) = partition(model.well.prodIndex) + 1;
%   partition(model.G.cells.num) = partition(model.G.cells.num) + 1;
  
  CG = generateCoarseGrid(model.G, partition);
  CG = coarsenGeometry(CG);
  
  % Define coarse rock
  rock_coarse = makeRock(CG, 30*milli*darcy, 0.3);
  
  weighting = accumarray(partition,1);
  
  coarse_model = initiateModel(CG, rock_coarse,model,weighting);

  % Add fields to the coarse grid to ensure that it passes as a
  % regular grid for our purposes.
  coarse_model.G.cartDims = coarse_dims;
  
  %% Restrict AD variables and defect
  
  
  coarse_p_init = accumarray(partition, p_ad.val)./weighting;
  coarse_sW_init = accumarray(partition,sW_ad.val)./weighting;
  
  % Until a better aproach is found, the ADI varaables is re-initiated 
  [p_ad_coarse, sW_ad_coarse] = initVariablesADI(coarse_p_init, coarse_sW_init);

  coarse_p_0 = accumarray(partition, p_ad_0.val)./weighting;
  coarse_sW_0 = accumarray(partition,sW_ad_0.val)./weighting;
  
  % Until a better aproach is found, the ADI varaables is re-initiated 
  [p_ad_0_coarse, sW_ad_0_coarse] = initVariablesADI(coarse_p_0, coarse_sW_0);
  
  coarse_defect = 0;
  if(~isempty(varargin))
    defect = varargin{1};
    [coarse_water_defect,coarse_oil_defect] ...
      = initVariablesADI(accumarray(partition,defect.water.val),accumarray(partition,defect.oil.val));

    coarse_defect = struct('water',coarse_water_defect, 'oil',coarse_oil_defect);
  end
  
  % Well conditions
  
  coarse_model.well.inRate = model.well.inRate;
  coarse_model.well.outRate = model.well.outRate;

end