function [p_approx, sW_approx,nit] ...
    = multigridV2(v1,v2,model,p_ad,sW_ad,tol,maxits,g,dt,pIx,sIx);
  %% Function description
  %
  % PARAMETERS:
  % model    - System model structure with grid, rock, phases and operator
  %            substructs
  %
  % v1_iter  - Number of presmoothing sweeps
  % v2_iter  - Number of postsmoothing sweeps
  % p_ad     - ADI struct for the pressure
  % s_ad     - ADI struct for the saturation
  % tol      - error tolerance for newton at the most coarse grid
  % maxits   - The maximum number of newton iterations perfomed on the most
  %            coarsed grid
  % g        - gravity constant 
  % dt       - current time step
  % pIx      - Index array for pressure values
  % sIx      - Index array for saturation values
  %
  % RETURNS:
  % p_ad     - Approximated values of the pressure stored in ADI structure
  % s_ad     - Approximated values of the saturation stored in ADI structure
  % nit      - Number of newton iterations performed at the most coarse grid
  % COMMENTS:
  %   This is currently a linear two-grid cycle
  %
  % SEE ALSO:
  %

 %% Presmoothing
 [p_ad,sW_ad,defect] = newtonTwoPhaseADV2(model,p_ad,sW_ad,tol,v1,g,dt,pIx,sIx);
 
  %% Set up of coarse grid
  [model_coarse,coarse_p_ad, coarse_sW_ad,coarse_defect_p_ad,coarse_defect_sW_ad, ...
    pIx_coarse, sIx_coarse] ...
   = coarseningV2(model, p_ad, sW_ad, defect,pIx,sIx);

  %% Multigrid core
  [correction_p,correction_sW,nit] ...
      = newtonTwoPhaseADV2(model_coarse,coarse_defect_p_ad,coarse_defect_sW_ad,tol,maxits,g,dt,pIx_coarse,sIx_coarse);

  
  %% Interpolating soluton from coarsed grid and compute ccorrected approximation
  [fine_correction_p, fine_correction_sW] = interpolate(model_coarse.G, correction_p, correction_sW);
  
  p_ad.val = p_ad.val + fine_correction_p;
  
  sW_ad.val = sW_ad.val + fine_correction_sW;
  
  % Postsmoothing
  [p_approx,sW_approx,nit] = newtonTwoPhaseADV2(model,p_ad,sW_ad,tol,v2,g,dt,pIx,sIx);
  
end