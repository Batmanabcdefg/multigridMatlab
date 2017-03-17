function [p_approx, sW_approx,nit] ...
    = multigridCycleV2(v1,v2,model,p_ad_0,sW_ad_0,tol,maxits,g,t,dt,pIx,sIx);
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
  [p_ad,sW_ad,defect] = newtonTwoPhaseADV2(model,p_ad_0,sW_ad_0,tol,v1,g,t,dt,pIx,sIx);
 
 
  %% Set up of coarse grid
  [coarse_model,coarse_p_ad, coarse_sW_ad,coarse_defect_p_ad,coarse_defect_sW_ad, ...
    pIx_coarse, sIx_coarse,coarse_p_ad_0, coarse_sW_ad_0] ...
   = coarseningV2(model, p_ad, sW_ad, defect,pIx,sIx,p_ad_0,sW_ad_0);
    
  %% Multigrid core
  [correction_p,correction_sW,nit] ...
      = newtonTwoPhaseADV2(coarse_model,coarse_p_ad,coarse_sW_ad,tol,maxits,g,t,dt,pIx_coarse,sIx_coarse,coarse_p_ad_0, coarse_sW_ad_0);

%     figure
%     subplot(6, 2, 1); plot(model.G.cells.indexMap,p_ad.val);
%     title('Pressure')
%     subplot(6, 2, 3); plot(1:coarse_model.G.cells.num,coarse_p_ad.val);
%     title('Coarse Pressure')
%     subplot(6, 2, 5); plot(1:coarse_model.G.cells.num,correction_p);
%     title('Corrected coarse Pressure')
%     
%     subplot(6, 2, 2); plot(model.G.cells.indexMap,sW_ad.val);
%     title('Saturation')
%     subplot(6, 2, 4); plot(1:coarse_model.G.cells.num,coarse_sW_ad.val);
%     title('Coarse Saturation')
%     subplot(6, 2, 6); plot(1:coarse_model.G.cells.num,correction_sW);
%     title('Corrected coarse Saturation')
%     drawnow
  
  %% Interpolating soluton from coarsed grid and compute ccorrected approximation

[fine_correction_p, fine_correction_sW] = interpolate(coarse_model.G, correction_p, correction_sW);
  
 % p_ad.val = p_ad.val + fine_correction_p;
  
 % sW_ad.val = sW_ad.val + fine_correction_sW;
%   
% figure(5)
% subplot(2,1,1)
% plotCellData(model.G, p_ad.val,'EdgeColor','k','EdgeAlpha',.2)
% title('Fine-scale solution')
% subplot(2,1,2)
% plotCellData(coarse_model.G, correction_p,'EdgeColor','none')
% %plotFaces(CG,(1:model_coarse.G.faces.num)','FaceColor','none');
% title('Coarse-scale solution')

  % Postsmoothing
  [p_approx,sW_approx,nit] = newtonTwoPhaseADV2(model,p_ad,sW_ad,tol,v2,g,t,dt,pIx,sIx);
  
end