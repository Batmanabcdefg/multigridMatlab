function [p_ad, sW_ad,nit] =  ...
    newtonTwoPhaseAD(model,p_ad,sW_ad,tol,maxits,dt,p_ad_0,sW_ad_0,varargin)
   %% Function description
   %
   % PARAMETERS:
   % model    - System model structure with grid, rock, phases and operator
   %            substructs
   %
   % p_ad     - ADI struct for the pressure
   % s_ad     - ADI struct for the saturation
   % tol      - error tolerance for newton
   % maxits   - The maximum number of newton iterations perfomed
   % g        - gravity constant 
   % dt       - current time step
   % p_ad_0   - Pressure ADI from previous timestep
   % sW_ad_0  - Saturation ADI from previous timestep
   % varargin - boundary condition from multigrid cycle
   %
   % RETURNS:
   % p_ad     - Approximated values of the pressure stored in ADI structure
   % s_ad     - Approximated values of the saturation stored in ADI structure
   % res      - Residual from the final newton iteration
   % nit      - Number of iterations performed
   % 
   % COMMENTS:
   % - The body of this function is simply copyed out and slightly modified from the original file twoPhaseAD.m
   %
   % SEE ALSO:
   %

%% Body

   p0 = double(p_ad_0);
   sW0 = double(sW_ad_0);
   
   nit = 0;
   resNorm = 1e99;
   old_res = resNorm;
   
  while (resNorm > tol) && (nit < maxits) && (old_res >= resNorm)
      old_res = resNorm;
     
      [water, oil] = computePhaseFlux(model,p_ad,sW_ad,dt,p0,sW0);
        
      % Check wether a defect have been passed or not.
      if(isempty(varargin) || isempty(varargin{1}))% || nit > 0)
          [water, oil] = computeBoundaryCondition(model,p_ad,sW_ad,water,oil);
      
      else
          boundaryCondition = varargin{1};
          water = water + boundaryCondition.water;
          oil = oil + boundaryCondition.oil;
          water_val = water(model.well.prodIndex).val;
          oil_val = oil(model.well.prodIndex).val;
          
          water(model.well.prodIndex) = water(model.well.prodIndex) ... %- water(model.well.prodIndex) ...
              + p_ad(model.well.prodIndex) - p_ad(model.well.prodIndex).val - water_val;
      
          oil(model.well.prodIndex) = oil(model.well.prodIndex) ...% - oil(model.well.prodIndex) ...
              + sW_ad(model.well.prodIndex)- sW_ad(model.well.prodIndex).val - oil_val;    
      end
     
      % Collect all equations
      eqs = {oil, water};
      % Concatenate equations and solve for update:
      eq  = cat(eqs{:});
      J   = eq.jac{1};  % Jacobian
      res = eq.val;     % residual
      upd = -(J \ res); % Newton update
      % Update variables
      p_ad.val  = p_ad.val  + upd(model.pIx);
      sW_ad.val = sW_ad.val + upd(model.sIx);
      sW_ad.val = min(sW_ad.val, 1);
      sW_ad.val = max(sW_ad.val, 0);
%     figure
%     plot(1:model.G.cells.num,sW_ad.val);
%     title('Saturation')
      
      resNorm = norm(res);
      nit     = nit + 1;
%       fprintf('  Iteration %3d:  Res = %.4e\n', nit, resNorm);
  end
   fprintf('  Iteration %3d:  Res = %.4e\n', nit, resNorm)
%    fprintf('Iterantions: %3d\n', nit)
end