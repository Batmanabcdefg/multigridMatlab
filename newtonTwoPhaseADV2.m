function [p_ad, sW_ad,nit] =  ...
    newtonTwoPhaseADV2(model,p_ad,sW_ad,tol,maxits,g,dt,varargin)
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
   if(isempty(varargin))
       p0  = double(p_ad); % Previous step pressure
       sW0 = double(sW_ad);
   else
       p0 = double(varargin{1});
       sW0 = double(varargin{2});
   end
   nit = 0;
   resNorm = 1e99;
   old_res = resNorm;
   
  while (resNorm > tol) && (nit < maxits) && (old_res >= resNorm)
      old_res = resNorm;
     
      [water, oil] = computePhaseFlux(model,p_ad,sW_ad,dt,g,p0,sW0);
      if(isempty(varargin) || (length(varargin) ~= 3) || nit > 0)
        [water, oil] = computeBoundaryCondition(model,p_ad,sW_ad,water,oil);
      else
        boundaryCondition = varargin{3};
        water = boundaryCondition.water;
        oil = boundaryCondition.oil;
%         %[water, oil] = boundaryConditions(model,p_ad,sW_ad,water,oil);
%         water.val = water.val - boundaryCondition.water.val;
%         oil.val = oil.val - boundaryCondition.oil.val;
       [water, oil] = initVariablesADI(water.val,oil.val);
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
      fprintf('  Iteration %3d:  Res = %.4e\n', nit, resNorm);
   end
end