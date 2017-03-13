function [p_ad, sW_ad,res,nit] =  ...
    newtonTwoPhaseADV2(model,p_ad,sW_ad,tol,maxits,g,dt,pIx,sIx)
   resNorm = 1e99;
   %% This function is simply copyed out from the original file twoPhaseAD.m
   %It is not indended to be a complete function on it's own due to the
   %massive need of input
   
   %%
   p0  = double(p_ad); % Previous step pressure
   sW0 = double(sW_ad);
   nit = 0;
 
  while (resNorm > tol) && (nit < maxits)
      % Evaluate properties
      rW = model.water.rhoW(p_ad);
      rW0 = model.water.rhoW(p0);
      rO = model.oil.rhoO(p_ad);
      rO0 = model.oil.rhoO(p0);
      
      % Define pressure drop over interface for both phases
      dp = model.operator.grad(p_ad);
      dpW = dp - g*model.operator.avg(rW).*model.operator.gradz; % Water
      dpO = dp - g*model.operator.avg(rO).*model.operator.gradz; % Oil
      % Pore volume of cells at current pressure and previous pressure 
      vol0 = model.rock.pv(p0);
      vol = model.rock.pv(p_ad);
      % Mobility: Relative permeability over constant viscosity
      mobW = model.water.krW(sW_ad)./model.water.muW;
      mobO = model.oil.krO(1-sW_ad)./model.oil.muO;
      % Define phases fluxes. Density and mobility is taken upwinded (value
      % on interface is defined as taken from the cell where the phase flow
      % is currently coming from). This gives more physical solutions than
      % averaging or downwinding.
      vW = -model.operator.upw(double(dpW) <= 0, rW.*mobW).*model.T.*dpW;
      vO = -model.operator.upw(double(dpO) <= 0, rO.*mobO).*model.T.*dpO;
      % Conservation of water and oil
      water = (1/dt).*(vol.*rW.*sW_ad - vol0.*rW0.*sW0) + model.operator.div(vW);
      oil   = (1/dt).*(vol.*rO.*(1-sW_ad) - vol0.*rO0.*(1-sW0)) + model.operator.div(vO);
      % Insert volumetric source term multiplied by density
      water(model.well.injIndex) = water(model.well.injIndex) - model.well.inRate.*model.water.rhoWS;
      % Set production cells to fixed pressure of 200 bar and zero water
      water(model.well.prodIndex) = p_ad(model.well.prodIndex) - 200*barsa;
      oil(model.well.prodIndex) = sW_ad(model.well.prodIndex);
      % Collect all equations
      eqs = {oil, water};
      % Concatenate equations and solve for update:
      eq  = cat(eqs{:});
      J   = eq.jac{1};  % Jacobian
      res = eq.val;     % residual
      upd = -(J \ res); % Newton update
      % Update variables
      p_ad.val   = p_ad.val   + upd(pIx);
      sW_ad.val = sW_ad.val + upd(sIx);
      sW_ad.val = min(sW_ad.val, 1);
      sW_ad.val = max(sW_ad.val, 0);
      
      resNorm = norm(res);
      nit     = nit + 1;
      fprintf('  Iteration %3d:  Res = %.4e\n', nit, resNorm);
   end
end