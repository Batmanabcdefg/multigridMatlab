function [p_ad, sW_ad,res,nit] =  ...
    newtonAD(p_ad,sW_ad,tol,maxits,rhoW,rhoO,grad,gradz,pv,krW,krO,muW,muO,T,g,avg,upw,dt,div,injIndex,inRate,rhoWS,prodIndex,pIx,sIx)
   resNorm = 1e99;
   %% This function is simply copyed out from the original file twoPhaseAD.m
   %It is not indended to be a complete function on it's own due to the
   %massive need of input
   
   %%
   p0  = double(p_ad); % Previous step pressure
   sW0 = double(sW_ad);
   nit = 1;
 
  while (resNorm > tol) && (nit <= maxits)
      % Evaluate properties
      rW = rhoW(p_ad);
      rW0 = rhoW(p0);
      rO = rhoO(p_ad);
      rO0 = rhoO(p0);
      
      % Define pressure drop over interface for both phases
      dp = grad(p_ad);
      dpW = dp - g*avg(rW).*gradz; % Water
      dpO = dp - g*avg(rO).*gradz; % Oil
      % Pore volume of cells at current pressure and previous pressure 
      vol0 = pv(p0);
      vol = pv(p_ad);
      % Mobility: Relative permeability over constant viscosity
      mobW = krW(sW_ad)./muW;
      mobO = krO(1-sW_ad)./muO;
      % Define phases fluxes. Density and mobility is taken upwinded (value
      % on interface is defined as taken from the cell where the phase flow
      % is currently coming from). This gives more physical solutions than
      % averaging or downwinding.
      vW = -upw(double(dpW) <= 0, rW.*mobW).*T.*dpW;
      vO = -upw(double(dpO) <= 0, rO.*mobO).*T.*dpO;
      % Conservation of water and oil
      water = (1/dt).*(vol.*rW.*sW_ad - vol0.*rW0.*sW0) + div(vW);
      oil   = (1/dt).*(vol.*rO.*(1-sW_ad) - vol0.*rO0.*(1-sW0)) + div(vO);
      % Insert volumetric source term multiplied by density
      water(injIndex) = water(injIndex) - inRate.*rhoWS;
      % Set production cells to fixed pressure of 200 bar and zero water
      water(prodIndex) = p_ad(prodIndex) - 200*barsa;
      oil(prodIndex) = sW_ad(prodIndex);
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