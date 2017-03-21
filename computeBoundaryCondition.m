function [water, oil] = computeBoundaryCondition(model,p_ad,sW_ad,water,oil)
      % Insert volumetric source term multiplied by density
      water(model.well.injIndex) = water(model.well.injIndex) - model.well.inRate.*model.water.rhoWS;
      % Set production cells to fixed pressure of 200 bar and zero water
      water(model.well.prodIndex) = p_ad(model.well.prodIndex) - 200*barsa;
      oil(model.well.prodIndex) = sW_ad(model.well.prodIndex);
      
end