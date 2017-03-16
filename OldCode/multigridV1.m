function [p_approx, sW_approx,nit] = multigridV1(v1_its,v2_its,G, rock,p_ad,sW_ad,tol,maxits,rhoW,rhoO,grad,gradz,pv,krW,krO,muW,muO,T,g,avg,upw,dt,div,injIndex,inRate,rhoWS,prodIndex,pIx,sIx,cr,p_r)
  %% Presmoothing
  [p_ad,sW_ad,defect] = newtonAD(p_ad,sW_ad,tol,v1_its,rhoW,rhoO,grad,gradz,pv,krW,krO,muW,muO,T,g,avg,upw,dt,div,injIndex,inRate,rhoWS,prodIndex,pIx,sIx);

  %% Set up of coarse grid
  [CG,coarse_rock,coarse_trans, coarse_p_ad, coarse_s_ad,coarse_defect_p_ad,coarse_defect_s_ad,N, ...
       pv_coarse, pIx_coarse, sIx_coarse, grad_coarse, div_coarse, avg_coarse, upw_coarse, prodIndex_coarse, gradz_coarse] ...
   = coarsening(G, rock, T, p_ad, sW_ad,defect,cr,p_r);

  %% Multigrid core
  [correction_p,correction_sW,nit] = newtonAD(coarse_defect_p_ad,coarse_defect_s_ad,tol,maxits,rhoW,rhoO, grad_coarse,gradz_coarse, ...
      pv_coarse,krW,krO,muW,muO,coarse_trans,g,avg_coarse,upw_coarse,dt,div_coarse,injIndex,inRate,rhoWS,prodIndex_coarse,pIx_coarse,sIx_coarse);

  
  %% Interpolating soluton from coarsed grid and compute ccorrected approximation
  [fine_correction_p, fine_correction_sW] = interpolate(CG, correction_p, correction_sW);
  
  p_ad.val = p_ad.val + fine_correction_p;
  
  sW_ad.val = sW_ad.val + fine_correction_sW;
  
  %%Postsmoothing
  [p_approx,sW_approx] = newtonAD(p_ad,sW_ad,tol,v2_its,rhoW,rhoO,grad,gradz,pv,krW,krO,muW,muO,T,g,avg,upw,dt,div,injIndex,inRate,rhoWS,prodIndex,pIx,sIx);
  
end