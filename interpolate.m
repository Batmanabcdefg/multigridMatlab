function [fine_p, fine_sW] = interpolate(CG, coarse_p_ad, coarse_sW_ad)
%%Apply the coarse grid partition field to interpolate the variables from the
%%coarse grid into the fine grid

% Need of: invertPartition(p) ? from invertPartition.m -> utils ->
% coarsegrid ->modules

  %fine_p_ad = coarseDataToFine(CG, coarse_p_ad);
  
  %fine_sW_ad = coarseDataToFine(CG, coarse_sW_ad);

  fine_p = coarse_p_ad.val(CG.partition);
  
  fine_sW = coarse_sW_ad.val(CG.partition);
end