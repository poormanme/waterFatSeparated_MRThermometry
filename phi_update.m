function phi = phi_update(imgs,Z,phi,tes,prc)
%% Multi-echo fat-suppressed MR thermometry using iterative separation 
%% of baseline water and fat images
% Function to fit DC shift due to Tx/Rx gain differences
%
% Creators: Megan Poorman, William Grissom
% Location: Vanderbilt University Institute of Imaging Science
% Created: 08/2014
% Updated: 05/2017
% Do not reproduce, distribute, or modify without proper citation according
% to license file
% 
%
% Inputs:
%   imgs: Dynamic multi-echo images (Nx x Ny x Necho)
%   Z:    Model image for each echo and each baseline (Nx x Ny x Necho x Nbaseline)
%   phi:  initial guess
%   tes:  echo times
%   prc:  direction of precession
%
% Output: 
%   phi: dim x dim fit

phi = phi(:);
innprod = Z.*conj(imgs);
innprod = permute(innprod,[3 1 2]); innprod = innprod(:,:).';

% loop over echoes to calc derivatives and curvatures
t1 = abs(innprod);
t2 = angle(innprod);
t3 = sin(t2);
t4 = t1.*t3;
hesssum = 0;gradsum = 0;
for ii = 1:length(tes)    
    hesssum = hesssum + t4(:,ii)./t2(:,ii);
    gradsum = gradsum - t4(:,ii);  
end

% update phi
hesssum(isnan(hesssum)) = 0;

dphi = (gradsum)./hesssum; 
dphi(isnan(dphi) | isinf(dphi)) = 0;

phi = phi + dphi;