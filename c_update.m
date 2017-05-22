function c = c_update(imgs,Z,A,c,tes,prc)
%% Real-time-compatible multiple echo fat-suppressed MR thermometry using
%% iterative separation of baseline water and fat images
% Function to find the coefficients on the backgrond polynomial fit to
% background phase
%
% Creators: Megan Poorman, William Grissom
% Location: Vanderbilt University Institute of Imaging Science
% Created: 08/2014
% Updated: 05/2017
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   imgs: Dynamic multi-echo images (Nx x Ny x Necho)
%   Z:    Model image for each echo and each baseline (Nx x Ny x Necho x Nbaseline)
%   A:    polynomial basis functions
%   c:    previous coefficient estimate
%   tes:  echo times
%   prc:  direction of precession
%
% Output: 
%   c:    polynomial coefficents


% reshape images and model into vectors
innprod = Z.*conj(imgs);
innprod = permute(innprod,[3 1 2]); innprod = innprod(:,:).';

% loop over echoes to calc derivatives and curvatures
t1 = abs(innprod);
t2 = angle(innprod);
t3 = sin(t2);
t4 = t1.*t3;
hesssum = 0;gradsum = 0;
for ii = 1:length(tes)    
    hesssum = hesssum + t4(:,ii)./t2(:,ii)*tes(ii)^2;
    gradsum = gradsum - (1-2*prc)*t4(:,ii)*tes(ii);        
end

% update c
hesssum(isnan(hesssum)) = 0;
dc = -(A'*bsxfun(@times,hesssum,A))\(A'*gradsum); 
dc(isnan(dc)) = 0;
c = c + dc;