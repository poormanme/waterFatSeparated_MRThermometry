function [w,libinds] = f_update(imgs,Z,usequadprog,nlib) 
%% Real-time-compatible multiple echo fat-suppressed MR thermometry using
%% iterative separation of baseline water and fat images
% Function to find the weights that best match the baseline library to the
% multi echo images
%
% Creators: Megan Poorman, William Grissom
% Location: Vanderbilt University Institute of Imaging Science
% Created: 08/2014
% Updated: 05/2017
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%
% Inputs:
%   imgs:        Dynamic multi-echo images (Nx x Ny x Necho)
%   Z:           Model image for each echo and each baseline (Nx x Ny x Necho x Nbaseline)
%   usequadprog: switch to use quadprog for constrained solution, or just
%   unconstrained least-squares
%   nlib:        number of library images
%
% Output: 
%   w:           baseline weights


if ~exist('usequadprog','var')
    usequadprog = 1;
end
if ~exist('nlib','var')
    nlib = size(Z,4); % if user doesn't specify nlib, don't compress
end

if size(Z,4) > 1 %---- if more than one library image
        
    % reshape Z into a matrix by collapsing space and echo dims
    Z = permute(Z,[4 1 2 3]); Z = Z(:,:).';
    
    % set up constraint: coefficients must add to 1
    Ceq = ones(1,size(Z,2));
    beq = 1;

    % Calculate inner product
    fprime = -real(imgs(:)'*Z);
    
    if nlib < size(Z,2)
        % find indices of largest nlib inner products (i.e. the nlib
        % library images that best correlate with the dynamic image).
        % these indices will be used for the fit.
        [~,libinds] = sort(abs(fprime),'descend');
        libinds = libinds(1:nlib);
        Z = Z(:,libinds);
        fprime = fprime(libinds);
        Ceq = Ceq(1:nlib);
    end

    % Calculate Hessian
    H = real(Z'*Z)';
    
    if usequadprog
    
        % solve for weights
        options = optimset;options.MaxIter = 100000;
        options.Display = 'off';
        warning('off','all')
        w = quadprog(H,fprime,[],[],Ceq,beq,zeros(size(Z,2),1),ones(size(Z,2),1),[],options);
        warning('on','all')
    
    else
        
        w = -H\fprime';
        
    end

else  %---- only one baseline, so weight vector = [1];

    w = 1;

end
