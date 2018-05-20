function c = c_update(imgs,Z,A,c,tes,prc,m)%Z2)
%% Real-time-compatible multiple echo fat-suppressed MR thermometry using
%% iterative separation of baseline water and fat images
% Function to find the coefficients on the backgrond polynomial fit to
% background phase
%
% keyboard
% Creators: Megan Poorman, William Grissom
% Location: Vanderbilt University Institute of Imaging Science
% Created: 08/2014
% Updated: 05/2017
% Do not reproduce, distribute, or modify without proper citation according
% to license file 
% NEWTON's METHOD GRADIENT DESCENT
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

% %%%ORIG
% reshape images and model into vectors
% keyboard
innprod = Z.*conj(imgs);
innprod = permute(innprod,[3 1 2]); 
innprod = innprod(:,:).';

% loop over echoes to calc derivatives and curvatures
% t1 = abs(Z2).*abs(imgs);
% t1 = permute(t1,[3 1 2]);
% t1 = t1(:,:).';%
t1 = (abs(innprod));
t2 = (angle(innprod)); %wrapped w/ heat

t3 = sin(t2);
t4 = t1.*t3;


% keyboard

% keyboard
hesssum = 0;gradsum = 0;
for ii = 1:length(tes)    
% aa(:,ii) = t4(:,ii)./t2(:,ii)*(tes(ii))^2;
% bb(:,ii) = - (1-2*prc)* t4(:,ii)*(tes(ii));
    hesssum = hesssum + t4(:,ii)./t2(:,ii)*(tes(ii))^2;
    gradsum = gradsum - (1-2*prc)* t4(:,ii)*(tes(ii));%

end
% disp(['Hess sum: ' num2str(hesssum(8150,1)) '------- Grad sum: ' num2str(gradsum(8150,1)) '--------c: ' num2str(c)]);
% aa = reshape(aa,[128 128 5]);
% bb = reshape(bb,[128 128 5]);
% figure;subplot(121);plot(squeeze(aa(64,64,:)));title('Hessian');xlabel('TE #');ylabel('Value');
% subplot(122);plot(squeeze(aa(64,64,:)));title('Gradient');xlabel('TE #');ylabel('Value');

%     figure(100);im(reshape(gradsum,[128 128]));drawnow;pause(0.1);
% % %POOJA
% dims = size(imgs);
% imgs = reshape(imgs,[prod(dims(1:2)) dims(3)]);
% Z = reshape(Z,[prod(dims(1:2)) dims(3)]);
% t1 = abs(imgs).*(Z);
% 
% for ii = 1:length(tes)
% %     keyboard
%     t2 = A*c*tes(ii)+m(:)*tes(ii)+angle(Z(:,ii))+angle(imgs(:,ii));
%     t3 = sin(t2);
% %     keyboard
%     t4 = mod(t2+pi,2*pi)-pi;
% %     keyboard
%     t5 = t1(:,ii).*t3;
%     gradsum = gradsum + (1-2*prc)* t5*tes(ii);
%     hesssum = hesssum + t5./t4*tes(ii)^2;
% 
% 
% end
% update c
hesssum(isnan(hesssum)) = 0;
% gradsum(isnan(gradsum)) = 0;
% keyboard
dc = -(A'*bsxfun(@times,hesssum,A))\(A'*gradsum);
% keyboard
dc(isnan(dc)) = 0;
c = c + dc;
% disp(['c is ' num2str(c) ' -grad sum is ' num2str(-A'*gradsum)]);pause(0.3);