%% Multi-echo fat-suppressed MR thermometry using iterative separation 
%% of baseline water and fat images
%
%  Sample function to run on simulated data
%
% Creators: Megan Poorman, William Grissom
% Location: Vanderbilt University Institute of Imaging Science
% Created: 08/2014
% Updated: 05/2017
% Do not reproduce, distribute, or modify without proper citation according
% to license file
%

clear all;close all;

%% %%%%%%%%%%%%%%%%%%%%%%%%
% algorithm parameters
%%%%%%%%%%%%%%%%%%%%%%%%%
libsource = 'true'; % Switch to use true Water/Fat baseline images or WF separation output ('toolbox')
% savedWFfile = 'WF_percentages_test/simulationDouble_mulit_WF';
nciter = 5; % # c iterations
nmiter = 5; % # m iterations
niter = 10; % # outer iterations for script version
lam = 10^-5; % l1 regularization parameter for script version
mthresh = 0.001; % heat-induced frequency threshold for 'significant' heat
stopthresh = 0.0001; % algorithm stopping threshold
initWithPreviousEstimate = 1;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scanner+sequence parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b0 = 3; % Tesla
dim = 128; % image size
fatppm = 3.4; % mean fat offset in ppm, used to calculate echo times
fatfrq = b0*42.5778*fatppm; % fat frequency in Hz for echo time calc
tecent = 0.01; % s, middle te time
tes = tecent + [-1/fatfrq -1/fatfrq/2 0 1/fatfrq/2 1/fatfrq]; % The TEs in sec
prc = 1; % direction of precession (scanner-dependent) +1 = heating induces negative apparent freq shift; -1 = heating induces positive apparent freq shift
order = 0; % polynomial background phase shift order 
hssig2 = 0.1; % hot spot sigma^2
dw0shift = 2*pi*10;%2*pi; % background polynomial frequency shift (will be replicated to all polynomial coeffs) (rad/sec)
phi0shift = 0;% Tx/Rx gain 
alpha = -0.01; % ppm/deg C;
dynmotind = 2; % index of motion vector to use for dynamic image
gamma = 42.5778; %MHz/T
noiselevel = 0; % add in noise (0.2 for SNR of 40)
noise = noiselevel*randn(dim,dim,length(tes));
fatPercentage = 0.5; %0.5 = 50% fat

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate baseline images 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--- generate the hot spot
disp('Generating hotspot');
[x,y] = meshgrid(-dim/2:dim/2-1);
hotspot = exp(-hssig2*(x.^2 + y.^2)); % normalized hot spot

maxtemps = 0:4:40; % KEEP 0 as first dyn to get lib with 0

for ii = 1:length(maxtemps)
    hotspotTrue(:,:,ii) = hotspot*alpha*maxtemps(ii)*b0*2*pi*gamma;
end


%--- generate the true W,F images
disp('Generating true Water and Fat images');
[x,y] = meshgrid(-dim/2:dim/2-1);
WTrue = x.^2 + 1.5*y.^2 <= (0.8*dim/2)^2; % water oval

f = (x.^2 + y.^2) <= (0.2*dim/2)^2; % small fat circle
FTrue = circshift(f,[0 -dim/4]) + f + circshift(f,[0 dim/4]); % replicate fat circle
R2starTrue = zeros(dim);
dw0True = zeros(dim)+dw0shift;

%apply WF ratio
WTrue = WTrue*(1-fatPercentage);
FTrue = FTrue*fatPercentage;

%---build polynomial matrix A
disp('Forming A polynomial Matrix');

[yc,xc] = meshgrid(linspace(-1/2,1/2,dim));
yc = yc(:);
xc = xc(:);
A = [];
for yp = 0:order
    for xp = 0:(order-yp)
        A = [A (xc.^xp).*(yc.^yp)];
    end
end

%--- get true polynomial coefficients and spatial phase shift
cTrue = zeros(size(A,2),1)+ dw0shift;  
AcTrue = reshape(A*cTrue,[dim dim]);% polynomial background phase shift

%--- create phi static shift
phiTrue = phi0shift; 

%--- generate library
disp('Generating image library');
motionvec = [-dim/8 0 dim/8]; % vector of object positions (integer) 
imgslib = calcmeimgs(WTrue,FTrue,tes,b0,dw0True,R2starTrue,0,prc)+noise;
imgslib = repmat(imgslib,[1 1 1 length(motionvec)]);

for ii = 1:length(motionvec) 
    
    % shift the images in first spatial dimension to simulate motion
    imgslib(:,:,:,ii) = circshift(imgslib(:,:,:,ii),motionvec(ii));
    
    % if we will not use Graphcut (usetruewf == 1), we will need the following:
    WTruebase(:,:,ii) = circshift(WTrue,motionvec(ii));
    FTruebase(:,:,ii) = circshift(FTrue,motionvec(ii));
    R2starTruebase(:,:,ii) = circshift(R2starTrue,motionvec(ii));
    dw0Truebase(:,:,ii) = circshift(dw0True,motionvec(ii));
    
    % to test the functions using the true values for all other variables,
    % we need the following:

    hotspotTrueshift(:,:,:,:,ii) = circshift(hotspotTrue,motionvec(ii));

        
    AcTrueshift(:,:,ii) = circshift(AcTrue,motionvec(ii));
    phiTrueshift(:,:,ii) = circshift(phiTrue,motionvec(ii));%MEP 20150727
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%
% generate dynamic images
%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Generating Dynamic images');

for ii = 1:length(maxtemps)
    imgsdyn(:,:,:,ii) = calcmeimgs(WTruebase(:,:,dynmotind).*exp(1i*angle(phiTrueshift(:,:,dynmotind))),FTruebase(:,:,dynmotind).*exp(1i*angle(phiTrueshift(:,:,dynmotind))),...
        tes,b0,dw0Truebase(:,:,dynmotind) + AcTrueshift(:,:,dynmotind),...%MEP 20150727
        R2starTruebase(:,:,dynmotind),hotspotTrueshift(:,:,ii,:,dynmotind),prc)+noise;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply water/fat separation to baseline images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Separating baseline WF');
switch libsource
    
    case 'toolbox'
        disp('using toolbox');
        %%%INSERT: code here to run data on ISMRM water/fat toolbox
        % use format as follows accounting for library of motion
        % ie perform separation on each position independently 
        warning('no toolbox outputs set');
        return
        % Wlib = outParams.species(1).amps;
        % Flib = outParams.species(2).amps;
        % dw0lib = outParams.fieldmap*2*pi;
        % R2starlib = outParams.r2starmap;
        
    case 'true'
        disp('using true water/fat');
        Wlib = WTruebase;
        Flib = FTruebase;
        dw0lib = dw0Truebase;
        R2starlib = R2starTruebase;
end



%% %%%%%%%%%%%%%%%%%%%% TEST COMBINED FUNCTION
disp('Testing combined function')

rad2degC = 1/(2*pi*b0*alpha*42.58);

% set up parameter structures

% algorithm parameters
algp.order = order;
algp.niter = niter;
algp.nmiter = nmiter;
algp.nciter = nciter;
algp.lam = lam;
algp.mthresh = mthresh;
algp.stopthresh = stopthresh;
algp.modeltest = 0;

% scan parameters
scanp.dim = [dim dim];
scanp.b0 = b0;
scanp.tes = tes;
scanp.prc = prc;

% baseline library
lib.Wlib = Wlib;
lib.Flib = Flib;
lib.dw0lib = dw0lib;
lib.R2starlib = R2starlib;

hsmask = (hotspotTrueshift(:,:,end,:,dynmotind)*rad2degC) > 0.5;
mfunc = [];

for ii = 1:length(maxtemps)
    if initWithPreviousEstimate && ii == 1
        [mfunc(:,:,ii),wts(:,ii),cfunc(:,ii),A,phi(:,:,ii)] = thermo_hybrid_waterfat(imgsdyn(:,:,:,ii), algp, scanp, lib,1,0,zeros(dim,dim));
    elseif initWithPreviousEstimate && ii > 1
        [mfunc(:,:,ii),wts(:,ii),cfunc(:,ii),A,phi(:,:,ii)] = thermo_hybrid_waterfat(imgsdyn(:,:,:,ii), algp, scanp, lib,1,0,mfunc(:,:,ii-1));
    else
        [mfunc(:,:,ii),wts(:,ii),cfunc(:,ii),A,phi(:,:,ii)] = thermo_hybrid_waterfat(imgsdyn(:,:,:,ii), algp, scanp, lib,1,0); %phi switch off, m switch on (solve for temperature without Tx/Rx gain)
    end
    RMSE_func(:,ii) = rmse(mfunc(:,:,ii),hotspotTrueshift(:,:,ii,:,dynmotind),hsmask);
end

figure;semilogy(maxtemps,RMSE_func); xlabel('Hotspot \Delta\circ C'); ylabel('RMSE \circ C');title('Temperature Error with Algorithm');ylim([0.01 10]);
figure;subplot(131);imagesc(mfunc(:,:,end));axis image;colorbar;title('Solved hotspot');
subplot(132);imagesc(hotspotTrueshift(:,:,end,:,dynmotind));axis image;colorbar;title('True hotspot');
subplot(133);imagesc(mfunc(:,:,end)-hotspotTrueshift(:,:,end,:,dynmotind));axis image;colorbar;title('Difference');

