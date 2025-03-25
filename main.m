% 3D forward monte carlo simulation for calculation of rrs
% Simulates photon streams in the forward direction
% Keeps track of physical depth, z, zenith and azimuthal angles


%% User Input

% Seawater or Freshwater
water = 'sea'
% water = 'fresh'

% Open ocean or water tank
setup = 'ocean'
% setup = 'tank'

% Choose type of expriment
exp = 0 % Run MC0 - black sky, for reflectance only
% exp = 1 % Run MC1 - MS tracking, full sky, pre-determined downwelling photons
% exp = 2 % Run MC2 - black tank, limited FOV uniform sky
% exp = 3 % Run MC3 - black tank, pre-determined downwelling photons

% wavelengths
wave = (400:10:850)';
nwave = size(wave,1);

% Incident photon zenith angle in air (degrees)
solzen = 30; % If you change this, you need to change the downwelling irradiance!

% Nomber of photons
nph = 30*10^6;

% Cutoff depth for infinitely deep water
cdep = 10000;

% Create struct and set IOPs
nsets = 100;
dsto = struct();

for iset = 1:nsets
    dsto(iset).G = 2*rand;
    dsto(iset).X = rand;
    dsto(iset).S = ((0.02-0.01)*rand)+0.01;
    dsto(iset).Y = 1.2*rand;
end

% Exiting photon nominal angles in air (degrees)
% zennom is w.r.t. upward direction (0 degrees)
zennom = [92.5, 100, 110, 120, 130, 140, 150, 160, 170, 180]; % rows 1-10
azinom = [0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180]; % cols 1-13
n_zen = size(zennom,2);
n_azi = size(azinom,2);

% Generate solid angle in sr for each quad
squad = sin(zennom'*pi/180).*ones(10,13)*(10*pi/180)*(15*pi/180); % General Quad solid angles
squad(1,:) = sin(92.5'*pi/180)*(5*pi/180)*(15*pi/180); % Equator solid angles
squad(10,:) = 0; % Polar Cap
squad(10,1) = 2*pi*(1 - cos(5*pi/180)); % Polar cap 
% Divide 0 and 180 quads by 2
squad(:,1) = squad(:,1)/2; 
squad(:,13) = squad(:,13)/2;

%% Preliminary calculations
% Incident photon zen angle after passing through surface
zen = asin(sin(solzen*pi/180)/1.33);
azi = 0*pi/180;

% Load absorption and scattering coefficients
if water == 'sea'
	watervals = xlsread('water_vals.xlsx','A2:C277');
	aw = interp1(watervals(:,1),watervals(:,2),wave);
	bbw = interp1(watervals(:,1),watervals(:,3),wave);
	bw = 2*bbw;
elseif water == 'fresh'
	watervals = xlsread('water_vals_fresh.xlsx','A2:C192');
	aw = interp1(watervals(:,1),watervals(:,2),wave);
	bbw = interp1(watervals(:,1),watervals(:,3),wave);
	bw = 2*bbw;

% Compute IOPs and define zero arrays
for iset = 1:nsets
    dsto(iset).ag = dsto(iset).G*exp(-dsto(iset).S*(wave-440));
    dsto(iset).bp = (dsto(iset).X*(550./wave).^dsto(iset).Y)/0.018;
    dsto(iset).bt = dsto(iset).bp + bw;
    
    % Allocate for output
    dsto(iset).outcount = zeros(n_zen,n_azi,nwave);
    dsto(iset).outrrs = zeros(n_zen,n_azi,nwave);
    dsto(iset).outR = zeros(nwave,1);
    dsto(iset).outscount = zeros(n_zen,n_azi,nwave); % number of scatters per exiting photon   
    dsto(iset).tscount = zeros(nwave,1); % number of scatters per (all photons)
end


% Load direct downwelling radiance PDF and scale no of photons (annual average)
if setup == 'ocean'

	%% For open ocean
	drad_ang = xlsread('downrad.xlsx','solzen30_pdf','C4:D124');
	drad_pdf = xlsread('downrad.xlsx','solzen30_pdf','E4:AX124');
	drad_nph = round(nph*drad_pdf);
	pcount = sum(drad_nph,1);

	% Load Fournier-Fourand CDF values
	ffgam = xlsread('phasefunc.xlsx','FF18','B3:B1803');
	ffcdf = xlsread('phasefunc.xlsx','FF18','F3:F1803');
	
elseif setup == 'tank'

	%% For water tank
	Load direct downwelling radiance PDF and scale no of photons (specific date, BB9 waves)
	drad_ang = xlsread('downrad_26072018.xlsx','solzen30_pdf','C4:D124');

	drad_pdf30 = xlsread('downrad_26072018.xlsx','solzen30_pdf','E4:M124');
	drad_nph30 = round(nph*drad_pdf30);
	pcount30 = sum(drad_nph30,1);

	drad_pdf20 = xlsread('downrad_26072018.xlsx','solzen20_pdf','E4:M124');
	drad_nph20 = round(nph*drad_pdf20);
	pcount20 = sum(drad_nph20,1);

	% Load Calcite CDF values
	ccgam = xlsread('phasefunc.xlsx','Calcite','B2:B1802');
	cccdf = xlsread('phasefunc.xlsx','Calcite','G2:G1802'); 


%% Forward MC 
tic
% parfor iset = 1:nsets
parfor iset = 1:nsets
    tic
        
    for iwave = 1:nwave 
	
		if exp == 0
			% Run MC0 - black sky, for reflectance only
			 [dsto(iset).outcount(:,:,iwave),dsto(iset).outR(iwave),dsto(iset).outrrs(:,:,iwave)] = fwd_mc_deep0(nph,n_zen,n_azi,zen,cdep,aw(iwave),bw(iwave),dsto(iset).ag(iwave),dsto(iset).bp(iwave),ffgam,ffcdf,squad);
		elseif exp == 1
			% Run MC1 - MS tracking, full sky, pre-determined downwelling photons
			[dsto(iset).outcount(:,:,iwave),dsto(iset).outR(iwave),dsto(iset).outrrs(:,:,iwave),dsto(iset).outscount(:,:,iwave),dsto(iset).tscount(iwave)] = fwd_mc_deep1(n_zen,n_azi,cdep,aw(iwave),bw(iwave),dsto(iset).ag(iwave),dsto(iset).bp(iwave),ffgam,ffcdf,squad,drad_ang,drad_nph(:,iwave),pcount(iwave));
			
		elseif exp == 2
			% Run MC2 - black tank, limited FOV uniform sky
			[dsto(iset).outcount(:,:,iwave),dsto(iset).outR(iwave),dsto(iset).outrrs(:,:,iwave)] = fwd_mc_deep2(nph,n_zen,n_azi,aw(iwave),bw(iwave),dsto(iset).ag(iwave),dsto(iset).bp(iwave),ffgam,ffcdf,squad);
		   
	   elseif exp == 3
			% Run MC3 - black tank, pre-determined downwelling photons
			[dsto(iset).outcount(:,:,iwave),dsto(iset).outR(iwave),dsto(iset).outrrs(:,:,iwave)] = fwd_mc_deep3(n_zen,n_azi,aw(iwave),bw(iwave),dsto(iset).ag(iwave),dsto(iset).bp(iwave),ffgam,ffcdf,squad,drad_ang,drad_nph(:,iwave),pcount(iwave));

    end
    toc
end
toc


