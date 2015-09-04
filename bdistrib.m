function bdistrib()

% This matlab script performs the same steps to the fortran program. The
% fortran and matlab versions are completely independent of each other. For
% identical inputs, they should give identical outputs to within roundoff
% error.

clear

basis_set_option = 2;
% 1 = sines only
% 2 = cosines only
% 3 = both sines and cosines

weight_option = 2;
% 1: w = 1/(nfp * |N|) so Fourier functions are basis functions.
% 2: w = 1/area so each basis function is a linear combination of Fourier functions.

% Resolution parameters:
% **********************************

nu_plasma = 63;
nu_middle = 51;
nu_outer  = 49;

nv_plasma = 48;
nv_middle = 47;
nv_outer  = 40;

mpol_plasma = 10;
mpol_middle = 9;
mpol_outer  = 8;

ntor_plasma = 7;
ntor_middle = 6;
ntor_outer  = 5;

%{
nu_plasma = 40;
nu_middle = 47;
nu_outer  = 48;

nv_plasma = 49;
nv_middle = 51;
nv_outer  = 63;

mpol_plasma = 5;
mpol_middle = 6;
mpol_outer  = 7;

ntor_plasma = 8;
ntor_middle = 9;
ntor_outer  = 10;
%}

% Options for the shape of the plasma surface:
% **********************************
geometry_option_plasma = 2;
woutFilename = 'C:\Users\landreman\Box Sync\MATLAB\20150601-01 Sfincs version 3\equilibria\wout_w7x_standardConfig.nc';
R0_plasma = 5.0;
a_plasma = 0.5;
nfp_imposed = 4;

% Options for the shape of the middle surface:
% **********************************
geometry_option_middle = 3;
R0_middle = 3.0;
a_middle = 1.0;
separation_middle = 0.35;
nescin_filename_middle = 'nescin.w7x_standardConfig_separation0.3';

% Options for the shape of the outer surface:
% **********************************
geometry_option_outer = 3;
R0_outer = 6.0;
a_outer = 1.5;
separation_outer = 0.7;
nescin_filename_outer = 'nescin.w7x_standardConfig_separation0.6';

% Options for the SVDs and pseudo-inverse:
% **********************************
pseudoinverse_thresholds = [1e-12];

n_singular_vectors_to_save = 12;

inverse_option = 0;
% 0 = Pseudo-inverse with customizable threshold, analogous to the implementation in the fortran version.
% 1 = Matlab "/"
% 2 = Matlab inv()
% 3 = Matlab pinv()

% Plotting options:
% **********************************

plot3DFigure = true;
%plot3DFigure = false;

%plotGrids = true;
plotGrids = false;

%plotVectors = true;
plotVectors = false;

%stopAfterInitialPlots = true;
stopAfterInitialPlots = false;

figureOffset = 0;

% Options related to checking fortran version
% *******************************************

compareToFortran = true;
%compareToFortran = false;

fortranNcFilename = 'C:\Users\landreman\Box Sync\MATLAB\bdistrib_out.compareToMatlab.nc';

fortranComparisonThreshhold = 1e-7;

% *************************************************************************
% *************************************************************************
% End of input parameters.
% *************************************************************************
% *************************************************************************

    function compareVariableToFortran(variableName, varargin)
        % Specify 'abs' as an argument to compare the absolute values.
        % This is useful for the singular vectors, which are only defined
        % up to a sign in practice.
        if ~ compareToFortran
            return
        end
        try
            fortranVariable = double(ncread(fortranNcFilename,variableName));
        catch
            fprintf(['*** Variable ',variableName,' does not exist in the fortran output file.\n'])
            return
        end
        matlabVariable = eval(variableName);
        assignin('base',[variableName,'_m'],matlabVariable)
        assignin('base',[variableName,'_f'],fortranVariable)
        if isrow(matlabVariable)
            matlabVariable = matlabVariable(:);
        end
        if isrow(fortranVariable)
            fortranVariable = fortranVariable(:);
        end
        try
            % Next lines will cause an exception if sizes are different:
            if nargin>1 && strcmp(varargin{1},'abs')
                differences = abs(abs(matlabVariable) - abs(fortranVariable)) > fortranComparisonThreshhold;
            else
                differences = abs(matlabVariable - fortranVariable) > fortranComparisonThreshhold;
            end
            if any(any(any(differences)))
                fprintf(['*** Variable ',variableName,' is the same size Matlab and fortran but differs in value. max|diff|=%g\n'],max(max(max(differences))))
            else
                fprintf(['    Variable ',variableName,' is the same in Matlab and fortran.\n'])
            end
        catch
            fprintf(['*** Variable ',variableName,' is a different size between Matlab and fortran.\n'])
        end
    end

compareVariableToFortran('nu_plasma')
compareVariableToFortran('nv_plasma')
compareVariableToFortran('nu_middle')
compareVariableToFortran('nv_middle')
compareVariableToFortran('nu_outer')
compareVariableToFortran('nv_outer')
compareVariableToFortran('geometry_option_plasma')
compareVariableToFortran('geometry_option_middle')
compareVariableToFortran('geometry_option_outer')
compareVariableToFortran('pseudoinverse_thresholds')

% *********************************************
% Set up Fourier arrays
% *********************************************

    function [mnmax, xm, xn] = setupFourierArrays(mpol,ntor)
        % xm is non-negative, while xn can be negative
        % xn is the rapidly increasing variable.
        
        % When xm=0, xn=1..ntor.
        % When xm>0, xn=-ntor..ntor.
        mnmax = ntor + mpol*(2*ntor+1);
        
        xm = zeros(mnmax,1);
        xn = zeros(mnmax,1);
        xn(1:ntor) = 1:ntor;
        nextIndex = ntor+1;
        for m = 1:mpol
            indices = nextIndex:(nextIndex+2*ntor);
            xm(indices) = m;
            xn(indices) = (-ntor):ntor;
            nextIndex = nextIndex + 2*ntor+1;
        end
    end

[mnmax_plasma, xm_plasma, xn_plasma] = setupFourierArrays(mpol_plasma, ntor_plasma);
[mnmax_middle, xm_middle, xn_middle] = setupFourierArrays(mpol_middle, ntor_middle);
[mnmax_outer, xm_outer, xn_outer]    = setupFourierArrays(mpol_outer, ntor_outer);

compareVariableToFortran('mpol_plasma')
compareVariableToFortran('ntor_plasma')
compareVariableToFortran('mpol_middle')
compareVariableToFortran('ntor_middle')
compareVariableToFortran('mpol_outer')
compareVariableToFortran('ntor_outer')
compareVariableToFortran('mnmax_plasma')
compareVariableToFortran('mnmax_middle')
compareVariableToFortran('mnmax_outer')
compareVariableToFortran('xn_plasma')
compareVariableToFortran('xm_plasma')
compareVariableToFortran('xn_plasma')
compareVariableToFortran('xm_middle')
compareVariableToFortran('xn_middle')
compareVariableToFortran('xm_outer')
compareVariableToFortran('xn_outer')
compareVariableToFortran('basis_set_option')
compareVariableToFortran('weight_option')

% *********************************************
% Initialize the plasma surface:
% *********************************************

switch geometry_option_plasma
    case {0,1}
        % Plain axisymmetric circular torus
        nfp = nfp_imposed;
        mnmax = 2;
        xm = [0,1];
        xn = [0,0];
        rmnc = [R0_plasma; a_plasma];
        zmns = [0; a_plasma];
        whichSurface = 2;
        Rmajor_p = R0_plasma;
        
    case {2}
        % Load flux surface info from VMEC
        filename = woutFilename;
        ns = double(ncread(filename,'ns'));
        Rmajor_p = double(ncread(filename,'Rmajor_p'));
        nfp = double(ncread(filename,'nfp'));
        xm = double(ncread(filename,'xm'));
        xn = double(ncread(filename,'xn'));
        xm_nyq = double(ncread(filename,'xm_nyq'));
        xn_nyq = double(ncread(filename,'xn_nyq'));
        rmnc = double(ncread(filename,'rmnc'));
        zmns = double(ncread(filename,'zmns'));
        bmnc = double(ncread(filename,'bmnc'));
        mnmax = double(ncread(filename,'mnmax'));
        mnmax_nyq = double(ncread(filename,'mnmax_nyq'));
        
        whichSurface = ns; % Choose the outermost surface
        % Discard the other surfaces:
        rmnc = rmnc(:,whichSurface);
        zmns = zmns(:,whichSurface);
        
    otherwise
        error('Invalid setting for geometry_option_plasma')
end

nvl_plasma = nv_plasma * nfp;
nvl_middle = nv_middle * nfp;
nvl_outer = nv_outer * nfp;

u_plasma = linspace(0,1,nu_plasma+1);
u_plasma(end) = [];
v_plasma = linspace(0,1,nv_plasma+1);
v_plasma(end) = [];
vl_plasma = linspace(0,nfp,nvl_plasma+1);
vl_plasma(end) = [];
[vl_plasma_2D, u_plasma_2D] = meshgrid(vl_plasma, u_plasma);

x = zeros(nu_plasma,nvl_plasma);
y = zeros(nu_plasma,nvl_plasma);
z = zeros(nu_plasma,nvl_plasma);

dxdu = zeros(nu_plasma,nvl_plasma);
dydu = zeros(nu_plasma,nvl_plasma);
dzdu = zeros(nu_plasma,nvl_plasma);

dxdv = zeros(nu_plasma,nvl_plasma);
dydv = zeros(nu_plasma,nvl_plasma);
dzdv = zeros(nu_plasma,nvl_plasma);

for i=1:mnmax
    angle = xm(i)*2*pi*u_plasma_2D-xn(i)*2*pi/nfp*vl_plasma_2D;
    angle2 = 2*pi*vl_plasma_2D/nfp;
    
    x = x + rmnc(i)*cos(angle).*cos(angle2);
    y = y + rmnc(i)*cos(angle).*sin(angle2);
    z = z + zmns(i)*sin(angle);
    
    dxdu = dxdu - 2*pi*xm(i)*rmnc(i)*sin(angle).*cos(angle2);
    dydu = dydu - 2*pi*xm(i)*rmnc(i)*sin(angle).*sin(angle2);
    dzdu = dzdu + 2*pi*xm(i)*zmns(i)*cos(angle);
    
    dxdv = dxdv + 2*pi/nfp*xn(i)*rmnc(i)*sin(angle).*cos(angle2) ...
        - 2*pi/nfp*rmnc(i)*cos(angle).*sin(angle2);
    dydv = dydv + 2*pi/nfp*xn(i)*rmnc(i)*sin(angle).*sin(angle2) ...
        + 2*pi/nfp*rmnc(i)*cos(angle).*cos(angle2);
    dzdv = dzdv - 2*pi/nfp*xn(i)*zmns(i)*cos(angle);
end

Nx = dydv .* dzdu - dzdv .* dydu;
Ny = dzdv .* dxdu - dxdv .* dzdu;
Nz = dxdv .* dydu - dydv .* dxdu;
norm_normal_plasma = sqrt(Nx.*Nx + Ny.*Ny + Nz.*Nz);
du_plasma = u_plasma(2)-u_plasma(1);
dv_plasma = v_plasma(2)-v_plasma(1);
area_plasma = sum(sum(norm_normal_plasma)) * du_plasma * dv_plasma;

r_plasma = zeros(3, nu_plasma, nvl_plasma);
drdu_plasma = zeros(3, nu_plasma, nvl_plasma);
drdv_plasma = zeros(3, nu_plasma, nvl_plasma);
normal_plasma = zeros(3, nu_plasma, nvl_plasma);

r_plasma(1,:,:) = x;
r_plasma(2,:,:) = y;
r_plasma(3,:,:) = z;
drdu_plasma(1,:,:) = dxdu;
drdu_plasma(2,:,:) = dydu;
drdu_plasma(3,:,:) = dzdu;
drdv_plasma(1,:,:) = dxdv;
drdv_plasma(2,:,:) = dydv;
drdv_plasma(3,:,:) = dzdv;
normal_plasma(1,:,:) = Nx;
normal_plasma(2,:,:) = Ny;
normal_plasma(3,:,:) = Nz;

compareVariableToFortran('nfp')
compareVariableToFortran('u_plasma')
compareVariableToFortran('v_plasma')
compareVariableToFortran('vl_plasma')

compareVariableToFortran('r_plasma')
compareVariableToFortran('drdu_plasma')
compareVariableToFortran('drdv_plasma')
compareVariableToFortran('normal_plasma')
compareVariableToFortran('norm_normal_plasma')


% *********************************************
% Initialize the middle and outer surfaces:
% *********************************************

    function [u, v, vl, u_2D, vl_2D, r, drdu, drdv, normal, norm_normal, area] ...
            = initSurface(nu, nv, geometry_option, R0, a, separation, nescin_filename)
        
        nvl = nv*nfp;
        u = linspace(0,1,nu+1);
        u(end) = [];
        v = linspace(0,1,nv+1);
        v(end) = [];
        vl = linspace(0,nfp,nvl+1);
        vl(end) = [];
        [vl_2D, u_2D] = meshgrid(vl, u);
        
        x = zeros(size(u_2D));
        y = zeros(size(u_2D));
        z = zeros(size(u_2D));
        dxdu = zeros(size(u_2D));
        dydu = zeros(size(u_2D));
        dzdu = zeros(size(u_2D));
        dxdv = zeros(size(u_2D));
        dydv = zeros(size(u_2D));
        dzdv = zeros(size(u_2D));
        
        switch(geometry_option)
            case {0,1}
                if geometry_option == 0
                    R0_to_use = Rmajor_p;
                else
                    R0_to_use = R0;
                end
                
                x = (R0_to_use + a * cos(2*pi*u_2D)) .* cos(2*pi*vl_2D/nfp);
                y = (R0_to_use + a * cos(2*pi*u_2D)) .* sin(2*pi*vl_2D/nfp);
                z = a * sin(2*pi*u_2D);
                
                dxdu = -a * 2 * pi * sin(2*pi*u_2D) .* cos(2*pi*vl_2D/nfp);
                dydu = -a * 2 * pi * sin(2*pi*u_2D) .* sin(2*pi*vl_2D/nfp);
                dzdu = a * 2 * pi * cos(2*pi*u_2D);
                
                dxdv = -2*pi/nfp*(R0_to_use + a * cos(2*pi*u_2D)) .* sin(2*pi*vl_2D/nfp);
                dydv =  2*pi/nfp*(R0_to_use + a * cos(2*pi*u_2D)) .* cos(2*pi*vl_2D/nfp);
                dzdv = zeros(size(u_2D));
                
            case 2
                error('geometry_option = 2 is not yet implemented for middle and outer surfaces.')
                
            case 3
                % Read coil surface from nescin file
                
                fid = fopen(nescin_filename,'r');
                search_string = '------ Current Surface';
                while true
                    line = fgetl(fid);
                    if strncmp(line,search_string,numel(search_string))
                        break
                    end
                end
                line = fgetl(fid); %Number of fourier modes in table
                line = fgetl(fid);
                mnmax_coil = sscanf(line,'%d');
                fprintf('  Reading %d modes from nescin file %s\n',mnmax_coil,nescin_filename)
                line = fgetl(fid); %Table of fourier coefficients
                line = fgetl(fid); %m,n,crc2,czs2,crs2,czc2
                xm_coil = zeros(mnmax_coil,1);
                xn_coil = zeros(mnmax_coil,1);
                rmnc_coil = zeros(mnmax_coil,1);
                zmns_coil = zeros(mnmax_coil,1);
                for i=1:mnmax_coil
                    line = fgetl(fid);
                    data = sscanf(line,'%d %d %g %g %g %g %g %g');
                    xm_coil(i) = data(1);
                    xn_coil(i) = data(2);
                    rmnc_coil(i) = data(3);
                    zmns_coil(i) = data(4);
                end
                fclose(fid);
                % Done reading nescin file.
                                
                for i = 1:mnmax_coil
                    angle = xm_coil(i)*2*pi*u_2D + xn_coil(i)*2*pi*vl_2D;
                    angle2 = 2*pi*vl_2D/nfp;
                    
                    x = x + rmnc_coil(i)*cos(angle).*cos(angle2);
                    y = y + rmnc_coil(i)*cos(angle).*sin(angle2);
                    z = z + zmns_coil(i)*sin(angle);
                    
                    dxdu = dxdu - 2*pi*xm_coil(i)*rmnc_coil(i)*sin(angle).*cos(angle2);
                    dydu = dydu - 2*pi*xm_coil(i)*rmnc_coil(i)*sin(angle).*sin(angle2);
                    dzdu = dzdu + 2*pi*xm_coil(i)*zmns_coil(i)*cos(angle);
                    
                    dxdv = dxdv - 2*pi*xn_coil(i)*rmnc_coil(i)*sin(angle).*cos(angle2) ...
                        - 2*pi/nfp*rmnc_coil(i)*cos(angle).*sin(angle2);
                    dydv = dydv - 2*pi*xn_coil(i)*rmnc_coil(i)*sin(angle).*sin(angle2) ...
                        + 2*pi/nfp*rmnc_coil(i)*cos(angle).*cos(angle2);
                    dzdv = dzdv + 2*pi*xn_coil(i)*zmns_coil(i)*cos(angle);
                    
                end
                
            otherwise
                error('Invalid geometry_option')
        end
        
        Nx = dydv .* dzdu - dzdv .* dydu;
        Ny = dzdv .* dxdu - dxdv .* dzdu;
        Nz = dxdv .* dydu - dydv .* dxdu;
        
        r = zeros(3, nu, nvl);
        drdu = zeros(3, nu, nvl);
        drdv = zeros(3, nu, nvl);
        normal = zeros(3, nu, nvl);
        
        r(1,:,:) = x;
        r(2,:,:) = y;
        r(3,:,:) = z;
        drdu(1,:,:) = dxdu;
        drdu(2,:,:) = dydu;
        drdu(3,:,:) = dzdu;
        drdv(1,:,:) = dxdv;
        drdv(2,:,:) = dydv;
        drdv(3,:,:) = dzdv;
        normal(1,:,:) = Nx;
        normal(2,:,:) = Ny;
        normal(3,:,:) = Nz;
        
        norm_normal = sqrt(Nx.*Nx + Ny.*Ny + Nz.*Nz);
        du = u(2)-u(1);
        dv = v(2)-v(1);
        area = sum(sum(norm_normal)) * du * dv;
    end

tic
fprintf('Initializing middle surface.\n')
[u_middle, v_middle, vl_middle, u_middle_2D, vl_middle_2D, r_middle, drdu_middle, drdv_middle, normal_middle, norm_normal_middle, area_middle] ...
    = initSurface(nu_middle, nv_middle, geometry_option_middle, R0_middle, a_middle, separation_middle, nescin_filename_middle);
fprintf('Done. Took %g seconds.\n',toc)

tic
fprintf('Initializing outer surface.\n')
[u_outer, v_outer, vl_outer, u_outer_2D, vl_outer_2D, r_outer, drdu_outer, drdv_outer, normal_outer, norm_normal_outer, area_outer] ...
    = initSurface(nu_outer, nv_outer, geometry_option_outer, R0_outer, a_outer, separation_outer, nescin_filename_outer);
fprintf('Done. Took %g seconds.\n',toc)

compareVariableToFortran('u_middle')
compareVariableToFortran('v_middle')
compareVariableToFortran('vl_middle')
compareVariableToFortran('r_middle')
compareVariableToFortran('drdu_middle')
compareVariableToFortran('drdv_middle')
compareVariableToFortran('normal_middle')
compareVariableToFortran('norm_normal_middle')

compareVariableToFortran('u_outer')
compareVariableToFortran('v_outer')
compareVariableToFortran('vl_outer')
compareVariableToFortran('r_outer')
compareVariableToFortran('drdu_outer')
compareVariableToFortran('drdv_outer')
compareVariableToFortran('normal_outer')
compareVariableToFortran('norm_normal_outer')
compareVariableToFortran('area_plasma')
compareVariableToFortran('area_middle')
compareVariableToFortran('area_outer')

% *********************************************
% Make 3D figure of surfaces
% *********************************************

%r_plasma = ncread(fortranNcFilename,'r_plasma');

if plot3DFigure
    r_plasma_toplot = r_plasma;
    r_middle_toplot = r_middle;
    r_outer_toplot = r_outer;
    
    % "Rotate" in theta so the seam in the plot is on the bottom
    nshift = round(nu_plasma*0.25);
    r_plasma_toplot = circshift(r_plasma_toplot, [0,nshift,0]);
    nshift = round(nu_middle*0.25);
    r_middle_toplot = circshift(r_middle_toplot, [0,nshift,0]);
    nshift = round(nu_outer*0.25);
    r_outer_toplot = circshift(r_outer_toplot, [0,nshift,0]);
    
    % Close surfaces for plotting:
    r_plasma_toplot(:,:,end+1) = r_plasma_toplot(:,:,1);
    r_plasma_toplot(:,end+1,:) = r_plasma_toplot(:,1,:);
    
    % For middle and outer surfaces, close in u, but don't bother closing
    % in v:
    r_middle_toplot(:,end+1,:) = r_middle_toplot(:,1,:);
    r_outer_toplot(:,end+1,:)  = r_outer_toplot(:,1,:);
    
    
    
    mask = vl_middle < 0.7*nfp;
    r_middle_toplot = r_middle_toplot(:,:,mask);
    
    mask = (vl_outer > 0.15*nfp) & (vl_outer < 0.55*nfp);
    r_outer_toplot = r_outer_toplot(:,:,mask);
    
    figure(1 + figureOffset)
    clf
    set(gcf,'Color','w')
    faceColor = [1,0,0];
    surf(squeeze(r_plasma_toplot(1,:,:)),squeeze(r_plasma_toplot(2,:,:)),squeeze(r_plasma_toplot(3,:,:)),'FaceColor',faceColor,'EdgeColor','none')
    hold on
    if plotGrids
        plot3(squeeze(r_plasma_toplot(1,:,:)),squeeze(r_plasma_toplot(2,:,:)),squeeze(r_plasma_toplot(3,:,:)),'.r')
    end
    daspect([1,1,1])
    %shading interp
    axis vis3d
    hold on
    
    faceColor = [1,0,1];
    %faceColor = [0,1,0];
    surf(squeeze(r_middle_toplot(1,:,:)),squeeze(r_middle_toplot(2,:,:)),squeeze(r_middle_toplot(3,:,:)),'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.75)
    
    faceColor = [0,0,1];
    surf(squeeze(r_outer_toplot(1,:,:)),squeeze(r_outer_toplot(2,:,:)),squeeze(r_outer_toplot(3,:,:)),'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.75)
    
    if plotGrids
        plot3(squeeze(r_middle_toplot(1,:,:)),squeeze(r_middle_toplot(2,:,:)),squeeze(r_middle_toplot(3,:,:)),'.m')
        plot3(squeeze(r_outer_toplot(1,:,:)), squeeze(r_outer_toplot(2,:,:)), squeeze(r_outer_toplot(3,:,:)),'.b')
    end
    %surf(X_coil,Y_coil,Z_coil,'FaceColor',faceColor,'EdgeColor','none')
    %surf(X_coil,Y_coil,Z_coil,'FaceColor',faceColor,'EdgeColor','none','FaceAlpha',0.75)
    light
    lighting gouraud
    zoom(1.6)
    campos([  574.9370 -457.0244  424.3304])
    camva(1.0271)
    axis off
    
end

if stopAfterInitialPlots
    return
end


% *********************************************
% If needed, compute conversions between Fourier basis and the basis
% functions.
% *********************************************

switch basis_set_option
    case {1,2}
        num_basis_functions_plasma = mnmax_plasma;
        num_basis_functions_middle = mnmax_middle;
        num_basis_functions_outer  = mnmax_outer;
    case {3}
        num_basis_functions_plasma = mnmax_plasma * 2;
        num_basis_functions_middle = mnmax_middle * 2;
        num_basis_functions_outer  = mnmax_outer  * 2;
    otherwise
        error('Invalid value for basis_set_option')
end

% *********************************************
% Compute mutual inductance matrices
% *********************************************

mu0 = 4*pi*(1e-7);
[v_outer_2D, u_outer_2D] = meshgrid(v_outer,u_outer);
du_outer = u_outer(2)-u_outer(1);
dv_outer = v_outer(2)-v_outer(1);

    function [inductanceMatrix, Fourier_to_basis, basis_to_Fourier] = ...
            computeInductanceMatrix(r, normal, norm_normal, u, v, mnmax, xm, xn, num_basis_functions, area)
        nu = size(r,2);
        nv = size(r,3)/nfp;
        if round(nv) ~= nv
            error('Something went wrong.')
        end
        
        [v2D, u2D] = meshgrid(v,u);
        du = u(2)-u(1);
        dv = v(2)-v(1);
        
        tic1 = tic;
        xToFourier = zeros(num_basis_functions,nu*nv);
        xToFourier_outer = zeros(nu_outer*nv_outer,num_basis_functions_outer);
        switch basis_set_option
            case {1}
                % sines only
                for imn = 1:mnmax
                    xToFourier(imn,:) = reshape(sin(2*pi*(xm(imn)*u2D + xn(imn)*v2D)), [nu*nv,1])';
                end
                for imn = 1:mnmax_outer
                    xToFourier_outer(:,imn) = reshape(sin(2*pi*(xm_outer(imn)*u_outer_2D + xn_outer(imn)*v_outer_2D)), [nu_outer*nv_outer,1]);
                end
            case {2}
                % cosines only
                for imn = 1:mnmax
                    xToFourier(imn,:) = reshape(cos(2*pi*(xm(imn)*u2D + xn(imn)*v2D)), [nu*nv,1])';
                end
                for imn = 1:mnmax_outer
                    xToFourier_outer(:,imn) = reshape(cos(2*pi*(xm_outer(imn)*u_outer_2D + xn_outer(imn)*v_outer_2D)), [nu_outer*nv_outer,1]);
                end
            case {3}
                % Both sines and cosines
                for imn = 1:mnmax
                    xToFourier(imn,:) = reshape(sin(2*pi*(xm(imn)*u2D + xn(imn)*v2D)), [nu*nv,1])';
                    xToFourier(imn+mnmax,:) = reshape(cos(2*pi*(xm(imn)*u2D + xn(imn)*v2D)), [nu*nv,1])';
                end
                for imn = 1:mnmax_outer
                    xToFourier_outer(:,imn) = reshape(sin(2*pi*(xm_outer(imn)*u_outer_2D + xn_outer(imn)*v_outer_2D)), [nu_outer*nv_outer,1]);
                    xToFourier_outer(:,imn+mnmax_outer) = reshape(cos(2*pi*(xm_outer(imn)*u_outer_2D + xn_outer(imn)*v_outer_2D)), [nu_outer*nv_outer,1]);
                end
        end
        fprintf('  xToFourier: %g\n',toc(tic1))
        
        xToFourier       = sqrt(2) * xToFourier;
        xToFourier_outer = sqrt(2) * xToFourier_outer;
        
        switch weight_option
            case 1
                Fourier_to_basis = eye(num_basis_functions);
                basis_to_Fourier = eye(num_basis_functions);
                 
            case 2
                %norm_normal = ones(nu,nv); % This line is for testing only!
                tic
                % Form area integral of f_i * f_j (product of 2 Fourier functions)
                if false
                    % Slow but transparent method
                    A = zeros(num_basis_functions);
                    for imn_row = 1:num_basis_functions
                        switch basis_set_option
                            case {1}
                                Fourier_row = sin(2*pi*(xm(imn_row)*u2D + xn(imn_row)*v2D));
                            case {2}
                                Fourier_row = cos(2*pi*(xm(imn_row)*u2D + xn(imn_row)*v2D));
                            case {3}
                                if imn_row <= mnmax
                                    Fourier_row = sin(2*pi*(xm(imn_row)*u2D + xn(imn_row)*v2D));
                                else
                                    Fourier_row = cos(2*pi*(xm(imn_row-mnmax)*u2D + xn(imn_row-mnmax)*v2D));
                                end
                        end
                        for imn_col = 1:num_basis_functions
                            switch basis_set_option
                                case {1}
                                    Fourier_col = sin(2*pi*(xm(imn_col)*u2D + xn(imn_col)*v2D));
                                case {2}
                                    Fourier_col = cos(2*pi*(xm(imn_col)*u2D + xn(imn_col)*v2D));
                                case {3}
                                    if imn_col <= mnmax
                                        Fourier_col = sin(2*pi*(xm(imn_col)*u2D + xn(imn_col)*v2D));
                                    else
                                        Fourier_col = cos(2*pi*(xm(imn_col-mnmax)*u2D + xn(imn_col-mnmax)*v2D));
                                    end
                            end
                            A(imn_row,imn_col) = sum(sum(Fourier_row .* Fourier_col .* norm_normal(:,1:nv))) * du * dv * nfp / area * 2; % 2 is from sqrt(2) factor in each Fourier function
                        end
                    end
                else
                    % Faster method
                    temp = reshape(norm_normal(:,1:nv), [nu*nv,1]);
                    A = xToFourier * diag(temp) * (xToFourier') * du * dv * nfp / area;
                end
                basis_to_Fourier = chol(A,'lower');
                Fourier_to_basis = inv(basis_to_Fourier);
                fprintf('  Fourier_to_basis: %g\n',toc(tic1))
                
            otherwise
                error('Invalid weight_option')
        end
        

        tic1 = tic;
        vIndices = 1:nv;
        inductanceMatrix_xBasis = zeros(nu*nv, nu_outer*nv_outer);
        for iu_outer = 1:nu_outer
            for iv_outer = 1:nv_outer
                index_outer = (iv_outer-1)*nu_outer + iu_outer;
                for l_outer = 0:(nfp-1)
                    ivl_outer = iv_outer + l_outer*nv_outer;
                    dx = r(1,:,vIndices) - r_outer(1,iu_outer,ivl_outer);
                    dy = r(2,:,vIndices) - r_outer(2,iu_outer,ivl_outer);
                    dz = r(3,:,vIndices) - r_outer(3,iu_outer,ivl_outer);
                    dr2 = dx.*dx + dy.*dy + dz.*dz;
                    denominator = dr2 .* sqrt(dr2);
                    temp = (normal(1,:,vIndices)*normal_outer(1,iu_outer,ivl_outer) ...
                        +   normal(2,:,vIndices)*normal_outer(2,iu_outer,ivl_outer) ...
                        +   normal(3,:,vIndices)*normal_outer(3,iu_outer,ivl_outer) ...
                        - (3./dr2) .* (dx .* normal(1,:,vIndices) + dy .* normal(2,:,vIndices) + dz .* normal(3,:,vIndices)) ...
                        .* (dx * normal_outer(1,iu_outer,ivl_outer) + dy * normal_outer(2,iu_outer,ivl_outer) + dz * normal_outer(3,iu_outer,ivl_outer))) ./ denominator;
                    inductanceMatrix_xBasis(:,index_outer) = inductanceMatrix_xBasis(:,index_outer) + ...
                        reshape(temp, [nu*nv,1]);
                end
            end
        end
        fprintf('  xBasis: %g\n',toc(tic1))
        
        tic1 = tic;
        inductanceMatrix = (nfp * mu0/(4*pi) * du * dv * du_outer * dv_outer) * xToFourier * inductanceMatrix_xBasis * xToFourier_outer;
        if weight_option>1
            inductanceMatrix = Fourier_to_basis * inductanceMatrix;
            %inductanceMatrix = basis_to_Fourier \ inductanceMatrix;
        end
        fprintf('  matmul: %g\n',toc(tic1))
    end

tic
fprintf('Building mutual inductance matrix between the plasma and outer surfaces.\n')
[inductance_plasma, Fourier_to_basis_plasma, basis_to_Fourier_plasma] = ...
    computeInductanceMatrix(r_plasma, normal_plasma, norm_normal_plasma, u_plasma, v_plasma, mnmax_plasma, xm_plasma, xn_plasma, num_basis_functions_plasma, area_plasma);
fprintf('Done. Took %g seconds.\n',toc)

tic
fprintf('Building mutual inductance matrix between the middle and outer surfaces.\n')
[inductance_middle, Fourier_to_basis_middle, basis_to_Fourier_middle] = ...
    computeInductanceMatrix(r_middle, normal_middle, norm_normal_middle, u_middle, v_middle, mnmax_middle, xm_middle, xn_middle, num_basis_functions_middle, area_middle);
fprintf('Done. Took %g seconds.\n',toc)

if weight_option>1
    compareVariableToFortran('basis_to_Fourier_plasma')
    compareVariableToFortran('basis_to_Fourier_middle')
end

compareVariableToFortran('inductance_plasma')
compareVariableToFortran('inductance_middle')
%return
% *********************************************
% Compute SVD of the two inductance matrices
% *********************************************

tic
fprintf('Computing SVD of the mutual inductance matrix between the plasma and outer surfaces.\n')
[svd_u_inductance_plasma, svd_s_inductance_plasma, svd_v_inductance_plasma] = svd(inductance_plasma);
svd_s_inductance_plasma = diag(svd_s_inductance_plasma);
fprintf('Done. Took %g seconds.\n',toc)

tic
fprintf('Computing SVD of the mutual inductance matrix between the middle and outer surfaces.\n')
[svd_u_inductance_middle, svd_s_inductance_middle, svd_v_inductance_middle] = svd(inductance_middle);
svd_s_inductance_middle_size = size(svd_s_inductance_middle);
svd_s_inductance_middle = diag(svd_s_inductance_middle);
fprintf('Done. Took %g seconds.\n',toc)

figure(figureOffset + 2)
clf
semilogy(svd_s_inductance_middle,'.m','DisplayName','Inductance matrix between middle and outer surfaces')
hold on
semilogy(svd_s_inductance_plasma,'.r','DisplayName','Inductance matrix between plasma and outer surfaces')
%set(gca,'xGrid','on','yGrid','on')
title('Singular values')

compareVariableToFortran('svd_s_inductance_plasma')
compareVariableToFortran('svd_s_inductance_middle')

% *********************************************
% Compute transfer matrix
% *********************************************

n_pseudoinverse_thresholds = numel(pseudoinverse_thresholds);
n_singular_values_retained = zeros(n_pseudoinverse_thresholds,1);
svd_s_transferMatrix = zeros(min([num_basis_functions_plasma,num_basis_functions_middle]), n_pseudoinverse_thresholds);
svd_u_transferMatrix = zeros(num_basis_functions_plasma,n_singular_vectors_to_save, n_pseudoinverse_thresholds);
svd_v_transferMatrix = zeros(num_basis_functions_middle,n_singular_vectors_to_save, n_pseudoinverse_thresholds);

compareVariableToFortran('n_pseudoinverse_thresholds')
compareVariableToFortran('pseudoinverse_thresholds')
compareVariableToFortran('n_singular_vectors_to_save')


colors = [0.9,0.6,0;
    0,0.7,0;
    0,0,1;
    0,0,0;];

if inverse_option==0
    Mpo_V = inductance_plasma * svd_v_inductance_middle;
end

for i = 1:n_pseudoinverse_thresholds
    % Decide how many singular values to keep:
    threshold = pseudoinverse_thresholds(i);
    keepMask = (svd_s_inductance_middle/svd_s_inductance_middle(1)) >= threshold;
    n_singular_values_retained(i) = sum(keepMask);
    fprintf('Forming transfer matrix for pinv threshold %g: retaining %d singular values.\n',threshold, n_singular_values_retained(i))
    tic
    
    switch inverse_option
        case 0
            fprintf('Using customizable pseudo-inverse to invert M_mo.\n')
            % Build the rectangular matrix of inverse singular values:
            inverseSingularValues = 1./svd_s_inductance_middle;
            inverseSingularValues(~keepMask) = 0;
            inverseSingularValues = inverseSingularValues(:);
            diagonalPart = spdiags(inverseSingularValues, 0, svd_s_inductance_middle_size(2), svd_s_inductance_middle_size(1));
            
            % Put the pieces together:
            transferMatrix = Mpo_V * diagonalPart * (svd_u_inductance_middle');
        case 1
            fprintf('Using Matlab''s "/" to invert M_mo.\n')
            transferMatrix = inductance_plasma / inductance_middle;
        case 2
            fprintf('Using Matlab''s "inv()" to invert M_mo.\n')
            transferMatrix = inductance_plasma * inv(inductance_middle);
        case 3
            fprintf('Using Matlab''s "pinv()" to invert M_mo.\n')
            transferMatrix = inductance_plasma * pinv(inductance_middle);
        otherwise
            error('Invalid option for inverse_option')
    end
    fprintf('  Assembled transfer matrix. Took %g seconds.\n',toc)
    
    % Carry out SVD
    tic
    [svd_u_transferMatrix_pre, svd_s_transferMatrix_pre, svd_v_transferMatrix_pre] = svd(transferMatrix);
    svd_s_transferMatrix(:,i) = diag(svd_s_transferMatrix_pre);
    fprintf('  Done with SVD. Took %g seconds.\n',toc)
    svd_u_transferMatrix(:,:,i) = basis_to_Fourier_plasma * svd_u_transferMatrix_pre(:,1:n_singular_vectors_to_save);
    svd_v_transferMatrix(:,:,i) = basis_to_Fourier_middle * svd_v_transferMatrix_pre(:,1:n_singular_vectors_to_save);
    
    colorindex = 1+mod(i-1,size(colors,1));
    plot(svd_s_transferMatrix(:,i),'.','DisplayName',['Transfer matrix, th=',num2str(threshold)],'Color',colors(colorindex,:))
end

switch basis_set_option
    case {1}
        svd_u_transferMatrix_sin = svd_u_transferMatrix;
        svd_u_transferMatrix_cos = zeros(size(svd_u_transferMatrix));
        svd_v_transferMatrix_sin = svd_v_transferMatrix;
        svd_v_transferMatrix_cos = zeros(size(svd_v_transferMatrix));
    case {2}
        svd_u_transferMatrix_cos = svd_u_transferMatrix;
        svd_u_transferMatrix_sin = zeros(size(svd_u_transferMatrix));
        svd_v_transferMatrix_cos = svd_v_transferMatrix;
        svd_v_transferMatrix_sin = zeros(size(svd_v_transferMatrix));
    case {3}
        svd_u_transferMatrix_sin = svd_u_transferMatrix(1:mnmax_plasma,:,:);
        svd_u_transferMatrix_cos = svd_u_transferMatrix(mnmax_plasma+1:end,:,:);
        svd_v_transferMatrix_sin = svd_v_transferMatrix(1:mnmax_middle,:,:);
        svd_v_transferMatrix_cos = svd_v_transferMatrix(mnmax_middle+1:end,:,:);
end

legend show
set(legend,'Box','off','Location','southwest')

compareVariableToFortran('svd_s_transferMatrix')
compareVariableToFortran('svd_u_transferMatrix_sin','abs')
compareVariableToFortran('svd_u_transferMatrix_cos','abs')
compareVariableToFortran('svd_v_transferMatrix_sin','abs')
compareVariableToFortran('svd_v_transferMatrix_cos','abs')

% *********************************************
% Done with the main calculation.
% Now plot singular vectors.
% *********************************************

numContours = 20;
nextFigure = 2;

nu_plot = 102;
nv_plot = 100;
u_plot = linspace(0,1,nu_plot);
v_plot = linspace(0,1,nv_plot);
[v_plot_2D, u_plot_2D] = meshgrid(v_plot,u_plot);
FourierToX = zeros(nu_plot*nv_plot,num_basis_functions_plasma);
switch basis_set_option
    case {1}
        for imn = 1:mnmax_plasma
            FourierToX(:,imn) = reshape(sin(2*pi*(xm_plasma(imn)*u_plot_2D + xn_plasma(imn)*v_plot_2D)), [nu_plot*nv_plot,1]);
        end
    case {2}
        for imn = 1:mnmax_plasma
            FourierToX(:,imn) = reshape(cos(2*pi*(xm_plasma(imn)*u_plot_2D + xn_plasma(imn)*v_plot_2D)), [nu_plot*nv_plot,1]);
        end
    case {3}
        for imn = 1:mnmax_plasma
            FourierToX(:,imn) = reshape(sin(2*pi*(xm_plasma(imn)*u_plot_2D + xn_plasma(imn)*v_plot_2D)), [nu_plot*nv_plot,1]);
            FourierToX(:,imn+mnmax_plasma) = reshape(cos(2*pi*(xm_plasma(imn)*u_plot_2D + xn_plasma(imn)*v_plot_2D)), [nu_plot*nv_plot,1]);
        end
end

for whichThreshold = 1:n_pseudoinverse_thresholds
    nextFigure = nextFigure+1;
    figure(nextFigure)
    clf
    numCols = ceil(sqrt(n_singular_vectors_to_save));
    numRows = ceil(n_singular_vectors_to_save / numCols);
    
    for whichPlot = 1:n_singular_vectors_to_save
        subplot(numRows,numCols,whichPlot)
        data = reshape(FourierToX * svd_u_transferMatrix(:,whichPlot,whichThreshold),[nu_plot,nv_plot]);
        contourf(v_plot, u_plot, data, numContours,'EdgeColor','none')
        colorbar
        xlabel('v')
        ylabel('u')
        title({['Singular vector u ',num2str(whichPlot)],['s=',num2str(svd_s_transferMatrix(whichPlot,whichThreshold))]})
    end
    stringForTop = ['Efficiency-ordered B_n distributions on plasma surface (threshold=',num2str(pseudoinverse_thresholds(whichThreshold)),')'];
    annotation('textbox',[0 0.96 1 .04],'HorizontalAlignment','center',...
        'Interpreter','none','VerticalAlignment','bottom',...
        'FontSize',11,'LineStyle','none','String',stringForTop);
    
end


FourierToX = zeros(nu_plot*nv_plot,num_basis_functions_middle);
switch basis_set_option
    case {1}
        for imn = 1:mnmax_middle
            FourierToX(:,imn) = reshape(sin(2*pi*(xm_middle(imn)*u_plot_2D + xn_middle(imn)*v_plot_2D)), [nu_plot*nv_plot,1]);
        end
    case {2}
        for imn = 1:mnmax_middle
            FourierToX(:,imn) = reshape(cos(2*pi*(xm_middle(imn)*u_plot_2D + xn_middle(imn)*v_plot_2D)), [nu_plot*nv_plot,1]);
        end
    case {3}
        for imn = 1:mnmax_middle
            FourierToX(:,imn) = reshape(sin(2*pi*(xm_middle(imn)*u_plot_2D + xn_middle(imn)*v_plot_2D)), [nu_plot*nv_plot,1]);
            FourierToX(:,imn+mnmax_middle) = reshape(cos(2*pi*(xm_middle(imn)*u_plot_2D + xn_middle(imn)*v_plot_2D)), [nu_plot*nv_plot,1]);
        end
end

for whichThreshold = 1:n_pseudoinverse_thresholds
    nextFigure = nextFigure+1;
    figure(nextFigure)
    clf
    numCols = ceil(sqrt(n_singular_vectors_to_save));
    numRows = ceil(n_singular_vectors_to_save / numCols);
    
    for whichPlot = 1:n_singular_vectors_to_save
        subplot(numRows,numCols,whichPlot)
        data = reshape(FourierToX * svd_v_transferMatrix(:,whichPlot,whichThreshold),[nu_plot,nv_plot]);
        contourf(v_plot, u_plot, data, numContours,'EdgeColor','none')
        colorbar
        xlabel('v')
        ylabel('u')
        title(['Singular vector v ',num2str(whichPlot)])
    end
    stringForTop = ['Efficiency-ordered B_n distributions on middle surface (threshold=',num2str(pseudoinverse_thresholds(whichThreshold)),')'];
    annotation('textbox',[0 0.96 1 .04],'HorizontalAlignment','center',...
        'Interpreter','none','VerticalAlignment','bottom',...
        'FontSize',11,'LineStyle','none','String',stringForTop);
    
end

end
