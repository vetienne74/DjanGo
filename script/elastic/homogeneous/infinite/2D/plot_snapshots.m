clear all ;
close all ;

% PARAMETERS
FONT_SIZE  = 12 ;
FLOAT_SIZE = 'float32' ;

figure('Position',[100 100 700 700])

% FIGURE 1: DG vs FDM - NO ABC
%-----------------------------
COMPONENT = "vz" ;
FIG_NAME  = 'HOMO-INF' ;
FILE1     = 'django.config.py.fem.continuous.O6.elastic.O1.sponge.xml.snapshot' ;
SCHEME1   = 'DGSE-P6-SPONGE' ;

FILE2     = 'django.config.py.fdm.staggered.O4.elastic.O1.pml.xml.snapshot';
SCHEME2   = 'FDM-O4-CPML' ;
ACQUI     = 'acquisition.config' ;
X0        = 0.0 ;

% FIGURE 1: DG vs FDM - ABC
%-----------------------------
%COMPONENT = "vx" ;
%FIG_NAME  = 'HOMO-FREE-SURF-1' ;
%FILE1     = 'REF.django.config.py.fem.discontinuous.O5.elastic.O1.acqui_ABC.xml.snapshot' ;
%SCHEME1   = 'DGM-P5-ABC' ;

%FILE2     = 'REF.django.config.py.fem.discontinuous.O5.elastic.O1.acqui_NO_ABC.xml.snapshot' ;
%SCHEME2   = 'DGM-P5-NO-ABC' ;
%ACQUI     = 'REF.acqui_ABC.config' ;
%X0        = 0.0 ;

%-----------------

NX       = 201 ;
NZ       = 201 ;
NSNAP    = 4 ;
DX       = 10.0 ;
DZ       = 10.0 ;
XOR      = 0. ;
ZOR      = 0. ;
PERC     = 0.4 ;
ISNAP    = 1 ;

X = (0:NX-1) * DX ;
Z = (0:NZ-1) * DZ ;

%==========================================================================
% 1st file
%==========================================================================

FILE     = FILE1 + "." + COMPONENT ;
TITLE    = SCHEME1 ;

% read wavefield
file1 = fopen(FILE, 'r', 'ieee-le');
for i=1:ISNAP
    val1 = fread(file1, [NX NZ], FLOAT_SIZE);
end
fclose(file1);

subplot(2,2,2)
%axes('XAxisLocation', 'top', 'FontSize', FONT_SIZE)
hold on
xlabel('x (m)', 'FontSize', FONT_SIZE)
ylabel('z (m)', 'FontSize', FONT_SIZE)
title(TITLE, 'FontSize', FONT_SIZE+1)

%colorbar('FontSize', FONT_SIZE, 'location','southoutside')  ;
colormap(gray) ; box on ; %colorbar ;

imagesc(X, Z, val1') ;

max_axis = max(max(abs(val1))) * PERC ;
caxis([-max_axis max_axis])

axis equal; axis tight; axis ij

%text(0, -30, 0, '(b)', 'FontSize', FONT_SIZE+6)

% plot acquisition
file1 = fopen(ACQUI, 'r') ;
acqui_type = fscanf(file1, '%g', 1) 
acqui_dim = fscanf(file1, '%g', 1) 
acqui_nsrc = fscanf(file1, '%g', 1)
for ii = 1:acqui_nsrc
    id = fscanf(file1, '%g', 1) ;
    zsrc  = fscanf(file1, '%g', 1) ;    
    xsrc  = fscanf(file1, '%g', 1) - X0 ;
    plot(xsrc,zsrc,'kp','MarkerFaceColor','w','MarkerSize',15)
end
acqui_nrec = fscanf(file1, '%g', 1)
for ii = 1:acqui_nrec
    id = fscanf(file1, '%g', 1) ;
    zrec  = fscanf(file1, '%g', 1) ;    
    xrec  = fscanf(file1, '%g', 1) - X0  ;
    plot(xrec,zrec,'wv','MarkerFaceColor','k','MarkerSize',7)
end

%==========================================================================
% 2nd file
%==========================================================================

FILE     = FILE2 + "." + COMPONENT ;
TITLE    = SCHEME2 ;

% read wavefield
file1 = fopen(FILE, 'r', 'ieee-le');
for i=1:ISNAP
    val2 = fread(file1, [NX NZ], FLOAT_SIZE);
end
fclose(file1);

subplot(2,2,3)
%axes('XAxisLocation', 'top', 'FontSize', FONT_SIZE)
hold on
xlabel('x (m)', 'FontSize', FONT_SIZE)
ylabel('z (m)', 'FontSize', FONT_SIZE)
title(TITLE, 'FontSize', FONT_SIZE+1)

%colorbar('FontSize', FONT_SIZE, 'location','southoutside')  ;
colormap(gray) ; box on ; %colorbar ;

imagesc(X, Z, val2') ;

%max_axis = max(max(abs(val))) * PERC ;
caxis([-max_axis max_axis])

axis equal; axis tight; axis ij
%text(0, -30, 0, '(c)', 'FontSize', FONT_SIZE+6)

% plot acquisition
file1 = fopen(ACQUI, 'r') ;
acqui_type = fscanf(file1, '%g', 1) 
acqui_dim = fscanf(file1, '%g', 1) 
acqui_nsrc = fscanf(file1, '%g', 1)
for ii = 1:acqui_nsrc
    id = fscanf(file1, '%g', 1) ;
    zsrc  = fscanf(file1, '%g', 1) ;    
    xsrc  = fscanf(file1, '%g', 1) - X0 ;
    plot(xsrc,zsrc,'kp','MarkerFaceColor','w','MarkerSize',15)
end
acqui_nrec = fscanf(file1, '%g', 1)
for ii = 1:acqui_nrec
    id = fscanf(file1, '%g', 1) ;
    zrec  = fscanf(file1, '%g', 1) ;    
    xrec  = fscanf(file1, '%g', 1) - X0  ;
    plot(xrec,zrec,'wv','MarkerFaceColor','k','MarkerSize',7)
end

%==========================================================================
% difference
%==========================================================================

TITLE    = 'Difference' ;

% compute diff
diff = val2 - val1 ;

% compute nrms
sum1=0.0;
sum2=0.0;
nval = numel(diff) ;
for ival=1:nval
    sum1=sum1+diff(ival)^2 ;
    sum2=sum2+val2(ival)^2 ;
end
nrms=sqrt(sum1/sum2)

subplot(2,2,4)
hold on
xlabel('x (m)', 'FontSize', FONT_SIZE)
ylabel('z (m)', 'FontSize', FONT_SIZE)
title(TITLE, 'FontSize', FONT_SIZE+1)

colormap(gray) ; box on ; 
imagesc(X, Z, diff') ;
caxis([-max_axis max_axis])

axis equal; axis tight; axis ij

error = sprintf('NRMSE %0.3g', nrms) ;
text(1000., 170. , 0, error, 'FontSize', FONT_SIZE+1,...
    'Color','w','HorizontalAlignment','center')
%text(0, -30, 0, '(d)', 'FontSize', FONT_SIZE+6)

% plot acquisition
file1 = fopen(ACQUI, 'r') ;
acqui_type = fscanf(file1, '%g', 1) 
acqui_dim = fscanf(file1, '%g', 1) 
acqui_nsrc = fscanf(file1, '%g', 1)
for ii = 1:acqui_nsrc
    id = fscanf(file1, '%g', 1) ;
    zsrc  = fscanf(file1, '%g', 1) ;    
    xsrc  = fscanf(file1, '%g', 1) - X0 ;
    plot(xsrc,zsrc,'kp','MarkerFaceColor','w','MarkerSize',15)
end
acqui_nrec = fscanf(file1, '%g', 1)
for ii = 1:acqui_nrec
    id = fscanf(file1, '%g', 1) ;
    zrec  = fscanf(file1, '%g', 1) ;    
    xrec  = fscanf(file1, '%g', 1) - X0  ;
    plot(xrec,zrec,'wv','MarkerFaceColor','k','MarkerSize',7)
end

%==========================================================================
% model
%==========================================================================

NX       = 300 ;
NZ       = 300 ;
NSNAP    = 1 ;
DX       = 2000.0 / (NX-1) ;
DZ       = 2000.0 / (NZ-1) ;
XOR      = 0. ;
ZOR      = 0. ;
PERC     = 1.0 ;
ISNAP    = NSNAP ;

X = (0:NX-1) * DX ;
Z = (0:NZ-1) * DZ ;

FILE     = 'mesh.vp.django.out.bin' ;
TITLE    = 'Wave speed' ;

% read wavefield
file1 = fopen(FILE, 'r', 'ieee-le');
for i=1:ISNAP
    val2 = fread(file1, [NX NZ], FLOAT_SIZE);
end
fclose(file1);

% invert z axis
val1 = val2 ;
for iz=1:NZ
    val1(iz,:) = val2(NZ-iz+1,:) ;
end

ax=subplot(2,2,1) ;
ax.LineWidth = 2 ;
hold on
xlabel('x (m)', 'FontSize', FONT_SIZE)
ylabel('z (m)', 'FontSize', FONT_SIZE)
title(TITLE, 'FontSize', FONT_SIZE+1)
colormap(gray) ; box on ;
colorbar('south') ;

max_val=max(max(val1))
min_val=min(min(val1))
if (max_val - min_val) < 10
    caxis([0.95*max_val 1.05*max_val])
end

imagesc(X, Z, val1) ;
axis equal; axis tight; axis ij
%text(0, -30, 0, '(a)', 'FontSize', FONT_SIZE+6)

filename = sprintf('%s-%s-%s-%s.tmp.jpg', FIG_NAME, SCHEME1, SCHEME2, COMPONENT) ;
print('-djpeg', filename)