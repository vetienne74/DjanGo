
close all
clear all

% input parameters
NZ = 1001 ;
DZ = 20. ;
NT = 401 ;
DT = 0.005 ;
VP = 4000 ;

PERC = 0.5 ;

% open file
%filename = 'pr.time.rec.cornflex.out' ;
filename = 'pr.time.rec.matlab.out' ;
file1 = fopen(filename, 'r', 'ieee-le');
wave_in = fread(file1, [NZ, NT], 'float32');
status = fclose(file1);

zz = (0:NZ-1) * DZ ;
tt = (0:NT-1) * DT ;

% plot input data
figure
hold on
title('Input wavefield')
xlabel('z (m)')
ylabel('time (s)')
imagesc(zz, tt, wave_in') 
axis ij
axis tight
colorbar 
colormap(gray)
max_val = PERC*max(max(abs(wave_in))) ;
caxis([-max_val max_val])

% inject wavefield into equation
wave_out = zeros(NZ, NT) ;
for it = 2:NT-1
    for iz = 2:NZ-2
        z_lap = (wave_in(iz-1, it) - 2*wave_in(iz, it) + wave_in(iz+1, it)) / DZ^2;
        t_lap = (wave_in(iz, it-1) - 2*wave_in(iz, it) + wave_in(iz, it+1)) / DT^2;
        wave_out(iz, it) = t_lap - VP^2 * z_lap ;
        %wave_out(iz, it) = z_lap ;
    end
end

% plot output data
figure
hold on
title('Output wavefield')
xlabel('z (m)')
ylabel('time (s)')
imagesc(zz, tt, wave_out') 
axis ij
axis tight
colorbar
colormap(gray)
max_val = PERC*max(max(abs(wave_out))) ;
caxis([-max_val max_val])

% plot retrieved source
figure
hold on
wavelet_out = wave_out(501,:) ;
plot(tt, wavelet_out,'--k', 'LineWidth', 2.)

% plot real source
filename = 'src.wavelet.django.out.bin' ;
file1 = fopen(filename, 'r', 'ieee-le');
wavelet_in = fread(file1, NT, 'float32');
status = fclose(file1);
plot(tt, wavelet_in,'-k', 'LineWidth', 1.)

ratio = max(abs(wavelet_out)) / max(abs(wavelet_in))
plot(tt, wavelet_out.*-ratio,'-r', 'LineWidth', 1.)

legend('retrieved', 'source', '-retrieved') 