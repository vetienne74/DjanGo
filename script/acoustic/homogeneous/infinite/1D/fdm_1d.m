
close all
clear all

NZ = 1001 ;
DZ = 20 ;
NT = 401 ;
DT = 0.005 ;
VP = 4000 ;

PERC = 0.5 ;

zsrc = 10000. ;
zrec_min = 0 ;
zrec_max = 20000 ;
nrec = (zrec_max-zrec_min) / DZ + 1

zz = (0:NZ-1) * DZ ;
tt = (0:NT-1) * DT ;

% read src wavelet
fwavelet = fopen('src.wavelet.django.out.bin', 'r', 'ieee-le') ;
wavelet = fread(fwavelet, NT, 'float32') ;
fclose(fwavelet) ;
figure
plot(wavelet)

% source
p = zeros(nrec, NT) ;
izsrc = zsrc / DZ + 1
izrec_min = zrec_min / DZ + 1 
izrec_max = zrec_max / DZ + 1 
for it = 1:NT
    %for iz = 1:NZ
        p(izsrc, it) = wavelet(it) ;
    %end
end

figure
hold on
title('input wavelet')
xlabel('z (m)')
ylabel('time (s)')
imagesc(zz, tt, p') 
axis ij
axis tight
colorbar 
colormap(gray)
max_val = PERC*max(max(abs(p))) ;
caxis([-max_val max_val])

% compute wavefield
p = zeros(nrec, NT) ;
izsrc = zsrc / DZ + 1
izrec_min = zrec_min / DZ + 1 
izrec_max = zrec_max / DZ + 1 
for it = 2:NT-1
    for iz = 2:NZ-1
        z_lap = (p(iz-1, it) - 2*p(iz, it) + p(iz+1, it)) / DZ^2;
        p(iz, it+1) = 2*p(iz, it) -p(iz, it-1) + (DT*VP)^2 * z_lap ;
    end
    p(izsrc, it+1) = p(izsrc, it+1) + DT*DT*wavelet(it) ;
end

figure
hold on
title('pressure wavefield')
xlabel('z (m)')
ylabel('time (s)')
imagesc(zz, tt, p') 
axis ij
axis tight
colorbar 
colormap(gray)
max_val = PERC*max(max(abs(p))) ;
caxis([-max_val max_val])

% write wavefield
fwave = fopen('pr.time.rec.matlab.out', 'w', 'ieee-le') ;
fwrite(fwave, p(izrec_min:izrec_max,:),'float32')
fclose(fwave) ;
