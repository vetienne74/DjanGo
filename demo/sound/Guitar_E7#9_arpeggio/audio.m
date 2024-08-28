
clear all ;
close all ;

% $$$ SELECT MODE $$$
MODE = 2 ; 
% 1 = write sound, 2 = play sound

if (MODE == 1)
    
    track = './speak_to_me_breath.mp3' ;
    
    [y Fs] = audioread(track) ;
    audioinfo(track)
    nt2 = size(y(:,:))
    
    tmin = 70. ;
    tmax = 83. ;
    
    imin = tmin * Fs
    imax = tmax * Fs
    nt = imax - imin + 1
    
    yy = y(imin:imax,1) ;
    tt = (0:nt-1) / Fs ;
    tmax = max (tt) ;
    
    % taper
    ttaper = 1 ;
    ntaper = ttaper * Fs
    
    i1 = 1 ;
    i2 = i1 + ntaper - 1;
    yy(i1:i2) = yy(i1:i2) .* ((1:ntaper) / ntaper)' ;
    i1 = nt - ntaper + 1 ;
    i2 = nt ;
    yy(i1:i2) = yy(i1:i2) .* ((ntaper:-1:1) / ntaper)' ;
    
    % play sound
    player=audioplayer(yy,Fs)
    play(player)
    isplaying(player)
    %stop(player)
    
    % write sound
    f1=fopen('speak_to_me.bin', 'w') ; 
    fwrite(f1, yy, 'float')
    fclose(f1) ; 
    
end

if (MODE == 2)
    
    Fs = 44100 ;    
    
    % sum the 6 strings    
    file_name = 'django.config.run1.xml.pr.out.bin' ;   
    fprintf('read file %s\n', file_name) ;      
    f1=fopen(file_name, 'r') ;
    yy = fread(f1, 'float') ;
    fclose(f1) ;    
    file_name = 'django.config.run2.xml.pr.out.bin' ;   
    fprintf('read file %s\n', file_name) ;      
    f1=fopen(file_name, 'r') ;
    yy = yy + fread(f1, 'float') ;
    fclose(f1) ;
    file_name = 'django.config.run3.xml.pr.out.bin' ;   
    fprintf('read file %s\n', file_name) ;      
    f1=fopen(file_name, 'r') ;
    yy = yy + fread(f1, 'float') ;
    fclose(f1) ;
    file_name = 'django.config.run4.xml.pr.out.bin' ;   
    fprintf('read file %s\n', file_name) ;      
    f1=fopen(file_name, 'r') ;
    yy = yy + fread(f1, 'float') ;
    fclose(f1) ;
    file_name = 'django.config.run5.xml.pr.out.bin' ;   
    fprintf('read file %s\n', file_name) ;      
    f1=fopen(file_name, 'r') ;
    yy = yy + fread(f1, 'float') ;
    fclose(f1) ;
    file_name = 'django.config.run6.xml.pr.out.bin' ;   
    fprintf('read file %s\n', file_name) ;      
    f1=fopen(file_name, 'r') ;
    yy = yy + fread(f1, 'float') ;
    fclose(f1) ;
    
    nt=size(yy,1) ;
    tt=(0:nt-1)/Fs ;
    fprintf('# time steps %d \n', nt) ;
    tmax = max(tt) ;
    fprintf('# tmax       %f \n', tmax) ;
    
    %player=audioplayer(yy,Fs) ;
    %play(player)
    %isplaying(player) ;
    %soundsc(yy,Fs) ;
    
    yy = yy ./ max(abs(yy(:))) ;
    
    player=audioplayer(yy,Fs) ;
    play(player)
    
    new_file = sprintf('django.E7#9.out.wav') ;
    audiowrite(new_file, yy, Fs) ;
    
    figure
    plot(tt,yy)
    title(file_name) 
end

%return

% AMPLITUDE

figure
plot(tt, yy)

%ntfft   = 2^12 ;
ntfft   = 2^12 ;
ntslice = floor (nt / ntfft)
ntfft2  = ntfft/2+1 ;
ff = Fs*(0:(ntfft/2))/ntfft;
df = Fs /ntfft ;

% SPECTROGRAM

spectro = zeros (ntslice, ntfft2) ;

for islice = 1 : ntslice
    
    istart = (islice - 1) * ntfft +1 ;
    %istart = islice ;
    iend   = istart + ntfft ;
    
    Y = fft (yy(istart:iend)) ;
    P = abs(Y/ntfft);
    
    spectro(islice, :) = P(1:ntfft2) ;
end

figure
title('Spectrogram')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
ax = gca;
hold on
colormap(hot) ;

imagesc(tt, ff, spectro')
caxis([0. max(max(spectro))*0.2])

axis([0 tmax 0 2000])

return

% ENVELOP

max_val = zeros(ntslice, 1) ;
for it = 1:ntslice
    max_val(it) = max(spectro(it,:)) ;
end
figure
title('Max. amp.')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
hold on
plot(max_val)

% CHROMATIC SCALE

A = 440. / 2 / 2 / 2

note(1) = A ;
for inote = 2: (12 * 4)
    note(inote) = note(inote-1) * 2^(1./12.) ;
end

scale = round(note / df) + 1 ;

max_scale = size(scale) ;
size_scale = max_scale(2)

chromatogram = zeros (ntslice, size_scale) ;
for islice = 1 : ntslice
    for iscale = 1 : size_scale
        chromatogram(islice, iscale) = spectro(islice, scale(iscale)) ;
    end
end

figure
title('Chromatogram')
xlabel('Time (s)')
ylabel('Scale degree')
ax = gca;
hold on
colormap(hot) ;

imagesc(tt, scale, chromatogram')
caxis([0. 0.02])
axis tight
%axis([0 tmax 0 1500])