

% build file with plucking force
% reference: Time domain simulation of a guitar: model and method
% Derveaux, Chaigne, Joly and Becache
% J. Acoust. Soc. Am. 2003

clear all
close all

% music scores
tmin    = 0.0 ;
tmax    = 167 ;
noire   = 0.3 ;
PLUCK   = 1 ;
HAMMER  = 2 ;
PULLOFF = 3 ;

%--------------------------------------------------------------------------
% STRING SETTINGS
%--------------------------------------------------------------------------
% string length (m)
length = 0.65 ;

% tunning E2 A2 D3 G3 B4 E4
% string 1 = high E
% string 6 = low E
A4 =440 ;
A3 = A4 / 2 ;
A2 = A3 / 2 ;

% std tunning 
A = A2 ;

freq(6) = A / 2^((5)/12) ; % E2
freq(5) = A ; % A2
freq(4) = A * 2^((5)/12) ; % D3
freq(3) = A * 2^((10)/12) ; % G3
freq(2) = A * 2^((12+2)/12) ; % B4
freq(1) = A * 2^((12+7)/12) ; % E4

disp(freq)

%--------------------------------------------------------------------------
% prelude C major
% Jean Sebastien Bach
%--------------------------------------------------------------------------

if (1)    
    prelude    
end

% rando sequence
if(0)
    n_note = 100 ;
    for inote = 1:n_note
        note(inote).string = randi(6) ;
        note(inote).time   = rand * tmax ;
    end
end

% define the pulse
t1   = 0.015 ;
t2   = 0.0154 ;
dt   = 0.0001 ;

nt = (tmax-tmin) / dt + 1 ;
hh = zeros(nt,1) ;
tt = zeros(nt,1) ;

for it = 1:nt
    tt(it) = (it-1) * dt ;
    if (tt(it) >= 0) && (tt(it) <= t1)
        hh(it) = 1-cos(pi*tt(it)/t1) ;
    elseif (tt(it)>t1) && (tt(it)<=t2)
        hh(it) = 1 + cos(pi*(tt(it)-t1)/t2) ;
    else
        hh(it) = 0.0 ;
    end
end

% plot figure
figure
plot(tt, hh)
xlabel('time (s)')
ylabel('amplitude')
title('Plucking force')
axis([0 0.05 0 2])

%--------------------------------------------------------------------------
% RIGHT HAND
%--------------------------------------------------------------------------

% build the source functions
% loop on the notes
force_func = zeros(6,nt) ;
for inote = 1:n_note
    if (note(inote).force == PLUCK)
        istring = note(inote).string ;
        delay   = round(note(inote).time / dt + 1) ;
        amp = 0.5 + (6 - istring) * 0.2 ;
        force_func(istring, delay+1:nt) = force_func(istring, delay+1:nt) + amp * hh(1:nt-delay)' ;
    end
end

% plot source functions
figure
xlabel('time (s)')
ylabel('amplitude')
title('Plucking force')
hold on
for istring = 1:6
    plot(tt, istring + force_func(istring,:))
end
axis ij

% write source functions
for istring = 1:6
    file_name = sprintf('force_string%d.out.bin', istring) ;
    f1=fopen(file_name, 'w') ;
    fwrite(f1, force_func(istring,:), 'float') ;
    fclose(f1) ;
end

%--------------------------------------------------------------------------
% LEFT HAND
%--------------------------------------------------------------------------

% build the left hand file
% loop on the notes
fret_left_hand = zeros(6,nt) ;
for inote = 1:n_note
    istring  = note(inote).string ;
    ifret    = note(inote).fret ;
    delay    = round(note(inote).time / dt + 1) ;
    %delay2   = min(round(noire / dt + 1), nt) ;
    %fret_left_hand(istring, delay+1:delay2) = ifret ; 
    fret_left_hand(istring, delay+1:nt) = ifret ; 
end

% plot fret 
figure
xlabel('time (s)')
ylabel('amplitude')
title('Left hand (fret)')
hold on
for istring = 1:6
    plot(tt, istring + fret_left_hand(istring,:))
end
axis ij

% convert fret to finger position
position_left_hand = zeros(6,nt) ;
for istring = 1:6
    for it = 1:nt
        position_left_hand(istring, it) = length * (1 - 2^(-fret_left_hand(istring, it)/12)) ;
    end
end

% plot positions
figure
xlabel('time (s)')
ylabel('amplitude')
title('Left hand (position)')
hold on
for istring = 1:6
    plot(tt, istring + position_left_hand(istring,:))
end
axis ij

% write source
for istring = 1:6
    file_name = sprintf('left_hand_string%d.out.bin', istring) ;
    f1=fopen(file_name, 'w') ;
    fwrite(f1, position_left_hand(istring,:), 'float') ;
    fclose(f1) ;
end

max_time = 0 ;
for inote =1:n_note
    max_time = max(note(inote).time, max_time) ;
end
max_time
