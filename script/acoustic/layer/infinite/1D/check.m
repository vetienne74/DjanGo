%======================================================================
%
% This program compare analytical and numerical solution in time
% and frequency domain
%
%-----------------------------------------------------------------------
%
% INPUT  : analytical and computed data
%
% OUTPUT : figures
%
%======================================================================

clear all ; close all

fprintf('---------------------------------------\n')
fprintf('   Plot seismograms \n')
fprintf('---------------------------------------\n\n')

%======================================================================
%
%                           INPUT PARAMETERS
%
%======================================================================

% DISPLAY PARAMETERS
NORMALIZE = 0 ;
FONT_SIZE = 13 ;

% ANALYTICAL SOLUTION
PLOT_ANALYTIC = 1 ;
NT_SOL   = 801 ;
DT_SOL   = .005 ;
REP_SOL  = './' ;
%FILE_SOL_F = 'freceiver' ;

% NUMERICAL SOLUTION 1
PLOT_NUM1 = 1 ;
NT_NUM    = 801 ;
DT_NUM    = 0.005 ;
REP_NUM   = './' ;
%FILE_NUM_F = 'P.freq' ;

% NUMERICAL SOLUTION 2
PLOT_NUM2 = 0 ;
NT_NUM2   = 801 ;
DT_NUM2   = 0.005 ;
REP_NUM2  = './' ;
%FILE_NUM_F = 'P.freq' ;

% VARIOUS
FIG_NAME      = 'MOD-HOMO-AC-INF-1D' ;
NT_PLOT       = 801 ;
TMAX          = 4.0 ;
NREC          = 11 ; % nb receivers
DELTA_REC     = 1 ;
GAP           = 100 ;
NORMALIZE     = 0 ;
DISPLAY_RMS   = 1 ;
PLOT_RESIDUAL = 1 ;

%=============================================================
%
%                  T I M E     D O M A I N
%
%=============================================================


dt        = TMAX/(NT_PLOT-1) ;
time_plot = (0:NT_PLOT-1) * dt ;
tot_rms   = 0. ;

for icomp = 3:3
    
    switch(icomp)
        
        case 1
            %FILE_SOL = 'VX_sol' ;
            FILE_NUM  = 'vx.time.rec.django.bin.out' ;
            FILE_NUM2 = 'vx.time.rec.specfem2d.out' ;
            alpha2    = 1 ;
            %shift2    = -0.002 ;
            shift2    = 0.0 ;
            TIT_FIG   = 'Vx component' ;
            %RATIO     = 1.e+7 ;
            
        case 2
            %FILE_SOL = 'VZ_sol' ;
            FILE_NUM  = 'vz.time.rec.django.bin.out' ;
            FILE_NUM2 = 'vz.time.rec.specfem2d.out' ;
            alpha2    = -1 ;
            %shift2    = -0.002 ;
            shift2    = 0.0 ;
            TIT_FIG  = 'Vz component' ;
            %RATIO    = 1.e+7 ;
            
        case 3
            FILE_SOL  = 'pr.time.rec.analytic_homo.out' ;            
            %FILE_NUM = 'django.config.py.femcontinuous.1st.O2.h10.no_pml.src_file1.xml.pr' ;
            %FILE_NUM = 'django.config.py.femcontinuous.1st.O5.h10.no_pml.src_file1.xml.pr' ;
            %FILE_NUM = 'django.config.py.fdm.2nd.O8.h10.no_pml.src_func.xml.pr' ;
            FILE_NUM = 'django.config.py.fdm.2nd.O2.magic_dt.no_pml.src_func.xml.pr' ; 
            FILE_NUM = 'django.config.py.fdm.1st.O2.magic_dt.no_pml.src_func.xml.pr' ; 
            RATIO_RES     = 10 ;
            %FILE_NUM2 = 'XXX' ;
            alpha2    = -1 ;
            shift2    = 0.0 ;            
            %shift2    = 0.005 ;
            TIT_FIG  = 'P component' ; 
            RATIO    = 1.e+4 ;           
            
        case 4
            FILE_SOL = 'P_sol' ;
            FILE_NUM = 'pr.time.rec.cornflex.out' ;
            TIT_FIG  = 'Seismograms (pressure component)' ;
            RATIO    = 1. ;
            
    end
    
    figure
    axes('FontSize', FONT_SIZE, 'LineWidth', 2)
    hold on
    
    FILE_NUM_TMP=replace(FILE_NUM,'.xml.pr','') ;
    FILE_NUM_TITLE=replace(FILE_NUM_TMP,'django.config.','') ;
    if (NORMALIZE == 1)
        title({FIG_NAME ; FILE_NUM_TITLE; 'Normalized amplitude'},'Interp','None','FontSize', FONT_SIZE, 'Color', 'k')
    else
        title({FIG_NAME ; FILE_NUM_TITLE; 'True amplitude'},'Interp','None', 'FontSize', FONT_SIZE, 'Color', 'k')
    end
    xlabel('Time (s)')
    ylabel('Receiver')
    
    %=============================================================
    %                    NUMERICAL SOLUTION 1
    %=============================================================
    
    if (PLOT_NUM1)
    
        NTMAX = NT_NUM;
        dt2   = DT_NUM;
        shift = shift2 / DT_NUM2 
        
        val_num  = zeros(NT_PLOT, NREC);
        val_numtmp  = zeros(NT_PLOT, NREC);
        
        % read seismograms
        filename = sprintf('%s%s', REP_NUM, FILE_NUM);
        fprintf('Read numerical solution 1 %s\n', filename)
        file1 = fopen(filename, 'r', 'ieee-le');
        val = fread(file1, [NREC,NTMAX], 'float32');
        fclose(file1);
        
        % interpolate seismograms
        for irec = 1:NREC
            
            val1 = val(irec, :);
            
            for ii = 1:NT_PLOT
                t = (ii - 1) * dt;
                n = floor(1 + (t / dt2));
                
                if ((n > 0) && (n < NTMAX-1))
                    
                    % linear interpolation between n and n+1
                    val_numtmp(ii, irec) = val1(n) + ((t - (n-1)*dt2)/dt2 * (val1(n+1) - val1(n)));
                end
            end
            
        end
        
        % shift
        if (shift == 0)
            val_num = val_numtmp ;
            
        else
            for irec = 1:NREC
                for ii = 1:NT_PLOT
                    iii = ii - shift ;
                    if (iii >= 1) && (iii <= NT_PLOT)
                        val_num(ii, irec) = val_numtmp(iii, irec) ;
                    end
                end
            end
        end
        
        % normalize
        max_val_num = max(max(abs(val_num))) ;
        if (NORMALIZE)            
            val_num = val_num ./ max_val_num ;
            RATIO = 1 ;
        end
        
        % plot seismograms
        for irec = 1:DELTA_REC:NREC
            if (irec == 1)
                p1 = plot(time_plot(1:NT_PLOT), (RATIO * val_num(1:NT_PLOT, irec)) + irec, ...
                    '-k','LineWidth', 1) ;
            else
                plot(time_plot(1:NT_PLOT), (RATIO * val_num(1:NT_PLOT, irec)) + irec, ...
                    '-k','LineWidth', 1) ;
            end
        end
        
    end
    
    %=============================================================
    %                    NUMERICAL SOLUTION 2
    %=============================================================
    
    if (PLOT_NUM2)
    
        NTMAX  = NT_NUM2;
        dt2    = DT_NUM2;
        shift  = shift2 / DT_NUM2 
        
        val_num2     = zeros(NT_PLOT, NREC);
        val_num2tmp  = zeros(NT_PLOT, NREC);
        
        % read seismograms
        filename = sprintf('%s%s', REP_NUM2, FILE_NUM2);
        fprintf('Read numerical solution 2 %s\n', filename)
        file1 = fopen(filename, 'r', 'ieee-le');
        %val = fread(file1, [NREC,NTMAX], 'float32');
        val = fread(file1, [NTMAX,NREC], 'float32');
        fclose(file1);
        
        % interpolate seismograms
        for irec = 1:NREC
            
            val1 = val(:, irec)';
            
            for ii = 1:NT_PLOT
                t = (ii - 1) * dt;
                n = floor(1 + (t / dt2));
                
                if ((n > 0) && (n < NTMAX-1))
                    
                    % linear interpolation between n and n+1
                    val_num2tmp(ii, irec) = val1(n) + ((t - (n-1)*dt2)/dt2 * (val1(n+1) - val1(n)));
                end
            end
            
        end
        
        % shift
        if (shift == 0)
            val_num2 = val_num2tmp ;
            
        else
            for irec = 1:NREC
                for ii = 1:NT_PLOT
                    iii = ii - shift ;
                    if (iii >= 1) && (iii <= NT_PLOT)
                        val_num2(ii, irec) = val_num2tmp(iii, irec) ;
                    end
                end
            end
        end
        
        % normalize
        max_val_num2 = max(max(abs(val_num2))) ;
        if (NORMALIZE)            
            val_num2 = alpha2 * val_num2 ./ max_val_num2 ;
            RATIO = 1 ;
        end
        
        % plot seismograms
        for irec = 1:DELTA_REC:NREC
            if (irec == 1)
                p2 = plot(time_plot(1:NT_PLOT), (RATIO * val_num2(1:NT_PLOT, irec)) + irec, ...
                    '--k','LineWidth', 2) ;
            else
                plot(time_plot(1:NT_PLOT), (RATIO * val_num2(1:NT_PLOT, irec)) + irec, ...
                    '--k','LineWidth', 2) ;
            end
        end
        
    end
    
    %=============================================================
    %                    ANALYTICAL solution
    %=============================================================
    
    if (PLOT_ANALYTIC)
        
        NTMAX = NT_SOL;
        dt2   = DT_SOL;
        
        val_sol  = zeros(NT_PLOT, NREC);
        
        % read seismograms
        filename = sprintf('%s%s', REP_SOL, FILE_SOL);
        fprintf('Read analytical solution %s\n', filename)
        file1 = fopen(filename, 'r', 'ieee-le');
        valtemp = fread(file1, [NTMAX,NREC], 'float32');
        val = transpose(valtemp) ;
        %val = valtemp ;
        status = fclose(file1);
        
        % interpolate seismograms
        for irec = 1:NREC
            
            val1 = val(irec, :);
            
            for ii = 1:NT_PLOT
                t = (ii - 1) * dt;
                n = floor(1 + (t / dt2));
                
                if ((n > 0) && (n < NTMAX-1))
                    
                    % linear interpolation between n and n+1
                    val_sol(ii, irec) = val1(n) + ((t - (n-1)*dt2)/dt2 * (val1(n+1) - val1(n)));
                end
            end
            
        end
        
        % normalize
        max_val_sol = max(max(abs(val_sol))) ;
        if (NORMALIZE)            
            val_sol = val_sol ./ max_val_sol ;
            RATIO = 1 ;
        end
        
        % plot seismograms
        for irec = 1:DELTA_REC:NREC
            if (irec == 1)
                p3 = plot(time_plot(1:NT_PLOT), (RATIO * val_sol(1:NT_PLOT, irec)) + irec, '--k','LineWidth', 2) ;
            else
                plot(time_plot(1:NT_PLOT), (RATIO * val_sol(1:NT_PLOT, irec)) + irec, '--k','LineWidth', 2) ;
            end
        end
        
    end
    
    %=============================================================
    %                    RESIDUALS
    %=============================================================
    
    if (1)
        
        % compute residuals
        val1 = val_num ; 
        val2 = val_sol ;        
        val_res = val1 - val2 ;
        
        % compute RMS
        sum_val_sol = 0. ;
        rms         = 0. ;
        
        for irec = 1:DELTA_REC:NREC     
            
            %if (irec ~= 6)
            for ii = 1:NT_PLOT 
                rms = rms + val_res(ii, irec)^2.; 
                sum_val_sol = sum_val_sol + val2(ii, irec)^2. ;
            end
            %end
            
            % plot residuals                        
            if (PLOT_RESIDUAL)
                if (irec == 1)
                    p4 = plot(time_plot(1:NT_PLOT), (RATIO * RATIO_RES * val_res(1:NT_PLOT, irec)) + irec, '-r','LineWidth', 1) ;
                else
                    plot(time_plot(1:NT_PLOT), (RATIO * RATIO_RES * val_res(1:NT_PLOT, irec)) + irec, '-r','LineWidth', 1) ;
                end
            end
        end
        
        if (DISPLAY_RMS)
            % initial RMS formula
            %nrms = 100*sqrt(rms / sum_val_sol) ;
            % alternative RMS formula
            %nrms = 100*sqrt(rms) / numel(val2) / (max(max(val2)) - min(min(val2))) ;
            %error = sprintf('%s / NRMS = %0.3g %%', date, nrms) ;            
            %fprintf('### %s / NRMS = %0.3g %% \n', date, nrms)
            
            nrms = sqrt(rms / sum_val_sol) ;
            error = sprintf('%s / NRMS = %0.3g', date, nrms) ;            
            fprintf('### %s / NRMS = %0.3g\n', date, nrms)
            
            text(0.6, +0.5 , 0, error, 'FontSize', FONT_SIZE)
        end
    end
    
    %===============================================================
    
    if (PLOT_RESIDUAL)
        diff_legend = sprintf("diff x%6.0f", RATIO_RES) ;
        legend([p1 p3 p4], 'django', 'analytical', diff_legend, 'Location', 'northwest')
    else
        legend([p1 p3], 'django', 'analytical', 'Location', 'northwest')
    end
    
    % save figure as .eps file
    %legend([p1 p2],'FDM O4','SEM O4')
    axis([0 max(time_plot) 0 NREC+1])
    %figure_name = sprintf('%s-TIME-%s.eps', FIG_NAME, FILE_NUM) ;
    %print('-depsc', figure_name)
    figure_name = sprintf('%s-TIME-%s.jpg', FIG_NAME, FILE_NUM_TITLE) ;
    print('-djpeg', figure_name)
    
    ratio_sol_num = max_val_sol / max_val_num
    ratio_num_sol = 1 / ratio_sol_num
    
end % for icomp=1:4

return



%=============================================================
%
%              F R E Q U E N C Y      D O M A I N
%
%=============================================================

figure
axes('FontSize', FONT_SIZE, 'LineWidth', 2)
hold on

title('Frequency solution (Pressure component)', 'FontSize', FONT_SIZE+2, 'Color', 'k')
xlabel('Receiver')
ylabel('Amplitude')

%=============================================================
%                    NUMERICAL solution
%=============================================================

if (1)
    
    % read solution
    filename = sprintf('%s%s', REP_NUM, FILE_NUM_F);
    fprintf('Read numerical solution %s\n', filename)
    file1 = fopen(filename, 'r', 'ieee-le');
    val_num = fread(file1, [2,NREC], 'float32');
    status = fclose(file1);
    
    max_amp = max(max(abs(val_num))) ;
    
    % plot real part
    plot(1:NREC, val_num(1, 1:NREC), '-k','LineWidth', 1)
    
    % plot imaginary part
    plot(1:NREC, val_num(2, 1:NREC) + 2*max_amp, '-r','LineWidth', 1)
    
end

%=============================================================
%                    ANALYTICAL solution
%=============================================================

if (1)
    
    % read solution
    filename = sprintf('%s%s', REP_SOL, FILE_SOL_F);
    fprintf('Read numerical solution %s\n', filename)
    file1 = fopen(filename, 'r', 'ieee-le');
    val_sol = fread(file1, [2,NREC], 'float32');
    status = fclose(file1);
    
    % normalize and invert imaginary part (due to opposite convetion of
    % the subroutine realft from Numerical recipes
    val_sol(1, :) = val_sol(1, :) / NT_SOL ;
    val_sol(2, :) = -val_sol(2, :) / NT_SOL ;
    
    % plot real part
    plot(1:NREC, val_sol(1, 1:NREC), '--k','LineWidth', 2)
    
    % plot imaginary part
    plot(1:NREC, val_sol(2, 1:NREC) + 2*max_amp, '--r','LineWidth', 2)
    
end

%=============================================================
%                    RESIDUALS
%=============================================================

if (1)
    
    % compute residuals
    val_res = val_num - val_sol ;
    
    tot_rms   = 0. ;
    
    % rms for real part
    rms = 0. ;
    for irec = 1:NREC; rms = rms + val_res(1, irec)^2.; end
    sum_val_sol = 0. ;
    for irec=1:NREC; sum_val_sol = sum_val_sol + val_sol(1, irec)^2. ; end
    rms = sqrt(rms / sum_val_sol) ;
    tot_rms = tot_rms + rms ;
    
    % rms for imaginary part
    rms = 0. ;
    for irec = 1:NREC; rms = rms + val_res(2, irec)^2.; end
    sum_val_sol = 0. ;
    for irec=1:NREC; sum_val_sol = sum_val_sol + val_sol(2, irec)^2. ; end
    rms = sqrt(rms / sum_val_sol) ;
    tot_rms = tot_rms + rms ;
    
    % plot real part
    plot(1:NREC, val_res(1,1:NREC), 'Color',[0.5,0.5,0.5],'LineWidth', 1)
    
    % plot imaginary part
    plot(1:NREC, val_res(2,1:NREC) + 2*max_amp, 'Color',[0.5,0.5,0.5],'LineWidth', 1)
    
    error = sprintf('%s / RMS = %0.3g', date, tot_rms);
    text(NREC/2, max_amp, 0, error, 'FontSize', FONT_SIZE)
    
end

%===============================================================

% save figure as .eps file
figure_name = sprintf('%s-FREQ.eps', FIG_NAME) ;
print('-depsc', figure_name)

return;