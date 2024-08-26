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
FONT_SIZE = 13 ;

% ANALYTICAL SOLUTION
PLOT_ANALYTIC = 1 ;
NT_SOL   = 500 ;
DT_SOL   = .004 ;
REP_SOL  = './' ;
%FILE_SOL_F = 'freceiver' ;

% NUMERICAL SOLUTION 1
PLOT_NUM1 = 1 ;
NT_NUM    = 1501 ;
DT_NUM    = 0.001 ;
REP_NUM   = './' ;
%FILE_NUM_F = 'P.freq' ;

% NUMERICAL SOLUTION 2
PLOT_NUM2 = 1 ;
NT_NUM2   = 1501 ;
DT_NUM2   = 0.001 ;
REP_NUM2  = './specfem2d/' ;
%REP_NUM2  = './' ;
%FILE_NUM_F = 'P.freq' ;

% VARIOUS
FIG_NAME  = 'MOD-HOMO-EL-INF-2D' ;
NT_PLOT   = 1501 ;
TMAX      = 1.5 ;
NREC      = 6 ; % nb receivers
DELTA_REC = 1 ;
GAP       = 100 ;
RATIO_RES = 10. ;
NORMALIZE = 0 ;
DISPLAY_RMS = 1 ;
PLOT_RESIDUAL = 1 ;

%=============================================================
%
%                  T I M E     D O M A I N
%
%=============================================================


dt        = TMAX/(NT_PLOT-1) ;
time_plot = (0:NT_PLOT-1) * dt ;
tot_rms   = 0. ;

for ifile = 1:4
    
    switch(ifile)
                
        case 1            
            FILE_NUM  = 'django.config.py.fdm.staggered.O4.elastic.O1.pml.xml.rec.vx' ;  
            SCHEME1   = 'FDM-O4' ;           
            FILE_NUM2 = 'vx.time.rec.specfem2d.out' ;
            SCHEME2   = 'SEM-O4' ;  
            alpha2    = 1 ;            
            shift2    = 0.0 ;
            COMP      = 'Vx' ;
            RATIO     = 1.e+7 ;           
            
        case 2            
            FILE_NUM  = 'django.config.py.fdm.staggered.O4.elastic.O1.pml.xml.rec.vz' ;  
            SCHEME1   = 'FDM-O4' ;            
            FILE_NUM2 = 'vz.time.rec.specfem2d.out' ;
            SCHEME2   = 'SEM-O4' ;   
            alpha2    = -1 ;            
            shift2    = 0.0 ;
            COMP      = 'Vz' ;
            RATIO     = 1.e+7 ;    
            
        case 3                       
            FILE_NUM  = 'django.config.py.fem.discontinuous.O6.elastic.O1.sponge.xml.rec.vx' ;
            SCHEME1   = 'DGM-P6' ;
            FILE_NUM2 = 'vx.time.rec.specfem2d.out' ;
            SCHEME2   = 'SEM-O4' ;  
            alpha2    = 1 ;            
            shift2    = 0.0 ;
            COMP      = 'Vx' ;
            RATIO     = 1.e+7 ;           
            
        case 4                       
            FILE_NUM  = 'django.config.py.fem.discontinuous.O6.elastic.O1.sponge.xml.rec.vz' ;
            SCHEME1   = 'DGM-P6' ;
            FILE_NUM2 = 'vz.time.rec.specfem2d.out' ;
            SCHEME2   = 'SEM-O4' ;   
            alpha2    = -1 ;            
            shift2    = 0.0 ;
            COMP      = 'Vz' ;
            RATIO     = 1.e+7 ;    
    end
    
    figure
    axes('FontSize', FONT_SIZE, 'LineWidth', 2)
    hold on
    
    if (NORMALIZE)
        sub_title = sprintf('%s - %s', COMP, 'Amplitude normalized') ;
    else
        sub_title = sprintf('%s - %s', COMP, 'True amplitude') ;
    end
    title({FIG_NAME ; sub_title}, 'FontSize', FONT_SIZE+2, 'Color', 'k')
    xlabel('Time (s)')
    ylabel('Receiver')
    
    %=============================================================
    %                    NUMERICAL SOLUTION 1
    %=============================================================
    
    if (PLOT_NUM1)
    
        NTMAX = NT_NUM;
        dt2   = DT_NUM;
        
        val_num  = zeros(NT_PLOT, NREC);
        
        % read seismograms
        filename = sprintf('%s%s', REP_NUM, FILE_NUM);
        fprintf('Read numerical solution 1 %s\n', filename)
        file1 = fopen(filename, 'r', 'ieee-le');
        val = fread(file1, [NREC,NTMAX], 'float32');
        fclose(file1);
        
        % interpolate seismograms
        for irec = 1:NREC
            
            val1 = 0.5 * val(irec, :);
            
            for ii = 1:NT_PLOT
                t = (ii - 1) * dt;
                n = floor(1 + (t / dt2));
                
                if ((n > 0) && (n < NTMAX-1))
                    
                    % linear interpolation between n and n+1
                    val_num(ii, irec) = val1(n) + ((t - (n-1)*dt2)/dt2 * (val1(n+1) - val1(n)));
                end
            end
            
        end
        
        % normalize
        max_val = max(max(abs(val_num))) ;
        if (NORMALIZE)            
            val_num = val_num ./ max_val ;
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
            
            val1 = 0.5 * val(:, irec)';
            
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
        
        max_val2 = max(max(abs(val_num2))) ;
        % normalize
        val_num2 = alpha2 * val_num2 ;
        if (NORMALIZE)            
            val_num2 = val_num2 ./ max_val2 ;
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
    
    if (0)
        
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
        
        % plot seismograms
        for irec = 1:DELTA_REC:NREC
            
            plot(time_plot(1:NT_PLOT), (RATIO * val_sol(1:NT_PLOT, irec)) + irec, '--k','LineWidth', 2)
            
        end
        
    end
    
    %=============================================================
    %                    RESIDUALS
    %=============================================================
    
    if (1)
        
        % compute residuals
        val_res = val_num - val_num2 ;
        
        % compute RMS
        sum_val_sol = 0. ;
        rms         = 0. ;
        
        for irec = 1:DELTA_REC:NREC            
            for ii = 1:NT_PLOT 
                rms = rms + val_res(ii, irec)^2.; 
                sum_val_sol = sum_val_sol + val_num2(ii, irec)^2. ;
            end
            
            % plot residuals
            if (PLOT_RESIDUAL)
                if (irec == 1)
                    p3 = plot(time_plot(1:NT_PLOT), (RATIO * RATIO_RES * val_res(1:NT_PLOT, irec)) + irec, ...
                        '-r','LineWidth', 1);
                else
                    plot(time_plot(1:NT_PLOT), (RATIO * RATIO_RES * val_res(1:NT_PLOT, irec)) + irec, ...
                        '-r','LineWidth', 1);
                end
            end
        end
        
        if (DISPLAY_RMS)
            error = sprintf('%s / NRMS = %0.3g', date, sqrt(rms / sum_val_sol)) ;
            text(0.6, +0.5 , 0, error, 'FontSize', FONT_SIZE)
        end
    end
    
    %===============================================================
    
    % save figure as .eps file
    if (PLOT_RESIDUAL)
        label3 = sprintf('Residual X%d', RATIO_RES) ;
        legend([p1 p2 p3],SCHEME1,SCHEME2, label3)
    else
        legend([p1 p2],SCHEME1,SCHEME2)
    end
    axis([0 max(time_plot) 0 NREC+1])
    
    figure_name = sprintf('%s-TIME-%s-%s.tmp.jpg', FIG_NAME, SCHEME1, COMP) ;
    print('-djpeg', figure_name)
    
    fprintf('Max val1 %d\n', max_val) ;
    fprintf('Max val2 %d\n', max_val2) ;
    fprintf('* ratio Max val1/val2 %d\n', max_val/max_val2) ;
    fprintf('* ratio Max val2/val1 %d\n', max_val2/max_val) ;
    
end % for ifile=1:4

return;