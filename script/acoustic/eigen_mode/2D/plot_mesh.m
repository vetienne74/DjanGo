
close all ;
clear all ;

%FILE_IN = "grid.point.django.out.ascii" ;
FILE_IN = "node.coord.django.out.ascii" ;

VAL = load(FILE_IN) ;
DIM = size(VAL) ;
NPOINT = DIM(1) ;

figure
hold on
if (DIM(2) == 1)
    % 1D MESH   
    for ip=1:NPOINT
        plot(VAL(ip,1), 1., '.k')
    end
    
elseif (DIM(2) == 2)
    % 2D MESH
    for ip=1:NPOINT
        plot(VAL(ip,2), VAL(ip,1), '.k')
    end   
end