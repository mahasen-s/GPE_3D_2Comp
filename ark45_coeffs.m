function [A,B,C] = ark45_coeffs(method)
% method is one of 'Dormand-Prince','Cash-Karp', or 'Fehlberg'
% Returns A,B,C, where A,B,C are the Butcher tableau matrix,
% weights, and nodes, respectively.

% Define Butcher tableaus
ARK45_DP_nodes  = [1/5, 3/10, 4/5, 8/9, 1, 1];
ARK45_DP_matrix = [
    1/5         3/40    44/45   19372/6561      9017/3168       35/384
    0           9/40    -56/15  -25360/2187     -355/33         0
    0           0       32/9    64448/6561      46732/5247      500/1113
    0           0       0       -212/729        49/176          125/192
    0           0       0       0               -5103/18656     -2187/6784
    0           0       0       0               0               11/84
    ];
ARK45_DP_weights = [
                    0           0   0           0           0               0           1
                    5179/57600 	0 	7571/16695 	393/640 	-92097/339200 	187/2100 	1/40
                    ];

ARK45_CK_nodes  = [1/5, 3/10, 3/5, 1, 7/8];
ARK45_CK_matrix = [
    1/5         3/40    3/10    -11/54      1631/55296
    0           9/40    -9/10   5/2         175/512
    0           0       6/5     -70/27      575/13824
    0           0       0       35/27       44275/110592
    0           0       0       0           253/4096
    ];
ARK45_CK_weights = [
                    37/378 	0 	250/621 	125/594 	0 	512/1771
                    2825/27648 	0 	18575/48384 	13525/55296 	277/14336 	1/4
                    ];
                
ARK45_F_nodes   = [1/4, 3/8, 12/13, 1, 1/2];
ARK45_F_matrix  = [
    1/4         3/32    1932/2197    439/216        -8/27
    0           9/32    -7200/2197   -8             2
    0           0       7296/2197    3680/513       -3544/2565
    0           0       0           -845/4104       1859/4104
    0           0       0           0               -11/40
    ];
ARK45_F_weights = [
                    16/135 	0 	6656/12825 	28561/56430 	-9/50 	2/55
                    25/216 	0 	1408/2565 	2197/4104       -1/5 	0
                    ];
                
                
switch lower(method)
    case 'dormand-prince'
        A   = ARK45_DP_matrix;
        B   = ARK45_DP_weights;
        C   = ARK45_DP_nodes;
        
    case 'cash-karp'
        A   = ARK45_CK_matrix;
        B   = ARK45_CK_weights;
        C   = ARK45_CK_nodes;
        
    case 'fehlberg'
        A   = ARK45_F_matrix;
        B   = ARK45_F_weights;
        C   = ARK45_F_nodes;
        
    otherwise
        error('Supplied method does not exist')
end
            


end