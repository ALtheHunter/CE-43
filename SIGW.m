function [A, G, theta] = SIGW(Y, W, A0, G0, theta0,  L_min, pruning, N_iter, lr, th)
%IR_FAR Y = WH * A * G + n
% L_min : minimun 
[Q, M] = size(Y);
[N, L] = size(A0);
% Ini
A = A0;
G = G0;
G_prev = G;
theta = theta0;
row = (-(N - 1)/2:(N - 1)/2)' ;
% lambda = 10;

Res = Y - W'*A*G;
Rnorm = norm(Res, 'fro')^2;


c1 = 0.1;

for iter = 1:N_iter
%     iter

    % calculate gradient
    dtheta = df_dtheta(Y, W, A);
    
    lr_iter = lr;
    while true && lr_iter > 1e-6
        % update theta(A) and gains(G)
        theta_lr = theta - lr_iter * dtheta;

        A_lr = exp( 1j*  pi * row * theta_lr )/sqrt(N);
        R_lr = W' * A_lr;
        G_lr = pinv( R_lr'*R_lr ) * R_lr' * Y;
        
        Res2 = Y - W'*A_lr*G_lr;
        Rnorm2 = norm(Res2, 'fro')^2;
        if Rnorm2 < Rnorm - c1 * lr_iter * norm(dtheta, 'fro')^2
            break;
        else
            lr_iter = lr_iter * 0.5;
        end
    end
    
    % update theta(A) and gains(G)
    theta = theta - lr_iter * dtheta;

    A = exp( 1j*  pi * row * theta )/sqrt(N);

    R = W' * A;
    G =pinv( R'*R ) * R' * Y;
     
    % obtain res and res norm
    Res = Y - W'*A*G;
    Rnorm = norm(Res, 'fro')^2;
    

    gamma = norm(G - G_prev, 'fro');
    % pruning
    if L > L_min && iter > 10
        Gnorm = sum(abs(G).^2, 2);
        index = find(Gnorm > pruning);
        theta = theta(index);
        A = A(:, index);
        G = G(index, :);
        L = numel(index);
    end
    
    
    grad_norm = norm(dtheta, 'fro')/numel(dtheta);

    % early stopping
    if gamma < th
       break; 
    end

    % update previous G
    G_prev = G;
end

end

% function l = loss(Y)

function df = df_dtheta(Y, W, A)
% f = - tr{ Y'R(R'R + D/lambda)^(-1)R'Y}
% R = W'A
    dA = dA_dtheta(A);
    R = W'*A;
    E = R'*R ;
    EI = pinv(E);
    [N, L] = size(A);
    M = size(Y, 2);
    df = zeros(1, L);
    for l = 1:L
        dA_dl = zeros(N, L);
        dA_dl(:, l) = dA(:, l);
        dR_dl = W'*dA_dl;
        dE_dl = dR_dl' * R + R' * dR_dl;
        dEI_dl = - EI * dE_dl * EI;
        dP_dl = dR_dl * EI * R' + R * dEI_dl * R' + R * EI * dR_dl';
        df_dl = - trace( Y' * dP_dl * Y);
        df(l) = real(df_dl)/M;
    end
end


function dA = dA_dtheta(A)
    N = size(A, 1);
    row = (-(N - 1)/2:(N - 1)/2)' ;
    dA = (1j*pi*row).*A;
end
