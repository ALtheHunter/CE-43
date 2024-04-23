clear all; 
clc
savefile = 1;
SNR_dB = 10;  SNR_linear=10.^(SNR_dB/10.);
sigma2=1/SNR_linear;
N_iter = 10; 
%%% system parameters
N = 256; % number of beams (transmit antennas)
K = 4; % number of users
N_RF = 4; % number of RF chains
M = 4; % number of subcarriers
L = 6; % number of paths per user


fc = 100e9; % carrier frequency
R = linspace(5,150,50);
sector = pi * 2 / 3;
fs = 100e6; % bandwidth
Q = 32;  % number of pilot blocks
tmax = 20e-9; % maximum of the path delay
f = zeros(1,M);
for m = 1:M
    f(m)=fc+fs/(M)*(m-1-(M-1)/2);
end
c = 3e8;
lambda_c = c/fc;
d = lambda_c / 2;
eps = 1e-3;

NMSE = zeros(1,length(R));
NMSE0 = zeros(1,length(R));
NMSE1 = zeros(1,length(R));
NMSE2 = zeros(1,length(R));
NMSE3 = zeros(1,length(R));
NMSE4 = zeros(1,length(R));
NMSE5 = zeros(1,length(R));


%DFT    
s = 2;
D = s*N; %字典规模
row = (-(N - 1)/2:(N - 1)/2)' ;
col = -1 + 2/D : 2/D : 1 ;
DFT  =  exp( 1j*  pi * row * col ) / sqrt(N);

% Quacode
rho = 3;
beta = 1.2;
rho_max = 64;
[Polarcodebook,label] = PolarCodeBook(N, s, d, lambda_c, beta, rho,rho_max);
S = size(Polarcodebook, 2);


R_len = length(R);
t0 = clock;
for i_r = 1:R_len
    Rmin = R(i_r);
    Rmax = Rmin;

    error = 0; error0 = 0; error1 = 0; error2 = 0; error3 = 0; error4 = 0;error5 = 0;
    for iter = 1:N_iter
        fprintf('R = %.4f[%d/%d] | iteration:[%d/%d] | run %.4f s\n', Rmin, i_r, R_len, iter, N_iter, etime(clock, t0)); 
        
        % Wideband spatial channel
        [H, hc, r, theta, G] = near_field_channel(N, K, L, d, fc, fs, M, Rmin, Rmax,sector,1);
        H = channel_norm(H);
        Phi = (2*randi(2,Q*N_RF,N) - 3)/sqrt(N);
        noise = sqrt(sigma2)*(randn(Q*N_RF,M)+1i*randn(Q*N_RF,M))/sqrt(2); 
        
        for k = 1 : K
            Hsf =  reshape(H(k, :, :), [N, M]);    % Hsf = F*Haf
            
            % adaptive selecting matrix
            Z = Phi*Hsf + noise;
            Znorm = norm(Z, 'fro');
            
            %% Oracle LS
             Hsf_hat0 = Oracle_LS(Z, Phi, N, M, r(k,:), theta(k,:), d, fc);
        
            %% Far field on-grid
            [Haf_hat1, sup1] = SOMP(Z, Phi*DFT, L*2, D, M);
            
            %% Far field off-grid
            [A2, G2] = SIGW( Z/Znorm, Phi'/Znorm, DFT(:, sup1), Haf_hat1(sup1,:), col(sup1),...
                L,  3, 20, 1, 1e-8);
            Hsf_hat2 = A2*G2;
            
            %% Near field on-grid
            [Haf_hat3,sup3] = SOMP(Z, Phi*Polarcodebook,L*2, S, M);
            
            %% Near fild off-grid
            [A4, G4, r4, theta4] =P_SIGW( Z/Znorm, Phi'/Znorm, Polarcodebook(:, sup3), Haf_hat3(sup3,:), label(2, sup3), label(1, sup3), ...
            fc, d, L, 3, 20, 10, 1, 1e-8);
            Hsf_hat4 = A4*G4;
            
            
            %% LS
             Hsf_hat5 = Hsf + ( sqrt(sigma2)*(randn(N,M)+1i*randn(N,M))/sqrt(2) );
            




            error0 = error0 + norm(Hsf - Hsf_hat0,'fro')^2/norm(Hsf,'fro')^2;
            error1 = error1 + norm(Hsf - DFT*Haf_hat1,'fro')^2/norm(Hsf,'fro')^2;
            error2 = error2 + norm(Hsf - Hsf_hat2,'fro')^2/norm(Hsf,'fro')^2;
            error3 = error3 + norm(Hsf - Polarcodebook*Haf_hat3,'fro')^2/norm(Hsf,'fro')^2;
            error4 = error4 + norm(Hsf - Hsf_hat4,'fro')^2/norm(Hsf,'fro')^2;
            error5 = error5 + norm(Hsf - Hsf_hat5,'fro')^2/norm(Hsf,'fro')^2;
        end     
    end
    NMSE0(i_r) = error0/K/N_iter;
    NMSE1(i_r) = error1/K/N_iter;
    NMSE2(i_r) = error2/K/N_iter;
    NMSE3(i_r) = error3/K/N_iter;
    NMSE4(i_r) = error4/K/N_iter;
    NMSE5(i_r) = error5/K/N_iter;
end

x = 1:2:size(R,2);
figure; hold on; box on; grid on;
plot(R(x),10*log10(NMSE1(x)),'-o','Linewidth',1.2,'markersize',5,'color',[0.4660 0.6740 0.1880],'MarkerFaceColor','w');
plot(R(x),10*log10(NMSE2(x)),'-d','Linewidth',1.2,'markersize',5,'color',[0.4660 0.6740 0.1880],'MarkerFaceColor','w');
plot(R(x),10*log10(NMSE3(x)),'-s','Linewidth',1.2,'markersize',5,'color',[0.8500 0.3250 0.0980],'MarkerFaceColor','w');
plot(R(x),10*log10(NMSE4(x)),'->','Linewidth',1.2,'markersize',5,'color',[0.8500 0.3250 0.0980],'MarkerFaceColor','w');
plot(R(x),10*log10(NMSE5(x)),'b-<','Linewidth',1.2,'markersize',5,'MarkerFaceColor','w');
plot(R(x),10*log10(NMSE0(x)) ,'m--','Linewidth',1.2,'markersize',5,'MarkerFaceColor','w');
xlabel('Distance (meters)');
ylabel('NMSE (dB)');
ylim([-30, 0]);
legend( 'SWOMP','SS-SIGW-OLS','Proposed P-SOMP','Proposed P-SIGW','LS', 'Oracle LS');

