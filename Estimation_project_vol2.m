clear all;
clear figure;
%Constants
N = 513;
n_0 = -256;
T = 10.^(-6);
f_0 = 10.^5;
F_s = 10^6;
w_0 = 2.*pi.*f_0;
A = 1;
k = 10:2:20;
phi = pi./8;
n = n_0:n_0+N-1;
P = (N.*(N-1))./(2);
Q = (N.*(N-1).*(2.*N-1))./(6);

%Making a vector of sigma based on snr
snr = -10:10:60;

sigma = zeros(1,8);

q = 1;

while q~= 9
    sigma(q)= sqrt(A.^2./(db2mag(snr(q)).*2));
    q = q + 1;
end 

%The signal of x
signal = A.*exp(1i.*(w_0.*n.*T+phi));

steps = 50;

w_fft = zeros(1, steps);
phi_fft = zeros(1, steps);
var_w_hat = zeros(1, 8);
var_phi_hat = zeros(1, 8);
phi_error = zeros(1,8);
w_error = zeros(1,8);
CRLB_w = zeros(1, 8);
CRLB_phi = zeros(1,8);
x_fft_matrix = zeros(8, N);

for j=1:6
   M = 2.^k(j);
   
    for j_1=1:8
        
           noise = normrnd(0, sigma(j_1), steps, N) + 1i.*normrnd(0, sigma(j_1), steps, N);
           
           for j_2=1:steps
               
                x = signal + noise(j_2,:);
    
                x_fft = fft(x,M);
    
                [argvalue, argmax] = max(abs(x_fft));
                m = argmax;
                
    
                w_fft(j_2) = (2.*pi.*m)./(M.*T);
                

                phi_fft(j_2) = angle(exp(-1i*w_fft(j_2)*n_0*T).*(x_fft(m)));
           end 
           
           CRLB_w(j_1) = (12.*(sigma(j_1)).^2)/(A.^2.*T.^2.*N.*(N.^2-1));
           CRLB_phi(j_1) = ((12.*(sigma(j_1)).^2).*(n_0.^2.*N+2.*n_0.*P+Q))./(A.^2.*N.^2.*(N.^2-1));
           
           %fprintf('SNR: %f \n', snr(j_1));
           %fprintf('M: %f \n', M);
           %fprintf('Mean of phi_fft: %f \n', mean(phi_fft));
           
           w_error = abs(w_0 - w_fft);
           
           var_w_hat(j_1) = var(w_error);
           
           fprintf('SNR: %f \n', snr(j_1));
           fprintf('M: %f \n', M);
           fprintf('omega fft: %f \n', mean(w_fft));
           fprintf('omega 0: %f \n', w_0);
           fprintf('omega error: %f \n', mean(w_error));
           fprintf('omega variance: %f \n', var_w_hat(j_1));
           

           %phi_error(j_1) = mean((phi - phi_fft).^2);
           %var_phi_hat(j_1) = var(phi_error);
           
          if k(j) == 10
            func = @(w, p) sum(abs(x - A*exp(1i*(w*n*T + phi))));

            [vals, ~, exitflag, output] = fminsearch(@(input) func(input(1), input(2)), [mean(w_fft), mean(phi_fft)]);

            w_mle = vals(1);
            phi_mle = vals(2);

            %fprintf('Estimated w: %.0f (%.3f%% off)\n', w_mle, 100*(w_mle - w_0)/w_0);
            %fprintf('Estimated phi: %f (%.3f%% off)\n', phi_mle,100*(phi_mle - phi)/phi)
            %fprintf('SNR: %f \n', snr(j_1)); 
            
          end 
    
    end
    %figure(1); 
    %semilogy(snr, CRLB_w, 'r'); hold on;
    %semilogy(snr, var_w_hat, 'g'); hold on;
    %legend({'CRLB','var({\omega}_{fft})'},'Location','northeast')
    %grid
    
    %figure(2); semilogy(snr, CRLB_phi, 'r'); hold on;
    %semilogy(snr, var_phi_hat, 'b'); hold on;
    %legend({'CRLB', 'var({\phi}_{fft})'}, 'Location', 'northeast'); 
    %grid
end

