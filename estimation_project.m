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

%while q~= 9
  %  sigma(q)= sqrt(A.^2./(db2mag(snr(q)).*2));
 %   q = q + 1;
%end 

%The signal of x
signal = A.*exp(1i.*(w_0.*n.*T+phi));

steps = 50;

w_fft = zeros(1, steps);
phi_fft = zeros(1, steps);
var_w_hat = zeros(1, 8);
var_phi_hat = zeros(1, 8);
phi_error = zeros(1, steps);
w_error = zeros(1, steps);
CRLB_w = zeros(1, 8);
CRLB_phi = zeros(1,8);
w_mle = zeros(1, 8);
phi_mle = zeros(1, 8);
var_w_min = zeros(1,8);


for j=1:6
   M = 2.^k(j);
   
    for j_1=1:8
        
           sigma(j_1)= sqrt(A.^2./(db2mag(snr(j_1)).*2));
           noise = normrnd(0, sigma(j_1), steps, N) + 1i.*normrnd(0, sigma(j_1), steps, N);
           
           
           for j_2=1:steps
               
                x = signal + noise(j_2,:);
    
                x_fft = fft(x,M);
    
                [argvalue, argmax] = max(abs(x_fft));
                m = argmax;
                
    
                w_fft(j_2) = (2.*pi.*m)./(M.*T);
                w_error(j_2) =(w_0 - w_fft(j_2)).^2;
                

                phi_fft(j_2) = angle(exp(-1i*w_fft(j_2)*n_0*T).*(x_fft(m)));
                phi_error(j_2) = (phi - phi_fft(j_2)).^2;
           end 
           
           CRLB_w(j_1) = (12.*(sigma(j_1)).^2)/(A.^2.*T.^2.*N.*(N.^2-1));
           CRLB_phi(j_1) = ((12.*(sigma(j_1)).^2).*(n_0.^2.*N+2.*n_0.*P+Q))./(A.^2.*N.^2.*(N.^2-1));
           
           var_error_w(j_1) =(1/steps)*(sum(w_error));
           var_error_phi(j_1) = (1/steps)*(sum(phi_error));
           
           
          if k(j) == 10
            func = @(w, p) sum(abs(x - A*exp(1i*(w*n*T + phi))));

            [vals, ~, exitflag, output] = fminsearch(@(input) func(input(1), input(2)), [mean(w_fft), mean(phi_fft)]);

            w_mle(j_1) = vals(1);
            
            phi_mle(j_1) = vals(2);
            
            var_w_min(j_1) = var(w_0 - w_mle(j_1));
            
            

            fprintf('Estimated w: %.0f (%.3f%% off)\n', w_mle(j_1), 100*(w_mle(j_1) - w_0)/w_0);
            fprintf('Estimated phi: %f (%.3f%% off)\n', phi_mle(j_1),100*(phi_mle(j_1) - phi)/phi)
            fprintf('SNR: %f \n', snr(j_1)); 
            
          end 
    
    end
    figure(1); 
    %semilogy((snr), CRLB_w ); hold on;
    %semilogy((snr), var_error_w ); hold on;
    grid
    
    %figure(2); 
    %semilogy((snr), CRLB_phi ); hold on;
    %semilogy((snr), var_error_phi ); hold on; 
    grid
    
    figure(3)
    semilogy((snr), CRLB_w, 'r'); hold on;
    semilogy((snr), var_w_min, 'g'); hold on;
    grid;
    
    figure(4)
    semilogy((snr), CRLB_phi, 'r'); hold on;
    semilogy((snr), phi_mle, 'b'); hold on;
    grid;
end

%plotting
figure(1)
legend('CRLB', '10', '12', '14', '16', '18', '20');
figure(2)
legend('CRLB', '10', '12', '14', '16', '18', '20');