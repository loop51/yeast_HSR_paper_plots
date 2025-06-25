clc;
clear;
l = 0.3048*2+ 1.2e-2;
load('detrend_data.mat');
load('processed_data.mat');


detrended = convert_to_real_imag(detrended);
s_filtered = convert_to_real_imag(s_filtered);

f = 10.5e9;

gamma_b = compute_gamma_l(detrended)/l;
gamma_p = compute_gamma_l(s_filtered)/l;

plot_alpha_beta(gamma_b, gamma_p, t);



%% 

% compute_complex_e(real(gamma_b(10)), imag(gamma_b(10)), f)
function plot_alpha_beta(base, peak, t)
    figure();
    subplot(2,1,1);
    plot(t, real(base));
    hold on; 
    plot(t, real(peak)); 
    grid on;
    ylabel('\alpha');
    
    subplot(2,1,2);
    plot(t, imag(base));
    hold on; 
    plot(t, imag(peak)); 
    grid on;
    ylabel('\beta');

    ax1 = subplot(2,1,1);
    ax2 = subplot(2,1,2);
    
    linkaxes([ax1,ax2],'x');

end

function [data] = compute_complex_e(alpha, beta, f)
    w = 2*pi*f;
    mu = 1.256637e-6;
    % syms epsilon sigma
    % 
    % x = (mu*epsilon/2);
    % % y = sqrt(1+(sigma/(w*epsilon))^2);
    % % equ1 = [alpha == w *sqrt(x*(y-1))];
    % % 
    % % equ2 = [beta == w *sqrt(x*(y+1))];
    % 
    % z= sigma/(w*epsilon);
    % y = 0.5*z^2 - 0.125*z^4 + (1/16)*z^6;
    % 
    % % equ1 = [alpha == w *sqrt(x*(y-1))];
    % equ1 = [alpha == w *sqrt(x*y)];
    % 
    % equ2 = [beta == w *sqrt(x*(y+2))];
    % 
    % 
    % s = vpasolve(equ1, equ2, [epsilon, sigma]);
    % 
    % e = s.epsilon - 1j*(s.sigma/w);
    
    % e = double(e);

    epsilon = (beta^2 - alpha^2)  / (mu * w^2);
    epsilon_double_prime = 2*alpha*beta/(mu * w^2);


    data = epsilon - 1j * epsilon_double_prime;
end
function [data] = compute_gamma_l(s)
    
    s11 = s.s11r + 1j * s.s11i;
    s21 = s.s21r+ 1j * s.s21i;
    for ii=1:length(s.s11r)    
        k(ii) = (((s11(ii)^2 - s21(ii)^2 +1)^2 -(2*s11(ii)))/(2*s21(ii))^2)^0.5;
    
        X(ii) = (1-s11(ii)^2+s21(ii)^2)/(2*s21(ii));
    end
    

    S_pos = log(X+k);
    S_neg = log(X-k);

    if real(S_pos)>=0
        data = S_pos;
    else
        data = S_neg;
    end
end

function [data] = convert_to_real_imag(s)
    s.s11m = 10.^(s.s11m/20);
    s.s11a  = deg2rad(s.s11a);

    s.s21m = 10.^(s.s21m/20);
    s.s21a  = deg2rad(s.s21a);

    [s.s11r, s.s11i] = pol2cart(s.s11a, s.s11m);
    [s.s21r, s.s21i] = pol2cart(s.s21a, s.s21m);
    data = s;
end


function [z] = compute_char_imp(s)
    s11 = s.s11r + 1j * s.s11i;
    s21 = s.s21r+ 1j * s.s21i;

    z0 = 50;

    a = ((1+s11)^2) - s21^2;
    b = ((1-s11)^2) -s21^2;
    c = z0^2 * a;
    z = sqrt(c/b);
end