% 2EiT_TO_P-17_2023-2024
% Temat projektu: Wyznaczanie charakterystyk częstotliwościowych czwórników RC i RL w oparciu o rozwiązywanie równań różniczkowych z zastosowaniem metody Rungego-Kutty
% Autorzy: Wojciech Kozub, Filip Żurek

clear; clc; close all;

% Wartości elementów
R = 1e6;            % [ohm]
L = 1e3;            % [H]
C = 1e-9;           % [F]

% Zakres częstotliwości
    % 0 Hz      fn = -1
    % 1 Hz      fn = 0
    % 1 MHz     fn = 6
fn_start = 0;       % Rząd wielkości początkowej częstotliwości (10^fn)
fn_stop = 3;        % Rząd wielkości końcowej częstotliwości (10^fn)

% Parametry symulacji
n = 30;             % Liczba częstotliwości testowych
ampl = 1;           % Amplituda sinusoidalnego napięcia wymuszającego [V]
fi = 0;             % Przesunięcie fazowe sinusoidalnego napięcia wymuszającego [rad]
u_c0 = 0;           % Warunek początkowy. Napięcie kondensatora [V]
i_l0 = 0;           % Warunek początkowy. Prąd cewki [A]
tau_n = 20;         % Liczba stałych czasowych (podstawa czasu)
min_tau = 10;       % Minimalna długość stanu nieustalonego. Jeśli mniej, ostrzega.
periods = 2;        % Liczba okresów stanu ustalonego
min_periods = 2*periods;        % Minimalna ilość okresów przebiegu wymuszającego (jeśli *tau_n stałych czasowych jest krótsze, niż długość *periods okresów)
p_const = 10000;     % Parametr "rozdzielczości/kroku". Krotność częstotliwości twierdzenia o próbkowaniu


% Dokładność ustalenia (z podręcznika);   Czas ustalenia;
% 63%                                       1   tau
% 90%                                       2.3 tau
% 99%                                       4.6 tau
% 99.9%                                     6.9 tau


% ---------------------RESZTA-KODU---------------------

% Menu
options = {'RC - dolnoprzepustowy (demo)', 'RC - górnoprzepustowy (demo)', 'RL - górnoprzepustowy (demo)', 'RC - dolnoprzepustowy', 'RC - górnoprzepustowy', 'RL - górnoprzepustowy'};
[choice, ~] = listdlg('Name','Menu','ListSize',[300,150],'OKString','Uruchom symulację','SelectionMode','single', 'PromptString', 'Wybierz czwórnik:', 'ListString', options);

if isempty(choice)
    disp('Anulowano wybór.');
    return;
end

switch choice
    % A, B - stałe równania różniczkowego
    case 1  % RC - całkujący demo
        R = 3.3e3;      % [ohm]
        C = 5.6e-12;    % [F]
        u_c0 = 1;       % [V]
        fn_start = 5;
        fn_stop = 9;
        A = -1 / (R * C);
        B = -A;
        tau = R * C;
        disp(['Rezystancja: ', num2str(R)]);
        disp(['Pojemność: ', num2str(C)]);

    case 2  % RC - różniczkujący demo
        R = 1.5e6;      % [ohm]
        C = 4.7e-9;     % [F]
        u_c0 = 0;       % [V]
        fn_start = 0;
        fn_stop = 3;
        A = -1 / (R * C);
        B = -A;
        tau = R * C;
        disp(['Rezystancja: ', num2str(R)]);
        disp(['Pojemność: ', num2str(C)]);

    case 3  % RL - różniczkujący demo
        R = 4.7e3;      % [ohm]
        L = 3.3e-3;     % [H]
        i_l0 = 0.5e-3;  % [A]
        fn_start = 4;
        fn_stop = 7;
        A = -1 * R / L;
        B = 1 / L;
        tau = L / R;
        disp(['Rezystancja: ', num2str(R)]);
        disp(['Indukcyjność: ', num2str(L)]);

    case 4  % RC - całkujący
        A = -1 / (R * C);
        B = -A;
        tau = R * C;
        disp(['Rezystancja: ', num2str(R)]);
        disp(['Pojemność: ', num2str(C)]);

    case 5  % RC - różniczkujący
        A = -1 / (R * C);
        B = -A;
        tau = R * C;
        disp(['Rezystancja: ', num2str(R)]);
        disp(['Pojemność: ', num2str(C)]);

    case 6  % RL - różniczkujący
        A = -1 * R / L;
        B = 1 / L;
        tau = L / R;
        disp(['Rezystancja: ', num2str(R)]);
        disp(['Indukcyjność: ', num2str(L)]);
end

f_span = logspace(fn_start, fn_stop, n);      % Zakres częstotliwości; generuje logarytmicznie rozmieszczone punkty
K = NaN(size(f_span));                        % Wzmocnienie
phase_shift = NaN(size(f_span));              % Przesunięcie fazowe
phase_shift_teor = NaN(size(f_span));         % Teoretyczne przesunięcie fazowe
stan_ustalony_po_n_tau = NaN(size(f_span));   % Liczba tau po jakiej przyjęto stan ustalony (wynika z wcześniejszych warunków)
tmax_t = tau * tau_n;                         % Końcowa wartość czasu 1
fg = 1 / (2 * pi * tau);                      % Częstotliwość 3dB
fig = figure('Name', options{choice}, 'NumberTitle', 'off', 'WindowState', 'maximized');  % Inicjalizacja okna dla wykresów

% Pętla główna
for i = 1:n
    dt = 1 / (p_const * f_span(i));     % Krok czasu
    tmax_T = min_periods/f_span(i);     % Końcowa wartość czasu 2
    tmax = max(tmax_t, tmax_T);         % Wybiera dłuższy czas 1 lub 2
    t = 0:dt:tmax;                      % Wektor czasu
    
    switch choice
        case {1,4}  % RC - całkujący
            u_out = rungego_kutty(A, B, ampl, fi, f_span(i), t, u_c0);   % Wyznaczenie przebiegu kondensatora
            correction_atan = -90;                                       % Korekta funkcji atan
    
        case {2,5}  % RC - różniczkujący
            u_c = rungego_kutty(A, B, ampl, fi, f_span(i), t, u_c0);    % Wyznaczenie napięcia kondensatora
            u_out = tau * diff(u_c) / dt;                               % Wyznaczenie napięcia rezystora
            t = t(1:end-1);                                             % Dostosowanie długości wektora czasu
            correction_atan = 0;                                      	% Korekta funkcji atan
                
        case {3,6}  % RL - różniczkujący
            i_l = rungego_kutty(A, B, ampl, fi, f_span(i), t, i_l0);    % Wyznaczenie prądu cewki
            u_out = L * diff(i_l) / dt;                                 % Wyznaczenie napięcia cewki
            t = t(1:end-1);                                             % Dostosowanie długości wektora czasu
            correction_atan = 0;                                        % Korekta funkcji atan
    end
    
    we = u_in(ampl, f_span(i), t, fi);                                  % Wektor przebiegu wejściowego

    % Wykres jednego wymuszenia ze stanem nieustalonym
    t_per_tau = t/tau;  % etytkiety osi x
    clf;                % Odświeżenie wykresów
    subplot(2,3,1);
        plot(t_per_tau, u_out, 'b-'); hold on; grid on;
        title('Rozwiązanie równania różniczkowego');
        xlabel('Liczba stałych czasowych (t/\tau)');
        ylabel('u_{wy}(t) [V]');
    
    T = 1/f_span(i);                % Okres
    ratio = 1 - periods/(tmax/T);   % Usuwana część początku przebiegu

    % Wycięcie końcówek przebiegów (stan ustalony)
    u_out = cut(u_out, ratio);
    we = cut(we, ratio);
    t = cut(t, ratio);
    
    % Wykres jednego wymuszenia (stan ustalony)
    t_per_tau = t/tau;  % etytkiety osi x
    subplot(2,3,4);
        plot(t_per_tau, u_out, 'b-'); hold on; grid on;
        plot(t_per_tau, we, 'r-');
        title('Stan ustalony');
        xlabel('Liczba stałych czasowych (t/\tau)');
        ylabel('u(t) [V]');
%        legend(['u_{wy}(t)'], ['u_{we}(t)']); % spowalnia symulacje
    
    K(i) = max(u_out) / ampl;                                                       		% Wyznaczenie wzmocnienia
    phase_shift(i) = shift(t, we, u_out, f_span(i), choice);                                % Wyznaczenie przesunięcia fazowego
    phase_shift_teor(i) = atan(1 / (2 * pi * f_span(i) * tau)) * 180 / pi + correction_atan;% Wyznaczenie teoretycznego przesunięcia fazowego

    % Rysowanie charakterystyk Bodego
    ax1 = subplot(2,3,[2,3]);
        semilogx(f_span, 20 * log10(K), 'b.-');  grid on;
        xl = xline(fg,'--','f_{3dB}');
        xl.LabelVerticalAlignment = 'middle';
        xl.LabelHorizontalAlignment = 'center';
        ylabel('Wzmocnienie [dB]');
        title('Charakterystyka Bodego');
        xlim([min(f_span),max(f_span)]);
    
    ax2 = subplot(2,3,[5,6]);
        semilogx(f_span, phase_shift, 'b.-'); grid on; hold on;
        semilogx(f_span, phase_shift_teor, 'r.');
        xl = xline(fg,'--','f_{3dB}');
        xl.LabelVerticalAlignment = 'bottom';
        xl.LabelHorizontalAlignment = 'center';
        ylabel('Faza [{\circ}]');
        xlabel('f [Hz]');
        xlim([min(f_span),max(f_span)]);
        legend(['Wyznaczona'], ['Teoretyczna']);    % spowalnia symulacje
    linkaxes([ax1, ax2], 'x');

    stan_ustalony_po_n_tau(i) = min(t_per_tau);
end

% Sprawdzenie minimalnej długości stanu nieustalonego
[SN, SN_i] = min(stan_ustalony_po_n_tau);
if SN <= min_tau
    disp(['Ostrzeżenie. Dla częstotliwości wejścia: ', num2str(round(f_span(SN_i))), 'Hz, stan nieustalony trwał tylko ', num2str(SN), ' tau.']);
    warndlg(['Dla częstotliwości wejścia: ', num2str(round(f_span(SN_i))), 'Hz, stan nieustalony trwał tylko ', num2str(SN), ' tau.'],'Ostrzeżenie','modal');
else
    msgbox('Ukończono obliczenia!');
end

% bode = figure('Name', options{choice});  % Inicjalizacja okna dla bodego
%     % Rysowanie charakterystyk Bodego
%     ax1 = subplot(2,1,1);
%         semilogx(f_span, 20 * log10(K), 'b.-');  grid on;
%         xl = xline(fg,'--','f_{3dB}');
%         xl.LabelVerticalAlignment = 'middle';
%         xl.LabelHorizontalAlignment = 'center';
%         ylabel('Wzmocnienie [dB]');
%         title('Charakterystyka Bodego');
%         xlim([min(f_span),max(f_span)]);
%     
%     ax2 = subplot(2,1,2);
%         semilogx(f_span, phase_shift, 'b.-'); grid on; hold on;
%         semilogx(f_span, phase_shift_teor, 'r.');
%         xl = xline(fg,'--','f_{3dB}');
%         xl.LabelVerticalAlignment = 'bottom';
%         xl.LabelHorizontalAlignment = 'center';
%         ylabel('Faza [{\circ}]');
%         xlabel('f [Hz]');
%         xlim([min(f_span),max(f_span)]);
%         legend(['Wyznaczona'], ['Teoretyczna']);    % spowalnia symulacje

% saveas(fig, ['wykres ', options{choice}], 'epsc');  % zapis grafiki wektorowej
% saveas(bode, ['bode ', options{choice}], 'epsc');

% -------------------------------------------------------------------------
% Funkcje

% Algorytm Rungego-Kutty czwartego rzędu
function x = rungego_kutty(A, B, a_gen, fi, f, t, x0)
    h = t(2) - t(1);
    x = NaN(size(t));
    x(1) = x0;    % Warunek początkowy (w chwili t0)
  
    for n = 1:length(t) - 1
        k1 = f_ode(t(n),          x(n),                A, B, a_gen, f, fi);
        k2 = f_ode(t(n) + h / 2,  x(n) + h * k1 / 2,   A, B, a_gen, f, fi);
        k3 = f_ode(t(n) + h / 2,  x(n) + h * k2 / 2,   A, B, a_gen, f, fi);
        k4 = f_ode(t(n) + h,      x(n) + h * k3,       A, B, a_gen, f, fi);
        x(n + 1) = x(n) + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    end
end

% Równanie różniczkowe
function dx = f_ode(t, x, A, B, a, f, fi)
    dx = A * x + B * u_in(a, f, t, fi);
end

% Napięcie wymuszające 
function u = u_in(a, f, t, fi)
    u = a * sin(2 * pi * f * t + fi);
end

% Zwraca przesunięcie fazowe między sygnałami x1 i x2 w stopniach
function phi = shift(t, x1, x2, f, choice)
    ref = islocalmax(x1);
    s = islocalmax(x2);
    
    % chwile czasowe szczytów
    p1 = t(ref);    
    p2 = t(s);

    % wyrównanie ilości elementów
    l = min(length(p1),length(p2));
    p1 = p1(1:l);
    p2 = p2(1:l);

    dif = p1 - p2;  % przesunięcia [s]
    
    phi = min(dif * 360 * f);   % przesunięcie [stopnie]
    
    % Niwelacja błędu (+/- 360)
    switch choice
        case {1,4}  % dolnoprzepustowy
            if (phi > 0)
                phi = phi - 360;
            elseif (phi < -90)
                phi = phi + 360;
            end                                         
                
        case {2,3,5,6}  % górnoprzepustowy
            if (phi > 90)
                phi = phi - 360;
            elseif (phi < 0)
                phi = phi + 360;
            end
    end
end

% Wycięcie końca przebiegu
function x_cut = cut(x, p)
    l = length(x);
    x_cut = x(round(p*l)+1 : end-1);
end
