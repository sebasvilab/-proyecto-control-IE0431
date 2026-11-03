%% IE-0431 Sistemas de Control - Proyecto final II-2026
% Equipo 2-10 | Planta: control de temperatura (TMP36 + TIP31C)
% Punto de operacion asignado: 79-27-60

% Puntos que se toman en esta código:

%   1) Carga y normaliza los datos experimentales (0-100 %)

%   2) Identifica 4 modelos: 3 continuos + 1 discreto (minimos cuadrados)
%        M1: System Identification Toolbox (P1D)
%        M2: Metodo de Alfaro 123c (POMTM)
%        M3: Metodo de Smith (dos puntos, POMTM)
%        M4: Minimos cuadrados, ARX discreto de orden 3

%   3) Compara los modelos (IAE, ISE y % de ajuste tipo toolbox)

%   4) Sintoniza 4 controladores PI (3 continuos + 1 discreto):
%        C1: PID Tuner de MATLAB (pidtune)
%        C2: SIMC de Skogestad (2003)
%        C3: AMIGO de Astrom y Hagglund
%        C4: PI discreto en sisotool con el modelo de minimos cuadrados

%   5) Simula servocontrol y regulatorio y mide Mp, ta2%, ep, IAE, TVu, Ms

clear; close all; clc;





% ---------------------------------------------
%% Configuración
% ---------------------------------------------


ARCHIVO_CSV = 'tclab_ident_data_Q1_v_79-27-60.txt' ;

% Archivo TCLab: 1 = tiempo [s], 2 = Heater 1 [%], 3 = Heater 2 [%],
%                4 = Temperature 1 [C], 5 = Temperature 2 [C]
% Se usa el lazo 1: u = Heater 1 (Q1), y = Temperature 1 (T1).
COL_T = 1; COL_U = 2; COL_Y = 4;
% Prueba realizada: Q1 = 0 -> 79 -> 27 -> 60 %. Se identifica con el ultimo
% escalon (27 -> 60 %), que parte de un estado casi estacionario.
N_PREVIO = 120;          % muestras antes del escalon que se conservan
% Normalizacion: el enunciado recomienda trabajar todo en 0-100 %.
Y_EN_CELSIUS = true;     % true si la salida viene en grados C
U_EN_PORCENTAJE = true;  % true si la senal de control ya viene en 0-100 %
U_MAX_FISICO = 255;      % solo se usa si U_EN_PORCENTAJE = false (ej. PWM 0-255)
T_MIN = -40; T_MAX = 150; % rango del sensor TMP36 [C]

% Escalones de prueba para las simulaciones de lazo cerrado [%]
DR = 5;   % cambio en el valor deseado (servocontrol)
DD = 5;   % cambio en la perturbacion (regulatorio)

% Especificaciones de diseno (punto 5 del enunciado). Son una propuesta:
% deben justificarse con la aplicacion industrial investigada.
ESP.Mp_max  = 10;    % sobrepaso maximo [%]
ESP.ta2_max = 500;   % tiempo de asentamiento al 2 % [s] (lazo abierto ~ 4*tau + L ~ 800 s)
ESP.ep_max  = 0;     % error permanente ante escalon [%] (servo y regulador)
ESP.Ms_max  = 1.8;   % robustez: sensibilidad maxima (1.2 - 2.0 es el rango usual)





%% ---------------------------------------------
% 1. Carga de datos
% ---------------------------------------------


% Busca el csv en la carpeta de este script. Si el nombre de ARCHIVO_CSV no
% existe y en la carpeta hay un solo .csv, usa ese automaticamente.
carpeta = fileparts(mfilename('fullpath'));
if isempty(carpeta), carpeta = pwd; end
ruta = fullfile(carpeta, ARCHIVO_CSV);
if ~isfile(ruta)
    lista = [dir(fullfile(carpeta, '*79-27-60*')); dir(fullfile(carpeta, '*.csv'))];
    if numel(lista) == 1
        ruta = fullfile(carpeta, lista(1).name);
    elseif isempty(lista)
        error(['No hay ningun archivo .csv en la carpeta:\n  %s\n' ...
               'Copie ahi el csv de su equipo (el que viene con la asignacion).'], carpeta);
    else
        error(['Hay varios .csv en la carpeta y ninguno se llama "%s":\n  %s\n' ...
               'Escriba el nombre correcto en ARCHIVO_CSV (seccion 0).'], ...
               ARCHIVO_CSV, strjoin({lista.name}, ', '));
    end
end

fprintf('Datos cargados de: %s\n', ruta);
datos = readmatrix(ruta);
datos = datos(all(~isnan(datos),2),:);   % quita filas vacias o encabezados
t = datos(:,COL_T); u = datos(:,COL_U); y = datos(:,COL_Y);
t = t - t(1);

if Y_EN_CELSIUS
    y = (y - T_MIN)/(T_MAX - T_MIN)*100;   % C -> %
end
if ~U_EN_PORCENTAJE
    u = u/U_MAX_FISICO*100;
end

Ts = median(diff(t));
if max(abs(diff(t) - Ts)) > 0.05*Ts      % remuestreo si el muestreo no es uniforme
    tn = (0:Ts:t(end))';
    u = interp1(t,u,tn,'previous'); y = interp1(t,y,tn,'linear');
    t = tn;
end
fprintf('Tiempo de muestreo de los datos: Ts = %.4f s\n', Ts);

figure('Name','Prueba completa');
subplot(2,1,1); plot(t, y*(T_MAX-T_MIN)/100 + T_MIN, 'LineWidth', 1.2); grid on;
ylabel('T_1 [C]'); title('Prueba experimental completa (Q1: 0 \rightarrow 79 \rightarrow 27 \rightarrow 60 %)');
subplot(2,1,2); stairs(t, u, 'LineWidth', 1.2); grid on;
ylabel('Q_1 [%]'); xlabel('Tiempo [s]');

% Tramo para identificar: desde N_PREVIO muestras antes del ultimo escalon
cambios = find(abs(diff(u)) > 0.5) + 1;
kult = cambios(end);
kini = max(1, kult - N_PREVIO);
t = t(kini:end) - t(kini); u = u(kini:end); y = y(kini:end);





% ---------------------------------------------
%% 2. Punto de operacion y variables de desviacion
% ---------------------------------------------


% Instante del escalon: primera muestra donde u cambia de forma apreciable
k0 = find(abs(u - u(1)) > 0.5, 1);
if isempty(k0), error('No se detecto el escalon en la senal de control.'); end
t0 = t(k0);

u0 = mean(u(1:max(k0-1,1)));
y0 = mean(y(1:max(k0-1,1)));
nfin = max(5, round(0.05*numel(y)));      % ultimo 5 % de datos = estado estacionario
uf = mean(u(end-nfin+1:end));
yf = mean(y(end-nfin+1:end));

du = u - u0;  dy = y - y0;
Du = uf - u0; Dy = yf - y0;
Kp_proc = Dy/Du;

fprintf('\nPunto de operacion (de los datos): u0 = %.2f %%, y0 = %.2f %% (%.2f C)\n', ...
        u0, y0, y0*(T_MAX-T_MIN)/100 + T_MIN);
fprintf('Escalon: Du = %.2f %%, Dy = %.2f %%, ganancia K = %.4f %%/%%\n', Du, Dy, Kp_proc);

figure('Name','Curva de reaccion');
subplot(2,1,1); plot(t,y,'LineWidth',1.2); grid on;
ylabel('y [%]'); title('Curva de reaccion del proceso (lazo abierto)');
subplot(2,1,2); stairs(t,u,'LineWidth',1.2); grid on;
ylabel('u [%]'); xlabel('Tiempo [s]');





% ---------------------------------------------
%% 3. Tiempos caracteristicos de la curva (sobre la salida suavizada)
% ---------------------------------------------


ys = movmean(dy, max(3, round(numel(dy)/200)));   % suaviza ruido de cuantizacion
yn = ys/Dy;                                        % respuesta normalizada 0 -> 1
tp = t - t0;                                       % tiempo desde el escalon
tx = @(p) cruce(tp(k0:end), yn(k0:end), p);        % tiempo al p (fraccion) de Dy

t25 = tx(0.25); t75 = tx(0.75);
t28 = tx(0.283); t63 = tx(0.632);





% ---------------------------------------------
%% 4. Identificacion de modelos
% ---------------------------------------------


s = tf('s');

% --- M1: System Identification Toolbox, process model P1D ---
dat = iddata(dy, du, Ts);
opt = procestOptions('InitialCondition','zero');   % SearchMethod 'auto' (interno)
M1id = procest(dat, 'P1D', opt);
P1 = tf(M1id);
P1 = tf(P1.Numerator, P1.Denominator, 'InputDelay', M1id.Td);

% --- M2: Alfaro 123c, primer orden mas tiempo muerto ---
tau2 = 0.9102*(t75 - t25);
L2   = max(0, 1.2620*t25 - 0.2620*t75);
P2   = Kp_proc/(tau2*s + 1); P2.InputDelay = L2;

% --- M3: Smith (28.3 % y 63.2 %), primer orden mas tiempo muerto ---
tau3 = 1.5*(t63 - t28);
L3   = max(0, t63 - tau3);
P3   = Kp_proc/(tau3*s + 1); P3.InputDelay = L3;

% --- M4: Minimos cuadrados, ARX de orden 3 ---
% Se submuestrea para que haya ~20 muestras por constante de tiempo y la
% regresion quede bien condicionada.
na = 3; nb = 3;
Nd = max(1, round((tau2/20)/Ts));
Tls = Nd*Ts;
idx = 1:Nd:numel(t);
yk = dy(idx); uk = du(idx); tk = t(idx);
nk = max(1, round(L2/Tls));                  % retardo discreto [muestras]
kini = max(na, nb+nk-1) + 1;
N = numel(yk);
Phi = zeros(N-kini+1, na+nb); Yr = yk(kini:N);
for k = kini:N
    Phi(k-kini+1,:) = [-yk(k-1:-1:k-na)', uk(k-nk:-1:k-nk-nb+1)'];
end
theta = Phi\Yr;                              % solucion de minimos cuadrados
a = [1, theta(1:na)'];
b = [zeros(1,nk), theta(na+1:end)'];
P4 = tf(b, a, Tls, 'Variable','z^-1');

fprintf('\n--- Modelos identificados (desviacion, %%/%%) ---\n');
fprintf('M1 Toolbox P1D : K=%.4f  tau=%.2f s  L=%.2f s\n', M1id.Kp, M1id.Tp1, M1id.Td);
fprintf('M2 Alfaro 123c : K=%.4f  tau=%.2f s  L=%.2f s\n', Kp_proc, tau2, L2);
fprintf('M3 Smith       : K=%.4f  tau=%.2f s  L=%.2f s\n', Kp_proc, tau3, L3);
fprintf('M4 MC ARX(3,3,%d), Ts=%.3f s\n', nk, Tls);

% Funciones de transferencia completas (quedan en el workspace como P1..P4)
disp('M1 - Toolbox P1D:'); P1 %#ok<NOPTS>
disp('M2 - Alfaro 123c:'); P2 %#ok<NOPTS>
disp('M3 - Smith:');       P3 %#ok<NOPTS>
disp('M4 - Minimos cuadrados (discreto):'); P4 %#ok<NOPTS>

% Comparacion con el toolbox (% de ajuste de cada modelo)
figure('Name','compare - continuos'); compare(dat, P1, P2, P3);
dat_ls = iddata(yk, uk, Tls);
figure('Name','compare - discreto');  compare(dat_ls, P4);





% ---------------------------------------------
%% 5. Comparacion de modelos
% ---------------------------------------------


ym1 = lsim(P1, du, t);
ym2 = lsim(P2, du, t);
ym3 = lsim(P3, du, t);
ym4 = interp1(tk, lsim(P4, uk, tk), t, 'previous', 'extrap');

Ym = [ym1 ym2 ym3 ym4];
nombres = {'M1 Toolbox P1D','M2 Alfaro 123c','M3 Smith','M4 MC ARX3'};
E   = dy - Ym;
IAE = trapz(t, abs(E))';
ISE = trapz(t, E.^2)';
FIT = (100*(1 - vecnorm(E)./norm(dy - mean(dy))))';   % misma formula (NRMSE) de compare()

Tmod = table(IAE, ISE, FIT, 'RowNames', nombres, ...
             'VariableNames', {'IAE','ISE','Ajuste_pct'});
disp(' '); disp('--- Comparacion de modelos ---'); disp(Tmod);

figure('Name','Comparacion de modelos');
plot(t, dy, 'k', 'LineWidth', 1.5); hold on;
plot(t, Ym, 'LineWidth', 1.2); grid on;
legend(['Datos reales', nombres], 'Location', 'southeast');
xlabel('Tiempo [s]'); ylabel('\Delta y [%]');
title('Respuesta de los modelos vs. proceso real');

% Mejor modelo continuo (los controladores continuos se sintonizan con el)
[~, ib] = min(IAE(1:3));
Pc = {P1, P2, P3}; P = Pc{ib};
[Kb, taub, Lb] = pomtm(P);
fprintf('Mejor modelo continuo (menor IAE): %s\n', nombres{ib});





% ---------------------------------------------
%% 6. Sintonizacion de controladores PI
% ---------------------------------------------


% Accion del controlador: al subir Q1 sube T1 (K > 0), entonces el
% controlador es de ACCION INVERSA: e = r - y  y  u = C(s)*e.
% Estructura: PI de un grado de libertad. Se escoge PI porque el proceso es
% de primer orden con tiempo muerto, la accion integral da ep = 0 ante
% escalones en r y d, y la derivativa amplificaria el ruido de cuantizacion
% del TMP36 sin mejora apreciable. Td = 0 entra en la estructura PID estandar.
if Kb > 0, accion = 'inversa'; else, accion = 'directa'; end
fprintf('\nGanancia del proceso K = %.4f > 0  ->  controlador de accion %s\n', Kb, accion);

% C1: PID Tuner de MATLAB, margen de fase de 60 grados
% pidtune es el mismo algoritmo de la app PID Tuner. Para las capturas del
% informe, abrir la app con:  pidTuner(P, 'PI')
optT = pidtuneOptions('PhaseMargin', 60, 'DesignFocus', 'balanced');
C1 = pidstd(pidtune(P, 'PI', optT));

% C2: SIMC de Skogestad, tau_c = L (sintonia "ajustada", buen compromiso
% desempeno-robustez). Si L es muy pequeno se limita tau_c >= 0.1*tau.
tauc = max(Lb, 0.1*taub);
Kc2 = taub/(Kb*(tauc + Lb));
Ti2 = min(taub, 4*(tauc + Lb));
C2 = pidstd(Kc2, Ti2);

% C3: AMIGO de Astrom y Hagglund (2004), PI para POMTM. Busca Ms ~ 1.4,
% es decir, prioriza robustez y suavidad de la senal de control.
Kc3 = 0.15/Kb + (0.35 - Lb*taub/(Lb + taub)^2)*taub/(Kb*Lb);
Ti3 = 0.35*Lb + 13*Lb*taub^2/(taub^2 + 12*Lb*taub + 7*Lb^2);
C3 = pidstd(Kc3, Ti3);

% C4: PI discreto con el modelo de minimos cuadrados (P4, Ts = Tls).
% Diseno en el Control System Designer:  sisotool(P4)
%  - El cero del PI se coloca sobre el polo dominante de P4 (cancelacion)
%  - La ganancia se ajusta en el LGR para respuesta sin sobrepaso y Ms ~ 1.5
% Los valores finales se dejan fijos aqui para que el script corra solo.
pz  = pole(P4); pz = real(pz(abs(imag(pz)) < 1e-9));
pd  = max(pz);                         % polo dominante (real, cercano a 1)
Ti4 = -Tls/log(pd);                    % cero del PI sobre ese polo
Kc4 = 8;                               % ganancia escogida en sisotool
C4  = pidstd(Kc4, Ti4, 'Ts', Tls, 'IFormula', 'BackwardEuler');
% sisotool(P4, C4)   % <- descomentar para ver/ajustar el diseno

Cs = {C1, C2, C3, C4}; cn = {'C1 PID Tuner','C2 SIMC','C3 AMIGO','C4 Discreto MC'};
fprintf('\n--- Controladores PI en forma estandar  C(s) = Kp(1 + 1/(Ti s) + Td s) ---\n');
for i = 1:4
    fprintf('%-15s Kp = %.4f %%/%%   Ti = %.2f s   Td = %.2f s\n', ...
            cn{i}, Cs{i}.Kp, Cs{i}.Ti, Cs{i}.Td);
end

% Ecuaciones de los controladores (se muestran igual que P1..P4)
% Forma estandar Kp*(1 + 1/(Ti*s)) y como funcion de transferencia.
disp('C1 - PID Tuner:');
C1 %#ok<NOPTS>
C1_tf = tf(C1) %#ok<NOPTS>
disp('C2 - SIMC:');
C2 %#ok<NOPTS>
C2_tf = tf(C2) %#ok<NOPTS>
disp('C3 - AMIGO:');
C3 %#ok<NOPTS>
C3_tf = tf(C3) %#ok<NOPTS>
disp('C4 - PI discreto (modelo de minimos cuadrados):');
C4 %#ok<NOPTS>
C4_tf = tf(C4) %#ok<NOPTS>





% ---------------------------------------------
%% 7. Simulacion: servocontrol y regulatorio en el punto de operacion
% ---------------------------------------------


Tf = 10*(taub + Lb);
res = zeros(4, 9);
figure('Name','Servocontrol'); figure('Name','Regulatorio');

for i = 1:4
    C = Cs{i};
    if i == 4            % el discreto se evalua con el modelo de minimos cuadrados
        Pi = P4; ts = (0:Tls:Tf)'; wmax = pi/Tls;
    else
        Pi = P;  ts = (0:Ts:Tf)';  wmax = 100/max(taub, 1e-3);
    end
    L = C*Pi;
    Myr = feedback(L, 1);   Mur = feedback(C, Pi);
    Myd = feedback(Pi, C);  Mud = -feedback(L, 1);
    S   = feedback(1, L);

    % Servocontrol (escalon DR en r)
    ysv = DR*step(Myr, ts);  usv = DR*step(Mur, ts);
    info = stepinfo(ysv, ts, DR, 'SettlingTimeThreshold', 0.02);
    Mp  = max(0, (max(ysv) - DR)/DR*100);
    ep  = (DR - ysv(end))/DR*100;
    IAEs = trapz(ts, abs(DR - ysv));
    TVs  = sum(abs(diff(usv)));

    % Regulatorio (escalon DD en d)
    yrg = DD*step(Myd, ts);  urg = DD*step(Mud, ts);
    IAEr = trapz(ts, abs(yrg));
    epr  = yrg(end)/DD*100;
    TVr  = sum(abs(diff(urg)));

    % Robustez
    w  = logspace(log10(wmax)-6, log10(wmax), 3000);
    Ms = max(abs(squeeze(freqresp(S, w))));

    res(i,:) = [Mp, info.SettlingTime, ep, IAEs, TVs, IAEr, epr, TVr, Ms];

    % Graficas en valores absolutos (punto de operacion + desviacion)
    figure(findobj('Name','Servocontrol'));
    subplot(3,1,1); plot(ts, y0 + ysv, 'LineWidth', 1.2); hold on; grid on;
    subplot(3,1,2); plot(ts, DR - ysv, 'LineWidth', 1.2); hold on; grid on;
    subplot(3,1,3); plot(ts, u0 + usv, 'LineWidth', 1.2); hold on; grid on;
    figure(findobj('Name','Regulatorio'));
    subplot(3,1,1); plot(ts, y0 + yrg, 'LineWidth', 1.2); hold on; grid on;
    subplot(3,1,2); plot(ts, -yrg, 'LineWidth', 1.2); hold on; grid on;
    subplot(3,1,3); plot(ts, u0 + urg, 'LineWidth', 1.2); hold on; grid on;

    if any(u0 + usv > 100 | u0 + usv < 0) || any(u0 + urg > 100 | u0 + urg < 0)
        warning('%s: la senal de control sale del rango 0-100 %%.', cn{i});
    end
end

figure(findobj('Name','Servocontrol'));
subplot(3,1,1); yline(y0 + DR, 'k--'); ylabel('y [%]');
title(sprintf('Servocontrol: escalon de %g %% en r', DR)); legend([cn, {'r'}]);
subplot(3,1,2); ylabel('e [%]');
subplot(3,1,3); ylabel('u [%]'); xlabel('Tiempo [s]');
figure(findobj('Name','Regulatorio'));
subplot(3,1,1); yline(y0, 'k--'); ylabel('y [%]');
title(sprintf('Regulatorio: escalon de %g %% en d', DD)); legend([cn, {'y_0'}]);
subplot(3,1,2); ylabel('e [%]');
subplot(3,1,3); ylabel('u [%]'); xlabel('Tiempo [s]');

Tctl = array2table(res, 'RowNames', cn, 'VariableNames', ...
    {'Mp_pct','ta2_s','ep_servo_pct','IAE_servo','TVu_servo', ...
     'IAE_reg','ep_reg_pct','TVu_reg','Ms'});
disp(' '); disp('--- Desempeno simulado de los controladores ---'); disp(Tctl);

% Verificacion de especificaciones
cumple = [res(:,1) <= ESP.Mp_max, res(:,2) <= ESP.ta2_max, ...
          abs(res(:,3)) <= ESP.ep_max + 0.5, abs(res(:,7)) <= ESP.ep_max + 0.5, ...
          res(:,9) <= ESP.Ms_max];
Tesp = array2table(cumple, 'RowNames', cn, 'VariableNames', ...
    {'Mp_ok','ta2_ok','ep_servo_ok','ep_reg_ok','Ms_ok'});
disp('--- Cumplimiento de especificaciones (1 = cumple) ---'); disp(Tesp);





% ---------------------------------------------
%% 8. Laboratorio
% ---------------------------------------------


% Si el Arduino recibe la temperatura en grados C en lugar de %, la ganancia
% se escala: Kp[%/C] = Kp[%/%] * 100/(T_MAX - T_MIN). Ti y Td no cambian.

for i = 1:4
    fprintf('%-15s Kp = %.4f %%/%%  (%.4f %%/C)   Ti = %.2f s   Td = 0\n', ...
        cn{i}, Cs{i}.Kp, Cs{i}.Kp*100/(T_MAX-T_MIN), Cs{i}.Ti);
end
fprintf('(C4 es discreto: periodo de muestreo Ts = %.0f s, integral por Euler hacia atras)\n', Tls);
fprintf('Punto de operacion: u0 = %.2f %%,  y0 = %.2f %% (%.2f C)\n', ...
        u0, y0, y0*(T_MAX-T_MIN)/100 + T_MIN);





% ---------------------------------------------
%% Funciones locales utilizadas
% ---------------------------------------------


function tc = cruce(t, yn, p)
% Primer instante (interpolado) en que la respuesta normalizada alcanza p
    k = find(yn >= p, 1);
    if isempty(k), error('La respuesta no alcanza el %.1f %% del cambio.', 100*p); end
    if k == 1, tc = t(1); return; end
    tc = t(k-1) + (p - yn(k-1))*(t(k) - t(k-1))/(yn(k) - yn(k-1));
end

function [K, tau, L] = pomtm(G)
% Extrae K, tau y L de un modelo K e^{-Ls}/(tau s + 1)
    [num, den] = tfdata(G, 'v');
    num = num(find(num, 1):end); den = den(find(den, 1):end);
    K   = num(end)/den(end);
    tau = den(1)/den(end);
    L   = G.InputDelay + G.OutputDelay + G.IODelay;
end