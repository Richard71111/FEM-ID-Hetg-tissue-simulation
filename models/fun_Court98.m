function [G, Iion, Ivec, dX, Iall] = fun_Court98(t, X, p, S)
%FUN_COURT98 Advance the Courtemanche 1998 ionic model by one time step.
% Inputs: time t, membrane voltage and states X, parameter structure p,
% and extracellular Na/K/Ca concentrations S.
% Outputs: updated states G, total ionic current Iion, species currents
% Ivec, state derivatives dX, and individual current components Iall.
if ~isfield(p,'dt')
    p.dt = 1;
end
dt = p.dt;
Nm = p.Npatches;


V = X(1:Nm);
G = X(Nm+1:end);
dX = zeros(length(V)+length(G),1);


Na_i = X(Nm+1:2*Nm);
m = X(2*Nm+1:3*Nm);
h = X(3*Nm+1:4*Nm);
j = X(4*Nm+1:5*Nm);
K_i = X(5*Nm+1:6*Nm);
oa = X(6*Nm+1:7*Nm);
oi = X(7*Nm+1:8*Nm);
ua = X(8*Nm+1:9*Nm);
ui = X(9*Nm+1:10*Nm);
xr = X(10*Nm+1:11*Nm);
xs = X(11*Nm+1:12*Nm);
Ca_i = X(12*Nm+1:13*Nm);
d = X(13*Nm+1:14*Nm);
f = X(14*Nm+1:15*Nm);
f_Ca = X(15*Nm+1:16*Nm);
Ca_rel = X(16*Nm+1:17*Nm);
u = X(17*Nm+1:18*Nm);
v = X(18*Nm+1:19*Nm);
w = X(19*Nm+1:20*Nm);
Ca_up = X(20*Nm+1:21*Nm);

if nargin < 4
    Na_o = 140;
    K_o = 5.4;
    Ca_o = 1.8;
else
    Na_o = S(1:Nm);   % extracellular Na, mM
    K_o = S(Nm+1:2*Nm); % extracellular Ko
    Ca_o = S(2*Nm+1:3*Nm); % extracellular Ca
end

% constants
R = 8.3143;
T = 310;
F = 96.4867;
Cm = 100; % pF
g_Na = 7.8;
g_K1 = 0.09;
K_Q10 = 3;
g_to = 0.1652;
g_Kr = 0.029411765;
g_Ks = 0.12941176;
g_Ca_L = 0.12375;
Km_Na_i = 10;
Km_K_o = 1.5;
i_NaK_max = 0.59933874;
g_B_Na = 0.0006744375;
g_B_Ca = 0.001131;
g_B_K = 0;
I_NaCa_max = 1600;
K_mNa = 87.5;
K_mCa = 1.38;
K_sat = 0.1;
gamma = 0.35;
i_CaP_max = 0.275;
K_rel = 30;
tau_tr = 180;
I_up_max = 0.005;
K_up = 0.00092;
Ca_up_max = 15;
CMDN_max = 0.05;
TRPN_max = 0.07;
CSQN_max = 10;
Km_CMDN = 0.00238;
Km_TRPN = 0.0005;
Km_CSQN = 0.8;
V_cell = 20100;
V_i =  V_cell.*0.680000;
tau_f_Ca = 2.00000;
sigma =  (1.00000./7.00000).*(exp(Na_o./67.3000) - 1.00000);
tau_u = 8.00000;
V_rel =  0.00480000.*V_cell;
V_up =  0.0552000.*V_cell;

f_Ca_infinity = power(1.00000+Ca_i./0.000350000,  - 1.00000);
diff_f_Ca = (f_Ca_infinity - f_Ca)./tau_f_Ca;

d_infinity = power(1.00000+exp((V+10.0000)./ - 8.00000),  - 1.00000);
tau_d = (1.00000 - exp((V+10.0000)./ - 6.24000))./( 0.0350000.*(V+10.0000).*(1.00000+exp((V+10.0000)./ - 6.24000)));
diff_d = (d_infinity - d)./tau_d;

f_infinity = exp( - (V+28.0000)./6.90000)./(1.00000+exp( - (V+28.0000)./6.90000));
tau_f =  9.00000.*power( 0.0197000.*exp(  - power(0.0337000, 2.00000).*power(V+10.0000, 2.00000))+0.0200000,  - 1.00000);
diff_f = (f_infinity - f)./tau_f;

tau_w =  ( 6.00000.*(1.00000 - exp( - (V - 7.90000)./5.00000)))./( (1.00000+ 0.300000.*exp( - (V - 7.90000)./5.00000)).*1.00000.*(V - 7.90000));
w_infinity = 1.00000 - power(1.00000+exp( - (V - 40.0000)./17.0000),  - 1.00000);
diff_w = (w_infinity - w)./tau_w;

alpha_m = ( 0.320000.*(V+47.1300))./(1.00000 - exp(  - 0.100000.*(V+47.1300)));
beta_m =  0.0800000.*exp( - V./11.0000);
m_inf = alpha_m./(alpha_m+beta_m);
tau_m = 1.00000./(alpha_m+beta_m);
diff_m = (m_inf - m)./tau_m;

a=1-1./(1+exp(-(V+40)/0.24));
alpha_h =  a.*(0.135000.*exp((V+80.0000)./ - 6.80000));
beta_h =   a.*(3.56000.*exp( 0.0790000.*V)+ 310000..*exp( 0.350000.*V)) + (1-a).*(1.00000./( 0.130000.*(1.00000+exp((V+10.6600)./ - 11.1000))));
alpha_j =  a.*(( (  - 127140..*exp( 0.244400.*V) -  3.47400e-05.*exp(  - 0.0439100.*V)).*(V+37.7800))./(1.00000+exp( 0.311000.*(V+79.2300))));
beta_j =   a.*(( 0.121200.*exp(  - 0.0105200.*V))./(1.00000+exp(  - 0.137800.*(V+40.1400)))) + (1-a).*((0.300000.*exp(  - 2.53500e-07.*V))./(1.00000+exp(- 0.100000.*(V+32.0000))));

h_inf = alpha_h./(alpha_h+beta_h);
tau_h = 1.00000./(alpha_h+beta_h);
diff_h = (h_inf - h)./tau_h;

j_inf = alpha_j./(alpha_j+beta_j);
tau_j = 1.00000./(alpha_j+beta_j);
diff_j = (j_inf - j)./tau_j;

alpha_oa =  0.650000.*power(exp((V -  - 10.0000)./ - 8.50000)+exp(((V -  - 10.0000) - 40.0000)./ - 59.0000),  - 1.00000);
beta_oa =  0.650000.*power(2.50000+exp(((V -  - 10.0000)+72.0000)./17.0000),  - 1.00000);
tau_oa = power(alpha_oa+beta_oa,  - 1.00000)./K_Q10;
oa_infinity = power(1.00000+exp(((V -  - 10.0000)+10.4700)./ - 17.5400),  - 1.00000);
diff_oa = (oa_infinity - oa)./tau_oa;

alpha_oi = power(18.5300+ 1.00000.*exp(((V -  - 10.0000)+103.700)./10.9500),  - 1.00000);
beta_oi = power(35.5600+ 1.00000.*exp(((V -  - 10.0000) - 8.74000)./ - 7.44000),  - 1.00000);
tau_oi = power(alpha_oi+beta_oi,  - 1.00000)./K_Q10;
oi_infinity = power(1.00000+exp(((V -  - 10.0000)+33.1000)./5.30000),  - 1.00000);
diff_oi = (oi_infinity - oi)./tau_oi;

alpha_ua =  0.650000.*power(exp((V -  - 10.0000)./ - 8.50000)+exp(((V -  - 10.0000) - 40.0000)./ - 59.0000),  - 1.00000);
beta_ua =  0.650000.*power(2.50000+exp(((V -  - 10.0000)+72.0000)./17.0000),  - 1.00000);
tau_ua = power(alpha_ua+beta_ua,  - 1.00000)./K_Q10;
ua_infinity = power(1.00000+exp(((V -  - 10.0000)+20.3000)./ - 9.60000),  - 1.00000);
diff_ua = (ua_infinity - ua)./tau_ua;

alpha_ui = power(21.0000+ 1.00000.*exp(((V -  - 10.0000) - 195.000)./ - 28.0000),  - 1.00000);
beta_ui = 1.00000./exp(((V -  - 10.0000) - 168.000)./ - 16.0000);
tau_ui = power(alpha_ui+beta_ui,  - 1.00000)./K_Q10;
ui_infinity = power(1.00000+exp(((V -  - 10.0000) - 109.450)./27.4800),  - 1.00000);
diff_ui = (ui_infinity - ui)./tau_ui;

alpha_xr =  ( 0.000300000.*(V+14.1000))./(1.00000 - exp((V+14.1000)./ - 5.00000));
beta_xr =  ( 7.38980e-05.*(V - 3.33280))./(exp((V - 3.33280)./5.12370) - 1.00000);
tau_xr = power(alpha_xr+beta_xr,  - 1.00000);
xr_infinity = power(1.00000+exp((V+14.1000)./ - 6.50000),  - 1.00000);
diff_xr = (xr_infinity - xr)./tau_xr;

alpha_xs =  ( 4.00000e-05.*(V - 19.9000))./(1.00000 - exp((V - 19.9000)./ - 17.0000));
beta_xs =  ( 3.50000e-05.*(V - 19.9000))./(exp((V - 19.9000)./9.00000) - 1.00000);
tau_xs =  0.500000.*power(alpha_xs+beta_xs,  - 1.00000);
xs_infinity = power(1.00000+exp((V - 19.9000)./ - 12.7000),  - 0.500000);
diff_xs = (xs_infinity - xs)./tau_xs;

E_K =  (( R.*T)./F).*log(K_o./K_i);
i_K1 = p.f_I(p.iik1,:)'.*( g_K1.*(V - E_K))./(1.00000+exp( 0.0700000.*(V+80.0000)));
i_to =  p.f_I(p.iito,:)'.*g_to.*power(oa, 3.00000).*oi.*(V - E_K);
g_Kur = 0.00500000+0.0500000./(1.00000+exp((V - 15.0000)./ - 13.0000));
i_Kur =  p.f_I(p.iikur,:)'.*g_Kur.*power(ua, 3.00000).*ui.*(V - E_K);
i_Kr = p.f_I(p.iikr,:)'.*( g_Kr.*xr.*(V - E_K))./(1.00000+exp((V+15.0000)./22.4000));
i_Ks =  p.f_I(p.iiks,:)'.*g_Ks.*power(xs, 2.00000).*(V - E_K);
f_NaK = power(1.00000+ 0.124500.*exp((  - 0.100000.*F.*V)./( R.*T))+ 0.0365000.*sigma.*exp((  - F.*V)./( R.*T)),  - 1.00000);
i_NaK = p.f_I(p.iinak,:)'.*( (( i_NaK_max.*f_NaK.*1.00000)./(1.00000+power(Km_Na_i./Na_i, 1.50000))).*K_o)./(K_o+Km_K_o);
i_B_K =  p.f_I(p.iibk,:)'.*g_B_K.*(V - E_K);
diff_K_i = ( 2.00000.*i_NaK - (i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_K))./( V_i.*F/Cm);

E_Na =  (( R.*T)./F).*log(Na_o./Na_i);
i_Na =  p.f_I(p.iina,:)'.*g_Na.*power(m, 3.00000).*h.*j.*(V - E_Na);
i_NaCa = p.f_I(p.iinaca,:)'.*( I_NaCa_max.*( exp(( gamma.*F.*V)./( R.*T)).*power(Na_i, 3.00000).*Ca_o -  exp(( (gamma - 1.00000).*F.*V)./( R.*T)).*power(Na_o, 3.00000).*Ca_i))./( (power(K_mNa, 3.00000)+power(Na_o, 3.00000)).*(K_mCa+Ca_o).*(1.00000+ K_sat.*exp(( (gamma - 1.00000).*V.*F)./( R.*T))));
i_B_Na =  p.f_I(p.iibna,:)'.*g_B_Na.*(V - E_Na);
diff_Na_i = (  - 3.00000.*i_NaK - ( 3.00000.*i_NaCa+i_B_Na+i_Na))./( V_i.*F/Cm);

% i_st = piecewise({t>=stim_start&t<=stim_end&(t - stim_start) -  floor((t - stim_start)./stim_period).*stim_period<=stim_duration, stim_amplitude }, 0.00000);

i_st = p.stim_amp*(mod(t,p.BCL)<p.stim_dur).*p.indstim;

% i_st = stim_amplitude*(mod(t,stim_period)>=stim_start & mod(t,stim_period)<(stim_start+stim_duration)).*p.indstim;

i_Ca_L =  p.f_I(p.iical,:)'.*g_Ca_L.*d.*f.*f_Ca.*(V - 65.0000);
i_CaP = p.f_I(p.iicap,:)'.*( i_CaP_max.*Ca_i)./(0.000500000+Ca_i);
E_Ca =  (( R.*T)./( 2.00000.*F)).*log(Ca_o./Ca_i);
i_B_Ca =  p.f_I(p.iibca,:)'.*g_B_Ca.*(V - E_Ca);

diff_V =  - (i_Na+i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_Na+i_B_K+i_B_Ca+i_NaK+i_CaP+i_NaCa+i_Ca_L+i_st);

% currents in pA/pF * uF = uA/uF * uF = uA
Iion =   p.Ctot*(i_Na+i_K1+i_to+i_Kur+i_Kr+i_Ks+i_B_Na+i_B_K+i_B_Ca+i_NaK+i_CaP+i_NaCa+i_Ca_L+i_st);

Ivec = nan(3*Nm,1);
Ivec(1:Nm) = p.Ctot*(i_Na + i_B_Na + 3*i_NaCa + 3*i_NaK);  % Na
Ivec(Nm+1:2*Nm) = p.Ctot*(i_K1 + i_to + i_Kur + i_Kr + i_Ks + i_B_K - 2*i_NaK); % K
Ivec(2*Nm+1:3*Nm) = p.Ctot*(i_B_Ca + i_CaP + i_Ca_L - 2*i_NaCa);

% all ionic currents
[Ncurrents,~] = size(p.f_I);
Iall = nan(Nm*Ncurrents,1);

Iall(p.iina:Ncurrents:end) = p.Ctot*i_Na;
Iall(p.iik1:Ncurrents:end) = p.Ctot*i_K1;
Iall(p.iito:Ncurrents:end) = p.Ctot*i_to;
Iall(p.iikur:Ncurrents:end) = p.Ctot*i_Kur;
Iall(p.iikr:Ncurrents:end) = p.Ctot*i_Kr;
Iall(p.iiks:Ncurrents:end) = p.Ctot*i_Ks;
Iall(p.iibna:Ncurrents:end) = p.Ctot*i_B_Na;
Iall(p.iibk:Ncurrents:end) = p.Ctot*i_B_K;
Iall(p.iibca:Ncurrents:end) = p.Ctot*i_B_Ca;
Iall(p.iinak:Ncurrents:end) = p.Ctot*i_NaK;
Iall(p.iicap:Ncurrents:end) = p.Ctot*i_CaP;
Iall(p.iinaca:Ncurrents:end) = p.Ctot*i_NaCa;
Iall(p.iical:Ncurrents:end) = p.Ctot*i_Ca_L;



i_rel =  K_rel.*power(u, 2.00000).*v.*w.*(Ca_rel - Ca_i);
i_tr = (Ca_up - Ca_rel)./tau_tr;
diff_Ca_rel =  (i_tr - i_rel).*power(1.00000+( CSQN_max.*Km_CSQN)./power(Ca_rel+Km_CSQN, 2.00000),  - 1.00000);
Fn =  1000.00.*( 1.00000e-15.*V_rel.*i_rel -  (1.00000e-15./( 2.00000.*F)).*( 0.500000.*Cm*i_Ca_L -  0.200000.*Cm*i_NaCa));
u_infinity = power(1.00000+exp( - (Fn - 3.41750e-13)./1.36700e-15),  - 1.00000);
diff_u = (u_infinity - u)./tau_u;

tau_v = 1.91000+ 2.09000.*power(1.00000+exp( - (Fn - 3.41750e-13)./1.36700e-15),  - 1.00000);
v_infinity = 1.00000 - power(1.00000+exp( - (Fn - 6.83500e-14)./1.36700e-15),  - 1.00000);
diff_v = (v_infinity - v)./tau_v;

i_up = I_up_max./(1.00000+K_up./Ca_i);
i_up_leak = ( I_up_max.*Ca_up)./Ca_up_max;
diff_Ca_up = i_up - (i_up_leak+( i_tr.*V_rel)./V_up);
B1 = Cm*( 2.00000.*i_NaCa - (i_CaP+i_Ca_L+i_B_Ca))./( 2.00000.*V_i.*F)+( V_up.*(i_up_leak - i_up)+ i_rel.*V_rel)./V_i;
B2 = 1.00000+( TRPN_max.*Km_TRPN)./power(Ca_i+Km_TRPN, 2.00000)+( CMDN_max.*Km_CMDN)./power(Ca_i+Km_CMDN, 2.00000);
diff_Ca_i = B1./B2;

G(1:Nm) = Na_i + dt*diff_Na_i;
G(Nm+1:2*Nm) = m_inf - (m_inf - m).*exp(-dt./tau_m);
G(2*Nm+1:3*Nm) = h_inf - (h_inf - h).*exp(-dt./tau_h);
G(3*Nm+1:4*Nm) = j_inf - (j_inf - j).*exp(-dt./tau_j);
G(4*Nm+1:5*Nm) = K_i + dt*diff_K_i;
G(5*Nm+1:6*Nm) = oa_infinity - (oa_infinity - oa).*exp(-dt./tau_oa);
G(6*Nm+1:7*Nm) = oi_infinity - (oi_infinity - oi).*exp(-dt./tau_oi);
G(7*Nm+1:8*Nm) = ua_infinity - (ua_infinity - ua).*exp(-dt./tau_ua);
G(8*Nm+1:9*Nm) = ui_infinity - (ui_infinity - ui).*exp(-dt./tau_ui);
G(9*Nm+1:10*Nm) = xr_infinity - (xr_infinity - xr).*exp(-dt./tau_xr);
G(10*Nm+1:11*Nm) = xs_infinity - (xs_infinity - xs).*exp(-dt./tau_xs);
G(11*Nm+1:12*Nm) = Ca_i + dt*diff_Ca_i;
G(12*Nm+1:13*Nm) = d_infinity - (d_infinity - d).*exp(-dt./tau_d);
G(13*Nm+1:14*Nm) = f_infinity - (f_infinity - f).*exp(-dt./tau_f);
G(14*Nm+1:15*Nm) = f_Ca_infinity - (f_Ca_infinity - f_Ca).*exp(-dt./tau_f_Ca);
G(15*Nm+1:16*Nm) = Ca_rel + dt*diff_Ca_rel;
G(16*Nm+1:17*Nm) = u_infinity - (u_infinity - u).*exp(-dt./tau_u);
G(17*Nm+1:18*Nm) = v_infinity - (v_infinity - v).*exp(-dt./tau_v);
G(18*Nm+1:19*Nm) = w_infinity - (w_infinity - w).*exp(-dt./tau_w);
G(19*Nm+1:20*Nm) = Ca_up + dt*diff_Ca_up;


if Nm == 1
    dX(1) = diff_V;
    dX(2) = diff_Na_i;
    dX(3) = diff_m;
    dX(4) = diff_h;
    dX(5) = diff_j;
    dX(6) = diff_K_i;
    dX(7) = diff_oa;
    dX(8) = diff_oi;
    dX(9) = diff_ua;
    dX(10) = diff_ui;
    dX(11) = diff_xr;
    dX(12) = diff_xs;
    dX(13) = diff_Ca_i;
    dX(14) = diff_d;
    dX(15) = diff_f;
    dX(16) = diff_f_Ca;
    dX(17) = diff_Ca_rel;
    dX(18) = diff_u;
    dX(19) = diff_v;
    dX(20) = diff_w;
    dX(21) = diff_Ca_up;
end
end

% function x = piecewise(cases, default)
%     set = [0];
%     for i = 1:2:length(cases)
%         if (length(cases{i+1}) == 1)
%             x(cases{i} & ~set,:) = cases{i+1};
%         else
%             x(cases{i} & ~set,:) = cases{i+1}(cases{i} & ~set);
%         end
%         set = set | cases{i};
%         if(set), break, end
%     end
%     if (length(default) == 1)
%         x(~set,:) = default;
%     else
%         x(~set,:) = default(~set);
%     end
% end
