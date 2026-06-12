function [G, Iion, Ivec, dX, Iall] = fun_LR1(t,X,p,S)
%FUN_LR1 Advance the Luo-Rudy 1991 ionic model by one time step.
% Inputs: time t, membrane voltage and states X, parameters p, and
% extracellular Na/K/Ca concentrations S.
% Outputs: updated states, total current, species currents, derivatives,
% and individual current components.
if ~isfield(p,'dt')
    p.dt = 1;
end
dt = p.dt;
Nm = p.Npatches;

Vm = X(1:Nm);
G = X(Nm+1:end);
dX = zeros(length(Vm)+length(G),1);

M = X(Nm+1:2*Nm);
H = X(2*Nm+1:3*Nm);
J = X(3*Nm+1:4*Nm);
D = X(4*Nm+1:5*Nm);
F = X(5*Nm+1:6*Nm);
Xi = X(6*Nm+1:7*Nm);
Ca_i = X(7*Nm+1:8*Nm);

Nao = S(1:Nm);   % extracellular Na
Ko = S(Nm+1:2*Nm); % extracellular Ko
Cao = S(2*Nm+1:3*Nm); % extracellular Ca


[INa,am,bm,ah,bh,aj,bj] = comp_INa(Vm,M,H,J,Nao, p);
[Isi,ad,bd,af,bf] = comp_Isi91(Vm,D,F,Ca_i,Cao,p);
[IK,ax,bx] = comp_IK(Vm,Xi,Nao,Ko,p);
IK1 = comp_IK1(Vm,Ko,p);
IKp = comp_IKp(Vm,Ko,p);
Ib = p.f_I(p.iib,:)'.*p.Ib_max.*(Vm+59.87);

Ivec = nan(Nm*3,1);    % ionic currents
Ivec(1:Nm) = INa;      % Na+
Ivec(Nm+1:2*Nm) = IK + IK1 + IKp; % K+
Ivec(2*Nm+1:3*Nm) = Isi;      % Ca2+


taum = 1./(am+bm); minf = am.*taum;
tauh = 1./(ah+bh); hinf = ah.*tauh;
tauj = 1./(aj+bj); jinf = aj.*tauj;
taud = 1./(ad+bd); dinf = ad.*taud;
tauf = 1./(af+bf); finf = af.*tauf;
taux = 1./(ax+bx); xinf = ax.*taux;

if p.mLR1   % modified LR1, as in Fenton-Karma 1998
    taud = taud/2;
    tauf = tauf/2;
end

% Rush-Larson for HH gates
G(1:Nm) = minf - (minf - M).*exp(-dt./taum);       % M
G(Nm+1:2*Nm) = hinf - (hinf - H).*exp(-dt./tauh);         % H
G(2*Nm+1:3*Nm) = jinf - (jinf - J).*exp(-dt./tauj);       % J
G(3*Nm+1:4*Nm) = dinf - (dinf - D).*exp(-dt./taud);       % D
G(4*Nm+1:5*Nm) = finf - (finf - F).*exp(-dt./tauf);       % F
G(5*Nm+1:6*Nm) = xinf - (xinf - Xi).*exp(-dt./taux);  

G(6*Nm+1:7*Nm) = Ca_i + dt*(-1e-4*Isi + 0.07*(1e-4-Ca_i)); % Ca_i

Istim = p.stim_amp*(mod(t,p.BCL)<p.stim_dur).*p.indstim;
Iion = -(Istim - INa - Isi - IK - IK1 - IKp - Ib);  % uA

% all ionic currents
[Ncurrents,~] = size(p.f_I);
Iall = nan(Nm*Ncurrents,1);

Iall(p.iina:Ncurrents:end) = INa;
Iall(p.iisi:Ncurrents:end) = Isi;
Iall(p.iik:Ncurrents:end) = IK;
Iall(p.iik1:Ncurrents:end) = IK1;
Iall(p.iikp:Ncurrents:end) = IKp;
Iall(p.iib:Ncurrents:end) = Ib;

if Nm == 1
    dX(1:Nm) = -Iion/p.Ctot;
    dX(Nm+1:2*Nm) = (minf - M)./taum;       % M
    dX(2*Nm+1:3*Nm) = (hinf - H)./tauh;         % H
    dX(3*Nm+1:4*Nm) = (jinf - J)./tauj;       % J
    dX(4*Nm+1:5*Nm) = (dinf - D)./taud;       % D
    dX(5*Nm+1:6*Nm) = (finf - F)./tauf;       % F
    dX(6*Nm+1:7*Nm) = (xinf - Xi)./taux;
    dX(7*Nm+1:8*Nm) = (-1e-4*Isi + 0.07*(1e-4-Ca_i)); % Ca_i
end
end

function [INa,alpham,betam,alphah,betah,alphaj,betaj] = comp_INa(V,m,H,J,Nao,data)
% INa    the fast sodium current in mammalian ventricular cells
%

ENa =data.RTF*log(Nao./data.Na_i);       % Nernst potential of Na, mV
INa = data.INa_max*data.f_I(data.iina,:)'.*m.*m.*m.*H.*J.*(V-ENa);

alpham = am(V); betam = bm(V);
alphah = ah(V); betah = bh(V);
alphaj = aj(V); betaj = bj(V);
end


function ah = ah(V)
a=1-1./(1+exp(-(V+40)/0.24));
ah= a.*0.135.*exp((80+V)./(-6.8));
%  ah= a.*0.135.*exp((80+V)./(-15.8));

end

function bh = bh(V)
a=1-1./(1+exp(-(V+40)/0.24));
bh= (1-a)./(0.13*(1+exp((V+10.66)/(-11.1)))) +(a).*(3.56*exp(0.079*V)+3.1*1e5*exp(0.35*V));
% bh= (1-a)./(0.13*(1+exp((V+10.66)/(-11.1)))) +(a).*(.1*3.56*exp(0.079*V)+3.1*1e5*exp(0.35*V));

end

function aj = aj(V)
a=1-1./(1+exp(-(V+40)/0.24));
aj =  a.*(-1.2714e5*exp(0.2444*V)-3.474e-5*exp(-0.04391*V)).*(V+37.78)./(1+exp(0.311*(V+79.23)));
end

function bj = bj(V)
a=1-1./(1+exp(-(V+40)/0.24));
bj= (1-a).*(0.3*exp(-2.535e-7*V)./(1+exp(-0.1*(V+32))))+(a).*(0.1212*exp(-0.01052*V)./(1+exp(-0.1378*(V+40.14))));
end

function am = am(V)
am = 0.32*(V+47.13)./(1-exp(-0.1*(V+47.13)));
end

function bm = bm(V)
bm = 0.08*exp(-V/11);
end

function [IK,alphax,betax] = comp_IK(V,x,Nao, Ko,data)
% IK    Time-dependent potassium current
Xi = xi(V);
EK = data.RTF.*log((Ko+data.PNa_K*Nao)/(data.K_i+data.PNa_K*data.Na_i));        % mV
sq=sqrt(Ko/5.4); 
GK = data.IK_max*sq;         % mS/um^2

IK = data.f_I(data.iik,:)'.*GK.*Xi.*x.*(V-EK);
alphax = ax(V); betax = bx(V);
end

function Xi = xi(V)
Xi = (V<=-100)+ ...
    (V > -100).*(2.837*(exp(0.04*(V+77))-1)./((V+77).*exp(0.04*(V+35))));
end

function ax = ax(V)
ax = 0.0005*exp(0.083*(V+50))./(1+exp(0.057*(V+50)));
end

function bx = bx(V)
bx = 0.0013*exp(-0.06*(V+20))./(1+exp(-0.04*(V+20)));
end 

function IK1 = comp_IK1(V,Ko,data)
% IK1    Time-independent potassium current

EK1 = data.RTF*log(Ko./data.K_i);
ak1 = 1.02./(1+exp(0.2385*(V-EK1-59.215)));
bk1 = (0.49124*exp(0.08032*(V-EK1+5.476))+exp(0.06175*(V-EK1-594.31)))./...
    (1+exp(-0.5143*(V-EK1+4.753)));
sq=sqrt(Ko/5.4); 

gK1 =data.IK1_max*sq.*ak1./(ak1+bk1);
IK1 = data.f_I(data.iik1,:)'.*gK1.*(V-EK1);
end

function [Is,alphad,betad,alphaf,betaf] = comp_Isi91(V,d,f,Ca_i,Cao,data)
% Isi    the slow inward current
% Esi = 7.7 - 13.0287*log(Ca_i);
Esi = data.RTF/2*log(Cao./Ca_i);
Is =  data.f_I(data.iisi,:)'.*data.Isi_max.*d.*f.*(V-Esi);

alphad = ad(V); betad = bd(V);
alphaf = af(V); betaf = bf(V);
end

function ad = ad(V)
ad = 0.095*exp(-0.01*(V-5))./(1+exp(-0.072*(V-5)));
end

function bd = bd(V)
bd = 0.07*exp(-0.017*(V+44))./(1+exp(0.05*(V+44)));
end

function af = af(V)
af = 0.012*exp(-0.008*(V+28))./(1+exp(0.15*(V+28)));
end

function bf = bf(V)
bf = 0.0065*exp(-0.02*(V+30))./(1+exp(-0.2*(V+30)));
end

function IKp = comp_IKp(V,Ko,data)
EK1 = data.RTF*log(Ko./data.K_i);
IKp = data.f_I(data.iikp,:)'.*data.IKp_max.*(V-EK1)./(1+exp((7.488-V)/5.98));
end

