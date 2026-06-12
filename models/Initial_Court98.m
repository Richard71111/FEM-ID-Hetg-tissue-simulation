function initial_state = Initial_Court98
%INITIAL_COURT98 Return the Court98 resting voltage and state variables.
% Input: none.
% Output: initial_state vector containing voltage followed by model states.

% initial conditions
V = -81.18;
Na_i = 1.117e+01;
m = 2.908e-3;
h = 9.649e-1;
j = 9.775e-1;
K_i = 1.39e+02;
oa = 3.043e-2;
oi = 9.992e-1;
ua = 4.966e-3;
ui = 9.986e-1;
xr = 3.296e-5;
xs = 1.869e-2;
Ca_i = 1.013e-4;
d = 1.367e-4;
f = 9.996e-1;
f_Ca = 7.755e-1;
Ca_rel = 1.488;
u = 2.35e-112;
v = 1;
w = 0.9992;
Ca_up = 1.488;

initial_state(1) = V;
initial_state(2) = Na_i;
initial_state(3) = m;
initial_state(4) = h;
initial_state(5) = j;
initial_state(6) = K_i;
initial_state(7) = oa;
initial_state(8) = oi;
initial_state(9) = ua;
initial_state(10) = ui;
initial_state(11) = xr;
initial_state(12) = xs;
initial_state(13) = Ca_i;
initial_state(14) = d;
initial_state(15) = f;
initial_state(16) = f_Ca;
initial_state(17) = Ca_rel;
initial_state(18) = u;
initial_state(19) = v;
initial_state(20) = w;
initial_state(21) = Ca_up;
