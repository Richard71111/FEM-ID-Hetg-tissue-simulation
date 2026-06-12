function model = setup_ionic_model(cfg)
%SETUP_IONIC_MODEL Select LR1, ORd11, or Court98 and set initial values.
% Input: cfg structure from default_config.
% Output: model structure with parameters, initial state, and current data.

Aax = 2 * pi * cfg.r * cfg.L;  % Axial membrane area, um^2.
Adisc = pi * cfg.r^2;  % Area of one intercalated disc, um^2.
Atot = Aax + 2 * Adisc;  % Total cell membrane area, um^2.
Ctot = Atot * cfg.Cm;  % Total cell capacitance, uF.

switch lower(cfg.model)
    case "lr1"
        [p, x0] = Initial_LR1(Atot);

        p.iina = 1;
        p.iisi = 2;
        p.iik = 3;
        p.iik1 = 4;
        p.iikp = 5;
        p.iib = 6;
        p.mLR1 = 1;

        Ncurrents = 6;
        Nstate = 7;
        ionic_fun = @fun_LR1;

        loc_vec = zeros(1, Ncurrents);
        loc_vec(p.iina) = cfg.locINa;
        loc_vec(p.iik1) = cfg.locIK1;
        loc_vec(p.iisi) = cfg.locICa;

        if isempty(cfg.stim_dur)
            p.stim_dur = 1;
        else
            p.stim_dur = cfg.stim_dur;
        end
        if isempty(cfg.stim_amp)
            p.stim_amp = 0.5 * 80e-8 * Atot;
        else
            p.stim_amp = cfg.stim_amp;
        end

    case "ord11"
        p.iina = 1;
        p.iinal = 2;
        p.iito = 3;
        p.iical = 4;
        p.iikr = 5;
        p.iiks = 6;
        p.iik1 = 7;
        p.iinaca_i = 8;
        p.iinaca_ss = 9;
        p.iinak = 10;
        p.iikb = 11;
        p.iinab = 12;
        p.iicab = 13;
        p.iipca = 14;

        p.fSERCA = 1;
        p.fRyR = 1;
        p.ftauhL = 1;
        p.fCaMKa = 1;
        p.fIleak = 1;
        p.fJrel = 1;

        x0 = Initial_ORd11();
        Ncurrents = 14;
        Nstate = 40;
        ionic_fun = @fun_ORd11;

        loc_vec = zeros(1, Ncurrents);
        loc_vec(p.iina) = cfg.locINa;
        loc_vec(p.iik1) = cfg.locIK1;
        loc_vec(p.iical) = cfg.locICa;
        loc_vec(p.iinak) = cfg.locINaK;
        loc_vec(p.iito) = 2 * Adisc / Atot;
        loc_vec(p.iikr) = 2 * Adisc / Atot;
        loc_vec(p.iiks) = 2 * Adisc / Atot;

        if isempty(cfg.stim_dur)
            p.stim_dur = 2;
        else
            p.stim_dur = cfg.stim_dur;
        end
        if isempty(cfg.stim_amp)
            p.stim_amp = 50;
        else
            p.stim_amp = cfg.stim_amp;
        end

    case "court98"
        p.iina = 1;
        p.iik1 = 2;
        p.iito = 3;
        p.iikur = 4;
        p.iikr = 5;
        p.iiks = 6;
        p.iibna = 7;
        p.iibk = 8;
        p.iibca = 9;
        p.iinak = 10;
        p.iicap = 11;
        p.iinaca = 12;
        p.iical = 13;

        x0 = Initial_Court98();
        Ncurrents = 13;
        Nstate = 20;
        ionic_fun = @fun_Court98;

        loc_vec = zeros(1, Ncurrents);
        loc_vec(p.iina) = cfg.locINa;
        loc_vec(p.iik1) = cfg.locIK1;
        loc_vec(p.iical) = cfg.locICa;
        loc_vec(p.iinak) = cfg.locINaK;
        loc_vec(p.iito) = 2 * Adisc / Atot;
        loc_vec(p.iikr) = 2 * Adisc / Atot;
        loc_vec(p.iiks) = 2 * Adisc / Atot;

        if isempty(cfg.stim_dur)
            p.stim_dur = 1;
        else
            p.stim_dur = cfg.stim_dur;
        end
        if isempty(cfg.stim_amp)
            p.stim_amp = -60;
        else
            p.stim_amp = cfg.stim_amp;
        end

    otherwise
        error("Unknown ionic model: %s", cfg.model);
end

p.L = cfg.L;
p.r = cfg.r;
p.Ctot = Ctot;
p.Na_o = cfg.Na_b;
p.K_o = cfg.K_b;
p.Ca_o = cfg.Ca_b;

model.p = p;
model.x0 = x0(:)';
model.Nstate = Nstate;
model.Ncurrents = Ncurrents;
model.ionic_fun = ionic_fun;
model.loc_vec = loc_vec;
model.scaleI = ones(1, Ncurrents);
end
