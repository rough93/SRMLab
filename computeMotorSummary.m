function s = computeMotorSummary(out, cfg)
%COMPUTEMOTORSUMMARY Compute OpenMotor-like readout metrics (stack-safe)

g0 = 9.80665;

% Pull histories
t  = out.t(:);
Pc = out.Pc(:);
F  = out.F(:);

Pa = cfg.noz.Pa;
At = cfg.noz.At;

% Histories required for flux calculations
mdot_noz = out.mdot_noz(:);
Ap_eff   = out.Ap(:);

% Optional per-segment port areas (stack)
if isfield(out,'Ap_seg') && ~isempty(out.Ap_seg)
    Ap_eff_seg = out.Ap_seg;
else
    Ap_eff_seg = [];
end

% Ensure all referenced histories have a common length
N = min([numel(t), numel(Pc), numel(F), numel(mdot_noz), numel(Ap_eff)]);
t        = t(1:N);
Pc       = Pc(1:N);
F        = F(1:N);
mdot_noz = mdot_noz(1:N);
Ap_eff   = Ap_eff(1:N);

if ~isempty(Ap_eff_seg)
    Ap_eff_seg = Ap_eff_seg(1:min(size(Ap_eff_seg,1), N), :);
    if size(Ap_eff_seg,1) < N
        % Pad with NaNs so indexing stays safe
        Ap_eff_seg(end+1:N, :) = NaN;
    end
end

% Burn time definition: from first time Pc exceeds 5% of peak to last time >5% of peak
Pc_pk = max(Pc);
th = 0.05 * Pc_pk;
i0 = find(Pc > th, 1, 'first');
i1 = find(Pc > th, 1, 'last');
if isempty(i0) || isempty(i1) || i1 < i0
    burnTime = 0;
else
    burnTime = t(i1) - t(i0);
end

% Impulse
Itot = trapz(t, F);

% Prop mass estimate
mprop = out.mprop_est;

% Delivered Isp
Isp_del = Itot / (max(mprop, eps) * g0);

% Avg pressure over burn window
if burnTime > 0
    t_b  = t(i0:i1);
    Pc_b = Pc(i0:i1);
    Pc_avg = trapz(t_b, Pc_b) / max(t_b(end) - t_b(1), eps);
else
    Pc_avg = mean(Pc);
end

% Peak chamber pressure
Pc_peak = Pc_pk;

% Kn = Ab/At (trim Ab to N too)
Ab = out.Ab(:);
Ab = Ab(1:min(numel(Ab), N));
if numel(Ab) < N
    Ab(end+1:N,1) = 0;
end

Kn_t   = Ab / At;
Kn_init = Kn_t(max(1, min(i0, N)));
Kn_peak = max(Kn_t);

% Port/throat ratio (use initial port area from geometry)
st0 = cfg.geo.geoFcn(0);
Ap0 = st0.Ap;
portThroat = Ap0 / At;

% Peak throat mass flux: Gt = mdot/At
Gt_pk = max(mdot_noz / At);

% ---- Peak port mass flux (burn-window max, OpenMotor-comparable) ----
KG_M2_TO_LB_IN2 = 2.20462262185 / (39.3700787402^2); % 0.00142233...

% Burn-window indices clamped to [1,N]
if ~isempty(i0) && ~isempty(i1) && i1 >= i0
    i0c = max(1, min(i0, N));
    i1c = max(1, min(i1, N));
    idx = (i0c:i1c).';
else
    idx = (1:N).';
end

mask = isfinite(mdot_noz(idx)) & mdot_noz(idx) > 0 & ...
    isfinite(Ap_eff(idx))   & Ap_eff(idx) > 0 & ...
    isfinite(Pc(idx))       & Pc(idx) > Pa;

valid = idx(mask);

% Effective (series) port mass flux peak
if ~isempty(valid)
    G_eff = mdot_noz(valid) ./ Ap_eff(valid); % kg/m^2/s
    s.peakPortMassFlux_kg_m2_s  = max(G_eff);
    s.peakPortMassFlux_lb_in2_s = s.peakPortMassFlux_kg_m2_s * KG_M2_TO_LB_IN2;
else
    s.peakPortMassFlux_kg_m2_s  = NaN;
    s.peakPortMassFlux_lb_in2_s = NaN;
end

% Per-segment peak fluxes over the same valid points
if ~isempty(Ap_eff_seg) && ~isempty(valid)
    nSeg = size(Ap_eff_seg, 2);
    s.peakPortMassFluxSeg_kg_m2_s  = nan(1, nSeg);
    s.peakPortMassFluxSeg_lb_in2_s = nan(1, nSeg);

    for j = 1:nSeg
        Apj = Ap_eff_seg(valid, j);
        vj = isfinite(Apj) & Apj > 0;

        if any(vj)
            Gj = mdot_noz(valid(vj)) ./ Apj(vj); % kg/m^2/s
            s.peakPortMassFluxSeg_kg_m2_s(j)  = max(Gj);
            s.peakPortMassFluxSeg_lb_in2_s(j) = s.peakPortMassFluxSeg_kg_m2_s(j) * KG_M2_TO_LB_IN2;
        end
    end
else
    s.peakPortMassFluxSeg_kg_m2_s  = [];
    s.peakPortMassFluxSeg_lb_in2_s = [];
end

% ---- Ideal and delivered thrust coefficients ----
Cf_ideal_pk = nozzleCfIdeal(cfg.noz.eps, cfg.prop.gamma, Pc_peak, cfg.noz.Pa);

den = trapz(t, Pc * At);
Cf_del = Itot / max(den, eps);

% ---- Volume loading (from geo.const) ----
Vcase  = cfg.geo.const.Vcase;
Vprop0 = cfg.geo.const.Vprop0;
volLoading = Vprop0 / Vcase;

% ---- Motor designation ----
s.designation = suggestDesignation(Itot, burnTime, 0.07);

% Bundle
s.Itot = Itot;
s.Isp_delivered = Isp_del;
s.burnTime = burnTime;

s.volLoading = volLoading;

s.Pc_avg = Pc_avg;
s.Pc_peak = Pc_peak;

s.Kn_init = Kn_init;
s.Kn_peak = Kn_peak;

s.Cf_ideal = Cf_ideal_pk;
s.Cf_delivered = Cf_del;

s.mprop = mprop;
s.propLength = cfg.geo.const.Lprop;

s.portThroat = portThroat;
s.peakMassFlux = Gt_pk;
end

function des = suggestDesignation(Itot, burnTime, tol)
% NAR/TRA-style motor designation: Letter = total impulse class, Number = avg thrust, (tol%)

if nargin < 3 || isempty(tol)
    tol = 0.07;
end

% Class upper bounds in N*s (TRA/NAR HPR convention)
classLetters = { ...
    "1/4A","1/2A","A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P" };
classMaxNs = [ ...
    0.625, 1.25, 2.50, 5, 10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20480, 40960, 81920 ];

% OpenMotor-style: choose class based on a conservative (lower) impulse with tolerance
I_nom = Itot * (1 - tol);

idx = find(I_nom <= classMaxNs, 1, 'first');
if isempty(idx)
    letter = "Above P";
else
    letter = classLetters{idx};
end

F_avg = Itot / max(burnTime, eps);

des = sprintf('%s%d (%.0f%%)', letter, round(F_avg), tol*100);
end
