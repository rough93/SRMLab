function geo = geo_stack(cfg)
%GEO_STACK Combine multiple grain segments into one equivalent geometry.
% Single regression depth x is applied to every segment.

segs = cfg.stack;

% ---- constants expected by computeMotorSummary ----
% Total propellant (stack) length
Ltot = 0;
Vprop0 = 0;

for i = 1:numel(segs)
    s = segs{i};
    Ltot = Ltot + s.L;

    % Initial prop volume for this segment (annulus)
    Vprop0 = Vprop0 + pi*(s.Ro^2 - s.Ri^2)*s.L;
end

% Case volume over the propellant length (needs case inner diameter)
Rcase = cfg.grain.case_inner_diameter_m/2;
Vcase = pi*Rcase^2 * Ltot;

geo.const.Vcase  = Vcase;
geo.const.Vprop0 = Vprop0;
geo.const.Lprop  = Ltot;


% Optional: automatically inhibit internal interfaces for BATES-like segments
% If you want OpenMotor-like behavior where touching faces do NOT burn:
segs = autoInhibitInternalFaces(segs);

V0 = 0;
if isfield(cfg, "grain") && isfield(cfg.grain, "free_volume_m3")
    V0 = cfg.grain.free_volume_m3;
end

geo.geoFcn = @(x) evalStack(max(x,0));

    function st = evalStack(x)
        Ab_sum   = 0;
        Vc_sum   = V0;
        dVdx_sum = 0;

        all_done = true;
        Ap_seg = nan(1, numel(segs));
        for i = 1:numel(segs)
            s = segs{i};

            switch lower(s.type)
                case "bates"
                    sti = geo_bates_segment(x, s);

                case "finocyl"
                    % Implement later
                    sti = geo_finocyl_segment(x, s);

                otherwise
                    error("geo_stack:UnknownType", "Unknown segment type: %s", s.type);
            end
            Ap_seg(i) = sti.Ap;
            Ab_sum   = Ab_sum   + sti.Ab;
            Vc_sum   = Vc_sum   + sti.Vc;
            dVdx_sum = dVdx_sum + sti.dVdx;

            all_done = all_done && sti.done;
        end

        st.Ab   = Ab_sum;
        st.Vc   = Vc_sum;
        st.dVdx = dVdx_sum;
        st.done = all_done;

        Ap_valid = Ap_seg(isfinite(Ap_seg) & Ap_seg > 0);
        if ~isempty(Ap_valid)
            st.Ap = min(Ap_valid);     % series restriction
        else
            st.Ap = NaN;
        end

        st.Ap_seg = Ap_seg;            % per-grain port areas (for per-grain G)
    end
end

function segs = autoInhibitInternalFaces(segs)
% Auto inhibit ends at internal interfaces for segments that burn ends.
% This is a simple default: touching faces do not burn.

if numel(segs) <= 1
    return;
end

for i = 1:numel(segs)
    if ~isfield(segs{i}, "inhibitEnds")
        segs{i}.inhibitEnds = true;
    end
end

% Only allow the outermost ends to burn if user explicitly wants them
% Keep whatever user set for first and last, but inhibit internal faces
% With this simple model, "inhibitEnds" is per segment and affects both ends,
% so we cannot selectively inhibit only one face without more fields.
% If you want that later, add inhibitHead / inhibitAft booleans.
end
