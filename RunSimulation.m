clc;clear;

cfg = makeMotorConfig();

% ===== Motor container =====
cfg.grain.case_inner_diameter_m = 0.0508;   % motor ID / case ID (m)
cfg.grain.free_volume_m3 = 3e-4;            % head-end free volume (m^3)

% ===== Grain stack =====
cfg.stack = { ...
    struct("type","bates", "L",0.1524, "Ro",0.0254, "Ri",0.0127, "inhibitEnds",true), ...
    struct("type","bates", "L",0.1524, "Ro",0.0254, "Ri",0.0127, "inhibitEnds",true) ...
};

cfg.noz.thrust_efficiency = 0.959; % default parameter


% ===== Nozzle =====
cfg.noz.throat_diameter_m = 0.01524;
cfg.noz.exit_diameter_m   = 0.0762;
cfg.noz.efficiency = 1;
cfg.noz.divergence_half_angle_deg = 15;
cfg.noz.convergence_half_angle_deg = 45;
cfg.noz.throat_length_m = 0.0;
cfg.noz.slag_buildup_coeff   = 0;
cfg.noz.throat_erosion_coeff = 0;
cfg.noz.Cd = 0.98;
cfg.noz.Pa = 101325;


% ===== Propellant =====
cfg.prop.name = "Blue Thunder";
cfg.prop.density_kgm3 = 1625.09;
cfg.prop.Pmin_Pa = 1e6;
cfg.prop.Pmax_Pa = 20e6; 
cfg.prop.a = 6.9947e-5;          % keep consistent with your burn-law units
cfg.prop.n = 0.321;
cfg.prop.gamma = 1.235;
cfg.prop.Tc_K = 2616.5;
cfg.prop.M_kgmol = 0.022959;

% ===== Numerics =====
cfg.num.dt = 5e-4;
cfg.num.tmax = 8;

[out, summary, derived] = runMotor(cfg);

disp(derived.noz)