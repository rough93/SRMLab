function [out, summary, derived] = runMotor(cfg)
%RUNMOTOR One-stop runner: validate -> derive -> geometry -> simulate -> summarize -> plot

% 1) Validate OpenMotor-like inputs
validateMotorConfig(cfg);

% 2) Compute derived fields (At, Ae, eps, c*)
derived = deriveMotorFields(cfg);

% 3) Attach derived nozzle fields into cfg (what the solver expects)
cfg.noz.At  = derived.noz.At;
cfg.noz.Ae  = derived.noz.Ae;
cfg.noz.eps = derived.noz.eps;

% If your nozzleMdotThrust uses Cd, Pa, eps, At, gamma, etc, this is enough.

% 4) Attach geometry object (adds cfg.geo.geoFcn and cfg.geo.const.*)
cfg.geo = geo_stack(cfg);

% 5) Run simulation
out = simulateMotor(cfg);

% ===== Propellant summary (derived c*) =====
R = 8.314462618;                    % J/mol-K
Rspec = R / cfg.prop.M_kgmol;       % J/kg-K

term = sqrt(cfg.prop.gamma/(Rspec*cfg.prop.Tc_K)) * ...
       (2/(cfg.prop.gamma+1))^((cfg.prop.gamma+1)/(2*(cfg.prop.gamma-1)));

cstar_implied = 1 / term;           % m/s

fprintf('\n--- Propellant Summary ---\n');
fprintf('c* (from γ,Tc,M) = %.1f m/s\n', cstar_implied);
fprintf('----------------------------------\n');



% 6) Summarize and display
summary = computeMotorSummary(out, cfg);
printMotorSummary(summary);
plotOpenMotorStyle(out, summary);

end
