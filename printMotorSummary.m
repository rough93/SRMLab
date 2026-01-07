function printMotorSummary(s)
fprintf('\n=== Motor Summary (%s) ===\n', s.designation);
fprintf('Total Impulse (N·s):          %.1f\n', s.Itot);
fprintf('Delivered Isp (s):            %.2f\n', s.Isp_delivered);
fprintf('Burn time (s):                %.3f\n', s.burnTime);
fprintf('Volume loading (Vprop/Vcase): %.3f\n', s.volLoading);
fprintf('Avg chamber pressure (psi):   %.3f\n', ((s.Pc_avg/1e6)*145.038));
fprintf('Peak chamber pressure (psi):  %.3f\n', ((s.Pc_peak/1e6)*145.038));
fprintf('Initial Kn:                   %.3f\n', s.Kn_init);
fprintf('Peak Kn:                      %.3f\n', s.Kn_peak);
fprintf('Ideal Cf (at peak Pc):        %.3f\n', s.Cf_ideal);
fprintf('Delivered Cf:                 %.3f\n', s.Cf_delivered);
fprintf('Prop mass (kg):               %.4f\n', s.mprop);
fprintf('Prop length (m):              %.3f\n', s.propLength);
fprintf('Port/Throat area ratio:       %.3f\n', s.portThroat);
fprintf('Peak port mass flux (kg/m2s): %.1f\n', s.peakPortMassFlux_kg_m2_s);
fprintf('Peak port mass flux (lb/in2*s): %.2f\n', s.peakPortMassFlux_lb_in2_s);

if isfield(s,'peakPortMassFlux_GN_tag') && s.peakPortMassFlux_GN_tag > 0
    fprintf('Peak mass flux (G:%d): %.2f lb/(in^2*s)\n', ...
        s.peakPortMassFlux_GN_tag, s.peakPortMassFlux_GN_lb_in2_s);
end

if isfield(s,'peakPortMassFluxSeg_lb_in2_s') && ~isempty(s.peakPortMassFluxSeg_lb_in2_s)
    for i = 1:numel(s.peakPortMassFluxSeg_lb_in2_s)
        fprintf('Peak mass flux seg %d (lb/in2*s): %.2f\n', ...
            i, s.peakPortMassFluxSeg_lb_in2_s(i));
    end
end
fprintf('===========================\n\n');
end
