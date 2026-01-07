function F = nozzleThrustFromCf(Pc, noz, gamma)
Pa = noz.Pa;

Pc_eff = max(Pc, Pa);
Cf = nozzleCfIdeal(noz.eps, gamma, Pc_eff, Pa);

eta = 1;
if isfield(noz,"efficiency") && ~isempty(noz.efficiency)
    eta = eta * noz.efficiency;
end
if isfield(noz,"thrust_efficiency") && ~isempty(noz.thrust_efficiency)
    eta = eta * noz.thrust_efficiency;
end

F = eta .* Cf .* Pc .* noz.At;

F(Pc <= Pa) = 0;
end
