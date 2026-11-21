function s = fresp_summary(f, G, G0)
%FRESP_SUMMARY Corner frequencies and low-frequency attenuation of one normalized response.
%
%   s = fresp_summary(f, G, G0)
%
%   f   1 x nf Hz, ascending;  G  1 x nf complex, normalized so G -> 1 at
%   high frequency;  G0  the value at f = 0 (real)
%
% Fields
%   G0        the gain at zero frequency, as a fraction of the frozen gain
%   dB0       the same in dB (20 log10)
%   f_half    where the attenuation in dB is half gone: |G| = sqrt(G0)
%   f_90      above this frequency |G| stays within 10% of the frozen
%             gain (|G| >= 0.9 from here up): "adaptation has let go"
%   f_3dB     above this frequency |G| stays within 3 dB of the frozen
%             gain (|G| >= 0.708)
%   phase_max the largest phase lead in degrees, and f_phase where it is
%   peak      the largest |G| (above 1 is a resonance of the neuron itself)
%
% A crossing is interpolated in log frequency. NaN when it is outside the
% grid (f_90 and f_3dB are the grid's first point when |G| never goes
% below the level).
%
% See also: fresp_run_one

a = abs(G);
s = struct('G0', G0, 'dB0', 20 * log10(abs(G0)));
s.f_half = last_crossing(f, a, sqrt(abs(G0)));
s.f_90 = last_crossing(f, a, 0.9);
s.f_3dB = last_crossing(f, a, 10^(-3 / 20));
ph = angle(G) * 180 / pi;
[s.phase_max, i] = max(ph);
s.f_phase = f(i);
s.peak = max(a);
end

function fc = last_crossing(f, a, level)
below = find(a < level, 1, 'last');
if isempty(below)
    fc = f(1);
elseif below == numel(f)
    fc = NaN;
else
    t = (level - a(below)) / (a(below + 1) - a(below));
    fc = exp(log(f(below)) + t * (log(f(below + 1)) - log(f(below))));
end
end
