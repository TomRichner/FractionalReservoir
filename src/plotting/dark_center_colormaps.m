function [diverging, sequential, info] = dark_center_colormaps(name, n)
% DARK_CENTER_COLORMAPS A diverging colormap with a DARK center, and the sequential map from white to that center color.
%
%   [diverging, sequential, info] = DARK_CENTER_COLORMAPS(name)
%   [diverging, sequential, info] = DARK_CENTER_COLORMAPS(name, n)   % n colors each (default 256)
%   names = DARK_CENTER_COLORMAPS()                                  % the names it knows
%
% Written for maps of a signed quantity whose ZERO is the thing to see (the
% largest Lyapunov exponent lambda_1, FD-029): on a white-centered diverging
% map the zero band reads as missing data, and a value near zero cannot be
% told from an empty cell. A dark center makes lambda_1 = 0 the darkest
% band of the map, and the two signs grow lighter away from it, each in its
% own hue. The second map is for a nonnegative quantity drawn beside it
% (memory capacity): white at zero, darkening gradually to the SAME color
% as the diverging map's center, so the two rows read as one palette.
%
% DIRECTION: the first row of `diverging` is the NEGATIVE end (cool: blue
% or cyan), the last row the positive end (warm: red, orange or pink), as
% in the blue-to-red maps used before. `sequential` runs from white (row 1)
% to the center color (row n).
%
% THE MAPS
%   'guppy'         CMasher's guppy, reversed (guppy_r): teal-blue, through
%                   blue and violet, to a dark purple-magenta center
%                   (CIELAB lightness about 25), then through crimson to
%                   light orange (lightness about 75). Perceptually uniform,
%                   the two halves at the same lightness.
%   'redshift'      CMasher's redshift: light cyan through blue to a BLACK
%                   center, then through red to light pink. Its sequential
%                   companion is therefore white to black, a gray scale:
%                   draw missing data in a color that is not gray when
%                   using it.
%   'dark_magenta'  constructed here, in CIELAB: a dark magenta center
%                   (L* 22, a* 48, b* -32), straight lines in L*a*b* to a
%                   light blue end (L* 78, a* -10, b* -26) and a light red
%                   end (L* 78, a* 28, b* 18). L* rises linearly from the
%                   center to each end, at the same rate on both sides, so
%                   lightness carries the distance from zero and the hue
%                   says the sign. The ends are as light and as saturated as
%                   the sRGB gamut allows at that lightness for a straight
%                   path; every color is inside it (checked below; nothing
%                   is clipped beyond round-off).
% The sequential maps are all made the same way: a straight line in CIELAB
% from white (L* 100, a* = b* = 0) to the diverging map's center color, so
% lightness falls linearly with the value.
%
% COLOR-VISION DEFICIENCY. The signal that matters, the distance from zero,
% is carried by lightness, which a dichromat sees as anyone does. The sign
% is carried by a cool against a warm hue (blue against red or orange),
% the axis that protan and deutan vision keeps; red against green, which it
% loses, is not used.
%
% SOURCES AND LICENSE. The 'guppy' and 'redshift' tables are from CMasher
% (Ellert van der Velden, https://github.com/1313e/CMasher, files
% src/cmasher/colormaps/<name>/<name>_hex.txt, 511 colors each; documented
% at https://cmasher.readthedocs.io/user/diverging.html), which is
% distributed under the BSD 3-Clause License, Copyright (c) 2019-2021,
% Ellert van der Velden; the license permits redistribution of the tables
% with this notice. Reference: van der Velden (2020), CMasher: Scientific
% colormaps for making accessible, informative and 'cmashing' plots,
% J. Open Source Softw. 5(46), 2004. Kept here are rows 1, 9, 17, ..., 505
% and 511 of each table (65 anchors, at their own positions along the map,
% (row - 1) / 510); interpolated linearly in sRGB between neighbors at most
% 8/510 of the range apart, they reproduce the full table to within a few
% 8-bit steps. The center color is the map at exactly 0.5.
% The CIELAB conversions are the standard sRGB (IEC 61966-2-1) and CIE
% 1976 L*a*b* formulas with the D65 white point, written out below so the
% function needs no toolbox.
%
% info: name, center (the center color, RGB in [0 1]), source (text).
%
% See also: fig_mc_gain_bias_maps

NAMES = {'guppy', 'redshift', 'dark_magenta'};
if nargin == 0
    diverging = NAMES;
    return
end
if nargin < 2; n = 256; end
validateattributes(n, {'numeric'}, {'scalar', 'integer', '>=', 3});
name = validatestring(name, NAMES, mfilename, 'name');

switch name
    case 'guppy'
        % CMasher guppy as published runs warm to cool; reversed so the
        % negative end is cool.
        anchors = flipud(hex_rows(GUPPY_HEX()));
        pos = 1 - flip(CMASHER_POSITIONS());
        source = 'CMasher guppy (reversed), BSD 3-Clause, (c) Ellert van der Velden';
    case 'redshift'
        anchors = hex_rows(REDSHIFT_HEX());
        pos = CMASHER_POSITIONS();
        source = 'CMasher redshift, BSD 3-Clause, (c) Ellert van der Velden';
    case 'dark_magenta'
        center = [22 48 -32]; cool = [78 -10 -26]; warm = [78 28 18];
        t = linspace(0, 1, 129)';
        lab = [cool + (center - cool) .* t; center + (warm - center) .* t(2:end)];
        anchors = lab2srgb(lab);
        assert(all(anchors(:) >= -1e-3 & anchors(:) <= 1 + 1e-3), ...
            'dark_center_colormaps:OutOfGamut', 'the constructed map left the sRGB gamut');
        anchors = min(max(anchors, 0), 1);
        pos = linspace(0, 1, size(anchors, 1));
        source = 'constructed in CIELAB (dark magenta center, light blue and light red ends)';
end
diverging = interp1(pos(:), anchors, linspace(0, 1, n)');
center_rgb = interp1(pos(:), anchors, 0.5);

% White to the center color, a straight line in CIELAB.
t = linspace(0, 1, n)';
lab = [100 0 0] + (srgb2lab(center_rgb) - [100 0 0]) .* t;
sequential = min(max(lab2srgb(lab), 0), 1);
% The two ends exactly, not to the round-off of the round trip through L*a*b*.
sequential([1 end], :) = [1 1 1; center_rgb];

info = struct('name', name, 'center', center_rgb, 'source', source);
end

%% ------------------------------------------------------------------------
function p = CMASHER_POSITIONS()
% Where the kept rows of a 511-color CMasher table sit along the map.
p = [0:8:504, 510] / 510;
end

function rgb = hex_rows(hex)
% Rows of RGB in [0 1] from a cell of 'RRGGBB' strings.
v = hex2dec(reshape(char(hex)', 2, [])');
rgb = reshape(v, 3, [])' / 255;
end

function lab = srgb2lab(rgb)
% sRGB in [0 1] (rows) to CIE L*a*b*, D65.
lin = rgb / 12.92;
hi = rgb > 0.04045;
lin(hi) = ((rgb(hi) + 0.055) / 1.055) .^ 2.4;
xyz = lin * M_RGB2XYZ()' ./ D65();
f = xyz .^ (1/3);
lo = xyz <= (6/29)^3;
f(lo) = xyz(lo) / (3 * (6/29)^2) + 4/29;
lab = [116 * f(:, 2) - 16, 500 * (f(:, 1) - f(:, 2)), 200 * (f(:, 2) - f(:, 3))];
end

function rgb = lab2srgb(lab)
% CIE L*a*b* (rows), D65, to sRGB; not clipped, so a caller can check the gamut.
fy = (lab(:, 1) + 16) / 116;
f = [fy + lab(:, 2) / 500, fy, fy - lab(:, 3) / 200];
xyz = f .^ 3;
lo = f <= 6/29;
xyz(lo) = 3 * (6/29)^2 * (f(lo) - 4/29);
lin = (xyz .* D65()) / M_RGB2XYZ()';
rgb = 12.92 * lin;
hi = lin > 0.0031308;
rgb(hi) = 1.055 * lin(hi) .^ (1/2.4) - 0.055;
end

function M = M_RGB2XYZ()
% Linear sRGB to CIE XYZ (IEC 61966-2-1).
M = [0.4124564 0.3575761 0.1804375; 0.2126729 0.7151522 0.0721750; 0.0193339 0.1191920 0.9503041];
end

function w = D65()
w = [0.95047 1.00000 1.08883];
end

function h = GUPPY_HEX()
% CMasher guppy, every eighth of its 511 colors and the last (see the header).
h = {'FA914F' 'FB8846' 'FC7E3D' 'FD7535' 'FE6B2C' 'FE6024' 'FE541C' 'FE4716' 'FD3913' ...
     'FC2915' 'F8181C' 'F30825' 'EC022E' 'E50336' 'DD083D' 'D50E43' 'CC1349' 'C4174E' ...
     'BC1A53' 'B41C57' 'AC1E5A' 'A41F5E' '9C1F61' '942063' '8C1F65' '851F67' '7E1D69' ...
     '771C6B' '701A6C' '69176D' '63136F' '5D0E71' '580874' '5B0C7F' '5D108A' '5F1396' ...
     '6116A2' '631AAF' '641DBC' '6422C8' '6427D4' '622DDF' '5F35E8' '5A3FEE' '5449F1' ...
     '4C54F1' '445DEF' '3C67EC' '356FE7' '2E77E3' '297EDE' '2584DA' '228AD6' '2290D3' ...
     '2296D0' '249BCD' '25A1CB' '27A6C9' '29ABC7' '2AB1C6' '2AB6C4' '29BCC3' '27C2C2' ...
     '23C7C1' '1ECCBF'};
end

function h = REDSHIFT_HEX()
% CMasher redshift, every eighth of its 511 colors and the last (see the header).
h = {'94F1F3' '89E8F0' '7DDFED' '72D6EA' '67CEE8' '5BC5E6' '4FBDE4' '43B4E3' '37ACE2' ...
     '2BA3E1' '209AE0' '1991DF' '1988DE' '217FDB' '2B75D4' '356CCA' '3A64BC' '3D5CAD' ...
     '3D559E' '3C4F8F' '3A4881' '364273' '333C66' '2F3559' '2A2F4C' '252940' '1F2235' ...
     '1A1C2A' '14151F' '0D0E15' '06060A' '020202' '000000' '040202' '0D0605' '180C0A' ...
     '23120F' '2D1713' '381B17' '42201B' '4D241F' '582722' '632B26' '6F2E29' '7A302C' ...
     '863230' '923434' '9E3639' 'AB373F' 'B63848' 'C03D55' 'C74562' 'CD4F70' 'D2597D' ...
     'D7638A' 'DB6E97' 'DF78A3' 'E383AF' 'E68EBB' 'EA99C7' 'EDA4D2' 'F1AEDE' 'F5B9E9' ...
     'F9C5F4' 'FCCDFC'};
end
