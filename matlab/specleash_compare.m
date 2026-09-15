function [d, info] = specleash_compare(P, Q, wlRange, options)
%SPECLEASH_COMPARE Frechet comparison of two spectra, after Gartrelle et al. (2021).
%   [d, info] = SPECLEASH_COMPARE(P, Q, [lo hi]) compares spectrum P (by
%   Gartrelle's convention the asteroid) with spectrum Q (the meteorite) over
%   lo-hi micrometres. P and Q are n-by-2 arrays [wavelength_um reflectance]
%   or names of two-column text files. Both must cover the whole range.
%
%   Name-value options:
%     NormaliseAt  wavelength (um) at which each spectrum is divided to unity,
%                  applied before cropping (default: none)
%     XScale       multiplier on wavelength, or "unit" to map the range onto [0 1] (default 1)
%     Plot         draw the leash plot (default false)
%     Labels       legend labels, {labelP, labelQ}
%
%   info.path, info.links, info.curveP, info.curveQ, info.bottleneck (wavelengths, um).
%   Same conventions and numbers as specleash.compare in Python.

arguments
    P
    Q
    wlRange (1, 2) double
    options.NormaliseAt double = []
    options.XScale = 1
    options.Plot (1, 1) logical = false
    options.Labels (1, 2) cell = {'P', 'Q'}
end

lo = wlRange(1);
hi = wlRange(2);
if ~(lo < hi)
    error('specleash:range', 'Range must have lo < hi.');
end
S1 = prepare(load_spectrum(P, 'P'), lo, hi, options.NormaliseAt, options.Labels{1});
S2 = prepare(load_spectrum(Q, 'Q'), lo, hi, options.NormaliseAt, options.Labels{2});
C1 = [scale_x(S1(:, 1), lo, hi, options.XScale), S1(:, 2)];
C2 = [scale_x(S2(:, 1), lo, hi, options.XScale), S2(:, 2)];
[d, path, links] = specleash_frechet(C1, C2);
[~, b] = max(links);
info = struct('path', path, 'links', links, 'curveP', C1, 'curveQ', C2, ...
    'bottleneck', [S1(path(b, 1), 1), S2(path(b, 2), 1)]);

if options.Plot
    figure;
    hold on;
    for k = 1:size(path, 1)
        plot([C1(path(k, 1), 1) C2(path(k, 2), 1)], [C1(path(k, 1), 2) C2(path(k, 2), 2)], ...
            'Color', [0.78 0.78 0.78], 'LineWidth', 0.4, 'HandleVisibility', 'off');
    end
    plot(C1(:, 1), C1(:, 2), '-', 'Color', [0.12 0.37 0.66], 'LineWidth', 1.5, 'DisplayName', options.Labels{1});
    plot(C2(:, 1), C2(:, 2), '-', 'Color', [0.69 0.23 0.18], 'LineWidth', 1.5, 'DisplayName', options.Labels{2});
    plot([C1(path(b, 1), 1) C2(path(b, 2), 1)], [C1(path(b, 1), 2) C2(path(b, 2), 2)], ...
        'Color', [0.88 0.56 0.04], 'LineWidth', 2.5, 'DisplayName', sprintf('d = %.4f', d));
    axis equal;
    xlabel('Wavelength (\mum)');
    ylabel('Normalised reflectance');
    legend('Location', 'best');
    title(sprintf('Discrete Frechet distance: %.6f', d));
    hold off;
end
end

function S = load_spectrum(X, name)
if ischar(X) || isstring(X)
    X = readmatrix(X, 'FileType', 'text');
end
validateattributes(X, {'numeric'}, {'2d', 'nonempty'}, 'specleash_compare', name);
S = sortrows(X(:, 1:2), 1);
S = S(all(isfinite(S), 2), :);
if any(diff(S(:, 1)) <= 0)
    error('specleash:wavelength', '%s: wavelength must be strictly increasing.', name);
end
end

function S = prepare(S, lo, hi, normaliseAt, name)
tol = 1e-9;
if ~isempty(normaliseAt)
    if normaliseAt < S(1, 1) - tol || normaliseAt > S(end, 1) + tol
        error('specleash:coverage', '%s covers %.4f-%.4f um; cannot normalise at %g um.', ...
            name, S(1, 1), S(end, 1), normaliseAt);
    end
    v = interp1(S(:, 1), S(:, 2), normaliseAt, 'linear', 'extrap');  % within tol of the data, as numpy.interp
    if ~(v > 0)
        error('specleash:normalise', '%s reflectance at %g um is not positive.', name, normaliseAt);
    end
    S(:, 2) = S(:, 2) / v;
end
if S(1, 1) > lo + tol || S(end, 1) < hi - tol
    error('specleash:coverage', '%s covers %.4f-%.4f um, not %g-%g um.', name, S(1, 1), S(end, 1), lo, hi);
end
inside = S(:, 1) > lo + tol & S(:, 1) < hi - tol;
% 'extrap' only ever reaches within tol of the data (coverage is checked above);
% without it interp1 returns NaN there, where numpy.interp returns the edge value.
S = [lo, interp1(S(:, 1), S(:, 2), lo, 'linear', 'extrap'); S(inside, :); ...
     hi, interp1(S(:, 1), S(:, 2), hi, 'linear', 'extrap')];
end

function x = scale_x(w, lo, hi, xs)
if (ischar(xs) || isstring(xs)) && strcmp(xs, 'unit')
    x = (w - lo) / (hi - lo);
elseif isnumeric(xs) && isscalar(xs) && xs > 0
    x = xs * w;
else
    error('specleash:xscale', 'XScale must be a positive number or "unit".');
end
end
