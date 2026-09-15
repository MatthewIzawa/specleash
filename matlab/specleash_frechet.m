function [d, path, links] = specleash_frechet(P, Q)
%SPECLEASH_FRECHET Exact discrete Frechet distance between two polylines.
%   d = SPECLEASH_FRECHET(P, Q) is the discrete Frechet distance between the
%   curves P (n-by-2) and Q (m-by-2).
%   [d, path, links] = SPECLEASH_FRECHET(P, Q) also returns the matching as a
%   k-by-2 array of index pairs [i j] from [1 1] to [n m], and the length of
%   each link along it. The longest link equals d.
%
%   Recurrence of Eiter and Mannila (1994), computed iteratively: curve length
%   is limited by memory, not by MATLAB's recursion limit.
%   Part of specleash (MIT licence).

validateattributes(P, {'numeric'}, {'2d', 'ncols', 2, 'nonempty', 'finite'}, mfilename, 'P', 1);
validateattributes(Q, {'numeric'}, {'2d', 'ncols', 2, 'nonempty', 'finite'}, mfilename, 'Q', 2);
n = size(P, 1);
m = size(Q, 1);
D = hypot(P(:, 1) - Q(:, 1).', P(:, 2) - Q(:, 2).');
C = zeros(n, m);
C(1, 1) = D(1, 1);
for j = 2:m
    C(1, j) = max(C(1, j - 1), D(1, j));
end
for i = 2:n
    C(i, 1) = max(C(i - 1, 1), D(i, 1));
    for j = 2:m
        best = C(i - 1, j - 1);
        if C(i - 1, j) < best, best = C(i - 1, j); end
        if C(i, j - 1) < best, best = C(i, j - 1); end
        C(i, j) = max(best, D(i, j));
    end
end
d = C(n, m);

if nargout > 1
    path = zeros(n + m - 1, 2);
    k = 1;
    i = n;
    j = m;
    path(k, :) = [i j];
    while i > 1 || j > 1
        if i == 1
            j = j - 1;
        elseif j == 1
            i = i - 1;
        else
            [~, which] = min([C(i - 1, j - 1), C(i - 1, j), C(i, j - 1)]);  % ties prefer the diagonal
            if which == 1
                i = i - 1; j = j - 1;
            elseif which == 2
                i = i - 1;
            else
                j = j - 1;
            end
        end
        k = k + 1;
        path(k, :) = [i j];
    end
    path = flipud(path(1:k, :));
    links = D(sub2ind([n m], path(:, 1), path(:, 2)));
end
end
