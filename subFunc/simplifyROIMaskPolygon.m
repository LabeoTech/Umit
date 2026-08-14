function [pgon, outerVerticesXY, info] = simplifyROIMaskPolygon(mask, varargin)
%SIMPLIFYROIMASKPOLYGON Convert a logical ROI mask to a compact polyshape.
%
%   [PGON, OUTERVERTICESXY, INFO] = SIMPLIFYROIMASKPOLYGON(MASK) converts
%   a two-dimensional logical ROI mask into a polyshape. Each boundary is
%   simplified with MATLAB's Ramer-Douglas-Peucker implementation,
%   REDUCEPOLY, using its normalized tolerance of 0.01. Candidate boundaries
%   are rasterized and must differ from their unsimplified representation by
%   no more than 0.5 percent of the ROI pixels (minimum one pixel).
%
%   MASK may contain disconnected components and internal holes. Both are
%   retained in PGON. OUTERVERTICESXY is the largest outer boundary and is
%   provided for existing callers whose serialized vertex field supports
%   one finite [x y] loop only. INFO reports the applied tolerance and
%   boundary vertex counts.
%
%   Name-Value options:
%       Tolerance - Normalized Ramer-Douglas-Peucker tolerance in [0, 1].
%                   Default: 0.01.

p = inputParser;
p.FunctionName = 'simplifyROIMaskPolygon';
addParameter(p, 'Tolerance', 0.01, @(x) isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x >= 0 && x <= 1);
parse(p, varargin{:});

if ~(islogical(mask) || isnumeric(mask)) || ndims(mask) ~= 2 %#ok<ISMAT>
    error('simplifyROIMaskPolygon:InvalidMask', ...
        'mask must be a two-dimensional logical or numeric array.');
end

mask = logical(mask);
if ~any(mask(:))
    error('simplifyROIMaskPolygon:EmptyMask', ...
        'mask must contain at least one true pixel.');
end

tolerance = double(p.Results.Tolerance);
[boundariesRC, ~, ~, adjacency] = bwboundaries(mask, 8, 'holes');

if isempty(boundariesRC)
    error('simplifyROIMaskPolygon:NoBoundary', ...
        'mask did not produce a valid boundary.');
end

nBoundaries = numel(boundariesRC);
isHole = false(nBoundaries, 1);
nestingDepth = zeros(nBoundaries, 1);
if ~isempty(adjacency)
    % BWBOUNDARIES records each boundary's immediate container in its row.
    % Odd nesting depth is a hole; even depth is an outer boundary (including
    % an island inside a hole).
    parentIdx = zeros(nBoundaries, 1);
    for iBoundary = 1:nBoundaries
        parent = find(adjacency(iBoundary, :), 1, 'first');
        if ~isempty(parent)
            parentIdx(iBoundary) = parent;
        end
    end
    nestingDepth = iGetNestingDepth(parentIdx);
    isHole = mod(nestingDepth, 2) == 1;
end

rawBoundariesXY = cell(nBoundaries, 1);
rawVertexCounts = zeros(nBoundaries, 1);
for iBoundary = 1:nBoundaries
    boundaryRC = double(boundariesRC{iBoundary});
    boundaryXY = [boundaryRC(:, 2), boundaryRC(:, 1)];
    boundaryXY = iRemoveClosingVertex(boundaryXY);

    if size(boundaryXY, 1) < 3
        error('simplifyROIMaskPolygon:InvalidBoundary', ...
            'Mask boundary %d has fewer than three vertices.', iBoundary);
    end

    rawVertexCounts(iBoundary) = size(boundaryXY, 1);
    rawBoundariesXY{iBoundary} = boundaryXY;
end

referenceMask = iRasterizeBoundaries(rawBoundariesXY, nestingDepth, size(mask));
maxPixelDifference = max(1, round(0.005 * nnz(referenceMask)));
[reducedBoundariesXY, appliedTolerance, pixelDifference, bUsedOriginalBoundaryFallback] = ...
    iSelectPixelFaithfulBoundaries(rawBoundariesXY, nestingDepth, size(mask), ...
    tolerance, maxPixelDifference);
reducedVertexCounts = cellfun(@(xy) size(xy, 1), reducedBoundariesXY);
pgon = iBuildPolyshape(reducedBoundariesXY, nBoundaries);

outerIdx = find(~isHole);
if isempty(outerIdx)
    error('simplifyROIMaskPolygon:InvalidBoundaryHierarchy', ...
        'No outer boundary was returned for the ROI mask.');
end

outerAreas = cellfun(@(xy) abs(polyarea(xy(:, 1), xy(:, 2))), ...
    reducedBoundariesXY(outerIdx));
[~, iLargest] = max(outerAreas);
outerVerticesXY = reducedBoundariesXY{outerIdx(iLargest)};

info = struct( ...
    'tolerance', tolerance, ...
    'nBoundaries', nBoundaries, ...
    'nHoles', sum(isHole), ...
    'nInputVertices', sum(rawVertexCounts), ...
    'nOutputVertices', sum(reducedVertexCounts), ...
    'appliedTolerance', appliedTolerance, ...
    'pixelDifference', pixelDifference, ...
    'maxPixelDifference', maxPixelDifference, ...
    'usedOriginalBoundaryFallback', bUsedOriginalBoundaryFallback);

end

function [boundariesOut, appliedTolerance, pixelDifference, bUsedFallback] = ...
        iSelectPixelFaithfulBoundaries(rawBoundariesXY, nestingDepth, imageSize, tolerance, maxPixelDifference)
%ISELECTPIXELFAITHFULBOUNDARIES Select the strongest pixel-faithful reduction.

referenceMask = iRasterizeBoundaries(rawBoundariesXY, nestingDepth, imageSize);
candidateTolerance = unique([tolerance, tolerance / 2, tolerance / 4, ...
    tolerance / 10, tolerance / 20, 0], 'stable');

for iTolerance = 1:numel(candidateTolerance)
    thisTolerance = candidateTolerance(iTolerance);
    candidate = cellfun(@(xy) iReduceBoundary(xy, thisTolerance), ...
        rawBoundariesXY, 'UniformOutput', false);

    if any(cellfun(@(xy) size(xy, 1) < 3, candidate))
        continue
    end

    candidateMask = iRasterizeBoundaries(candidate, nestingDepth, imageSize);
    pixelDifference = nnz(xor(candidateMask, referenceMask));
    if pixelDifference <= maxPixelDifference
        boundariesOut = candidate;
        appliedTolerance = thisTolerance;
        bUsedFallback = thisTolerance == 0;
        return
    end
end

boundariesOut = rawBoundariesXY;
appliedTolerance = 0;
pixelDifference = 0;
bUsedFallback = true;

end

function rasterMask = iRasterizeBoundaries(boundariesXY, nestingDepth, imageSize)
%IRASTERIZEBOUNDARIES Rasterize nested boundaries with their original topology.

rasterMask = false(imageSize);
[~, order] = sort(nestingDepth, 'ascend');
for iOrder = 1:numel(order)
    iBoundary = order(iOrder);
    xy = boundariesXY{iBoundary};
    thisMask = poly2mask(xy(:, 1), xy(:, 2), imageSize(1), imageSize(2));
    if mod(nestingDepth(iBoundary), 2) == 0
        rasterMask = rasterMask | thisMask;
    else
        rasterMask = rasterMask & ~thisMask;
    end
end

end

function nestingDepth = iGetNestingDepth(parentIdx)
%IGETNESTINGDEPTH Calculate boundary nesting depth from immediate parents.

nBoundaries = numel(parentIdx);
nestingDepth = zeros(nBoundaries, 1);

for iBoundary = 1:nBoundaries
    visited = false(nBoundaries, 1);
    current = iBoundary;

    while parentIdx(current) ~= 0
        if visited(current)
            error('simplifyROIMaskPolygon:InvalidBoundaryHierarchy', ...
                'Mask boundary hierarchy contains a cycle.');
        end
        visited(current) = true;
        nestingDepth(iBoundary) = nestingDepth(iBoundary) + 1;
        current = parentIdx(current);
    end
end

end

function boundaryXY = iRemoveClosingVertex(boundaryXY)
%IREMOVECLOSINGVERTEX Convert a closed boundary to one finite vertex loop.

boundaryXY = double(boundaryXY);
boundaryXY = boundaryXY(all(isfinite(boundaryXY), 2), :);

if size(boundaryXY, 1) > 1 && isequal(boundaryXY(1, :), boundaryXY(end, :))
    boundaryXY(end, :) = [];
end

if size(boundaryXY, 1) > 1
    duplicate = [false; all(diff(boundaryXY, 1, 1) == 0, 2)];
    boundaryXY(duplicate, :) = [];
end

end

function reducedXY = iReduceBoundary(boundaryXY, tolerance)
%IREDUCEBOUNDARY Simplify a closed loop without changing its endpoints.

closedXY = [boundaryXY; boundaryXY(1, :)];
try
    reducedXY = reducepoly(closedXY, tolerance);
catch
    % REDUCEPOLY is available with Image Processing Toolbox. If a valid
    % boundary happens to trigger a toolbox edge case, preserve it exactly.
    reducedXY = closedXY;
end

reducedXY = iRemoveClosingVertex(reducedXY);

end

function pgon = iBuildPolyshape(boundariesXY, expectedBoundaryCount)
%IBUILDPOLYSHAPE Assemble NaN-separated outer and hole boundaries.

[pgon, bValid] = iCreatePolyshape(boundariesXY, expectedBoundaryCount);

if bValid
    return
end

error('simplifyROIMaskPolygon:InvalidPolyshape', ...
    'Mask boundaries could not be represented as a valid polyshape.');

end

function [pgon, bValid] = iCreatePolyshape(boundariesXY, expectedBoundaryCount)
%ICREATEPOLYSHAPE Construct and validate a polyshape from finite loops.

x = [];
y = [];
for iBoundary = 1:numel(boundariesXY)
    xy = boundariesXY{iBoundary};
    x = [x; xy(:, 1); NaN]; %#ok<AGROW>
    y = [y; xy(:, 2); NaN]; %#ok<AGROW>
end

pgon = polyshape(x, y, 'Simplify', true);
bValid = ~isempty(pgon.Vertices) && numboundaries(pgon) == expectedBoundaryCount;

end
