function qc = getTformQCMetrics(tform)
%GETTFORMQCMETRICS Extract convention-safe QC metrics from a 2-D transform.
%
%   qc = GETTFORMQCMETRICS(tform) returns a scalar struct describing the
%   geometric content of a 2-D image transform:
%
%       qc.translationXY - [tx ty] translation in pixels
%       qc.rotationDeg   - rotation angle in degrees, counter-clockwise
%       qc.scaleXY       - [sx sy] column norms of the linear part
%       qc.determinant   - determinant of the linear part
%       qc.matrixA       - the transform as a premultiply 3-by-3 matrix
%       qc.convention    - 'premultiply', 'postmultiply' or 'numeric', i.e.
%                          the convention detected on the input
%       qc.transformType - geometric family of the transform
%                          ('translation', 'rigid', 'similarity', 'affine',
%                          'projective' or 'unknown')
%
%   MATLAB exposes two incompatible conventions for 2-D transforms:
%
%       * the legacy AFFINE2D / PROJECTIVE2D classes store a POSTMULTIPLY
%         matrix T, where [x' y' 1] = [x y 1] * T, so the translation lives
%         in T(3,1:2) and the linear part is T(1:2,1:2);
%       * the modern TRANSLTFORM2D / RIGIDTFORM2D / SIMTFORM2D /
%         AFFINETFORM2D / PROJTFORM2D classes store a PREMULTIPLY matrix A,
%         where [x'; y'; 1] = A * [x; y; 1], so the translation lives in
%         A(1:2,3) and the linear part is A(1:2,1:2).
%
%   The two linear parts are transposes of one another, so reading the same
%   element index from both yields a different translation and an opposite
%   rotation sign. This function normalises every supported input to the
%   premultiply convention first, so callers obtain the same metrics
%   whichever class IMREGTFORM or a saved resource happens to return.
%
%   A numeric 3-by-3 input is treated as a premultiply matrix unless its
%   last COLUMN is [0 0 1]' and its last ROW is not, which unambiguously
%   identifies a postmultiply matrix.
%
%   Non-finite or unsupported inputs produce NaN metrics rather than an
%   error, so QC reporting never aborts a registration run.
%
%   See also IMREGTFORM, AFFINE2D, AFFINETFORM2D.

qc = struct( ...
    'translationXY', [NaN NaN], ...
    'rotationDeg', NaN, ...
    'scaleXY', [NaN NaN], ...
    'determinant', NaN, ...
    'matrixA', nan(3), ...
    'convention', 'unknown', ...
    'transformType', iTransformType(tform));

[A, convention] = iResolveMatrixA(tform);
if isempty(A)
    return
end

qc.matrixA = A;
qc.convention = convention;

linearPart = A(1:2,1:2);
qc.translationXY = A(1:2,3)';
qc.scaleXY = [norm(linearPart(:,1)), norm(linearPart(:,2))];
qc.determinant = det(linearPart);
qc.rotationDeg = atan2d(linearPart(2,1), linearPart(1,1));

end

% =========================================================================
% Local helpers
% =========================================================================
function transformType = iTransformType(tform)
%ITRANSFORMTYPE Map a transform object to its geometric family.

classMap = { ...
    'transltform2d',  'translation'; ...
    'rigidtform2d',   'rigid'; ...
    'simtform2d',     'similarity'; ...
    'affinetform2d',  'affine'; ...
    'projtform2d',    'projective'; ...
    'affine2d',       'affine'; ...
    'projective2d',   'projective'};

transformType = 'unknown';
for k = 1:size(classMap,1)
    if isa(tform, classMap{k,1})
        transformType = classMap{k,2};
        return
    end
end

end

function [A, convention] = iResolveMatrixA(tform)
%IRESOLVEMATRIXA Return the transform as a premultiply 3-by-3 matrix.

A = [];
convention = 'unknown';

if isnumeric(tform)
    if ~isequal(size(tform), [3 3]) || any(~isfinite(tform), 'all')
        return
    end
    candidate = double(tform);
    if isequal(candidate(1:3,3)', [0 0 1]) && ~isequal(candidate(3,1:3), [0 0 1])
        % Last column is [0 0 1]': unambiguously a postmultiply matrix.
        A = candidate';
        convention = 'postmultiply';
    else
        A = candidate;
        convention = 'numeric';
    end
    return
end

% Modern premultiply classes expose A. They also still expose a deprecated
% postmultiply T for backwards compatibility, so A must be tested first;
% AFFINE2D and PROJECTIVE2D expose only T.
if isprop(tform, 'A')
    candidate = double(tform.A);
    convention = 'premultiply';
    if isequal(size(candidate), [3 3]) && all(isfinite(candidate), 'all')
        A = candidate;
    end
    return
end

if isprop(tform, 'T')
    candidate = double(tform.T);
    convention = 'postmultiply';
    if isequal(size(candidate), [3 3]) && all(isfinite(candidate), 'all')
        A = candidate';
    end
end

end
