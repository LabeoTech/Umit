function varargout = rng(varargin)
%RNG Test-only shadow of the built-in `rng` function.
%
%   When the environment variable UMIT_TEST_MOCK_RNG_MODE is 'fixed', always
%   reseeds the global random stream with a fixed seed, regardless of the
%   requested mode (e.g. the 'shuffle' argument NormalisationFiltering
%   passes in). This makes fminsearch's randomized initial guess
%   reproducible across repeated calls within a single test, which is
%   otherwise not guaranteed (tracked separately as DFR-20260819-001).
%
%   In any other state (variable unset or different), this folder is on the
%   path outside a test (e.g. addpath(genpath(pwd))) and the call is passed
%   through to the real RNG, with its arguments and outputs unchanged.
%
%   A test enables the fixed seed with:
%       setenv('UMIT_TEST_MOCK_RNG_MODE', 'fixed')
%   and should restore it afterwards (e.g. via onCleanup or
%   matlab.unittest.fixtures.EnvironmentVariableFixture).

if strcmp(getenv('UMIT_TEST_MOCK_RNG_MODE'), 'fixed')
    s = RandStream('mt19937ar', 'Seed', 12345);
    RandStream.setGlobalStream(s);
    return
end

% Pass through to the real function. It is an M-file, so run it from its own
% folder, where it takes precedence over this shadow.
paths = which('rng', '-all');
paths = paths(~startsWith(paths, fileparts(mfilename('fullpath'))));
here = pwd;
cd(fileparts(paths{1}));
cleanup = onCleanup(@() cd(here)); %#ok<NASGU>
[varargout{1:nargout}] = rng(varargin{:});

end
