function results = run_testIO()
% Runs the UMT I/O unit-test suite.
%
%   results = run_testIO()

results = runtests('testIO');
disp(table(results))
end