function tests = TestValidatorThrows %#ok<STOUT>
%TESTVALIDATORTHROWS Fixture that fails while constructing a test suite.

error('TestValidatorThrows:setupFailure', ...
    'Intentional validator infrastructure failure.');
tests = functiontests(localfunctions); %#ok<UNRCH>
end
