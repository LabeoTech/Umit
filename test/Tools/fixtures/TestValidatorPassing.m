function tests = TestValidatorPassing
%TESTVALIDATORPASSING Passing function-based test fixture for the validator.

tests = functiontests(localfunctions);
end

function testSampleFunction(testCase)
testCase.verifyEqual(validatorSampleFunction(2), 3);
end
