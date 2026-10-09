function tests = TestValidatorFailing
%TESTVALIDATORFAILING Controlled failing fixture for validator reporting.

tests = functiontests(localfunctions);
end

function testControlledFailure(testCase)
testCase.verifyEqual(1, 2, 'Intentional validator fixture failure.');
end
