#!/bin/bash

# Run tests
R -e "x <- devtools::test(quiet=T)" > /test_result.txt

# Read test results from file
res=$(cat /test_result.txt)

# # return test result as an output (will be deprecated)
echo ::set-output name=test_result::$res

# return exit status
failed_tests_count=$(grep -i -c "Failed tests" /test_result.txt)

if (( failed_tests_count > 0 )); then
    echo "::error:: There are failed tests that failed, please download artifacts for more details."
else
    echo "All tests passed."
    exit 0
fi