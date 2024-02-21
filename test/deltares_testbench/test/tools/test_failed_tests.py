import pytest
from tools.failed_tests import get_failed_tests


class TestFailedTests:
    def test_get_failed_tests(self, tmpdir):
        # Arrange
        expected_result = ['test1', 'test3']
        test_log = tmpdir.join("test_log.log")
        with open(str(test_log), 'a') as f:
            f.write("|test1|INFO|ERROR\n")
            f.write("|test2|INFO|SUCCESS\n")
            f.write("|test3|INFO|ERROR\n")

        # Act
        result = get_failed_tests(str(test_log))

        # Assert
        assert result == expected_result
