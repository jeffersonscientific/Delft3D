"""Tests for pin_and_tag_builds.py."""

from unittest.mock import Mock, call, patch

import pytest

from ci_tools.dimrset_delivery.dimr_context import DimrAutomationContext
from ci_tools.dimrset_delivery.lib.git_client import GitClient
from ci_tools.dimrset_delivery.lib.teamcity import TeamCity
from ci_tools.dimrset_delivery.pin_and_tag_builds import PinAndTagHelper
from ci_tools.dimrset_delivery.services import Services
from ci_tools.dimrset_delivery.settings.teamcity_settings import Settings, TeamcityIds


class TestPinAndTagBuilds:
    """Test cases for pin_and_tag_builds function."""

    def setup_method(self) -> None:
        """Set up test fixtures."""
        self.mock_context = Mock(spec=DimrAutomationContext)
        self.mock_context.build_id = "12345"
        self.mock_context.dry_run = False
        self.mock_context.settings = Mock(spec=Settings)
        self.mock_context.settings.dry_run_prefix = "[TEST]"
        self.mock_context.log = Mock()
        self.mock_context.kernel_versions = {"build.vcs.number": "abc123def"}
        self.mock_context.dimr_version = "1.2.3"

        self.mock_services = Mock(spec=Services)
        self.mock_services.git_client = Mock(spec=GitClient)
        self.mock_services.teamcity = Mock(spec=TeamCity)
        self.mock_services.teamcity.teamcity_ids = Mock(spec=TeamcityIds)
        self.mock_services.teamcity.teamcity_ids.delft3d_windows_collect_build_type_id = "windows_build_type"
        self.mock_services.teamcity.teamcity_ids.delft3d_linux_collect_build_type_id = "linux_build_type"

        self.mock_context.settings.teamcity_ids = self.mock_services.teamcity.teamcity_ids

    def test_pin_and_tag_builds_success(self) -> None:
        """Test successful pinning and tagging of builds."""
        # Arrange
        self.mock_services.teamcity.get_dependent_build_ids_with_filter.return_value = ["id1", "id2"]
        helper = PinAndTagHelper(self.mock_context, self.mock_services)

        # Act
        helper.pin_and_tag_builds()

        # Assert
        self.mock_context.log.assert_has_calls(
            [call("Pinning and tagging builds..."), call("Build pinning and tagging completed successfully!")]
        )
        self.mock_services.teamcity.add_tag_to_build_with_dependencies.assert_called_once_with(
            "12345", tag="DIMRset_1.2.3"
        )
        teamcity_ids_list = list(vars(self.mock_context.settings.teamcity_ids).values())
        self.mock_services.teamcity.get_dependent_build_ids_with_filter.assert_called_once_with(
            "12345", teamcity_ids_list
        )
        self.mock_services.teamcity.pin_build.assert_any_call(build_id="id1")
        self.mock_services.teamcity.pin_build.assert_any_call(build_id="id2")
        self.mock_services.teamcity.pin_build.assert_any_call(build_id="12345")  # Also pins the original build
        self.mock_services.git.tag_commit.assert_called_once_with("abc123def", "DIMRset_1.2.3")

    def test_pin_and_tag_builds_dry_run(self) -> None:
        """Test dry run mode - should only print what would be done."""
        # Arrange
        self.mock_context.dry_run = True
        helper = PinAndTagHelper(self.mock_context, self.mock_services)

        # Act
        with patch("builtins.print"):
            helper.pin_and_tag_builds()

        # Assert
        self.mock_context.log.assert_any_call("Pinning and tagging builds...")
        self.mock_context.log.assert_any_call("Pin and tag TC builds")
        # Check for the combined log call with both arguments
        expected_log_call = ("[TEST] Would tag git commit with:", "commit=abc123def, tag=DIMRset_1.2.3")
        assert expected_log_call in [call.args for call in self.mock_context.log.call_args_list]

    def test_pin_and_tag_builds_missing_teamcity_client(self) -> None:
        """Test error when TeamCity client is None."""
        # Arrange
        self.mock_services.teamcity = None
        helper = PinAndTagHelper(self.mock_context, self.mock_services)

        # Act & Assert
        with pytest.raises(ValueError, match="TeamCity client is required but not initialized"):
            helper.pin_and_tag_builds()

    def test_pin_and_tag_builds_missing_git_client(self) -> None:
        """Test error when Git client is None."""
        # Arrange
        self.mock_services.git = None
        helper = PinAndTagHelper(self.mock_context, self.mock_services)

        # Act & Assert
        with pytest.raises(ValueError, match="Git client is required but not initialized"):
            helper.pin_and_tag_builds()

    def test_pin_and_tag_builds_with_different_versions(self) -> None:
        """Test with different kernel and DIMR versions."""
        # Arrange
        self.mock_context.kernel_versions = {"build.vcs.number": "xyz789abc"}
        self.mock_context.dimr_version = "2.0.0"
        self.mock_services.teamcity.get_dependent_build_ids_with_filter.return_value = ["id3"]
        helper = PinAndTagHelper(self.mock_context, self.mock_services)

        # Act
        helper.pin_and_tag_builds()

        # Assert
        self.mock_services.teamcity.add_tag_to_build_with_dependencies.assert_called_once_with(
            "12345", tag="DIMRset_2.0.0"
        )
        teamcity_ids_list = list(vars(self.mock_context.settings.teamcity_ids).values())
        self.mock_services.teamcity.get_dependent_build_ids_with_filter.assert_called_once_with(
            "12345", teamcity_ids_list
        )
        self.mock_services.teamcity.pin_build.assert_any_call(build_id="id3")
        self.mock_services.teamcity.pin_build.assert_any_call(build_id="12345")  # Also pins the original build
        self.mock_services.git.tag_commit.assert_called_once_with("xyz789abc", "DIMRset_2.0.0")

    def test_pin_and_tag_builds_success_message(self) -> None:
        """Test that success message is printed after completion."""
        # Arrange
        self.mock_services.teamcity.get_dependent_build_ids_with_filter.return_value = []
        helper = PinAndTagHelper(self.mock_context, self.mock_services)

        # Act
        helper.pin_and_tag_builds()

        # Assert
        self.mock_context.log.assert_called_with("Build pinning and tagging completed successfully!")

    def test_pin_and_tag_builds_teamcity_exception(self) -> None:
        """Test when TeamCity client raises an exception."""
        # Arrange
        self.mock_services.teamcity.add_tag_to_build_with_dependencies.side_effect = Exception("TeamCity error")
        helper = PinAndTagHelper(self.mock_context, self.mock_services)

        # Act & Assert
        with pytest.raises(Exception, match="TeamCity error"):
            helper.pin_and_tag_builds()

    def test_pin_and_tag_builds_git_client_exception(self) -> None:
        """Test when Git client raises an exception."""
        # Arrange
        self.mock_services.teamcity.get_dependent_build_ids_with_filter.return_value = []
        self.mock_services.git.tag_commit.side_effect = Exception("Git error")
        helper = PinAndTagHelper(self.mock_context, self.mock_services)

        # Act & Assert
        with pytest.raises(Exception, match="Git error"):
            helper.pin_and_tag_builds()

    def test_pin_and_tag_builds_dry_run_teamcity_client_required(self) -> None:
        """Test that dry run without teamcity client."""
        # Arrange
        self.mock_context.dry_run = True
        self.mock_services.teamcity = None
        helper = PinAndTagHelper(self.mock_context, self.mock_services)

        # Act & Assert - should raise an exception
        with pytest.raises(ValueError, match="TeamCity client is required but not initialized"):
            helper.pin_and_tag_builds()

    def test_pin_and_tag_builds_dry_run_git_client_required(self) -> None:
        """Test that dry run without git client."""
        # Arrange
        self.mock_context.dry_run = True
        self.mock_services.git = None
        helper = PinAndTagHelper(self.mock_context, self.mock_services)

        # Act & Assert - should raise an exception
        with pytest.raises(ValueError, match="Git client is required but not initialized"):
            helper.pin_and_tag_builds()

    def test_pin_and_tag_builds_context_method_calls(self) -> None:
        """Test that all expected context methods are called in the correct order."""
        # Arrange
        self.mock_services.teamcity.get_dependent_build_ids_with_filter.return_value = []
        helper = PinAndTagHelper(self.mock_context, self.mock_services)

        # Act
        helper.pin_and_tag_builds()

        # Assert method call order
        assert self.mock_context.log.call_count == 2
        self.mock_context.log.assert_has_calls(
            [call("Pinning and tagging builds..."), call("Build pinning and tagging completed successfully!")]
        )
