"""Tests for prepare_email.py."""

from unittest.mock import Mock, patch

import pytest

from ci_tools.dimrset_delivery.dimr_context import DimrAutomationContext
from ci_tools.dimrset_delivery.helpers.result_testbank_parser import ResultTestBankParser
from ci_tools.dimrset_delivery.prepare_email import (
    get_previous_testbank_result_parser,
    get_tag_from_build_info,
    parse_version,
    prepare_email,
)


class TestPrepareEmail:
    """Test cases for prepare_email function."""

    @patch("ci_tools.dimrset_delivery.prepare_email.EmailHelper")
    @patch("ci_tools.dimrset_delivery.prepare_email.get_testbank_result_parser")
    @patch("ci_tools.dimrset_delivery.prepare_email.get_previous_testbank_result_parser")
    def test_prepare_email_success(
        self,
        mock_get_previous_parser: Mock,
        mock_get_parser: Mock,
        mock_email_helper: Mock,
    ) -> None:
        """Test successful email preparation."""
        # Arrange
        mock_context = Mock(spec=DimrAutomationContext)
        mock_context.dry_run = False
        mock_context.get_kernel_versions.return_value = {"kernel1": "1.0.0", "kernel2": "2.0.0"}
        mock_context.get_dimr_version.return_value = "1.2.3"

        mock_current_parser = Mock(spec=ResultTestBankParser)
        mock_previous_parser = Mock(spec=ResultTestBankParser)
        mock_get_parser.return_value = mock_current_parser
        mock_get_previous_parser.return_value = mock_previous_parser

        mock_helper_instance = Mock()
        mock_email_helper.return_value = mock_helper_instance

        # Act
        prepare_email(mock_context)

        # Assert
        mock_context.print_status.assert_called_once_with("Preparing email template...")
        mock_context.get_kernel_versions.assert_called_once()
        mock_context.get_dimr_version.assert_called_once()
        mock_get_parser.assert_called_once()
        mock_get_previous_parser.assert_called_once_with(mock_context)

        mock_email_helper.assert_called_once_with(
            dimr_version="1.2.3",
            kernel_versions={"kernel1": "1.0.0", "kernel2": "2.0.0"},
            current_parser=mock_current_parser,
            previous_parser=mock_previous_parser,
        )
        mock_helper_instance.generate_template.assert_called_once()

    @patch("ci_tools.dimrset_delivery.prepare_email.get_testbank_result_parser")
    @patch("builtins.print")
    def test_prepare_email_dry_run(
        self,
        mock_print: Mock,
        mock_get_parser: Mock,
    ) -> None:
        """Test prepare_email in dry run mode."""
        # Arrange
        mock_context = Mock(spec=DimrAutomationContext)
        mock_context.dry_run = True
        mock_context.get_kernel_versions.return_value = {"kernel1": "1.0.0"}
        mock_context.get_dimr_version.return_value = "1.2.3"

        # Act
        prepare_email(mock_context)

        # Assert
        mock_context.print_status.assert_called_once_with("Preparing email template...")
        mock_context.get_kernel_versions.assert_called_once()
        mock_context.get_dimr_version.assert_called_once()
        mock_get_parser.assert_not_called()

        # Check that print was called with the correct arguments
        mock_print.assert_called_once_with("[DRY-RUN] Would prepare email template for DIMR version:", "1.2.3")


class TestGetPreviousTestbankResultParser:
    """Test cases for get_previous_testbank_result_parser function."""

    def test_no_teamcity_client_raises_error(self) -> None:
        """Test that missing TeamCity client raises ValueError."""
        # Arrange
        mock_context = Mock(spec=DimrAutomationContext)
        mock_context.teamcity = None

        # Act & Assert
        with pytest.raises(ValueError, match="TeamCity client is required but not initialized"):
            get_previous_testbank_result_parser(mock_context)

    def test_no_current_build_info_returns_none(self) -> None:
        """Test that missing current build info returns None."""
        # Arrange
        mock_context = Mock(spec=DimrAutomationContext)
        mock_context.build_id = "12345"
        mock_context.teamcity = Mock()
        mock_context.teamcity.get_full_build_info_for_build_id.return_value = None

        # Act
        result = get_previous_testbank_result_parser(mock_context)

        # Assert
        assert result is None
        mock_context.teamcity.get_full_build_info_for_build_id.assert_called_once_with("12345")

    def test_no_build_type_id_returns_none(self) -> None:
        """Test that missing buildTypeId returns None."""
        # Arrange
        mock_context = Mock(spec=DimrAutomationContext)
        mock_context.build_id = "12345"
        mock_context.teamcity = Mock()
        mock_context.teamcity.get_full_build_info_for_build_id.return_value = {}

        # Act
        result = get_previous_testbank_result_parser(mock_context)

        # Assert
        assert result is None

    def test_no_builds_returns_none(self) -> None:
        """Test that no builds available returns None."""
        # Arrange
        mock_context = Mock(spec=DimrAutomationContext)
        mock_context.build_id = "12345"
        mock_context.teamcity = Mock()
        mock_context.teamcity.get_full_build_info_for_build_id.return_value = {
            "buildTypeId": "bt123",
            "tags": {"tag": [{"name": "DIMRset_1.2.3"}]},
        }
        mock_context.teamcity.get_builds_for_build_type_id.return_value = None

        # Act
        result = get_previous_testbank_result_parser(mock_context)

        # Assert
        assert result is None

    @patch("ci_tools.dimrset_delivery.prepare_email.ResultTestBankParser")
    def test_successful_previous_parser_retrieval(self, mock_parser_class: Mock) -> None:
        """Test successful retrieval of previous testbank result parser."""
        # Arrange
        mock_context = Mock(spec=DimrAutomationContext)
        mock_context.build_id = "12345"
        mock_context.teamcity = Mock()

        current_build_info = {"buildTypeId": "bt123", "tags": {"tag": [{"name": "DIMRset_1.2.3"}]}}

        builds_response = {
            "build": [
                {"id": 12346},  # Different build
                {"id": 12344},  # Previous build
            ]
        }

        previous_build_info = {"buildTypeId": "bt123", "tags": {"tag": [{"name": "DIMRset_1.2.2"}]}}

        mock_artifact_content = b"artifact content"

        mock_context.teamcity.get_full_build_info_for_build_id.side_effect = [
            current_build_info,  # For current build
            None,  # For first loop build (12346)
            previous_build_info,  # For second loop build (12344)
        ]
        mock_context.teamcity.get_builds_for_build_type_id.return_value = builds_response
        mock_context.teamcity.get_build_artifact.return_value = mock_artifact_content

        mock_parser_instance = Mock(spec=ResultTestBankParser)
        mock_parser_class.return_value = mock_parser_instance

        # Act
        result = get_previous_testbank_result_parser(mock_context)

        # Assert
        assert result == mock_parser_instance
        mock_parser_class.assert_called_once_with("artifact content")
        mock_context.teamcity.get_build_artifact.assert_called_once()

    def test_no_previous_version_found_returns_none(self) -> None:
        """Test that no previous version found returns None."""
        # Arrange
        mock_context = Mock(spec=DimrAutomationContext)
        mock_context.build_id = "12345"
        mock_context.teamcity = Mock()

        current_build_info = {"buildTypeId": "bt123", "tags": {"tag": [{"name": "DIMRset_1.2.3"}]}}

        builds_response = {
            "build": [
                {"id": 12346},  # Different build
            ]
        }

        other_build_info = {
            "buildTypeId": "bt123",
            "tags": {"tag": [{"name": "DIMRset_1.2.4"}]},  # Higher version
        }

        mock_context.teamcity.get_full_build_info_for_build_id.side_effect = [
            current_build_info,  # For current build
            other_build_info,  # For other build
        ]
        mock_context.teamcity.get_builds_for_build_type_id.return_value = builds_response

        # Act
        result = get_previous_testbank_result_parser(mock_context)

        # Assert
        assert result is None

    def test_no_artifact_returns_none(self) -> None:
        """Test that missing artifact returns None."""
        # Arrange
        mock_context = Mock(spec=DimrAutomationContext)
        mock_context.build_id = "12345"
        mock_context.teamcity = Mock()

        current_build_info = {"buildTypeId": "bt123", "tags": {"tag": [{"name": "DIMRset_1.2.3"}]}}

        builds_response = {
            "build": [
                {"id": 12344},  # Previous build
            ]
        }

        previous_build_info = {"buildTypeId": "bt123", "tags": {"tag": [{"name": "DIMRset_1.2.2"}]}}

        mock_context.teamcity.get_full_build_info_for_build_id.side_effect = [
            current_build_info,  # For current build
            previous_build_info,  # For previous build
        ]
        mock_context.teamcity.get_builds_for_build_type_id.return_value = builds_response
        mock_context.teamcity.get_build_artifact.return_value = None

        # Act
        result = get_previous_testbank_result_parser(mock_context)

        # Assert
        assert result is None


class TestGetTagFromBuildInfo:
    """Test cases for get_tag_from_build_info function."""

    def test_no_tags_returns_default(self) -> None:
        """Test that build info without tags returns default tuple."""
        # Arrange
        build_info = {}

        # Act
        result = get_tag_from_build_info(build_info)

        # Assert
        assert result == (0, 0, 0)

    def test_empty_tags_returns_default(self) -> None:
        """Test that empty tags returns default tuple."""
        # Arrange
        build_info = {"tags": {"tag": []}}

        # Act
        result = get_tag_from_build_info(build_info)

        # Assert
        assert result == (0, 0, 0)

    def test_no_dimrset_tag_returns_default(self) -> None:
        """Test that tags without DIMRset prefix return default tuple."""
        # Arrange
        build_info = {"tags": {"tag": [{"name": "some_other_tag"}, {"name": "another_tag"}]}}

        # Act
        result = get_tag_from_build_info(build_info)

        # Assert
        assert result == (0, 0, 0)

    def test_valid_dimrset_tag_returns_version(self) -> None:
        """Test that valid DIMRset tag returns parsed version."""
        # Arrange
        build_info = {"tags": {"tag": [{"name": "some_other_tag"}, {"name": "DIMRset_1.2.3"}, {"name": "another_tag"}]}}

        # Act
        result = get_tag_from_build_info(build_info)

        # Assert
        assert result == (1, 2, 3)

    def test_multiple_dimrset_tags_returns_last_valid(self) -> None:
        """Test that multiple DIMRset tags returns the last valid one."""
        # Arrange
        build_info = {
            "tags": {"tag": [{"name": "DIMRset_1.0.0"}, {"name": "DIMRset_2.3.4"}, {"name": "some_other_tag"}]}
        }

        # Act
        result = get_tag_from_build_info(build_info)

        # Assert
        assert result == (2, 3, 4)

    def test_valid_non_standard_dimrset_tag_returns_version(self) -> None:
        """Test that valid DIMRset tag with non-standard format returns parsed version."""
        # Arrange
        build_info = {
            "tags": {
                "tag": [
                    {"name": "some_other_tag"},
                    {"name": "DIMRset_1.2"},  # Valid but only major.minor
                    {"name": "another_tag"},
                ]
            }
        }

        # Act
        result = get_tag_from_build_info(build_info)

        # Assert
        assert result == (1, 2)

    def test_invalid_dimrset_tag_returns_default(self) -> None:
        """Test that invalid DIMRset tag format returns default tuple."""
        # Arrange
        build_info = {
            "tags": {
                "tag": [
                    {"name": "DIMRset_invalid"},
                    {"name": "DIMRset_a.b.c"},  # Completely invalid
                ]
            }
        }

        # Act
        result = get_tag_from_build_info(build_info)

        # Assert
        assert result == (0, 0, 0)


class TestParseVersion:
    """Test cases for parse_version function."""

    def test_valid_dimrset_tag_returns_tuple(self) -> None:
        """Test that valid DIMRset tag returns version tuple."""
        # Act & Assert
        assert parse_version("DIMRset_1.2.3") == (1, 2, 3)
        assert parse_version("DIMRset_10.20.30") == (10, 20, 30)
        assert parse_version("DIMRset_0.1.0") == (0, 1, 0)
        # These are also valid even if not standard semantic versions
        assert parse_version("DIMRset_1.2") == (1, 2)
        assert parse_version("DIMRset_1.2.3.4") == (1, 2, 3, 4)

    def test_invalid_prefix_returns_none(self) -> None:
        """Test that tag without DIMRset prefix returns None."""
        # Act & Assert
        assert parse_version("NotDIMRset_1.2.3") is None
        assert parse_version("1.2.3") is None
        assert parse_version("") is None

    def test_invalid_version_format_returns_none(self) -> None:
        """Test that invalid version format returns None."""
        # Act & Assert
        # Note: "1.2" is actually valid and returns (1, 2)
        # Note: "1.2.3.4" is actually valid and returns (1, 2, 3, 4)
        # Only truly invalid formats return None
        assert parse_version("DIMRset_a.b.c") is None
        assert parse_version("DIMRset_") is None
        assert parse_version("DIMRset_1.") is None
        assert parse_version("DIMRset_.1") is None

    def test_none_input_returns_none(self) -> None:
        """Test that None input returns None."""
        # Act & Assert
        assert parse_version(None) is None

    def test_empty_string_returns_none(self) -> None:
        """Test that empty string returns None."""
        # Act & Assert
        assert parse_version("") is None

    def test_dimrset_prefix_only_returns_none(self) -> None:
        """Test that DIMRset prefix without version returns None."""
        # Act & Assert
        assert parse_version("DIMRset_") is None


class TestIntegration:
    """Integration tests for prepare_email functionality."""

    @patch("ci_tools.dimrset_delivery.prepare_email.EmailHelper")
    @patch("ci_tools.dimrset_delivery.prepare_email.get_testbank_result_parser")
    @patch("ci_tools.dimrset_delivery.prepare_email.get_previous_testbank_result_parser")
    def test_prepare_email_with_no_previous_parser(
        self,
        mock_get_previous_parser: Mock,
        mock_get_parser: Mock,
        mock_email_helper: Mock,
    ) -> None:
        """Test email preparation when no previous parser is available."""
        # Arrange
        mock_context = Mock(spec=DimrAutomationContext)
        mock_context.dry_run = False
        mock_context.get_kernel_versions.return_value = {"kernel1": "1.0.0"}
        mock_context.get_dimr_version.return_value = "1.2.3"

        mock_current_parser = Mock(spec=ResultTestBankParser)
        mock_get_parser.return_value = mock_current_parser
        mock_get_previous_parser.return_value = None  # No previous parser available

        mock_helper_instance = Mock()
        mock_email_helper.return_value = mock_helper_instance

        # Act
        prepare_email(mock_context)

        # Assert
        mock_email_helper.assert_called_once_with(
            dimr_version="1.2.3",
            kernel_versions={"kernel1": "1.0.0"},
            current_parser=mock_current_parser,
            previous_parser=None,
        )
        mock_helper_instance.generate_template.assert_called_once()
