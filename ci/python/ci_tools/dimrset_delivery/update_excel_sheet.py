#!/usr/bin/env python3
"""Update the Excel sheet with this week's release information."""

from datetime import datetime, timezone
from typing import List

from openpyxl import load_workbook
from openpyxl.worksheet.worksheet import Worksheet

from ci_tools.dimrset_delivery.common_utils import ResultTestBankParser, get_testbank_result_parser
from ci_tools.dimrset_delivery.dimr_context import (
    DimrAutomationContext,
    create_context_from_args,
    parse_common_arguments,
)
from ci_tools.dimrset_delivery.lib.ssh_client import Direction
from ci_tools.dimrset_delivery.services import Services


class ExcelHelper:
    """Object responsible for updating the Excel sheet."""

    def __init__(
        self,
        context: DimrAutomationContext,
        parser: ResultTestBankParser,
    ) -> None:
        """
        Create a new instance of ExcelHelper.

        Args:
            teamcity (TeamCity): A wrapper for the TeamCity REST API.
            filepath (str): Path to the Excel file.
            dimr_version (str): The DIMR version to update the Excel for.
            kernel_versions (Dict[str, str]): A dictionary mapping kernel names to their version.
            parser (ResultTestBankParser): A parser for the latest test bench results.
        """
        self.__filepath = context.settings.versions_excel_filename
        self.__dimr_version = context.dimr_version
        self.__kernel_versions = context.kernel_versions
        self.__parser = parser
        self.__sheet_name = context.settings.sheet_name
        self.__name_column = context.settings.name_column

    def append_row(self) -> None:
        """Append a new row to the Excel sheet with this week's DIMR information."""
        row = self.__prepare_row_to_insert()
        print(row)

        try:
            workbook = load_workbook(filename=self.__filepath)
            worksheet = workbook[self.__sheet_name]

            if self.__worksheet_already_contains_row(worksheet):
                print("The Excel sheet already contains a row for this DIMRset value and revision number.")
                return

            worksheet.append(row)
            workbook.save(filename=self.__filepath)
        except Exception as e:
            print("Could not update the excel: \n")
            print(e)
        finally:
            workbook.close()

    def __prepare_row_to_insert(self) -> List[str]:
        """Prepare a row to be inserted in the Excel sheet."""
        row = []

        row.append("")  # Column A (empty column)
        row.append(str(datetime.now(tz=timezone.utc).date()))  # Column B (Date)
        row.append(f"DIMRset {self.__dimr_version}")  # Column C (DIMR version)
        row.append("")  # Column D (Revision)
        row.append("FLOW1D2D now in GitHub")  # Column E (Flow1D)
        row.append("OSS")  # Column F (FlowFM)
        row.append(self.__kernel_versions["build.vcs.number"])  # Column G (OSS)
        row.append("DRR now in GitHub")  # Column H (RR)
        row.append("FBC now in GitHub")  # Column I (FBC)
        row.append(self.__parser.get_percentage_total_passing())  # Column J (Percentage passing)
        row.append(self.__parser.get_total_tests())  # Column K (Total Number of          cases)
        row.append(self.__parser.get_total_passing())  # Column L (      Number of green    cases)
        row.append(self.__parser.get_total_failing())  # Column M (      Number of red      cases)
        row.append(self.__parser.get_total_exceptions())  # Column N (      Number of crashing cases)
        row.append("")  # Column O (Docker hub)
        row.append("Flow1D and RR: only Windows")  # Column P (Remarks)

        return row

    def __worksheet_already_contains_row(self, worksheet: Worksheet) -> bool:
        """
        Check if the Excel sheet already contains a row for the given DIMRset.

        Return True if the Excel sheet already contains such a row.
        """
        name_already_exists = False

        name_column = worksheet[self.__name_column]

        for cell in name_column:
            if cell.value == f"DIMRset {self.__dimr_version}":
                name_already_exists = True

        return name_already_exists


def update_excel_sheet(context: DimrAutomationContext, services: Services) -> None:
    """Update the Excel sheet with this week's release information.

    Parameters
    ----------
    context : DimrAutomationContext
        The automation context containing necessary clients and configuration.
    """
    context.log("Updating Excel sheet...")

    if context.dry_run:
        context.log(f"Would update Excel sheet with DIMR version: {context.dimr_version}")
        context.log("Would download Excel from network drive")
        context.log("Would append new row with release information")
        context.log("Would upload updated Excel back to network drive")
        return

    parser = get_testbank_result_parser(context.settings.path_to_release_test_results_artifact)
    path_to_excel_file = f"/p/d-hydro/dimrset/{context.settings.versions_excel_filename}"

    if services.ssh is None:
        raise ValueError("SSH client is required but not initialized")
    if services.teamcity is None:
        raise ValueError("TeamCity client is required but not initialized")

    services.ssh.secure_copy(
        context.settings.versions_excel_filename,
        path_to_excel_file,
        Direction.FROM,
    )
    helper = ExcelHelper(
        context=context,
        parser=parser,
    )
    helper.append_row()
    services.ssh.secure_copy(
        context.settings.versions_excel_filename,
        path_to_excel_file,
        Direction.TO,
    )

    context.log("Excel sheet update completed successfully!")


if __name__ == "__main__":
    args = parse_common_arguments()
    context = create_context_from_args(args, require_atlassian=False, require_git=False)
    services = Services(context)

    context.log("Starting Excel sheet update...")
    update_excel_sheet(context, services)
    context.log("Finished")
