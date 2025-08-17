#!/usr/bin/env python3
"""Pin and tag the appropriate builds in TeamCity."""

from ci_tools.dimrset_delivery.dimr_context import (
    DimrAutomationContext,
    create_context_from_args,
    parse_common_arguments,
)
from ci_tools.dimrset_delivery.services import Services


class PinAndTagHelper:
    """Helper class to pin and tag builds in TeamCity and Git.

    This class uses the provided automation context to interact with TeamCity and Git clients,
    pinning and tagging builds as required for DIMRset delivery.

    Attributes
    ----------
    context : DimrAutomationContext
        The automation context containing necessary clients and configuration.

    Methods
    -------
    pin_and_tag_builds() -> None
        Pin and tag the appropriate builds.
    """

    def __init__(self, context: DimrAutomationContext, services: Services) -> None:
        self.__context = context
        self.__dry_run = context.dry_run
        self.__dry_run_prefix = context.settings.dry_run_prefix
        self.__kernel_versions = context.kernel_versions
        self.__dimr_version = context.dimr_version
        self.__build_id = context.build_id
        self.__teamcity_ids = context.settings.teamcity_ids

        self.__git_client = services.git
        self.__teamcity = services.teamcity

    def pin_and_tag_builds(self) -> None:
        """Pin and tag the appropriate builds.

        Parameters
        ----------
        context : DimrAutomationContext
            The automation context containing necessary clients and configuration.
        """
        self.__context.log("Pinning and tagging builds...")

        if self.__teamcity is None:
            raise ValueError("TeamCity client is required but not initialized")
        if self.__git_client is None:
            raise ValueError("Git client is required but not initialized")

        if self.__dry_run:
            self.__context.log("Pin and tag TC builds")
            self.__context.log(
                f"{self.__dry_run_prefix} Would tag git commit with:",
                f"commit={self.__kernel_versions['build.vcs.number']}, tag=DIMRset_{self.__dimr_version}",
            )
        else:
            self.__pin_and_tag_builds_teamcity()
            self.__git_client.tag_commit(self.__kernel_versions["build.vcs.number"], f"DIMRset_{self.__dimr_version}")

        self.__context.log("Build pinning and tagging completed successfully!")

    def __pin_and_tag_builds_teamcity(self) -> None:
        """Tag all builds and pin the appropriate builds in TeamCity."""
        if self.__teamcity is None:
            raise ValueError("TeamCity client is required but not initialized")

        tag = f"DIMRset_{self.__dimr_version}"
        self.__teamcity.add_tag_to_build_with_dependencies(self.__build_id, tag=tag)
        # Only pin specific builds
        teamcity_ids_list = list(vars(self.__teamcity_ids).values())
        build_ids_to_pin = self.__teamcity.get_dependent_build_ids_with_filter(self.__build_id, teamcity_ids_list)
        build_ids_to_pin.append(self.__build_id)
        for build_id in build_ids_to_pin:
            self.__teamcity.pin_build(build_id=build_id)


if __name__ == "__main__":
    args = parse_common_arguments()
    context = create_context_from_args(args, require_atlassian=False, require_ssh=False)
    services = Services(context)

    context.log("Starting build pinning and tagging...")
    helper = PinAndTagHelper(context=context, services=services)
    helper.pin_and_tag_builds()
    context.log("Finished")
