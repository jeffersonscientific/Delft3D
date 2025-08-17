from typing import Optional

from ci_tools.dimrset_delivery.dimr_context import DimrAutomationContext
from ci_tools.dimrset_delivery.lib.atlassian import Atlassian
from ci_tools.dimrset_delivery.lib.git_client import GitClient
from ci_tools.dimrset_delivery.lib.ssh_client import SshClient
from ci_tools.dimrset_delivery.lib.teamcity import TeamCity


class Services:
    """Dataclass for storing initialized service clients."""

    atlassian: Optional[Atlassian] = None
    teamcity: Optional[TeamCity] = None
    ssh: Optional[SshClient] = None
    git: Optional[GitClient] = None

    def __init__(self, context: DimrAutomationContext) -> None:
        """Initialize services based on requirements."""
        if context.requirements.atlassian:
            if not context.credentials.atlassian_username or not context.credentials.atlassian_password:
                raise ValueError("Atlassian credentials are required but not provided")
            self.atlassian = Atlassian(
                username=context.credentials.atlassian_username,
                password=context.credentials.atlassian_password,
                context=context,
            )
        else:
            self.atlassian = None

        if context.requirements.teamcity:
            if not context.credentials.teamcity_username or not context.credentials.teamcity_password:
                raise ValueError("TeamCity credentials are required but not provided")
            self.teamcity = TeamCity(
                username=context.credentials.teamcity_username,
                password=context.credentials.teamcity_password,
                context=context,
            )
            context.dimr_version = self.teamcity.get_dimr_version_from_context()
            context.kernel_versions = self.teamcity.get_kernel_versions_from_context()
            context.branch_name = self.teamcity.get_branch_name_from_context()
        else:
            self.teamcity = None

        if context.requirements.ssh:
            if not context.credentials.ssh_username or not context.credentials.ssh_password:
                raise ValueError("SSH credentials are required but not provided")
            self.ssh = SshClient(
                username=context.credentials.ssh_username,
                password=context.credentials.ssh_password,
                context=context,
            )
        else:
            self.ssh = None

        if context.requirements.git:
            if not context.credentials.git_username or not context.credentials.git_pat:
                raise ValueError("Git credentials are required but not provided")
            self.git = GitClient(
                username=context.credentials.git_username,
                password=context.credentials.git_pat,
                context=context,
            )
        else:
            self.git = None
