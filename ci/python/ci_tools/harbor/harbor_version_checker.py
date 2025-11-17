import argparse
import re
from dataclasses import dataclass

import requests

# Configuration
BASE_URL = "https://containers.deltares.nl/api/v2.0/projects/delft3d/repositories/delft3dfm/artifacts"
PAGE_SIZE = 100  # Maximum page size for Harbor API


@dataclass(order=True)
class YearVersion:
    """Represents a year.version tag (e.g., 2026.01, release-2025.02)."""

    year: int
    version: int
    tag: str


@dataclass(order=True)
class SemVer:
    """Represents a semantic version tag (e.g., 2.30.06-development)."""

    major: int
    minor: int
    patch: int
    tag: str


def parse_year_version_tags(tags: list[str]) -> list[YearVersion]:
    """
    Parse tags matching year.version format (e.g., 2026.01, release-2025.02).

    Args:
        tags: List of tag strings

    Returns
    -------
        List of YearVersion objects sorted by version descending
    """
    year_version_pattern = re.compile(r"^(\d{4})\.(\d{2})(?:-release)?$")
    year_versions = []

    for tag in tags:
        # Handle both "release-YYYY.MM" and "YYYY.MM-release" formats
        clean_tag = tag.replace("release-", "").replace("-release", "")
        match = year_version_pattern.match(clean_tag)
        if match:
            year = int(match.group(1))
            version = int(match.group(2))
            year_versions.append(YearVersion(year=year, version=version, tag=tag))

    year_versions.sort(reverse=True)
    return year_versions


def parse_semver_tags(tags: list[str]) -> list[SemVer]:
    """
    Parse tags matching semantic versioning format (e.g., 2.30.06-development).

    Args:
        tags: List of tag strings

    Returns
    -------
        List of SemVer objects sorted by version descending
    """
    semver_pattern = re.compile(r"^(\d+)\.(\d+)\.(\d+)(?:-[\w]+)?$")
    semvers = []

    for tag in tags:
        match = semver_pattern.match(tag)
        if match:
            major = int(match.group(1))
            minor = int(match.group(2))
            patch = int(match.group(3))
            semvers.append(SemVer(major=major, minor=minor, patch=patch, tag=tag))

    semvers.sort(reverse=True)
    return semvers


def fetch_all_tags(username: str, password: str) -> list[str]:
    """
    Fetch all tags from Harbor API artifacts.

    Args:
        username: Harbor API username
        password: Harbor API password

    Returns
    -------
        List of all tag names from all artifacts
    """
    # Fetch all artifacts with pagination
    all_artifacts = []
    page = 1

    while True:
        url = f"{BASE_URL}?page={page}&page_size={PAGE_SIZE}"
        print(f"Fetching page {page}...")
        response = requests.get(url, auth=(username, password))

        if response.status_code != 200:
            print(f"Error: Failed to query API. Status code: {response.status_code}")
            print(f"Response: {response.text}")
            break

        artifacts = response.json()
        if not artifacts:  # No more results
            break

        all_artifacts.extend(artifacts)
        page += 1

    print(f"\nTotal artifacts fetched: {len(all_artifacts)}")

    # Extract all tags from all artifacts
    all_tags = []
    for artifact in all_artifacts:
        if artifact.get("tags"):
            for tag in artifact["tags"]:
                tag_name = tag.get("name")
                if tag_name:
                    all_tags.append(tag_name)

    return all_tags


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.

    Returns
    -------
        Parsed command line arguments
    """
    parser = argparse.ArgumentParser(description="Check Harbor repository for latest version tags")
    parser.add_argument(
        "--harbor-username",
        type=str,
        required=True,
        help="Harbor API username",
    )
    parser.add_argument(
        "--harbor-password",
        type=str,
        required=True,
        help="Harbor API password",
    )
    return parser.parse_args()


def main() -> None:
    """Fetch artifacts from Harbor API and display latest version tags."""
    args = parse_arguments()
    all_tags = fetch_all_tags(args.harbor_username, args.harbor_password)

    # Print the list of tags
    print("\nAll tags found:")
    for tag in sorted(all_tags):
        print(f"  - {tag}")

    print(f"\nTotal tags: {len(all_tags)}")

    # Find latest versions
    year_versions = parse_year_version_tags(all_tags)
    semvers = parse_semver_tags(all_tags)

    # Create latest_versions dictionary
    latest_versions = {
        "latest_year_version": year_versions[0].tag if year_versions else None,
        "latest_semver": semvers[0].tag if semvers else None,
        "all_year_versions": [yv.tag for yv in year_versions],
        "all_semvers": [sv.tag for sv in semvers],
    }

    print("\n" + "=" * 60)
    print("LATEST VERSIONS")
    print("=" * 60)
    print(f"Latest Year.Version tag: {latest_versions['latest_year_version']}")
    print(f"Latest Semver tag:       {latest_versions['latest_semver']}")

    print("\nAll Year.Version tags found:")
    if latest_versions["all_year_versions"]:
        for tag in latest_versions["all_year_versions"]:
            print(f"  - {tag}")
    else:
        print("  None found")

    print("\nAll Semver tags found:")
    if latest_versions["all_semvers"]:
        for tag in latest_versions["all_semvers"]:
            print(f"  - {tag}")
    else:
        print("  None found")


if __name__ == "__main__":
    main()
