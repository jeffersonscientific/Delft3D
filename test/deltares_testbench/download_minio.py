import argparse

from src.config.credentials import Credentials
from src.utils.handlers.minio_handler import MinIOHandler
from src.utils.logging.i_logger import ILogger


class PrintLogger(ILogger):
    def debug(self, message: str) -> None:
        print(f"DEBUG: {message}")

    def error(self, message: str) -> None:
        print(f"ERROR: {message}")

    def exception(self, message: str) -> None:
        print(f"EXCEPTION: {message}")

    def info(self, message: str) -> None:
        print(f"INFO: {message}")

    def log(self, message: str) -> None:
        print(f"LOG: {message}")

    def warning(self, message: str) -> None:
        print(f"WARNING: {message}")


def create_argument_parser() -> argparse.ArgumentParser:
    """Create custom argument parser."""
    parser = argparse.ArgumentParser(description="Retrieve status of a testbench running on TeamCity")

    parser.add_argument("-u", "--username", help="Username for accessing TeamCity.", dest="username")
    parser.add_argument(
        "-p",
        "--password",
        help="Password belonging to username for accessing TeamCity.",
        dest="password",
    )
    parser.add_argument(
        "-e",
        "--engine_dir",
        help="Engine directory to checkout.",
        dest="engine_dir",
    )
    return parser


parser = create_argument_parser()
args = parser.parse_args()
credentials = Credentials()
credentials.username = args.username
credentials.password = args.password
logger = PrintLogger()

handler = MinIOHandler()
handler.download(
    from_path=f"https://s3.deltares.nl/dsc-testbench/cases/{args.engine_dir}",
    to_path=".",
    credentials=credentials,
    version=None,
    logger=logger,
)
