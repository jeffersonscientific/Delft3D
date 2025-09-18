import io
import logging

from pytest import LogCaptureFixture

from ci_tools.teamcity.logging import TeamCityFormatter, enter_test_context


class TestTeamCityServiceMessageFormatter:
    def test_format(self, caplog: LogCaptureFixture) -> None:
        formatter = TeamCityFormatter()

        logger = logging.getLogger()
        logger.setLevel(logging.DEBUG)

        stream = io.StringIO()
        handler = logging.StreamHandler(stream)
        handler.setFormatter(formatter)
        logger.addHandler(handler)

        logger.info("foo")
        with enter_test_context("test_1", logger):
            logger.warning("bar")
        logger.error("baz")
        try:
            with enter_test_context("test_2", logger):
                logger.debug("qux")
                raise ValueError("Kaboom!")
        except ValueError:
            pass
        logger.info("quux")

        assert logger
