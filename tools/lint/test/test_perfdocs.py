import contextlib
import mock
import os
import pytest
import shutil
import tempfile

import mozunit

LINTER = "perfdocs"


class PerfDocsLoggerMock:
    LOGGER = None
    PATHS = []
    FAILED = True


"""
This is a sample mozperftest test that we use for testing
the verification process.
"""
SAMPLE_TEST = """
"use strict";

async function setUp(context) {
  context.log.info("setUp example!");
}

async function test(context, commands) {
  context.log.info("Test with setUp/tearDown example!");
  await commands.measure.start("https://www.sitespeed.io/");
  await commands.measure.start("https://www.mozilla.org/en-US/");
}

async function tearDown(context) {
  context.log.info("tearDown example!");
}

module.noexport = {};

module.exports = {
  setUp,
  tearDown,
  test,
  owner: "Performance Testing Team",
  name: "Example",
  description: "The description of the example test.",
  longDescription: `
  This is a longer description of the test perhaps including information
  about how it should be run locally or links to relevant information.
  `
};
"""


SAMPLE_CONFIG = """
name: mozperftest
manifest: None
static-only: False
suites:
    suite:
        description: "Performance tests from the 'suite' folder."
        tests:
            Example: ""
"""


DYNAMIC_SAMPLE_CONFIG = """
name: {}
manifest: None
static-only: False
suites:
    suite:
        description: "Performance tests from the 'suite' folder."
        tests:
            Example: ""
"""


@contextlib.contextmanager
def temp_file(name="temp", tempdir=None, content=None):
    if tempdir is None:
        tempdir = tempfile.mkdtemp()
    path = os.path.join(tempdir, name)
    if content is not None:
        with open(path, "w") as f:
            f.write(content)
    try:
        yield path
    finally:
        try:
            shutil.rmtree(tempdir)
        except FileNotFoundError:
            pass


@contextlib.contextmanager
def temp_dir():
    tempdir = tempfile.mkdtemp()
    try:
        yield tempdir
    finally:
        try:
            shutil.rmtree(tempdir)
        except FileNotFoundError:
            pass


def setup_sample_logger(logger, structured_logger, top_dir):
    from perfdocs.logger import PerfDocLogger

    PerfDocLogger.LOGGER = structured_logger
    PerfDocLogger.PATHS = ["perfdocs"]
    PerfDocLogger.TOP_DIR = top_dir

    import perfdocs.verifier as vf
    import perfdocs.gatherer as gt
    import perfdocs.generator as gn

    gt.logger = logger
    vf.logger = logger
    gn.logger = logger


@mock.patch("perfdocs.generator.Generator")
@mock.patch("perfdocs.verifier.Verifier")
@mock.patch("perfdocs.logger.PerfDocLogger", new=PerfDocsLoggerMock)
def test_perfdocs_start_and_fail(generator, verifier, structured_logger, config, paths):
    from perfdocs.perfdocs import run_perfdocs

    with temp_file("bad", content="foo") as temp:
        run_perfdocs(config, logger=structured_logger, paths=[temp], generate=False)
        assert PerfDocsLoggerMock.LOGGER == structured_logger
        assert PerfDocsLoggerMock.PATHS == [temp]
        assert PerfDocsLoggerMock.FAILED

    assert verifier.validate_tree.assert_called_once()
    assert generator.generate_perfdocs.assert_not_called()


@mock.patch("perfdocs.generator.Generator")
@mock.patch("perfdocs.verifier.Verifier")
@mock.patch("perfdocs.logger.PerfDocLogger", new=PerfDocsLoggerMock)
def test_perfdocs_start_and_pass(generator, verifier, structured_logger, config, paths):
    from perfdocs.perfdocs import run_perfdocs

    PerfDocsLoggerMock.FAILED = False
    with temp_file("bad", content="foo") as temp:
        run_perfdocs(config, logger=structured_logger, paths=[temp], generate=False)
        assert PerfDocsLoggerMock.LOGGER == structured_logger
        assert PerfDocsLoggerMock.PATHS == [temp]
        assert not PerfDocsLoggerMock.FAILED

    assert verifier.validate_tree.assert_called_once()
    assert generator.generate_perfdocs.assert_called_once()


@mock.patch("perfdocs.logger.PerfDocLogger", new=PerfDocsLoggerMock)
def test_perfdocs_bad_paths(structured_logger, config, paths):
    from perfdocs.perfdocs import run_perfdocs

    with pytest.raises(Exception):
        run_perfdocs(config, logger=structured_logger, paths=["bad"], generate=False)


@mock.patch("perfdocs.logger.PerfDocLogger")
def test_perfdocs_verification(logger, structured_logger, perfdocs_sample):
    top_dir = perfdocs_sample["top_dir"]
    setup_sample_logger(logger, structured_logger, top_dir)

    from perfdocs.verifier import Verifier

    verifier = Verifier(top_dir)
    verifier.validate_tree()

    # Make sure that we had no warnings
    assert logger.warning.call_count == 0
    assert logger.log.call_count == 1
    assert len(logger.mock_calls) == 1


@mock.patch("perfdocs.logger.PerfDocLogger")
def test_perfdocs_framework_gatherers(logger, structured_logger, perfdocs_sample):
    top_dir = perfdocs_sample["top_dir"]
    setup_sample_logger(logger, structured_logger, top_dir)

    # Check to make sure that every single framework
    # gatherer that has been implemented produces a test list
    # in every suite that contains a test with an associated
    # manifest.
    from perfdocs.gatherer import frameworks

    for framework, gatherer in frameworks.items():
        with open(perfdocs_sample["config"], "w") as f:
            f.write(DYNAMIC_SAMPLE_CONFIG.format(framework))

        fg = gatherer(perfdocs_sample["config"], top_dir)
        if getattr(fg, "get_test_list", None) is None:
            # Skip framework gatherers that have not
            # implemented a method to build a test list.
            continue

        # Setup some framework-specific things here if needed
        if framework == "raptor":
            fg._manifest_path = perfdocs_sample["manifest"]
            fg._get_subtests_from_ini = mock.Mock()
            fg._get_subtests_from_ini.return_value = {
                "Example": perfdocs_sample["manifest"]
            }

        for suite, suitetests in fg.get_test_list().items():
            assert suite == "suite"
            for test, manifest in suitetests.items():
                assert test == "Example"
                assert manifest == perfdocs_sample["manifest"]


def test_perfdocs_logger_failure(config, paths):
    from perfdocs.logger import PerfDocLogger

    PerfDocLogger.LOGGER = None
    with pytest.raises(Exception):
        PerfDocLogger()

    PerfDocLogger.PATHS = []
    with pytest.raises(Exception):
        PerfDocLogger()


if __name__ == "__main__":
    mozunit.main()
