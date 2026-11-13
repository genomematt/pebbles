import os
import pathlib

import pytest

from pebbles.pebbles import cli, main

DATA_DIR = pathlib.Path(__file__).parent.joinpath("data")
TEST_SAM = str(DATA_DIR.joinpath("map.sam"))
TEST_BAM = str(DATA_DIR.joinpath("map.bam"))

EXPECTED_COUNT = (
    'variant\tcount\n'
    'AY286018:g.16_18delGAC\t1\n'
    'AY286018:g.18_19insATG\t1\n'
    'AY286018:g.19_20delinsAG\t2\n'
    'AY286018:g.19_21delinsATG\t1\n'
    'AY286018:g.59A>T\t2\n'
)


@pytest.mark.parametrize("arguments", [
    f"count --max 2 --min_quality 5 {TEST_BAM}",
    ["count", "--max", "2", "--min_quality", "5", TEST_BAM],
])
def test_cli_uses_supplied_arguments(arguments):
    args = cli(arguments)
    assert args.command == 'count'
    assert args.max == 2
    assert args.min_quality == 5
    assert args.infile == TEST_BAM


def test_main_count(capsys):
    main(['count', TEST_BAM])
    assert capsys.readouterr().out == EXPECTED_COUNT


def test_missing_subcommand_is_a_usage_error(capsys):
    with pytest.raises(SystemExit) as exit_info:
        cli([])
    assert exit_info.value.code == 2
    assert 'required' in capsys.readouterr().err


def test_version(capsys):
    from pebbles import VERSION
    with pytest.raises(SystemExit) as exit_info:
        cli(['--version'])
    assert exit_info.value.code == 0
    assert capsys.readouterr().out == f'pebbles v{VERSION}\n'


@pytest.mark.parametrize("path", [TEST_SAM, TEST_BAM])
def test_format_detected_from_content(tmp_path, capsys, path):
    # extension does not match the content, so the format must come from the file itself
    renamed = tmp_path / "alignments.data"
    renamed.write_bytes(pathlib.Path(path).read_bytes())
    main(['count', str(renamed)])
    assert capsys.readouterr().out == EXPECTED_COUNT


@pytest.mark.parametrize("path", [TEST_SAM, TEST_BAM])
def test_reads_stdin(capfd, path):
    with open(path, 'rb') as handle:
        # htslib reads file descriptor 0 directly, so point it at the test file
        saved = os.dup(0)
        os.dup2(handle.fileno(), 0)
        try:
            main(['count', '-'])
        finally:
            os.dup2(saved, 0)
            os.close(saved)
    assert capfd.readouterr().out == EXPECTED_COUNT


@pytest.mark.parametrize("make_path", [
    lambda tmp_path: tmp_path / "missing.bam",
    lambda tmp_path: tmp_path,
    lambda tmp_path: (tmp_path / "junk.txt", (tmp_path / "junk.txt").write_text("not alignments\n"))[0],
])
def test_unreadable_input_is_a_clean_error(tmp_path, make_path):
    path = str(make_path(tmp_path))
    with pytest.raises(SystemExit) as exit_info:
        main(['count', path])
    assert f'cannot read alignments from {path}' in str(exit_info.value.code)
