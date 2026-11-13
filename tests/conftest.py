import pytest


@pytest.fixture
def sam_factory(tmp_path):
    """Writes SAM records against a 10bp reference R (ACGTACGTAC) to a file and returns its path."""

    def _create_sam(records: list, filename: str = "records.sam"):
        path = tmp_path / filename
        lines = ["@SQ\tSN:R\tLN:10"] + ["\t".join(fields) for fields in records]
        path.write_text("\n".join(lines) + "\n")
        return str(path)

    return _create_sam
