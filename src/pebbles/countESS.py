from collections.abc import Iterable, Mapping
from itertools import islice
from typing import Generator, Optional

import dask.dataframe as dd
import numpy as np
import pandas as pd  # type: ignore
from more_itertools import ichunked

from pebbles import VERSION
from pysam import AlignmentFile

from countess.core.parameters import (
    ArrayParam,
    BooleanParam,
    FileArrayParam,
    FileParam,
    FloatParam,
    MultiParam,
    StringParam,
)
from countess.core.plugins import DaskInputPlugin
from countess.utils.dask import concat_dataframes, merge_dataframes



class CountSAMPlugin(DaskInputPlugin):
    """Counts occurences of alleles in a SAM file"""

    name = "SAM to Counts"
    title = "Load from FastQ"
    description = "Uses the Pebbles package to call variants from a SAM file to HGVS g. strings"
    version = VERSION

    file_types = [("SAM", "*.sam"),]

    parameters = {
        "max": BooleanParam("Maximum number of variants in a valid allele (read/alignment)", 1),
    }

    _file_permissions = 'r'

    def count_file_to_dataframe(self, file_param, column_suffix="", row_limit=None, _file_permissions='r'):
        records = []
        count_column_name = "count"
        if column_suffix:
            count_column_name += "_" + column_suffix

        with pysam.AlignmentFile(file_param["filename"].value, "_file_permissions") as fh:
            records = pebbles.count_dict(fh, self.parameters["max"].value)

        return pd.DataFrame.from_records(
            records, columns=("allele", count_column_name)
        )

    def combine_dfs(self, dfs):
        """first concatenate the count dataframes, then (optionally) group them by sequence"""

        combined_df = concat_dask_dataframes(dfs)

        if len(combined_df) and self.parameters["group"].value:
            combined_df = combined_df.groupby(by=["allele"]).sum()

        return combined_df

class CountBAMPlugin(CountSAMPlugin):
    """Counts occurences of alleles in a SAM file"""

    name = "BAM to Counts"
    title = "Count from BAM"
    description = "Uses the Pebbles package to call variants from a BAM file to HGVS g. strings"
    version = VERSION

    file_types = [("BAM", "*.bam"),]

    parameters = {
        "max": BooleanParam("Maximum number of variants in a valid allele (read/alignment)", 1),
    }

    _file_permissions = 'rb'

