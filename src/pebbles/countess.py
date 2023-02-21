import pandas as pd  # type: ignore

import pysam
from pebbles import VERSION
from pebbles.pebbles import count_dict

from countess.core.parameters import (
    ArrayParam,
    BooleanParam,
    FileArrayParam,
    FileParam,
    FloatParam,
    MultiParam,
    StringParam,
    IntegerParam,
)
from countess.core.plugins import DaskInputPlugin


class CountSAMPlugin(DaskInputPlugin):
    """Counts occurences of alleles in a SAM file"""

    name = "SAM to Counts"
    title = "Count alleles in a SAM file (pebbles)"
    description = "Uses the Pebbles package to call variants from a SAM file to HGVS g. strings"
    version = VERSION

    file_types = [("SAM", "*.sam"),]

    parameters = {
        "max": IntegerParam("Maximum number of variants in a valid allele (read/alignment)", 1),
    }

    _file_permissions = 'r'

    def read_file_to_dataframe(self, file_param, logger, column_suffix="", row_limit=None):
        records = []
        count_column_name = "count"
        if column_suffix:
            count_column_name += "_" + str(column_suffix)

        with pysam.AlignmentFile(file_param["filename"].value, _file_permissions) as fh:
            records = count_dict(fh, self.parameters["max"].value)

        return pd.DataFrame.from_dict(
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

