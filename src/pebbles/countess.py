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
from countess.utils.dask import concat_dataframes as concat_dask_dataframes


class CountSAMPlugin(DaskInputPlugin):
    """Counts occurrences of alleles in a SAM file"""

    name = "SAM to Counts"
    title = "Count alleles in a SAM file (pebbles)"
    description = "Uses the Pebbles package to call variants from a SAM file to HGVS g. strings"
    version = VERSION

    file_types = [("SAM", "*.sam"),]

    parameters = {
        "max": IntegerParam("Maximum number of variants in a valid allele (read/alignment)", 1),
    }

    _file_permissions = 'r'

    def read_file_to_dataframe(self, file_param, column_suffix="", row_limit=None):
        records = {}
        count_column_name = "count"
        if column_suffix:
            count_column_name += "_" + str(column_suffix)

        with pysam.AlignmentFile(file_param["filename"].value, self._file_permissions) as fh:
            records = count_dict(fh, self.parameters["max"].value, row_limit)

        if records:
            return pd.DataFrame.from_records(
                list(records.items()),  columns=("allele", count_column_name)
            )
        else:
            # records may be an empty dictionary if no variants in row_limit
            # return a useful shaped object to use in preview
            return pd.DataFrame.from_records(
                [(f'Warning:no_mutations_first_{row_limit}_alignments',0)],  columns=("allele", count_column_name)
            )


    def combine_dfs(self, dfs):
        """first concatenate the count dataframes, then group them by allele"""

        combined_df = concat_dask_dataframes(dfs)

        if len(combined_df):
            combined_df = combined_df.groupby(by=["allele"]).sum()

        return combined_df

class CountBAMPlugin(CountSAMPlugin):
    """Counts occurrences of alleles in a BAM file"""

    name = "BAM to Counts"
    title = "Count alleles in a BAM file (pebbles)"
    description = "Uses the Pebbles package to call variants from a BAM file to HGVS g. strings"
    version = VERSION

    file_types = [("BAM", "*.bam"),]

    parameters = {
        "max": IntegerParam("Maximum number of variants in a valid allele (read/alignment)", 1),
    }

    _file_permissions = 'rb'

