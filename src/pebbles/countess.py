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
from countess.core.plugins import PandasInputPlugin


class CountSAMPlugin(PandasInputPlugin):
    """Counts occurrences of alleles in a SAM file"""

    name = "SAM to Counts"
    title = "Count alleles in a SAM file (pebbles)"
    description = "Uses the Pebbles package to call variants from a SAM file to HGVS g. strings"
    version = VERSION

    file_types = [("SAM", "*.sam"),]
    file_mode = 'r'

    parameters = {
        "max": IntegerParam("Maximum number of variants in a valid allele (read/alignment)", 1),
        "min_quality": IntegerParam("Minimum quality score of alignment for a valid allele", 0),
    }


    def read_file_to_dataframe(self, file_param, logger, row_limit=None):
        records = {}
        count_column_name = "count"

        with pysam.AlignmentFile(file_param["filename"].value, self.file_mode) as fh:
            records = count_dict(fh,
                                 max_variants=self.parameters["max"].value,
                                 min_quality=self.parameters["min_quality"].value,
                                 row_limit=row_limit,
                                 logger=logger)

        return pd.DataFrame.from_records(
            list(records.items()),  columns=("allele", count_column_name)
            )


    def combine_dfs(self, dfs):
        """first concatenate the count dataframes, then group them by allele"""

        combined_df = pd.concat(dfs)

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
    file_mode = 'rb'

    parameters = {
        "max": IntegerParam("Maximum number of variants in a valid allele (read/alignment)", 1),
        "min_quality": IntegerParam("Minimum quality score of alignment for a valid allele", 0),
    }
