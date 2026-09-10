import json
import logging
from pathlib import Path

import numpy as np
from tiled.client import from_uri
from tpx3awkward.processing import Tpx3Config, convert_tpx3_files_parallel

logger = logging.getLogger(__name__)

chx_catalog = from_uri("https://tiled.nsls2.bnl.gov")["chx"]["raw"]


def process_raw_files_sids(
    sids: list[str],
    output_dir: Path | None = None,
    tpx3_config_fpath: Path | None = None,
) -> None:
    """Convert the raw Timepix3 files associated with each scan ID.

    When ``output_dir`` is omitted, each scan is written beside its source
    data with the ``assets`` path component replaced by ``Compressed_Data``.
    """
    config = None
    if tpx3_config_fpath is not None:
        with tpx3_config_fpath.open() as config_file:
            config = Tpx3Config(**json.load(config_file))

    for sid in sids:
        try:
            raw_filepaths = chx_catalog[sid]["primary"]["data"][
                "tpx3_files_raw_filepaths"
            ].read()[0]
        except (KeyError, IndexError):
            logger.warning("Scan %s does not have Timepix3 files", sid)
            continue

        if not isinstance(raw_filepaths, np.ndarray):
            logger.warning("Scan %s has invalid Timepix3 file metadata", sid)
            continue

        filepaths = [str(filepath).removeprefix("file:") for filepath in raw_filepaths]
        if not filepaths:
            logger.warning("Scan %s does not have Timepix3 files", sid)
            continue

        scan_output_dir = output_dir
        if scan_output_dir is None:
            source_dir = Path(filepaths[0]).parent
            scan_output_dir = Path(str(source_dir).replace("assets", "Compressed_Data"))

        convert_tpx3_files_parallel(
            filepaths,
            output_dir=scan_output_dir,
            config=config,
        )
