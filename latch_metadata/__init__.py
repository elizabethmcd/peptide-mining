
from latch.types.metadata import (
    NextflowMetadata,
    LatchAuthor,
    NextflowRuntimeResources
)
from latch.types.directory import LatchDir

from .parameters import generated_parameters

NextflowMetadata(
    display_name='peptide-mining',
    author=LatchAuthor(
        name="Elizabeth McDaniel, Taylor Reiter",
    ),
    parameters=generated_parameters,
    runtime_resources=NextflowRuntimeResources(
        cpus=12,
        memory=36,
        storage_gib=100,
    ),
    log_dir=LatchDir("latch:///your_log_dir"),
)
