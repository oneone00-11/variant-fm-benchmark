"""E7: the pipeline from a given stage on, in the reference image, to write the commit's
outputs without running the stages before it, whose outputs the commit leaves unchanged
(a copy of scripts/stage_E6/pipeline_from.py for this folder's run records). E7 runs it
from stage 17, the predictors' training provenance, whose fusion row's text it corrects,
and from stage 41, the layout. Only a full run decides: the commit is then checked by
scripts/stage_E7/check_scope.py against the previous commit and by
`container/run.sh verify`, which runs every stage.

    docker run ... IMAGE bash scripts/stage_E7/in_image.sh scripts/stage_E7/pipeline_from.py 17
"""
from __future__ import annotations

import runpy
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _runlog import REPO, RunLog  # noqa: E402


def main() -> None:
    start = int(sys.argv[1])
    with RunLog(__file__) as log:
        log.note("from_stage", start)
        sys.argv = [str(REPO / "scripts/reproduce_evidence.py"), "--from", str(start)]
        try:
            runpy.run_path(sys.argv[0], run_name="__main__")
        except SystemExit as e:
            if e.code not in (0, None):
                raise


if __name__ == "__main__":
    main()
