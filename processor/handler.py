"""
The Lambda runtime. A lambda step's payload carries inputDir, outputDir and the step's parameters
as fields (where a container gets them as environment variables); the handler sets those
variables and runs the same main(), so both runtimes write the same outputs. The thread count is
fixed when the function starts (processor/threads.py: the function's vCPUs, about 2 at 3008 MB);
an event's THREADS does not change it.
"""
import json
import os

from processor import main as app


def handler(event, context):
    os.environ["INPUT_DIR"] = event.get("inputDir") or ""
    os.environ["OUTPUT_DIR"] = event.get("outputDir") or ""
    # A warm function keeps its environment between invocations, so a parameter this run doesn't set is
    # cleared, not left from the last one.
    for name in app.PARAMETERS:
        value = event.get(name)
        if value is None:
            os.environ.pop(name, None)
        elif isinstance(value, bool):
            os.environ[name] = str(value).lower()
        elif isinstance(value, (dict, list)):
            os.environ[name] = json.dumps(value) if isinstance(value, dict) else ",".join(map(str, value))
        else:
            os.environ[name] = str(value)
    code = app.main()
    if code != 0:
        raise RuntimeError(f"sc-prep exited with {code}")      # an exception fails the step; a return value would not
    return {"status": "ok", "outputDir": os.environ["OUTPUT_DIR"]}
