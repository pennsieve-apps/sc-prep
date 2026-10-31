"""
One image, two runtimes. Fixes the thread count (processor/threads.py) before numpy is imported;
then, as a Lambda function (AWS_LAMBDA_RUNTIME_API is set), hands over to the Lambda runtime
interface client with processor.handler, or, as a container, runs main() once.
"""
import os
import sys

from processor import threads  # noqa: F401

if os.environ.get("AWS_LAMBDA_RUNTIME_API"):
    os.execvp(sys.executable, [sys.executable, "-m", "awslambdaric", "processor.handler.handler"])

from processor.main import main  # noqa: E402

sys.exit(main())
