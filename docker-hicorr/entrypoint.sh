#!/bin/bash
set -e

# Setup any runtime environment variables
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

# Execute the command passed to docker run
exec "$@"
