# The image of CI's gcc job (.github/workflows/ci.yml), so a Codespace builds with the compiler CI tests.
FROM gcc:16

RUN apt-get update \
    && apt-get install -y --no-install-recommends cmake ninja-build ccache \
    && rm -rf /var/lib/apt/lists/*
