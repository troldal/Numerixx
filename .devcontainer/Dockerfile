# The image of CI's gcc job (.github/workflows/ci.yml), so a Codespace builds with the compiler CI tests.
FROM gcc:16

# bash for the substitution in Doxygen's URL; pipefail so that a failed download fails its pipe.
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

RUN apt-get update \
    && apt-get install -y --no-install-recommends cmake ninja-build ccache \
    && rm -rf /var/lib/apt/lists/*

# Clang 22 and libc++, as CI's clang job, from the LLVM project's packages for this Debian release (trixie has 19).
# clang++.cfg makes libc++ the default, as CI's flags do, so that `cmake --preset clang` needs no flags here; without
# it Clang would take Debian's libstdc++ 14, not GCC 16's.
RUN curl -fsSL https://apt.llvm.org/llvm-snapshot.gpg.key | gpg --dearmor -o /usr/share/keyrings/llvm.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/llvm.gpg] https://apt.llvm.org/trixie/ llvm-toolchain-trixie-22 main" \
        > /etc/apt/sources.list.d/llvm.list \
    && apt-get update \
    && apt-get install -y --no-install-recommends clang-22 libc++-22-dev libc++abi-22-dev \
    && rm -rf /var/lib/apt/lists/* \
    && ln -s /usr/bin/clang-22 /usr/local/bin/clang \
    && ln -s /usr/bin/clang++-22 /usr/local/bin/clang++ \
    && echo "-stdlib=libc++" > /usr/lib/llvm-22/bin/clang++.cfg

# The documentation (docs/CMakeLists.txt): Doxygen 1.18 from its release, as CI's docs job (Debian has 1.9), and uv,
# which fetches Python and the pinned packages of docs/sphinx/requirements.txt on the first docs build.
ARG DOXYGEN_VERSION=1.18.0
ARG DOXYGEN_SHA256=14fa81bdc34171edb5f1f02b1d60e74802f0439b77fa44e592565d517d72df90
RUN curl -fsSL -o /tmp/doxygen.tar.gz \
        "https://github.com/doxygen/doxygen/releases/download/Release_${DOXYGEN_VERSION//./_}/doxygen-${DOXYGEN_VERSION}.linux.bin.tar.gz" \
    && echo "$DOXYGEN_SHA256  /tmp/doxygen.tar.gz" | sha256sum -c - \
    && tar xzf /tmp/doxygen.tar.gz -C /opt \
    && ln -s /opt/doxygen-${DOXYGEN_VERSION}/bin/doxygen /usr/local/bin/doxygen \
    && rm /tmp/doxygen.tar.gz
COPY --from=ghcr.io/astral-sh/uv:0.13.5 /uv /uvx /usr/local/bin/

# OpenAI's Codex CLI, for the independent review at the end of a phase (CLAUDE.md), with the bubblewrap from the same
# release, in which it sandboxes the commands it runs (`-s read-only`); devcontainer.json lets bubblewrap make its
# namespaces.
ARG CODEX_VERSION=0.161.0
ARG CODEX_SHA256=b1efb95097660d7f2e5a3887618a23f2ea1b0d548078bf92b0f7a5d229a0cef2
ARG BWRAP_SHA256=68f3b912c9d5f113343f6e054b5bb2be3cbbc09d21fcbc67ea6e997ed4b1cbfd
RUN cd /tmp \
    && release="https://github.com/openai/codex/releases/download/rust-v${CODEX_VERSION}" \
    && curl -fsSLO "$release/codex-x86_64-unknown-linux-musl.tar.gz" \
    && curl -fsSLO "$release/bwrap-x86_64-unknown-linux-musl.tar.gz" \
    && echo "$CODEX_SHA256  codex-x86_64-unknown-linux-musl.tar.gz" | sha256sum -c - \
    && echo "$BWRAP_SHA256  bwrap-x86_64-unknown-linux-musl.tar.gz" | sha256sum -c - \
    && tar xzf codex-x86_64-unknown-linux-musl.tar.gz \
    && tar xzf bwrap-x86_64-unknown-linux-musl.tar.gz \
    && install codex-x86_64-unknown-linux-musl /usr/local/bin/codex \
    && install bwrap-x86_64-unknown-linux-musl /usr/local/bin/bwrap \
    && rm /tmp/codex-x86_64-unknown-linux-musl* /tmp/bwrap-x86_64-unknown-linux-musl*

# The PDF manual (numerixx_docs_pdf): a minimal TeX Live with only what the manual's LaTeX loads, from TeX Live 2025's
# frozen archive, so the image does not change with the live repository. tlmgr checks each package's checksum against
# the archive's database, whose signature is valid but made with a key that has expired since: tlmgr warns, and goes on.
# The packages are the ones the manual reads (latexmk's .fls), found by compiling it; pdflatex embeds one URW Times
# font (utmr8a), hence times.
ARG TEXLIVE_REPOSITORY=https://ftp.math.utah.edu/pub/tex/historic/systems/texlive/2025/tlnet-final
ARG TEXLIVE_INSTALLER_SHA256=311df9f1477fd90c520159d1feddc2d6270f010d8349d1f6bdb9461a93b48a5c
COPY texlive.profile /tmp/texlive.profile
RUN mkdir /tmp/install-tl \
    && curl -fsSL -o /tmp/install-tl.tar.gz "$TEXLIVE_REPOSITORY/install-tl-unx.tar.gz" \
    && echo "$TEXLIVE_INSTALLER_SHA256  /tmp/install-tl.tar.gz" | sha256sum -c - \
    && tar xzf /tmp/install-tl.tar.gz -C /tmp/install-tl --strip-components=1 \
    && /tmp/install-tl/install-tl -profile /tmp/texlive.profile -repository "$TEXLIVE_REPOSITORY" \
    && /opt/texlive/bin/x86_64-linux/tlmgr install \
        latexmk \
        plex fontaxes mweights txfonts times \
        titlesec tabulary varwidth framed fancyvrb float wrapfig capt-of needspace upquote parskip \
        cmap etoolbox pict2e ellipse booktabs xcolor xkeyval \
    && rm -rf /tmp/install-tl /tmp/install-tl.tar.gz /tmp/texlive.profile
ENV PATH="/opt/texlive/bin/x86_64-linux:${PATH}"
