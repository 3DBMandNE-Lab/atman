# syntax=docker/dockerfile:1.6
#
# Atman reproducible build.
#
# Stage 1 builds the Rust workspace against a pinned toolchain.
# Stage 2 is a slim runtime carrying only the compiled binary and the
# example data required to reproduce the Dube pipeline end-to-end.
#
# To run the bundled Dube example:
#   docker build -t atman:1.0.0 .
#   docker run --rm -v "$(pwd)/out:/out" atman:1.0.0 \
#     atman ingest \
#       --platform olink-explore-ngs --parser dube \
#       --output-dir /out \
#       /data/dube_heat_2023/20212016_Dube_NPX_2021-11-30.csv \
#       /data/dube_heat_2023/20212017_Dube_NPX_2021-12-13_OID30253_corrected.csv

FROM rust:1.75-slim-bookworm AS builder

WORKDIR /src

# Cache-friendly dependency layer: resolve the dependency graph with minimal
# temporary source files, then overlay real sources.
COPY Cargo.toml Cargo.lock ./
COPY crates/atman-core/Cargo.toml crates/atman-core/Cargo.toml
COPY crates/atman/Cargo.toml crates/atman/Cargo.toml
RUN mkdir -p crates/atman-core/src crates/atman/src \
    && echo 'fn main() {}' > crates/atman/src/main.rs \
    && echo '' > crates/atman-core/src/lib.rs \
    && cargo build --workspace --release \
    && rm -rf crates/atman-core/src crates/atman/src

COPY crates crates
RUN cargo build --workspace --release \
    && strip target/release/atman

FROM debian:bookworm-slim AS runtime

RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /src/target/release/atman /usr/local/bin/atman
COPY example_data /data

RUN useradd --create-home --shell /bin/bash runner
USER runner
WORKDIR /home/runner

ENTRYPOINT ["atman"]
CMD ["--help"]
