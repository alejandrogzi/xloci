# syntax=docker/dockerfile:1.7

FROM rust:1.93.0-bookworm as builder
WORKDIR /app

COPY xloci/Cargo.toml xloci/Cargo.lock ./
COPY xloci/src ./src

RUN cargo build --release

FROM debian:bookworm-slim
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /app/target/release/xloci /usr/local/bin/xloci

ENTRYPOINT ["xloci"]
