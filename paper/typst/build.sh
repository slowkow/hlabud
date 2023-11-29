#!/usr/bin/env bash
#
# Install rustup:
#
#     curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
#
# Then install rust and cargo with rustup.
#
# Then install typst:
#
#     cargo install --git https://github.com/typst/typst
#

# Build once
# typst compile main.typ

# Continously rebuild each time we make an edit
typst watch main.typ

# For diagrams, install d2:
#
#     go install oss.terrastruct.com/d2@latest
