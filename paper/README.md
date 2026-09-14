# How to reproduce the figures and paper

## Prerequisites

- [`Julia`](https://julialang.org/) version `1.13`
- [Docker](https://www.docker.com/)

## Figures

```bash
bash paper/script/make_figures.sh
```

This instantiates the project environment and runs every script in `paper/script/`, writing the figures to `paper/image/`.

## Paper

Start Docker, then run

```bash
bash paper/script/make_paper.sh
```
