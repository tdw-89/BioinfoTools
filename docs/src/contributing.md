```@meta
CurrentModule = BioinfoTools
```

```@id contributing
```

# Contributing & support

Whether you are opening an issue, contributing code, or extending the documentation, the guidelines
below will help keep the project consistent and easy to maintain.

## Getting help

- Open a GitHub issue for bugs, feature requests, or documentation gaps.
- Include the Julia version, BioinfoTools commit, relevant snippets, and (if possible) a reduced test
  case using public data.
- Discussions about new APIs or design changes are welcome—draft your proposal with motivation,
  alternatives considered, and expected inputs/outputs.

## Development workflow

1. Fork the repository and create a feature branch.
2. Activate the project, develop locally, and run the tests:

   ```bash
   julia --project test/runtests.jl
   ```

3. Update or add docstrings for all public APIs touched by your change.
4. Update this documentation when the user-facing behavior changes.

!!! tip "Style"
    Follow the Julia style guide, keep functions short, and prefer explicit type annotations where
    they communicate biological intent (e.g., genomic coordinates, strands, or counts).

## Documentation

- Keep `docs/Project.toml` and `docs/Manifest.toml` in sync with the code base.
- Use short sections, callouts (`!!! note` / `!!! warning`), and doctested snippets to keep examples
  trustworthy.
- Run `julia --project=docs -e 'import Pkg; Pkg.instantiate(); include("docs/make.jl")'` before
  opening a pull request to ensure the static site still builds cleanly.

## Release checklist

1. Bump the version in `Project.toml`.
2. Tag the release and push it upstream.
3. Build the docs and deploy them to your chosen static host (GitHub Pages, Netlify, etc.).
4. Announce notable changes in the README and release notes.
