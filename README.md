# Modern Techniques in Modelling Github Repo

This repo contains the material to create the MTM course page, at: <https://nicholasdavies.github.io/mtm>

All the raw material is in the top level folder and is written in Quarto.

To modify the published web site, you have to first run

``` r
quarto::quarto_render()
```

after updating the .qmd files, then commit and push changes to the git repo.

To make non-Quarto files available, just put them in a subdirectory (such as `slides/` or `code/`, but NOT `docs/`) and link to them from a qmd file. Quarto will detect the link and move it to the appropriate place automatically. Link to it just as e.g. `slides/Lecture_01_Intro.pdf`.

The actual rendered website is in `docs/`. Do not change anything in there manually.
