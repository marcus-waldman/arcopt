# arcopt Hex Sticker

## Concept

A stylized contour plot showing ARC escaping a saddle point -- the signature capability of cubic regularization. The visual story: "Newton gets stuck at saddle points; ARC detects negative curvature and escapes to a minimum."

## Test Function

Use a double-well potential with a saddle point at the origin:

```
f(x, y) = x^4/4 - x^2/2 + y^2/2
```

**Properties:**

- Gradient: `(x^3 - x, y)`
- Hessian: `diag(3x^2 - 1, 1)`
- **Saddle point at (0, 0)** with Hessian `diag(-1, 1)` -- genuinely indefinite
- Two minima at `(-1, 0)` and `(1, 0)`
- Contour structure: two basins separated by a saddle ridge

## Design Specification

### Layout

- **Background**: contour-filled plot of `f(x, y)` over the region around the saddle and both minima
- **Optimization path**: bold curve starting near the saddle (e.g., `(0.05, 0.3)`) that escapes along the negative curvature direction and converges to one of the minima, with small dots at each iterate
- **Saddle marker**: small marker (star, diamond, or x) at the origin to mark the saddle point
- **Text**: "arcopt" in white, clean sans-serif font at the bottom of the hex

### Color Palette (Colorblind-Friendly)

- Use `viridis` or `cividis` for contour fill (`cividis` is specifically optimized for deuteranopia/protanopia)
- Hex border: deep navy or charcoal
- Optimization path: bold white or bright yellow
- Iterate dots: white with thin dark outline
- Text: white

### Hex Sticker Parameters

- Standard hex sticker dimensions (use `hexSticker` package defaults)
- Dark border color matching the palette
- Minimal/no subplot padding so the contour fills the hex nicely

## Implementation Tasks

- [ ] Install required packages (`hexSticker`, `ggplot2`, `viridis`, `showtext`)
- [ ] Define `f(x, y)`, gradient, and Hessian in R
- [ ] Generate the optimization path by running `arcopt()` from a starting point near the saddle (e.g., `x0 = c(0.05, 0.3)`) and recording iterates
- [ ] If `arcopt()` iteration recording isn't available, manually simulate the path using the ARC algorithm or hard-code a representative path
- [ ] Create the contour subplot in `ggplot2`:
  - [ ] Evaluate `f` on a grid (e.g., `x in [-1.8, 1.8]`, `y in [-1.2, 1.2]`)
  - [ ] Use `geom_contour_filled()` with `viridis` or `cividis` scale
  - [ ] Overlay path with `geom_path()` (white/yellow, linewidth ~1)
  - [ ] Add iterate points with `geom_point()`
  - [ ] Add saddle marker at origin
  - [ ] Strip all axes, labels, legends, and background (theme_void)
- [ ] Build the hex sticker with `hexSticker::sticker()`:
  - [ ] Feed in the ggplot subplot
  - [ ] Set package name `"arcopt"`, white text, clean font
  - [ ] Set hex fill/border colors (dark navy/charcoal)
  - [ ] Adjust `s_x`, `s_y`, `s_width`, `s_height` so the contour fills the hex well
  - [ ] Adjust `p_x`, `p_y`, `p_size` for text placement
- [ ] Save output as both PNG (300 dpi) and SVG
- [ ] Iterate on aesthetics:
  - [ ] Try `viridis` vs `cividis` vs `mako` and pick the best looking one
  - [ ] Adjust number of contour levels for visual clarity at small size
  - [ ] Ensure readability when printed at ~2 inch (5 cm) width
  - [ ] Consider reducing contour detail if it looks too busy at sticker scale

## Reference Files

- Existing saddle test: `tests/testthat/test-arcopt.R` (lines 151-170)
- Benchmark functions: `benchmarks/` directory
- Main optimizer: `R/arcopt.R`

## R Snippet to Get Started

```r
# install.packages(c("hexSticker", "ggplot2", "viridis", "showtext"))
library(ggplot2)
library(viridis)
library(hexSticker)

# --- Test function ---
f <- function(x, y) x^4 / 4 - x^2 / 2 + y^2 / 2

# --- Grid for contour ---
grid <- expand.grid(
  x = seq(-1.8, 1.8, length.out = 300),
  y = seq(-1.2, 1.2, length.out = 300)
)
grid$z <- f(grid$x, grid$y)

# --- Placeholder path (replace with actual arcopt iterates) ---
# This is a representative trajectory escaping the saddle at (0,0)
# toward the minimum at (1, 0)
path <- data.frame(
  x = c(0.05, 0.12, 0.35, 0.62, 0.85, 0.95, 0.99, 1.0),
  y = c(0.30, 0.22, 0.14, 0.08, 0.03, 0.01, 0.00, 0.0)
)

# --- Contour subplot ---
p <- ggplot(grid, aes(x = x, y = y, z = z)) +
  geom_contour_filled(bins = 15) +
  scale_fill_viridis_d(option = "cividis") +
  geom_path(
    data = path, aes(x = x, y = y, z = NULL),
    color = "white", linewidth = 1, lineend = "round"
  ) +
  geom_point(
    data = path, aes(x = x, y = y, z = NULL),
    color = "white", size = 1.5
  ) +
  annotate("point", x = 0, y = 0, shape = 4, color = "white", size = 3, stroke = 1.5) +
  theme_void() +
  theme(legend.position = "none")

# --- Hex sticker ---
sticker(
  p,
  package = "arcopt",
  p_size = 20,
  p_color = "white",
  p_y = 1.4,
  s_x = 1,
  s_y = 0.85,
  s_width = 1.6,
  s_height = 1.2,
  h_fill = "#1a1a2e",
  h_color = "#16213e",
  filename = "hexsticker/arcopt_sticker.png",
  dpi = 300
)
```

## Notes

- The placeholder path above is approximate. Ideally, run `arcopt()` on this function and record the actual iterates for authenticity.
- If the contour looks too busy at sticker size, reduce `bins` in `geom_contour_filled()` to 8-10.
- The function `f(x,y) = x^4/4 - x^2/2 + y^2/2` was chosen because it has a genuinely indefinite Hessian at the saddle (eigenvalues -1 and +1), making it a real test of ARC's negative curvature detection -- not just a cosmetic choice.
- For font customization, `showtext` + Google Fonts (e.g., "Fira Sans", "Roboto") can give a cleaner look than R's default fonts.
