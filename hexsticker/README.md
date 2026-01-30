# arcopt Hex Sticker

Final CU Boulder branded hex stickers for the arcopt R package.

## Files

### Stickers (300 dpi PNG)

**Black Background:**
- `arcopt_final_gold_black.png` - Standard CU Gold (#CFB87C)
- `arcopt_final_darkgold_black.png` - Accessible CU Gold (#8D7334)

**Transparent Background:**
- `arcopt_final_gold_transparent.png` - Standard CU Gold (#CFB87C) ⭐ Recommended
- `arcopt_final_darkgold_transparent.png` - Accessible CU Gold (#8D7334)

### Scripts

- `create_final_sticker.R` - Main script to generate all sticker variants
- `saddle_function.R` - Double-well potential test function definition

### Documentation

- `FINAL_README.md` - Complete design documentation
- `TASKLIST.md` - Original task list and design concept

## Design Features

- **2 symmetric paths** escaping from near saddle point at (0, 0)
- **3 greyscale contours** (white=minima, grey=saddle)
- **Thin gold lines** (0.8pt, no arrows, no points)
- **CU Gold branding** (text and hex border)
- **Black or transparent background** for versatility

## Quick Start

To regenerate the stickers:

```r
source("hexsticker/create_final_sticker.R")
```

This will create all 4 variants (2 black background, 2 transparent background).

## Recommendations

**Primary choice**: `arcopt_final_gold_transparent.png` ⭐

The transparent background version is the most versatile - it works on any background color and is ideal for websites, slides, READMEs, and documentation.

For full design details, see [FINAL_README.md](FINAL_README.md).
