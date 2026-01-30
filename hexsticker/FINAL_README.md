# Final CU Boulder Hex Stickers

## Overview

**Professional, clean design with CU Boulder branding in black or transparent backgrounds.**

This is the final recommended design showcasing ARC's saddle point escape capability with:
- **Black or transparent background** for versatility
- **Greyscale contours** showing the optimization landscape
- **CU Gold thin trajectories** showing escape paths
- **CU Gold branding** (text and border)

## Design Specifications

### Visual Elements

✓ **Black or transparent background** - Modern, professional look
✓ **Greyscale contours** - White (smallest loss) to grey (saddle point)
✓ **Thin gold lines** - Clean trajectory visualization (0.8pt width)
✓ **No arrows, no points** - Minimalist, elegant appearance
✓ **2 symmetric paths** - Left and right of saddle
✓ **3 contour levels** - Saddle + 2 around minima

### Color Scheme

**Background**: Black (#000000)

**Contours** (greyscale gradient):
- **White** (f = -0.24) - Innermost contour, closest to minima
- **Grey 70%** (f = -0.22) - Middle contour around minima
- **Grey 40%** (f = 0) - Saddle point level

The greyscale gradient intuitively shows the loss landscape: brighter = better (lower loss).

**Branding**:
- **CU Gold** (#CFB87C) - Hex border, text, trajectories
- **Accessible CU Gold** (#8D7334) - Darker variant for increased contrast

### Trajectories

**Starting points** (symmetric, left and right of saddle):
1. **Path 1**: (-0.05, 0.1) → left minimum (-1, 0) in 4 iterations
2. **Path 2**: (0.05, 0.1) → right minimum (1, 0) in 4 iterations

**Style**:
- Line width: 0.8pt (thin, clean lines)
- No arrows, no points
- Minimalist appearance
- CU Gold color

## Available Variants

Four variants with different gold shades and backgrounds:

### Black Background Variants

#### 1. CU Gold (Primary) ⭐
**File**: `arcopt_final_gold_black.png` (50KB)
- Background: Black
- Trajectories, text, border: CU Gold (#CFB87C)
- **Recommended**: Bold, modern look for most uses

#### 2. Accessible CU Gold (High Contrast)
**File**: `arcopt_final_darkgold_black.png` (50KB)
- Background: Black
- Trajectories, text, border: Accessible CU Gold (#8D7334)
- Darker gold for maximum contrast
- Better for small sizes or low-contrast displays

### Transparent Background Variants

#### 3. CU Gold (Transparent) ⭐
**File**: `arcopt_final_gold_transparent.png` (50KB)
- Background: Transparent
- Trajectories, text, border: CU Gold (#CFB87C)
- **Recommended**: Adapts to any background color
- Ideal for websites, slides, and flexible usage

#### 4. Accessible CU Gold (Transparent)
**File**: `arcopt_final_darkgold_transparent.png` (50KB)
- Background: Transparent
- Trajectories, text, border: Accessible CU Gold (#8D7334)
- Darker gold with transparent background
- Maximum flexibility and contrast

## Design Philosophy

**"Less is more, but strategic"**

The final design strips away unnecessary elements while keeping everything needed to tell the story:

1. **Black or transparent background** - Versatile and focuses attention on the gold paths
2. **Greyscale contours** - Show landscape without overwhelming the visualization
3. **Thin lines** - Clean, elegant appearance without clutter
4. **No arrows, no points** - Minimalist aesthetic, path continuity is clear
5. **Symmetric layout** - Visually balanced, shows both directions from saddle
6. **CU Gold branding** - Consistent brand identity throughout

## Comparison with Previous Versions

| Feature | Original | Minimal (White) | Final |
|---------|----------|-----------------|-------|
| Background | Colorful fill | White | Black or transparent |
| Contours | Viridis (filled) | Gray/black lines | Greyscale lines |
| Paths | White | CU Gold | CU Gold |
| Path style | Thick w/ points | Thick w/ points | Thin, clean |
| Direction indicator | None | None | None |
| Text color | White | Black | CU Gold |
| Border color | Dark navy | CU Gold | CU Gold |
| Visual style | Academic | Clean & bright | Bold & modern |
| File size | 61KB | 50KB | 50KB |

## Visual Impact

The final design:
- **Stands out** with bold black or adapts with transparent background
- **Focuses attention** on the gold trajectories with clean lines
- **Looks modern** and professional with minimalist aesthetic
- **Works well** at small sizes (hex sticker scale ~2 inches)
- **High contrast** makes paths immediately visible
- **Intuitive** - brighter contours = better regions
- **Versatile** - transparent version works anywhere

The greyscale contours provide just enough context without competing with the main story (the escape trajectories). The thin lines without arrows create a clean, elegant look that's both professional and eye-catching.

## Recommendations

### Background Choice

**Black background** (`*_black.png`):
- Bold, modern, striking appearance
- High contrast makes paths stand out
- Best for light-colored surfaces (white pages, light slides)
- Makes a strong visual statement

**Transparent background** (`*_transparent.png`): ⭐
- Adapts to any background color
- More versatile and flexible
- Works on both light and dark backgrounds
- Recommended for most uses (websites, slides, READMEs)

### Color Choice

**Primary choice**: CU Gold (`*_gold_*.png`) ⭐

Use the standard CU Gold for most applications. The color is vibrant and matches official CU Boulder branding.

**Alternative**: Use Accessible CU Gold (`*_darkgold_*.png`) if:
- The sticker will be printed very small
- Viewing conditions have poor contrast
- You want maximum readability
- Using on a lighter background where regular gold doesn't pop enough

## Usage Guidelines

**Ideal for**:
- R package hex sticker (primary use case)
- README badges and shields
- Slide presentations (especially dark-themed)
- Academic posters
- Website headers
- GitHub social preview cards
- Stickers and swag

**Black background variants** work best on:
- White or light backgrounds (maximum contrast)
- Light grey backgrounds
- Clean, minimal surfaces

**Transparent background variants** work best on:
- Any solid color background (adapts automatically)
- Websites with varying themes
- Slides with different color schemes
- Documents where background might change

**All variants** - Avoid:
- Busy patterned backgrounds (competes for attention)
- Colors too similar to gold (low contrast)

## Technical Details

- **Resolution**: 300 dpi
- **Format**: PNG with transparency
- **Dimensions**: Standard hexSticker package dimensions
- **File size**: ~50KB
- **Color mode**: RGB
- **Contour levels**: 3 (f = -0.24, -0.22, 0)
- **Trajectory width**: 0.8pt
- **Style**: Thin lines, no arrows, no points

## Reproduction

```r
source("hexsticker/saddle_function.R")
source("hexsticker/create_final_sticker.R")
```

The script:
1. Loads arcopt package
2. Generates 2 symmetric paths from left/right of saddle
3. Creates greyscale contours (white=minima, grey=saddle)
4. Adds thin gold paths (no arrows, no points)
5. Sets black or transparent background with CU Gold branding
6. Saves as 300 dpi PNG files

## Mathematical Story

The visualization communicates:

1. **Saddle point topology** - Grey contour at f=0 passes through origin
2. **Negative curvature** - Paths diverge from near saddle in opposite directions
3. **Basin of attraction** - White/light grey contours show two minima basins
4. **Efficient convergence** - Clean paths reach minima in just 4 iterations each
5. **Robustness** - Symmetric starting points both succeed

This is the complete narrative of cubic regularization's advantage: detecting indefiniteness and efficiently escaping saddle points.
