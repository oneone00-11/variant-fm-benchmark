# PNG figures, pixel by pixel

Written by `scripts/stage_B6/compare_png_pixels.py`; do not edit by hand.

The 13 PNG figures the comparator called "differs (figure)" when every stage ran in the reference image (commit 7125246) against the commit before it (e12220c), each decoded and compared pixel by pixel with its committed version. `docs/rebaseline-report.md`, which lists them, is left as it was written.

PNG figures, pixel by pixel: old = commit e12220c, new = commit 7125246.

bytes identical: 0; bytes differ, pixels identical: 13; pixels differ: 0

| file | status | bytes old -> new | image | where pixels differ |
|---|---|---|---|---|
| `phase1/reports/evidence/fig_data/fig1_draft.png` | bytes differ, pixels identical | 30,423 -> 35,063 | 990x440 RGBA |  |
| `phase1/reports/evidence/fig_data/fig2_draft.png` | bytes differ, pixels identical | 75,014 -> 78,999 | 1430x396 RGBA |  |
| `phase1/reports/evidence/fig_data/fig3_draft.png` | bytes differ, pixels identical | 40,347 -> 46,463 | 990x440 RGBA |  |
| `phase1/reports/evidence/fig_data/fig4_draft.png` | bytes differ, pixels identical | 39,926 -> 48,461 | 990x440 RGBA |  |
| `phase1/reports/evidence/figures/figure1.png` | bytes differ, pixels identical | 175,428 -> 197,346 | 2161x1393 RGBA |  |
| `phase1/reports/evidence/figures/figure2.png` | bytes differ, pixels identical | 351,525 -> 372,748 | 2161x1771 RGBA |  |
| `phase1/reports/evidence/figures/figure3.png` | bytes differ, pixels identical | 338,290 -> 378,129 | 2161x1771 RGBA |  |
| `phase1/reports/evidence/figures/figure4.png` | bytes differ, pixels identical | 240,110 -> 270,564 | 2161x1866 RGBA |  |
| `phase1/reports/evidence/figures/figureS1.png` | bytes differ, pixels identical | 589,774 -> 628,632 | 2161x2952 RGBA |  |
| `phase1/reports/evidence/figures/figureS2.png` | bytes differ, pixels identical | 166,707 -> 177,746 | 2161x1240 RGBA |  |
| `phase1/reports/evidence/figures/figureS3.png` | bytes differ, pixels identical | 148,711 -> 164,328 | 2161x1122 RGBA |  |
| `phase1/reports/evidence/figures/figure_practical.png` | bytes differ, pixels identical | 330,780 -> 349,609 | 2161x1511 RGBA |  |
| `phase1/reports/evidence/figures/figure_practical_slides.png` | bytes differ, pixels identical | 293,996 -> 316,791 | 2000x1125 RGBA |  |
