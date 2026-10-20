"""Build a descriptive W2-to-W3 matching flow from aggregate counts only."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from xml.sax.saxutils import escape

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "source_data/figS5_followup_flow.csv"
OUT = ROOT / "figures/FigureS5_followup_support_flow"
SCALE = 3

def font_file(bold: bool = False) -> str:
    names = (
        ["/System/Library/Fonts/Supplemental/Arial Bold.ttf",
         "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"]
        if bold else
        ["/System/Library/Fonts/Supplemental/Arial.ttf",
         "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"]
    )
    return next((name for name in names if Path(name).is_file()), "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf")


def main() -> None:
    with SOURCE.open(newline="", encoding="utf-8") as f:
        rows = {r["stage"]: r for r in csv.DictReader(f)}
    n = {k: int(v["n"]) for k, v in rows.items()}
    assert n["wave2_baseline"] == n["baseline_grades_without_match"] + n["continuing_grades"]
    assert n["continuing_grades"] == n["one_unsupported_school"] + n["supported_cells"]
    assert n["supported_cells"] == n["within_support_unmatched"] + n["matched_wave3"]

    boxes = [
        (350, 35, 300, 90, "Wave-2 baseline: 8,902", "11 school codes", "#e8efff", "#374151"),
        (40, 195, 360, 90, "2,961: grades 6, 9, 12", "No matched wave-3 outcome", "#f2f2f2", "#374151"),
        (510, 195, 450, 90, "5,941: continuing grades", "Grades 4, 5, 7, 8, 10, 11", "#e8efff", "#374151"),
        (40, 355, 360, 90, "1,422: one high school", "Grades 10–11; no matched outcome", "#f2f2f2", "#374151"),
        (510, 355, 450, 90, "4,519: observed-support domain", "Ten schools with any matched outcome", "#e8efff", "#374151"),
        (500, 500, 210, 80, "4,219 matched", "93.4% of 4,519", "#d9f1e7", "#267355"),
        (745, 500, 210, 80, "300 unmatched", "6.6% of 4,519", "#f2f2f2", "#374151"),
    ]
    arrows = [((500, 125), (220, 195)), ((500, 125), (735, 195)),
              ((735, 285), (220, 355)), ((735, 285), (735, 355)),
              ((735, 445), (605, 500)), ((735, 445), (850, 500))]
    notes = [
        (610, "Absence of a match does not identify why wave-3 data are missing."),
        (635, "The 93.4% denominator is defined by observed school × grade support."),
    ]

    img = Image.new("RGB", (1000 * SCALE, 650 * SCALE), "white")
    draw = ImageDraw.Draw(img)
    regular = ImageFont.truetype(font_file(), 16 * SCALE)
    small = ImageFont.truetype(font_file(), 14 * SCALE)
    bold = ImageFont.truetype(font_file(bold=True), 18 * SCALE)
    for (sx, sy), (ex, ey) in arrows:
        sx, sy, ex, ey = (v * SCALE for v in (sx, sy, ex, ey))
        draw.line((sx, sy, ex, ey), fill="#64748b", width=3 * SCALE)
        theta = math.atan2(ey - sy, ex - sx)
        length, half = 11 * SCALE, 6 * SCALE
        bx, by = ex - length * math.cos(theta), ey - length * math.sin(theta)
        left = (bx + half * math.sin(theta), by - half * math.cos(theta))
        right = (bx - half * math.sin(theta), by + half * math.cos(theta))
        draw.polygon([(ex, ey), left, right], fill="#64748b")
    for x, y, w, h, title, subtitle, fill, edge in boxes:
        draw.rounded_rectangle(tuple(v * SCALE for v in (x, y, x+w, y+h)), radius=13*SCALE, fill=fill, outline=edge, width=2*SCALE)
        draw.text(((x+w/2)*SCALE, (y+31)*SCALE), title, font=bold, fill="#111827", anchor="mm")
        draw.text(((x+w/2)*SCALE, (y+65)*SCALE), subtitle, font=regular, fill="#374151", anchor="mm")
    for y, note in notes:
        draw.text((500*SCALE, y*SCALE), note, font=small, fill="#4b5563", anchor="mm")
    img.save(OUT.with_suffix(".png"), dpi=(300, 300))
    img.save(OUT.with_suffix(".pdf"), resolution=300)

    svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="1000" height="650" viewBox="0 0 1000 650">', '<rect width="1000" height="650" fill="white"/>']
    for (sx, sy), (ex, ey) in arrows:
        theta = math.atan2(ey - sy, ex - sx)
        bx, by = ex - 11*math.cos(theta), ey - 11*math.sin(theta)
        left = (bx + 6*math.sin(theta), by - 6*math.cos(theta))
        right = (bx - 6*math.sin(theta), by + 6*math.cos(theta))
        svg.append(f'<path d="M {sx} {sy} L {ex} {ey}" stroke="#64748b" stroke-width="3"/>')
        svg.append(f'<polygon points="{ex},{ey} {left[0]:.2f},{left[1]:.2f} {right[0]:.2f},{right[1]:.2f}" fill="#64748b"/>')
    for x, y, w, h, title, subtitle, fill, edge in boxes:
        svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="13" fill="{fill}" stroke="{edge}" stroke-width="2"/>')
        svg.append(f'<text x="{x+w/2}" y="{y+37}" text-anchor="middle" font-family="Arial" font-size="18" font-weight="bold" fill="#111827">{escape(title)}</text>')
        svg.append(f'<text x="{x+w/2}" y="{y+70}" text-anchor="middle" font-family="Arial" font-size="16" fill="#374151">{escape(subtitle)}</text>')
    for y, note in notes:
        svg.append(f'<text x="500" y="{y+5}" text-anchor="middle" font-family="Arial" font-size="14" fill="#4b5563">{escape(note)}</text>')
    OUT.with_suffix(".svg").write_text("\n".join(svg + ["</svg>"]) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
