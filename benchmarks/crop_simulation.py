#!/usr/bin/env python3
"""Split the native-vector simulation g/h SVG into two standalone panels."""

from __future__ import annotations

import argparse
import re
import xml.etree.ElementTree as ET
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    raw = args.source.read_text(encoding="utf-8")
    match = re.search(r'<svg\b[^>]*\bviewBox="0 0 ([0-9.]+) ([0-9.]+)"', raw)
    if not match:
        raise ValueError("Expected a full-size Matplotlib SVG viewBox")
    source_width, source_height = map(float, match.groups())
    panels = (("simulation_time.svg", 0, 34, 575, 382),
              ("simulation_memory.svg", 645, 34, min(534, int(source_width - 645)), 382))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    original = re.compile(r'(<svg\b[^>]*?)width="[^"]+" height="[^"]+" viewBox="[^"]+"')
    for name, x, y, width, height in panels:
        if not (width > 0 and x + width <= source_width and y + height <= source_height):
            raise ValueError(f"Crop is outside source SVG: {name}")
        replacement = f'width="{width}pt" height="{height}pt" viewBox="{x} {y} {width} {height}"'
        cropped, count = original.subn(lambda found: found.group(1) + replacement, raw, count=1)
        if count != 1:
            raise ValueError(f"Could not rewrite SVG viewBox: {name}")
        output = args.output_dir / name
        output.write_text(cropped, encoding="utf-8")
        tree = ET.parse(output)
        if tree.findall('.//{http://www.w3.org/2000/svg}image'):
            raise ValueError(f"Raster image embedded in {name}")
        print(output)


if __name__ == "__main__":
    main()
