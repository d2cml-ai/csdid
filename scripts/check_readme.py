#!/usr/bin/env python3
"""Reject known README generation regressions before the bot opens a PR."""

import argparse
from html import unescape
from html.parser import HTMLParser
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import parse_qs, urlsplit


REQUIRED_SECTIONS = (
    "Getting Started",
    "Basic Example",
    "Event studies",
    "Overall effect of participating in treatment",
    "Additional options",
    "How to cite",
)
REQUIRED_BADGES = {
    "PyPI version": "https://img.shields.io/pypi/v/csdid.svg",
    "Downloads": "https://static.pepy.tech/personalized-badge/csdid?",
    "Last commit": "https://img.shields.io/github/last-commit/d2cml-ai/csdid.svg",
    "GitHub stars": "https://img.shields.io/github/stars/d2cml-ai/csdid.svg",
    "GitHub issues": "https://img.shields.io/github/issues/d2cml-ai/csdid.svg",
    "License": "https://img.shields.io/github/license/d2cml-ai/csdid.svg",
}


class _ImageParser(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.images = []

    def handle_starttag(self, tag, attrs):
        if tag.casefold() == "img":
            attributes = dict(attrs)
            self.images.append((attributes.get("alt", ""), attributes.get("src", "")))


def image_references(content):
    """Read Markdown images and HTML images, including escaped query strings."""
    parser = _ImageParser()
    parser.feed(content)
    references = re.findall(r"!\[([^\]]*)\]\(([^\s)]+)\)", content)
    references.extend(parser.images)
    return [(" ".join((alt or "").split()), unescape(source or "")) for alt, source in references]


def check_document(content, name):
    """Check authored content that the old bot source silently discarded."""
    errors = []
    normalized = " ".join(content.split())
    if not re.search(r"\*\*csdid\*\* package implements", normalized):
        errors.append(f"{name}: the Python package must be identified as csdid")
    if re.search(r"\*\*did\*\* package", normalized):
        errors.append(f"{name}: the Python package was renamed to did")

    headings = {
        " ".join(match.split()).casefold()
        for match in re.findall(r"^#{1,6}\s+(.+?)\s*#*\s*$", content, re.MULTILINE)
    }
    for heading in REQUIRED_SECTIONS:
        if heading.casefold() not in headings:
            errors.append(f"{name}: missing section '{heading}'")
    for aggregation in ("dynamic", "group"):
        if not re.search(r"\.aggte\(\s*typec\s*=\s*['\"]" + aggregation + r"['\"]", content):
            errors.append(f"{name}: missing {aggregation} aggregation example")
    images = image_references(content)
    for badge, image_url in REQUIRED_BADGES.items():
        sources = [source for alt, source in images if alt == badge and source.startswith(image_url)]
        if not sources:
            errors.append(f"{name}: missing {badge} badge for csdid")
        if badge == "Downloads":
            for source in sources:
                if parse_qs(urlsplit(source).query).get("left_text") != ["Downloads"]:
                    errors.append(f"{name}: Downloads badge left_text must be exactly Downloads")
    if not re.search(r"@software\s*\{csdid\s*,", content):
        errors.append(f"{name}: missing csdid software citation")
    return errors


def check_readme_paths(paths):
    """Git can track case aliases that a macOS checkout cannot distinguish."""
    aliases = sorted(path for path in paths if path.casefold() == "readme.md")
    if aliases != ["README.md"]:
        return [f"keep one canonical README.md; tracked README paths: {aliases}"]
    return []


def check_images(content, root):
    errors = []
    references = [source for _, source in image_references(content) if source.startswith("README_files/")]
    if not references:
        errors.append("README.md: no generated example figures")
    for reference in references:
        path = (root / reference).resolve()
        try:
            path.relative_to(root.resolve())
        except ValueError:
            errors.append(f"README.md: figure path escapes the repository: {reference}")
            continue
        if not path.is_file():
            errors.append(f"README.md: missing generated figure {reference}")
    return errors


def check_repository(root, source_only=False):
    errors = []
    tracked = subprocess.run(
        ["git", "-C", str(root), "ls-files", "-z"],
        check=True, capture_output=True, text=True,
    ).stdout.split("\0")
    errors.extend(check_readme_paths(tracked))
    names = ["README.qmd"] if source_only else ["README.qmd", "README.md"]
    for name in names:
        path = root / name
        if not path.is_file():
            errors.append(f"missing {name}")
            continue
        content = path.read_text(encoding="utf-8")
        errors.extend(check_document(content, name))
        if name == "README.md":
            errors.extend(check_images(content, root))
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--source-only", action="store_true")
    args = parser.parse_args()
    errors = check_repository(args.root, args.source_only)
    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        return 1
    print("README content, canonical path, and generated figures are valid.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
