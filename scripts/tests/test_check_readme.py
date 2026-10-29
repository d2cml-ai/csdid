"""Exercise the content losses observed in the generated README bot PR."""

import importlib.util
from pathlib import Path
import re
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("check_readme", ROOT / "scripts/check_readme.py")
guard = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(guard)


class ReadmeGuardTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = (ROOT / "README.qmd").read_text(encoding="utf-8")
        cls.rendered = (ROOT / "README.md").read_text(encoding="utf-8")

    def test_repaired_source_and_rendered_readme_are_accepted(self):
        self.assertEqual(guard.check_document(self.source, "README.qmd"), [])
        self.assertEqual(guard.check_document(self.rendered, "README.md"), [])
        self.assertEqual(guard.check_images(self.rendered, ROOT), [])

    def test_old_bot_package_rename_is_rejected(self):
        damaged = self.rendered.replace("**csdid** package implements", "**did** package implements")
        self.assertTrue(damaged != self.rendered)
        self.assertTrue(any("renamed to did" in error for error in guard.check_document(damaged, "bot")))

    def test_each_lost_badge_is_rejected(self):
        for badge in guard.REQUIRED_BADGES:
            with self.subTest(badge=badge):
                alt = r"\s+".join(re.escape(word) for word in badge.split())
                damaged = re.sub(r"\[!\[" + alt + r"\]\([^)]+\)\]\([^)]+\)\s*", "", self.rendered)
                damaged = re.sub(r"<img\b(?=[^>]*\balt=['\"]" + alt + r"['\"])[^>]*>", "", damaged,
                                 flags=re.IGNORECASE)
                self.assertTrue(damaged != self.rendered)
                self.assertTrue(any(badge + " badge" in error for error in guard.check_document(damaged, "bot")))

    def test_html_download_badge_accepts_escaped_query(self):
        for content in (self.source, self.rendered):
            with self.subTest(document="source" if content == self.source else "rendered"):
                images = dict(guard.image_references(content))
                original = images["Downloads"]
                escaped = original.replace("&", "&amp;")
                html_badge = '<img src="' + escaped + '" alt="Downloads" />'
                candidate = re.sub(r"<img\b(?=[^>]*\balt=['\"]Downloads['\"])[^>]*>", html_badge, content,
                                   flags=re.IGNORECASE)
                candidate = re.sub(r"!\[Downloads\]\([^)]+\)", html_badge, candidate)
                self.assertEqual(guard.check_document(candidate, "HTML badge"), [])

    def test_real_quarto_download_label_corruption_is_rejected(self):
        # Quarto appended .png to this exact extensionless SVG URL during regeneration.
        corrupt_url = (
            "https://static.pepy.tech/personalized-badge/csdid?period=total"
            "&units=international_system&left_color=blue&right_color=grey"
            "&left_text=Downloads.png"
        )
        for image in ('![Downloads](' + corrupt_url + ')',
                      '<img alt="Downloads" src="' + corrupt_url.replace("&", "&amp;") + '">'):
            with self.subTest(markup=image[:15]):
                errors = guard.check_document(self.rendered + "\n" + image, "corrupted badge")
                self.assertTrue(any("left_text must be exactly Downloads" in error for error in errors))

    def test_lost_examples_and_citation_are_rejected(self):
        for title in ("Event studies", "Overall effect of participating in treatment", "How to cite"):
            with self.subTest(section=title):
                damaged = re.sub(r"^## " + re.escape(title) + r"\n.*?(?=^## |\Z)", "", self.rendered,
                                 flags=re.MULTILINE | re.DOTALL | re.IGNORECASE)
                self.assertTrue(damaged != self.rendered)
                self.assertTrue(any("missing section '" + title + "'" in error
                                    for error in guard.check_document(damaged, "bot")))

    def test_citation_text_is_required_even_when_heading_remains(self):
        damaged = self.rendered.replace("@software{csdid,", "@software{other,")
        self.assertTrue(damaged != self.rendered)
        self.assertTrue(any("software citation" in error for error in guard.check_document(damaged, "bot")))

    def test_example_code_is_required_even_when_headings_remain(self):
        for aggregation in ("dynamic", "group"):
            with self.subTest(aggregation=aggregation):
                damaged = re.sub(r"\.aggte\(\s*typec\s*=\s*['\"]" + aggregation + r"['\"]", ".removed()", self.rendered)
                self.assertTrue(damaged != self.rendered)
                self.assertTrue(any(aggregation + " aggregation example" in error
                                    for error in guard.check_document(damaged, "bot")))

    def test_case_variant_readmes_are_rejected(self):
        self.assertEqual(guard.check_readme_paths(["README.md", "README.qmd"]), [])
        self.assertTrue(guard.check_readme_paths(["README.md", "readme.md"]))
        self.assertTrue(guard.check_readme_paths(["readme.md"]))

    def test_missing_generated_image_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            errors = guard.check_images(self.rendered, Path(directory))
        self.assertTrue(any("missing generated figure" in error for error in errors))


if __name__ == "__main__":
    unittest.main()
