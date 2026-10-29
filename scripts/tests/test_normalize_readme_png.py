"""PNG compression-only counterexamples; standard library and Git only."""
import importlib.util
from pathlib import Path
import struct
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import zlib


SPEC = importlib.util.spec_from_file_location(
    "normalize_readme_png", Path(__file__).parents[1] / "normalize_readme_png.py")
png = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(png)
RAW = (b"\x00" + b"\x10\x20\x30\xff" * 2) * 2


def chunk(kind, data):
    return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data) & 0xffffffff)


def image(raw=RAW, level=9, text=b"Software\x00test", split=False, compressed=None, text_after=False):
    header = chunk(b"IHDR", struct.pack(">IIBBBBB", 2, 2, 8, 6, 0, 0, 0))
    value = zlib.compress(raw, level) if compressed is None else compressed
    parts = [value[:len(value) // 2], value[len(value) // 2:]] if split else [value]
    idat = b"".join(chunk(b"IDAT", part) for part in parts)
    metadata = chunk(b"tEXt", text)
    return png.PNG_SIGNATURE + header + (idat + metadata if text_after else metadata + idat) + chunk(b"IEND", b"")


class PNGFingerprints(unittest.TestCase):
    def test_compression_levels_and_idat_segmentation_match_exact_scanline_stream(self):
        original, generated = image(), image(level=0, split=True)
        self.assertNotEqual(original, generated)
        self.assertEqual(png.fingerprint(original), png.fingerprint(generated))

    def test_changed_scanlines_metadata_or_metadata_position_remain_distinct(self):
        for generated in (image(raw=RAW[:-1] + b"\x00"), image(text=b"Software\x00changed"), image(text_after=True)):
            with self.subTest(generated=generated):
                self.assertNotEqual(png.fingerprint(image()), png.fingerprint(generated))

    def test_signature_truncation_crc_headers_and_trailing_bytes_fail(self):
        valid = image()
        bad_crc = valid[:-1] + bytes([valid[-1] ^ 1])
        no_idat = png.PNG_SIGNATURE + chunk(b"IHDR", struct.pack(">IIBBBBB", 2, 2, 8, 6, 0, 0, 0)) + chunk(b"IEND", b"")
        for invalid in (b"not a PNG", valid[:10], valid[:-1], bad_crc, valid + b"trailing",
                        png.PNG_SIGNATURE + chunk(b"IDAT", zlib.compress(RAW)) + chunk(b"IEND", b""),
                        no_idat, png.PNG_SIGNATURE + chunk(b"IHDR", b"short") + valid[33:]):
            with self.subTest(invalid=invalid), self.assertRaises(png.PNGError):
                png.fingerprint(invalid)

    def test_nonconsecutive_idat_and_duplicate_header_fail(self):
        header = chunk(b"IHDR", struct.pack(">IIBBBBB", 2, 2, 8, 6, 0, 0, 0))
        value = zlib.compress(RAW)
        invalid = (png.PNG_SIGNATURE + header + chunk(b"IDAT", value[:3]) + chunk(b"tEXt", b"a\x00b")
                   + chunk(b"IDAT", value[3:]) + chunk(b"IEND", b""))
        for candidate in (invalid, png.PNG_SIGNATURE + header + image()[8:]):
            with self.subTest(candidate=candidate), self.assertRaises(png.PNGError):
                png.fingerprint(candidate)

    def test_invalid_truncated_or_concatenated_zlib_streams_fail(self):
        value = zlib.compress(RAW)
        for compressed in (b"invalid", value[:-1], value + b"extra", value + value, zlib.compress(b"")):
            with self.subTest(compressed=compressed), self.assertRaises(png.PNGError):
                png.fingerprint(image(compressed=compressed))

    def test_file_and_inflated_limits_are_enforced_without_unbounded_flush(self):
        with self.assertRaisesRegex(png.PNGError, "exceeds"):
            png.fingerprint(image(raw=b"\x00" * 100000), max_inflated=32)
        with self.assertRaisesRegex(png.PNGError, "positive"):
            png.fingerprint(image(), max_inflated=-1)
        with patch.object(png, "MAX_FILE_BYTES", 10), self.assertRaisesRegex(png.PNGError, "oversized"):
            png.fingerprint(image())
        self.assertEqual(png.fingerprint(image(), max_inflated=len(RAW)), png.fingerprint(image()))


class PNGRestoration(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.folder = self.root / png.FIGURES
        self.folder.mkdir(parents=True)
        self.original = self.folder / "figure.png"
        self.original.write_bytes(image())
        self.git("init", "-q")
        self.git("config", "user.email", "png-test@example.invalid")
        self.git("config", "user.name", "PNG test")
        self.git("config", "commit.gpgsign", "false")
        self.git("add", ".")
        self.git("commit", "-qm", "tracked figure")

    def git(self, *args):
        return subprocess.check_output(["git", "-C", str(self.root), *args], stderr=subprocess.DEVNULL)

    def test_compression_only_changes_restore_exact_head_bytes_and_are_idempotent(self):
        self.original.write_bytes(image(level=0, split=True))
        self.assertEqual(png.normalize(self.root), [(png.FIGURES / "figure.png").as_posix()])
        self.assertEqual(self.original.read_bytes(), image())
        self.assertEqual(self.git("diff", "--name-only"), b"")
        self.assertEqual(png.normalize(self.root), [])

    def test_changed_pixels_metadata_new_figures_and_removals_are_preserved(self):
        for changed in (image(raw=RAW[:-1] + b"\x00"), image(text=b"Software\x00changed")):
            self.original.write_bytes(changed)
            self.assertEqual(png.normalize(self.root), [])
            self.assertEqual(self.original.read_bytes(), changed)
        added = self.folder / "new.png"
        added.write_bytes(image(level=0))
        self.original.unlink()
        self.assertEqual(png.normalize(self.root), [])
        self.assertFalse(self.original.exists())
        self.assertEqual(added.read_bytes(), image(level=0))

    def test_malformed_new_figure_blocks_before_any_other_file_is_restored(self):
        changed = image(level=0)
        self.original.write_bytes(changed)
        (self.folder / "later.png").write_bytes(image()[:-1])
        with self.assertRaises(png.PNGError):
            png.normalize(self.root)
        self.assertEqual(self.original.read_bytes(), changed)

    def test_bad_crc_or_invalid_zlib_in_new_figures_is_not_ignored(self):
        for invalid in (image()[:-1] + b"\xff", image(compressed=b"invalid")):
            added = self.folder / "new.png"
            added.write_bytes(invalid)
            with self.subTest(invalid=invalid), self.assertRaises(png.PNGError):
                png.normalize(self.root)

    def test_symlink_figure_cannot_write_outside_the_repository(self):
        outside = self.root / "outside.png"
        outside.write_bytes(image(level=0))
        self.original.unlink()
        self.original.symlink_to(outside)
        with self.assertRaises(png.PNGError):
            png.normalize(self.root)
        self.assertEqual(outside.read_bytes(), image(level=0))


if __name__ == "__main__":
    unittest.main()
