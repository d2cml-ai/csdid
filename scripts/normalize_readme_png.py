"""Restore HEAD figure bytes only for PNG compression/chunking-only changes."""
import argparse
import os
from pathlib import Path
import struct
import subprocess
import sys
import zlib


FIGURES = Path("README_files/figure-commonmark")
PNG_SIGNATURE = b"\x89PNG\r\n\x1a\n"
MAX_FILE_BYTES = 32 * 1024 * 1024
MAX_INFLATED_BYTES = 64 * 1024 * 1024


class PNGError(ValueError):
    pass


def fingerprint(data, max_inflated=MAX_INFLATED_BYTES):
    """Compare containers and scanline streams, without decoding pixels."""
    if max_inflated < 1:
        raise PNGError("PNG inflated-data limit must be positive")
    if len(data) > MAX_FILE_BYTES or not data.startswith(PNG_SIGNATURE):
        raise PNGError("Invalid or oversized PNG signature/file")
    offset, chunks, compressed = len(PNG_SIGNATURE), [], []
    idat_started = idat_ended = ended = False
    while offset < len(data):
        start = offset
        if len(data) - offset < 12:
            raise PNGError("Truncated PNG chunk")
        length = struct.unpack_from(">I", data, offset)[0]
        kind = data[offset + 4:offset + 8]
        if (len(kind) != 4 or any(not (65 <= value <= 90 or 97 <= value <= 122) for value in kind)
                or not 65 <= kind[2] <= 90):
            raise PNGError("Invalid PNG chunk type")
        offset += 8
        if length > len(data) - offset - 4:
            raise PNGError("Truncated PNG chunk payload")
        payload = data[offset:offset + length]
        crc = struct.unpack_from(">I", data, offset + length)[0]
        if zlib.crc32(kind + payload) & 0xffffffff != crc:
            raise PNGError("Invalid PNG chunk CRC")
        offset += length + 4
        if not chunks:
            if kind != b"IHDR" or length != 13:
                raise PNGError("PNG must start with one IHDR")
            width, height, depth, color, compression, filtering, interlace = struct.unpack(">IIBBBBB", payload)
            depths = {0: (1, 2, 4, 8, 16), 2: (8, 16), 3: (1, 2, 4, 8), 4: (8, 16), 6: (8, 16)}
            if (not 0 < width < 2 ** 31 or not 0 < height < 2 ** 31
                    or depth not in depths.get(color, ()) or compression or filtering or interlace not in (0, 1)):
                raise PNGError("Invalid PNG IHDR")
        elif kind == b"IHDR":
            raise PNGError("Duplicate PNG IHDR")
        if kind == b"IDAT":
            if idat_ended:
                raise PNGError("PNG IDAT chunks must be consecutive")
            if not idat_started:
                chunks.append(None)  # Preserve metadata position around the stream.
            idat_started = True
            compressed.append(payload)
        else:
            if idat_started:
                idat_ended = True
            chunks.append(data[start:offset])
        if kind == b"IEND":
            if length or not idat_started or offset != len(data):
                raise PNGError("Invalid PNG IEND or trailing bytes")
            ended = True
            break
    if not ended:
        raise PNGError("PNG is missing IEND")
    inflater = zlib.decompressobj()
    try:
        inflated = inflater.decompress(b"".join(compressed), max_inflated + 1)
    except zlib.error as error:
        raise PNGError("Invalid PNG zlib stream") from error
    if len(inflated) > max_inflated or inflater.unconsumed_tail:
        raise PNGError("PNG inflated data exceeds the limit")
    if not inflater.eof or inflater.unused_data:
        raise PNGError("Incomplete PNG zlib stream or trailing compressed data")
    if not inflated:
        raise PNGError("PNG scanline stream is empty")
    return tuple(chunks), inflated


def normalize(repo):
    """Validate every generated figure before restoring any proven equivalent."""
    repo = Path(repo).resolve()
    def git(*args):
        return subprocess.check_output(["git", "-C", str(repo), *args])
    head = git("rev-parse", "HEAD").decode().strip()
    tracked = {}
    for record in git("ls-tree", "-r", "-z", head, "--", FIGURES.as_posix()).split(b"\0"):
        if record:
            metadata, name = record.split(b"\t", 1)
            mode, kind, oid = metadata.decode().split()
            path = Path(os.fsdecode(name))
            if path.suffix == ".png":
                tracked[path] = (mode, kind, oid)
    folder = repo / FIGURES
    for parent in (repo / "README_files", folder):
        if parent.is_symlink():
            raise PNGError("Generated figure directories must not be symlinks")
    generated = []
    for directory, folders, files in os.walk(folder, followlinks=False):
        if any((Path(directory) / name).is_symlink() for name in folders):
            raise PNGError("Generated figure directories must not be symlinks")
        generated.extend(Path(directory) / name for name in files if name.endswith(".png"))
    restore = []
    for path in sorted(generated):
        relative = path.relative_to(repo)
        if path.is_symlink() or not path.is_file():
            raise PNGError("Generated PNG must be a regular file: " + str(relative))
        with path.open("rb") as stream:
            current = stream.read(MAX_FILE_BYTES + 1)
        try:
            current_fingerprint = fingerprint(current)
            if relative not in tracked:
                continue  # Validate new figures, but never replace them.
            mode, kind, oid = tracked[relative]
            if kind != "blob" or mode not in ("100644", "100755"):
                raise PNGError("Tracked PNG must be a regular file")
            if int(git("cat-file", "-s", oid)) > MAX_FILE_BYTES:
                raise PNGError("Tracked PNG exceeds the file limit")
            original = git("cat-file", "blob", oid)
            if current != original and current_fingerprint == fingerprint(original):
                restore.append((path, original))
        except PNGError as error:
            raise PNGError(str(relative) + ": " + str(error)) from error
    # A malformed later figure cannot cause earlier candidates to be modified.
    for path, original in restore:
        path.write_bytes(original)
    return [path.relative_to(repo).as_posix() for path, _ in restore]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", default=".")
    args = parser.parse_args()
    try:
        restored = normalize(args.repo)
    except (PNGError, OSError, subprocess.CalledProcessError) as error:
        print("PNG normalization blocked: " + str(error), file=sys.stderr)
        return 1
    print("Restored %d compression-only PNG changes." % len(restored))
    for path in restored:
        print(path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
