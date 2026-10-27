"""Small file primitives shared at private cache and input boundaries."""

import os
import stat
import tempfile
from contextlib import contextmanager, suppress
from pathlib import Path


def stat_fingerprint(status):
    return status.st_dev, status.st_ino, status.st_size, status.st_mtime_ns, status.st_ctime_ns


def private_directory(path: Path, *, create=True):
    if not isinstance(path, Path):
        raise TypeError("cache directory must be a Path")
    try:
        status = os.lstat(path)
    except FileNotFoundError:
        if not create:
            raise
        try:
            os.mkdir(path, 0o700)
        except FileExistsError:
            pass
        status = os.lstat(path)
    if not stat.S_ISDIR(status.st_mode):
        raise OSError("private cache directory required")
    if os.name == "posix" and (status.st_uid != os.getuid() or status.st_mode & 0o077):
        raise OSError("owner-private cache directory required")
    return path


@contextmanager
def open_regular(path: Path, *, follow_symlinks=True):
    flags = os.O_RDONLY | getattr(os, "O_NONBLOCK", 0) | getattr(os, "O_BINARY", 0)
    if not follow_symlinks:
        flags |= getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags)
    try:
        status = os.fstat(descriptor)
        if not stat.S_ISREG(status.st_mode):
            raise OSError("regular file required")
        stream = os.fdopen(descriptor, "rb")
        descriptor = -1
        with stream:
            yield stream, status
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def read_regular_bounded(path: Path, maximum: int, *, follow_symlinks=True, check_size=False):
    """Read one regular file with explicit symlink and stat-size policies."""
    if type(maximum) is not int or maximum < 0:
        raise ValueError("maximum must be a nonnegative integer")
    with open_regular(path, follow_symlinks=follow_symlinks) as (stream, before):
        if check_size and before.st_size > maximum:
            raise OverflowError("regular file exceeds maximum")
        content = stream.read(maximum + 1)
        after = os.fstat(stream.fileno())
    if len(content) > maximum:
        raise OverflowError("regular file exceeds maximum")
    return content, before, after


def atomic_replace_bytes(path: Path, content: bytes, *, prefix=None, sync=True, mode=None):
    """Replace a file atomically; the caller owns directory and error policy."""
    if not isinstance(path, Path) or not isinstance(content, bytes):
        raise TypeError("Path and bytes required")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=f".{path.name}." if prefix is None else prefix,
            suffix=".tmp",
            dir=path.parent,
            delete=False,
        ) as stream:
            temporary = Path(stream.name)
            stream.write(content)
            if sync:
                stream.flush()
                os.fsync(stream.fileno())
        os.replace(temporary, path)
        temporary = None
        if mode is not None:
            os.chmod(path, mode)
    finally:
        if temporary is not None:
            with suppress(OSError):
                temporary.unlink(missing_ok=True)
