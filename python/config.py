"""
Central configuration loader (mirrors the generatoroutputanalysis convention).

All bulk-data and report paths live in a JSON file.  Resolution order:
    1. $FASERCAL_SETTINGS               (explicit override)
    2. json/ev.json                    (personal, git-ignored)
    3. json/settings_template.json     (committed fallback template)

Import ``get(key)`` to resolve a path.  Keys beginning with ``_`` are comments.
Shell scripts share the resolution via:  python3 python/config.py <key>
"""
import json
import os
import sys
from pathlib import Path


def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def settings_path() -> Path:
    env = os.environ.get("FASERCAL_SETTINGS")
    if env:
        return Path(env)
    ev = repo_root() / "json" / "ev.json"
    if ev.exists():
        return ev
    return repo_root() / "json" / "settings_template.json"


def load_settings() -> dict:
    with open(settings_path()) as f:
        return {k: v for k, v in json.load(f).items() if not k.startswith("_")}


SETTINGS = load_settings()


def get(key: str, default=None):
    """Return a configured path (str) by key, or ``default`` if absent."""
    return SETTINGS.get(key, default)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: config.py <settings-key>")
    val = get(sys.argv[1])
    if val is None:
        sys.exit(f"config.py: no key '{sys.argv[1]}' in {settings_path()}")
    print(val)
