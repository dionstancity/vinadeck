from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("STREAMLIT_GLOBAL_DEVELOPMENT_MODE", "false")
os.environ.setdefault("STREAMLIT_BROWSER_GATHER_USAGE_STATS", "false")
os.environ.setdefault("STREAMLIT_SERVER_FILE_WATCHER_TYPE", "none")
os.environ.setdefault("STREAMLIT_SERVER_HEADLESS", "false")
os.environ.setdefault("STREAMLIT_LOGGER_LEVEL", "info")

from streamlit import config as st_config
from streamlit.web import bootstrap

st_config.set_option("global.developmentMode", False)


def resolve_app_path() -> Path:
    candidates: list[Path] = []
    if hasattr(sys, "_MEIPASS"):
        meipass = Path(sys._MEIPASS)
        candidates.extend([meipass / "app.py", meipass / "_internal" / "app.py"])

    executable_dir = Path(sys.executable).resolve().parent
    candidates.extend([executable_dir / "app.py", executable_dir / "_internal" / "app.py"])
    candidates.append(Path(__file__).resolve().with_name("app.py"))

    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError("Could not locate bundled app.py")


def main() -> None:
    bootstrap.run(
        str(resolve_app_path()),
        False,
        [],
        {
            "browser.gatherUsageStats": False,
            "logger.level": "info",
            "server.fileWatcherType": "none",
            "server.headless": False,
            "global.developmentMode": False,
        },
    )


if __name__ == "__main__":
    main()
