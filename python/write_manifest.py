import hashlib
import subprocess
from pathlib import Path


def main():
    root = Path(__file__).resolve().parents[1]
    listed = subprocess.run(["git", "ls-files"], cwd=root, check=True,
                            capture_output=True, text=True).stdout.splitlines()
    lines = []
    for name in sorted(listed):
        if name == "MANIFEST.sha256":
            continue
        digest = hashlib.sha256((root / name).read_bytes()).hexdigest()
        lines.append(f"{digest}  {name}")
    (root / "MANIFEST.sha256").write_text("\n".join(lines) + "\n")
    print(f"MANIFEST.sha256 lists {len(lines)} tracked files")


if __name__ == "__main__":
    main()
