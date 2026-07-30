import fcntl
import hashlib
import json
import os
import ssl
import urllib.error
import urllib.request
from typing import Any


BASE_DIR = os.path.abspath(os.path.dirname(__file__))
DATA_DIR = os.environ.get("DEPMAP_DATA_DIR", os.path.join(BASE_DIR, "data"))
USER_AGENT = "DepMap-Dependency-Explorer/1.0"
DOWNLOAD_CHUNK_BYTES = 4 * 1024 * 1024
CRISPR_RUNTIME_RELEASE = "DepMap Public 24Q4"
RNAI_RUNTIME_RELEASE = "DEMETER2 Data v6"

RUNTIME_FILES = (
    {
        "release": CRISPR_RUNTIME_RELEASE,
        "filename": "Model.csv",
        "minimum_bytes": 100_000,
        "url": "https://ndownloader.figshare.com/files/51065297",
        "md5_hash": "675210d17675f3517b0ce39a3c274f16",
    },
    {
        "release": CRISPR_RUNTIME_RELEASE,
        "filename": "CRISPRGeneEffect.csv",
        "minimum_bytes": 10_000_000,
        "url": "https://ndownloader.figshare.com/files/51064667",
        "md5_hash": "6edf7ade09b9b34199210b559d4745d3",
    },
    {
        "release": RNAI_RUNTIME_RELEASE,
        "filename": "D2_combined_gene_dep_scores.csv",
        "minimum_bytes": 10_000_000,
        "url": "https://ndownloader.figshare.com/files/13515395",
        "md5_hash": "e82566a48993753c20ea5f480f1867ec",
    },
)


def _manifest_path() -> str:
    return os.path.join(DATA_DIR, ".runtime-data.json")


def _load_manifest() -> dict[str, Any]:
    try:
        with open(_manifest_path()) as manifest_file:
            payload = json.load(manifest_file)
        return payload if isinstance(payload, dict) else {}
    except (FileNotFoundError, json.JSONDecodeError, OSError):
        return {}


def _mark_managed_file(spec: dict[str, Any]) -> None:
    path = _manifest_path()
    partial_path = f"{path}.part-{os.getpid()}"
    payload = _load_manifest()
    payload.setdefault("files", {})[spec["filename"]] = {
        "release": spec["release"],
        "md5_hash": spec["md5_hash"],
    }
    try:
        with open(partial_path, "w") as manifest_file:
            json.dump(payload, manifest_file, separators=(",", ":"))
        os.replace(partial_path, path)
    finally:
        if os.path.exists(partial_path):
            os.remove(partial_path)


def runtime_data_status() -> dict[str, Any]:
    manifest = _load_manifest().get("files", {})
    files = {}
    for spec in RUNTIME_FILES:
        path = os.path.join(DATA_DIR, spec["filename"])
        size = os.path.getsize(path) if os.path.isfile(path) else 0
        managed_file = manifest.get(spec["filename"], {})
        files[spec["filename"]] = {
            "ready": size >= spec["minimum_bytes"],
            "bytes": size,
            "release": managed_file.get("release") or "Existing server data",
        }
    return {
        "ready": all(item["ready"] for item in files.values()),
        "directory": DATA_DIR,
        "files": files,
    }


def _open_url(url: str):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        return urllib.request.urlopen(request, timeout=120)
    except urllib.error.URLError as exc:
        certificate_error = isinstance(exc.reason, ssl.SSLCertVerificationError)
        if (
            not certificate_error
            or os.environ.get("DEPMAP_ALLOW_UNVERIFIED_TLS", "false").lower() != "true"
        ):
            raise
        return urllib.request.urlopen(
            request,
            context=ssl._create_unverified_context(),
            timeout=120,
        )


def _download_file(record: dict[str, Any], destination: str, minimum_bytes: int) -> None:
    partial_path = f"{destination}.part-{os.getpid()}"
    expected_md5 = (record.get("md5_hash") or "").strip().lower()
    digest = hashlib.md5(usedforsecurity=False)
    downloaded = 0
    try:
        with _open_url(record["url"]) as response, open(partial_path, "wb") as output:
            while True:
                chunk = response.read(DOWNLOAD_CHUNK_BYTES)
                if not chunk:
                    break
                output.write(chunk)
                digest.update(chunk)
                downloaded += len(chunk)
        if downloaded < minimum_bytes:
            raise RuntimeError(
                f"Downloaded {record['filename']} is unexpectedly small ({downloaded:,} bytes)."
            )
        if expected_md5 and digest.hexdigest().lower() != expected_md5:
            raise RuntimeError(f"Checksum verification failed for {record['filename']}.")
        with open(partial_path, "rb") as downloaded_file:
            prefix = downloaded_file.read(64).lstrip().lower()
        if prefix.startswith((b"<!doctype html", b"<html")):
            raise RuntimeError(f"Download for {record['filename']} returned HTML instead of CSV data.")
        os.replace(partial_path, destination)
        print(f"Provisioned {record['filename']} ({downloaded / 1048576:.1f} MB)", flush=True)
    finally:
        if os.path.exists(partial_path):
            os.remove(partial_path)


def ensure_depmap_runtime_data() -> dict[str, Any]:
    status = runtime_data_status()
    if status["ready"]:
        return status

    os.makedirs(DATA_DIR, exist_ok=True)
    lock_path = os.path.join(DATA_DIR, ".runtime-data.lock")
    with open(lock_path, "w") as lock_file:
        fcntl.flock(lock_file, fcntl.LOCK_EX)
        status = runtime_data_status()
        if status["ready"]:
            return status

        print(f"Provisioning raw DepMap data in {DATA_DIR}", flush=True)
        for spec in RUNTIME_FILES:
            destination = os.path.join(DATA_DIR, spec["filename"])
            if os.path.isfile(destination) and os.path.getsize(destination) >= spec["minimum_bytes"]:
                continue
            _download_file(spec, destination, spec["minimum_bytes"])
            _mark_managed_file(spec)

    status = runtime_data_status()
    if not status["ready"]:
        raise RuntimeError("Raw DepMap runtime data provisioning did not complete.")
    return status
