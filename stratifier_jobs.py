"""Persisted jobs for the single-host Render service, shared by all web workers."""

import fcntl
import logging
import os
import threading
import uuid
from datetime import datetime, timezone
from contextlib import closing


def timestamp():
    return datetime.now(timezone.utc).isoformat()


class StratifierJobs:
    def __init__(self, connect, execute, compute, save, directory):
        self.connect = lambda: closing(connect())
        self.execute = execute
        self.compute = compute
        self.save = save
        self.directory = directory
        self.thread = None
        self.guard = threading.Lock()
        with self.connect() as db:
            self.execute(db, """
                CREATE TABLE IF NOT EXISTS stratifier_jobs (
                    id TEXT PRIMARY KEY, request_key TEXT NOT NULL UNIQUE,
                    prompt TEXT NOT NULL, status TEXT NOT NULL, stage TEXT NOT NULL,
                    error TEXT NOT NULL DEFAULT '', result_id INTEGER,
                    attempts INTEGER NOT NULL DEFAULT 0,
                    created_at TEXT NOT NULL, updated_at TEXT NOT NULL
                )
            """)
            db.commit()

    def enqueue(self, prompt, request_key):
        job_id = uuid.uuid4().hex
        now = timestamp()
        with self.connect() as db:
            self.execute(db, """
                INSERT INTO stratifier_jobs
                    (id, request_key, prompt, status, stage, created_at, updated_at)
                VALUES (?, ?, ?, 'queued', 'Waiting to start', ?, ?)
                ON CONFLICT(request_key) DO NOTHING
            """, (job_id, request_key, prompt, now, now))
            db.commit()
            job = dict(self.execute(db, "SELECT * FROM stratifier_jobs WHERE request_key=?", (request_key,)).fetchone())
        self.start()
        return job

    def get(self, job_id):
        with self.connect() as db:
            row = self.execute(db, "SELECT * FROM stratifier_jobs WHERE id=?", (job_id,)).fetchone()
            return dict(row) if row else None

    def recent(self):
        with self.connect() as db:
            return [dict(row) for row in self.execute(db,
                "SELECT * FROM stratifier_jobs ORDER BY created_at DESC LIMIT 10").fetchall()]

    def progress(self, job_id, stage):
        with self.connect() as db:
            self.execute(db, "UPDATE stratifier_jobs SET stage=?, updated_at=? WHERE id=?",
                         (stage, timestamp(), job_id))
            db.commit()

    def start(self):
        with self.guard:
            if self.thread and self.thread.is_alive():
                return
            self.thread = threading.Thread(target=self.run, daemon=True, name="stratifier-builder")
            self.thread.start()

    def run(self):
        os.makedirs(self.directory, exist_ok=True)
        # flock also serializes workers in other processes and releases on a crash.
        with open(os.path.join(self.directory, ".stratifier-jobs.lock"), "a") as lock:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                return
            with self.connect() as db:
                self.execute(db, """UPDATE stratifier_jobs SET status='queued',
                    stage='Resuming after server restart' WHERE status='running' AND attempts < 3""")
                self.execute(db, """UPDATE stratifier_jobs SET status='failed',
                    error='The server interrupted this build repeatedly. Please retry.',
                    stage='Build interrupted' WHERE status='running'""")
                db.commit()
            while True:
                with self.connect() as db:
                    row = self.execute(db, """SELECT * FROM stratifier_jobs
                        WHERE status='queued' ORDER BY created_at LIMIT 1""").fetchone()
                    if not row:
                        return
                    job = dict(row)
                    self.execute(db, """UPDATE stratifier_jobs SET status='running',
                        attempts=attempts+1, updated_at=? WHERE id=?""", (timestamp(), job["id"]))
                    db.commit()
                try:
                    analysis, source, quality = self.compute(job["prompt"],
                        progress=lambda stage: self.progress(job["id"], stage))
                    # Save the analysis and completion together: retries cannot duplicate a result.
                    with self.connect() as db:
                        result_id = self.save(db, job["prompt"], analysis, source, quality)
                        self.execute(db, """UPDATE stratifier_jobs SET status='complete',
                            stage='Ready to explore', result_id=?, updated_at=? WHERE id=?""",
                            (result_id, timestamp(), job["id"]))
                        db.commit()
                except Exception as exc:
                    logging.getLogger(__name__).exception("Stratifier job %s failed", job["id"])
                    with self.connect() as db:
                        self.execute(db, """UPDATE stratifier_jobs SET status='failed',
                            stage='Could not complete this comparison', error=?, updated_at=? WHERE id=?""",
                            (str(exc)[:1500], timestamp(), job["id"]))
                        db.commit()
