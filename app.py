import csv
import json
import math
import os
import sqlite3
from datetime import datetime, timezone
from functools import lru_cache, wraps
from typing import Any

from flask import (
    Flask,
    flash,
    g,
    redirect,
    render_template,
    request,
    session,
    url_for,
)
from werkzeug.security import check_password_hash, generate_password_hash

from dependency_stratifiers import MODEL_PATH, compute_stratifier


BASE_DIR = os.path.abspath(os.path.dirname(__file__))
DB_PATH = os.environ.get("DB_PATH", os.path.join(BASE_DIR, "app.db"))
DATABASE_URL = os.environ.get("DATABASE_URL", "").strip()
USING_POSTGRES = bool(DATABASE_URL)
DEPENDENCY_SUMMARY_PATH = os.path.join(BASE_DIR, "static", "data", "hpv_dependency_summary.json")
UTC = timezone.utc
MIN_DEPENDENCY_GROUP_VALUES = 3
MIN_DEPENDENCY_GROUP_COVERAGE = 0.5


app = Flask(__name__)
app.config["SECRET_KEY"] = os.environ.get("SECRET_KEY", "dev-secret-change-me")
app.config["ALLOW_SELF_SIGNUP"] = os.environ.get("ALLOW_SELF_SIGNUP", "true").lower() == "true"


def now_utc() -> datetime:
    return datetime.now(tz=UTC)


def to_iso(dt: datetime) -> str:
    return dt.astimezone(UTC).replace(microsecond=0).isoformat()


def from_iso(value: str) -> datetime:
    return datetime.fromisoformat(value).astimezone(UTC)


def connect_db() -> Any:
    if USING_POSTGRES:
        try:
            import psycopg
            from psycopg.rows import dict_row
        except ImportError as exc:
            raise RuntimeError(
                "DATABASE_URL is set, but PostgreSQL support is not installed. "
                "Run `pip install -r requirements.txt`."
            ) from exc
        return psycopg.connect(DATABASE_URL, row_factory=dict_row)

    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def db_execute(db: Any, statement: str, params: tuple[Any, ...] = ()) -> Any:
    if USING_POSTGRES:
        statement = statement.replace("?", "%s")
    return db.execute(statement, params)


def insert_row(db: Any, statement: str, params: tuple[Any, ...]) -> int:
    if USING_POSTGRES:
        cursor = db_execute(db, f"{statement.rstrip().rstrip(';')} RETURNING id", params)
        return int(cursor.fetchone()["id"])
    cursor = db_execute(db, statement, params)
    return int(cursor.lastrowid)


def get_db() -> Any:
    if "db" not in g:
        g.db = connect_db()
    return g.db


@app.teardown_appcontext
def close_db(exception: Exception | None) -> None:
    db = g.pop("db", None)
    if db is not None:
        db.close()


def init_db() -> None:
    if not USING_POSTGRES:
        db_dir = os.path.dirname(DB_PATH)
        if db_dir:
            os.makedirs(db_dir, exist_ok=True)

    sqlite_schema = (
        """
        CREATE TABLE IF NOT EXISTS users (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            email TEXT NOT NULL UNIQUE,
            first_name TEXT NOT NULL,
            password_hash TEXT NOT NULL,
            role TEXT NOT NULL CHECK(role IN ('admin', 'forecaster', 'viewer')) DEFAULT 'forecaster',
            is_active INTEGER NOT NULL DEFAULT 1,
            created_at TEXT NOT NULL
        );

        CREATE TABLE IF NOT EXISTS audit_logs (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            actor_id INTEGER,
            event_type TEXT NOT NULL,
            entity_type TEXT NOT NULL,
            entity_id INTEGER,
            payload_json TEXT,
            created_at TEXT NOT NULL,
            FOREIGN KEY (actor_id) REFERENCES users(id)
        );

        CREATE TABLE IF NOT EXISTS dependency_stratifiers (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            label TEXT NOT NULL,
            prompt TEXT NOT NULL,
            analysis_json TEXT NOT NULL,
            source_json TEXT NOT NULL,
            quality_json TEXT NOT NULL,
            created_at TEXT NOT NULL,
            updated_at TEXT NOT NULL
        )
        """,
    )
    postgres_schema = (
        """
        CREATE TABLE IF NOT EXISTS users (
            id BIGSERIAL PRIMARY KEY,
            email TEXT NOT NULL UNIQUE,
            first_name TEXT NOT NULL,
            password_hash TEXT NOT NULL,
            role TEXT NOT NULL CHECK(role IN ('admin', 'forecaster', 'viewer')) DEFAULT 'forecaster',
            is_active INTEGER NOT NULL DEFAULT 1,
            created_at TEXT NOT NULL
        )
        """,
        """
        CREATE TABLE IF NOT EXISTS audit_logs (
            id BIGSERIAL PRIMARY KEY,
            actor_id BIGINT,
            event_type TEXT NOT NULL,
            entity_type TEXT NOT NULL,
            entity_id BIGINT,
            payload_json TEXT,
            created_at TEXT NOT NULL,
            FOREIGN KEY (actor_id) REFERENCES users(id)
        )
        """,
        """
        CREATE TABLE IF NOT EXISTS dependency_stratifiers (
            id BIGSERIAL PRIMARY KEY,
            label TEXT NOT NULL,
            prompt TEXT NOT NULL,
            analysis_json TEXT NOT NULL,
            source_json TEXT NOT NULL,
            quality_json TEXT NOT NULL,
            created_at TEXT NOT NULL,
            updated_at TEXT NOT NULL
        )
        """,
    )

    db = connect_db()
    for statement in postgres_schema if USING_POSTGRES else sqlite_schema:
        if USING_POSTGRES:
            db_execute(db, statement)
        else:
            db.executescript(statement)
    db.commit()

    admin_email = os.environ.get("ADMIN_EMAIL")
    admin_password = os.environ.get("ADMIN_PASSWORD")
    admin_name = os.environ.get("ADMIN_FIRST_NAME", "Admin")

    existing_admin = db_execute(
        db, "SELECT id FROM users WHERE role='admin' AND is_active=1 LIMIT 1"
    ).fetchone()
    if not existing_admin and admin_email and admin_password:
        db_execute(
            db,
            """
            INSERT INTO users (email, first_name, password_hash, role, is_active, created_at)
            VALUES (?, ?, ?, 'admin', 1, ?)
            """,
            (admin_email.lower().strip(), admin_name.strip(), generate_password_hash(admin_password), to_iso(now_utc())),
        )
        db.commit()

    db.close()


def log_event(actor_id: int | None, event_type: str, entity_type: str, entity_id: int | None, payload: dict[str, Any]) -> None:
    db = get_db()
    db_execute(
        db,
        """
        INSERT INTO audit_logs (actor_id, event_type, entity_type, entity_id, payload_json, created_at)
        VALUES (?, ?, ?, ?, ?, ?)
        """,
        (actor_id, event_type, entity_type, entity_id, json.dumps(payload), to_iso(now_utc())),
    )
    db.commit()


def get_user(user_id: int) -> Any | None:
    db = get_db()
    return db_execute(db, "SELECT * FROM users WHERE id=? AND is_active=1", (user_id,)).fetchone()


def current_user() -> Any | None:
    uid = session.get("user_id")
    if not uid:
        return None
    return get_user(int(uid))


def login_required(fn):
    @wraps(fn)
    def wrapper(*args, **kwargs):
        if not current_user():
            return redirect(url_for("login"))
        return fn(*args, **kwargs)

    return wrapper


def roles_required(*roles: str):
    def decorator(fn):
        @wraps(fn)
        def wrapper(*args, **kwargs):
            user = current_user()
            if not user or user["role"] not in roles:
                flash("You do not have permission to access this page.", "error")
                return redirect(url_for("dashboard"))
            return fn(*args, **kwargs)

        return wrapper

    return decorator


def load_builtin_dependency_summary() -> dict[str, Any]:
    with open(DEPENDENCY_SUMMARY_PATH) as f:
        return json.load(f)


def fetch_custom_stratifiers() -> list[Any]:
    return db_execute(
        get_db(),
        """
        SELECT *
        FROM dependency_stratifiers
        ORDER BY updated_at DESC
        """
    ).fetchall()


def custom_stratifier_analysis(row: Any) -> dict[str, Any]:
    analysis = json.loads(row["analysis_json"])
    analysis["id"] = f"custom-{row['id']}"
    analysis["label"] = f"{analysis['label']} *"
    analysis["category"] = "Custom stratifiers"
    analysis["custom_id"] = row["id"]
    filter_low_coverage_rows(analysis)
    return analysis


def filter_low_coverage_rows(analysis: dict[str, Any]) -> None:
    for dataset_key, rows in analysis.get("datasets", {}).items():
        included = analysis.get("included_models", {}).get(dataset_key, {})
        min_positive_n = int(
            included.get("min_positive_n")
            or min_dependency_values(len(included.get("positive") or []))
        )
        min_negative_n = int(
            included.get("min_negative_n")
            or min_dependency_values(int(included.get("negative_n") or 0))
        )
        filtered_rows = sorted(
            [
                row
                for row in rows
                if len(row) >= 6 and row[4] >= min_positive_n and row[5] >= min_negative_n
            ],
            key=lambda x: (float(x[1]), str(x[0])),
        )
        for rank, row in enumerate(filtered_rows, start=1):
            if len(row) >= 7:
                row[6] = rank
        analysis["datasets"][dataset_key] = filtered_rows


def min_dependency_values(group_size: int) -> int:
    if group_size <= 0:
        return 1
    return min(
        group_size,
        max(MIN_DEPENDENCY_GROUP_VALUES, math.ceil(group_size * MIN_DEPENDENCY_GROUP_COVERAGE)),
    )


@lru_cache(maxsize=1)
def depmap_cancer_model_count() -> int:
    if not os.path.exists(MODEL_PATH):
        return 0
    with open(MODEL_PATH, newline="") as handle:
        return sum(1 for _ in csv.DictReader(handle))


def add_stratifier_prevalence(analysis: dict[str, Any]) -> None:
    positive_count = len(analysis.get("positive_models", []))
    total_count = int(analysis.get("prevalence_total") or depmap_cancer_model_count())
    analysis["stratifier_prevalence"] = {
        "positive": positive_count,
        "total": total_count,
        "frequency": positive_count / total_count if total_count else 0,
        "denominator": analysis.get("prevalence_denominator")
        or ("all DepMap cancer models" if total_count else "positive cohort only"),
    }


def build_dependency_summary() -> dict[str, Any]:
    summary = load_builtin_dependency_summary()
    summary["analyses"].extend(custom_stratifier_analysis(row) for row in fetch_custom_stratifiers())
    for analysis in summary["analyses"]:
        add_stratifier_prevalence(analysis)
    return summary


def quality_status(analysis: dict[str, Any], quality: dict[str, Any]) -> str:
    positives = len(analysis.get("positive_models", []))
    crispr_n = len(analysis.get("included_models", {}).get("crispr", {}).get("positive", []))
    rnai_n = len(analysis.get("included_models", {}).get("rnai", {}).get("positive", []))
    if positives < 5 or crispr_n < 5 or rnai_n < 3:
        return "Needs review"
    if quality.get("weaknesses"):
        return "Review advised"
    return "Looks usable"


@app.before_request
def hydrate_user() -> None:
    g.user = current_user()


@app.context_processor
def inject_globals() -> dict[str, Any]:
    return {"current_user": g.get("user"), "now_utc": now_utc()}


@app.template_filter("dt")
def format_dt(value: str) -> str:
    if not value:
        return ""
    dt = from_iso(value).astimezone()
    return dt.strftime("%b %d, %Y %I:%M %p")


@app.route("/")
def home():
    if not current_user():
        return redirect(url_for("login"))
    return redirect(url_for("dashboard"))


@app.route("/register", methods=["GET", "POST"])
def register():
    if not app.config["ALLOW_SELF_SIGNUP"]:
        flash("Self-signup is disabled. Ask your admin to create your account.", "error")
        return redirect(url_for("login"))

    if request.method == "POST":
        first_name = request.form.get("first_name", "").strip()
        email = request.form.get("email", "").strip().lower()
        password = request.form.get("password", "")

        if not first_name or not email or not password:
            flash("First name, email, and password are required.", "error")
            return render_template("register.html")

        if len(password) < 8:
            flash("Password must be at least 8 characters.", "error")
            return render_template("register.html")

        db = get_db()
        exists = db_execute(db, "SELECT id FROM users WHERE email=?", (email,)).fetchone()
        if exists:
            flash("Email is already registered.", "error")
            return render_template("register.html")

        db_execute(
            db,
            """
            INSERT INTO users (email, first_name, password_hash, role, is_active, created_at)
            VALUES (?, ?, ?, 'forecaster', 1, ?)
            """,
            (email, first_name, generate_password_hash(password), to_iso(now_utc())),
        )
        db.commit()
        flash("Account created. Please sign in.", "success")
        return redirect(url_for("login"))

    return render_template("register.html")


@app.route("/login", methods=["GET", "POST"])
def login():
    if request.method == "POST":
        email = request.form.get("email", "").strip().lower()
        password = request.form.get("password", "")

        db = get_db()
        user = db_execute(
            db, "SELECT * FROM users WHERE email=? AND is_active=1", (email,)
        ).fetchone()

        if not user or not check_password_hash(user["password_hash"], password):
            flash("Invalid credentials.", "error")
            return render_template("login.html")

        session.clear()
        session["user_id"] = user["id"]
        return redirect(url_for("dashboard"))

    return render_template("login.html")


@app.route("/logout", methods=["POST"])
@login_required
def logout():
    session.clear()
    return redirect(url_for("login"))


@app.route("/dashboard")
@login_required
def dashboard():
    return render_template("dashboard.html")


@app.route("/leaderboard")
@app.route("/bets/new", methods=["GET", "POST"])
@app.route("/bets/<int:bet_id>")
@app.route("/bets/<int:bet_id>/approve", methods=["POST"])
@app.route("/bets/<int:bet_id>/discard", methods=["POST"])
@app.route("/bets/<int:bet_id>/forecast", methods=["POST"])
@app.route("/bets/<int:bet_id>/resolve", methods=["POST"])
def legacy_depmap_redirect(bet_id=None):
    return redirect(url_for("hpv_dependencies"))


@app.route("/admin/users", methods=["GET", "POST"])
@login_required
@roles_required("admin")
def admin_users():
    db = get_db()

    if request.method == "POST":
        action = request.form.get("action", "").strip()

        if action == "create":
            first_name = request.form.get("first_name", "").strip()
            email = request.form.get("email", "").strip().lower()
            password = request.form.get("password", "")
            role = request.form.get("role", "forecaster").strip()
            if role not in {"admin", "forecaster", "viewer"}:
                role = "forecaster"

            if not first_name or not email or len(password) < 8:
                flash("User create failed. Check first name, email, and password length (>=8).", "error")
                return redirect(url_for("admin_users"))

            exists = db_execute(db, "SELECT id FROM users WHERE email=?", (email,)).fetchone()
            if exists:
                flash("Email already exists.", "error")
                return redirect(url_for("admin_users"))

            created_user_id = insert_row(
                db,
                """
                INSERT INTO users (email, first_name, password_hash, role, is_active, created_at)
                VALUES (?, ?, ?, ?, 1, ?)
                """,
                (email, first_name, generate_password_hash(password), role, to_iso(now_utc())),
            )
            db.commit()
            log_event(current_user()["id"], "user_created", "user", created_user_id, {"role": role})
            flash("User created.", "success")

        elif action == "toggle_active":
            target_id_raw = request.form.get("user_id", "")
            try:
                target_id = int(target_id_raw)
            except ValueError:
                flash("Invalid user id.", "error")
                return redirect(url_for("admin_users"))

            target = db_execute(db, "SELECT * FROM users WHERE id=?", (target_id,)).fetchone()
            if not target:
                flash("User not found.", "error")
                return redirect(url_for("admin_users"))

            if target["id"] == current_user()["id"]:
                flash("You cannot deactivate your own account.", "error")
                return redirect(url_for("admin_users"))

            new_state = 0 if target["is_active"] else 1
            db_execute(db, "UPDATE users SET is_active=? WHERE id=?", (new_state, target_id))
            db.commit()
            log_event(current_user()["id"], "user_toggled", "user", target_id, {"is_active": new_state})
            flash("User updated.", "success")

        return redirect(url_for("admin_users"))

    users = db_execute(
        db,
        "SELECT id, first_name, email, role, is_active, created_at FROM users ORDER BY created_at DESC"
    ).fetchall()
    return render_template("admin_users.html", users=users)


@app.route("/scoring")
@login_required
def scoring_explainer():
    return redirect(url_for("hpv_dependencies"))


@app.route("/hpv-dependencies")
@login_required
def hpv_dependencies():
    return render_template("hpv_dependencies.html")


@app.route("/api/dependency-summary")
@login_required
def dependency_summary():
    return build_dependency_summary()


@app.route("/stratifiers")
@login_required
@roles_required("forecaster", "admin")
def manage_stratifiers():
    rows = []
    for row in fetch_custom_stratifiers():
        analysis = json.loads(row["analysis_json"])
        quality = json.loads(row["quality_json"])
        rows.append(
            {
                "id": row["id"],
                "label": row["label"],
                "prompt": row["prompt"],
                "analysis": analysis,
                "source": json.loads(row["source_json"]),
                "quality": quality,
                "quality_status": quality_status(analysis, quality),
                "created_at": row["created_at"],
                "updated_at": row["updated_at"],
            }
        )
    return render_template(
        "stratifiers.html",
        stratifiers=rows,
        openai_ready=bool(os.environ.get("OPENAI_API_KEY")),
        openai_model=os.environ.get("OPENAI_MODEL", "gpt-5.5"),
    )


@app.route("/stratifiers", methods=["POST"])
@login_required
@roles_required("forecaster", "admin")
def create_stratifier():
    prompt = request.form.get("prompt", "").strip()
    if not prompt:
        flash("Describe the cell-line difference you want to stratify.", "error")
        return redirect(url_for("manage_stratifiers"))

    try:
        analysis, source, quality = compute_stratifier(prompt)
    except Exception as exc:
        flash(f"Could not create stratifier: {exc}", "error")
        return redirect(url_for("manage_stratifiers"))

    db = get_db()
    now = to_iso(now_utc())
    db_execute(
        db,
        """
        INSERT INTO dependency_stratifiers
            (label, prompt, analysis_json, source_json, quality_json, created_at, updated_at)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        (
            analysis["label"],
            prompt,
            json.dumps(analysis, separators=(",", ":")),
            json.dumps(source, separators=(",", ":")),
            json.dumps(quality, separators=(",", ":")),
            now,
            now,
        ),
    )
    db.commit()
    flash(f"Added stratifier: {analysis['label']}", "success")
    return redirect(url_for("manage_stratifiers"))


@app.route("/stratifiers/<int:stratifier_id>/delete", methods=["POST"])
@login_required
@roles_required("admin")
def delete_stratifier(stratifier_id: int):
    db = get_db()
    db_execute(db, "DELETE FROM dependency_stratifiers WHERE id=?", (stratifier_id,))
    db.commit()
    flash("Stratifier deleted.", "success")
    return redirect(url_for("manage_stratifiers"))


@app.route("/ret-allostery")
def ret_allostery():
    return redirect(url_for("hpv_dependencies"))


@app.route("/ddx3-selectivity")
def ddx3_selectivity():
    return redirect(url_for("hpv_dependencies"))


init_db()


if __name__ == "__main__":
    app.run(debug=True)
