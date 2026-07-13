import json
import os
import sqlite3
from datetime import datetime, timedelta, timezone
from functools import wraps
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

from dependency_stratifiers import compute_stratifier


BASE_DIR = os.path.abspath(os.path.dirname(__file__))
DB_PATH = os.environ.get("DB_PATH", os.path.join(BASE_DIR, "app.db"))
DEPENDENCY_SUMMARY_PATH = os.path.join(BASE_DIR, "static", "data", "hpv_dependency_summary.json")
UTC = timezone.utc


app = Flask(__name__)
app.config["SECRET_KEY"] = os.environ.get("SECRET_KEY", "dev-secret-change-me")
app.config["ALLOW_SELF_SIGNUP"] = os.environ.get("ALLOW_SELF_SIGNUP", "true").lower() == "true"


def now_utc() -> datetime:
    return datetime.now(tz=UTC)


def to_iso(dt: datetime) -> str:
    return dt.astimezone(UTC).replace(microsecond=0).isoformat()


def from_iso(value: str) -> datetime:
    return datetime.fromisoformat(value).astimezone(UTC)


def parse_local_datetime(value: str) -> datetime:
    # HTML datetime-local has no timezone; interpret as local server time and convert to UTC.
    local_dt = datetime.strptime(value, "%Y-%m-%dT%H:%M")
    return local_dt.astimezone()


def clamp(num: float, low: float = 0.0, high: float = 1.0) -> float:
    return max(low, min(high, num))


def get_db() -> sqlite3.Connection:
    if "db" not in g:
        conn = sqlite3.connect(DB_PATH)
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA foreign_keys = ON")
        g.db = conn
    return g.db


@app.teardown_appcontext
def close_db(exception: Exception | None) -> None:
    db = g.pop("db", None)
    if db is not None:
        db.close()


def init_db() -> None:
    db_dir = os.path.dirname(DB_PATH)
    if db_dir:
        os.makedirs(db_dir, exist_ok=True)

    db = sqlite3.connect(DB_PATH)
    db.execute("PRAGMA foreign_keys = ON")
    db.executescript(
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

        CREATE TABLE IF NOT EXISTS bets (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            title TEXT NOT NULL,
            description TEXT,
            category TEXT,
            market_type TEXT NOT NULL CHECK(market_type IN ('binary', 'numeric', 'multiple_choice')),
            options_json TEXT,
            creator_id INTEGER NOT NULL,
            status TEXT NOT NULL CHECK(status IN ('proposed', 'active', 'resolved', 'discarded', 'expired')) DEFAULT 'proposed',
            close_at TEXT NOT NULL,
            resolve_at TEXT NOT NULL,
            approval_deadline TEXT NOT NULL,
            allow_multiple_resolutions INTEGER NOT NULL DEFAULT 1,
            discarded_reason TEXT,
            created_at TEXT NOT NULL,
            FOREIGN KEY (creator_id) REFERENCES users(id)
        );

        CREATE TABLE IF NOT EXISTS bet_approvals (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            bet_id INTEGER NOT NULL,
            user_id INTEGER NOT NULL,
            approved_at TEXT NOT NULL,
            UNIQUE (bet_id, user_id),
            FOREIGN KEY (bet_id) REFERENCES bets(id) ON DELETE CASCADE,
            FOREIGN KEY (user_id) REFERENCES users(id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS forecasts (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            bet_id INTEGER NOT NULL,
            user_id INTEGER NOT NULL,
            choice TEXT,
            probability_yes REAL,
            point_estimate REAL,
            confidence REAL NOT NULL,
            created_at TEXT NOT NULL,
            updated_at TEXT NOT NULL,
            locked_at TEXT NOT NULL,
            UNIQUE (bet_id, user_id),
            FOREIGN KEY (bet_id) REFERENCES bets(id) ON DELETE CASCADE,
            FOREIGN KEY (user_id) REFERENCES users(id) ON DELETE CASCADE
        );

        CREATE TABLE IF NOT EXISTS resolutions (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            bet_id INTEGER NOT NULL,
            label TEXT NOT NULL,
            outcome_binary INTEGER,
            outcome_numeric REAL,
            outcome_option TEXT,
            weight REAL NOT NULL DEFAULT 1.0,
            source_url TEXT,
            notes TEXT,
            is_void INTEGER NOT NULL DEFAULT 0,
            resolver_id INTEGER NOT NULL,
            resolved_at TEXT NOT NULL,
            FOREIGN KEY (bet_id) REFERENCES bets(id) ON DELETE CASCADE,
            FOREIGN KEY (resolver_id) REFERENCES users(id)
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
        );
        """
    )

    admin_email = os.environ.get("ADMIN_EMAIL")
    admin_password = os.environ.get("ADMIN_PASSWORD")
    admin_name = os.environ.get("ADMIN_FIRST_NAME", "Admin")

    existing_admin = db.execute("SELECT id FROM users WHERE role='admin' AND is_active=1 LIMIT 1").fetchone()
    if not existing_admin and admin_email and admin_password:
        db.execute(
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
    db.execute(
        """
        INSERT INTO audit_logs (actor_id, event_type, entity_type, entity_id, payload_json, created_at)
        VALUES (?, ?, ?, ?, ?, ?)
        """,
        (actor_id, event_type, entity_type, entity_id, json.dumps(payload), to_iso(now_utc())),
    )
    db.commit()


def refresh_bet_statuses() -> None:
    db = get_db()
    ts = to_iso(now_utc())

    db.execute(
        """
        UPDATE bets
        SET status='expired'
        WHERE status='proposed' AND approval_deadline < ?
        """,
        (ts,),
    )

    db.execute(
        """
        UPDATE bets
        SET status='active'
        WHERE status='proposed'
          AND approval_deadline >= ?
          AND EXISTS (
            SELECT 1
            FROM bet_approvals ba
            WHERE ba.bet_id = bets.id AND ba.user_id != bets.creator_id
          )
        """,
        (ts,),
    )

    db.commit()


def get_user(user_id: int) -> sqlite3.Row | None:
    db = get_db()
    return db.execute("SELECT * FROM users WHERE id=? AND is_active=1", (user_id,)).fetchone()


def current_user() -> sqlite3.Row | None:
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


def can_approve_bet(user: sqlite3.Row, bet: sqlite3.Row) -> bool:
    if user["role"] not in {"forecaster", "admin"}:
        return False
    if user["id"] == bet["creator_id"]:
        return False
    if bet["status"] != "proposed":
        return False
    return from_iso(bet["approval_deadline"]) >= now_utc()


def get_forecast_for_user(bet_id: int, user_id: int) -> sqlite3.Row | None:
    db = get_db()
    return db.execute("SELECT * FROM forecasts WHERE bet_id=? AND user_id=?", (bet_id, user_id)).fetchone()


def can_view_forecasts(user: sqlite3.Row, bet_id: int) -> bool:
    if user["role"] in {"admin", "viewer"}:
        return True
    own = get_forecast_for_user(bet_id, user["id"])
    return own is not None


def score_binary(prob_yes: float, outcome: int) -> float:
    brier = (prob_yes - float(outcome)) ** 2
    return clamp(1.0 - brier) * 100.0


def score_numeric(point_estimate: float, confidence: float, actual: float) -> float:
    scale = max(abs(actual), 1.0)
    relative_error = abs(point_estimate - actual) / scale
    distance_component = clamp(1.0 - relative_error)

    # Confidence is treated as probability that estimate lands in +/-10% band.
    hit = 1.0 if relative_error <= 0.10 else 0.0
    confidence_component = clamp(1.0 - ((confidence - hit) ** 2))

    score = (0.65 * distance_component) + (0.35 * confidence_component)
    return clamp(score) * 100.0


def score_multiple(choice: str, confidence: float, outcome_option: str) -> float:
    outcome = 1 if choice == outcome_option else 0
    return score_binary(confidence, outcome)


def resolve_score_for_bet(bet: sqlite3.Row, forecast: sqlite3.Row, resolutions: list[sqlite3.Row]) -> float | None:
    active_resolutions = [r for r in resolutions if r["is_void"] == 0]
    if not active_resolutions:
        return None

    weighted_total = 0.0
    weight_sum = 0.0

    for res in active_resolutions:
        weight = float(res["weight"] or 1.0)
        market_type = bet["market_type"]

        if market_type == "binary":
            if res["outcome_binary"] is None or forecast["probability_yes"] is None:
                continue
            sc = score_binary(float(forecast["probability_yes"]), int(res["outcome_binary"]))
        elif market_type == "numeric":
            if res["outcome_numeric"] is None or forecast["point_estimate"] is None:
                continue
            sc = score_numeric(float(forecast["point_estimate"]), float(forecast["confidence"]), float(res["outcome_numeric"]))
        else:
            if not res["outcome_option"] or not forecast["choice"]:
                continue
            sc = score_multiple(str(forecast["choice"]), float(forecast["confidence"]), str(res["outcome_option"]))

        weighted_total += (sc * weight)
        weight_sum += weight

    if weight_sum == 0:
        return None
    return round(weighted_total / weight_sum, 2)


def fetch_leaderboard() -> list[dict[str, Any]]:
    db = get_db()

    users = db.execute(
        """
        SELECT id, first_name, role
        FROM users
        WHERE is_active=1 AND role IN ('forecaster', 'admin')
        ORDER BY first_name
        """
    ).fetchall()

    resolved_bets = db.execute("SELECT * FROM bets WHERE status='resolved'").fetchall()
    resolutions_by_bet: dict[int, list[sqlite3.Row]] = {}
    for bet in resolved_bets:
        rows = db.execute("SELECT * FROM resolutions WHERE bet_id=?", (bet["id"],)).fetchall()
        resolutions_by_bet[bet["id"]] = rows

    leaderboard: list[dict[str, Any]] = []

    for user in users:
        scores: list[float] = []
        for bet in resolved_bets:
            forecast = db.execute(
                "SELECT * FROM forecasts WHERE bet_id=? AND user_id=?",
                (bet["id"], user["id"]),
            ).fetchone()
            if not forecast:
                continue
            score = resolve_score_for_bet(bet, forecast, resolutions_by_bet.get(bet["id"], []))
            if score is not None:
                scores.append(score)

        avg_score = round(sum(scores) / len(scores), 2) if scores else 0.0
        leaderboard.append(
            {
                "user_id": user["id"],
                "first_name": user["first_name"],
                "role": user["role"],
                "avg_score": avg_score,
                "resolved_count": len(scores),
            }
        )

    leaderboard.sort(key=lambda x: (x["avg_score"], x["resolved_count"]), reverse=True)
    for idx, entry in enumerate(leaderboard, start=1):
        entry["rank"] = idx
    return leaderboard


def load_builtin_dependency_summary() -> dict[str, Any]:
    with open(DEPENDENCY_SUMMARY_PATH) as f:
        return json.load(f)


def fetch_custom_stratifiers() -> list[sqlite3.Row]:
    return get_db().execute(
        """
        SELECT *
        FROM dependency_stratifiers
        ORDER BY datetime(updated_at) DESC
        """
    ).fetchall()


def custom_stratifier_analysis(row: sqlite3.Row) -> dict[str, Any]:
    analysis = json.loads(row["analysis_json"])
    analysis["id"] = f"custom-{row['id']}"
    analysis["label"] = f"{analysis['label']} *"
    analysis["custom_id"] = row["id"]
    return analysis


def build_dependency_summary() -> dict[str, Any]:
    summary = load_builtin_dependency_summary()
    summary["analyses"].extend(custom_stratifier_analysis(row) for row in fetch_custom_stratifiers())
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
    if g.user:
        refresh_bet_statuses()


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
        exists = db.execute("SELECT id FROM users WHERE email=?", (email,)).fetchone()
        if exists:
            flash("Email is already registered.", "error")
            return render_template("register.html")

        db.execute(
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
        user = db.execute("SELECT * FROM users WHERE email=? AND is_active=1", (email,)).fetchone()

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
    db = get_db()

    proposed = db.execute(
        """
        SELECT b.*, u.first_name AS creator_name,
               (SELECT COUNT(*) FROM bet_approvals ba WHERE ba.bet_id = b.id) AS approval_count
        FROM bets b
        JOIN users u ON u.id = b.creator_id
        WHERE b.status='proposed'
        ORDER BY b.created_at DESC
        """
    ).fetchall()

    active = db.execute(
        """
        SELECT b.*, u.first_name AS creator_name
        FROM bets b
        JOIN users u ON u.id = b.creator_id
        WHERE b.status='active'
        ORDER BY b.close_at ASC
        """
    ).fetchall()

    resolved = db.execute(
        """
        SELECT b.*, u.first_name AS creator_name
        FROM bets b
        JOIN users u ON u.id = b.creator_id
        WHERE b.status='resolved'
        ORDER BY b.resolve_at DESC
        LIMIT 20
        """
    ).fetchall()

    leaderboard = fetch_leaderboard()

    return render_template(
        "dashboard.html",
        proposed=proposed,
        active=active,
        resolved=resolved,
        leaderboard=leaderboard,
    )


@app.route("/leaderboard")
@login_required
def leaderboard():
    leaderboard_data = fetch_leaderboard()
    return render_template("leaderboard.html", leaderboard=leaderboard_data)


@app.route("/bets/new", methods=["GET", "POST"])
@login_required
@roles_required("forecaster", "admin")
def create_bet():
    if request.method == "POST":
        title = request.form.get("title", "").strip()
        description = request.form.get("description", "").strip()
        category = request.form.get("category", "").strip()
        market_type = request.form.get("market_type", "").strip()
        close_at_raw = request.form.get("close_at", "")
        resolve_at_raw = request.form.get("resolve_at", "")
        options_raw = request.form.get("options", "").strip()

        if not title or market_type not in {"binary", "numeric", "multiple_choice"}:
            flash("Title and a valid market type are required.", "error")
            return render_template("bet_form.html")

        try:
            close_at = parse_local_datetime(close_at_raw)
            resolve_at = parse_local_datetime(resolve_at_raw)
        except ValueError:
            flash("Invalid close or resolve date.", "error")
            return render_template("bet_form.html")

        now = now_utc()
        if close_at <= now:
            flash("Close time must be in the future.", "error")
            return render_template("bet_form.html")
        if resolve_at <= close_at:
            flash("Resolve time must be after close time.", "error")
            return render_template("bet_form.html")

        options_json = None
        if market_type == "multiple_choice":
            options = [line.strip() for line in options_raw.splitlines() if line.strip()]
            if len(options) < 2:
                flash("Multiple choice bets need at least two options.", "error")
                return render_template("bet_form.html")
            options_json = json.dumps(options)

        approval_deadline = min(now + timedelta(days=5), close_at)

        db = get_db()
        cur = db.execute(
            """
            INSERT INTO bets (
                title, description, category, market_type, options_json,
                creator_id, status, close_at, resolve_at, approval_deadline,
                allow_multiple_resolutions, created_at
            )
            VALUES (?, ?, ?, ?, ?, ?, 'proposed', ?, ?, ?, 1, ?)
            """,
            (
                title,
                description,
                category,
                market_type,
                options_json,
                current_user()["id"],
                to_iso(close_at),
                to_iso(resolve_at),
                to_iso(approval_deadline),
                to_iso(now),
            ),
        )
        bet_id = cur.lastrowid
        db.commit()

        log_event(current_user()["id"], "bet_created", "bet", bet_id, {"title": title, "market_type": market_type})

        flash("Bet proposed. It becomes active after one non-creator approval.", "success")
        return redirect(url_for("bet_detail", bet_id=bet_id))

    return render_template("bet_form.html")


@app.route("/bets/<int:bet_id>")
@login_required
def bet_detail(bet_id: int):
    db = get_db()
    bet = db.execute(
        """
        SELECT b.*, u.first_name AS creator_name
        FROM bets b
        JOIN users u ON u.id = b.creator_id
        WHERE b.id=?
        """,
        (bet_id,),
    ).fetchone()

    if not bet:
        flash("Bet not found.", "error")
        return redirect(url_for("dashboard"))

    approvals = db.execute(
        """
        SELECT ba.*, u.first_name
        FROM bet_approvals ba
        JOIN users u ON u.id = ba.user_id
        WHERE ba.bet_id=?
        ORDER BY ba.approved_at ASC
        """,
        (bet_id,),
    ).fetchall()

    own_forecast = get_forecast_for_user(bet_id, current_user()["id"])

    forecasts = []
    if can_view_forecasts(current_user(), bet_id):
        forecasts = db.execute(
            """
            SELECT f.*, u.first_name
            FROM forecasts f
            JOIN users u ON u.id = f.user_id
            WHERE f.bet_id=?
            ORDER BY f.updated_at DESC
            """,
            (bet_id,),
        ).fetchall()

    resolutions = db.execute(
        """
        SELECT r.*, u.first_name AS resolver_name
        FROM resolutions r
        JOIN users u ON u.id = r.resolver_id
        WHERE r.bet_id=?
        ORDER BY r.resolved_at DESC
        """,
        (bet_id,),
    ).fetchall()

    options = json.loads(bet["options_json"]) if bet["options_json"] else []
    now = now_utc()
    close_dt = from_iso(bet["close_at"])

    can_submit_forecast = (
        current_user()["role"] in {"forecaster", "admin"}
        and bet["status"] == "active"
        and now <= close_dt
    )

    can_edit = False
    if own_forecast:
        can_edit = now <= from_iso(own_forecast["locked_at"]) and now <= close_dt and bet["status"] == "active"

    return render_template(
        "bet_detail.html",
        bet=bet,
        approvals=approvals,
        own_forecast=own_forecast,
        forecasts=forecasts,
        resolutions=resolutions,
        options=options,
        can_approve=can_approve_bet(current_user(), bet),
        can_submit_forecast=can_submit_forecast,
        can_edit=can_edit,
        can_view_forecast_list=can_view_forecasts(current_user(), bet_id),
    )


@app.route("/bets/<int:bet_id>/approve", methods=["POST"])
@login_required
@roles_required("forecaster", "admin")
def approve_bet(bet_id: int):
    db = get_db()
    bet = db.execute("SELECT * FROM bets WHERE id=?", (bet_id,)).fetchone()
    if not bet:
        flash("Bet not found.", "error")
        return redirect(url_for("dashboard"))

    if not can_approve_bet(current_user(), bet):
        flash("You cannot approve this bet.", "error")
        return redirect(url_for("bet_detail", bet_id=bet_id))

    exists = db.execute(
        "SELECT id FROM bet_approvals WHERE bet_id=? AND user_id=?",
        (bet_id, current_user()["id"]),
    ).fetchone()
    if exists:
        flash("You already approved this bet.", "error")
        return redirect(url_for("bet_detail", bet_id=bet_id))

    db.execute(
        "INSERT INTO bet_approvals (bet_id, user_id, approved_at) VALUES (?, ?, ?)",
        (bet_id, current_user()["id"], to_iso(now_utc())),
    )

    db.execute(
        """
        UPDATE bets
        SET status='active'
        WHERE id=? AND status='proposed'
        """,
        (bet_id,),
    )

    db.commit()

    log_event(current_user()["id"], "bet_approved", "bet", bet_id, {})

    flash("Bet approved and activated.", "success")
    return redirect(url_for("bet_detail", bet_id=bet_id))


@app.route("/bets/<int:bet_id>/discard", methods=["POST"])
@login_required
@roles_required("admin")
def discard_bet(bet_id: int):
    reason = request.form.get("reason", "").strip()
    db = get_db()
    bet = db.execute("SELECT * FROM bets WHERE id=?", (bet_id,)).fetchone()
    if not bet:
        flash("Bet not found.", "error")
        return redirect(url_for("dashboard"))

    if bet["status"] in {"resolved", "discarded"}:
        flash("Bet cannot be discarded in its current state.", "error")
        return redirect(url_for("bet_detail", bet_id=bet_id))

    db.execute("UPDATE bets SET status='discarded', discarded_reason=? WHERE id=?", (reason, bet_id))
    db.commit()

    log_event(current_user()["id"], "bet_discarded", "bet", bet_id, {"reason": reason})

    flash("Bet discarded.", "success")
    return redirect(url_for("bet_detail", bet_id=bet_id))


@app.route("/bets/<int:bet_id>/forecast", methods=["POST"])
@login_required
@roles_required("forecaster", "admin")
def submit_forecast(bet_id: int):
    db = get_db()
    bet = db.execute("SELECT * FROM bets WHERE id=?", (bet_id,)).fetchone()
    if not bet:
        flash("Bet not found.", "error")
        return redirect(url_for("dashboard"))

    if bet["status"] != "active":
        flash("This bet is not active.", "error")
        return redirect(url_for("bet_detail", bet_id=bet_id))

    now = now_utc()
    if now > from_iso(bet["close_at"]):
        flash("Forecasting window has closed.", "error")
        return redirect(url_for("bet_detail", bet_id=bet_id))

    confidence_raw = request.form.get("confidence", "")
    try:
        confidence = float(confidence_raw) / 100.0
    except ValueError:
        flash("Invalid confidence value.", "error")
        return redirect(url_for("bet_detail", bet_id=bet_id))

    confidence = clamp(confidence, 0.5, 1.0)

    choice = None
    probability_yes = None
    point_estimate = None

    if bet["market_type"] == "binary":
        choice = request.form.get("choice", "").strip().lower()
        if choice not in {"yes", "no"}:
            flash("Binary bets require Yes or No.", "error")
            return redirect(url_for("bet_detail", bet_id=bet_id))
        probability_yes = confidence if choice == "yes" else (1.0 - confidence)

    elif bet["market_type"] == "numeric":
        point_raw = request.form.get("point_estimate", "")
        try:
            point_estimate = float(point_raw)
        except ValueError:
            flash("Numeric bets require a valid point estimate.", "error")
            return redirect(url_for("bet_detail", bet_id=bet_id))

    else:
        choice = request.form.get("choice", "").strip()
        options = json.loads(bet["options_json"] or "[]")
        if choice not in options:
            flash("Invalid option selected.", "error")
            return redirect(url_for("bet_detail", bet_id=bet_id))

    existing = get_forecast_for_user(bet_id, current_user()["id"])

    if existing:
        if now > from_iso(existing["locked_at"]):
            flash("Your 10-minute edit window has ended.", "error")
            return redirect(url_for("bet_detail", bet_id=bet_id))

        db.execute(
            """
            UPDATE forecasts
            SET choice=?, probability_yes=?, point_estimate=?, confidence=?, updated_at=?
            WHERE id=?
            """,
            (choice, probability_yes, point_estimate, confidence, to_iso(now), existing["id"]),
        )
        db.commit()
        log_event(current_user()["id"], "forecast_updated", "bet", bet_id, {"user_id": current_user()["id"]})
        flash("Forecast updated.", "success")
    else:
        lock_time = min(now + timedelta(minutes=10), from_iso(bet["close_at"]))
        db.execute(
            """
            INSERT INTO forecasts (
                bet_id, user_id, choice, probability_yes, point_estimate,
                confidence, created_at, updated_at, locked_at
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                bet_id,
                current_user()["id"],
                choice,
                probability_yes,
                point_estimate,
                confidence,
                to_iso(now),
                to_iso(now),
                to_iso(lock_time),
            ),
        )
        db.commit()
        log_event(current_user()["id"], "forecast_submitted", "bet", bet_id, {"user_id": current_user()["id"]})
        flash("Forecast submitted.", "success")

    return redirect(url_for("bet_detail", bet_id=bet_id))


@app.route("/bets/<int:bet_id>/resolve", methods=["POST"])
@login_required
@roles_required("admin")
def resolve_bet(bet_id: int):
    db = get_db()
    bet = db.execute("SELECT * FROM bets WHERE id=?", (bet_id,)).fetchone()
    if not bet:
        flash("Bet not found.", "error")
        return redirect(url_for("dashboard"))

    label = request.form.get("label", "").strip() or "Primary"
    source_url = request.form.get("source_url", "").strip()
    notes = request.form.get("notes", "").strip()
    weight_raw = request.form.get("weight", "1")
    void_raw = request.form.get("is_void", "0")

    try:
        weight = max(0.1, float(weight_raw))
    except ValueError:
        weight = 1.0

    is_void = 1 if void_raw == "1" else 0

    outcome_binary = None
    outcome_numeric = None
    outcome_option = None

    if bet["market_type"] == "binary":
        resolved_value = request.form.get("outcome_binary", "").strip().lower()
        if resolved_value not in {"yes", "no"}:
            flash("Binary resolution requires Yes or No.", "error")
            return redirect(url_for("bet_detail", bet_id=bet_id))
        outcome_binary = 1 if resolved_value == "yes" else 0

    elif bet["market_type"] == "numeric":
        numeric_raw = request.form.get("outcome_numeric", "").strip()
        try:
            outcome_numeric = float(numeric_raw)
        except ValueError:
            flash("Numeric resolution requires a valid number.", "error")
            return redirect(url_for("bet_detail", bet_id=bet_id))

    else:
        selected = request.form.get("outcome_option", "").strip()
        options = json.loads(bet["options_json"] or "[]")
        if selected not in options:
            flash("Resolution option must match one of the configured options.", "error")
            return redirect(url_for("bet_detail", bet_id=bet_id))
        outcome_option = selected

    db.execute(
        """
        INSERT INTO resolutions (
            bet_id, label, outcome_binary, outcome_numeric, outcome_option, weight,
            source_url, notes, is_void, resolver_id, resolved_at
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            bet_id,
            label,
            outcome_binary,
            outcome_numeric,
            outcome_option,
            weight,
            source_url,
            notes,
            is_void,
            current_user()["id"],
            to_iso(now_utc()),
        ),
    )

    db.execute("UPDATE bets SET status='resolved' WHERE id=?", (bet_id,))
    db.commit()

    log_event(
        current_user()["id"],
        "bet_resolved",
        "bet",
        bet_id,
        {
            "label": label,
            "is_void": is_void,
            "weight": weight,
        },
    )

    flash("Resolution recorded. You can add more resolutions if needed.", "success")
    return redirect(url_for("bet_detail", bet_id=bet_id))


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

            exists = db.execute("SELECT id FROM users WHERE email=?", (email,)).fetchone()
            if exists:
                flash("Email already exists.", "error")
                return redirect(url_for("admin_users"))

            cur = db.execute(
                """
                INSERT INTO users (email, first_name, password_hash, role, is_active, created_at)
                VALUES (?, ?, ?, ?, 1, ?)
                """,
                (email, first_name, generate_password_hash(password), role, to_iso(now_utc())),
            )
            db.commit()
            log_event(current_user()["id"], "user_created", "user", cur.lastrowid, {"role": role})
            flash("User created.", "success")

        elif action == "toggle_active":
            target_id_raw = request.form.get("user_id", "")
            try:
                target_id = int(target_id_raw)
            except ValueError:
                flash("Invalid user id.", "error")
                return redirect(url_for("admin_users"))

            target = db.execute("SELECT * FROM users WHERE id=?", (target_id,)).fetchone()
            if not target:
                flash("User not found.", "error")
                return redirect(url_for("admin_users"))

            if target["id"] == current_user()["id"]:
                flash("You cannot deactivate your own account.", "error")
                return redirect(url_for("admin_users"))

            new_state = 0 if target["is_active"] else 1
            db.execute("UPDATE users SET is_active=? WHERE id=?", (new_state, target_id))
            db.commit()
            log_event(current_user()["id"], "user_toggled", "user", target_id, {"is_active": new_state})
            flash("User updated.", "success")

        return redirect(url_for("admin_users"))

    users = db.execute(
        "SELECT id, first_name, email, role, is_active, created_at FROM users ORDER BY created_at DESC"
    ).fetchall()
    return render_template("admin_users.html", users=users)


@app.route("/scoring")
@login_required
def scoring_explainer():
    return render_template("scoring.html")


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
    return render_template("stratifiers.html", stratifiers=rows)


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
    db.execute(
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
    db.execute("DELETE FROM dependency_stratifiers WHERE id=?", (stratifier_id,))
    db.commit()
    flash("Stratifier deleted.", "success")
    return redirect(url_for("manage_stratifiers"))


@app.route("/ret-allostery")
def ret_allostery():
    return render_template("ret_allostery.html")


@app.route("/ddx3-selectivity")
def ddx3_selectivity():
    return render_template("ddx3_selectivity.html")


init_db()


if __name__ == "__main__":
    app.run(debug=True)
