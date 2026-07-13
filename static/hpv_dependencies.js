const charts = {
  crispr: {
    canvas: document.getElementById("crispr-chart"),
    count: document.getElementById("crispr-count"),
    subtitle: document.getElementById("crispr-subtitle"),
    title: "CRISPR",
    color: "#111111",
    points: [],
  },
  rnai: {
    canvas: document.getElementById("rnai-chart"),
    count: document.getElementById("rnai-count"),
    subtitle: document.getElementById("rnai-subtitle"),
    title: "siRNA",
    color: "#d4a900",
    points: [],
  },
};

const analysisSelect = document.getElementById("analysis-select");
const searchInput = document.getElementById("gene-search");
const resultBody = document.querySelector("#gene-results tbody");
const metaPanel = document.getElementById("hpv-meta");
const tooltip = document.getElementById("chart-tooltip");

let summary = null;
let activeAnalysis = null;
let activeGene = "";

function escapeHtml(value) {
  return String(value || "")
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}

function normalizeRows(rows) {
  return rows.map((row) => ({
    gene: row[0],
    score: row[1],
    positive_average: row[2],
    negative_average: row[3],
    positive_n: row[4],
    negative_n: row[5],
    rank: row[6],
  }));
}

function normalizeSummary(data) {
  data.analyses.forEach((analysis) => {
    analysis.datasets.crispr = normalizeRows(analysis.datasets.crispr);
    analysis.datasets.rnai = normalizeRows(analysis.datasets.rnai);
  });
  return data;
}

function formatScore(value) {
  return Number.isFinite(value) ? value.toFixed(4) : "";
}

function formatPair(entry) {
  if (!entry) {
    return "";
  }
  return `${formatScore(entry.positive_average)} / ${formatScore(entry.negative_average)}`;
}

function resizeCanvas(canvas) {
  const rect = canvas.getBoundingClientRect();
  const ratio = window.devicePixelRatio || 1;
  canvas.width = Math.max(1, Math.round(rect.width * ratio));
  canvas.height = Math.max(1, Math.round(rect.height * ratio));
  const ctx = canvas.getContext("2d");
  ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
  return { width: rect.width, height: rect.height, ctx };
}

function niceTicks(min, max, count = 5) {
  if (min === max) {
    return [min];
  }
  const span = max - min;
  const step = Math.pow(10, Math.floor(Math.log10(span / count)));
  const err = span / count / step;
  const niceStep = err >= 7.5 ? step * 10 : err >= 3.5 ? step * 5 : err >= 1.5 ? step * 2 : step;
  const ticks = [];
  for (let value = Math.ceil(min / niceStep) * niceStep; value <= max + niceStep / 2; value += niceStep) {
    ticks.push(value);
  }
  return ticks;
}

function yLabel() {
  return `${activeAnalysis.positive_label} minus ${activeAnalysis.negative_label}`;
}

function drawChart(key) {
  const chart = charts[key];
  const data = activeAnalysis.datasets[key];
  const { width, height, ctx } = resizeCanvas(chart.canvas);
  const margin = { top: 18, right: 18, bottom: 40, left: 56 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const scores = data.map((d) => d.score);
  const minScore = Math.min(...scores);
  const maxScore = Math.max(...scores);
  const yMin = Math.floor((minScore - 0.05) * 10) / 10;
  const yMax = Math.ceil((maxScore + 0.05) * 10) / 10;
  const highlighted = activeGene ? data.filter((d) => d.gene.toLowerCase().includes(activeGene)) : [];

  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);
  ctx.strokeStyle = "#d9ddd8";
  ctx.lineWidth = 1;
  ctx.strokeRect(margin.left, margin.top, plotWidth, plotHeight);

  ctx.font = "12px IBM Plex Sans, sans-serif";
  ctx.fillStyle = "#5e6d74";
  ctx.textAlign = "right";
  ctx.textBaseline = "middle";
  for (const tick of niceTicks(yMin, yMax)) {
    const y = margin.top + ((yMax - tick) / (yMax - yMin)) * plotHeight;
    ctx.strokeStyle = "#eef1ed";
    ctx.beginPath();
    ctx.moveTo(margin.left, y);
    ctx.lineTo(margin.left + plotWidth, y);
    ctx.stroke();
    ctx.fillText(tick.toFixed(1), margin.left - 9, y);
  }

  if (yMin < 0 && yMax > 0) {
    const zeroY = margin.top + (yMax / (yMax - yMin)) * plotHeight;
    ctx.strokeStyle = "#1f2a30";
    ctx.lineWidth = 1.2;
    ctx.beginPath();
    ctx.moveTo(margin.left, zeroY);
    ctx.lineTo(margin.left + plotWidth, zeroY);
    ctx.stroke();
  }

  ctx.save();
  ctx.translate(15, margin.top + plotHeight / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = "center";
  ctx.fillStyle = "#1f2a30";
  ctx.fillText("Differential dependency", 0, 0);
  ctx.restore();

  ctx.textAlign = "center";
  ctx.fillStyle = "#5e6d74";
  ctx.fillText("Genes ranked by differential score", margin.left + plotWidth / 2, height - 14);

  const sorted = [...data].sort((a, b) => a.rank - b.rank);
  const positions = [];
  ctx.fillStyle = chart.color;
  ctx.globalAlpha = activeGene ? 0.16 : 0.58;
  sorted.forEach((d, i) => {
    const x = margin.left + (i / Math.max(1, sorted.length - 1)) * plotWidth;
    const y = margin.top + ((yMax - d.score) / (yMax - yMin)) * plotHeight;
    positions.push({ x, y, datum: d });
    ctx.beginPath();
    ctx.arc(x, y, 1.55, 0, Math.PI * 2);
    ctx.fill();
  });
  ctx.globalAlpha = 1;

  if (highlighted.length) {
    ctx.fillStyle = "#f2a83b";
    ctx.strokeStyle = "#1f2a30";
    ctx.lineWidth = 1.5;
    highlighted.forEach((d) => {
      const p = positions.find((item) => item.datum.gene === d.gene);
      if (!p) {
        return;
      }
      ctx.beginPath();
      ctx.arc(p.x, p.y, 4.4, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();
    });
  }

  chart.points = positions;
  chart.count.value = `${data.length.toLocaleString()} genes`;
  chart.subtitle.textContent = `${activeAnalysis.positive_label} average minus ${activeAnalysis.negative_label} average`;
}

function renderMeta() {
  const crisprCounts = activeAnalysis.included_models.crispr;
  const rnaiCounts = activeAnalysis.included_models.rnai;
  const positiveModels = activeAnalysis.positive_models || [];
  const modelRows = positiveModels
    .map((m) => {
      const note = m.grouping_note ? ` - ${escapeHtml(m.grouping_note)}` : "";
      return `<li>${escapeHtml(m.cell_line)} <span class="muted">(${escapeHtml(m.model_id)}; ${escapeHtml(m.disease || "unknown disease")}${note})</span></li>`;
    })
    .join("");

  metaPanel.innerHTML = `
    <div>
      <strong>${escapeHtml(activeAnalysis.label)}</strong>
      <span class="small-line">${escapeHtml(activeAnalysis.source)}. Score is ${escapeHtml(activeAnalysis.positive_label)} average minus ${escapeHtml(activeAnalysis.negative_label)} average.</span>
    </div>
    <div>
      <strong>Screened comparison sets</strong>
      <span class="small-line">CRISPR: ${crisprCounts.positive.length} ${escapeHtml(activeAnalysis.positive_label)} vs ${crisprCounts.negative_n} ${escapeHtml(activeAnalysis.negative_label)}. siRNA: ${rnaiCounts.positive.length} ${escapeHtml(activeAnalysis.positive_label)} vs ${rnaiCounts.negative_n} ${escapeHtml(activeAnalysis.negative_label)}.</span>
    </div>
    <details class="model-details">
      <summary>${positiveModels.length.toLocaleString()} positive-group DepMap models</summary>
      <ul>${modelRows || "<li>No positive models matched this analysis.</li>"}</ul>
    </details>
  `;
}

function renderResults() {
  const rnai = new Map(activeAnalysis.datasets.rnai.map((d) => [d.gene, d]));
  const crispr = new Map(activeAnalysis.datasets.crispr.map((d) => [d.gene, d]));
  const genes = new Set([...rnai.keys(), ...crispr.keys()]);
  const rows = [...genes]
    .filter((gene) => (activeGene ? gene.toLowerCase().includes(activeGene) : false))
    .sort()
    .slice(0, 50)
    .map((gene) => {
      const r = rnai.get(gene);
      const c = crispr.get(gene);
      return `
        <tr>
          <td><strong>${escapeHtml(gene)}</strong></td>
          <td>${c ? formatScore(c.score) + ` (#${c.rank})` : ""}</td>
          <td>${formatPair(c)}</td>
          <td>${r ? formatScore(r.score) + ` (#${r.rank})` : ""}</td>
          <td>${formatPair(r)}</td>
        </tr>
      `;
    });

  resultBody.innerHTML = rows.length
    ? rows.join("")
    : '<tr><td colspan="5" class="muted">Type a gene symbol to highlight it in both charts.</td></tr>';
}

function renderAll() {
  renderMeta();
  drawChart("crispr");
  drawChart("rnai");
  renderResults();
}

function populateAnalyses() {
  analysisSelect.innerHTML = summary.analyses
    .map((analysis) => `<option value="${escapeHtml(analysis.id)}">${escapeHtml(analysis.label)}</option>`)
    .join("");
  activeAnalysis = summary.analyses[0];
}

function setAnalysis(analysisId) {
  activeAnalysis = summary.analyses.find((analysis) => analysis.id === analysisId) || summary.analyses[0];
  hideTooltip();
  renderAll();
}

function nearestPoint(chart, x, y) {
  let nearest = null;
  let best = Infinity;
  for (const point of chart.points) {
    const dist = Math.hypot(point.x - x, point.y - y);
    if (dist < best) {
      best = dist;
      nearest = point;
    }
  }
  return best <= 8 ? nearest : null;
}

function showTooltip(event, key, point) {
  const d = point.datum;
  tooltip.innerHTML = `
    <strong>${escapeHtml(d.gene)}</strong>
    <span>${charts[key].title} rank #${d.rank.toLocaleString()}</span>
    <span>Diff: ${formatScore(d.score)}</span>
    <span>${escapeHtml(activeAnalysis.positive_label)}: ${formatScore(d.positive_average)} (n=${d.positive_n})</span>
    <span>${escapeHtml(activeAnalysis.negative_label)}: ${formatScore(d.negative_average)} (n=${d.negative_n})</span>
  `;
  tooltip.classList.remove("hidden");
  tooltip.style.left = `${event.clientX + 14}px`;
  tooltip.style.top = `${event.clientY + 14}px`;
}

function hideTooltip() {
  tooltip.classList.add("hidden");
}

for (const [key, chart] of Object.entries(charts)) {
  chart.canvas.addEventListener("mousemove", (event) => {
    if (!activeAnalysis) {
      return;
    }
    const rect = chart.canvas.getBoundingClientRect();
    const point = nearestPoint(chart, event.clientX - rect.left, event.clientY - rect.top);
    if (point) {
      showTooltip(event, key, point);
    } else {
      hideTooltip();
    }
  });
  chart.canvas.addEventListener("mouseleave", hideTooltip);
}

fetch("/api/dependency-summary")
  .then((res) => {
    if (!res.ok) {
      throw new Error("Missing data summary");
    }
    return res.json();
  })
  .then((data) => {
    summary = normalizeSummary(data);
    populateAnalyses();
    renderAll();
  })
  .catch(() => {
    metaPanel.innerHTML = '<div><strong>Data summary not found</strong><span class="small-line">Run scripts/build_hpv_dependency_data.py.</span></div>';
  });

analysisSelect.addEventListener("change", (event) => {
  setAnalysis(event.target.value);
});

searchInput.addEventListener("input", (event) => {
  activeGene = event.target.value.trim().toLowerCase();
  if (summary) {
    renderAll();
  }
});

window.addEventListener("resize", () => {
  if (summary) {
    renderAll();
  }
});
