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
const geneOptions = document.getElementById("gene-options");
const negativeDisplayMode = document.getElementById("negative-display-mode");
const cutoffInput = document.getElementById("negative-cutoff");
const cutoffValue = document.getElementById("negative-cutoff-value");
const cutoffSettings = document.getElementById("cutoff-settings");
const negativeColorHelp = document.getElementById("negative-color-help");
const resultBody = document.querySelector("#gene-results tbody");
const metaPanel = document.getElementById("hpv-meta");
const stratifierPrevalencePanel = document.getElementById("stratifier-prevalence-panel");
const tooltip = document.getElementById("chart-tooltip");

let summary = null;
let activeAnalysis = null;
let activeGene = "";
let geneQuery = "";
let analysisLoadToken = 0;
let negativeEssentialityCutoff = Number(cutoffInput.value);

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
    raw_difference: Number.isFinite(row[7]) ? row[7] : row[2] - row[3],
  }));
}

function normalizeAnalysis(analysis) {
  if (analysis.datasets && !analysis.datasets_normalized) {
    analysis.datasets.crispr = normalizeRows(analysis.datasets.crispr);
    analysis.datasets.rnai = normalizeRows(analysis.datasets.rnai);
    analysis.datasets_normalized = true;
  }
  return analysis;
}

function normalizeSummary(data) {
  data.analyses.forEach(normalizeAnalysis);
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

function visibleRows(key) {
  const rows = activeAnalysis.datasets[key];
  if (negativeDisplayMode.value === "color") {
    return rows;
  }
  return rows.filter((row) => row.negative_average > negativeEssentialityCutoff);
}

function negativeDependencyColor(value) {
  const normalized = Math.max(0, Math.min(1, (value + 1.5) / 1.5));
  const hue = 4 + normalized * 136;
  const lightness = 43 + normalized * 5;
  return `hsl(${hue}, 66%, ${lightness}%)`;
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
  const error = span / count / step;
  const niceStep = error >= 7.5 ? step * 10 : error >= 3.5 ? step * 5 : error >= 1.5 ? step * 2 : step;
  const ticks = [];
  for (let value = Math.ceil(min / niceStep) * niceStep; value <= max + niceStep / 2; value += niceStep) {
    ticks.push(value);
  }
  return ticks;
}

function paddedExtent(values, fallbackMin = -1, fallbackMax = 1) {
  if (!values.length) {
    return [fallbackMin, fallbackMax];
  }
  const min = Math.min(...values);
  const max = Math.max(...values);
  const padding = Math.max(0.05, (max - min) * 0.06);
  return [min - padding, max + padding];
}

function xLabel() {
  return `Absolute dependency in ${activeAnalysis.positive_label}`;
}

function yLabel() {
  if (activeAnalysis.effect_metric === "hedges_g") {
    return `Hedges' g: ${activeAnalysis.positive_label} vs ${activeAnalysis.negative_label}`;
  }
  return `${activeAnalysis.positive_label} minus ${activeAnalysis.negative_label}`;
}

function drawQuadrants(ctx, margin, plotWidth, plotHeight, cutoffX, zeroY) {
  ctx.save();
  ctx.beginPath();
  ctx.rect(margin.left, margin.top, plotWidth, plotHeight);
  ctx.clip();
  ctx.fillStyle = "rgba(245, 196, 0, 0.07)";
  ctx.fillRect(margin.left, margin.top, Math.max(0, cutoffX - margin.left), Math.max(0, zeroY - margin.top));
  ctx.fillStyle = "rgba(35, 122, 87, 0.06)";
  ctx.fillRect(margin.left, zeroY, Math.max(0, cutoffX - margin.left), Math.max(0, margin.top + plotHeight - zeroY));
  ctx.fillStyle = "rgba(180, 58, 58, 0.04)";
  ctx.fillRect(cutoffX, margin.top, Math.max(0, margin.left + plotWidth - cutoffX), Math.max(0, zeroY - margin.top));
  ctx.restore();
}

function drawNegativeDependencyLegend(ctx, width) {
  const left = Math.max(70, width - 214);
  const top = 16;
  const gradient = ctx.createLinearGradient(left, top, left + 112, top);
  gradient.addColorStop(0, negativeDependencyColor(-1.5));
  gradient.addColorStop(0.67, negativeDependencyColor(-0.5));
  gradient.addColorStop(1, negativeDependencyColor(0));
  ctx.font = "11px IBM Plex Sans, sans-serif";
  ctx.fillStyle = "#5e6d74";
  ctx.textAlign = "left";
  ctx.textBaseline = "top";
  ctx.fillText("Negative cohort dependency", left, 2);
  ctx.fillStyle = gradient;
  ctx.fillRect(left, top, 112, 8);
  ctx.strokeStyle = "#b8bbb6";
  ctx.strokeRect(left, top, 112, 8);
  ctx.fillStyle = "#5e6d74";
  ctx.textAlign = "left";
  ctx.textBaseline = "top";
  ctx.fillText("<= -1.5", left, top + 11);
  ctx.textAlign = "right";
  ctx.fillText(">= 0", left + 112, top + 11);
}

function drawChart(key) {
  const chart = charts[key];
  const allData = activeAnalysis.datasets[key];
  const data = visibleRows(key);
  const { width, height, ctx } = resizeCanvas(chart.canvas);
  const margin = { top: negativeDisplayMode.value === "color" ? 44 : 28, right: 18, bottom: 54, left: 58 };
  const plotWidth = width - margin.left - margin.right;
  const plotHeight = height - margin.top - margin.bottom;
  const [xMin, xMax] = paddedExtent(
    [...allData.map((d) => d.positive_average), negativeEssentialityCutoff],
    -1,
    0
  );
  const [yMin, yMax] = paddedExtent(allData.map((d) => d.score), -1, 1);
  const xAt = (value) => margin.left + ((value - xMin) / (xMax - xMin)) * plotWidth;
  const yAt = (value) => margin.top + ((yMax - value) / (yMax - yMin)) * plotHeight;
  const cutoffX = xAt(negativeEssentialityCutoff);
  const zeroY = yAt(0);
  const highlighted = activeGene
    ? data.filter((d) => d.gene.toLowerCase() === activeGene.toLowerCase())
    : [];

  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = "#ffffff";
  ctx.fillRect(0, 0, width, height);
  if (negativeDisplayMode.value === "cutoff") {
    drawQuadrants(ctx, margin, plotWidth, plotHeight, cutoffX, zeroY);
  }
  ctx.strokeStyle = "#d9ddd8";
  ctx.lineWidth = 1;
  ctx.strokeRect(margin.left, margin.top, plotWidth, plotHeight);

  ctx.font = "12px IBM Plex Sans, sans-serif";
  ctx.fillStyle = "#5e6d74";
  ctx.textBaseline = "middle";
  ctx.textAlign = "right";
  for (const tick of niceTicks(yMin, yMax)) {
    const y = yAt(tick);
    ctx.strokeStyle = "#eef1ed";
    ctx.beginPath();
    ctx.moveTo(margin.left, y);
    ctx.lineTo(margin.left + plotWidth, y);
    ctx.stroke();
    ctx.fillText(tick.toFixed(1), margin.left - 9, y);
  }

  ctx.textAlign = "center";
  ctx.textBaseline = "top";
  for (const tick of niceTicks(xMin, xMax)) {
    const x = xAt(tick);
    ctx.strokeStyle = "#eef1ed";
    ctx.beginPath();
    ctx.moveTo(x, margin.top);
    ctx.lineTo(x, margin.top + plotHeight);
    ctx.stroke();
    ctx.fillStyle = "#5e6d74";
    ctx.fillText(tick.toFixed(1), x, margin.top + plotHeight + 8);
  }

  if (yMin < 0 && yMax > 0) {
    ctx.strokeStyle = "#1f2a30";
    ctx.lineWidth = 1.2;
    ctx.beginPath();
    ctx.moveTo(margin.left, zeroY);
    ctx.lineTo(margin.left + plotWidth, zeroY);
    ctx.stroke();
  }

  if (negativeDisplayMode.value === "cutoff" && cutoffX >= margin.left && cutoffX <= margin.left + plotWidth) {
    ctx.strokeStyle = "#b07d00";
    ctx.lineWidth = 1.4;
    ctx.setLineDash([5, 4]);
    ctx.beginPath();
    ctx.moveTo(cutoffX, margin.top);
    ctx.lineTo(cutoffX, margin.top + plotHeight);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = "#7a5900";
    ctx.textAlign = cutoffX > width - 145 ? "right" : "left";
    ctx.textBaseline = "top";
    ctx.fillText(`Negative cutoff ${negativeEssentialityCutoff.toFixed(2)}`, cutoffX + (cutoffX > width - 145 ? -5 : 5), 7);
  }

  ctx.save();
  ctx.translate(15, margin.top + plotHeight / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.textAlign = "center";
  ctx.textBaseline = "middle";
  ctx.fillStyle = "#1f2a30";
  ctx.fillText(activeAnalysis.effect_metric === "hedges_g" ? "Relative essentiality (Hedges' g)" : `Differential dependency: ${yLabel()}`, 0, 0);
  ctx.restore();

  ctx.textAlign = "center";
  ctx.textBaseline = "bottom";
  ctx.fillStyle = "#5e6d74";
  ctx.fillText(xLabel(), margin.left + plotWidth / 2, height - 5);

  if (negativeDisplayMode.value === "color") {
    drawNegativeDependencyLegend(ctx, width);
  }

  const positions = [];
  ctx.globalAlpha = activeGene ? 0.16 : 0.58;
  data.forEach((d) => {
    const x = xAt(d.positive_average);
    const y = yAt(d.score);
    positions.push({ x, y, datum: d });
    ctx.fillStyle = negativeDisplayMode.value === "color" ? negativeDependencyColor(d.negative_average) : chart.color;
    ctx.beginPath();
    ctx.arc(x, y, 1.7, 0, Math.PI * 2);
    ctx.fill();
  });
  ctx.globalAlpha = 1;

  if (highlighted.length) {
    ctx.fillStyle = "#f2a83b";
    ctx.strokeStyle = "#1f2a30";
    ctx.lineWidth = 1.5;
    highlighted.forEach((d) => {
      const point = positions.find((item) => item.datum.gene === d.gene);
      if (!point) {
        return;
      }
      ctx.beginPath();
      ctx.arc(point.x, point.y, 4.5, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();
    });
  }

  chart.points = positions;
  chart.count.value = negativeDisplayMode.value === "cutoff"
    ? `${data.length.toLocaleString()} / ${allData.length.toLocaleString()} genes`
    : `${allData.length.toLocaleString()} genes`;
  chart.subtitle.textContent = `${activeAnalysis.positive_label} absolute dependency (X) and ${activeAnalysis.effect_metric === "hedges_g" ? "standardized positive-vs-negative effect (Y)" : `${yLabel()} (Y)`}${negativeDisplayMode.value === "color" ? "; dot color is negative-cohort dependency" : ""}`;
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
      <span class="small-line">${escapeHtml(activeAnalysis.source)}. ${activeAnalysis.effect_metric === "hedges_g" ? `Y-axis is Hedges' g standardized by pooled within-cohort variability; negative values indicate greater essentiality in ${escapeHtml(activeAnalysis.positive_label)}.` : `Y-axis is ${escapeHtml(yLabel())}.`} X-axis is absolute dependency in ${escapeHtml(activeAnalysis.positive_label)}.</span>
    </div>
    <div>
      <strong>Screened comparison sets</strong>
      <span class="small-line">CRISPR: ${crisprCounts.positive.length} ${escapeHtml(activeAnalysis.positive_label)} vs ${crisprCounts.negative_n} ${escapeHtml(activeAnalysis.negative_label)}; genes require ${crisprCounts.min_positive_n || 1}+ / ${crisprCounts.min_negative_n || 1}+ values. siRNA: ${rnaiCounts.positive.length} ${escapeHtml(activeAnalysis.positive_label)} vs ${rnaiCounts.negative_n} ${escapeHtml(activeAnalysis.negative_label)}; genes require ${rnaiCounts.min_positive_n || 1}+ / ${rnaiCounts.min_negative_n || 1}+ values.</span>
    </div>
    <details class="model-details">
      <summary>${positiveModels.length.toLocaleString()} positive-group DepMap models</summary>
      <ul>${modelRows || "<li>No positive models matched this analysis.</li>"}</ul>
    </details>
  `;
}

function formatPercent(value) {
  const percent = value * 100;
  return `${percent >= 10 ? percent.toFixed(0) : percent.toFixed(1)}%`;
}

function renderStratifierPrevalence() {
  const prevalence = activeAnalysis.stratifier_prevalence;
  if (!prevalence) {
    stratifierPrevalencePanel.innerHTML = "";
    return;
  }
  stratifierPrevalencePanel.innerHTML = `
    <div>
      <h2>Stratifier Prevalence</h2>
      <p class="small-line">${escapeHtml(activeAnalysis.positive_label)} among ${escapeHtml(prevalence.denominator)}.</p>
    </div>
    <div class="prevalence-stat">
      <strong>${formatPercent(prevalence.frequency)}</strong>
      <span>${prevalence.positive.toLocaleString()} of ${prevalence.total.toLocaleString()} ${escapeHtml(prevalence.denominator)}</span>
    </div>
  `;
}

function renderResults() {
  const rnai = new Map(visibleRows("rnai").map((d) => [d.gene, d]));
  const crispr = new Map(visibleRows("crispr").map((d) => [d.gene, d]));
  const genes = new Set([...rnai.keys(), ...crispr.keys()]);
  const rows = [...genes]
    .filter((gene) => (geneQuery ? gene.toLowerCase().includes(geneQuery.toLowerCase()) : false))
    .sort()
    .slice(0, 50)
    .map((gene) => {
      const r = rnai.get(gene);
      const c = crispr.get(gene);
      const isSelected = activeGene && gene.toLowerCase() === activeGene.toLowerCase();
      return `
        <tr>
          <td><button class="gene-result-button${isSelected ? " is-active" : ""}" type="button" data-gene="${escapeHtml(gene)}" aria-pressed="${isSelected ? "true" : "false"}">${escapeHtml(gene)}</button></td>
          <td>${c ? formatScore(c.positive_average) : ""}</td>
          <td>${c ? formatScore(c.score) + ` (#${c.rank})` : ""}</td>
          <td>${formatPair(c)}</td>
          <td>${r ? formatScore(r.positive_average) : ""}</td>
          <td>${r ? formatScore(r.score) + ` (#${r.rank})` : ""}</td>
          <td>${formatPair(r)}</td>
        </tr>
      `;
    });

  resultBody.innerHTML = rows.length
    ? rows.join("")
    : `<tr><td colspan="7" class="muted">${geneQuery && negativeDisplayMode.value === "cutoff" ? "No matched genes pass the active negative-cohort cutoff." : "Search for a gene, then select an exact symbol to highlight it."}</td></tr>`;
}

function renderAll() {
  renderMeta();
  renderStratifierPrevalence();
  drawChart("crispr");
  drawChart("rnai");
  renderResults();
}

function populateAnalyses() {
  const groups = new Map();
  summary.analyses.forEach((analysis) => {
    const group = analysis.category || "Additional analyses";
    if (!groups.has(group)) {
      groups.set(group, []);
    }
    groups.get(group).push(analysis);
  });
  analysisSelect.innerHTML = [...groups.entries()]
    .map(([group, analyses]) => `
      <optgroup label="${escapeHtml(group)}">
        ${analyses.map((analysis) => `<option value="${escapeHtml(analysis.id)}">${escapeHtml(analysis.label)}</option>`).join("")}
      </optgroup>
    `)
    .join("");
}

function populateGeneOptions() {
  const genes = new Set([
    ...activeAnalysis.datasets.crispr.map((row) => row.gene),
    ...activeAnalysis.datasets.rnai.map((row) => row.gene),
  ]);
  geneOptions.innerHTML = [...genes]
    .sort()
    .map((gene) => `<option value="${escapeHtml(gene)}"></option>`)
    .join("");
}

function selectExactGene(value) {
  const query = value.trim();
  geneQuery = query;
  if (!activeAnalysis || !activeAnalysis.datasets) {
    activeGene = "";
    return;
  }
  const genes = new Map(
    [
      ...activeAnalysis.datasets.crispr,
      ...activeAnalysis.datasets.rnai,
    ].map((row) => [row.gene.toLowerCase(), row.gene])
  );
  activeGene = genes.get(query.toLowerCase()) || "";
  if (activeGene) {
    geneQuery = activeGene;
    searchInput.value = activeGene;
  }
}

async function setAnalysis(analysisId) {
  const token = ++analysisLoadToken;
  activeAnalysis = summary.analyses.find((analysis) => analysis.id === analysisId) || summary.analyses[0];
  analysisSelect.value = activeAnalysis.id;
  hideTooltip();
  if (!activeAnalysis.datasets) {
    metaPanel.innerHTML = '<div><strong>Loading analysis</strong><span class="small-line">Fetching dependency data for the selected comparison...</span></div>';
    try {
      const response = await fetch(activeAnalysis.data_url);
      if (!response.ok) {
        throw new Error("Analysis data unavailable");
      }
      Object.assign(activeAnalysis, normalizeAnalysis(await response.json()));
    } catch (error) {
      if (token === analysisLoadToken) {
        metaPanel.innerHTML = '<div><strong>Analysis unavailable</strong><span class="small-line">Rebuild the dependency analysis data and try again.</span></div>';
      }
      return;
    }
  }
  if (token !== analysisLoadToken) {
    return;
  }
  populateGeneOptions();
  selectExactGene(searchInput.value);
  renderAll();
}

function nearestPoint(chart, x, y) {
  let nearest = null;
  let best = Infinity;
  for (const point of chart.points) {
    const distance = Math.hypot(point.x - x, point.y - y);
    if (distance < best) {
      best = distance;
      nearest = point;
    }
  }
  return best <= 8 ? nearest : null;
}

function showTooltip(event, key, point) {
  const d = point.datum;
  tooltip.innerHTML = `
    <strong>${escapeHtml(d.gene)}</strong>
    <span>${charts[key].title} positive-cohort selectivity rank #${d.rank.toLocaleString()}</span>
    <span>Positive absolute: ${formatScore(d.positive_average)}</span>
    <span>${activeAnalysis.effect_metric === "hedges_g" ? "Standardized effect (Hedges' g)" : "Differential"}: ${formatScore(d.score)}</span>
    <span>Raw mean difference: ${formatScore(d.raw_difference)}</span>
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

function updateCutoffControl() {
  negativeEssentialityCutoff = Number(cutoffInput.value);
  cutoffValue.value = negativeEssentialityCutoff.toFixed(2);
  const showingCutoff = negativeDisplayMode.value === "cutoff";
  cutoffSettings.classList.toggle("hidden", !showingCutoff);
  negativeColorHelp.classList.toggle("hidden", showingCutoff);
  cutoffInput.disabled = !showingCutoff;
  cutoffValue.classList.toggle("is-disabled", cutoffInput.disabled);
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
    updateCutoffControl();
    return setAnalysis(summary.analyses[0].id);
  })
  .catch(() => {
    metaPanel.innerHTML = '<div><strong>Data summary not found</strong><span class="small-line">Run scripts/build_hpv_dependency_data.py.</span></div>';
  });

analysisSelect.addEventListener("change", async (event) => {
  await setAnalysis(event.target.value);
});

searchInput.addEventListener("input", (event) => {
  selectExactGene(event.target.value);
  if (summary && activeAnalysis && activeAnalysis.datasets) {
    renderAll();
  }
});

resultBody.addEventListener("click", (event) => {
  const button = event.target.closest(".gene-result-button");
  if (!button) {
    return;
  }
  selectExactGene(button.dataset.gene || "");
  renderAll();
});

cutoffInput.addEventListener("input", () => {
  updateCutoffControl();
  hideTooltip();
  if (summary && activeAnalysis?.datasets && activeAnalysis?.included_models) {
    renderAll();
  }
});

negativeDisplayMode.addEventListener("change", () => {
  updateCutoffControl();
  hideTooltip();
  if (summary && activeAnalysis?.datasets && activeAnalysis?.included_models) {
    renderAll();
  }
});

window.addEventListener("resize", () => {
  if (summary && activeAnalysis?.datasets && activeAnalysis?.included_models) {
    renderAll();
  }
});
