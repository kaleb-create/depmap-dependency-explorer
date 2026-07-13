(async function () {
  const data = await fetch("/static/data/ddx3/ddx3_selectivity.json").then((response) => response.json());
  const viewerEl = document.getElementById("ddx3Viewer");
  const state = { model: "x", labels: false };
  const coreSelector = { chain: "A", resi: "132-607" };
  const colors = {
    "ATP clamp": "#d94f4f",
    "RNA rim": "#1c8f75",
    "Allosteric bridge": "#4f6fd9"
  };

  const viewer = $3Dmol.createViewer(viewerEl, { backgroundColor: "#fffdf7" });
  const ddx3x = await fetch(data.models.ddx3x.pdbPath).then((response) => response.text());
  const ddx3y = await fetch(data.models.ddx3y.pdbPath).then((response) => response.text());
  const xModel = viewer.addModel(ddx3x, "pdb");
  const yModel = viewer.addModel(ddx3y, "pdb");

  function residueOr(residues) {
    return residues.map((resi) => ({ resi }));
  }

  function activeModelSelector() {
    return state.model === "overlay" ? coreSelector : { ...coreSelector, model: state.model === "x" ? xModel.getID() : yModel.getID() };
  }

  function paintModel() {
    viewer.removeAllLabels();
    viewer.setStyle({}, {});
    yModel.setStyle(coreSelector, { cartoon: { color: "#c7d1d2", opacity: state.model === "overlay" ? 0.38 : 1 } });
    xModel.setStyle(coreSelector, { cartoon: { colorscheme: { prop: "b", gradient: "roygb", min: 40, max: 100 } } });

    if (state.model === "y") {
      yModel.setStyle(coreSelector, { cartoon: { color: "#6d7d86" } });
      xModel.setStyle(coreSelector, {});
    }
    if (state.model === "x") {
      yModel.setStyle(coreSelector, {});
    }

    data.motifGroups.forEach((group) => {
      viewer.setStyle(
        { ...activeModelSelector(), chain: "A", or: residueOr(group.residues) },
        { cartoon: { color: colors[group.name], opacity: 0.95 }, stick: { color: colors[group.name], radius: 0.14 } }
      );
    });

    data.hotspotSubstitutions.forEach((sub) => {
      viewer.addSphere({
        center: { x: sub.coord[0], y: sub.coord[1], z: sub.coord[2] },
        radius: sub.minMotifDistance <= 14 ? 1.35 : 0.9,
        color: sub.minMotifDistance <= 14 ? "#ffb000" : "#8a3ffc",
        alpha: 0.9
      });
      if (state.labels && sub.minMotifDistance <= 15) {
        viewer.addLabel(`${sub.to}${sub.x}`, {
          position: { x: sub.coord[0], y: sub.coord[1], z: sub.coord[2] },
          backgroundColor: "#1f2a30",
          fontColor: "#ffffff",
          fontSize: 11,
          showBackground: true
        });
      }
    });

    viewer.zoomTo(coreSelector);
    viewer.render();
  }

  document.querySelectorAll("[data-ddx3-model]").forEach((button) => {
    button.addEventListener("click", () => {
      state.model = button.getAttribute("data-ddx3-model");
      document.querySelectorAll("[data-ddx3-model]").forEach((modelButton) => {
        modelButton.classList.toggle("btn-ghost", modelButton !== button);
      });
      paintModel();
    });
  });

  document.getElementById("ddx3Labels").addEventListener("change", (event) => {
    state.labels = event.target.checked;
    paintModel();
  });

  document.getElementById("ddx3Identity").textContent = `${data.summary.fullIdentity}%`;
  document.getElementById("ddx3CoreSubs").textContent = data.summary.coreSubstitutions;
  document.getElementById("ddx3NearMotif").textContent = data.summary.nearMotifCoreSubstitutions;

  const readoutEl = document.getElementById("ddx3Readouts");
  data.readouts.forEach((readout) => {
    const card = document.createElement("article");
    card.className = "ddx3-readout";
    card.innerHTML = `
      <div class="ret-kpi-value">${readout.value}</div>
      <h3>${readout.label}</h3>
      <p>${readout.interpretation}</p>
    `;
    readoutEl.appendChild(card);
  });

  const tableBody = document.getElementById("ddx3Hotspots");
  data.hotspotSubstitutions.forEach((sub) => {
    const row = document.createElement("tr");
    const closest = Object.entries(sub.motifDistances).sort((a, b) => a[1] - b[1])[0];
    row.innerHTML = `
      <td><strong>${sub.from}${sub.y} -> ${sub.to}${sub.x}</strong></td>
      <td>${sub.domain}</td>
      <td>${closest[0]}</td>
      <td>${closest[1].toFixed(1)} A</td>
      <td>${sub.plddt.toFixed(1)}</td>
    `;
    tableBody.appendChild(row);
  });

  const sourceList = document.getElementById("ddx3SourceList");
  data.sources.forEach((source) => {
    const anchor = document.createElement("a");
    anchor.className = "ret-source";
    anchor.href = source.href;
    anchor.target = "_blank";
    anchor.rel = "noreferrer";
    anchor.textContent = source.label;
    sourceList.appendChild(anchor);
  });

  paintModel();
})();
