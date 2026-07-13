(async function () {
  const data = await fetch("/static/data/ret/ret_allostery.json").then((r) => r.json());

  const pocketColors = Object.fromEntries(data.pockets.map((p) => [p.slug, p.color]));
  const viewerEl = document.getElementById("retViewer");
  const viewer = $3Dmol.createViewer(viewerEl, { backgroundColor: "#fffdf7" });
  const pdbText = await fetch(data.alphafold.pdbPath).then((r) => r.text());
  viewer.addModel(pdbText, "pdb");

  const colorConfidence = () => {
    viewer.setStyle({}, {});
    viewer.setStyle({}, { cartoon: { colorscheme: { prop: "b", gradient: "roygb", min: 0, max: 100 } } });
  };

  const colorCartoon = () => {
    viewer.setStyle({}, {});
    viewer.setStyle({}, { cartoon: { color: "#6d7d86" } });
    data.pockets.forEach((p) => {
      const residueSpec = p.residues.map((resi) => ({ resi }));
      viewer.setStyle({ chain: "A", or: residueSpec }, { stick: { colorscheme: p.color, radius: 0.18 } });
    });
  };

  const colorSurface = () => {
    viewer.setStyle({}, {});
    viewer.setStyle({}, { cartoon: { color: "#d8dedc", opacity: 0.55 } });
    viewer.addSurface($3Dmol.SurfaceType.VDW, { color: "#9ad8d5", opacity: 0.32 }, {});
  };

  colorConfidence();
  data.pockets.forEach((p) => {
    const [x, y, z] = p.rank === 1 ? [29, 10, -15] : p.rank === 2 ? [25, 15, -31] : [22, 3, -25];
    viewer.addSphere({
      center: { x, y, z },
      radius: 2.4,
      color: p.color,
      alpha: 0.78,
      wireframe: false
    });
  });

  viewer.zoomTo({ chain: "A", resi: [705, 1013] });
  viewer.render();

  document.querySelectorAll("[data-style]").forEach((button) => {
    button.addEventListener("click", () => {
      document.querySelectorAll("[data-style]").forEach((b) => b.classList.add("btn-ghost"));
      button.classList.remove("btn-ghost");
      const style = button.getAttribute("data-style");
      if (style === "confidence") colorConfidence();
      if (style === "cartoon") colorCartoon();
      if (style === "surface") colorSurface();
      viewer.zoomTo({ chain: "A", resi: [705, 1013] });
      viewer.render();
    });
  });

  const switches = document.getElementById("retPocketSwitches");
  data.pockets.forEach((p) => {
    const label = document.createElement("label");
    label.className = "ret-switch";
    label.innerHTML = `<input type="checkbox" checked data-pocket="${p.slug}" /><span class="ret-swatch" style="background:${p.color}"></span>${p.name}`;
    switches.appendChild(label);
  });

  switches.addEventListener("change", (event) => {
    if (!(event.target instanceof HTMLInputElement)) return;
    viewer.removeAllShapes();
    data.pockets
      .filter((p) => switches.querySelector(`input[data-pocket="${p.slug}"]`)?.checked)
      .forEach((p) => {
        const [x, y, z] = p.rank === 1 ? [29, 10, -15] : p.rank === 2 ? [25, 15, -31] : [22, 3, -25];
        viewer.addSphere({ center: { x, y, z }, radius: 2.4, color: p.color, alpha: 0.78 });
      });
    viewer.render();
  });

  document.getElementById("retMeanPlddt").textContent = data.alphafold.meanPLDDT.toFixed(1);
  document.getElementById("retAfSummary").textContent = data.alphafold.summary;
  document.getElementById("retAfLink").href = data.alphafold.sourceHref;

  const ranksEl = document.getElementById("retPocketRanks");
  data.pockets
    .slice()
    .sort((a, b) => b.score - a.score)
    .forEach((p) => {
      const card = document.createElement("article");
      card.className = "ret-pocket-card";
      card.innerHTML = `
        <div class="ret-pocket-row">
          <div>
            <div class="eyebrow">Rank ${p.rank}</div>
            <h3>${p.name}</h3>
          </div>
          <div class="ret-score" style="--score:${p.score}; --pocket:${p.color}">
            <span>${p.score.toFixed(1)}</span>
          </div>
        </div>
        <p>${p.description}</p>
        <div class="ret-chip-row">${p.residues.map((r) => `<span class="chip">R${r}</span>`).join("")}</div>
      `;
      ranksEl.appendChild(card);
    });

  const ligandGrid = document.getElementById("retLigandGrid");
  data.ligands.forEach((ligand) => {
    const card = document.createElement("article");
    card.className = "ret-ligand-card";
    card.innerHTML = `
      <img src="${ligand.imagePath}" alt="${ligand.name} structure" />
      <div class="ret-ligand-body">
        <div class="ret-ligand-top">
          <div>
            <h3>${ligand.name}</h3>
            <p>${ligand.aliases.join(" · ")}</p>
          </div>
          <span class="chip">${ligand.selectivity}</span>
        </div>
        <p>${ligand.notes}</p>
        <div class="ret-chip-row">${ligand.primaryPockets.map((slug) => {
          const pocket = data.pockets.find((p) => p.slug === slug);
          return `<span class="chip" style="border-color:${pocket.color}; color:${pocket.color}">${pocket.name}</span>`;
        }).join("")}</div>
        <dl class="ret-props">
          <div><dt>MW</dt><dd>${ligand.molecularWeight}</dd></div>
          <div><dt>xLogP</dt><dd>${ligand.xlogp}</dd></div>
          <div><dt>TPSA</dt><dd>${ligand.tpsa}</dd></div>
          <div><dt>PDB</dt><dd>${ligand.evidencePdb.length ? ligand.evidencePdb.join(", ") : "Comparator"}</dd></div>
        </dl>
        <a class="text-link" href="${ligand.sourceHref}" target="_blank" rel="noreferrer">PubChem record</a>
      </div>
    `;
    ligandGrid.appendChild(card);
  });

  const sourceList = document.getElementById("retSourceList");
  data.sources.forEach((src) => {
    const a = document.createElement("a");
    a.className = "ret-source";
    a.href = src.href;
    a.target = "_blank";
    a.rel = "noreferrer";
    a.textContent = src.label;
    sourceList.appendChild(a);
  });
})();
