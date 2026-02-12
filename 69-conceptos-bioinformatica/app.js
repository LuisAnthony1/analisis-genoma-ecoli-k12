(function () {
  const rawData = window.BIOINFO_GRAPH;
  if (!rawData) {
    console.error("BIOINFO_GRAPH no está definido. Asegúrate de cargar data.js antes de app.js.");
    return;
  }

  const { bloques, nodes, links } = rawData;

  const filtersContainer = document.getElementById("filters-container");
  const infoPanel = document.getElementById("info-panel");
  const infoTitle = document.getElementById("info-title");
  const infoBlock = document.getElementById("info-block");
  const infoDesc = document.getElementById("info-desc");
  const infoCloseBtn = document.getElementById("info-close");
  const searchInput = document.getElementById("concept-search");
  const conceptsDatalist = document.getElementById("concepts-datalist");
  const navPanel = document.getElementById("nav-panel");
  const navBlockName = document.getElementById("nav-block-name");
  const navConceptName = document.getElementById("nav-concept-name");
  const navPrevBtn = document.getElementById("nav-prev");
  const navNextBtn = document.getElementById("nav-next");
  const navCounter = document.getElementById("nav-counter");

  // Estado de filtros y flujo
  // Iniciamos con NINGÚN bloque seleccionado para evitar mostrar el grafo por defecto
  const activeBlocks = new Set();

  // Mapas de ayuda
  const bloqueColorById = new Map(bloques.map((b) => [b.id, b.color]));
  const bloqueNombreById = new Map(bloques.map((b) => [b.id, b.nombre]));
  // Mapa de nodos por ID para acceso rápido
  const nodeById = new Map(nodes.map((n) => [n.id, n]));
  const checkboxByBlock = new Map();
  let selectedNavBlockId = null;
  let conceptList = [];
  let currentIndex = 0;
  let currentActiveNodeId = null;
  let _blockFocusTimer = null;

  // Helper: extraer ID de source/target (puede ser string o objeto)
  function getId(nodeOrId) {
    return typeof nodeOrId === "object" && nodeOrId !== null ? nodeOrId.id : nodeOrId;
  }

  // Construir UI de filtros
  function buildBlockFilters() {
    filtersContainer.innerHTML = "";
    bloques.forEach((bloque) => {
      const pill = document.createElement("label");
      pill.className = "filter-pill active";

      const checkbox = document.createElement("input");
      checkbox.type = "checkbox";
      checkbox.checked = false;
      checkbox.dataset.blockId = bloque.id;
      checkboxByBlock.set(bloque.id, checkbox);

      const dot = document.createElement("span");
      dot.className = "dot";
      dot.style.backgroundColor = bloque.color;

      const labelSpan = document.createElement("span");
      labelSpan.className = "label";
      labelSpan.textContent = bloque.nombre;
      labelSpan.title = bloque.nombre;

      pill.appendChild(checkbox);
      pill.appendChild(dot);
      pill.appendChild(labelSpan);
      filtersContainer.appendChild(pill);

      checkbox.addEventListener("change", () => {
        if (checkbox.checked) {
          activeBlocks.add(bloque.id);
          pill.classList.add("active");
          setSelectedBlock(bloque.id, true);
        } else {
          activeBlocks.delete(bloque.id);
          pill.classList.remove("active");
          if (selectedNavBlockId === bloque.id) {
            selectedNavBlockId = null;
            conceptList = [];
            currentIndex = 0;
            currentActiveNodeId = null;
            updateNavUI();
          }
        }
        updateFilteredGraphData();
      });
    });
  }

  // Construir datos filtrados según bloques activos
  function getFilteredData() {
    // Filtrar nodos: solo los que pertenecen a bloques activos
    const filteredNodes = nodes.filter((n) => activeBlocks.has(n.bloque));
    const allowedIds = new Set(filteredNodes.map((n) => n.id));

    // Filtrar links: solo los que conectan dos nodos visibles
    // Maneja tanto strings (IDs) como objetos (nodos ya procesados por la librería)
    const filteredLinks = links.filter((l) => {
      const sourceId = getId(l.source);
      const targetId = getId(l.target);
      // Un link es visible solo si AMBOS nodos están en la lista de permitidos
      return allowedIds.has(sourceId) && allowedIds.has(targetId);
    });

    return {
      nodes: filteredNodes,
      links: filteredLinks
    };
  }

  // Crear instancia de ForceGraph3D (usa THREE global del HTML)
  const Graph = ForceGraph3D()(document.getElementById("3d-graph"));

  Graph
    .backgroundColor("rgba(0,0,0,0)") // se ve el fondo "space" del body
    // Tooltip al pasar el mouse: solo nombre del concepto
    .nodeLabel((node) => node.nombre || "")
    .nodeColor((node) => bloqueColorById.get(node.bloque) || "#ffffff")
    .nodeOpacity(0.95)
    .nodeRelSize(6)
    .linkColor(() => "rgba(200,220,255,0.45)")
    .linkOpacity(0.6)
    .linkWidth(0.5)
    .linkDirectionalArrowLength(3)
    .linkDirectionalArrowRelPos(1)
    .linkDirectionalParticles(0)
    .linkDirectionalParticleWidth(1.5)
    .linkDirectionalParticleColor(() => "#ffffff")
    .linkLabel((link) => link.relacion || "");

  // Configuración de verbos permanentes en aristas usando SpriteText
  if (typeof SpriteText !== "undefined") {
    // Texto en aristas con el verbo de relación, centrado en la línea
    Graph.linkThreeObject((link) => {
      const sprite = new SpriteText(link.relacion || "");
      sprite.color = "#f5f5ff"; // visible sobre fondo negro
      sprite.textHeight = 2.8;
      return sprite;
    }).linkThreeObjectExtend(true);

    Graph.linkPositionUpdate((sprite, { start, end }) => {
      // Colocar el SpriteText en el punto medio de la arista
      const middlePos = {
        x: (start.x + end.x) / 2,
        y: (start.y + end.y) / 2,
        z: (start.z + end.z) / 2
      };
      Object.assign(sprite.position, middlePos);
    });
  }

  // Mostrar detalles del nodo en el panel lateral
  function mostrarNodoEnPanel(node) {
    if (!node || !infoPanel) return;

    const bloqueNombre = bloqueNombreById.get(node.bloque) || "";
    infoTitle.textContent = node.nombre || "Concepto";
    infoBlock.textContent = bloqueNombre ? `Bloque: ${bloqueNombre}` : "";
    if (infoDesc) {
      infoDesc.textContent = node.descripcion || "";
    }

    infoPanel.classList.remove("hidden");
    infoPanel.classList.add("visible");
  }

  if (infoCloseBtn && infoPanel) {
    infoCloseBtn.addEventListener("click", () => {
      infoPanel.classList.remove("visible");
      infoPanel.classList.add("hidden");
    });
  }

  function buildNodeObject(node) {
    const group = new THREE.Group();
    const label = new SpriteText(node.nombre || "");
    label.color = "#ffd86b";
    label.textHeight = node.id === currentActiveNodeId ? 4.5 : 3.2;
    label.position.y = 10;
    group.add(label);
    if (node.id === currentActiveNodeId) {
      const color = bloqueColorById.get(node.bloque) || "#ffffff";
      const geom = new THREE.SphereGeometry(9, 16, 16);
      const mat = new THREE.MeshBasicMaterial({ color, transparent: true, opacity: 0.25 });
      group.add(new THREE.Mesh(geom, mat));
    }
    return group;
  }

  Graph.nodeThreeObject(buildNodeObject).nodeThreeObjectExtend(true);

  function refreshObjects() {
    try {
      Graph.refresh();
    } catch (e) {
      try {
        Graph.graphData(Graph.graphData());
      } catch {}
    }
  }

  function updateNavUI() {
    if (!navPanel) return;
    const bName = selectedNavBlockId ? bloqueNombreById.get(selectedNavBlockId) || "" : "";
    navBlockName.textContent = bName;
    const concept = conceptList.length ? conceptList[currentIndex] : null;
    navConceptName.textContent = concept ? (concept.nombre || "") : "";
    navCounter.textContent = conceptList.length ? `${currentIndex + 1} / ${conceptList.length}` : "0 / 0";
  }

  function setSelectedBlock(blockId, fromCheckbox) {
    selectedNavBlockId = blockId;
    conceptList = nodes.filter((n) => n.bloque === blockId).slice().sort((a, b) => (a.nombre || "").localeCompare(b.nombre || ""));
    currentIndex = 0;
    updateNavUI();
    if (!activeBlocks.has(blockId)) {
      activeBlocks.add(blockId);
      const cb = checkboxByBlock.get(blockId);
      if (cb) cb.checked = true;
      updateFilteredGraphData();
    }
    if (!fromCheckbox) {
      updateFilteredGraphData();
    }
    if (_blockFocusTimer) clearTimeout(_blockFocusTimer);
    _blockFocusTimer = setTimeout(() => {
      const first = conceptList[0];
      if (first) {
        focusConcept(first.id);
      }
    }, 120);
  }

  function focusConcept(conceptId) {
    const node = nodeById.get(conceptId);
    if (!node) return;
    currentActiveNodeId = conceptId;
    mostrarNodoEnPanel(node);
    updateNavUI();
    refreshObjects();
  }

  function navigate(delta) {
    if (!conceptList.length) return;
    currentIndex = (currentIndex + delta + conceptList.length) % conceptList.length;
    const id = conceptList[currentIndex].id;
    focusConcept(id);
  }

  if (navPrevBtn) {
    navPrevBtn.addEventListener("click", () => navigate(-1));
  }
  if (navNextBtn) {
    navNextBtn.addEventListener("click", () => navigate(1));
  }

  // Click en nodo sin mover la cámara
  Graph.onNodeClick((node) => {
    mostrarNodoEnPanel(node);
    selectedNavBlockId = node.bloque;
    conceptList = nodes.filter((n) => n.bloque === node.bloque).slice().sort((a, b) => (a.nombre || "").localeCompare(b.nombre || ""));
    currentIndex = Math.max(0, conceptList.findIndex((c) => c.id === node.id));
    currentActiveNodeId = node.id;
    updateNavUI();
    refreshObjects();
  });

  // Alternar flujo (resalta rutas "dogma central" y similares con partículas)
  function setFlujoActivo(active) {
    flujoActivo = active;
    if (flujoActivo) {
      Graph.linkDirectionalParticles(2);
      Graph.linkDirectionalParticleSpeed(0.006);
      Graph.linkColor((link) => {
        // enlaces que recorren el flujo principal o técnicas asociadas
        const flujoClave = new Set([
          "adn",
          "arn",
          "arn_mensajero",
          "traduccion",
          "proteina",
          "dogma_central",
          "pcr",
          "clonacion_molecular",
          "bioinformatica"
        ]);

        const srcId = typeof link.source === "object" ? link.source.id : link.source;
        const tgtId = typeof link.target === "object" ? link.target.id : link.target;

        const isFlujo = flujoClave.has(srcId) || flujoClave.has(tgtId);

        return isFlujo ? "rgba(255,230,128,0.95)" : "rgba(170,190,240,0.25)";
      });
    } else {
      Graph.linkDirectionalParticles(0);
      Graph.linkColor(() => "rgba(200,220,255,0.45)");
    }
  }

  function populateConceptsDatalist() {
    if (!conceptsDatalist) return;
    conceptsDatalist.innerHTML = "";
    nodes.forEach((n) => {
      const opt = document.createElement("option");
      opt.value = n.nombre || "";
      conceptsDatalist.appendChild(opt);
    });
  }

  function findAndSelectByName(name) {
    if (!name) return;
    const value = name.toLowerCase();
    let found = nodes.find((n) => (n.nombre || "").toLowerCase() === value);
    if (!found) {
      found = nodes.find((n) => (n.nombre || "").toLowerCase().includes(value));
    }
    if (!found) return;
    setSelectedBlock(found.bloque, false);
    // Cancelar el auto-focus al primer concepto del bloque
    if (_blockFocusTimer) { clearTimeout(_blockFocusTimer); _blockFocusTimer = null; }
    currentIndex = Math.max(0, conceptList.findIndex((c) => c.id === found.id));
    focusConcept(found.id);
  }

  if (searchInput) {
    searchInput.addEventListener("change", () => findAndSelectByName(searchInput.value.trim()));
    searchInput.addEventListener("keydown", (e) => {
      if (e.key === "Enter") {
        findAndSelectByName(searchInput.value.trim());
      }
    });
  }

  // Actualizar datos según filtros
  function updateFilteredGraphData() {
    const data = getFilteredData();
    Graph.graphData(data);
    setTimeout(() => {
      if (data.nodes.length > 0) {
        try {
          Graph.d3ReheatSimulation();
        } catch (e) {
          console.warn("No se pudo re-calentar la simulación:", e);
        }
      }
    }, 100);
  }

  // Inicializar
  buildBlockFilters();
  updateFilteredGraphData();
  populateConceptsDatalist();
  // Sin selección inicial de bloque para no mostrar grafo hasta que el usuario elija
})();

