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
  const conceptsListContainer = document.getElementById("concepts-list");
  let selectedNavBlockId = null;
  let conceptList = [];
  let currentIndex = 0;
  let currentActiveNodeId = null;
  let _blockFocusTimer = null;
  let isGlobalNav = false; // true = navegar los 69, false = solo el bloque
  // Map to track concept list item DOM elements by node id
  const conceptItemByNodeId = new Map();

  // Importancia por conexiones (para escalar texto)
  const connectionCount = new Map();
  nodes.forEach((n) => connectionCount.set(n.id, 0));
  links.forEach((l) => {
    const s = typeof l.source === "object" ? l.source.id : l.source;
    const t = typeof l.target === "object" ? l.target.id : l.target;
    connectionCount.set(s, (connectionCount.get(s) || 0) + 1);
    connectionCount.set(t, (connectionCount.get(t) || 0) + 1);
  });
  const maxConns = Math.max(1, ...connectionCount.values());
  // Tracking de sprites para zoom dinámico
  const labelDataByNodeId = new Map();

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
        // Sincronizar el checkbox de "Seleccionar todo"
        const selectAllCb = document.getElementById("select-all-check");
        if (selectAllCb) selectAllCb.checked = activeBlocks.size === bloques.length;
      });
    });

    // Botón "Seleccionar todo"
    const selectAllCheck = document.getElementById("select-all-check");
    if (selectAllCheck) {
      selectAllCheck.addEventListener("change", () => {
        if (selectAllCheck.checked) {
          bloques.forEach((b) => {
            activeBlocks.add(b.id);
            const cb = checkboxByBlock.get(b.id);
            if (cb) { cb.checked = true; cb.closest(".filter-pill").classList.add("active"); }
          });
        } else {
          bloques.forEach((b) => {
            activeBlocks.delete(b.id);
            const cb = checkboxByBlock.get(b.id);
            if (cb) { cb.checked = false; cb.closest(".filter-pill").classList.remove("active"); }
          });
          selectedNavBlockId = null;
          conceptList = [];
          currentIndex = 0;
          currentActiveNodeId = null;
          updateNavUI();
        }
        updateFilteredGraphData();
      });
    }
  }

  // Construir lista de todos los conceptos organizada por bloques
  function buildConceptsList() {
    if (!conceptsListContainer) return;
    conceptsListContainer.innerHTML = "";
    conceptItemByNodeId.clear();
    let globalNum = 0;

    bloques.forEach((bloque) => {
      const bloqueNodes = nodes.filter((n) => n.bloque === bloque.id);
      if (!bloqueNodes.length) return;

      const group = document.createElement("div");
      group.className = "concept-block-group";

      // Header del bloque
      const header = document.createElement("div");
      header.className = "concept-block-header";
      const dot = document.createElement("span");
      dot.className = "dot";
      dot.style.backgroundColor = bloque.color;
      const nameSpan = document.createElement("span");
      nameSpan.textContent = bloque.nombre;
      const countSpan = document.createElement("span");
      countSpan.className = "count";
      countSpan.textContent = `(${bloqueNodes.length})`;
      header.appendChild(dot);
      header.appendChild(nameSpan);
      header.appendChild(countSpan);
      group.appendChild(header);

      // Items de conceptos
      bloqueNodes.forEach((node) => {
        globalNum++;
        const item = document.createElement("div");
        item.className = "concept-item";
        item.dataset.nodeId = node.id;
        item.dataset.blockId = bloque.id;

        const numSpan = document.createElement("span");
        numSpan.className = "concept-num";
        numSpan.textContent = globalNum + ".";

        const nameEl = document.createElement("span");
        nameEl.className = "concept-name";
        nameEl.textContent = node.nombre;
        nameEl.title = node.nombre;

        item.appendChild(numSpan);
        item.appendChild(nameEl);
        group.appendChild(item);

        conceptItemByNodeId.set(node.id, item);

        item.addEventListener("click", () => {
          // Navegación global: los 69 conceptos en orden
          isGlobalNav = true;
          conceptList = nodes.slice(); // todos los 69
          currentIndex = Math.max(0, conceptList.findIndex((c) => c.id === node.id));
          selectedNavBlockId = null;
          // Activar bloque del concepto si no está activo
          if (!activeBlocks.has(bloque.id)) {
            activeBlocks.add(bloque.id);
            const cb = checkboxByBlock.get(bloque.id);
            if (cb) { cb.checked = true; cb.closest(".filter-pill").classList.add("active"); }
            updateFilteredGraphData();
          }
          setTimeout(() => focusConcept(node.id), 200);
        });
      });

      conceptsListContainer.appendChild(group);
    });
  }

  // Actualizar highlight del item activo en la lista
  function updateActiveListItem(nodeId) {
    conceptItemByNodeId.forEach((el, id) => {
      el.classList.toggle("active", id === nodeId);
    });
    // Scroll el item activo a la vista
    const activeEl = conceptItemByNodeId.get(nodeId);
    if (activeEl) {
      activeEl.scrollIntoView({ behavior: "smooth", block: "nearest" });
    }
  }

  // Animar cámara hacia un nodo del grafo 3D
  // Vuelo 3D pero al llegar queda de frente (vista 2D centrada en TODA la pantalla)
  function flyToNode(nodeId) {
    const graphNode = Graph.graphData().nodes.find((n) => n.id === nodeId);
    if (!graphNode || graphNode.x === undefined) return;
    const distance = 120;

    // Compensar el panel izquierdo para centrar en el viewport completo
    const panel = document.getElementById("glass-panel");
    const panelW = panel ? panel.offsetWidth + 32 : 0; // 32px = margins
    const graphEl = document.getElementById("3d-graph");
    const graphW = graphEl ? graphEl.offsetWidth : window.innerWidth;
    const cam = Graph.camera();
    const fov = cam.fov * Math.PI / 180;
    const visibleW = 2 * distance * Math.tan(fov / 2);
    // Offset en unidades 3D: desplazar la vista a la derecha para que el nodo
    // quede al centro de la pantalla completa (no solo del canvas)
    const xOffset = (panelW / 2) * (visibleW / graphW);

    Graph.cameraPosition(
      { x: graphNode.x + xOffset, y: graphNode.y, z: graphNode.z + distance },
      { x: graphNode.x + xOffset, y: graphNode.y, z: graphNode.z },
      2500
    );
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
    .nodeRelSize(3)
    .linkColor(() => "rgba(200,220,255,0.45)")
    .linkOpacity(0.6)
    .linkWidth(0.4)
    .linkDirectionalArrowLength(2.5)
    .linkDirectionalArrowRelPos(1)
    .linkDirectionalParticles(0)
    .linkDirectionalParticleWidth(1.5)
    .linkDirectionalParticleColor(() => "#ffffff")
    .linkLabel((link) => link.relacion || "")
    // Estirar links y forzar distribución 3D real
    .d3Force("link", Graph.d3Force("link").distance(120))
    .d3Force("charge", Graph.d3Force("charge").strength(-200));

  // Configuración de verbos permanentes en aristas usando SpriteText
  if (typeof SpriteText !== "undefined") {
    // Texto en aristas con el verbo de relación, centrado en la línea
    Graph.linkThreeObject((link) => {
      const sprite = new SpriteText(link.relacion || "");
      sprite.color = "#f5f5ff"; // visible sobre fondo negro
      sprite.textHeight = 1.8;
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
  const infoConnections = document.getElementById("info-connections");
  const infoAnalogy = document.getElementById("info-analogy");

  function mostrarNodoEnPanel(node) {
    if (!node || !infoPanel) return;

    const bloqueNombre = bloqueNombreById.get(node.bloque) || "";
    infoTitle.textContent = node.nombre || "Concepto";
    infoBlock.textContent = bloqueNombre ? `Bloque: ${bloqueNombre}` : "";
    if (infoDesc) {
      infoDesc.textContent = node.descripcion || "";
    }

    // Conexiones: buscar links donde este nodo es source o target
    if (infoConnections) {
      const conns = [];
      links.forEach((l) => {
        const srcId = getId(l.source);
        const tgtId = getId(l.target);
        if (srcId === node.id) {
          const tgt = nodeById.get(tgtId);
          if (tgt) conns.push({ dir: "→", name: tgt.nombre, rel: l.relacion || "" });
        } else if (tgtId === node.id) {
          const src = nodeById.get(srcId);
          if (src) conns.push({ dir: "←", name: src.nombre, rel: l.relacion || "" });
        }
      });
      if (conns.length) {
        infoConnections.innerHTML = conns.map((c) =>
          `<div class="conn-item"><span class="conn-arrow">${c.dir}</span><span class="conn-target">${c.name}</span><span class="conn-relation">${c.rel}</span></div>`
        ).join("");
      } else {
        infoConnections.textContent = "Sin conexiones";
      }
    }

    // Analogía
    if (infoAnalogy) {
      const analogias = window.BIOINFO_ANALOGIAS;
      const txt = analogias && analogias[node.id] ? analogias[node.id] : "Analogia no disponible";
      infoAnalogy.textContent = txt;
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
    const conns = connectionCount.get(node.id) || 0;
    const importance = 0.7 + (conns / maxConns) * 0.6; // 0.7 a 1.3
    const isActive = node.id === currentActiveNodeId;
    const baseHeight = (isActive ? 3 : 2) * importance;

    const label = new SpriteText(node.nombre || "");
    label.color = "#ffd86b";
    label.textHeight = baseHeight;
    label.position.y = 7;
    group.add(label);

    // Guardar referencia para zoom dinámico
    labelDataByNodeId.set(node.id, {
      sprite: label,
      baseScaleX: label.scale.x,
      baseScaleY: label.scale.y
    });

    if (isActive) {
      const color = bloqueColorById.get(node.bloque) || "#ffffff";
      const geom = new THREE.SphereGeometry(6, 16, 16);
      const mat = new THREE.MeshBasicMaterial({ color, transparent: true, opacity: 0.25 });
      group.add(new THREE.Mesh(geom, mat));
    }
    return group;
  }

  Graph.nodeThreeObject(buildNodeObject).nodeThreeObjectExtend(true);

  // Escalar texto dinámicamente según zoom de la cámara
  let _zoomTimer = null;
  function adjustLabelsForZoom() {
    if (_zoomTimer) return;
    _zoomTimer = setTimeout(() => {
      _zoomTimer = null;
      const cam = Graph.camera();
      const dist = cam.position.length();
      // Lejos (~500+) → texto grande, Cerca (~120 flyTo) → texto chico
      const factor = Math.max(0.3, Math.min(1.8, dist / 450));
      labelDataByNodeId.forEach(({ sprite, baseScaleX, baseScaleY }) => {
        sprite.scale.set(baseScaleX * factor, baseScaleY * factor, 1);
      });
    }, 60);
  }

  // Conectar al control de cámara después de inicializar
  setTimeout(() => {
    try {
      Graph.controls().addEventListener("change", adjustLabelsForZoom);
    } catch (e) { /* controls no disponibles aún */ }
  }, 500);

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
    const concept = conceptList.length ? conceptList[currentIndex] : null;
    // En modo global: mostrar bloque del concepto actual; en modo bloque: el bloque seleccionado
    const blockId = isGlobalNav && concept ? concept.bloque : selectedNavBlockId;
    const bName = blockId ? bloqueNombreById.get(blockId) || "" : "";
    navBlockName.textContent = isGlobalNav ? (bName ? `${bName} · Global` : "Global") : bName;
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
    updateActiveListItem(conceptId);
    refreshObjects();
    flyToNode(conceptId);
  }

  function navigate(delta) {
    if (!conceptList.length) return;
    currentIndex = (currentIndex + delta + conceptList.length) % conceptList.length;
    const concept = conceptList[currentIndex];
    // En modo global, activar el bloque del siguiente concepto si no está visible
    if (isGlobalNav && !activeBlocks.has(concept.bloque)) {
      activeBlocks.add(concept.bloque);
      const cb = checkboxByBlock.get(concept.bloque);
      if (cb) { cb.checked = true; cb.closest(".filter-pill").classList.add("active"); }
      updateFilteredGraphData();
      setTimeout(() => focusConcept(concept.id), 200);
      return;
    }
    focusConcept(concept.id);
  }

  if (navPrevBtn) {
    navPrevBtn.addEventListener("click", () => navigate(-1));
  }
  if (navNextBtn) {
    navNextBtn.addEventListener("click", () => navigate(1));
  }

  // Click en nodo del grafo 3D: navegación solo por bloque
  Graph.onNodeClick((node) => {
    isGlobalNav = false;
    selectedNavBlockId = node.bloque;
    conceptList = nodes.filter((n) => n.bloque === node.bloque);
    currentIndex = Math.max(0, conceptList.findIndex((c) => c.id === node.id));
    focusConcept(node.id);
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
  buildConceptsList();
  updateFilteredGraphData();
  populateConceptsDatalist();
  // Sin selección inicial de bloque para no mostrar grafo hasta que el usuario elija
})();

