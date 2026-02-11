/**
 * GenomeHub - Core Application
 *
 * Maneja sidebar con magic curve, navegacion, submenus,
 * estado global, temas y utilidades comunes.
 */

// =============================================================================
// ESTADO GLOBAL
// =============================================================================

const AppState = {
    currentSection: 'search',
    currentTab: null,
    selectedGenomes: new Set(),
    libraryGenomes: [],
    analysisGenome: null,
    resultsGenome: null,
    theme: localStorage.getItem('theme') || 'dark',
    apiBase: ''
};

// =============================================================================
// SIDEBAR - Toggle, Indicador, Submenus
// =============================================================================

function toggleSidebar() {
    document.getElementById('sidebar').classList.toggle('open');
    // Re-posicionar indicador despues de la transicion
    setTimeout(() => {
        const active = document.querySelector('.nav-item.active');
        if (active) moveIndicator(active);
    }, 420);
}

function moveIndicator(element) {
    const indicator = document.getElementById('nav-indicator');
    if (!indicator || !element) return;
    indicator.style.top = element.offsetTop + 'px';
}

function toggleSubmenu(navItem) {
    const wasOpen = navItem.classList.contains('submenu-open');

    // Cerrar todos los submenus
    document.querySelectorAll('.nav-item.submenu-open').forEach(item => {
        item.classList.remove('submenu-open');
    });

    // Abrir/cerrar este
    if (!wasOpen) {
        navItem.classList.add('submenu-open');
        // Si el sidebar esta cerrado, abrirlo para ver submenu
        const sidebar = document.getElementById('sidebar');
        if (!sidebar.classList.contains('open')) {
            sidebar.classList.add('open');
        }
    }
}

// =============================================================================
// NAVEGACION
// =============================================================================

function showSection(sectionName) {
    console.log(`[APP] Navegando a: ${sectionName}`);

    // Ocultar todas las secciones
    document.querySelectorAll('.section-panel').forEach(s => {
        s.classList.remove('active');
        s.style.display = 'none';
    });

    // Mostrar seccion
    const target = document.getElementById(`section-${sectionName}`);
    if (target) {
        target.classList.add('active');
        target.style.display = 'block';
    }

    AppState.currentSection = sectionName;

    // Sincronizar sidebar: buscar nav-item que corresponda a esta seccion y activarlo
    // 'results' se muestra bajo el nav-item 'analysis'
    const sectionToNav = { results: 'analysis', 'evo-results': 'evolucion' };
    const navSection = sectionToNav[sectionName] || sectionName;
    const navItem = document.querySelector(`.nav-item[data-section="${navSection}"]`);
    if (navItem && !navItem.classList.contains('active')) {
        document.querySelectorAll('.nav-item').forEach(i => i.classList.remove('active'));
        navItem.classList.add('active');
        moveIndicator(navItem);
    }

    // Cargar datos segun seccion
    if (sectionName === 'library') {
        loadLibrary();
    } else if (sectionName === 'analysis') {
        AnalysisRunner.loadGenomeSelector();
    } else if (sectionName === 'compare') {
        loadCompareSelectors();
    } else if (sectionName === 'results') {
        loadResultsGenomeSelector();
    } else if (sectionName === 'visor3d') {
        loadVisor3DGenomeSelector();
    } else if (sectionName === 'informe') {
        loadInformeGenomeSelector();
    } else if (sectionName === 'chat') {
        loadChatGenomeSelector();
    } else if (sectionName === 'evolucion') {
        loadEvolucionGenomeSelector();
    }
}

/**
 * Navegar a un nav-item del sidebar.
 * Para items sin submenu: activa la seccion.
 * Para items con submenu: abre el submenu y muestra la seccion principal.
 */
function navigateToItem(navItem) {
    const section = navItem.getAttribute('data-section');
    const hasSubmenu = navItem.classList.contains('has-submenu');

    // Actualizar activo
    document.querySelectorAll('.nav-item').forEach(i => i.classList.remove('active'));
    navItem.classList.add('active');
    moveIndicator(navItem);

    if (hasSubmenu) {
        toggleSubmenu(navItem);

        // Mostrar la seccion principal del grupo
        if (section === 'analysis') {
            showSection('analysis');
        } else if (section === 'compare') {
            showSection('compare');
        } else if (section === 'visor3d') {
            showSection('visor3d');
        } else if (section === 'evolucion') {
            showSection('evolucion');
        }
    } else {
        // Cerrar submenus abiertos
        document.querySelectorAll('.nav-item.submenu-open').forEach(i => {
            i.classList.remove('submenu-open');
        });
        showSection(section);
    }

    // Limpiar submenu-items activos
    document.querySelectorAll('.submenu-item.active').forEach(si => si.classList.remove('active'));
    AppState.currentTab = null;
}

/**
 * Click en un sub-item del sidebar.
 * Carga el dashboard o vista correspondiente.
 */
function navigateToSubItem(subItem) {
    const tab = subItem.getAttribute('data-tab');
    const parentItem = subItem.closest('.nav-item');

    // Marcar sub-item activo
    document.querySelectorAll('.submenu-item.active').forEach(si => si.classList.remove('active'));
    subItem.classList.add('active');

    // Asegurarse que el padre este activo
    document.querySelectorAll('.nav-item').forEach(i => i.classList.remove('active'));
    parentItem.classList.add('active');
    moveIndicator(parentItem);

    AppState.currentTab = tab;

    const section = parentItem.getAttribute('data-section');

    if (section === 'analysis') {
        // Sub-items de analisis → cargar dashboard
        showSection('results');
        loadResultsGenomeSelector(tab);
    } else if (section === 'compare') {
        if (tab === 'comp-config') {
            showSection('compare');
        } else if (tab === 'comp-resultados') {
            showSection('compare');
            // El dashboard de comparacion ya se muestra en section-compare
        }
    } else if (section === 'visor3d') {
        showSection('visor3d');
        loadVisor3DTab(tab);
    } else if (section === 'evolucion') {
        // Sub-items de evolucion → cargar dashboard de resultados
        showSection('evo-results');
        loadEvolucionTab(tab);
    }
}

// =============================================================================
// CARGADORES DE SELECTORES DE GENOMA (nuevas secciones)
// =============================================================================

async function loadVisor3DGenomeSelector() {
    const select = document.getElementById('visor3d-genome');
    if (!select) return;

    try {
        const resp = await fetch('/api/results_genomes');
        const data = await resp.json();
        if (data.success && data.genomes && data.genomes.length > 0) {
            select.innerHTML = data.genomes.map(g =>
                `<option value="${g.basename}">${g.label}</option>`
            ).join('');
        } else {
            select.innerHTML = '<option value="">No hay resultados</option>';
        }
    } catch {
        select.innerHTML = '<option value="">Error</option>';
    }
}

async function loadInformeGenomeSelector() {
    const select = document.getElementById('informe-genome');
    if (!select) return;

    try {
        const resp = await fetch('/api/results_genomes');
        const data = await resp.json();
        if (data.success && data.genomes && data.genomes.length > 0) {
            select.innerHTML = data.genomes.map(g =>
                `<option value="${g.basename}">${g.label}</option>`
            ).join('');
            const preselect = (typeof currentResultsGenome !== 'undefined' && currentResultsGenome) || '';
            if (preselect) select.value = preselect;
        } else {
            select.innerHTML = '<option value="">No hay resultados</option>';
        }
    } catch {
        select.innerHTML = '<option value="">Error</option>';
    }
}

async function loadChatGenomeSelector() {
    const select = document.getElementById('chat-genome');
    if (!select) return;

    try {
        const resp = await fetch('/api/results_genomes');
        const data = await resp.json();
        if (data.success && data.genomes && data.genomes.length > 0) {
            select.innerHTML = '<option value="">Sin genoma</option>' +
                data.genomes.map(g =>
                    `<option value="${g.basename}">${g.label}</option>`
                ).join('');
            const preselect = (typeof currentResultsGenome !== 'undefined' && currentResultsGenome) || '';
            if (preselect) select.value = preselect;
        }
    } catch {
        // silencioso
    }
}

function loadVisor3DTab(tab) {
    const container = document.getElementById('visor3d-container');
    const genome = document.getElementById('visor3d-genome')?.value;

    if (!genome) {
        container.innerHTML = `
            <div class="text-center py-12 text-secondary">
                <div class="text-6xl mb-3">🔬</div>
                <p>Selecciona un genoma con analisis de proteinas</p>
            </div>
        `;
        return;
    }

    // Cargar datos de proteinas y renderizar segun tab
    container.innerHTML = `
        <div class="flex items-center justify-center py-12">
            <div class="w-8 h-8 border-4 border-emerald-500 border-t-transparent rounded-full animate-spin mr-4"></div>
            <span class="text-secondary">Cargando datos de proteinas...</span>
        </div>
    `;

    DashboardRenderer.fetchResultData(genome, `analisis_proteinas_${genome}.json`).then(data => {
        if (!data) {
            container.innerHTML = `
                <div class="text-center py-12 text-secondary">
                    <div class="text-6xl mb-3">📭</div>
                    <p>No hay datos de proteinas para este genoma</p>
                    <p class="text-sm mt-2">Ejecuta el analisis de proteinas primero</p>
                </div>
            `;
            return;
        }

        if (tab === '3d-primaria') {
            renderVisor3DPrimaria(data, container);
        } else if (tab === '3d-secundaria') {
            renderVisor3DSecundaria(data, container);
        } else if (tab === '3d-terciaria') {
            renderVisor3DTerciaria(data, container);
        } else if (tab === '3d-cuaternaria') {
            renderVisor3DCuaternaria(data, container);
        }
    });
}

// Renders basicos para el visor 3D (se pueden expandir despues)
function renderVisor3DPrimaria(data, container) {
    const primaria = data.estructura_primaria || {};
    const stats = primaria.estadisticas_generales || {};
    const comp = primaria.composicion_aminoacidos || {};
    const top10 = primaria.top_10_mas_grandes || [];
    const mutaciones = primaria.mutaciones_patogenicas || {};
    const mutList = mutaciones.genes_analizados || mutaciones.mutaciones || [];

    // Calcular stats adicionales de composicion
    const compEntries = Object.entries(comp).sort((a, b) => b[1].porcentaje - a[1].porcentaje);
    const hidrof = compEntries.filter(([k]) => ['Ala', 'Val', 'Ile', 'Leu', 'Met', 'Phe', 'Trp', 'Pro', 'A', 'V', 'I', 'L', 'M', 'F', 'W', 'P'].includes(k));
    const hidrofPct = hidrof.reduce((s, [, v]) => s + (v.porcentaje || 0), 0);

    container.innerHTML = `
        <h3 class="text-xl font-bold text-primary mb-4">Estructura Primaria - Secuencia de Aminoacidos</h3>
        <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
            ${DashboardRenderer.statsCard('Total Proteinas', DashboardRenderer.fmt(stats.total_analizadas || stats.total_proteinas || 0), 'proteinas codificadas')}
            ${DashboardRenderer.statsCard('Longitud Promedio', DashboardRenderer.fmt(Math.round(stats.longitud_promedio_aa || 0)) + ' aa', 'aminoacidos', 'cyan')}
            ${DashboardRenderer.statsCard('Peso Molecular', DashboardRenderer.fmt(Math.round((stats.peso_molecular_promedio_da || 0) / 1000)) + ' kDa', 'promedio', 'violet')}
            ${DashboardRenderer.statsCard('Punto Isoelectrico', (stats.pi_promedio || 0).toFixed(1), 'pI promedio', 'amber')}
        </div>

        <!-- Composicion de Aminoacidos -->
        ${compEntries.length > 0 ? `
        <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
            <h4 class="text-sm font-semibold text-primary mb-3">Composicion de Aminoacidos del Proteoma</h4>
            <div class="grid grid-cols-2 gap-4 mb-4">
                ${DashboardRenderer.statsCard('Hidrofobicos', hidrofPct.toFixed(1) + '%', 'Ala, Val, Ile, Leu, Met, Phe, Trp, Pro', 'amber')}
                ${DashboardRenderer.statsCard('Hidrofilicos', (100 - hidrofPct).toFixed(1) + '%', 'Cargados + polares + otros', 'cyan')}
            </div>
            <canvas id="chart-3d-aa-comp" height="280"></canvas>
        </div>` : ''}

        <!-- Mutaciones Patogenicas -->
        ${mutList.length > 0 ? `
        <div class="bg-card rounded-xl p-5 border border-red-200 mb-6">
            <h4 class="text-sm font-semibold text-red-600 mb-3">Mutaciones Patogenicas - Genes de Resistencia a Antibioticos</h4>
            <p class="text-xs text-secondary mb-3">Genes asociados a resistencia antibiotica y virulencia detectados en el genoma.</p>
            <div class="overflow-x-auto max-h-[350px] overflow-y-auto">
                <table class="w-full text-xs">
                    <thead class="bg-red-50 dark:bg-red-900/20 sticky top-0">
                        <tr>
                            <th class="px-2 py-2 text-left text-secondary">Gen</th>
                            <th class="px-2 py-2 text-left text-secondary">Producto</th>
                            <th class="px-2 py-2 text-left text-secondary">Categoria</th>
                            <th class="px-2 py-2 text-right text-secondary">Long. (aa)</th>
                            <th class="px-2 py-2 text-right text-secondary">MW (kDa)</th>
                        </tr>
                    </thead>
                    <tbody>
                        ${mutList.map(m => `
                            <tr class="border-t border-slate-100 dark:border-slate-800">
                                <td class="px-2 py-2 font-mono font-bold text-red-600">${m.nombre_gen || m.locus_tag || 'N/A'}</td>
                                <td class="px-2 py-2 text-secondary">${(m.producto || '').substring(0, 50)}</td>
                                <td class="px-2 py-2"><span class="px-1.5 py-0.5 bg-red-100 dark:bg-red-900/30 text-red-700 rounded text-[10px]">${m.categoria || m.tipo || 'Resistencia'}</span></td>
                                <td class="px-2 py-2 text-right text-primary">${DashboardRenderer.fmt(m.longitud_aa || 0)}</td>
                                <td class="px-2 py-2 text-right text-primary">${m.peso_molecular_kda ? m.peso_molecular_kda.toFixed(1) : ((m.peso_molecular_da || 0) / 1000).toFixed(1)}</td>
                            </tr>
                        `).join('')}
                    </tbody>
                </table>
            </div>
        </div>` : ''}

        <!-- Top 10 Proteinas -->
        <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
            <h4 class="text-sm font-semibold text-primary mb-3">Top 10 Proteinas mas Grandes</h4>
            <div class="overflow-x-auto">
                <table class="w-full text-xs">
                    <thead><tr>
                        <th class="text-left px-2 py-2 border-b border-slate-200 text-secondary">#</th>
                        <th class="text-left px-2 py-2 border-b border-slate-200 text-secondary">Gen</th>
                        <th class="text-left px-2 py-2 border-b border-slate-200 text-secondary">Producto</th>
                        <th class="text-right px-2 py-2 border-b border-slate-200 text-secondary">Long. (aa)</th>
                        <th class="text-right px-2 py-2 border-b border-slate-200 text-secondary">MW (kDa)</th>
                        <th class="text-right px-2 py-2 border-b border-slate-200 text-secondary">pI</th>
                    </tr></thead>
                    <tbody>
                        ${top10.map((p, i) => `
                            <tr class="border-t border-slate-100 hover:bg-emerald-50 dark:hover:bg-emerald-900/10">
                                <td class="px-2 py-2 text-secondary">${i + 1}</td>
                                <td class="px-2 py-2 font-medium text-primary">${p.nombre_gen || p.locus_tag || 'N/A'}</td>
                                <td class="px-2 py-2 text-secondary line-clamp-1">${p.producto || 'N/A'}</td>
                                <td class="px-2 py-2 text-right text-primary font-mono">${DashboardRenderer.fmt(p.longitud_aa || 0)}</td>
                                <td class="px-2 py-2 text-right text-primary">${p.peso_molecular_kda ? p.peso_molecular_kda.toFixed(1) : ((p.peso_molecular_da || 0) / 1000).toFixed(1)}</td>
                                <td class="px-2 py-2 text-right text-primary">${(p.punto_isoelectrico || 0).toFixed(1)}</td>
                            </tr>
                        `).join('')}
                    </tbody>
                </table>
            </div>
        </div>
    `;

    // Chart de composicion de aminoacidos
    if (compEntries.length > 0) {
        setTimeout(() => {
            const paleta = ['#10b981', '#06b6d4', '#8b5cf6', '#f59e0b', '#ef4444', '#6366f1', '#ec4899', '#14b8a6', '#f97316', '#a855f7',
                            '#84cc16', '#22d3ee', '#c084fc', '#fb923c', '#f87171', '#818cf8', '#f472b6', '#2dd4bf', '#fbbf24', '#a78bfa'];
            DashboardRenderer.createChart('chart-3d-aa-comp', {
                type: 'bar',
                data: {
                    labels: compEntries.map(([k]) => {
                        const nombre = AA_NOMBRE_COMPLETO[k];
                        return nombre ? k + ' (' + nombre.substring(0, 3) + ')' : k;
                    }),
                    datasets: [{
                        label: 'Frecuencia (%)',
                        data: compEntries.map(([, v]) => v.porcentaje || 0),
                        backgroundColor: compEntries.map((_, i) => paleta[i % paleta.length]),
                        borderRadius: 4
                    }]
                },
                options: {
                    responsive: true,
                    plugins: {
                        legend: { display: false },
                        tooltip: {
                            callbacks: {
                                afterLabel: (ctx) => {
                                    const key = compEntries[ctx.dataIndex][0];
                                    return AA_NOMBRE_COMPLETO[key] || '';
                                }
                            }
                        }
                    },
                    scales: { x: { ticks: { font: { size: 9 }, maxRotation: 45 } }, y: { beginAtZero: true, title: { display: true, text: '%' } } }
                }
            });
        }, 100);
    }
}

function renderVisor3DSecundaria(data, container) {
    const sec = data.estructura_secundaria || {};
    const prom = sec.promedio_proteoma || {};
    const topHelice = sec.top_10_mas_helice || [];
    const topLamina = sec.top_10_mas_lamina || [];

    container.innerHTML = `
        <h3 class="text-xl font-bold text-primary mb-4">Estructura Secundaria - Helices, Laminas y Giros</h3>
        <div class="grid grid-cols-3 gap-4 mb-6">
            ${DashboardRenderer.statsCard('Helice alfa', (prom.helix || 0).toFixed(1) + '%', 'promedio del proteoma', 'emerald')}
            ${DashboardRenderer.statsCard('Lamina beta', (prom.sheet || 0).toFixed(1) + '%', 'promedio del proteoma', 'cyan')}
            ${DashboardRenderer.statsCard('Giros/bucles', (prom.turn || prom.coil || 0).toFixed(1) + '%', 'promedio del proteoma', 'amber')}
        </div>

        <!-- Pie chart distribucion -->
        <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
            <div class="bg-card rounded-xl p-5 border border-slate-200">
                <h4 class="text-sm font-semibold text-primary mb-3">Distribucion Promedio</h4>
                <canvas id="chart-3d-sec-pie" height="250"></canvas>
            </div>
            <div class="bg-card rounded-xl p-5 border border-slate-200">
                <h4 class="text-sm font-semibold text-primary mb-3">Interpretacion</h4>
                <div class="space-y-3 text-sm text-secondary">
                    <p><strong class="text-emerald-500">Helice alfa (${(prom.helix || 0).toFixed(1)}%):</strong> Estructura enrollada estabilizada por puentes de hidrogeno intracadena. Comun en proteinas transmembrana y factores de transcripcion.</p>
                    <p><strong class="text-cyan-500">Lamina beta (${(prom.sheet || 0).toFixed(1)}%):</strong> Cadenas extendidas unidas lateralmente. Predomina en enzimas y proteinas de union a sustrato.</p>
                    <p><strong class="text-amber-500">Giros/Bucles (${(prom.turn || prom.coil || 0).toFixed(1)}%):</strong> Regiones flexibles que conectan elementos de estructura secundaria. Importantes para la funcion y especificidad.</p>
                    <p class="text-xs mt-2">Prediccion: metodo Chou-Fasman (BioPython). Para mayor precision, usar herramientas como DSSP sobre estructuras 3D experimentales.</p>
                </div>
            </div>
        </div>

        <!-- Top proteinas con mas helice -->
        ${topHelice.length > 0 ? `
        <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
            <h4 class="text-sm font-semibold text-emerald-600 mb-3">Top 10 Proteinas con Mayor Contenido de Helice Alfa</h4>
            <div class="overflow-x-auto">
                <table class="w-full text-xs">
                    <thead class="bg-emerald-50 dark:bg-emerald-900/20">
                        <tr>
                            <th class="px-2 py-2 text-left text-secondary">#</th>
                            <th class="px-2 py-2 text-left text-secondary">Gen</th>
                            <th class="px-2 py-2 text-left text-secondary">Producto</th>
                            <th class="px-2 py-2 text-right text-secondary">Helice %</th>
                            <th class="px-2 py-2 text-right text-secondary">Long. (aa)</th>
                        </tr>
                    </thead>
                    <tbody>
                        ${topHelice.map((p, i) => `
                            <tr class="border-t border-slate-100">
                                <td class="px-2 py-2 text-secondary">${i + 1}</td>
                                <td class="px-2 py-2 font-mono font-bold text-primary">${p.nombre_gen || p.locus_tag || 'N/A'}</td>
                                <td class="px-2 py-2 text-secondary">${(p.producto || '').substring(0, 45)}</td>
                                <td class="px-2 py-2 text-right font-bold text-emerald-500">${(p.helix || p.porcentaje_helice || 0).toFixed(1)}%</td>
                                <td class="px-2 py-2 text-right text-primary">${DashboardRenderer.fmt(p.longitud_aa || 0)}</td>
                            </tr>
                        `).join('')}
                    </tbody>
                </table>
            </div>
        </div>` : ''}

        <!-- Top proteinas con mas lamina -->
        ${topLamina.length > 0 ? `
        <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
            <h4 class="text-sm font-semibold text-cyan-600 mb-3">Top 10 Proteinas con Mayor Contenido de Lamina Beta</h4>
            <div class="overflow-x-auto">
                <table class="w-full text-xs">
                    <thead class="bg-cyan-50 dark:bg-cyan-900/20">
                        <tr>
                            <th class="px-2 py-2 text-left text-secondary">#</th>
                            <th class="px-2 py-2 text-left text-secondary">Gen</th>
                            <th class="px-2 py-2 text-left text-secondary">Producto</th>
                            <th class="px-2 py-2 text-right text-secondary">Lamina %</th>
                            <th class="px-2 py-2 text-right text-secondary">Long. (aa)</th>
                        </tr>
                    </thead>
                    <tbody>
                        ${topLamina.map((p, i) => `
                            <tr class="border-t border-slate-100">
                                <td class="px-2 py-2 text-secondary">${i + 1}</td>
                                <td class="px-2 py-2 font-mono font-bold text-primary">${p.nombre_gen || p.locus_tag || 'N/A'}</td>
                                <td class="px-2 py-2 text-secondary">${(p.producto || '').substring(0, 45)}</td>
                                <td class="px-2 py-2 text-right font-bold text-cyan-500">${(p.sheet || p.porcentaje_lamina || 0).toFixed(1)}%</td>
                                <td class="px-2 py-2 text-right text-primary">${DashboardRenderer.fmt(p.longitud_aa || 0)}</td>
                            </tr>
                        `).join('')}
                    </tbody>
                </table>
            </div>
        </div>` : ''}
    `;

    // Pie chart
    setTimeout(() => {
        DashboardRenderer.createChart('chart-3d-sec-pie', {
            type: 'doughnut',
            data: {
                labels: ['Helice alfa', 'Lamina beta', 'Giros/Bucles'],
                datasets: [{
                    data: [prom.helix || 0, prom.sheet || 0, prom.turn || prom.coil || 0],
                    backgroundColor: ['#10b981', '#06b6d4', '#f59e0b']
                }]
            },
            options: { responsive: true, plugins: { legend: { position: 'bottom' } } }
        });
    }, 100);
}

function renderVisor3DTerciaria(data, container) {
    const terciaria = data.estructura_terciaria || {};
    const pdb = terciaria.proteinas_con_pdb || [];
    const af = terciaria.proteinas_alphafold || [];
    const todas = [...pdb, ...af];
    const cisteinas = terciaria.analisis_cisteinas || {};

    const options = todas.length > 0
        ? todas.map((p, i) => `<option value="${i}">${p.nombre_gen || p.locus_tag || 'Proteina'} - ${p.fuente || 'PDB'}: ${p.pdb_id || p.uniprot_id || 'N/A'}</option>`).join('')
        : '<option value="">No hay estructuras disponibles</option>';

    container.innerHTML = `
        <h3 class="text-xl font-bold text-primary mb-4">Estructura Terciaria - Visor 3D</h3>

        <!-- Stats de estructuras encontradas -->
        <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
            ${DashboardRenderer.statsCard('PDB (RCSB)', pdb.length, 'estructuras experimentales', 'emerald')}
            ${DashboardRenderer.statsCard('AlphaFold', af.length, 'predicciones IA', 'cyan')}
            ${DashboardRenderer.statsCard('Total 3D', todas.length, 'proteinas con estructura', 'violet')}
            ${DashboardRenderer.statsCard('Puentes S-S', cisteinas.total_puentes_disulfuro || cisteinas.total_cisteinas || 0, cisteinas.total_cisteinas ? cisteinas.total_cisteinas + ' cisteinas' : 'puentes disulfuro', 'amber')}
        </div>

        <!-- Visor 3D -->
        <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
            <div class="flex flex-wrap gap-3 items-center mb-4">
                <select id="protein-3d-select" class="flex-1 px-3 py-2 bg-slate-50 dark:bg-slate-800 rounded-lg border border-slate-200 dark:border-slate-700 text-sm text-primary" onchange="DashboardRenderer._load3DProtein()">
                    ${options}
                </select>
                <div class="flex gap-1">
                    <button onclick="DashboardRenderer._setStyle3D('cartoon')" class="px-3 py-1.5 bg-emerald-500 text-white text-xs rounded-lg">Cintas</button>
                    <button onclick="DashboardRenderer._setStyle3D('stick')" class="px-3 py-1.5 bg-slate-200 dark:bg-slate-700 text-xs rounded-lg">Varillas</button>
                    <button onclick="DashboardRenderer._setStyle3D('sphere')" class="px-3 py-1.5 bg-slate-200 dark:bg-slate-700 text-xs rounded-lg">Esferas</button>
                    <button onclick="DashboardRenderer._setStyle3D('line')" class="px-3 py-1.5 bg-slate-200 dark:bg-slate-700 text-xs rounded-lg">Lineas</button>
                    <button onclick="DashboardRenderer._setStyle3D('surface')" class="px-3 py-1.5 bg-slate-200 dark:bg-slate-700 text-xs rounded-lg">Superficie</button>
                </div>
            </div>
            <div id="protein-3d-viewer" style="width:100%; height:500px; position:relative; border-radius:12px; overflow:hidden; background:#1a1a2e;">
                <p class="text-center text-secondary text-sm pt-12">Selecciona una proteina y espera a que cargue la estructura</p>
            </div>
            <div id="protein-3d-loading" class="hidden mt-2 text-center">
                <div class="inline-block w-5 h-5 border-2 border-emerald-500 border-t-transparent rounded-full animate-spin"></div>
                <span class="text-xs text-secondary ml-2">Descargando estructura 3D...</span>
            </div>
            <div id="protein-3d-info" class="mt-3 p-4 bg-slate-50 dark:bg-slate-800 rounded-lg">
                <p class="text-xs text-secondary" id="protein-3d-info-text">Selecciona una proteina del desplegable para ver su estructura tridimensional. Las fuentes incluyen RCSB PDB (experimental) y AlphaFold (prediccion).</p>
            </div>
        </div>

        <!-- Tabla resumen de proteinas con estructura -->
        ${todas.length > 0 ? `
        <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
            <h4 class="text-sm font-semibold text-primary mb-3">Proteinas con Estructura 3D Disponible</h4>
            <div class="overflow-x-auto max-h-[350px] overflow-y-auto">
                <table class="w-full text-xs">
                    <thead class="bg-slate-50 dark:bg-slate-800 sticky top-0">
                        <tr>
                            <th class="px-2 py-2 text-left text-secondary">Gen</th>
                            <th class="px-2 py-2 text-left text-secondary">Producto</th>
                            <th class="px-2 py-2 text-center text-secondary">Fuente</th>
                            <th class="px-2 py-2 text-left text-secondary">ID</th>
                            <th class="px-2 py-2 text-right text-secondary">Long. (aa)</th>
                        </tr>
                    </thead>
                    <tbody>
                        ${todas.map(p => `
                            <tr class="border-t border-slate-100">
                                <td class="px-2 py-2 font-mono font-bold text-primary">${p.nombre_gen || p.locus_tag || 'N/A'}</td>
                                <td class="px-2 py-2 text-secondary">${(p.producto || '').substring(0, 40)}</td>
                                <td class="px-2 py-2 text-center"><span class="px-1.5 py-0.5 ${p.fuente === 'PDB' || p.pdb_id ? 'bg-emerald-100 text-emerald-700' : 'bg-cyan-100 text-cyan-700'} rounded text-[10px]">${p.fuente || (p.pdb_id ? 'PDB' : 'AlphaFold')}</span></td>
                                <td class="px-2 py-2 font-mono text-primary">${p.pdb_id || p.uniprot_id || 'N/A'}</td>
                                <td class="px-2 py-2 text-right text-primary">${DashboardRenderer.fmt(p.longitud_aa || 0)}</td>
                            </tr>
                        `).join('')}
                    </tbody>
                </table>
            </div>
        </div>` : ''}
    `;

    // Guardar datos para el visor
    DashboardRenderer._proteinas3DData = todas;
    DashboardRenderer._viewer3D = null;
    DashboardRenderer._current3DStyle = 'cartoon';
}

function renderVisor3DCuaternaria(data, container) {
    const cuat = data.estructura_cuaternaria || {};
    const complejos = cuat.complejos_detectados || [];
    const distTamano = {};
    complejos.forEach(c => {
        const n = c.num_subunidades || 0;
        const key = n <= 2 ? 'Dimero (2)' : n <= 4 ? 'Tetramero (3-4)' : n <= 8 ? 'Octamero (5-8)' : 'Grande (>8)';
        distTamano[key] = (distTamano[key] || 0) + 1;
    });

    container.innerHTML = `
        <h3 class="text-xl font-bold text-primary mb-4">Estructura Cuaternaria - Complejos Proteicos</h3>
        <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
            ${DashboardRenderer.statsCard('Complejos', DashboardRenderer.fmt(cuat.total_complejos || 0), 'detectados', 'emerald')}
            ${DashboardRenderer.statsCard('Subunidades', DashboardRenderer.fmt(cuat.total_subunidades || 0), 'proteinas en complejos', 'cyan')}
            ${DashboardRenderer.statsCard('Complejo Mayor', (cuat.complejo_mas_grande?.nombre_complejo || 'N/A').substring(0, 20), (cuat.complejo_mas_grande?.num_subunidades || 0) + ' subunidades', 'violet')}
            ${DashboardRenderer.statsCard('Promedio', (complejos.length > 0 ? (complejos.reduce((s, c) => s + (c.num_subunidades || 0), 0) / complejos.length).toFixed(1) : '0') + ' sub.', 'subunidades por complejo', 'amber')}
        </div>

        <!-- Distribucion de tamanos -->
        ${Object.keys(distTamano).length > 0 ? `
        <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
            <div class="bg-card rounded-xl p-5 border border-slate-200">
                <h4 class="text-sm font-semibold text-primary mb-3">Distribucion por Tamano de Complejo</h4>
                <canvas id="chart-3d-cuat-dist" height="250"></canvas>
            </div>
            <div class="bg-card rounded-xl p-5 border border-slate-200">
                <h4 class="text-sm font-semibold text-primary mb-3">Interpretacion</h4>
                <div class="space-y-2 text-sm text-secondary">
                    <p>Los complejos proteicos son ensamblajes de multiples subunidades que realizan funciones biologicas coordinadas.</p>
                    <p><strong class="text-emerald-500">Dimeros:</strong> Forma mas comun. Regulacion alosterica y cooperatividad.</p>
                    <p><strong class="text-cyan-500">Tetrameros:</strong> Tipicos de enzimas metabolicas y canales ionicos.</p>
                    <p><strong class="text-violet-500">Complejos grandes:</strong> Ribosomas, proteasomas, ATP sintasa.</p>
                </div>
            </div>
        </div>` : ''}

        <!-- Complejos detectados - expandibles -->
        ${complejos.length > 0 ? `
        <div class="bg-card rounded-xl p-5 border border-slate-200">
            <h4 class="text-sm font-semibold text-primary mb-3">Complejos Detectados (${complejos.length})</h4>
            <div class="space-y-3 max-h-[500px] overflow-y-auto">
                ${complejos.slice(0, 40).map((c, i) => `
                    <div class="bg-slate-50 dark:bg-slate-800 rounded-lg p-3 cursor-pointer hover:bg-slate-100 dark:hover:bg-slate-700 transition" onclick="DashboardRenderer._toggleExpandRow('cuat-detail-${i}')">
                        <div class="flex items-center justify-between mb-1">
                            <span class="font-medium text-sm text-primary">${c.nombre_complejo || 'Complejo'}</span>
                            <div class="flex items-center gap-2">
                                ${c.pdb_id ? `<span class="text-[10px] px-1.5 py-0.5 bg-emerald-100 text-emerald-700 rounded">PDB: ${c.pdb_id}</span>` : ''}
                                <span class="text-xs px-2 py-0.5 bg-emerald-500/10 text-emerald-500 rounded-full">${c.num_subunidades || 0} subunidades</span>
                            </div>
                        </div>
                        <p class="text-xs text-secondary">${(c.subunidades || []).map(s => s.nombre_gen || s.locus_tag).join(', ')}</p>
                    </div>
                    <div id="cuat-detail-${i}" class="hidden px-3 pb-3">
                        <div class="bg-white dark:bg-slate-900 rounded-lg p-3 border border-slate-200 text-xs">
                            <table class="w-full">
                                <thead><tr>
                                    <th class="text-left text-secondary px-1 py-1">Subunidad</th>
                                    <th class="text-left text-secondary px-1 py-1">Producto</th>
                                    <th class="text-right text-secondary px-1 py-1">Long. (aa)</th>
                                </tr></thead>
                                <tbody>
                                    ${(c.subunidades || []).map(s => `
                                        <tr class="border-t border-slate-100">
                                            <td class="px-1 py-1 font-mono text-primary">${s.nombre_gen || s.locus_tag || 'N/A'}</td>
                                            <td class="px-1 py-1 text-secondary">${(s.producto || '').substring(0, 40)}</td>
                                            <td class="px-1 py-1 text-right text-primary">${DashboardRenderer.fmt(s.longitud_aa || 0)}</td>
                                        </tr>
                                    `).join('')}
                                </tbody>
                            </table>
                        </div>
                    </div>
                `).join('')}
            </div>
        </div>` : '<p class="text-secondary text-center py-8">No se detectaron complejos proteicos</p>'}
    `;

    // Chart distribucion
    if (Object.keys(distTamano).length > 0) {
        setTimeout(() => {
            const entries = Object.entries(distTamano);
            DashboardRenderer.createChart('chart-3d-cuat-dist', {
                type: 'doughnut',
                data: {
                    labels: entries.map(([k]) => k),
                    datasets: [{
                        data: entries.map(([, v]) => v),
                        backgroundColor: ['#10b981', '#06b6d4', '#8b5cf6', '#f59e0b']
                    }]
                },
                options: { responsive: true, plugins: { legend: { position: 'bottom' } } }
            });
        }, 100);
    }
}

// =============================================================================
// EVOLUCION / PANGENOMA
// =============================================================================

async function loadEvolucionGenomeSelector() {
    const container = document.getElementById('evo-genome-list');
    if (!container) return;

    try {
        const resp = await fetch('/api/genomes');
        const data = await resp.json();
        if (data.success && data.genomes && data.genomes.length > 0) {
            container.innerHTML = data.genomes.map(g => {
                const sizeMb = g.length_mb || g.size_mb || 0;
                const acc = g.accession_id && g.accession_id !== 'N/A' ? g.accession_id : '';
                const desc = g.description || g.organism || g.basename.replace(/_/g, ' ');
                return `
                <label class="flex items-start gap-2 cursor-pointer bg-slate-50 dark:bg-slate-800 px-3 py-2.5 rounded-lg border border-slate-200 dark:border-slate-700 hover:border-emerald-500 transition text-xs group">
                    <input type="checkbox" class="evo-genome-check w-3.5 h-3.5 text-emerald-500 rounded mt-0.5 shrink-0" value="${g.basename}">
                    <div class="min-w-0">
                        <p class="text-primary font-medium text-xs leading-tight group-hover:text-emerald-500 transition">${desc.substring(0, 60)}</p>
                        <p class="text-[10px] text-secondary mt-0.5">${acc ? acc + ' | ' : ''}${sizeMb ? sizeMb + ' Mb | ' : ''}${g.basename.substring(0, 35)}</p>
                    </div>
                </label>`;
            }).join('');

            container.querySelectorAll('.evo-genome-check').forEach(cb => {
                cb.addEventListener('change', updateEvoCount);
            });
        } else {
            container.innerHTML = '<p class="text-sm text-secondary col-span-full text-center py-4">No hay genomas descargados. Ve a "Buscar Genomas" para descargar.</p>';
        }
    } catch {
        container.innerHTML = '<p class="text-sm text-red-500 col-span-full text-center py-4">Error al cargar genomas</p>';
    }
}

function updateEvoCount() {
    const checked = document.querySelectorAll('.evo-genome-check:checked');
    const countEl = document.getElementById('evo-selected-count');
    if (countEl) countEl.textContent = checked.length;
}

function evoSelectAll() {
    document.querySelectorAll('.evo-genome-check').forEach(cb => { cb.checked = true; });
    updateEvoCount();
}

function evoSelectNone() {
    document.querySelectorAll('.evo-genome-check').forEach(cb => { cb.checked = false; });
    updateEvoCount();
}

async function runEvolucion() {
    const checked = document.querySelectorAll('.evo-genome-check:checked');
    if (checked.length < 2) {
        showNotification('Selecciona al menos 2 bacterias para comparar', 'warning');
        return;
    }

    const genomesList = Array.from(checked).map(cb => cb.value).join(',');
    const btn = document.getElementById('run-evolucion-btn');
    const progress = document.getElementById('evo-progress');
    const consoleDiv = document.getElementById('evo-console');
    const consoleOutput = document.getElementById('evo-console-output');

    btn.disabled = true;
    btn.textContent = 'Analizando...';
    progress.classList.remove('hidden');
    consoleDiv.classList.remove('hidden');
    consoleOutput.textContent = 'Iniciando analisis de comunidad bacteriana...\n';

    try {
        const resp = await fetch('/api/run_analysis', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                script: 'analisis_evolucion',
                genome_basename: genomesList
            })
        });

        const result = await resp.json();
        progress.classList.add('hidden');
        btn.disabled = false;
        btn.textContent = 'Comparar Bacterias';

        // Show console output
        if (result.output) {
            consoleOutput.textContent += result.output;
        }
        // Auto scroll console to bottom
        consoleOutput.scrollTop = consoleOutput.scrollHeight;

        if (result.success && result.return_code === 0) {
            consoleOutput.textContent += '\n[OK] Analisis completado exitosamente.\n';
            consoleOutput.textContent += 'Haz click en las sub-pestanas del menu (Genes, Arbol, Parentesco) para ver los resultados.\n';
            showNotification('Analisis completado. Ve a las sub-pestanas para ver resultados.', 'success');
        } else {
            consoleOutput.textContent += '\n[ERROR] El analisis fallo.\n';
            showNotification('Error en el analisis. Revisa la consola.', 'error');
        }
        consoleOutput.scrollTop = consoleOutput.scrollHeight;
    } catch (error) {
        progress.classList.add('hidden');
        btn.disabled = false;
        btn.textContent = 'Comparar Bacterias';
        consoleOutput.textContent += '\n[ERROR] Error de conexion: ' + error.message + '\n';
        showNotification('Error de conexion: ' + error.message, 'error');
    }
}

async function loadEvolucionResults(tab) {
    const dashboard = document.getElementById('evo-dashboard');
    if (!dashboard) return;

    const activeTab = tab || AppState.currentTab || 'evo-pangenoma';

    dashboard.innerHTML = `
        <div class="flex items-center justify-center py-12">
            <div class="w-8 h-8 border-4 border-emerald-500 border-t-transparent rounded-full animate-spin mr-4"></div>
            <span class="text-secondary">Cargando resultados...</span>
        </div>
    `;

    try {
        const resp = await fetch('/api/evolucion_resultado');
        const result = await resp.json();

        if (!result.success || !result.data) {
            dashboard.innerHTML = `
                <div class="text-center py-12 text-secondary">
                    <div class="text-6xl mb-3">🧬</div>
                    <p class="text-lg font-medium">No hay resultados todavia</p>
                    <p class="text-sm mt-2">Primero ve a la pestana principal "Pangenoma" en el menu,<br>selecciona tus bacterias y presiona "Comparar Bacterias"</p>
                    <button onclick="navigateToItem(document.querySelector('.nav-item[data-section=evolucion]'))" class="mt-4 px-6 py-2 bg-emerald-500 hover:bg-emerald-600 text-white rounded-lg text-sm transition">
                        Ir a configurar analisis
                    </button>
                </div>
            `;
            return;
        }

        const data = result.data;
        DashboardRenderer.destroyCharts();

        if (activeTab === 'evo-pangenoma') {
            DashboardRenderer.renderEvolucionPangenoma(data, dashboard);
        } else if (activeTab === 'evo-arbol') {
            DashboardRenderer.renderEvolucionArbol(data, dashboard);
        } else if (activeTab === 'evo-matriz') {
            DashboardRenderer.renderEvolucionMatriz(data, dashboard);
        } else {
            DashboardRenderer.renderEvolucionPangenoma(data, dashboard);
        }
    } catch (error) {
        dashboard.innerHTML = `
            <div class="text-center py-12 text-red-500">
                <p>Error al cargar resultados: ${error.message}</p>
            </div>
        `;
    }
}

function loadEvolucionTab(tab) {
    AppState.currentTab = tab;
    loadEvolucionResults(tab);
}

// =============================================================================
// CHAT (seccion completa en vez de modal)
// =============================================================================

async function sendChatMessage() {
    const input = document.getElementById('chat-input');
    const messagesDiv = document.getElementById('chat-messages');
    const genome = document.getElementById('chat-genome')?.value || '';
    const message = input.value.trim();

    if (!message) return;
    input.value = '';

    // Agregar mensaje del usuario
    messagesDiv.innerHTML += `
        <div class="flex justify-end">
            <div class="bg-emerald-500 text-white rounded-xl rounded-tr-sm px-4 py-2 max-w-[75%] text-sm">${message}</div>
        </div>
    `;

    // Loading
    const loadingId = 'chat-loading-' + Date.now();
    messagesDiv.innerHTML += `
        <div id="${loadingId}" class="flex justify-start">
            <div class="bg-slate-100 dark:bg-slate-800 rounded-xl rounded-tl-sm px-4 py-2 text-sm text-secondary">
                <div class="flex items-center gap-2">
                    <div class="w-4 h-4 border-2 border-violet-500 border-t-transparent rounded-full animate-spin"></div>
                    Pensando...
                </div>
            </div>
        </div>
    `;
    messagesDiv.scrollTop = messagesDiv.scrollHeight;

    try {
        const resp = await fetch('/api/chat', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ message, context: genome ? `genome:${genome}` : 'general' })
        });

        const data = await resp.json();
        const reply = data.response || 'No se pudo obtener respuesta';

        // Reemplazar loading con respuesta
        const loadingEl = document.getElementById(loadingId);
        if (loadingEl) {
            loadingEl.innerHTML = `
                <div class="bg-slate-100 dark:bg-slate-800 rounded-xl rounded-tl-sm px-4 py-3 max-w-[85%] text-sm text-primary leading-relaxed">
                    ${reply.replace(/\n/g, '<br>')}
                </div>
            `;
        }
    } catch (error) {
        const loadingEl = document.getElementById(loadingId);
        if (loadingEl) {
            loadingEl.innerHTML = `
                <div class="bg-red-500/10 border border-red-500/20 rounded-xl px-4 py-2 text-sm text-red-500">
                    Error: ${error.message}
                </div>
            `;
        }
    }

    messagesDiv.scrollTop = messagesDiv.scrollHeight;
}

// =============================================================================
// SISTEMA DE TEMAS
// =============================================================================

function toggleTheme() {
    AppState.theme = AppState.theme === 'dark' ? 'light' : 'dark';
    localStorage.setItem('theme', AppState.theme);
    applyTheme(AppState.theme);
}

function applyTheme(theme) {
    // Actualizar clase en <html> (para Tailwind dark: prefixes) y <body>
    document.documentElement.classList.remove('light', 'dark');
    document.documentElement.classList.add(theme);
    document.body.classList.remove('light', 'dark');
    document.body.classList.add(theme);

    // Actualizar icono y texto del boton de tema
    const btn = document.getElementById('theme-toggle-btn');
    if (btn) {
        const iconEl = btn.querySelector('#theme-icon');
        if (iconEl) {
            iconEl.setAttribute('data-lucide', theme === 'dark' ? 'moon' : 'sun');
        }
        const textEl = btn.querySelector('.link-text');
        if (textEl) textEl.textContent = theme === 'dark' ? 'Tema Oscuro' : 'Tema Claro';
    }

    // Actualizar colores del indicador para que coincida con el fondo
    const indicator = document.getElementById('nav-indicator');
    if (indicator) {
        const bgColor = theme === 'dark' ? '#000000' : '#e2e8f0';
        indicator.style.background = bgColor;
        indicator.style.setProperty('--indicator-bg', bgColor);
    }

    // Re-crear iconos lucide despues del cambio
    if (typeof lucide !== 'undefined') {
        lucide.createIcons();
    }
}

// =============================================================================
// NOTIFICACIONES TOAST
// =============================================================================

function showNotification(message, type = 'info') {
    console.log(`[${type.toUpperCase()}] ${message}`);

    const container = document.getElementById('notifications-container') || createNotificationContainer();
    const notification = document.createElement('div');
    notification.className = 'notification toast animate-slide-in';

    const colors = { success: 'bg-emerald-500', error: 'bg-red-500', warning: 'bg-amber-500', info: 'bg-cyan-500' };
    const icons = { success: '✓', error: '✕', warning: '!', info: 'ⓘ' };

    notification.innerHTML = `
        <div class="flex items-center gap-3 px-4 py-3 rounded-lg shadow-lg ${colors[type]} text-white">
            <span class="text-xl font-bold">${icons[type]}</span>
            <span class="flex-1">${message}</span>
            <button onclick="this.parentElement.parentElement.remove()" class="text-white hover:text-gray-200">✕</button>
        </div>
    `;

    container.appendChild(notification);
    setTimeout(() => {
        notification.style.opacity = '0';
        setTimeout(() => notification.remove(), 300);
    }, 5000);
}

function createNotificationContainer() {
    const container = document.createElement('div');
    container.id = 'notifications-container';
    container.className = 'fixed top-4 right-4 z-[9999] space-y-2 max-w-md';
    document.body.appendChild(container);
    return container;
}

// =============================================================================
// UTILIDADES
// =============================================================================

function formatFileSize(bytes) {
    if (bytes === 0) return '0 Bytes';
    const k = 1024;
    const sizes = ['Bytes', 'KB', 'MB', 'GB'];
    const i = Math.floor(Math.log(bytes) / Math.log(k));
    return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i];
}

function formatDate(isoString) {
    if (!isoString) return 'N/A';
    const date = new Date(isoString);
    return date.toLocaleDateString('es-ES', { year: 'numeric', month: 'short', day: 'numeric', hour: '2-digit', minute: '2-digit' });
}

function formatNumber(num) {
    return new Intl.NumberFormat('es-ES').format(num);
}

function showLoading(elementId) {
    const el = document.getElementById(elementId);
    if (el) el.classList.remove('hidden');
}

function hideLoading(elementId) {
    const el = document.getElementById(elementId);
    if (el) el.classList.add('hidden');
}

// =============================================================================
// INICIALIZACION
// =============================================================================

document.addEventListener('DOMContentLoaded', () => {
    console.log('[APP] Inicializando GenomeHub...');

    // Inicializar iconos Lucide
    if (typeof lucide !== 'undefined') {
        lucide.createIcons();
    }

    // --- Event listeners del sidebar ---

    // Nav items (click en el link principal)
    document.querySelectorAll('.nav-item').forEach(item => {
        const link = item.querySelector('.nav-link');
        if (link) {
            link.addEventListener('click', (e) => {
                e.preventDefault();
                navigateToItem(item);
            });
        }
    });

    // Sub-items
    document.querySelectorAll('.submenu-item').forEach(subItem => {
        subItem.addEventListener('click', (e) => {
            e.stopPropagation();
            navigateToSubItem(subItem);
        });
    });

    // Visor 3D genome selector change
    const visor3dGenome = document.getElementById('visor3d-genome');
    if (visor3dGenome) {
        visor3dGenome.addEventListener('change', () => {
            if (AppState.currentTab && AppState.currentTab.startsWith('3d-')) {
                loadVisor3DTab(AppState.currentTab);
            }
        });
    }

    // Aplicar tema
    applyTheme(AppState.theme);

    // Posicionar indicador en el item activo inicial
    const activeItem = document.querySelector('.nav-item.active');
    if (activeItem) {
        moveIndicator(activeItem);
    }

    // Mostrar seccion inicial
    showSection('search');

    console.log('[APP] GenomeHub inicializado');
});

// Estilos CSS para animaciones
const style = document.createElement('style');
style.textContent = `
    .animate-slide-in {
        animation: slideIn 0.3s ease-out;
    }
    @keyframes slideIn {
        from { transform: translateX(100%); opacity: 0; }
        to { transform: translateX(0); opacity: 1; }
    }
    .notification { transition: opacity 0.3s ease; }
`;
document.head.appendChild(style);
