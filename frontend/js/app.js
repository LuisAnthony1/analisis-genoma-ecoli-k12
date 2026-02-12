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
    const top10peq = primaria.top_10_mas_pequenas || [];
    const categorias = primaria.categorias_funcionales || {};
    const mutaciones = primaria.mutaciones_patogenicas || {};
    const mutList = mutaciones.genes_analizados || mutaciones.mutaciones || [];

    const compEntries = Object.entries(comp).sort((a, b) => b[1].porcentaje - a[1].porcentaje);
    const hidrof = compEntries.filter(([k]) => ['Ala', 'Val', 'Ile', 'Leu', 'Met', 'Phe', 'Trp', 'Pro', 'A', 'V', 'I', 'L', 'M', 'F', 'W', 'P'].includes(k));
    const hidrofPct = hidrof.reduce((s, [, v]) => s + (v.porcentaje || 0), 0);

    // Categorias funcionales para grafico
    const catEntries = Object.entries(categorias).filter(([, v]) => (v.total || v.cantidad || v) > 0).sort((a, b) => (b[1].total || b[1].cantidad || b[1]) - (a[1].total || a[1].cantidad || a[1]));
    const catNombresES = {
        'enzimas': 'Enzimas', 'transportadores': 'Transportadores', 'reguladores': 'Reguladores',
        'estructurales': 'Estructurales', 'virulencia': 'Virulencia', 'chaperonas': 'Chaperonas',
        'replicacion': 'Replicacion', 'hipoteticas': 'Hipoteticas', 'metabolismo': 'Metabolismo',
        'enzymes': 'Enzimas', 'transporters': 'Transportadores', 'regulators': 'Reguladores',
        'structural': 'Estructurales', 'chaperones': 'Chaperonas', 'replication': 'Replicacion',
        'hypothetical': 'Hipoteticas', 'metabolism': 'Metabolismo'
    };

    container.innerHTML = `
        <h3 class="text-xl font-bold text-primary mb-2">Estructura Primaria - Secuencia de Aminoacidos</h3>
        <p class="text-sm text-secondary mb-5">La estructura primaria es la secuencia lineal de aminoacidos que forma cada proteina. Aqui analizamos la composicion del proteoma completo del genoma.</p>

        <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
            ${DashboardRenderer.statsCard('Total Proteinas', DashboardRenderer.fmt(stats.total_analizadas || stats.total_proteinas || 0), 'proteinas codificadas')}
            ${DashboardRenderer.statsCard('Longitud Promedio', DashboardRenderer.fmt(Math.round(stats.longitud_promedio_aa || 0)) + ' aa', 'aminoacidos por proteina', 'cyan')}
            ${DashboardRenderer.statsCard('Peso Molecular', DashboardRenderer.fmt(Math.round((stats.peso_molecular_promedio_da || 0) / 1000)) + ' kDa', 'promedio por proteina', 'violet')}
            ${DashboardRenderer.statsCard('Punto Isoelectrico', (stats.pi_promedio || 0).toFixed(1), 'pI promedio (' + (stats.pi_promedio < 7 ? 'acido' : 'basico') + ')', 'amber')}
        </div>

        <!-- Resumen visual: que tipo de proteinas tiene este genoma -->
        <div class="grid grid-cols-2 md:grid-cols-4 gap-3 mb-6">
            ${DashboardRenderer.statsCard('Hidrofobicas', DashboardRenderer.fmt(stats.proteinas_hidrofobicas || 0), `${((stats.proteinas_hidrofobicas || 0) / (stats.total_analizadas || 1) * 100).toFixed(0)}% del total`, 'amber')}
            ${DashboardRenderer.statsCard('Hidrofilicas', DashboardRenderer.fmt(stats.proteinas_hidrofilicas || 0), `${((stats.proteinas_hidrofilicas || 0) / (stats.total_analizadas || 1) * 100).toFixed(0)}% del total`, 'cyan')}
            ${DashboardRenderer.statsCard('Estables', DashboardRenderer.fmt(stats.proteinas_estables || 0), 'indice inestabilidad < 40', 'emerald')}
            ${DashboardRenderer.statsCard('Inestables', DashboardRenderer.fmt(stats.proteinas_inestables || 0), 'indice inestabilidad >= 40', 'red')}
        </div>

        <!-- Composicion y Categorias lado a lado -->
        <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
            <!-- Composicion de Aminoacidos -->
            ${compEntries.length > 0 ? `
            <div class="bg-card rounded-xl p-5 border border-slate-200">
                <h4 class="text-sm font-semibold text-primary mb-1">Composicion de Aminoacidos</h4>
                <p class="text-[10px] text-secondary mb-3">Frecuencia de cada aminoacido en todas las proteinas del genoma. <span class="text-amber-500 font-bold">${hidrofPct.toFixed(1)}% hidrofobicos</span>, <span class="text-cyan-500 font-bold">${(100 - hidrofPct).toFixed(1)}% hidrofilicos</span>.</p>
                <canvas id="chart-3d-aa-comp" height="250"></canvas>
            </div>` : ''}

            <!-- Categorias funcionales -->
            ${catEntries.length > 0 ? `
            <div class="bg-card rounded-xl p-5 border border-slate-200">
                <h4 class="text-sm font-semibold text-primary mb-1">Categorias Funcionales</h4>
                <p class="text-[10px] text-secondary mb-3">Clasificacion de proteinas segun su funcion en la celula.</p>
                <canvas id="chart-3d-cat-func" height="250"></canvas>
            </div>` : ''}
        </div>

        <!-- Aminoacidos visual grid -->
        ${compEntries.length > 0 ? `
        <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
            <h4 class="text-sm font-semibold text-primary mb-3">Los 20 Aminoacidos - Composicion del Proteoma</h4>
            <div class="grid grid-cols-4 md:grid-cols-5 lg:grid-cols-10 gap-2">
                ${compEntries.map(([k, v], i) => {
                    const nombre = AA_NOMBRE_COMPLETO[k] || k;
                    const pct = v.porcentaje || 0;
                    const isHidrof = ['Ala', 'Val', 'Ile', 'Leu', 'Met', 'Phe', 'Trp', 'Pro', 'A', 'V', 'I', 'L', 'M', 'F', 'W', 'P'].includes(k);
                    return `
                    <div class="text-center p-2 rounded-lg ${isHidrof ? 'bg-amber-50 dark:bg-amber-900/15 border border-amber-200' : 'bg-cyan-50 dark:bg-cyan-900/15 border border-cyan-200'}" title="${nombre} (${k}) - ${pct.toFixed(2)}% - ${isHidrof ? 'Hidrofobico' : 'Hidrofilico'}">
                        <p class="text-lg font-bold ${isHidrof ? 'text-amber-600' : 'text-cyan-600'}">${k}</p>
                        <p class="text-[9px] text-secondary truncate">${nombre}</p>
                        <p class="text-xs font-bold text-primary">${pct.toFixed(1)}%</p>
                    </div>`;
                }).join('')}
            </div>
            <div class="flex gap-4 mt-3 text-[10px] text-secondary">
                <span class="flex items-center gap-1"><span class="w-3 h-3 rounded bg-amber-100 border border-amber-200"></span> Hidrofobico (no polar)</span>
                <span class="flex items-center gap-1"><span class="w-3 h-3 rounded bg-cyan-100 border border-cyan-200"></span> Hidrofilico (polar/cargado)</span>
            </div>
        </div>` : ''}

        <!-- Mutaciones Patogenicas -->
        ${mutList.length > 0 ? `
        <div class="bg-card rounded-xl p-5 border border-red-200 mb-6">
            <h4 class="text-sm font-semibold text-red-600 mb-1">Genes de Resistencia a Antibioticos</h4>
            <p class="text-xs text-secondary mb-3">Genes asociados a resistencia antibiotica detectados en el genoma. Estos pueden conferir resistencia a farmacos.</p>
            <div class="overflow-x-auto max-h-[350px] overflow-y-auto">
                <table class="w-full text-xs">
                    <thead class="bg-red-50 dark:bg-red-900/20 sticky top-0">
                        <tr>
                            <th class="px-2 py-2 text-left text-secondary">Gen</th>
                            <th class="px-2 py-2 text-left text-secondary">Funcion</th>
                            <th class="px-2 py-2 text-left text-secondary">Categoria</th>
                            <th class="px-2 py-2 text-right text-secondary">Tamano (aa)</th>
                            <th class="px-2 py-2 text-right text-secondary">Peso (kDa)</th>
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

        <!-- Top 10 proteinas grandes y pequenas -->
        <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
            <div class="bg-card rounded-xl p-5 border border-emerald-200">
                <h4 class="text-sm font-semibold text-emerald-600 mb-3">Top 10 Proteinas Mas Grandes</h4>
                <div class="space-y-1.5">
                    ${top10.map((p, i) => {
                        const kda = p.peso_molecular_kda ? p.peso_molecular_kda : ((p.peso_molecular_da || 0) / 1000);
                        const maxAA = top10[0]?.longitud_aa || 1;
                        return `
                        <div class="group cursor-default" title="${p.producto || 'Sin anotacion'} | ${DashboardRenderer.fmt(p.longitud_aa)} aa | ${kda.toFixed(1)} kDa | pI ${(p.punto_isoelectrico || 0).toFixed(1)}">
                            <div class="flex items-center gap-2 text-xs">
                                <span class="text-secondary w-4 text-right">${i + 1}.</span>
                                <span class="font-mono font-bold text-primary w-14 truncate">${p.nombre_gen || p.locus_tag || 'N/A'}</span>
                                <div class="flex-1 bg-slate-100 dark:bg-slate-800 rounded-full h-2.5">
                                    <div class="h-2.5 rounded-full bg-emerald-500" style="width:${(p.longitud_aa / maxAA * 100).toFixed(0)}%"></div>
                                </div>
                                <span class="font-bold text-emerald-500 w-16 text-right">${DashboardRenderer.fmt(p.longitud_aa)} aa</span>
                            </div>
                            <p class="text-[10px] text-secondary pl-6 truncate hidden group-hover:block">${p.producto || 'Proteina hipotetica'}</p>
                        </div>`;
                    }).join('')}
                </div>
            </div>
            ${top10peq.length > 0 ? `
            <div class="bg-card rounded-xl p-5 border border-amber-200">
                <h4 class="text-sm font-semibold text-amber-600 mb-3">Top 10 Proteinas Mas Pequenas</h4>
                <div class="space-y-1.5">
                    ${top10peq.map((p, i) => {
                        const kda = p.peso_molecular_kda ? p.peso_molecular_kda : ((p.peso_molecular_da || 0) / 1000);
                        return `
                        <div class="group cursor-default" title="${p.producto || 'Sin anotacion'} | ${DashboardRenderer.fmt(p.longitud_aa)} aa | ${kda.toFixed(1)} kDa">
                            <div class="flex items-center gap-2 text-xs">
                                <span class="text-secondary w-4 text-right">${i + 1}.</span>
                                <span class="font-mono font-bold text-primary w-14 truncate">${p.nombre_gen || p.locus_tag || 'N/A'}</span>
                                <span class="text-amber-500 font-bold w-14 text-right">${DashboardRenderer.fmt(p.longitud_aa)} aa</span>
                                <span class="text-secondary flex-1 truncate text-[10px]">${(p.producto || 'Hipotetica').substring(0, 35)}</span>
                            </div>
                        </div>`;
                    }).join('')}
                </div>
            </div>` : ''}
        </div>
    `;

    // Charts
    setTimeout(() => {
        if (compEntries.length > 0) {
            const paleta = ['#10b981', '#06b6d4', '#8b5cf6', '#f59e0b', '#ef4444', '#6366f1', '#ec4899', '#14b8a6', '#f97316', '#a855f7',
                            '#84cc16', '#22d3ee', '#c084fc', '#fb923c', '#f87171', '#818cf8', '#f472b6', '#2dd4bf', '#fbbf24', '#a78bfa'];
            DashboardRenderer.createChart('chart-3d-aa-comp', {
                type: 'bar',
                data: {
                    labels: compEntries.map(([k]) => {
                        const nombre = AA_NOMBRE_COMPLETO[k];
                        return nombre ? nombre.substring(0, 4) : k;
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
                                title: (ctx) => {
                                    const key = compEntries[ctx[0].dataIndex][0];
                                    return `${AA_NOMBRE_COMPLETO[key] || key} (${key})`;
                                },
                                afterLabel: (ctx) => {
                                    const [k, v] = compEntries[ctx.dataIndex];
                                    return `Total: ${DashboardRenderer.fmt(v.total || v.cantidad || 0)} residuos`;
                                }
                            }
                        }
                    },
                    scales: { x: { ticks: { font: { size: 8 }, maxRotation: 60 } }, y: { beginAtZero: true, title: { display: true, text: '%' } } }
                }
            });
        }

        if (catEntries.length > 0) {
            const catColors = ['#10b981', '#06b6d4', '#8b5cf6', '#f59e0b', '#ef4444', '#6366f1', '#ec4899', '#14b8a6'];
            DashboardRenderer.createChart('chart-3d-cat-func', {
                type: 'doughnut',
                data: {
                    labels: catEntries.map(([k]) => catNombresES[k.toLowerCase()] || k),
                    datasets: [{
                        data: catEntries.map(([, v]) => v.total || v.cantidad || v),
                        backgroundColor: catColors
                    }]
                },
                options: { responsive: true, plugins: { legend: { position: 'right', labels: { font: { size: 10 } } } } }
            });
        }
    }, 100);
}

function renderVisor3DSecundaria(data, container) {
    const sec = data.estructura_secundaria || {};
    const prom = sec.promedio_proteoma || {};
    const topHelice = sec.top_10_mas_helice || [];
    const topLamina = sec.top_10_mas_lamina || [];

    // Combinar ambas listas para obtener proteinas con datos de estructura secundaria
    const allProteins = [];
    const seen = new Set();
    [...topHelice, ...topLamina].forEach(p => {
        const key = p.locus_tag || p.nombre_gen;
        if (!seen.has(key)) {
            seen.add(key);
            allProteins.push(p);
        }
    });

    container.innerHTML = `
        <h3 class="text-xl font-bold text-primary mb-2">Estructura Secundaria - Helices Alfa, Laminas Beta y Bucles</h3>
        <p class="text-sm text-secondary mb-5">La estructura secundaria son los patrones locales que forma la cadena de aminoacidos: espirales (helices alfa), flechas planas (laminas beta) y lazos flexibles (bucles).</p>

        <div class="grid grid-cols-3 gap-4 mb-6">
            ${DashboardRenderer.statsCard('Helice alfa', (prom.helix || 0).toFixed(1) + '%', 'promedio del proteoma', 'emerald')}
            ${DashboardRenderer.statsCard('Lamina beta', (prom.sheet || 0).toFixed(1) + '%', 'promedio del proteoma', 'cyan')}
            ${DashboardRenderer.statsCard('Giros/bucles', (prom.turn || prom.coil || 0).toFixed(1) + '%', 'promedio del proteoma', 'amber')}
        </div>

        <!-- Ilustracion visual animada de estructuras secundarias -->
        <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
            <h4 class="text-sm font-semibold text-primary mb-3">Como se ve cada tipo de estructura</h4>
            <div class="grid grid-cols-1 md:grid-cols-3 gap-6">
                <!-- Helice alfa -->
                <div class="text-center">
                    <div class="bg-emerald-50 dark:bg-emerald-900/15 rounded-xl p-4 mb-3" style="min-height:180px;">
                        <svg viewBox="0 0 200 160" width="200" height="160" class="mx-auto">
                            <defs>
                                <linearGradient id="helixGrad" x1="0" y1="0" x2="1" y2="0">
                                    <stop offset="0%" stop-color="#10b981"/>
                                    <stop offset="100%" stop-color="#059669"/>
                                </linearGradient>
                            </defs>
                            <!-- Helice alfa - espiral animada -->
                            <path d="M30,140 C50,140 50,120 70,120 C90,120 90,100 110,100 C130,100 130,80 150,80 C170,80 170,60 180,50"
                                  fill="none" stroke="url(#helixGrad)" stroke-width="4" stroke-linecap="round">
                                <animate attributeName="stroke-dashoffset" from="300" to="0" dur="2s" fill="freeze"/>
                                <animate attributeName="stroke-dasharray" from="0,300" to="300,0" dur="2s" fill="freeze"/>
                            </path>
                            <!-- Loops de la espiral -->
                            ${[0,1,2,3,4].map(i => `
                                <ellipse cx="${40 + i * 35}" cy="${130 - i * 20}" rx="18" ry="8"
                                    fill="none" stroke="#10b981" stroke-width="2.5" opacity="0.6"
                                    transform="rotate(-15, ${40 + i * 35}, ${130 - i * 20})">
                                    <animate attributeName="opacity" from="0" to="0.6" begin="${i * 0.3}s" dur="0.5s" fill="freeze"/>
                                </ellipse>
                            `).join('')}
                            <!-- Puentes de hidrogeno (lineas punteadas) -->
                            ${[0,1,2,3].map(i => `
                                <line x1="${45 + i * 35}" y1="${125 - i * 20}" x2="${55 + i * 35}" y2="${118 - i * 20}"
                                    stroke="#f59e0b" stroke-width="1" stroke-dasharray="2,2" opacity="0.5">
                                    <animate attributeName="opacity" from="0" to="0.5" begin="${0.5 + i * 0.3}s" dur="0.5s" fill="freeze"/>
                                </line>
                            `).join('')}
                            <text x="100" y="155" text-anchor="middle" fill="#10b981" font-size="10" font-weight="bold">3.6 aa por giro</text>
                        </svg>
                    </div>
                    <p class="text-sm font-bold text-emerald-600">Helice Alfa</p>
                    <p class="text-xs text-secondary mt-1">Espiral estabilizada por puentes de hidrogeno. Cada giro tiene 3.6 aminoacidos.</p>
                </div>

                <!-- Lamina beta -->
                <div class="text-center">
                    <div class="bg-cyan-50 dark:bg-cyan-900/15 rounded-xl p-4 mb-3" style="min-height:180px;">
                        <svg viewBox="0 0 200 160" width="200" height="160" class="mx-auto">
                            <!-- Flechas paralelas (laminas beta) -->
                            ${[0,1,2,3].map(i => `
                                <g>
                                    <rect x="20" y="${20 + i * 35}" width="140" height="14" rx="2"
                                        fill="${i % 2 === 0 ? '#06b6d4' : '#0891b2'}" opacity="0.7">
                                        <animate attributeName="width" from="0" to="140" begin="${i * 0.25}s" dur="0.6s" fill="freeze"/>
                                    </rect>
                                    <polygon points="${i % 2 === 0 ? '160,' + (20 + i * 35) + ' 175,' + (27 + i * 35) + ' 160,' + (34 + i * 35) : '20,' + (20 + i * 35) + ' 5,' + (27 + i * 35) + ' 20,' + (34 + i * 35)}"
                                        fill="${i % 2 === 0 ? '#06b6d4' : '#0891b2'}" opacity="0.7">
                                        <animate attributeName="opacity" from="0" to="0.7" begin="${0.3 + i * 0.25}s" dur="0.3s" fill="freeze"/>
                                    </polygon>
                                </g>
                            `).join('')}
                            <!-- Puentes H entre laminas -->
                            ${[0,1,2].map(i => `
                                ${[0,1,2,3,4].map(j => `
                                    <line x1="${35 + j * 30}" y1="${34 + i * 35}" x2="${35 + j * 30}" y2="${55 + i * 35}"
                                        stroke="#f59e0b" stroke-width="1" stroke-dasharray="2,2" opacity="0.4">
                                        <animate attributeName="opacity" from="0" to="0.4" begin="${1 + (i * 5 + j) * 0.05}s" dur="0.3s" fill="freeze"/>
                                    </line>
                                `).join('')}
                            `).join('')}
                            <text x="100" y="155" text-anchor="middle" fill="#06b6d4" font-size="10" font-weight="bold">Cadenas paralelas/antiparalelas</text>
                        </svg>
                    </div>
                    <p class="text-sm font-bold text-cyan-600">Lamina Beta</p>
                    <p class="text-xs text-secondary mt-1">Cadenas extendidas unidas por puentes de hidrogeno entre hebras adyacentes.</p>
                </div>

                <!-- Giros/Bucles -->
                <div class="text-center">
                    <div class="bg-amber-50 dark:bg-amber-900/15 rounded-xl p-4 mb-3" style="min-height:180px;">
                        <svg viewBox="0 0 200 160" width="200" height="160" class="mx-auto">
                            <!-- Curva irregular (bucle) -->
                            <path d="M20,130 Q40,100 60,110 T100,70 T140,90 T180,40"
                                  fill="none" stroke="#f59e0b" stroke-width="4" stroke-linecap="round">
                                <animate attributeName="stroke-dashoffset" from="400" to="0" dur="2s" fill="freeze"/>
                                <animate attributeName="stroke-dasharray" from="0,400" to="400,0" dur="2s" fill="freeze"/>
                            </path>
                            <!-- Aminoacidos como puntos -->
                            ${[[30,118],[50,105],[75,95],[100,70],[120,82],[140,90],[160,65],[180,40]].map(([x,y], i) => `
                                <circle cx="${x}" cy="${y}" r="4" fill="#f59e0b" opacity="0">
                                    <animate attributeName="opacity" from="0" to="0.8" begin="${i * 0.2}s" dur="0.3s" fill="freeze"/>
                                </circle>
                            `).join('')}
                            <text x="100" y="155" text-anchor="middle" fill="#f59e0b" font-size="10" font-weight="bold">Region flexible</text>
                        </svg>
                    </div>
                    <p class="text-sm font-bold text-amber-600">Giros y Bucles</p>
                    <p class="text-xs text-secondary mt-1">Regiones flexibles que conectan helices y laminas. Claves para la funcion enzimatica.</p>
                </div>
            </div>
        </div>

        <!-- Distribucion visual y grafico -->
        <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
            <div class="bg-card rounded-xl p-5 border border-slate-200">
                <h4 class="text-sm font-semibold text-primary mb-3">Distribucion Promedio del Proteoma</h4>
                <canvas id="chart-3d-sec-pie" height="250"></canvas>
            </div>
            <div class="bg-card rounded-xl p-5 border border-slate-200">
                <h4 class="text-sm font-semibold text-primary mb-3">Composicion por Proteina</h4>
                <p class="text-xs text-secondary mb-3">Barra visual de cada proteina analizada. <span class="text-emerald-500 font-bold">Verde</span> = helice, <span class="text-cyan-500 font-bold">azul</span> = lamina, <span class="text-amber-500 font-bold">amarillo</span> = bucle.</p>
                <div class="space-y-2 max-h-[280px] overflow-y-auto pr-1">
                    ${allProteins.slice(0, 20).map(p => {
                        const h = p.helix || p.porcentaje_helice || 0;
                        const s = p.sheet || p.porcentaje_lamina || 0;
                        const t = Math.max(0, 100 - h - s);
                        return `
                        <div class="group">
                            <div class="flex items-center gap-2 text-[10px] mb-0.5">
                                <span class="font-mono font-bold text-primary w-16 truncate">${p.nombre_gen || p.locus_tag || '?'}</span>
                                <span class="text-secondary truncate flex-1">${(p.producto || '').substring(0, 30)}</span>
                            </div>
                            <div class="flex h-3 rounded-full overflow-hidden bg-slate-100 dark:bg-slate-800" title="${(p.nombre_gen || p.locus_tag || '')} - Helice: ${h.toFixed(1)}%, Lamina: ${s.toFixed(1)}%, Bucle: ${t.toFixed(1)}%">
                                <div style="width:${h}%; background:#10b981;" class="transition-all duration-500"></div>
                                <div style="width:${s}%; background:#06b6d4;" class="transition-all duration-500"></div>
                                <div style="width:${t}%; background:#f59e0b;" class="transition-all duration-500"></div>
                            </div>
                        </div>`;
                    }).join('')}
                </div>
            </div>
        </div>

        <!-- Top proteinas por tipo -->
        <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
            ${topHelice.length > 0 ? `
            <div class="bg-card rounded-xl p-5 border border-emerald-200">
                <h4 class="text-sm font-semibold text-emerald-600 mb-3">Top 10 - Mayor Helice Alfa</h4>
                <div class="space-y-1.5">
                    ${topHelice.map((p, i) => {
                        const h = p.helix || p.porcentaje_helice || 0;
                        return `
                        <div class="flex items-center gap-2 text-xs">
                            <span class="text-secondary w-4 text-right">${i + 1}.</span>
                            <span class="font-mono font-bold text-primary w-14 truncate">${p.nombre_gen || p.locus_tag || 'N/A'}</span>
                            <div class="flex-1 bg-slate-100 dark:bg-slate-800 rounded-full h-2.5">
                                <div class="h-2.5 rounded-full bg-emerald-500" style="width:${h}%"></div>
                            </div>
                            <span class="font-bold text-emerald-500 w-12 text-right">${h.toFixed(1)}%</span>
                        </div>`;
                    }).join('')}
                </div>
            </div>` : ''}
            ${topLamina.length > 0 ? `
            <div class="bg-card rounded-xl p-5 border border-cyan-200">
                <h4 class="text-sm font-semibold text-cyan-600 mb-3">Top 10 - Mayor Lamina Beta</h4>
                <div class="space-y-1.5">
                    ${topLamina.map((p, i) => {
                        const s = p.sheet || p.porcentaje_lamina || 0;
                        return `
                        <div class="flex items-center gap-2 text-xs">
                            <span class="text-secondary w-4 text-right">${i + 1}.</span>
                            <span class="font-mono font-bold text-primary w-14 truncate">${p.nombre_gen || p.locus_tag || 'N/A'}</span>
                            <div class="flex-1 bg-slate-100 dark:bg-slate-800 rounded-full h-2.5">
                                <div class="h-2.5 rounded-full bg-cyan-500" style="width:${s}%"></div>
                            </div>
                            <span class="font-bold text-cyan-500 w-12 text-right">${s.toFixed(1)}%</span>
                        </div>`;
                    }).join('')}
                </div>
            </div>` : ''}
        </div>
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
        ? todas.map((p, i) => `<option value="${i}">${p.nombre_gen || p.locus_tag || 'Proteina'} - ${(p.producto || '').substring(0, 40)} [${p.fuente || (p.pdb_id ? 'PDB' : 'AlphaFold')}]</option>`).join('')
        : '<option value="">No hay estructuras disponibles</option>';

    container.innerHTML = `
        <h3 class="text-xl font-bold text-primary mb-2">Estructura Terciaria - Visor Molecular 3D</h3>
        <p class="text-sm text-secondary mb-5">Visualiza la forma tridimensional de las proteinas. Usa el mouse para rotar (click + arrastrar), zoom (rueda del mouse) y mover (click derecho + arrastrar).</p>

        <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
            ${DashboardRenderer.statsCard('PDB (RCSB)', pdb.length, 'estructuras experimentales', 'emerald')}
            ${DashboardRenderer.statsCard('AlphaFold', af.length, 'predicciones IA', 'cyan')}
            ${DashboardRenderer.statsCard('Total 3D', todas.length, 'proteinas con estructura', 'violet')}
            ${DashboardRenderer.statsCard('Puentes S-S', cisteinas.total_puentes_disulfuro || cisteinas.total_cisteinas || 0, cisteinas.total_cisteinas ? cisteinas.total_cisteinas + ' cisteinas' : 'puentes disulfuro', 'amber')}
        </div>

        <!-- Visor 3D principal -->
        <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
            <div class="flex flex-wrap gap-3 items-center mb-4">
                <select id="protein-3d-select" class="flex-1 min-w-[200px] px-3 py-2 bg-slate-50 dark:bg-slate-800 rounded-lg border border-slate-200 dark:border-slate-700 text-sm text-primary" onchange="DashboardRenderer._load3DProtein()">
                    ${options}
                </select>
            </div>
            <div class="flex flex-wrap gap-1 mb-3">
                <button id="btn-style-cartoon" onclick="DashboardRenderer._setStyle3D('cartoon')" class="px-3 py-1.5 bg-emerald-500 text-white text-xs rounded-lg transition">Cintas</button>
                <button id="btn-style-stick" onclick="DashboardRenderer._setStyle3D('stick')" class="px-3 py-1.5 bg-slate-200 dark:bg-slate-700 text-xs rounded-lg hover:bg-slate-300 transition">Varillas</button>
                <button id="btn-style-sphere" onclick="DashboardRenderer._setStyle3D('sphere')" class="px-3 py-1.5 bg-slate-200 dark:bg-slate-700 text-xs rounded-lg hover:bg-slate-300 transition">Esferas</button>
                <button id="btn-style-line" onclick="DashboardRenderer._setStyle3D('line')" class="px-3 py-1.5 bg-slate-200 dark:bg-slate-700 text-xs rounded-lg hover:bg-slate-300 transition">Lineas</button>
                <button id="btn-style-surface" onclick="DashboardRenderer._setStyle3D('surface')" class="px-3 py-1.5 bg-slate-200 dark:bg-slate-700 text-xs rounded-lg hover:bg-slate-300 transition">Superficie</button>
                <button id="btn-style-ss" onclick="DashboardRenderer._setStyle3D('ss')" class="px-3 py-1.5 bg-slate-200 dark:bg-slate-700 text-xs rounded-lg hover:bg-slate-300 transition">Estructura 2&deg;</button>
            </div>
            <div style="position:relative;">
                <div id="protein-3d-viewer" style="width:100%; height:550px; position:relative; border-radius:12px; overflow:hidden; background:#1a1a2e;">
                    <p class="text-slate-400 text-sm text-center pt-12">Selecciona una proteina del desplegable para ver su estructura 3D</p>
                </div>
                <!-- Play/Pause button - esquina superior derecha -->
                <button id="btn-3d-spin" onclick="DashboardRenderer._toggle3DSpin()"
                    class="absolute top-3 right-3 px-3 py-1.5 bg-black/50 hover:bg-black/70 text-white text-xs rounded-lg backdrop-blur-sm transition z-20 border border-white/20"
                    style="font-size:12px;">
                    &#9654; Play
                </button>
                <!-- Leyenda de colores -->
                <div class="absolute bottom-3 left-3 px-3 py-2 bg-black/50 backdrop-blur-sm rounded-lg text-[10px] text-white/80 z-20 border border-white/20">
                    <span class="inline-block w-2 h-2 rounded-full mr-1" style="background:#ff0000"></span>N-terminal &nbsp;
                    <span class="inline-block w-2 h-2 rounded-full mr-1" style="background:#00ff00"></span>Centro &nbsp;
                    <span class="inline-block w-2 h-2 rounded-full mr-1" style="background:#0000ff"></span>C-terminal
                </div>
                <!-- Controles info -->
                <div class="absolute top-3 left-3 px-3 py-2 bg-black/50 backdrop-blur-sm rounded-lg text-[10px] text-white/80 z-20 border border-white/20">
                    Click + arrastrar = rotar &nbsp; | &nbsp; Rueda = zoom &nbsp; | &nbsp; Click derecho = mover
                </div>
            </div>
            <div id="protein-3d-loading" class="hidden mt-2 text-center">
                <div class="inline-block w-5 h-5 border-2 border-emerald-500 border-t-transparent rounded-full animate-spin"></div>
                <span class="text-xs text-secondary ml-2">Descargando estructura 3D...</span>
            </div>
            <div class="mt-3 p-4 bg-slate-50 dark:bg-slate-800 rounded-lg">
                <p class="text-xs" id="protein-3d-info-text">Selecciona una proteina para ver su estructura tridimensional. Fuentes: RCSB PDB (experimental) y AlphaFold (prediccion IA).</p>
            </div>
        </div>

        <!-- Tabla de proteinas con estructura -->
        ${todas.length > 0 ? `
        <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
            <h4 class="text-sm font-semibold text-primary mb-1">Proteinas con Estructura 3D Disponible (${todas.length})</h4>
            <p class="text-xs text-secondary mb-3">Click en una fila para cargarla en el visor 3D.</p>
            <div class="overflow-x-auto max-h-[350px] overflow-y-auto">
                <table class="w-full text-xs">
                    <thead class="bg-slate-50 dark:bg-slate-800 sticky top-0">
                        <tr>
                            <th class="px-2 py-2 text-left text-secondary">Gen</th>
                            <th class="px-2 py-2 text-left text-secondary">Funcion</th>
                            <th class="px-2 py-2 text-center text-secondary">Fuente</th>
                            <th class="px-2 py-2 text-left text-secondary">ID</th>
                            <th class="px-2 py-2 text-right text-secondary">Tamano (aa)</th>
                        </tr>
                    </thead>
                    <tbody>
                        ${todas.map((p, i) => `
                            <tr class="border-t border-slate-100 cursor-pointer hover:bg-emerald-50 dark:hover:bg-emerald-900/10 transition" onclick="document.getElementById('protein-3d-select').value='${i}'; DashboardRenderer._load3DProtein();">
                                <td class="px-2 py-2 font-mono font-bold text-primary">${p.nombre_gen || p.locus_tag || 'N/A'}</td>
                                <td class="px-2 py-2 text-secondary">${(p.producto || 'Proteina hipotetica').substring(0, 45)}</td>
                                <td class="px-2 py-2 text-center"><span class="px-1.5 py-0.5 ${p.fuente === 'PDB' || p.pdb_id ? 'bg-emerald-100 dark:bg-emerald-900/30 text-emerald-700' : 'bg-cyan-100 dark:bg-cyan-900/30 text-cyan-700'} rounded text-[10px]">${p.fuente || (p.pdb_id ? 'PDB' : 'AlphaFold')}</span></td>
                                <td class="px-2 py-2 font-mono text-primary">${p.pdb_id || p.uniprot_id || 'N/A'}</td>
                                <td class="px-2 py-2 text-right text-primary">${DashboardRenderer.fmt(p.longitud_aa || 0)}</td>
                            </tr>
                        `).join('')}
                    </tbody>
                </table>
            </div>
        </div>` : ''}

        <!-- Buscar cualquier proteina -->
        <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
            <h4 class="text-sm font-semibold text-primary mb-1">Buscar cualquier proteina</h4>
            <p class="text-xs text-secondary mb-3">Si una proteina no aparece en la lista, puedes buscarla por su ID de PDB (ej: <code>1L2Y</code>) o de UniProt/AlphaFold (ej: <code>P0A870</code>). Tambien puedes buscar por nombre de gen en UniProt.</p>
            <div class="flex gap-2">
                <input id="manual-protein-id" type="text" placeholder="ID de PDB o UniProt (ej: 1L2Y, P0A870)" class="flex-1 px-3 py-2 bg-slate-50 dark:bg-slate-800 rounded-lg border border-slate-200 dark:border-slate-700 text-sm text-primary">
                <select id="manual-protein-source" class="px-3 py-2 bg-slate-50 dark:bg-slate-800 rounded-lg border border-slate-200 dark:border-slate-700 text-sm">
                    <option value="rcsb">PDB (RCSB)</option>
                    <option value="alphafold">AlphaFold (UniProt ID)</option>
                </select>
                <button onclick="_loadManualProtein()" class="px-4 py-2 bg-emerald-500 hover:bg-emerald-600 text-white text-sm rounded-lg transition">Cargar</button>
            </div>
            <p id="manual-protein-status" class="text-xs text-secondary mt-2 hidden"></p>
        </div>

        <!-- Info educativa -->
        <div class="bg-violet-50 dark:bg-violet-900/15 border border-violet-200 rounded-xl p-5">
            <h4 class="text-sm font-semibold text-violet-700 mb-2">Sobre la estructura terciaria</h4>
            <div class="grid grid-cols-1 md:grid-cols-2 gap-4 text-xs text-secondary">
                <div>
                    <p class="mb-2">La <strong class="text-primary">estructura terciaria</strong> es la forma 3D completa de la proteina. Determina su funcion biologica.</p>
                    <p class="mb-2"><strong>Fuerzas que la estabilizan:</strong></p>
                    <ul class="list-disc pl-4 space-y-1">
                        <li><strong class="text-amber-600">Puentes disulfuro (S-S)</strong> - enlaces covalentes entre cisteinas</li>
                        <li><strong class="text-emerald-600">Interacciones hidrofobicas</strong> - aminoacidos no polares se agrupan en el interior</li>
                        <li><strong class="text-cyan-600">Puentes de hidrogeno</strong> - entre cadenas laterales polares</li>
                        <li><strong class="text-violet-600">Interacciones ionicas</strong> - entre cargas opuestas</li>
                    </ul>
                </div>
                <div>
                    <p class="mb-2"><strong>Metodos de determinacion:</strong></p>
                    <ul class="list-disc pl-4 space-y-1">
                        <li><strong class="text-emerald-600">PDB (RCSB)</strong> - Rayos X, RMN, Cryo-EM (estructura real)</li>
                        <li><strong class="text-cyan-600">AlphaFold</strong> - Prediccion por IA de Google DeepMind (muy precisa)</li>
                    </ul>
                    <p class="mt-2"><strong>Modos de visualizacion:</strong></p>
                    <ul class="list-disc pl-4 space-y-1">
                        <li><strong>Cintas</strong> - muestra helices y laminas</li>
                        <li><strong>Varillas</strong> - muestra cada enlace quimico</li>
                        <li><strong>Esferas</strong> - muestra el tamano real de los atomos</li>
                        <li><strong>Superficie</strong> - muestra la forma externa de la proteina</li>
                        <li><strong>Estructura 2&deg;</strong> - colorea por tipo: <span class="text-emerald-500">helice</span>, <span class="text-blue-500">lamina</span>, <span class="text-amber-500">bucle</span></li>
                    </ul>
                </div>
            </div>
        </div>
    `;

    // Guardar datos para el visor
    DashboardRenderer._proteinas3DData = todas;
    DashboardRenderer._viewer3D = null;
    DashboardRenderer._current3DStyle = 'cartoon';
    DashboardRenderer._spinning3D = false;
}

async function _loadManualProtein() {
    const idInput = document.getElementById('manual-protein-id');
    const sourceSelect = document.getElementById('manual-protein-source');
    const status = document.getElementById('manual-protein-status');
    if (!idInput || !sourceSelect) return;

    const id = idInput.value.trim();
    if (!id) { status.textContent = 'Escribe un ID'; status.classList.remove('hidden'); return; }

    const source = sourceSelect.value;
    status.textContent = 'Buscando estructura...';
    status.className = 'text-xs text-amber-500 mt-2';

    // Add as a manual entry to the protein list
    const manualProtein = {
        nombre_gen: id,
        producto: `Busqueda manual (${source === 'alphafold' ? 'AlphaFold' : 'PDB'})`,
        pdb_id: id,
        uniprot_id: source === 'alphafold' ? id : '',
        source: source,
        fuente: source === 'alphafold' ? 'AlphaFold' : 'PDB',
        url: source === 'alphafold' ? `https://alphafold.ebi.ac.uk/files/AF-${id}-F1-model_v4.pdb` : ''
    };

    DashboardRenderer._proteinas3DData.push(manualProtein);
    const newIdx = DashboardRenderer._proteinas3DData.length - 1;

    // Set the select to the new index and load
    const select = document.getElementById('protein-3d-select');
    if (select) {
        const opt = document.createElement('option');
        opt.value = newIdx;
        opt.textContent = `${id} - Busqueda manual [${source === 'alphafold' ? 'AlphaFold' : 'PDB'}]`;
        select.appendChild(opt);
        select.value = newIdx;
    }

    try {
        await DashboardRenderer._load3DProtein();
        status.textContent = 'Estructura cargada correctamente';
        status.className = 'text-xs text-emerald-500 mt-2';
    } catch (e) {
        status.textContent = 'No se encontro la estructura. Verifica el ID.';
        status.className = 'text-xs text-red-500 mt-2';
    }
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
            // Si el link tiene href real (externo), dejar que el navegador lo abra
            if (link.hasAttribute('href') && link.getAttribute('href') !== '#') return;
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
