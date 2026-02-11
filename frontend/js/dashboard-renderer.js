/**
 * GenomeHub - Dashboard Renderer
 *
 * Renderiza dashboards interactivos desde datos JSON usando Chart.js.
 * Cada análisis tiene su propio dashboard con stats cards, gráficos y tablas.
 */

// Mapeo global: abreviatura de aminoácido → nombre completo en español
const AA_NOMBRE_COMPLETO = {
    'Ala': 'Alanina', 'Arg': 'Arginina', 'Asn': 'Asparagina', 'Asp': 'Aspartato',
    'Cys': 'Cisteina', 'Glu': 'Glutamato', 'Gln': 'Glutamina', 'Gly': 'Glicina',
    'His': 'Histidina', 'Ile': 'Isoleucina', 'Leu': 'Leucina', 'Lys': 'Lisina',
    'Met': 'Metionina', 'Phe': 'Fenilalanina', 'Pro': 'Prolina', 'Ser': 'Serina',
    'Thr': 'Treonina', 'Trp': 'Triptofano', 'Tyr': 'Tirosina', 'Val': 'Valina',
    'Stop': 'Codon de parada',
    // Con letra entre parentesis (formato del backend)
    'Ala (A)': 'Alanina', 'Arg (R)': 'Arginina', 'Asn (N)': 'Asparagina', 'Asp (D)': 'Aspartato',
    'Cys (C)': 'Cisteina', 'Glu (E)': 'Glutamato', 'Gln (Q)': 'Glutamina', 'Gly (G)': 'Glicina',
    'His (H)': 'Histidina', 'Ile (I)': 'Isoleucina', 'Leu (L)': 'Leucina', 'Lys (K)': 'Lisina',
    'Met (M)': 'Metionina', 'Phe (F)': 'Fenilalanina', 'Pro (P)': 'Prolina', 'Ser (S)': 'Serina',
    'Thr (T)': 'Treonina', 'Trp (W)': 'Triptofano', 'Tyr (Y)': 'Tirosina', 'Val (V)': 'Valina',
    // Por letra sola
    'A': 'Alanina', 'R': 'Arginina', 'N': 'Asparagina', 'D': 'Aspartato',
    'C': 'Cisteina', 'E': 'Glutamato', 'Q': 'Glutamina', 'G': 'Glicina',
    'H': 'Histidina', 'I': 'Isoleucina', 'L': 'Leucina', 'K': 'Lisina',
    'M': 'Metionina', 'F': 'Fenilalanina', 'P': 'Prolina', 'S': 'Serina',
    'T': 'Treonina', 'W': 'Triptofano', 'Y': 'Tirosina', 'V': 'Valina',
};

const DashboardRenderer = {

    // Charts activos para destruir al cambiar de vista
    activeCharts: [],

    /**
     * Destruir todos los charts activos
     */
    destroyCharts() {
        this.activeCharts.forEach(chart => chart.destroy());
        this.activeCharts = [];
    },

    /**
     * Crear una tarjeta de estadística
     */
    statsCard(title, value, subtitle = '', color = 'emerald') {
        return `
            <div class="bg-card rounded-xl p-5 border border-slate-200">
                <p class="text-xs font-medium text-secondary uppercase tracking-wide">${title}</p>
                <p class="text-2xl font-bold text-${color}-500 mt-1">${value}</p>
                ${subtitle ? `<p class="text-xs text-secondary mt-1">${subtitle}</p>` : ''}
            </div>
        `;
    },

    /**
     * Crear un chart y añadirlo al registro
     */
    createChart(canvasId, config) {
        const canvas = document.getElementById(canvasId);
        if (!canvas) return null;
        const chart = new Chart(canvas.getContext('2d'), config);
        this.activeCharts.push(chart);
        return chart;
    },

    /**
     * Formatear números grandes
     */
    fmt(num) {
        if (num === undefined || num === null) return 'N/A';
        return new Intl.NumberFormat('es-ES').format(num);
    },

    // =========================================================================
    // DASHBOARD DE GENES
    // =========================================================================

    renderGenesDashboard(data, container) {
        const stats = data.estadisticas_generales || {};
        const dist = data.distribucion_tamanos || {};
        const lit = data.comparacion_literatura || {};
        const extremos = data.genes_extremos || {};
        const refs = data.referencias_bibliograficas || {};
        const genesVsCds = data.analisis_genes_vs_cds || {};

        container.innerHTML = `
            <!-- Stats Cards -->
            <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
                ${this.statsCard('Total Genes (CDS)', this.fmt(stats.total_genes), 'Secuencias codificantes')}
                ${this.statsCard('Densidad Génica', (stats.densidad_genica_porcentaje || 0) + '%', 'del genoma codifica proteínas', 'cyan')}
                ${this.statsCard('GC Promedio (CDS)', (stats.contenido_gc_cds?.promedio || 0) + '%', 'Contenido G+C', 'amber')}
                ${this.statsCard('Tamaño Promedio', this.fmt(Math.round(stats.tamano_gen?.promedio_pb || 0)) + ' pb', (Math.round((stats.tamano_gen?.promedio_pb || 0) / 3)) + ' aminoácidos', 'violet')}
            </div>

            <!-- Gráficos -->
            <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
                <!-- Distribución de tamaños -->
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Distribución de Tamaños de Genes</h3>
                    <canvas id="chart-genes-sizes" height="250"></canvas>
                </div>

                <!-- Distribución por hebra -->
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Distribución por Hebra</h3>
                    <canvas id="chart-genes-strands" height="250"></canvas>
                </div>
            </div>

            <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
                <!-- Comparación con literatura -->
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Comparación con Literatura</h3>
                    ${Object.keys(lit).length > 0 ? '<canvas id="chart-genes-literature" height="250"></canvas>' : '<p class="text-secondary text-sm">Sin datos de literatura para este organismo</p>'}
                </div>

                <!-- Genes extremos -->
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Genes Extremos</h3>
                    <canvas id="chart-genes-extremes" height="250"></canvas>
                    ${extremos.gen_mas_largo ? `
                    <div class="mt-3 space-y-2 text-xs">
                        <div class="p-2 bg-emerald-50 dark:bg-emerald-900/20 rounded-lg">
                            <span class="font-bold text-emerald-600">Mas largo:</span>
                            <span class="text-primary">${extremos.gen_mas_largo.nombre || extremos.gen_mas_largo.locus_tag}</span> -
                            <span class="text-secondary">${extremos.gen_mas_largo.producto || 'Sin anotacion'}</span>
                            <span class="text-emerald-500 font-bold">(${this.fmt(extremos.gen_mas_largo.longitud_pb)} pb)</span>
                        </div>
                        <div class="p-2 bg-amber-50 dark:bg-amber-900/20 rounded-lg">
                            <span class="font-bold text-amber-600">Mas corto:</span>
                            <span class="text-primary">${extremos.gen_mas_corto.nombre || extremos.gen_mas_corto.locus_tag}</span> -
                            <span class="text-secondary">${extremos.gen_mas_corto.producto || 'Sin anotacion'}</span>
                            <span class="text-amber-500 font-bold">(${this.fmt(extremos.gen_mas_corto.longitud_pb)} pb)</span>
                        </div>
                    </div>` : ''}
                </div>
            </div>

            <!-- Top 10 genes mas largos -->
            ${extremos.top_10_mas_largos ? `
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h3 class="text-sm font-semibold text-primary mb-2">Top 10 Genes Mas Largos</h3>
                <p class="text-xs text-secondary mb-3">Los genes mas grandes del genoma. Haz click en una fila para ver mas detalles.</p>
                <div class="overflow-x-auto">
                    <table class="w-full text-sm">
                        <thead class="bg-slate-50 dark:bg-slate-800">
                            <tr>
                                <th class="px-3 py-2 text-left text-secondary font-medium">#</th>
                                <th class="px-3 py-2 text-left text-secondary font-medium">Gen</th>
                                <th class="px-3 py-2 text-left text-secondary font-medium">Proteina que codifica</th>
                                <th class="px-3 py-2 text-right text-secondary font-medium">Tamano (pb)</th>
                                <th class="px-3 py-2 text-right text-secondary font-medium">Aminoacidos</th>
                                <th class="px-3 py-2 text-center text-secondary font-medium">Hebra</th>
                            </tr>
                        </thead>
                        <tbody>
                            ${extremos.top_10_mas_largos.map((g, i) => `
                                <tr class="border-t border-slate-100 dark:border-slate-800 cursor-pointer hover:bg-emerald-50 dark:hover:bg-emerald-900/10 transition" onclick="DashboardRenderer._toggleExpandRow('gene-largo-${i}')">
                                    <td class="px-3 py-2 text-secondary">${i + 1}</td>
                                    <td class="px-3 py-2 font-mono font-bold text-primary">${g.nombre || g.locus_tag}</td>
                                    <td class="px-3 py-2 text-secondary">${g.producto || 'Sin anotacion'}</td>
                                    <td class="px-3 py-2 text-right font-bold text-emerald-500">${this.fmt(g.longitud_pb)}</td>
                                    <td class="px-3 py-2 text-right text-primary">${this.fmt(g.num_aminoacidos)}</td>
                                    <td class="px-3 py-2 text-center ${g.hebra === '+' ? 'text-emerald-500' : 'text-cyan-500'} font-bold">${g.hebra}</td>
                                </tr>
                                <tr id="gene-largo-${i}" class="hidden">
                                    <td colspan="6" class="px-4 py-3 bg-slate-50 dark:bg-slate-800/50 border-b border-slate-200">
                                        <div class="grid grid-cols-2 md:grid-cols-4 gap-3 text-xs">
                                            <div><span class="text-secondary">Locus:</span> <span class="font-mono text-primary">${g.locus_tag}</span></div>
                                            <div><span class="text-secondary">Posicion:</span> <span class="text-primary">${this.fmt(g.inicio)} - ${this.fmt(g.fin)}</span></div>
                                            <div><span class="text-secondary">GC:</span> <span class="text-primary">${g.contenido_gc}%</span></div>
                                            <div><span class="text-secondary">ID Proteina:</span> <span class="font-mono text-primary">${g.proteina_id || 'N/A'}</span></div>
                                        </div>
                                        <p class="text-xs text-secondary mt-2"><strong>Funcion:</strong> ${g.producto || 'Proteina hipotetica sin funcion conocida'}</p>
                                    </td>
                                </tr>
                            `).join('')}
                        </tbody>
                    </table>
                </div>
            </div>` : ''}

            <!-- Top 10 genes mas cortos -->
            ${extremos.top_10_mas_cortos ? `
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h3 class="text-sm font-semibold text-primary mb-2">Top 10 Genes Mas Cortos</h3>
                <p class="text-xs text-secondary mb-3">Los genes mas pequenos del genoma. Haz click en una fila para ver mas detalles.</p>
                <div class="overflow-x-auto">
                    <table class="w-full text-sm">
                        <thead class="bg-slate-50 dark:bg-slate-800">
                            <tr>
                                <th class="px-3 py-2 text-left text-secondary font-medium">#</th>
                                <th class="px-3 py-2 text-left text-secondary font-medium">Gen</th>
                                <th class="px-3 py-2 text-left text-secondary font-medium">Proteina que codifica</th>
                                <th class="px-3 py-2 text-right text-secondary font-medium">Tamano (pb)</th>
                                <th class="px-3 py-2 text-right text-secondary font-medium">Aminoacidos</th>
                                <th class="px-3 py-2 text-center text-secondary font-medium">Hebra</th>
                            </tr>
                        </thead>
                        <tbody>
                            ${extremos.top_10_mas_cortos.map((g, i) => `
                                <tr class="border-t border-slate-100 dark:border-slate-800 cursor-pointer hover:bg-amber-50 dark:hover:bg-amber-900/10 transition" onclick="DashboardRenderer._toggleExpandRow('gene-corto-${i}')">
                                    <td class="px-3 py-2 text-secondary">${i + 1}</td>
                                    <td class="px-3 py-2 font-mono font-bold text-primary">${g.nombre || g.locus_tag}</td>
                                    <td class="px-3 py-2 text-secondary">${g.producto || 'Sin anotacion'}</td>
                                    <td class="px-3 py-2 text-right font-bold text-amber-500">${this.fmt(g.longitud_pb)}</td>
                                    <td class="px-3 py-2 text-right text-primary">${this.fmt(g.num_aminoacidos)}</td>
                                    <td class="px-3 py-2 text-center ${g.hebra === '+' ? 'text-emerald-500' : 'text-cyan-500'} font-bold">${g.hebra}</td>
                                </tr>
                                <tr id="gene-corto-${i}" class="hidden">
                                    <td colspan="6" class="px-4 py-3 bg-slate-50 dark:bg-slate-800/50 border-b border-slate-200">
                                        <div class="grid grid-cols-2 md:grid-cols-4 gap-3 text-xs">
                                            <div><span class="text-secondary">Locus:</span> <span class="font-mono text-primary">${g.locus_tag}</span></div>
                                            <div><span class="text-secondary">Posicion:</span> <span class="text-primary">${this.fmt(g.inicio)} - ${this.fmt(g.fin)}</span></div>
                                            <div><span class="text-secondary">GC:</span> <span class="text-primary">${g.contenido_gc}%</span></div>
                                            <div><span class="text-secondary">ID Proteina:</span> <span class="font-mono text-primary">${g.proteina_id || 'N/A'}</span></div>
                                        </div>
                                        <p class="text-xs text-secondary mt-2"><strong>Funcion:</strong> ${g.producto || 'Proteina hipotetica sin funcion conocida'}</p>
                                    </td>
                                </tr>
                            `).join('')}
                        </tbody>
                    </table>
                </div>
            </div>` : ''}

            <!-- Referencias bibliográficas -->
            ${Object.keys(refs).length > 0 ? `
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h3 class="text-sm font-semibold text-primary mb-3">Referencias Bibliograficas</h3>
                <div class="space-y-2 text-sm text-secondary">
                    ${Object.entries(refs).map(([k, v]) => `<p><span class="font-medium">${k}:</span> ${v}</p>`).join('')}
                </div>
            </div>` : ''}

            <!-- GC por Posicion de Codon -->
            ${data.gc_por_posicion_codon ? (() => {
                const gcp = data.gc_por_posicion_codon;
                return `
                <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                    <h3 class="text-sm font-semibold text-primary mb-4">Contenido GC por Posicion del Codon</h3>
                    <div class="grid grid-cols-3 gap-4 mb-4">
                        ${this.statsCard('GC1 (1ra pos.)', gcp.gc1_porcentaje + '%', 'Mas conservada', 'emerald')}
                        ${this.statsCard('GC2 (2da pos.)', gcp.gc2_porcentaje + '%', 'Afecta aminoacido', 'cyan')}
                        ${this.statsCard('GC3 (3ra pos.)', gcp.gc3_porcentaje + '%', 'Refleja sesgo codonico', 'amber')}
                    </div>
                    <canvas id="chart-genes-gc-pos" height="200"></canvas>
                    <p class="text-xs text-secondary mt-3">${gcp.interpretacion || ''}</p>
                </div>`;
            })() : ''}

            <!-- Genes Esenciales (Keio Collection) -->
            ${data.genes_esenciales ? (() => {
                const esc = data.genes_esenciales;
                const cats = esc.categorias_esenciales || {};
                const encontrados = esc.genes_encontrados || [];
                const noEncontrados = esc.genes_no_encontrados || [];
                const statsEsc = esc.estadisticas || {};
                return `
                <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                    <h3 class="text-sm font-semibold text-primary mb-4">Genes Esenciales (Keio Collection)</h3>
                    <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-4">
                        ${this.statsCard('Encontrados', esc.total_encontrados + '/' + esc.total_referencia, 'genes esenciales', 'emerald')}
                        ${this.statsCard('Cobertura', esc.porcentaje_encontrados + '%', 'del set de referencia', 'cyan')}
                        ${this.statsCard('GC Promedio', (statsEsc.gc_promedio_esenciales || 'N/A') + '%', 'en genes esenciales', 'amber')}
                        ${this.statsCard('Long. Promedio', this.fmt(Math.round(statsEsc.longitud_promedio_esenciales || 0)) + ' pb', 'genes esenciales', 'violet')}
                    </div>
                    <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-4">
                        <div>
                            <h4 class="text-xs font-semibold text-secondary mb-2">Categorias Funcionales</h4>
                            <canvas id="chart-genes-esenciales-cat" height="250"></canvas>
                        </div>
                        <div>
                            <h4 class="text-xs font-semibold text-secondary mb-2">Genes Encontrados (${encontrados.length})</h4>
                            <div class="max-h-64 overflow-y-auto">
                                <table class="w-full text-xs">
                                    <thead class="bg-slate-50 dark:bg-slate-800 sticky top-0">
                                        <tr>
                                            <th class="px-2 py-1 text-left text-secondary">Gen</th>
                                            <th class="px-2 py-1 text-left text-secondary">Producto</th>
                                            <th class="px-2 py-1 text-right text-secondary">pb</th>
                                            <th class="px-2 py-1 text-right text-secondary">GC%</th>
                                        </tr>
                                    </thead>
                                    <tbody>
                                        ${encontrados.slice(0, 40).map(g => `
                                            <tr class="border-t border-slate-100">
                                                <td class="px-2 py-1 font-mono text-primary">${g.nombre_gen}</td>
                                                <td class="px-2 py-1 text-secondary line-clamp-1">${(g.producto || '').substring(0, 40)}</td>
                                                <td class="px-2 py-1 text-right text-primary">${this.fmt(g.longitud_pb)}</td>
                                                <td class="px-2 py-1 text-right text-primary">${g.contenido_gc}%</td>
                                            </tr>
                                        `).join('')}
                                    </tbody>
                                </table>
                            </div>
                        </div>
                    </div>
                    ${noEncontrados.length > 0 ? `
                    <div class="p-3 bg-amber-50 dark:bg-amber-900/10 rounded-lg">
                        <p class="text-xs text-secondary"><strong class="text-amber-600">No encontrados (${noEncontrados.length}):</strong>
                        <span class="font-mono">${noEncontrados.slice(0, 15).join(', ')}${noEncontrados.length > 15 ? '...' : ''}</span></p>
                    </div>` : ''}
                    <p class="text-xs text-secondary mt-2">Fuente: ${esc.fuente || 'Keio collection'}</p>
                </div>`;
            })() : ''}

            <!-- Densidad Genica por Ventana -->
            ${data.densidad_por_ventana ? (() => {
                const dens = data.densidad_por_ventana;
                return `
                <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                    <h3 class="text-sm font-semibold text-primary mb-4">Densidad Genica a lo Largo del Genoma (ventana ${Math.round(dens.ventana_pb / 1000)} kb)</h3>
                    <div class="grid grid-cols-3 gap-4 mb-4">
                        ${this.statsCard('Promedio', dens.promedio_densidad + '%', 'densidad genica', 'emerald')}
                        ${this.statsCard('Maxima', dens.max_densidad + '%', 'region mas densa', 'cyan')}
                        ${this.statsCard('Minima', dens.min_densidad + '%', 'region menos densa', 'amber')}
                    </div>
                    <canvas id="chart-genes-densidad" height="200"></canvas>
                </div>`;
            })() : ''}

            <!-- Descargar datos -->
            <div class="flex gap-2">
                <button onclick="DashboardRenderer.downloadJSON('genes')" class="px-4 py-2 bg-emerald-500 hover:bg-emerald-600 text-white text-sm rounded-lg transition">Descargar JSON</button>
            </div>
        `;

        // Renderizar charts
        setTimeout(() => {
            // Distribución de tamaños
            if (Object.keys(dist).length > 0) {
                const labels = Object.keys(dist);
                const values = labels.map(k => dist[k].cantidad || 0);
                this.createChart('chart-genes-sizes', {
                    type: 'bar',
                    data: {
                        labels: labels.map(l => l.replace(' pares de bases', ' pb')),
                        datasets: [{
                            label: 'Genes',
                            data: values,
                            backgroundColor: ['#10b981', '#06b6d4', '#8b5cf6', '#f59e0b', '#ef4444', '#6366f1'],
                            borderRadius: 6
                        }]
                    },
                    options: { responsive: true, plugins: { legend: { display: false } }, scales: { y: { beginAtZero: true } } }
                });
            }

            // Hebras
            if (stats.distribucion_hebras) {
                this.createChart('chart-genes-strands', {
                    type: 'doughnut',
                    data: {
                        labels: ['Hebra directa (+)', 'Hebra complementaria (-)'],
                        datasets: [{
                            data: [stats.distribucion_hebras.forward, stats.distribucion_hebras.reverse],
                            backgroundColor: ['#10b981', '#06b6d4']
                        }]
                    },
                    options: { responsive: true, plugins: { legend: { position: 'bottom' } } }
                });
            }

            // Literatura
            if (Object.keys(lit).length > 0) {
                const litLabels = Object.keys(lit);
                const observados = litLabels.map(k => lit[k].observado);
                const literatura = litLabels.map(k => lit[k].literatura);
                this.createChart('chart-genes-literature', {
                    type: 'bar',
                    data: {
                        labels: litLabels,
                        datasets: [
                            { label: 'Observado', data: observados, backgroundColor: '#10b981', borderRadius: 4 },
                            { label: 'Literatura', data: literatura, backgroundColor: '#64748b', borderRadius: 4 }
                        ]
                    },
                    options: { responsive: true, plugins: { legend: { position: 'top' } }, scales: { y: { beginAtZero: true } } }
                });
            }

            // Genes extremos
            if (extremos.gen_mas_largo && extremos.gen_mas_corto) {
                const largo = extremos.gen_mas_largo;
                const corto = extremos.gen_mas_corto;
                this.createChart('chart-genes-extremes', {
                    type: 'bar',
                    data: {
                        labels: [largo.nombre_gen || largo.locus_tag || 'Más largo', corto.nombre_gen || corto.locus_tag || 'Más corto'],
                        datasets: [{
                            label: 'Longitud (pb)',
                            data: [largo.longitud_pb, corto.longitud_pb],
                            backgroundColor: ['#10b981', '#f59e0b'],
                            borderRadius: 6
                        }]
                    },
                    options: { indexAxis: 'y', responsive: true, plugins: { legend: { display: false } } }
                });
            }

            // GC por posicion de codon
            if (data.gc_por_posicion_codon) {
                const gcp = data.gc_por_posicion_codon;
                this.createChart('chart-genes-gc-pos', {
                    type: 'bar',
                    data: {
                        labels: ['GC1 (1ra posicion)', 'GC2 (2da posicion)', 'GC3 (3ra posicion)'],
                        datasets: [{
                            label: 'Contenido GC (%)',
                            data: [gcp.gc1_porcentaje, gcp.gc2_porcentaje, gcp.gc3_porcentaje],
                            backgroundColor: ['#10b981', '#06b6d4', '#f59e0b'],
                            borderRadius: 6
                        }]
                    },
                    options: {
                        responsive: true,
                        plugins: { legend: { display: false } },
                        scales: { y: { beginAtZero: true, max: 100, title: { display: true, text: 'GC (%)' } } }
                    }
                });
            }

            // Genes esenciales - categorias funcionales
            if (data.genes_esenciales) {
                const cats = data.genes_esenciales.categorias_esenciales || {};
                const catEntries = Object.entries(cats).sort((a, b) => b[1] - a[1]);
                if (catEntries.length > 0) {
                    const paleta = ['#10b981', '#06b6d4', '#8b5cf6', '#f59e0b', '#ef4444', '#6366f1', '#ec4899', '#14b8a6', '#f97316', '#a855f7'];
                    this.createChart('chart-genes-esenciales-cat', {
                        type: 'doughnut',
                        data: {
                            labels: catEntries.map(([k]) => k),
                            datasets: [{
                                data: catEntries.map(([, v]) => v),
                                backgroundColor: catEntries.map((_, i) => paleta[i % paleta.length])
                            }]
                        },
                        options: { responsive: true, plugins: { legend: { position: 'right', labels: { font: { size: 10 }, boxWidth: 12 } } } }
                    });
                }
            }

            // Densidad genica por ventana (line/area chart)
            if (data.densidad_por_ventana) {
                const dens = data.densidad_por_ventana;
                const ventanas = dens.ventanas || [];
                if (ventanas.length > 0) {
                    this.createChart('chart-genes-densidad', {
                        type: 'line',
                        data: {
                            labels: ventanas.map(v => Math.round(v.inicio_pb / 1000) + ' kb'),
                            datasets: [{
                                label: 'Densidad genica (%)',
                                data: ventanas.map(v => v.densidad_porcentaje),
                                borderColor: '#10b981',
                                backgroundColor: 'rgba(16, 185, 129, 0.1)',
                                fill: true,
                                tension: 0.3,
                                pointRadius: 0,
                                borderWidth: 2
                            }]
                        },
                        options: {
                            responsive: true,
                            plugins: { legend: { display: false } },
                            scales: {
                                x: { title: { display: true, text: 'Posicion en el genoma' }, ticks: { maxTicksLimit: 15, font: { size: 9 } } },
                                y: { beginAtZero: true, title: { display: true, text: 'Densidad (%)' } }
                            }
                        }
                    });
                }
            }
        }, 100);
    },

    // =========================================================================
    // DASHBOARD DE CODONES
    // =========================================================================

    renderCodonesDashboard(data, container) {
        const inicio = data.codones_inicio || {};
        const parada = data.codones_parada || {};
        const gc = data.contenido_gc || {};
        const codones64 = data.conteo_64_codones || {};

        container.innerHTML = `
            <!-- Stats Cards -->
            <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
                ${this.statsCard('Total ATG', this.fmt(inicio.conteo_total), inicio.densidad_por_kb + ' /kb')}
                ${this.statsCard('Codones de Parada', this.fmt(parada.total_codones_parada), 'TAA + TAG + TGA', 'red')}
                ${this.statsCard('Contenido GC', (gc.contenido_gc_porcentaje || 0) + '%', 'G + C', 'cyan')}
                ${this.statsCard('Codones Únicos', codones64.codones_detalle ? Object.keys(codones64.codones_detalle).filter(c => codones64.codones_detalle[c].conteo > 0).length + '/64' : 'N/A', 'de los 64 posibles', 'violet')}
            </div>

            <!-- Gráficos -->
            <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
                <!-- Codones de parada (pie) -->
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Codones de Parada</h3>
                    <canvas id="chart-codones-stop" height="250"></canvas>
                </div>

                <!-- Composición nucleotídica -->
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Composición Nucleotídica</h3>
                    <canvas id="chart-codones-nucleotides" height="250"></canvas>
                </div>
            </div>

            <!-- Gráfico de 64 codones -->
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h3 class="text-sm font-semibold text-primary mb-4">Frecuencia de los 64 Codones</h3>
                <div style="height: 400px; overflow-x: auto;">
                    <canvas id="chart-codones-64" height="380"></canvas>
                </div>
            </div>

            <!-- Tabla de 64 codones -->
            ${codones64.codones_detalle ? `
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h3 class="text-sm font-semibold text-primary mb-4">Tabla Completa de 64 Codones</h3>
                <div class="overflow-x-auto max-h-96 overflow-y-auto">
                    <table class="w-full text-sm">
                        <thead class="bg-slate-50 dark:bg-slate-800 sticky top-0">
                            <tr>
                                <th class="px-3 py-2 text-left text-secondary font-medium">Codón</th>
                                <th class="px-3 py-2 text-left text-secondary font-medium">Aminoácido</th>
                                <th class="px-3 py-2 text-right text-secondary font-medium">Conteo</th>
                                <th class="px-3 py-2 text-right text-secondary font-medium">Frecuencia %</th>
                                <th class="px-3 py-2 text-right text-secondary font-medium">Densidad /kb</th>
                            </tr>
                        </thead>
                        <tbody>
                            ${Object.entries(codones64.codones_detalle).sort((a, b) => b[1].conteo - a[1].conteo).map(([codon, info]) => `
                                <tr class="border-t border-slate-100 dark:border-slate-800">
                                    <td class="px-3 py-2 font-mono font-bold text-primary">${codon}</td>
                                    <td class="px-3 py-2 text-secondary"><span title="${AA_NOMBRE_COMPLETO[info.aminoacido] || info.aminoacido}" class="cursor-help border-b border-dotted border-slate-400">${info.aminoacido}</span></td>
                                    <td class="px-3 py-2 text-right text-primary">${this.fmt(info.conteo)}</td>
                                    <td class="px-3 py-2 text-right text-primary">${info.frecuencia_porcentaje}%</td>
                                    <td class="px-3 py-2 text-right text-secondary">${info.densidad_por_kb}</td>
                                </tr>
                            `).join('')}
                        </tbody>
                    </table>
                </div>
            </div>` : ''}

            <!-- RSCU (Uso Relativo de Codones Sinonimos) -->
            ${data.rscu ? (() => {
                const rscu = data.rscu;
                const preferidos = rscu.codones_preferidos || {};
                const evitados = rscu.codones_evitados || {};
                const raros = rscu.codones_raros || {};
                const rscuPorCodon = rscu.rscu_por_codon || {};
                const topRscu = Object.entries(rscuPorCodon).sort((a, b) => b[1] - a[1]).slice(0, 25);
                return `
                <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                    <h3 class="text-sm font-semibold text-primary mb-4">RSCU - Uso Relativo de Codones Sinonimos</h3>
                    <div class="grid grid-cols-3 gap-4 mb-4">
                        ${this.statsCard('Preferidos', rscu.total_preferidos || Object.keys(preferidos).length, 'RSCU > 1.3', 'emerald')}
                        ${this.statsCard('Evitados', rscu.total_evitados || Object.keys(evitados).length, 'RSCU < 0.7', 'amber')}
                        ${this.statsCard('Raros', rscu.total_raros || Object.keys(raros).length, 'RSCU < 0.5', 'red')}
                    </div>
                    <div class="mb-4" style="height: 400px;">
                        <canvas id="chart-codones-rscu" height="380"></canvas>
                    </div>
                    ${Object.keys(raros).length > 0 ? `
                    <div class="p-3 bg-red-50 dark:bg-red-900/10 rounded-lg mb-3">
                        <h4 class="text-xs font-semibold text-red-600 mb-2">Codones Raros (RSCU < 0.5)</h4>
                        <div class="flex flex-wrap gap-2">
                            ${Object.entries(raros).map(([codon, val]) => `
                                <span class="px-2 py-1 bg-red-100 dark:bg-red-900/30 rounded text-xs font-mono text-red-700">${codon}: ${typeof val === 'number' ? val.toFixed(2) : val}</span>
                            `).join('')}
                        </div>
                    </div>` : ''}
                    <p class="text-xs text-secondary">RSCU = 1.0 indica uso sin sesgo. Valores > 1.3 indican preferencia, < 0.7 indican evitamiento. El sesgo codonico refleja presion selectiva por eficiencia traduccional.</p>
                </div>`;
            })() : ''}

            <!-- Numero Efectivo de Codones (Nc) -->
            ${data.numero_efectivo_codones ? (() => {
                const nc = data.numero_efectivo_codones;
                const ncVal = nc.nc || nc.valor || 0;
                const ncColor = ncVal < 35 ? 'amber' : ncVal < 45 ? 'cyan' : 'emerald';
                return `
                <div class="bg-gradient-to-br from-${ncColor}-500/5 to-${ncColor}-500/10 rounded-xl p-6 border border-${ncColor}-500/20 mb-6">
                    <h3 class="text-sm font-semibold text-primary mb-4">Numero Efectivo de Codones (Nc)</h3>
                    <div class="flex items-center gap-6">
                        <div class="text-center">
                            <p class="text-5xl font-bold text-${ncColor}-500">${typeof ncVal === 'number' ? ncVal.toFixed(1) : ncVal}</p>
                            <p class="text-xs text-secondary mt-1">Nc</p>
                        </div>
                        <div class="flex-1">
                            <div class="w-full bg-slate-200 dark:bg-slate-700 rounded-full h-3 mb-2">
                                <div class="h-3 rounded-full bg-gradient-to-r from-amber-500 via-cyan-500 to-emerald-500" style="width: ${Math.min(100, ((ncVal - 20) / 41) * 100)}%"></div>
                            </div>
                            <div class="flex justify-between text-xs text-secondary">
                                <span>20 (sesgo maximo)</span>
                                <span>40</span>
                                <span>61 (sin sesgo)</span>
                            </div>
                        </div>
                    </div>
                    <p class="text-sm text-primary mt-4">${nc.interpretacion || (ncVal < 35 ? 'Sesgo codonico fuerte: el organismo usa un subconjunto restringido de codones sinonimos.' : ncVal < 45 ? 'Sesgo codonico moderado: preferencia parcial por ciertos codones.' : 'Sesgo codonico bajo: uso relativamente uniforme de codones sinonimos.')}</p>
                    <p class="text-xs text-secondary mt-2">Referencia: ${nc.referencia || 'Wright F. (1990) Gene 87:23-29'}</p>
                </div>`;
            })() : ''}

            <!-- Sesgo Codonico por Aminoacido -->
            ${data.sesgo_por_aminoacido ? (() => {
                const sesgo = data.sesgo_por_aminoacido;
                const aaEntries = Object.entries(sesgo).filter(([, v]) => v.codones && v.codones.length > 1);
                return `
                <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                    <h3 class="text-sm font-semibold text-primary mb-4">Sesgo Codonico por Aminoacido</h3>
                    <p class="text-xs text-secondary mb-3">Para cada aminoacido con codones sinonimos: codon preferido (verde) vs evitado (rojo).</p>
                    <div class="overflow-x-auto max-h-[500px] overflow-y-auto">
                        <table class="w-full text-sm">
                            <thead class="bg-slate-50 dark:bg-slate-800 sticky top-0">
                                <tr>
                                    <th class="px-3 py-2 text-left text-secondary font-medium">Aminoacido</th>
                                    <th class="px-3 py-2 text-center text-secondary font-medium">Codon Preferido</th>
                                    <th class="px-3 py-2 text-right text-secondary font-medium">% Uso</th>
                                    <th class="px-3 py-2 text-center text-secondary font-medium">Codon Evitado</th>
                                    <th class="px-3 py-2 text-right text-secondary font-medium">% Uso</th>
                                    <th class="px-3 py-2 text-right text-secondary font-medium">Num. Codones</th>
                                </tr>
                            </thead>
                            <tbody>
                                ${aaEntries.map(([aa, info]) => `
                                    <tr class="border-t border-slate-100 dark:border-slate-800">
                                        <td class="px-3 py-2 font-medium text-primary"><span title="${AA_NOMBRE_COMPLETO[aa] || aa}" class="cursor-help border-b border-dotted border-slate-400">${aa}</span></td>
                                        <td class="px-3 py-2 text-center font-mono font-bold text-emerald-500">${info.preferido || '-'}</td>
                                        <td class="px-3 py-2 text-right text-emerald-600 font-medium">${info.preferido_porcentaje ? info.preferido_porcentaje + '%' : '-'}</td>
                                        <td class="px-3 py-2 text-center font-mono font-bold text-red-500">${info.evitado || '-'}</td>
                                        <td class="px-3 py-2 text-right text-red-600 font-medium">${info.evitado_porcentaje ? info.evitado_porcentaje + '%' : '-'}</td>
                                        <td class="px-3 py-2 text-right text-secondary">${info.codones ? info.codones.length : '-'}</td>
                                    </tr>
                                `).join('')}
                            </tbody>
                        </table>
                    </div>
                </div>`;
            })() : ''}

            <div class="flex gap-2">
                <button onclick="DashboardRenderer.downloadJSON('codones')" class="px-4 py-2 bg-emerald-500 hover:bg-emerald-600 text-white text-sm rounded-lg transition">Descargar JSON</button>
            </div>
        `;

        setTimeout(() => {
            // Codones de parada (doughnut)
            if (parada.conteos) {
                this.createChart('chart-codones-stop', {
                    type: 'doughnut',
                    data: {
                        labels: ['TAA', 'TAG', 'TGA'],
                        datasets: [{
                            data: [parada.conteos.TAA, parada.conteos.TAG, parada.conteos.TGA],
                            backgroundColor: ['#10b981', '#f59e0b', '#ef4444']
                        }]
                    },
                    options: { responsive: true, plugins: { legend: { position: 'bottom' } } }
                });
            }

            // Nucleótidos
            if (gc.conteo_nucleotidos) {
                const n = gc.conteo_nucleotidos;
                this.createChart('chart-codones-nucleotides', {
                    type: 'bar',
                    data: {
                        labels: ['Adenina (A)', 'Timina (T)', 'Guanina (G)', 'Citosina (C)'],
                        datasets: [{
                            label: 'Conteo',
                            data: [n.A, n.T, n.G, n.C],
                            backgroundColor: ['#ef4444', '#f59e0b', '#10b981', '#06b6d4'],
                            borderRadius: 6
                        }]
                    },
                    options: { responsive: true, plugins: { legend: { display: false } }, scales: { y: { beginAtZero: true } } }
                });
            }

            // 64 codones (bar chart)
            if (codones64.codones_detalle) {
                const sorted = Object.entries(codones64.codones_detalle).sort((a, b) => b[1].conteo - a[1].conteo);
                const labels = sorted.map(([c]) => c);
                const values = sorted.map(([, i]) => i.conteo);
                const aaInfo = sorted.map(([, i]) => i.aminoacido);
                const colors = sorted.map(([, i]) => {
                    if (i.aminoacido === 'Stop') return '#ef4444';
                    if (i.aminoacido === 'Met (M)') return '#f59e0b';
                    return '#10b981';
                });

                this.createChart('chart-codones-64', {
                    type: 'bar',
                    data: {
                        labels,
                        datasets: [{
                            label: 'Conteo',
                            data: values,
                            backgroundColor: colors,
                            borderRadius: 2
                        }]
                    },
                    options: {
                        responsive: true,
                        maintainAspectRatio: false,
                        plugins: {
                            legend: { display: false },
                            tooltip: {
                                callbacks: {
                                    afterLabel: (ctx) => {
                                        const aa = aaInfo[ctx.dataIndex];
                                        const nombre = AA_NOMBRE_COMPLETO[aa] || aa;
                                        return aa + (nombre !== aa ? ' - ' + nombre : '');
                                    }
                                }
                            }
                        },
                        scales: {
                            x: { ticks: { font: { size: 9, family: 'monospace' }, maxRotation: 90 } },
                            y: { beginAtZero: true }
                        }
                    }
                });
            }

            // RSCU horizontal bar chart
            if (data.rscu) {
                const rscuPorCodon = data.rscu.rscu_por_codon || {};
                const topRscu = Object.entries(rscuPorCodon).sort((a, b) => b[1] - a[1]).slice(0, 25);
                if (topRscu.length > 0) {
                    this.createChart('chart-codones-rscu', {
                        type: 'bar',
                        data: {
                            labels: topRscu.map(([c]) => c),
                            datasets: [{
                                label: 'RSCU',
                                data: topRscu.map(([, v]) => v),
                                backgroundColor: topRscu.map(([, v]) => v > 1.3 ? '#10b981' : v < 0.7 ? '#ef4444' : '#94a3b8'),
                                borderRadius: 4
                            }]
                        },
                        options: {
                            indexAxis: 'y',
                            responsive: true,
                            maintainAspectRatio: false,
                            plugins: { legend: { display: false } },
                            scales: {
                                x: { beginAtZero: true, title: { display: true, text: 'RSCU' } },
                                y: { ticks: { font: { size: 10, family: 'monospace' } } }
                            }
                        }
                    });
                }
            }
        }, 100);
    },

    // =========================================================================
    // DASHBOARD DE DISTANCIAS INTERGÉNICAS
    // =========================================================================

    renderDistanciasDashboard(data, container) {
        const stats = data.estadisticas_generales || {};
        const hebra = data.estadisticas_por_hebra || {};
        const dist = data.distribucion_tipos || {};
        const top20 = data.top_20_regiones_grandes || [];

        container.innerHTML = `
            <!-- Stats Cards -->
            <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
                ${this.statsCard('Total Pares', this.fmt(stats.total_pares_genes), 'pares de genes analizados')}
                ${this.statsCard('Solapados', this.fmt(stats.genes_solapados), stats.porcentaje_solapados + '% del total', 'red')}
                ${this.statsCard('Distancia Promedio', this.fmt(Math.round(stats.distancia_promedio || 0)) + ' pb', 'entre genes espaciados', 'cyan')}
                ${this.statsCard('Regiones Grandes', this.fmt(stats.regiones_grandes), '> 500 pb', 'amber')}
            </div>

            <!-- Gráficos -->
            <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
                <!-- Distribución de tipos -->
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Distribución de Tipos de Distancias</h3>
                    <canvas id="chart-dist-types" height="250"></canvas>
                </div>

                <!-- Por hebra -->
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Análisis por Hebra</h3>
                    <canvas id="chart-dist-strand" height="250"></canvas>
                </div>
            </div>

            <!-- Tabla Top 20 regiones grandes -->
            ${top20.length > 0 ? `
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h3 class="text-sm font-semibold text-primary mb-4">Top 20 Regiones Intergenicas Grandes (posibles islas de patogenicidad)</h3>
                <p class="text-xs text-secondary mb-3">Haz click en una fila para ver mas detalles sobre los genes flanqueantes.</p>
                <div class="overflow-x-auto max-h-[600px] overflow-y-auto">
                    <table class="w-full text-sm">
                        <thead class="bg-slate-50 dark:bg-slate-800 sticky top-0">
                            <tr>
                                <th class="px-3 py-2 text-left text-secondary font-medium">#</th>
                                <th class="px-3 py-2 text-left text-secondary font-medium">Gen 1</th>
                                <th class="px-3 py-2 text-left text-secondary font-medium">Codifica (Gen 1)</th>
                                <th class="px-3 py-2 text-left text-secondary font-medium">Gen 2</th>
                                <th class="px-3 py-2 text-left text-secondary font-medium">Codifica (Gen 2)</th>
                                <th class="px-3 py-2 text-right text-secondary font-medium">Distancia (pb)</th>
                                <th class="px-3 py-2 text-center text-secondary font-medium">Hebras</th>
                            </tr>
                        </thead>
                        <tbody>
                            ${top20.map((r, i) => `
                                <tr class="border-t border-slate-100 dark:border-slate-800 cursor-pointer hover:bg-emerald-50 dark:hover:bg-emerald-900/10 transition" onclick="DashboardRenderer._toggleExpandRow('dist-detail-${i}')">
                                    <td class="px-3 py-2 text-secondary">${i + 1}</td>
                                    <td class="px-3 py-2 font-mono font-bold text-primary">${r.gen1}</td>
                                    <td class="px-3 py-2 text-secondary text-xs">${r.gen1_producto || 'Sin anotacion'}</td>
                                    <td class="px-3 py-2 font-mono font-bold text-primary">${r.gen2}</td>
                                    <td class="px-3 py-2 text-secondary text-xs">${r.gen2_producto || 'Sin anotacion'}</td>
                                    <td class="px-3 py-2 text-right font-bold text-emerald-500">${this.fmt(r.distancia_pb)}</td>
                                    <td class="px-3 py-2 text-center text-secondary">${r.hebra_gen1} -> ${r.hebra_gen2}</td>
                                </tr>
                                <tr id="dist-detail-${i}" class="hidden">
                                    <td colspan="7" class="px-4 py-3 bg-slate-50 dark:bg-slate-800/50 border-b border-slate-200">
                                        <div class="grid grid-cols-1 md:grid-cols-2 gap-4 text-xs">
                                            <div class="p-3 bg-white dark:bg-slate-900 rounded-lg border border-slate-200">
                                                <p class="font-bold text-emerald-500 mb-1">Gen 1: ${r.gen1}</p>
                                                <p><span class="text-secondary">Proteina:</span> <span class="text-primary font-medium">${r.gen1_producto || 'No anotado'}</span></p>
                                                <p><span class="text-secondary">Hebra:</span> ${r.hebra_gen1 === '+' ? 'Directa (+)' : 'Complementaria (-)'}</p>
                                            </div>
                                            <div class="p-3 bg-white dark:bg-slate-900 rounded-lg border border-slate-200">
                                                <p class="font-bold text-cyan-500 mb-1">Gen 2: ${r.gen2}</p>
                                                <p><span class="text-secondary">Proteina:</span> <span class="text-primary font-medium">${r.gen2_producto || 'No anotado'}</span></p>
                                                <p><span class="text-secondary">Hebra:</span> ${r.hebra_gen2 === '+' ? 'Directa (+)' : 'Complementaria (-)'}</p>
                                            </div>
                                        </div>
                                        <p class="text-xs text-secondary mt-2">Distancia intergenica: <strong class="text-amber-500">${this.fmt(r.distancia_pb)} pb</strong> - Region grande que podria contener elementos regulatorios, islas genomicas o genes no anotados.</p>
                                    </td>
                                </tr>
                            `).join('')}
                        </tbody>
                    </table>
                </div>
            </div>` : ''}

            <div class="flex gap-2">
                <button onclick="DashboardRenderer.downloadJSON('distancias')" class="px-4 py-2 bg-emerald-500 hover:bg-emerald-600 text-white text-sm rounded-lg transition">Descargar JSON</button>
            </div>
        `;

        setTimeout(() => {
            // Distribución de tipos
            if (Object.keys(dist).length > 0) {
                const entries = Object.entries(dist).sort((a, b) => b[1] - a[1]);
                this.createChart('chart-dist-types', {
                    type: 'doughnut',
                    data: {
                        labels: entries.map(([k]) => k),
                        datasets: [{
                            data: entries.map(([, v]) => v),
                            backgroundColor: ['#10b981', '#06b6d4', '#8b5cf6', '#f59e0b', '#ef4444', '#6366f1', '#ec4899']
                        }]
                    },
                    options: { responsive: true, plugins: { legend: { position: 'right', labels: { font: { size: 11 } } } } }
                });
            }

            // Por hebra
            this.createChart('chart-dist-strand', {
                type: 'bar',
                data: {
                    labels: ['Misma hebra', 'Diferente hebra'],
                    datasets: [
                        {
                            label: 'Total de pares',
                            data: [hebra.misma_hebra_total || 0, hebra.diferente_hebra_total || 0],
                            backgroundColor: ['#10b981', '#06b6d4'],
                            borderRadius: 6
                        }
                    ]
                },
                options: { responsive: true, plugins: { legend: { display: false } }, scales: { y: { beginAtZero: true } } }
            });
        }, 100);
    },

    // =========================================================================
    // DASHBOARD DE COMPARACIÓN
    // =========================================================================

    renderComparacionDashboard(data, container) {
        const metricas = data.metricas_generales || {};
        const virulencia = data.genes_virulencia || {};
        const organismos = data.organismos_comparados || {};
        const usoCodones = data.uso_codones || {};
        const distTamanos = data.distribucion_tamanos || {};
        const distGc = data.distribucion_gc || {};
        const interpretacionIA = data.interpretacion_ia || '';

        const ecoli = metricas.ecoli || {};
        const salmonella = metricas.salmonella || {};

        const nombre1 = ecoli.nombre || organismos.organismo_1?.nombre || 'Genoma 1';
        const nombre2 = salmonella.nombre || organismos.organismo_2?.nombre || 'Genoma 2';

        const tamStats1 = distTamanos.ecoli?.estadisticas || {};
        const tamStats2 = distTamanos.salmonella?.estadisticas || {};
        const gcStats1 = distGc.ecoli?.estadisticas || {};
        const gcStats2 = distGc.salmonella?.estadisticas || {};

        container.innerHTML = `
            <!-- Stats lado a lado -->
            <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
                <div class="bg-card rounded-xl p-5 border border-emerald-500/30">
                    <h3 class="text-lg font-bold text-emerald-500 mb-3">${nombre1}</h3>
                    <div class="space-y-2 text-sm">
                        <p><span class="text-secondary">Longitud:</span> <span class="font-bold text-primary">${this.fmt(ecoli.longitud_genoma_pb)} pb</span></p>
                        <p><span class="text-secondary">Genes CDS:</span> <span class="font-bold text-primary">${this.fmt(ecoli.total_genes_cds)}</span></p>
                        <p><span class="text-secondary">GC:</span> <span class="font-bold text-primary">${ecoli.contenido_gc_porcentaje}%</span></p>
                        <p><span class="text-secondary">Densidad:</span> <span class="font-bold text-primary">${ecoli.densidad_genica_porcentaje}%</span></p>
                        <p><span class="text-secondary">Tam. promedio gen:</span> <span class="font-bold text-primary">${this.fmt(tamStats1.promedio || ecoli.tamano_promedio_gen_pb || 0)} pb</span></p>
                    </div>
                </div>
                <div class="bg-card rounded-xl p-5 border border-cyan-500/30">
                    <h3 class="text-lg font-bold text-cyan-500 mb-3">${nombre2}</h3>
                    <div class="space-y-2 text-sm">
                        <p><span class="text-secondary">Longitud:</span> <span class="font-bold text-primary">${this.fmt(salmonella.longitud_genoma_pb)} pb</span></p>
                        <p><span class="text-secondary">Genes CDS:</span> <span class="font-bold text-primary">${this.fmt(salmonella.total_genes_cds)}</span></p>
                        <p><span class="text-secondary">GC:</span> <span class="font-bold text-primary">${salmonella.contenido_gc_porcentaje}%</span></p>
                        <p><span class="text-secondary">Densidad:</span> <span class="font-bold text-primary">${salmonella.densidad_genica_porcentaje}%</span></p>
                        <p><span class="text-secondary">Tam. promedio gen:</span> <span class="font-bold text-primary">${this.fmt(tamStats2.promedio || salmonella.tamano_promedio_gen_pb || 0)} pb</span></p>
                    </div>
                </div>
            </div>

            <!-- Graficos comparativos - fila 1 -->
            <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Comparacion de Metricas Generales</h3>
                    <canvas id="chart-compare-metrics" height="250"></canvas>
                </div>
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Genes de Virulencia</h3>
                    <canvas id="chart-compare-virulence" height="250"></canvas>
                </div>
            </div>

            <!-- Graficos comparativos - fila 2 -->
            <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Distribucion de Tamanos de Genes</h3>
                    <canvas id="chart-compare-sizes" height="250"></canvas>
                </div>
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Distribucion de GC por Gen</h3>
                    <canvas id="chart-compare-gc" height="250"></canvas>
                </div>
            </div>

            <!-- Uso de codones -->
            ${usoCodones.ecoli ? `
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h3 class="text-sm font-semibold text-primary mb-4">Comparacion de Uso de Codones (Top 10)</h3>
                <canvas id="chart-compare-codons" height="300"></canvas>
            </div>` : ''}

            <!-- Categorias de virulencia -->
            ${virulencia.ecoli && virulencia.salmonella ? `
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h3 class="text-sm font-semibold text-primary mb-4">Desglose de Genes de Virulencia por Categoria</h3>
                <canvas id="chart-compare-vir-cats" height="250"></canvas>
            </div>` : ''}

            <!-- Interpretacion IA -->
            ${interpretacionIA ? `
            <div class="bg-gradient-to-br from-violet-500/5 to-fuchsia-500/5 rounded-xl p-6 border border-violet-500/20 mb-6">
                <div class="flex items-center gap-3 mb-4">
                    <div class="w-10 h-10 rounded-xl bg-gradient-to-br from-violet-500 to-fuchsia-500 flex items-center justify-center text-xl shadow-lg">✨</div>
                    <div>
                        <h3 class="font-bold text-primary">Interpretacion IA</h3>
                        <p class="text-xs text-secondary">Analisis biologico generado por Gemini AI</p>
                    </div>
                </div>
                <div class="text-sm text-primary leading-relaxed space-y-3">
                    ${interpretacionIA.split('\n').filter(p => p.trim()).map(p =>
                        `<p>${p.replace(/\*\*(.*?)\*\*/g, '<strong>$1</strong>').replace(/\*(.*?)\*/g, '<em>$1</em>')}</p>`
                    ).join('')}
                </div>
            </div>` : ''}
        `;

        setTimeout(() => {
            // Metricas comparativas
            this.createChart('chart-compare-metrics', {
                type: 'bar',
                data: {
                    labels: ['Longitud (Mb)', 'Genes CDS', 'GC %', 'Densidad %'],
                    datasets: [
                        {
                            label: nombre1,
                            data: [
                                (ecoli.longitud_genoma_pb || 0) / 1e6,
                                ecoli.total_genes_cds || 0,
                                ecoli.contenido_gc_porcentaje || 0,
                                ecoli.densidad_genica_porcentaje || 0
                            ],
                            backgroundColor: '#10b981',
                            borderRadius: 4
                        },
                        {
                            label: nombre2,
                            data: [
                                (salmonella.longitud_genoma_pb || 0) / 1e6,
                                salmonella.total_genes_cds || 0,
                                salmonella.contenido_gc_porcentaje || 0,
                                salmonella.densidad_genica_porcentaje || 0
                            ],
                            backgroundColor: '#06b6d4',
                            borderRadius: 4
                        }
                    ]
                },
                options: { responsive: true, plugins: { legend: { position: 'top' } }, scales: { y: { beginAtZero: true } } }
            });

            // Virulencia total
            if (virulencia.ecoli && virulencia.salmonella) {
                this.createChart('chart-compare-virulence', {
                    type: 'bar',
                    data: {
                        labels: [nombre1, nombre2],
                        datasets: [{
                            label: 'Genes de virulencia',
                            data: [virulencia.ecoli.total || 0, virulencia.salmonella.total || 0],
                            backgroundColor: ['#10b981', '#06b6d4'],
                            borderRadius: 6
                        }]
                    },
                    options: { responsive: true, plugins: { legend: { display: false } }, scales: { y: { beginAtZero: true } } }
                });
            }

            // Distribucion de tamanos
            const rangos1 = distTamanos.ecoli?.rangos || {};
            const rangos2 = distTamanos.salmonella?.rangos || {};
            if (Object.keys(rangos1).length > 0) {
                const labels = Object.keys(rangos1);
                this.createChart('chart-compare-sizes', {
                    type: 'bar',
                    data: {
                        labels: labels.map(l => l + ' pb'),
                        datasets: [
                            { label: nombre1, data: labels.map(l => rangos1[l] || 0), backgroundColor: '#10b981', borderRadius: 4 },
                            { label: nombre2, data: labels.map(l => rangos2[l] || 0), backgroundColor: '#06b6d4', borderRadius: 4 }
                        ]
                    },
                    options: { responsive: true, plugins: { legend: { position: 'top' } }, scales: { y: { beginAtZero: true } } }
                });
            }

            // Distribucion GC
            const gcRangos1 = distGc.ecoli?.rangos || {};
            const gcRangos2 = distGc.salmonella?.rangos || {};
            if (Object.keys(gcRangos1).length > 0) {
                const labels = Object.keys(gcRangos1);
                this.createChart('chart-compare-gc', {
                    type: 'bar',
                    data: {
                        labels,
                        datasets: [
                            { label: nombre1, data: labels.map(l => gcRangos1[l] || 0), backgroundColor: '#10b981', borderRadius: 4 },
                            { label: nombre2, data: labels.map(l => gcRangos2[l] || 0), backgroundColor: '#06b6d4', borderRadius: 4 }
                        ]
                    },
                    options: { responsive: true, plugins: { legend: { position: 'top' } }, scales: { y: { beginAtZero: true } } }
                });
            }

            // Uso de codones (top 10)
            if (usoCodones.ecoli && usoCodones.salmonella) {
                const top1 = usoCodones.ecoli.top_10 || [];
                const top2 = usoCodones.salmonella.top_10 || [];
                const allCodons = [...new Set([...top1.map(c => c.codon), ...top2.map(c => c.codon)])].slice(0, 12);
                const freq1Map = {};
                const freq2Map = {};
                top1.forEach(c => freq1Map[c.codon] = c.frecuencia);
                top2.forEach(c => freq2Map[c.codon] = c.frecuencia);

                this.createChart('chart-compare-codons', {
                    type: 'bar',
                    data: {
                        labels: allCodons,
                        datasets: [
                            { label: nombre1, data: allCodons.map(c => freq1Map[c] || 0), backgroundColor: '#10b981', borderRadius: 3 },
                            { label: nombre2, data: allCodons.map(c => freq2Map[c] || 0), backgroundColor: '#06b6d4', borderRadius: 3 }
                        ]
                    },
                    options: {
                        responsive: true,
                        plugins: { legend: { position: 'top' } },
                        scales: {
                            x: { ticks: { font: { family: 'monospace', weight: 'bold' } } },
                            y: { beginAtZero: true, title: { display: true, text: 'Frecuencia %' } }
                        }
                    }
                });
            }

            // Categorias de virulencia
            if (virulencia.ecoli?.categorias && virulencia.salmonella?.categorias) {
                const cats1 = virulencia.ecoli.categorias || {};
                const cats2 = virulencia.salmonella.categorias || {};
                const allCats = [...new Set([...Object.keys(cats1), ...Object.keys(cats2)])];
                this.createChart('chart-compare-vir-cats', {
                    type: 'bar',
                    data: {
                        labels: allCats,
                        datasets: [
                            { label: nombre1, data: allCats.map(c => cats1[c] || 0), backgroundColor: '#10b981', borderRadius: 4 },
                            { label: nombre2, data: allCats.map(c => cats2[c] || 0), backgroundColor: '#06b6d4', borderRadius: 4 }
                        ]
                    },
                    options: {
                        indexAxis: 'y',
                        responsive: true,
                        plugins: { legend: { position: 'top' } },
                        scales: { x: { beginAtZero: true } }
                    }
                });
            }
        }, 100);
    },

    // =========================================================================
    // DASHBOARD DE PROTEINAS (con visor 3D)
    // =========================================================================

    renderProteinasDashboard(data, container) {
        const meta = data.metadata || {};
        const primaria = data.estructura_primaria || {};
        const stats = primaria.estadisticas_generales || {};
        const comp = primaria.composicion_aminoacidos || {};
        const categorias = primaria.categorias_funcionales || {};
        const peptidos = primaria.peptidos_senal || {};
        const mutaciones = primaria.mutaciones_patogenicas || [];
        const secundaria = data.estructura_secundaria || {};
        const terciaria = data.estructura_terciaria || {};
        const cuaternaria = data.estructura_cuaternaria || {};
        const pdbList = terciaria.proteinas_con_pdb || [];
        const afList = terciaria.proteinas_alphafold || [];
        const cisteinas = terciaria.analisis_cisteinas || {};
        const eduInfo = terciaria.info_educativa || {};
        const topGrandes = primaria.top_10_mas_grandes || [];
        const topPequenas = primaria.top_10_mas_pequenas || [];
        const promSec = secundaria.promedio_proteoma || {};

        // Unir proteinas con estructura disponible para el visor 3D
        const proteinasConEstructura = [];
        pdbList.forEach(p => {
            if (p.encontrado && p.pdb_ids && p.pdb_ids.length > 0) {
                proteinasConEstructura.push({
                    label: `${p.nombre_gen || p.locus_tag} - ${(p.producto || '').substring(0, 50)}`,
                    pdb_id: p.pdb_ids[0],
                    source: 'rcsb',
                    nombre_gen: p.nombre_gen,
                    producto: p.producto
                });
            }
        });
        afList.forEach(p => {
            if (p.encontrado && p.alphafold_url) {
                proteinasConEstructura.push({
                    label: `${p.nombre_gen || p.locus_tag} - ${(p.producto || '').substring(0, 50)} (AlphaFold)`,
                    pdb_id: p.uniprot_id,
                    source: 'alphafold',
                    url: p.alphafold_url,
                    nombre_gen: p.nombre_gen,
                    producto: p.producto
                });
            }
        });

        // Generar opciones del selector 3D
        const viewer3dOptions = proteinasConEstructura.length > 0
            ? proteinasConEstructura.map((p, i) => `<option value="${i}">${p.label}</option>`).join('')
            : '<option value="">No hay estructuras disponibles</option>';

        container.innerHTML = `
            <!-- Stats Cards -->
            <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
                ${this.statsCard('Total Proteinas', this.fmt(stats.total_analizadas || meta.total_proteinas), 'Proteoma completo')}
                ${this.statsCard('Longitud Promedio', (stats.longitud_promedio_aa || 0) + ' aa', 'Mediana: ' + (stats.longitud_mediana_aa || 0) + ' aa', 'cyan')}
                ${this.statsCard('Peso Molecular', this.fmt(Math.round((stats.peso_molecular_promedio_da || 0) / 1000)) + ' kDa', 'Promedio', 'amber')}
                ${this.statsCard('pI Promedio', stats.pi_promedio || 0, 'Acidas: ' + this.fmt(stats.proteinas_acidas || 0) + ' / Basicas: ' + this.fmt(stats.proteinas_basicas || 0), 'violet')}
            </div>

            <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
                ${this.statsCard('GRAVY', stats.gravy_promedio || 0, (stats.proteinas_hidrofobicas || 0) + ' hidrofobicas / ' + (stats.proteinas_hidrofilicas || 0) + ' hidrofilicas', 'emerald')}
                ${this.statsCard('Estabilidad', this.fmt(stats.proteinas_estables || 0) + ' estables', (stats.proteinas_inestables || 0) + ' inestables (indice>40)', 'cyan')}
                ${this.statsCard('Peptidos Senal', this.fmt(peptidos.total_detectados || 0), (peptidos.porcentaje || 0) + '% del proteoma', 'amber')}
                ${this.statsCard('Complejos', this.fmt(cuaternaria.total_complejos || 0), (cuaternaria.total_subunidades || 0) + ' subunidades totales', 'violet')}
            </div>

            <!-- VISOR 3D DE PROTEINAS -->
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h3 class="text-sm font-semibold text-primary mb-2">Visor 3D de Estructura Proteica</h3>
                <p class="text-xs text-secondary mb-3">Selecciona una proteina con estructura conocida (base de datos PDB o prediccion AlphaFold) para visualizarla en 3D. Arrastra para rotar, rueda del mouse para acercar/alejar.</p>

                <div class="flex flex-col md:flex-row gap-3 mb-3">
                    <select id="protein-3d-select" class="flex-1 px-3 py-2 bg-slate-50 dark:bg-slate-800 rounded-lg border border-slate-200 dark:border-slate-700 text-sm text-primary" onchange="DashboardRenderer._load3DProtein()">
                        <option value="">-- Seleccionar proteina --</option>
                        ${viewer3dOptions}
                    </select>
                    <div class="flex gap-2">
                        <button onclick="DashboardRenderer._setStyle3D('cartoon')" class="px-3 py-2 bg-emerald-500/10 hover:bg-emerald-500/20 text-emerald-500 text-xs font-medium rounded-lg transition" title="Cintas (helices y laminas)">Cintas</button>
                        <button onclick="DashboardRenderer._setStyle3D('stick')" class="px-3 py-2 bg-cyan-500/10 hover:bg-cyan-500/20 text-cyan-500 text-xs font-medium rounded-lg transition" title="Varillas (enlaces atomicos)">Varillas</button>
                        <button onclick="DashboardRenderer._setStyle3D('sphere')" class="px-3 py-2 bg-violet-500/10 hover:bg-violet-500/20 text-violet-500 text-xs font-medium rounded-lg transition" title="Esferas (radio de van der Waals)">Esferas</button>
                        <button onclick="DashboardRenderer._setStyle3D('line')" class="px-3 py-2 bg-amber-500/10 hover:bg-amber-500/20 text-amber-500 text-xs font-medium rounded-lg transition" title="Lineas (esqueleto de la proteina)">Lineas</button>
                        <button onclick="DashboardRenderer._setStyle3D('surface')" class="px-3 py-2 bg-pink-500/10 hover:bg-pink-500/20 text-pink-500 text-xs font-medium rounded-lg transition" title="Superficie molecular">Superficie</button>
                    </div>
                </div>

                <div id="protein-3d-container" style="width:100%; height:450px; position:relative; border-radius:12px; overflow:hidden; background: #1a1a2e;">
                    <div id="protein-3d-placeholder" class="flex items-center justify-center h-full">
                        <div class="text-center">
                            <div class="text-5xl mb-3 opacity-50">🧬</div>
                            <p class="text-slate-400 text-sm">${proteinasConEstructura.length > 0 ? 'Selecciona una proteina de la lista para ver su estructura en 3D' : 'Ejecuta el analisis con conexion a internet para obtener estructuras de PDB y AlphaFold'}</p>
                            ${proteinasConEstructura.length > 0 ? `<p class="text-slate-500 text-xs mt-1">${proteinasConEstructura.length} estructuras disponibles para visualizar</p>` : ''}
                        </div>
                    </div>
                    <div id="protein-3d-viewer" style="width:100%; height:100%; position:absolute; top:0; left:0; display:none;"></div>
                </div>
                <div id="protein-3d-info" class="hidden mt-3 p-3 bg-slate-50 dark:bg-slate-800 rounded-lg">
                    <p id="protein-3d-name" class="text-sm font-semibold text-primary"></p>
                    <p id="protein-3d-detail" class="text-xs text-secondary mt-1"></p>
                </div>
                <div id="protein-3d-loading" class="hidden mt-2 text-center">
                    <div class="inline-block w-5 h-5 border-2 border-emerald-500 border-t-transparent rounded-full animate-spin"></div>
                    <span class="text-xs text-secondary ml-2">Descargando y renderizando estructura 3D...</span>
                </div>
            </div>

            <!-- Graficos principales -->
            <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
                <!-- Composicion de aminoacidos -->
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Composicion de Aminoacidos del Proteoma</h3>
                    <canvas id="chart-prot-aminoacidos" height="280"></canvas>
                </div>

                <!-- Categorias funcionales -->
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Categorias Funcionales</h3>
                    <canvas id="chart-prot-categorias" height="280"></canvas>
                </div>
            </div>

            <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
                <!-- Estructura secundaria -->
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Estructura Secundaria Promedio</h3>
                    <canvas id="chart-prot-secundaria" height="280"></canvas>
                </div>

                <!-- Distribucion de pI -->
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Acidas vs Basicas (pI)</h3>
                    <canvas id="chart-prot-pi" height="280"></canvas>
                </div>
            </div>

            <!-- Top 10 proteinas mas grandes -->
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h3 class="text-sm font-semibold text-primary mb-4">Top 10 Proteinas mas Grandes</h3>
                <div class="overflow-x-auto">
                    <table class="w-full text-xs">
                        <thead>
                            <tr>
                                <th class="text-left px-2 py-2 border-b border-slate-200 text-secondary">#</th>
                                <th class="text-left px-2 py-2 border-b border-slate-200 text-secondary">Gen</th>
                                <th class="text-left px-2 py-2 border-b border-slate-200 text-secondary">Producto</th>
                                <th class="text-right px-2 py-2 border-b border-slate-200 text-secondary">Longitud (aa)</th>
                                <th class="text-right px-2 py-2 border-b border-slate-200 text-secondary">MW (kDa)</th>
                                <th class="text-right px-2 py-2 border-b border-slate-200 text-secondary">pI</th>
                                <th class="text-right px-2 py-2 border-b border-slate-200 text-secondary">GRAVY</th>
                            </tr>
                        </thead>
                        <tbody>
                            ${topGrandes.map((p, i) => `
                                <tr class="hover:bg-slate-50 dark:hover:bg-slate-800">
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700">${i+1}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 font-mono font-medium text-primary">${p.nombre_gen || p.locus_tag}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-secondary max-w-[250px] truncate">${p.producto || ''}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-right font-bold text-emerald-500">${this.fmt(p.longitud_aa)}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-right">${(p.peso_molecular_da / 1000).toFixed(1)}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-right">${p.punto_isoelectrico || '-'}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-right">${p.gravy || '-'}</td>
                                </tr>
                            `).join('')}
                        </tbody>
                    </table>
                </div>
            </div>

            <!-- Mutaciones Patogenicas -->
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h3 class="text-sm font-semibold text-primary mb-2">Genes Criticos - Resistencia a Antibioticos</h3>
                <p class="text-xs text-secondary mb-4">Genes involucrados en resistencia a antibioticos. Se analiza su presencia y posiciones clave de mutacion.</p>
                <div class="overflow-x-auto">
                    <table class="w-full text-xs">
                        <thead>
                            <tr>
                                <th class="text-left px-2 py-2 border-b border-slate-200 text-secondary">Gen</th>
                                <th class="text-left px-2 py-2 border-b border-slate-200 text-secondary">Descripcion</th>
                                <th class="text-center px-2 py-2 border-b border-slate-200 text-secondary">Estado</th>
                                <th class="text-right px-2 py-2 border-b border-slate-200 text-secondary">Longitud (aa)</th>
                                <th class="text-left px-2 py-2 border-b border-slate-200 text-secondary">Posiciones Clave</th>
                            </tr>
                        </thead>
                        <tbody>
                            ${mutaciones.map(m => {
                                const estado = m.encontrado
                                    ? (m.diferencia_longitud !== null && Math.abs(m.diferencia_longitud) < 20
                                        ? '<span class="px-2 py-0.5 bg-emerald-500/10 text-emerald-500 rounded-full text-[10px]">Normal</span>'
                                        : '<span class="px-2 py-0.5 bg-amber-500/10 text-amber-500 rounded-full text-[10px]">Longitud atipica</span>')
                                    : '<span class="px-2 py-0.5 bg-red-500/10 text-red-500 rounded-full text-[10px]">No encontrado</span>';
                                const posiciones = (m.posiciones_clave || []).map(pos =>
                                    `<span class="inline-block px-1.5 py-0.5 bg-slate-100 dark:bg-slate-700 rounded text-[10px] mr-1 mb-1" title="${pos.descripcion_mutacion}">pos${pos.posicion}: ${pos.aminoacido_actual}</span>`
                                ).join('');
                                return `
                                <tr class="hover:bg-slate-50 dark:hover:bg-slate-800">
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 font-mono font-bold text-primary">${m.gen}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-secondary text-[11px] max-w-[250px]">${m.descripcion || ''}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-center">${estado}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-right">${m.encontrado ? (m.longitud_observada_aa + ' / ' + m.longitud_esperada_aa) : '-'}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700">${posiciones || '-'}</td>
                                </tr>`;
                            }).join('')}
                        </tbody>
                    </table>
                </div>
            </div>

            <!-- Complejos Cuaternarios -->
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h3 class="text-sm font-semibold text-primary mb-2">Complejos Proteicos (Estructura Cuaternaria)</h3>
                <p class="text-xs text-secondary mb-4">${cuaternaria.total_complejos || 0} complejos multi-subunidad detectados. Click para expandir.</p>
                <div class="overflow-x-auto max-h-[400px] overflow-y-auto">
                    <table class="w-full text-xs">
                        <thead class="sticky top-0 bg-card">
                            <tr>
                                <th class="text-left px-2 py-2 border-b border-slate-200 text-secondary">Complejo</th>
                                <th class="text-center px-2 py-2 border-b border-slate-200 text-secondary">Subunidades</th>
                                <th class="text-right px-2 py-2 border-b border-slate-200 text-secondary">Peso Total (aa)</th>
                            </tr>
                        </thead>
                        <tbody>
                            ${(cuaternaria.complejos_detectados || []).slice(0, 25).map((c, i) => `
                                <tr class="hover:bg-slate-50 dark:hover:bg-slate-800 cursor-pointer transition" onclick="DashboardRenderer._toggleExpandRow('complejo-detail-${i}')">
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 font-medium text-primary">${(c.nombre_complejo || '').substring(0, 60)}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-center font-bold text-violet-500">${c.num_subunidades}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-right">${this.fmt(c.peso_total_estimado_aa || 0)}</td>
                                </tr>
                                <tr id="complejo-detail-${i}" class="hidden">
                                    <td colspan="3" class="px-3 py-3 bg-slate-50 dark:bg-slate-800/50 border-b border-slate-200">
                                        <p class="text-xs font-semibold text-primary mb-2">Subunidades:</p>
                                        <div class="space-y-1">
                                            ${(c.subunidades || []).map(s => `
                                                <div class="flex items-center gap-2">
                                                    <span class="font-mono text-emerald-500 text-[10px]">${s.nombre_gen || s.locus_tag}</span>
                                                    <span class="text-secondary text-[10px]">${s.producto || ''}</span>
                                                    <span class="text-slate-400 text-[10px]">(${s.longitud_aa} aa)</span>
                                                </div>
                                            `).join('')}
                                        </div>
                                    </td>
                                </tr>
                            `).join('')}
                        </tbody>
                    </table>
                </div>
            </div>

            <!-- Info Educativa - Estructura Terciaria -->
            <div class="bg-card rounded-xl p-6 border border-slate-200 mb-6">
                <h3 class="text-lg font-bold text-primary mb-3">${eduInfo.titulo || 'Estructura Terciaria'}</h3>
                <p class="text-sm text-secondary mb-4">${eduInfo.descripcion || ''}</p>

                ${eduInfo.fuerzas_estabilizadoras ? `
                <h4 class="text-sm font-semibold text-primary mb-2">Fuerzas Estabilizadoras</h4>
                <div class="grid grid-cols-1 md:grid-cols-2 gap-2 mb-4">
                    ${(eduInfo.fuerzas_estabilizadoras || []).map(f => `
                        <div class="p-3 bg-slate-50 dark:bg-slate-800 rounded-lg">
                            <p class="text-xs font-semibold text-emerald-500">${f.nombre}</p>
                            <p class="text-xs text-secondary mt-1">${f.descripcion}</p>
                        </div>
                    `).join('')}
                </div>` : ''}

                ${eduInfo.metodos_experimentales ? `
                <h4 class="text-sm font-semibold text-primary mb-2">Metodos de Determinacion Estructural</h4>
                <div class="grid grid-cols-1 md:grid-cols-2 gap-2 mb-4">
                    ${(eduInfo.metodos_experimentales || []).map(m => `
                        <div class="p-3 bg-slate-50 dark:bg-slate-800 rounded-lg">
                            <p class="text-xs font-semibold text-cyan-500">${m.nombre}</p>
                            <p class="text-xs text-secondary mt-1">${m.descripcion}</p>
                        </div>
                    `).join('')}
                </div>` : ''}

                ${eduInfo.nota_bacterias ? `
                <div class="p-3 bg-amber-50 dark:bg-amber-900/20 rounded-lg border border-amber-200 dark:border-amber-800">
                    <p class="text-xs text-amber-700 dark:text-amber-400">${eduInfo.nota_bacterias}</p>
                </div>` : ''}
            </div>

            <!-- Cisteinas y Puentes Disulfuro -->
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h3 class="text-sm font-semibold text-primary mb-2">Analisis de Cisteinas (Potenciales Puentes Disulfuro)</h3>
                <div class="grid grid-cols-3 gap-4 mb-4">
                    ${this.statsCard('Con 2+ Cys', this.fmt(cisteinas.total_con_multiples_cys || 0), 'proteinas', 'amber')}
                    ${this.statsCard('Puentes Potenciales', this.fmt(cisteinas.potenciales_puentes_disulfuro || 0), 'Cys/2 estimados', 'violet')}
                    ${this.statsCard('Top Cisteinas', (cisteinas.proteinas && cisteinas.proteinas[0]) ? cisteinas.proteinas[0].num_cisteinas + ' Cys' : '-', (cisteinas.proteinas && cisteinas.proteinas[0]) ? (cisteinas.proteinas[0].nombre_gen || cisteinas.proteinas[0].locus_tag) : '', 'cyan')}
                </div>
            </div>
        `;

        // Guardar datos de proteinas 3D para el visor
        this._proteinas3DData = proteinasConEstructura;
        this._viewer3D = null;
        this._current3DStyle = 'cartoon';

        // Renderizar charts
        setTimeout(() => {
            this._renderProteinCharts(comp, categorias, promSec, stats);
        }, 100);
    },

    // --- Charts de Proteinas ---
    _renderProteinCharts(comp, categorias, promSec, stats) {
        // 1. Composicion de aminoacidos (barras)
        const compData = comp.composicion_porcentaje || comp;
        if (compData && typeof compData === 'object') {
            const aas = Object.keys(compData).filter(k => k.length === 1).sort();
            const aaNames = {'A':'Ala','R':'Arg','N':'Asn','D':'Asp','C':'Cys','E':'Glu','Q':'Gln','G':'Gly','H':'His','I':'Ile','L':'Leu','K':'Lys','M':'Met','F':'Phe','P':'Pro','S':'Ser','T':'Thr','W':'Trp','Y':'Tyr','V':'Val'};
            const valores = aas.map(aa => compData[aa] || 0);

            this.createChart('chart-prot-aminoacidos', {
                type: 'bar',
                data: {
                    labels: aas.map(aa => aa + ' (' + (aaNames[aa] || '') + ')'),
                    datasets: [{
                        label: 'Porcentaje (%)',
                        data: valores,
                        backgroundColor: valores.map(v => v > 7 ? '#10b981' : v > 4 ? '#3b82f6' : '#94a3b8'),
                        borderRadius: 3
                    }]
                },
                options: {
                    responsive: true,
                    plugins: {
                        legend: { display: false },
                        tooltip: {
                            callbacks: {
                                title: (items) => {
                                    const aa = aas[items[0].dataIndex];
                                    const abr = aaNames[aa] || '';
                                    return abr + ' - ' + (AA_NOMBRE_COMPLETO[aa] || '');
                                }
                            }
                        }
                    },
                    scales: {
                        y: { beginAtZero: true, title: { display: true, text: '%' } },
                        x: { ticks: { font: { size: 9, family: 'monospace' } } }
                    }
                }
            });
        }

        // 2. Categorias funcionales (doughnut)
        if (categorias && typeof categorias === 'object') {
            const catNames = Object.keys(categorias).filter(k => categorias[k].total > 0);
            const catValues = catNames.map(k => categorias[k].total);
            const catColors = ['#10b981', '#3b82f6', '#f59e0b', '#8b5cf6', '#ef4444', '#06b6d4', '#22c55e', '#64748b', '#ec4899'];

            this.createChart('chart-prot-categorias', {
                type: 'doughnut',
                data: {
                    labels: catNames.map((n, i) => n + ' (' + catValues[i] + ')'),
                    datasets: [{
                        data: catValues,
                        backgroundColor: catColors.slice(0, catNames.length),
                        borderWidth: 2,
                        borderColor: '#fff'
                    }]
                },
                options: {
                    responsive: true,
                    plugins: { legend: { position: 'right', labels: { font: { size: 10 } } } }
                }
            });
        }

        // 3. Estructura secundaria (pie)
        if (promSec && promSec.helix !== undefined) {
            this.createChart('chart-prot-secundaria', {
                type: 'pie',
                data: {
                    labels: [
                        'Helice alfa (' + promSec.helix + '%)',
                        'Lamina beta (' + promSec.sheet + '%)',
                        'Giros/Coil (' + promSec.turn + '%)'
                    ],
                    datasets: [{
                        data: [promSec.helix, promSec.sheet, promSec.turn],
                        backgroundColor: ['#10b981', '#3b82f6', '#f59e0b'],
                        borderWidth: 2,
                        borderColor: '#fff'
                    }]
                },
                options: {
                    responsive: true,
                    plugins: { legend: { position: 'bottom', labels: { font: { size: 11 } } } }
                }
            });
        }

        // 4. Distribucion acidas vs basicas
        if (stats.proteinas_acidas !== undefined) {
            this.createChart('chart-prot-pi', {
                type: 'doughnut',
                data: {
                    labels: [
                        'Acidas pI<7 (' + this.fmt(stats.proteinas_acidas) + ')',
                        'Basicas pI>7 (' + this.fmt(stats.proteinas_basicas) + ')'
                    ],
                    datasets: [{
                        data: [stats.proteinas_acidas, stats.proteinas_basicas],
                        backgroundColor: ['#ef4444', '#3b82f6'],
                        borderWidth: 2,
                        borderColor: '#fff'
                    }]
                },
                options: {
                    responsive: true,
                    plugins: { legend: { position: 'bottom', labels: { font: { size: 11 } } } }
                }
            });
        }
    },

    // --- Visor 3D con 3Dmol.js ---
    _proteinas3DData: [],
    _viewer3D: null,
    _current3DStyle: 'cartoon',

    async _load3DProtein() {
        const select = document.getElementById('protein-3d-select');
        const viewerDiv = document.getElementById('protein-3d-viewer');
        const placeholder = document.getElementById('protein-3d-placeholder');
        const loadingDiv = document.getElementById('protein-3d-loading');
        const infoDiv = document.getElementById('protein-3d-info');
        const nameEl = document.getElementById('protein-3d-name');
        const detailEl = document.getElementById('protein-3d-detail');

        if (!select || select.value === '') return;

        const idx = parseInt(select.value);
        const proteinData = this._proteinas3DData[idx];
        if (!proteinData) return;

        // Mostrar loading
        if (loadingDiv) loadingDiv.classList.remove('hidden');
        if (placeholder) placeholder.style.display = 'none';
        if (viewerDiv) viewerDiv.style.display = 'block';

        try {
            // Construir URL del proxy
            let proxyUrl;
            if (proteinData.source === 'alphafold' && proteinData.url) {
                proxyUrl = `/api/proxy_pdb?pdb_id=${proteinData.pdb_id}&source=alphafold&url=${encodeURIComponent(proteinData.url)}`;
            } else {
                proxyUrl = `/api/proxy_pdb?pdb_id=${proteinData.pdb_id}&source=rcsb`;
            }

            const resp = await fetch(proxyUrl);
            if (!resp.ok) {
                throw new Error('No se pudo descargar la estructura');
            }

            const pdbData = await resp.text();

            // Detectar formato
            const format = pdbData.includes('_atom_site') ? 'cif' : 'pdb';

            // Crear o limpiar visor
            if (this._viewer3D) {
                this._viewer3D.clear();
            } else {
                if (typeof $3Dmol !== 'undefined') {
                    viewerDiv.innerHTML = '';
                    this._viewer3D = $3Dmol.createViewer(viewerDiv, {
                        backgroundColor: '#1a1a2e',
                        antialias: true
                    });
                } else {
                    throw new Error('3Dmol.js no cargado');
                }
            }

            // Cargar estructura
            this._viewer3D.addModel(pdbData, format);
            this._applyStyle3D(this._current3DStyle);
            this._viewer3D.zoomTo();
            this._viewer3D.render();
            this._viewer3D.spin('y', 0.5);

            // Mostrar info
            if (infoDiv) infoDiv.classList.remove('hidden');
            if (nameEl) nameEl.textContent = `${proteinData.nombre_gen || ''} - ${proteinData.producto || ''}`;
            if (detailEl) detailEl.textContent = `Fuente: ${proteinData.source === 'alphafold' ? 'AlphaFold (prediccion por inteligencia artificial)' : 'RCSB PDB (estructura experimental)'} | Identificador: ${proteinData.pdb_id} | Formato: ${format.toUpperCase()}`;

        } catch (err) {
            if (viewerDiv) viewerDiv.innerHTML = `<div class="flex items-center justify-center h-full"><p class="text-red-400 text-sm text-center px-4">${err.message || 'Error cargando estructura 3D'}</p></div>`;
            if (infoDiv) infoDiv.classList.add('hidden');
        }

        if (loadingDiv) loadingDiv.classList.add('hidden');
    },

    _setStyle3D(style) {
        this._current3DStyle = style;
        if (this._viewer3D) {
            this._applyStyle3D(style);
            this._viewer3D.render();
        }
    },

    _applyStyle3D(style) {
        if (!this._viewer3D) return;
        this._viewer3D.setStyle({}, {});
        const colorScheme = { prop: 'ss', map: { h: '#10b981', s: '#3b82f6', c: '#94a3b8' } };

        switch (style) {
            case 'cartoon':
                this._viewer3D.setStyle({}, { cartoon: { color: 'spectrum' } });
                break;
            case 'stick':
                this._viewer3D.setStyle({}, { stick: { colorscheme: 'Jmol', radius: 0.15 } });
                break;
            case 'sphere':
                this._viewer3D.setStyle({}, { sphere: { colorscheme: 'Jmol', scale: 0.3 } });
                break;
            case 'line':
                this._viewer3D.setStyle({}, { line: { colorscheme: 'Jmol' } });
                break;
            case 'surface':
                this._viewer3D.setStyle({}, { cartoon: { color: 'spectrum', opacity: 0.5 } });
                this._viewer3D.addSurface($3Dmol.SurfaceType.VDW, {
                    opacity: 0.7,
                    color: 'white',
                    colorscheme: { prop: 'b', gradient: new $3Dmol.Gradient.RWB(0, 100) }
                });
                break;
        }
    },

    // =========================================================================
    // VISTA DE ARCHIVOS (lista simple con descargar/eliminar)
    // =========================================================================

    renderArchivosView(tablas, figuras, genome, container) {
        const allFiles = [...(tablas || []), ...(figuras || [])];

        if (allFiles.length === 0) {
            container.innerHTML = `
                <div class="text-center py-12 text-secondary">
                    <div class="text-6xl mb-3">📄</div>
                    <p>No hay archivos de resultados para este genoma</p>
                </div>
            `;
            return;
        }

        const tipoArchivo = (file) => {
            if (file.extension === 'json') return { icon: '📋', label: 'JSON', color: 'emerald' };
            if (file.extension === 'csv') return { icon: '📊', label: 'CSV', color: 'cyan' };
            if (file.extension === 'png') return { icon: '🖼️', label: 'Imagen', color: 'violet' };
            if (file.extension === 'pdf') return { icon: '📑', label: 'PDF', color: 'red' };
            return { icon: '📄', label: 'Archivo', color: 'slate' };
        };

        container.innerHTML = `
            <!-- Resumen -->
            <div class="grid grid-cols-3 gap-4 mb-6">
                ${this.statsCard('Tablas', (tablas || []).length + ' archivos', 'JSON y CSV', 'emerald')}
                ${this.statsCard('Graficos', (figuras || []).length + ' imagenes', 'PNG generados', 'violet')}
                ${this.statsCard('Total', allFiles.length + ' archivos', this._calcTotalSize(allFiles), 'cyan')}
            </div>

            <!-- Grid de archivos -->
            <div class="grid grid-cols-1 md:grid-cols-2 gap-4">
                ${(tablas || []).map(file => {
                    const tipo = tipoArchivo(file);
                    return `
                    <div class="bg-card rounded-xl border border-slate-200 hover:border-${tipo.color}-500/30 transition overflow-hidden" id="file-card-${file.filename.replace(/\./g, '-')}">
                        <div class="px-5 py-4">
                            <div class="flex items-start justify-between mb-3">
                                <div class="flex items-center gap-3 min-w-0">
                                    <span class="text-2xl">${tipo.icon}</span>
                                    <div class="min-w-0">
                                        <p class="text-sm font-semibold text-primary truncate">${file.filename}</p>
                                        <p class="text-xs text-secondary">${file.size_kb} KB · ${tipo.label}</p>
                                    </div>
                                </div>
                                <span class="px-2 py-0.5 bg-${tipo.color}-500/10 text-${tipo.color}-500 text-xs font-medium rounded-full">${tipo.label}</span>
                            </div>

                            <!-- Preview area -->
                            <div class="bg-slate-50 dark:bg-slate-800 rounded-lg p-3 mb-3 max-h-48 overflow-y-auto">
                                <div id="preview-${file.filename.replace(/\./g, '-')}" class="text-xs font-mono text-secondary">
                                    <p class="text-center text-slate-400">Cargando vista previa...</p>
                                </div>
                            </div>

                            <!-- Acciones -->
                            <div class="flex gap-2">
                                <a href="${file.path}" download class="flex-1 px-3 py-2 bg-emerald-500 hover:bg-emerald-600 text-white text-xs font-medium rounded-lg transition text-center">
                                    Descargar
                                </a>
                                <button onclick="deleteResult('${file.filename}', 'tablas')" class="px-3 py-2 bg-red-500/10 hover:bg-red-500/20 text-red-500 text-xs font-medium rounded-lg transition">
                                    Eliminar
                                </button>
                            </div>
                        </div>
                    </div>`;
                }).join('')}

                ${(figuras || []).map(file => {
                    const tipo = tipoArchivo(file);
                    return `
                    <div class="bg-card rounded-xl border border-slate-200 hover:border-violet-500/30 transition overflow-hidden">
                        <div class="px-5 py-4">
                            <div class="flex items-start justify-between mb-3">
                                <div class="flex items-center gap-3 min-w-0">
                                    <span class="text-2xl">${tipo.icon}</span>
                                    <div class="min-w-0">
                                        <p class="text-sm font-semibold text-primary truncate">${file.filename}</p>
                                        <p class="text-xs text-secondary">${file.size_kb} KB · ${tipo.label}</p>
                                    </div>
                                </div>
                                <span class="px-2 py-0.5 bg-violet-500/10 text-violet-500 text-xs font-medium rounded-full">Imagen</span>
                            </div>

                            <!-- Preview imagen -->
                            <div class="bg-slate-50 dark:bg-slate-800 rounded-lg p-2 mb-3">
                                <img src="${file.path}" alt="${file.filename}" class="w-full rounded-lg" loading="lazy" onerror="this.parentElement.innerHTML='<p class=\\'text-xs text-center text-slate-400 py-4\\'>No se pudo cargar la imagen</p>'">
                            </div>

                            <!-- Acciones -->
                            <div class="flex gap-2">
                                <a href="${file.path}" download class="flex-1 px-3 py-2 bg-emerald-500 hover:bg-emerald-600 text-white text-xs font-medium rounded-lg transition text-center">
                                    Descargar
                                </a>
                                <button onclick="deleteResult('${file.filename}', 'figuras')" class="px-3 py-2 bg-red-500/10 hover:bg-red-500/20 text-red-500 text-xs font-medium rounded-lg transition">
                                    Eliminar
                                </button>
                            </div>
                        </div>
                    </div>`;
                }).join('')}
            </div>
        `;

        // Cargar previews de JSON y CSV
        setTimeout(() => {
            (tablas || []).forEach(file => {
                if (file.extension === 'json') {
                    this._cargarPreviewJSON(file, genome);
                } else if (file.extension === 'csv') {
                    this._cargarPreviewCSV(file, genome);
                }
            });
        }, 100);
    },

    _calcTotalSize(files) {
        const total = files.reduce((acc, f) => acc + (f.size_kb || 0), 0);
        return total > 1024 ? (total / 1024).toFixed(1) + ' MB' : total.toFixed(0) + ' KB';
    },

    async _cargarPreviewJSON(file, genome) {
        const previewId = `preview-${file.filename.replace(/\./g, '-')}`;
        const el = document.getElementById(previewId);
        if (!el) return;

        try {
            const data = await this.fetchResultData(genome, file.filename);
            if (!data) {
                el.innerHTML = '<p class="text-slate-400">Sin preview</p>';
                return;
            }

            // Mostrar primeros 10 campos clave
            const campos = Object.keys(data).slice(0, 10);
            let html = '<div class="space-y-1">';
            for (const campo of campos) {
                const valor = data[campo];
                let valorStr = '';
                if (typeof valor === 'object' && valor !== null) {
                    const subKeys = Object.keys(valor).length;
                    valorStr = `<span class="text-cyan-500">{${subKeys} campos}</span>`;
                } else if (typeof valor === 'number') {
                    valorStr = `<span class="text-emerald-500">${this.fmt(valor)}</span>`;
                } else if (typeof valor === 'string') {
                    valorStr = `<span class="text-amber-500">"${valor.substring(0, 40)}${valor.length > 40 ? '...' : ''}"</span>`;
                } else {
                    valorStr = String(valor);
                }
                html += `<p><span class="text-violet-400">${campo}:</span> ${valorStr}</p>`;
            }
            if (Object.keys(data).length > 10) {
                html += `<p class="text-slate-400 mt-1">... y ${Object.keys(data).length - 10} campos mas</p>`;
            }
            html += '</div>';
            el.innerHTML = html;
        } catch {
            el.innerHTML = '<p class="text-slate-400">Error al cargar preview</p>';
        }
    },

    async _cargarPreviewCSV(file, genome) {
        const previewId = `preview-${file.filename.replace(/\./g, '-')}`;
        const el = document.getElementById(previewId);
        if (!el) return;

        try {
            const resp = await fetch(`/api/result_data?genome=${genome}&file=${file.filename}`);
            const result = await resp.json();
            const rows = result.data || result.rows || [];
            if (!result.success || rows.length === 0) {
                el.innerHTML = '<p class="text-slate-400">Sin datos</p>';
                return;
            }

            const headers = result.headers || Object.keys(rows[0]);
            const maxRows = Math.min(rows.length, 8);
            const maxCols = Math.min(headers.length, 6);

            let html = '<div class="overflow-x-auto"><table class="w-full text-xs">';
            html += '<thead><tr>';
            for (let i = 0; i < maxCols; i++) {
                html += `<th class="text-left px-2 py-1 text-violet-400 border-b border-slate-300 dark:border-slate-600">${headers[i]}</th>`;
            }
            if (headers.length > maxCols) html += `<th class="px-2 py-1 text-slate-400 border-b border-slate-300">...</th>`;
            html += '</tr></thead><tbody>';

            for (let r = 0; r < maxRows; r++) {
                const row = rows[r];
                html += '<tr>';
                for (let i = 0; i < maxCols; i++) {
                    const val = row[headers[i]] || '';
                    const display = String(val).substring(0, 30) + (String(val).length > 30 ? '...' : '');
                    html += `<td class="px-2 py-1 text-secondary border-b border-slate-100 dark:border-slate-700">${display}</td>`;
                }
                if (headers.length > maxCols) html += '<td class="px-2 py-1 text-slate-400">...</td>';
                html += '</tr>';
            }
            html += '</tbody></table></div>';

            if (rows.length > maxRows) {
                html += `<p class="text-slate-400 text-center mt-2">... y ${rows.length - maxRows} filas mas</p>`;
            }

            el.innerHTML = html;
        } catch {
            el.innerHTML = '<p class="text-slate-400">Error al cargar preview</p>';
        }
    },

    // =========================================================================
    // DASHBOARD DE ESTRUCTURA GENOMICA
    // =========================================================================

    async renderEstructuraDashboard(data, genome, container) {
        const comp = data.composicion_genomica || {};
        const stats = data.estadisticas || {};
        const operones = data.operones_putativos || {};
        const edu = data.educativo || {};
        const gcReg = data.gc_por_region || {};

        container.innerHTML = `
            <!-- Stats Cards -->
            <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
                ${this.statsCard('Codificante (CDS)', (comp.codificante_cds?.porcentaje || 0) + '%', this.fmt(comp.codificante_cds?.bases || 0) + ' pb')}
                ${this.statsCard('Intergenico', (comp.intergenico?.porcentaje || 0) + '%', this.fmt(comp.intergenico?.bases || 0) + ' pb', 'amber')}
                ${this.statsCard('Operones Putativos', this.fmt(operones.total || 0), (operones.porcentaje_genes_en_operones || 0) + '% de genes', 'violet')}
                ${this.statsCard('GC Codificante', (gcReg.gc_codificante || 0) + '%', 'No-cod: ' + (gcReg.gc_no_codificante || 0) + '%', 'cyan')}
            </div>

            <!-- Graficos -->
            <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
                <!-- Composicion genomica -->
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Composicion del Genoma</h3>
                    <canvas id="chart-estructura-comp" height="280"></canvas>
                </div>

                <!-- Features anotados -->
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Features Anotados</h3>
                    <canvas id="chart-estructura-features" height="280"></canvas>
                </div>
            </div>

            <!-- Mapa del genoma -->
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h3 class="text-sm font-semibold text-primary mb-2">Mapa Lineal del Genoma</h3>
                <p class="text-xs text-secondary mb-3">Vista lineal del cromosoma. Genes en hebra directa (+) arriba, hebra complementaria (-) abajo. Pasa el cursor sobre un gen para ver detalles.</p>
                <div class="flex gap-2 mb-3">
                    <button onclick="DashboardRenderer._zoomGenomeMap(1.5)" class="px-3 py-1 bg-slate-100 dark:bg-slate-800 rounded text-xs hover:bg-slate-200 transition">Zoom +</button>
                    <button onclick="DashboardRenderer._zoomGenomeMap(0.67)" class="px-3 py-1 bg-slate-100 dark:bg-slate-800 rounded text-xs hover:bg-slate-200 transition">Zoom -</button>
                    <button onclick="DashboardRenderer._zoomGenomeMap(0)" class="px-3 py-1 bg-slate-100 dark:bg-slate-800 rounded text-xs hover:bg-slate-200 transition">Reiniciar</button>
                    <span class="text-xs text-secondary ml-2 self-center" id="genome-map-zoom-label">100%</span>
                </div>
                <div id="genome-map-container" class="overflow-x-auto border border-slate-100 dark:border-slate-700 rounded-lg" style="min-height: 180px;">
                    <p class="text-center text-secondary text-sm py-8">Cargando mapa del genoma...</p>
                </div>
            </div>

            <!-- Seccion educativa -->
            <div class="bg-card rounded-xl p-6 border border-slate-200 mb-6">
                <h3 class="text-lg font-bold text-primary mb-4">${edu.por_que_no_intrones?.titulo || 'Estructura del Gen Bacteriano'}</h3>
                <div class="space-y-3 text-sm text-secondary leading-relaxed">
                    ${(edu.por_que_no_intrones?.explicacion || []).map(p => `<p>${p}</p>`).join('')}
                    ${edu.por_que_no_intrones?.nota ? `<p class="text-xs italic text-amber-600 dark:text-amber-400 mt-2 p-3 bg-amber-50 dark:bg-amber-900/20 rounded-lg">${edu.por_que_no_intrones.nota}</p>` : ''}
                </div>

                <!-- Diagrama de estructura del gen bacteriano -->
                <h4 class="text-sm font-semibold text-primary mt-6 mb-3">${edu.estructura_gen_bacteriano?.titulo || 'Estructura del Gen'}</h4>
                <div class="flex items-center gap-1 overflow-x-auto pb-2">
                    ${(edu.estructura_gen_bacteriano?.elementos || []).map(el => `
                        <div class="flex flex-col items-center min-w-[100px]">
                            <div class="w-full h-12 rounded-lg flex items-center justify-center text-white text-xs font-bold shadow-md" style="background: ${el.color}">${el.nombre}</div>
                            <p class="text-[10px] text-secondary mt-1 text-center">${el.posicion}</p>
                        </div>
                        <div class="text-slate-300 text-lg flex-shrink-0">→</div>
                    `).join('').replace(/→<\/div>$/, '</div>')}
                </div>
                <div class="mt-3 grid grid-cols-1 md:grid-cols-2 gap-2">
                    ${(edu.estructura_gen_bacteriano?.elementos || []).map(el => `
                        <div class="flex items-start gap-2 p-2 rounded-lg bg-slate-50 dark:bg-slate-800">
                            <div class="w-3 h-3 rounded-sm flex-shrink-0 mt-0.5" style="background: ${el.color}"></div>
                            <div>
                                <span class="text-xs font-semibold text-primary">${el.nombre}:</span>
                                <span class="text-xs text-secondary ml-1">${el.descripcion}</span>
                            </div>
                        </div>
                    `).join('')}
                </div>
            </div>

            <!-- Comparacion eucariota vs procariota -->
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h3 class="text-sm font-semibold text-primary mb-4">${edu.comparacion_eucariota?.titulo || 'Comparacion'}</h3>
                <div class="overflow-x-auto">
                    <table class="w-full text-sm">
                        <thead>
                            <tr>
                                <th class="text-left px-3 py-2 border-b border-slate-200 text-secondary">Aspecto</th>
                                <th class="text-left px-3 py-2 border-b border-slate-200 text-emerald-500">Bacteria</th>
                                <th class="text-left px-3 py-2 border-b border-slate-200 text-violet-500">Eucariota</th>
                            </tr>
                        </thead>
                        <tbody>
                            ${(edu.comparacion_eucariota?.diferencias || []).map(d => `
                                <tr>
                                    <td class="px-3 py-2 border-b border-slate-100 dark:border-slate-700 font-medium text-primary">${d.aspecto}</td>
                                    <td class="px-3 py-2 border-b border-slate-100 dark:border-slate-700 text-secondary">${d.bacteria}</td>
                                    <td class="px-3 py-2 border-b border-slate-100 dark:border-slate-700 text-secondary">${d.eucariota}</td>
                                </tr>
                            `).join('')}
                        </tbody>
                    </table>
                </div>
            </div>

            <!-- Tabla de operones -->
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h3 class="text-sm font-semibold text-primary mb-4">Operones Putativos (Top 20)</h3>
                <p class="text-xs text-secondary mb-3">Grupos de genes adyacentes en la misma hebra con distancia intergenica < 50 pb. Haz click para ver todos los genes.</p>
                <div class="overflow-x-auto max-h-[500px] overflow-y-auto">
                    <table class="w-full text-xs">
                        <thead class="sticky top-0 bg-card">
                            <tr>
                                <th class="text-left px-2 py-2 border-b border-slate-200 text-secondary">#</th>
                                <th class="text-left px-2 py-2 border-b border-slate-200 text-secondary">Genes</th>
                                <th class="text-center px-2 py-2 border-b border-slate-200 text-secondary">N</th>
                                <th class="text-center px-2 py-2 border-b border-slate-200 text-secondary">Hebra</th>
                                <th class="text-right px-2 py-2 border-b border-slate-200 text-secondary">Inicio</th>
                                <th class="text-right px-2 py-2 border-b border-slate-200 text-secondary">Fin</th>
                                <th class="text-right px-2 py-2 border-b border-slate-200 text-secondary">Longitud</th>
                            </tr>
                        </thead>
                        <tbody>
                            ${(operones.operones || []).sort((a, b) => b.num_genes - a.num_genes).slice(0, 20).map((op, idx) => `
                                <tr class="hover:bg-slate-50 dark:hover:bg-slate-800 cursor-pointer transition" onclick="DashboardRenderer._toggleExpandRow('operon-detail-${idx}')">
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700">${op.numero}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 font-mono text-primary">${op.genes.slice(0, 6).join(', ')}${op.genes.length > 6 ? ' ...' : ''}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-center font-bold text-emerald-500">${op.num_genes}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-center">${op.hebra}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-right">${this.fmt(op.inicio)}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-right">${this.fmt(op.fin)}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-right">${this.fmt(op.longitud_total)} pb</td>
                                </tr>
                                <tr id="operon-detail-${idx}" class="hidden">
                                    <td colspan="7" class="px-3 py-3 bg-slate-50 dark:bg-slate-800/50 border-b border-slate-200">
                                        <p class="text-xs font-semibold text-primary mb-2">Todos los genes del operon ${op.numero}:</p>
                                        <div class="flex flex-wrap gap-1">
                                            ${op.genes.map(g => `<span class="px-2 py-1 bg-emerald-100 dark:bg-emerald-900/30 text-emerald-700 dark:text-emerald-400 rounded font-mono text-[10px]">${g}</span>`).join('')}
                                        </div>
                                        <p class="text-xs text-secondary mt-2">Posicion: ${this.fmt(op.inicio)} - ${this.fmt(op.fin)} | Hebra: ${op.hebra === '+' ? 'Directa' : 'Complementaria'} | ${op.num_genes} genes en ${this.fmt(op.longitud_total)} pb</p>
                                    </td>
                                </tr>
                            `).join('')}
                        </tbody>
                    </table>
                </div>
            </div>
        `;

        // Renderizar charts
        this._renderChartEstructuraComp(comp);
        this._renderChartEstructuraFeatures(data.features_anotados || {});

        // Cargar mapa del genoma (necesita CSV de genes)
        this._cargarMapaGenoma(genome, data.longitud_genoma);
    },

    _renderChartEstructuraComp(comp) {
        const labels = [];
        const values = [];
        const colors = ['#3b82f6', '#f59e0b', '#8b5cf6', '#ef4444', '#06b6d4', '#6b7280'];

        const tipos = [
            ['codificante_cds', 'CDS'],
            ['intergenico', 'Intergenico'],
            ['rRNA', 'rRNA'],
            ['tRNA', 'tRNA'],
            ['otro_RNA_funcional', 'Otro RNA'],
            ['regiones_repetitivas', 'Repetitivo']
        ];

        for (const [key, label] of tipos) {
            if (comp[key] && comp[key].porcentaje > 0) {
                labels.push(label + ' (' + comp[key].porcentaje + '%)');
                values.push(comp[key].porcentaje);
            }
        }

        this.createChart('chart-estructura-comp', {
            type: 'doughnut',
            data: {
                labels,
                datasets: [{
                    data: values,
                    backgroundColor: colors.slice(0, values.length),
                    borderWidth: 2,
                    borderColor: '#fff'
                }]
            },
            options: {
                responsive: true,
                plugins: {
                    legend: { position: 'bottom', labels: { font: { size: 11 } } }
                }
            }
        });
    },

    _renderChartEstructuraFeatures(features) {
        const exclude = ['source', 'gene'];
        const labels = [];
        const values = [];
        const colors = ['#3b82f6', '#10b981', '#f59e0b', '#8b5cf6', '#ef4444', '#06b6d4', '#ec4899', '#6b7280'];

        for (const [tipo, count] of Object.entries(features)) {
            if (!exclude.includes(tipo) && count > 0) {
                labels.push(tipo);
                values.push(count);
            }
        }

        // Ordenar por count descendente
        const sorted = labels.map((l, i) => ({ l, v: values[i] })).sort((a, b) => b.v - a.v).slice(0, 8);

        this.createChart('chart-estructura-features', {
            type: 'bar',
            data: {
                labels: sorted.map(x => x.l),
                datasets: [{
                    label: 'Cantidad',
                    data: sorted.map(x => x.v),
                    backgroundColor: colors.slice(0, sorted.length)
                }]
            },
            options: {
                responsive: true,
                indexAxis: 'y',
                plugins: { legend: { display: false } },
                scales: {
                    x: { title: { display: true, text: 'Cantidad' } }
                }
            }
        });
    },

    _genomeMapZoom: 1,
    _lastGenomeMapGenes: null,
    _lastGenomeMapLength: null,

    async _cargarMapaGenoma(genome, genomeLength) {
        const mapContainer = document.getElementById('genome-map-container');
        if (!mapContainer) return;

        try {
            const csvFile = `lista_genes_${genome}.csv`;
            const resp = await fetch(`/api/result_data?genome=${genome}&file=${csvFile}`);
            const result = await resp.json();

            if (!result.success || !result.data || result.data.length === 0) {
                mapContainer.innerHTML = '<p class="text-center text-secondary text-sm py-4">No hay datos de genes para el mapa. Ejecuta el analisis de genes primero.</p>';
                return;
            }

            const genes = result.data;
            this._genomeMapZoom = 1;
            this._lastGenomeMapGenes = genes;
            this._lastGenomeMapLength = genomeLength;
            this._renderGenomeMapSVG(genes, genomeLength, mapContainer);
        } catch (e) {
            mapContainer.innerHTML = '<p class="text-center text-red-400 text-sm py-4">Error cargando mapa del genoma</p>';
        }
    },

    _renderGenomeMapSVG(genes, genomeLength, container) {
        const width = Math.max(container.clientWidth - 20, 800) * this._genomeMapZoom;
        const height = 160;
        const margin = 30;
        const trackFwd = 35;
        const trackCenter = 75;
        const trackRev = 115;
        const geneH = 28;

        let svg = `<svg width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" xmlns="http://www.w3.org/2000/svg" style="min-width: ${width}px">`;

        // Fondo
        svg += `<rect width="${width}" height="${height}" fill="transparent"/>`;

        // Linea central del cromosoma
        svg += `<line x1="${margin}" y1="${trackCenter}" x2="${width - margin}" y2="${trackCenter}" stroke="#475569" stroke-width="2" stroke-dasharray="4,2"/>`;

        // Marcas de escala
        const scaleStep = Math.pow(10, Math.floor(Math.log10(genomeLength / 10)));
        for (let pos = 0; pos <= genomeLength; pos += scaleStep) {
            const x = margin + (pos / genomeLength) * (width - 2 * margin);
            svg += `<line x1="${x}" y1="${trackCenter - 5}" x2="${x}" y2="${trackCenter + 5}" stroke="#64748b" stroke-width="1"/>`;
            if (pos % (scaleStep * 2) === 0) {
                const label = pos >= 1000000 ? (pos / 1000000).toFixed(1) + ' Mb' : pos >= 1000 ? (pos / 1000).toFixed(0) + ' kb' : pos;
                svg += `<text x="${x}" y="${height - 5}" text-anchor="middle" fill="#94a3b8" font-size="9">${label}</text>`;
            }
        }

        // Genes como rectangulos
        for (const gene of genes) {
            const inicio = parseInt(gene.inicio || gene.start || 0);
            const fin = parseInt(gene.fin || gene.end || 0);
            const hebra = gene.hebra || gene.strand || '+';
            const nombre = gene.nombre_gen || gene.gene || gene.locus_tag || '';

            const x = margin + (inicio / genomeLength) * (width - 2 * margin);
            const w = Math.max(1, ((fin - inicio) / genomeLength) * (width - 2 * margin));
            const y = hebra === '+' || hebra === '1' ? trackFwd : trackRev;
            const color = hebra === '+' || hebra === '1' ? '#10b981' : '#06b6d4';

            svg += `<rect x="${x}" y="${y}" width="${w}" height="${geneH}" fill="${color}" opacity="0.6" rx="1">`;
            svg += `<title>${nombre}\n${gene.producto || ''}\nPos: ${this.fmt(inicio)}-${this.fmt(fin)} (${hebra})\n${this.fmt(fin - inicio)} pb</title>`;
            svg += `</rect>`;
        }

        // Labels de hebras
        svg += `<text x="${margin}" y="25" fill="#10b981" font-size="11" font-weight="bold">Directa (+) - ${genes.filter(g => (g.hebra || g.strand) === '+' || (g.hebra || g.strand) === '1').length} genes</text>`;
        svg += `<text x="${margin}" y="${height - 15}" fill="#06b6d4" font-size="11" font-weight="bold">Complementaria (-) - ${genes.filter(g => (g.hebra || g.strand) === '-' || (g.hebra || g.strand) === '-1').length} genes</text>`;

        svg += '</svg>';
        container.innerHTML = svg;
    },

    _zoomGenomeMap(factor) {
        const container = document.getElementById('genome-map-container');
        const label = document.getElementById('genome-map-zoom-label');
        if (!container) return;

        if (factor === 0) {
            this._genomeMapZoom = 1;
        } else {
            this._genomeMapZoom = Math.max(0.5, Math.min(10, this._genomeMapZoom * factor));
        }

        if (label) label.textContent = Math.round(this._genomeMapZoom * 100) + '%';

        // Re-render usando los datos almacenados para que todo se recalule con el nuevo zoom
        if (this._lastGenomeMapGenes && this._lastGenomeMapLength) {
            this._renderGenomeMapSVG(this._lastGenomeMapGenes, this._lastGenomeMapLength, container);
        }
    },

    // =========================================================================
    // BUSQUEDA DE SECUENCIA
    // =========================================================================

    renderBusquedaSecuencia(genome, container) {
        container.innerHTML = `
            <div class="max-w-4xl mx-auto">
                <!-- Input de busqueda -->
                <div class="bg-card rounded-xl p-6 border border-slate-200 mb-6">
                    <h3 class="text-lg font-bold text-primary mb-2">Localizar Secuencia en el Genoma</h3>
                    <p class="text-sm text-secondary mb-4">Busca una secuencia de nucleotidos y encuentra su posicion exacta, en que gen cae y en que hebra se encuentra.</p>

                    <div class="flex gap-3">
                        <input
                            type="text"
                            id="buscar-secuencia-input"
                            placeholder="Ingresa secuencia (ej: ATGCGATCGA, min 4 nucleotidos)"
                            class="flex-1 px-4 py-3 bg-slate-50 dark:bg-slate-800 rounded-xl border border-slate-200 dark:border-slate-700 focus:border-emerald-500 focus:ring-2 focus:ring-emerald-500/20 outline-none text-sm text-primary font-mono uppercase"
                            onkeydown="if(event.key==='Enter') DashboardRenderer._buscarSecuencia('${genome}')"
                        >
                        <button
                            onclick="DashboardRenderer._buscarSecuencia('${genome}')"
                            class="px-6 py-3 bg-gradient-to-r from-emerald-500 to-cyan-500 hover:from-emerald-600 hover:to-cyan-600 text-white rounded-xl font-medium text-sm transition shadow-lg"
                        >
                            Buscar
                        </button>
                    </div>
                    <p class="text-xs text-secondary mt-2">Se busca en ambas hebras (directa y complementaria inversa). Maximo 200 resultados.</p>
                </div>

                <!-- Resultados -->
                <div id="buscar-secuencia-resultados"></div>
            </div>
        `;
    },

    async _buscarSecuencia(genome) {
        const input = document.getElementById('buscar-secuencia-input');
        const resultContainer = document.getElementById('buscar-secuencia-resultados');
        const secuencia = (input?.value || '').trim().toUpperCase().replace(/[^ATGCN]/g, '');

        if (secuencia.length < 4) {
            resultContainer.innerHTML = '<p class="text-amber-500 text-sm text-center py-4">La secuencia debe tener al menos 4 nucleotidos (A, T, G, C)</p>';
            return;
        }

        resultContainer.innerHTML = `
            <div class="text-center py-8">
                <div class="w-8 h-8 border-2 border-emerald-500 border-t-transparent rounded-full animate-spin mx-auto mb-3"></div>
                <p class="text-secondary text-sm">Buscando "${secuencia}" en el genoma...</p>
            </div>
        `;

        try {
            const resp = await fetch(`/api/buscar_secuencia?genome=${genome}&secuencia=${secuencia}`);
            const data = await resp.json();

            if (!data.success) {
                resultContainer.innerHTML = `<p class="text-red-500 text-sm text-center py-4">${data.error}</p>`;
                return;
            }

            if (data.total_matches === 0) {
                resultContainer.innerHTML = `
                    <div class="bg-card rounded-xl p-6 border border-slate-200 text-center">
                        <div class="text-4xl mb-2">🔍</div>
                        <p class="text-primary font-medium">No se encontro la secuencia</p>
                        <p class="text-sm text-secondary mt-1">"${secuencia}" no fue encontrada en ninguna hebra del genoma (${this.fmt(data.longitud_genoma)} pb)</p>
                    </div>
                `;
                return;
            }

            let html = `
                <div class="bg-card rounded-xl p-5 border border-slate-200 mb-4">
                    <div class="flex items-center justify-between mb-3">
                        <h4 class="text-sm font-semibold text-primary">Resultados de Busqueda</h4>
                        <span class="px-3 py-1 bg-emerald-500/10 text-emerald-500 rounded-full text-xs font-medium">
                            ${data.total_matches} coincidencia${data.total_matches > 1 ? 's' : ''}
                        </span>
                    </div>
                    <p class="text-xs text-secondary mb-3">Secuencia: <span class="font-mono text-primary">${secuencia}</span> (${secuencia.length} nt) en genoma de ${this.fmt(data.longitud_genoma)} pb</p>

                    <div class="overflow-x-auto max-h-[500px] overflow-y-auto">
                        <table class="w-full text-xs">
                            <thead class="sticky top-0 bg-card">
                                <tr>
                                    <th class="text-left px-2 py-2 border-b border-slate-200 text-secondary">#</th>
                                    <th class="text-left px-2 py-2 border-b border-slate-200 text-secondary">Posicion</th>
                                    <th class="text-center px-2 py-2 border-b border-slate-200 text-secondary">Hebra</th>
                                    <th class="text-left px-2 py-2 border-b border-slate-200 text-secondary">Gen</th>
                                    <th class="text-left px-2 py-2 border-b border-slate-200 text-secondary">Producto</th>
                                    <th class="text-left px-2 py-2 border-b border-slate-200 text-secondary">Contexto</th>
                                </tr>
                            </thead>
                            <tbody>
            `;

            for (let i = 0; i < data.matches.length; i++) {
                const m = data.matches[i];
                const genNombre = m.gen?.nombre || 'Intergenico';
                const genProducto = m.gen?.producto || '';
                const hebraColor = m.hebra === '+' ? 'text-emerald-500' : 'text-cyan-500';

                // Resaltar match en contexto
                let contexto = m.contexto || '';
                const matchIdx = contexto.toUpperCase().indexOf(secuencia);
                if (matchIdx >= 0) {
                    contexto = contexto.substring(0, matchIdx) +
                        '<span class="bg-yellow-200 dark:bg-yellow-800 font-bold">' +
                        contexto.substring(matchIdx, matchIdx + secuencia.length) +
                        '</span>' + contexto.substring(matchIdx + secuencia.length);
                }

                html += `
                    <tr class="hover:bg-slate-50 dark:hover:bg-slate-800">
                        <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700">${i + 1}</td>
                        <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 font-mono">${this.fmt(m.posicion)}-${this.fmt(m.fin)}</td>
                        <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-center font-bold ${hebraColor}">${m.hebra}</td>
                        <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 font-medium text-primary">${genNombre}</td>
                        <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-secondary max-w-[200px] truncate">${genProducto}</td>
                        <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 font-mono text-[10px]">${contexto}</td>
                    </tr>
                `;
            }

            html += '</tbody></table></div></div>';
            resultContainer.innerHTML = html;
        } catch (e) {
            resultContainer.innerHTML = `<p class="text-red-500 text-sm text-center py-4">Error de conexion: ${e.message}</p>`;
        }
    },

    // =========================================================================
    // EXTRACTOR DE SECUENCIAS
    // =========================================================================

    renderExtractorSecuencia(genome, container) {
        // Obtener total de genes del JSON de genes
        this.fetchResultData(genome, `analisis_genes_${genome}.json`).then(genesData => {
            const totalGenes = genesData?.estadisticas_generales?.total_genes || '?';
            container.innerHTML = `
                <div class="max-w-5xl mx-auto">
                    <div class="bg-card rounded-xl p-6 border border-slate-200 mb-6">
                        <h3 class="text-lg font-bold text-primary mb-2">Extractor de Secuencias Genicas</h3>
                        <p class="text-sm text-secondary mb-4">Extrae secuencias de ADN de un rango de genes. Los genes estan ordenados por posicion en el genoma.</p>

                        <div class="grid grid-cols-1 md:grid-cols-4 gap-4 mb-4">
                            <div>
                                <label class="block text-xs font-medium text-secondary mb-1">Gen Inicio (1-based)</label>
                                <input type="number" id="extractor-start" min="1" max="${totalGenes}" value="1"
                                    class="w-full px-3 py-2 bg-slate-50 dark:bg-slate-800 rounded-lg border border-slate-200 dark:border-slate-700 text-sm text-primary focus:border-emerald-500 outline-none">
                            </div>
                            <div>
                                <label class="block text-xs font-medium text-secondary mb-1">Gen Fin (1-based)</label>
                                <input type="number" id="extractor-end" min="1" max="${totalGenes}" value="5"
                                    class="w-full px-3 py-2 bg-slate-50 dark:bg-slate-800 rounded-lg border border-slate-200 dark:border-slate-700 text-sm text-primary focus:border-emerald-500 outline-none">
                            </div>
                            <div>
                                <label class="block text-xs font-medium text-secondary mb-1">Modo</label>
                                <select id="extractor-mode"
                                    class="w-full px-3 py-2 bg-slate-50 dark:bg-slate-800 rounded-lg border border-slate-200 dark:border-slate-700 text-sm text-primary focus:border-emerald-500 outline-none">
                                    <option value="cds">Solo CDS (secuencias codificantes)</option>
                                    <option value="todo">Region completa (con ADN intergenico)</option>
                                </select>
                            </div>
                            <div class="flex items-end">
                                <button onclick="DashboardRenderer._ejecutarExtraccion('${genome}')"
                                    class="w-full px-4 py-2 bg-gradient-to-r from-emerald-500 to-cyan-500 hover:from-emerald-600 hover:to-cyan-600 text-white rounded-lg font-medium text-sm transition shadow-lg">
                                    Extraer
                                </button>
                            </div>
                        </div>
                        <p class="text-xs text-secondary">Total de genes en el genoma: <strong class="text-primary">${totalGenes}</strong>. E. coli no tiene intrones: "Solo CDS" = secuencias codificantes concatenadas, "Region completa" = incluye ADN intergenico.</p>
                    </div>

                    <div id="extractor-resultados"></div>
                </div>
            `;
        });
    },

    async _ejecutarExtraccion(genome) {
        const startInput = document.getElementById('extractor-start');
        const endInput = document.getElementById('extractor-end');
        const modeSelect = document.getElementById('extractor-mode');
        const resultContainer = document.getElementById('extractor-resultados');

        const geneStart = parseInt(startInput?.value) || 1;
        const geneEnd = parseInt(endInput?.value) || 5;
        const mode = modeSelect?.value || 'cds';

        if (geneStart < 1 || geneEnd < geneStart) {
            resultContainer.innerHTML = '<p class="text-amber-500 text-sm text-center py-4">Rango invalido: Gen inicio debe ser >= 1 y Gen fin >= Gen inicio</p>';
            return;
        }

        resultContainer.innerHTML = `
            <div class="text-center py-8">
                <div class="w-8 h-8 border-2 border-emerald-500 border-t-transparent rounded-full animate-spin mx-auto mb-3"></div>
                <p class="text-secondary text-sm">Extrayendo genes ${geneStart} a ${geneEnd} (modo: ${mode})...</p>
            </div>
        `;

        try {
            const resp = await fetch(`/api/extraer_secuencia?genome=${genome}&gene_start=${geneStart}&gene_end=${geneEnd}&mode=${mode}`);
            const data = await resp.json();

            if (!data.success) {
                resultContainer.innerHTML = `<p class="text-red-500 text-sm text-center py-4">${data.error}</p>`;
                return;
            }

            const modoTexto = mode === 'cds' ? 'Solo CDS (codificantes)' : 'Region completa (con intergenico)';

            // Formatear secuencia en formato FASTA (60 chars por linea con numeros de posicion)
            const seq = data.secuencia_total || '';
            let fastaLines = `>${genome}_genes_${geneStart}-${geneEnd}_${mode} | ${data.num_genes} genes | ${this.fmt(data.longitud_secuencia)} pb\n`;
            for (let i = 0; i < seq.length; i += 60) {
                fastaLines += seq.substring(i, i + 60) + '\n';
            }

            let html = `
                <!-- Stats -->
                <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
                    ${this.statsCard('Genes Extraidos', data.num_genes, 'de ' + data.total_genes_genome + ' totales')}
                    ${this.statsCard('Longitud', this.fmt(data.longitud_secuencia) + ' pb', modoTexto, 'cyan')}
                    ${this.statsCard('Rango', geneStart + ' - ' + geneEnd, 'posicion en genoma', 'amber')}
                    ${this.statsCard('Modo', mode === 'cds' ? 'CDS' : 'Completo', mode === 'cds' ? 'Solo codificante' : 'Con intergenico', 'violet')}
                </div>

                <!-- Tabla de genes -->
                <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                    <h4 class="text-sm font-semibold text-primary mb-3">Genes en el Rango</h4>
                    <div class="overflow-x-auto max-h-[400px] overflow-y-auto">
                        <table class="w-full text-xs">
                            <thead class="bg-slate-50 dark:bg-slate-800 sticky top-0">
                                <tr>
                                    <th class="px-2 py-2 text-left text-secondary">#</th>
                                    <th class="px-2 py-2 text-left text-secondary">Gen</th>
                                    <th class="px-2 py-2 text-left text-secondary">Locus</th>
                                    <th class="px-2 py-2 text-left text-secondary">Producto</th>
                                    <th class="px-2 py-2 text-right text-secondary">Inicio</th>
                                    <th class="px-2 py-2 text-right text-secondary">Fin</th>
                                    <th class="px-2 py-2 text-center text-secondary">Hebra</th>
                                    <th class="px-2 py-2 text-right text-secondary">Longitud</th>
                                </tr>
                            </thead>
                            <tbody>
            `;

            for (const g of data.genes_info) {
                html += `
                    <tr class="border-t border-slate-100 dark:border-slate-800 hover:bg-slate-50 dark:hover:bg-slate-800/50">
                        <td class="px-2 py-2 text-secondary">${g.index}</td>
                        <td class="px-2 py-2 font-mono font-bold text-primary">${g.nombre_gen || '-'}</td>
                        <td class="px-2 py-2 font-mono text-secondary">${g.locus_tag}</td>
                        <td class="px-2 py-2 text-secondary max-w-[250px] truncate">${g.producto || 'Sin anotacion'}</td>
                        <td class="px-2 py-2 text-right text-primary">${this.fmt(g.inicio)}</td>
                        <td class="px-2 py-2 text-right text-primary">${this.fmt(g.fin)}</td>
                        <td class="px-2 py-2 text-center font-bold ${g.hebra === '+' ? 'text-emerald-500' : 'text-cyan-500'}">${g.hebra}</td>
                        <td class="px-2 py-2 text-right font-medium text-primary">${this.fmt(g.longitud_pb)} pb</td>
                    </tr>
                `;
            }

            html += `
                            </tbody>
                        </table>
                    </div>
                </div>

                <!-- Secuencia FASTA -->
                <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                    <div class="flex items-center justify-between mb-3">
                        <h4 class="text-sm font-semibold text-primary">Secuencia (formato FASTA)</h4>
                        <div class="flex gap-2">
                            <button onclick="DashboardRenderer._copiarSecuencia()" class="px-3 py-1.5 bg-emerald-500 hover:bg-emerald-600 text-white text-xs rounded-lg transition">
                                Copiar al Portapapeles
                            </button>
                            <button onclick="DashboardRenderer._descargarFASTA('${genome}', ${geneStart}, ${geneEnd}, '${mode}')" class="px-3 py-1.5 bg-cyan-500 hover:bg-cyan-600 text-white text-xs rounded-lg transition">
                                Descargar FASTA
                            </button>
                        </div>
                    </div>
                    <pre id="extractor-secuencia-pre" class="bg-slate-900 text-emerald-400 p-4 rounded-lg text-xs font-mono overflow-x-auto max-h-[400px] overflow-y-auto leading-relaxed select-all">${fastaLines}</pre>
                    <p class="text-xs text-secondary mt-2">${this.fmt(data.longitud_secuencia)} pares de bases | ${data.num_genes} genes | Modo: ${modoTexto}</p>
                </div>
            `;

            resultContainer.innerHTML = html;

            // Guardar secuencia para copiar
            DashboardRenderer._extractedFasta = fastaLines;

        } catch (e) {
            resultContainer.innerHTML = `<p class="text-red-500 text-sm text-center py-4">Error de conexion: ${e.message}</p>`;
        }
    },

    _extractedFasta: '',

    async _copiarSecuencia() {
        try {
            await navigator.clipboard.writeText(DashboardRenderer._extractedFasta);
            const btn = event.target;
            const originalText = btn.textContent;
            btn.textContent = 'Copiado!';
            btn.classList.replace('bg-emerald-500', 'bg-green-600');
            setTimeout(() => {
                btn.textContent = originalText;
                btn.classList.replace('bg-green-600', 'bg-emerald-500');
            }, 2000);
        } catch (e) {
            alert('Error al copiar: ' + e.message);
        }
    },

    _descargarFASTA(genome, start, end, mode) {
        const blob = new Blob([DashboardRenderer._extractedFasta], { type: 'text/plain' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `${genome}_genes_${start}-${end}_${mode}.fasta`;
        a.click();
        URL.revokeObjectURL(url);
    },

    // =========================================================================
    // UTILIDADES
    // =========================================================================

    /**
     * Toggle expand/collapse de una fila de detalle
     */
    _toggleExpandRow(rowId) {
        const row = document.getElementById(rowId);
        if (!row) return;
        if (row.classList.contains('hidden')) {
            row.classList.remove('hidden');
            row.style.animation = 'fadeIn 0.2s ease-out';
        } else {
            row.classList.add('hidden');
        }
    },

    // Cache de datos cargados
    _cache: {},

    async fetchResultData(genome, filename) {
        const key = `${genome}/${filename}`;
        if (this._cache[key]) return this._cache[key];

        try {
            const resp = await fetch(`/api/result_data?genome=${genome}&file=${filename}`);
            const result = await resp.json();
            if (result.success) {
                this._cache[key] = result.data;
                return result.data;
            }
        } catch (e) {
            console.error('Error fetching result data:', e);
        }
        return null;
    },

    clearCache() {
        this._cache = {};
    },

    downloadJSON(type) {
        const genome = document.getElementById('results-genome')?.value;
        if (!genome) return;

        const fileMap = {
            genes: `analisis_genes_${genome}.json`,
            codones: `analisis_codones_${genome}.json`,
            distancias: `analisis_distancias_${genome}.json`
        };

        const filename = fileMap[type];
        if (filename) {
            const a = document.createElement('a');
            a.href = `backend/resultados/${genome}/tablas/${filename}`;
            a.download = filename;
            a.click();
        }
    }
};
