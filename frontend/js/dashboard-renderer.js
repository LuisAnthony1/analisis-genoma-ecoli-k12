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
     * Etiqueta de hebra biologica
     */
    hebraLabel(hebra, format = 'short') {
        if (hebra === '+' || hebra === '1' || hebra === 1) {
            return format === 'full' ? "Hebra 5'\u21923' (+)" : "5'\u21923'";
        }
        return format === 'full' ? "Hebra 3'\u21925' (-)" : "3'\u21925'";
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
                    <h3 class="text-sm font-semibold text-primary mb-1">Distribucion por Hebra</h3>
                    <p class="text-[10px] text-secondary mb-3"><span class="inline-block w-2.5 h-2.5 rounded-full mr-1" style="background:#10b981"></span>5'&rarr;3' = hebra directa (+) &nbsp; <span class="inline-block w-2.5 h-2.5 rounded-full mr-1" style="background:#06b6d4"></span>3'&rarr;5' = hebra complementaria (-)</p>
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
                                    <td class="px-3 py-2 text-center ${g.hebra === '+' ? 'text-emerald-500' : 'text-cyan-500'} font-bold">${g.hebra === '+' ? "5'\u21923'" : "3'\u21925'"}</td>
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
                                    <td class="px-3 py-2 text-center ${g.hebra === '+' ? 'text-emerald-500' : 'text-cyan-500'} font-bold">${g.hebra === '+' ? "5'\u21923'" : "3'\u21925'"}</td>
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
                        labels: ["Hebra 5'\u21923' (+)", "Hebra 3'\u21925' (-)"],
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
                    <h3 class="text-sm font-semibold text-primary mb-1">Analisis por Hebra</h3>
                    <p class="text-[10px] text-secondary mb-3">Misma hebra = ambos genes en 5'&rarr;3' o ambos en 3'&rarr;5'. Diferente hebra = uno en cada direccion.</p>
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
                                    <td class="px-3 py-2 text-center text-secondary">${r.hebra_gen1 === '+' ? "5'\u21923'" : "3'\u21925'"} &rarr; ${r.hebra_gen2 === '+' ? "5'\u21923'" : "3'\u21925'"}</td>
                                </tr>
                                <tr id="dist-detail-${i}" class="hidden">
                                    <td colspan="7" class="px-4 py-3 bg-slate-50 dark:bg-slate-800/50 border-b border-slate-200">
                                        <div class="grid grid-cols-1 md:grid-cols-2 gap-4 text-xs">
                                            <div class="p-3 bg-white dark:bg-slate-900 rounded-lg border border-slate-200">
                                                <p class="font-bold text-emerald-500 mb-1">Gen 1: ${r.gen1}</p>
                                                <p><span class="text-secondary">Proteina:</span> <span class="text-primary font-medium">${r.gen1_producto || 'No anotado'}</span></p>
                                                <p><span class="text-secondary">Hebra:</span> ${r.hebra_gen1 === '+' ? "5'\u21923' (+)" : "3'\u21925' (-)"}</p>
                                            </div>
                                            <div class="p-3 bg-white dark:bg-slate-900 rounded-lg border border-slate-200">
                                                <p class="font-bold text-cyan-500 mb-1">Gen 2: ${r.gen2}</p>
                                                <p><span class="text-secondary">Proteina:</span> <span class="text-primary font-medium">${r.gen2_producto || 'No anotado'}</span></p>
                                                <p><span class="text-secondary">Hebra:</span> ${r.hebra_gen2 === '+' ? "5'\u21923' (+)" : "3'\u21925' (-)"}</p>
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
        const ortologos = data.ortologos || {};
        const mutaciones = data.mutaciones || {};
        const sintenia = data.sintenia || {};
        const islas = data.islas_genomicas || {};
        const funcional = data.categorias_funcionales || {};
        const trnaRrna = data.trna_rrna || {};
        const resistencia = data.resistencia_antibiotica || {};
        const gcRegional = data.gc_regional || {};
        const aniData = data.ani || {};

        const ec = metricas.ecoli || {};
        const sal = metricas.salmonella || {};
        const nombre1 = ec.nombre || organismos.organismo_1?.nombre || 'Genoma 1';
        const nombre2 = sal.nombre || organismos.organismo_2?.nombre || 'Genoma 2';
        const n1 = nombre1.length > 25 ? nombre1.substring(0, 25) + '...' : nombre1;
        const n2 = nombre2.length > 25 ? nombre2.substring(0, 25) + '...' : nombre2;
        const ortoStats = ortologos.estadisticas || {};
        const mutStats = mutaciones.estadisticas || {};
        const sinStats = sintenia.estadisticas || {};
        const islStats = islas.estadisticas || {};
        const rnaStats = trnaRrna.estadisticas || {};
        const resStats = resistencia.estadisticas || {};
        const gcStats = gcRegional.estadisticas || {};

        // Helper for ANI badge color
        const ani = aniData.ani_estimado || mutStats.ani_ortologos || 0;
        const aniColor = ani >= 96 ? 'emerald' : ani >= 90 ? 'amber' : 'red';
        const kaKs = mutStats.ratio_ka_ks || 0;
        const kaKsLabel = kaKs < 1 ? 'Seleccion purificadora' : kaKs > 1 ? 'Seleccion positiva' : 'Neutral';

        container.innerHTML = `
            <!-- Mini-nav -->
            <div class="sticky top-0 z-20 bg-card/90 backdrop-blur-sm rounded-xl p-3 border border-slate-200 mb-6 flex flex-wrap gap-2">
                <a href="#sec-resumen" class="px-3 py-1 bg-emerald-500/10 text-emerald-500 text-xs rounded-full hover:bg-emerald-500/20 transition">Resumen</a>
                <a href="#sec-ortologos" class="px-3 py-1 bg-violet-500/10 text-violet-500 text-xs rounded-full hover:bg-violet-500/20 transition">Ortologos</a>
                <a href="#sec-mutaciones" class="px-3 py-1 bg-rose-500/10 text-rose-500 text-xs rounded-full hover:bg-rose-500/20 transition">Mutaciones</a>
                <a href="#sec-sintenia" class="px-3 py-1 bg-blue-500/10 text-blue-500 text-xs rounded-full hover:bg-blue-500/20 transition">Sintenia</a>
                <a href="#sec-islas" class="px-3 py-1 bg-amber-500/10 text-amber-500 text-xs rounded-full hover:bg-amber-500/20 transition">Islas Genomicas</a>
                <a href="#sec-funcional" class="px-3 py-1 bg-cyan-500/10 text-cyan-500 text-xs rounded-full hover:bg-cyan-500/20 transition">Funcional</a>
                <a href="#sec-resistencia" class="px-3 py-1 bg-red-500/10 text-red-500 text-xs rounded-full hover:bg-red-500/20 transition">Resistencia</a>
                <a href="#sec-arn" class="px-3 py-1 bg-indigo-500/10 text-indigo-500 text-xs rounded-full hover:bg-indigo-500/20 transition">tRNA/rRNA</a>
                <a href="#sec-gc" class="px-3 py-1 bg-teal-500/10 text-teal-500 text-xs rounded-full hover:bg-teal-500/20 transition">GC Regional</a>
                <a href="#sec-codones" class="px-3 py-1 bg-pink-500/10 text-pink-500 text-xs rounded-full hover:bg-pink-500/20 transition">Codones</a>
            </div>

            <!-- =============== SECCION: RESUMEN =============== -->
            <div id="sec-resumen" class="mb-8">
                <h2 class="text-xl font-bold text-primary mb-4 flex items-center gap-2"><span class="w-1.5 h-6 bg-emerald-500 rounded-full"></span> Resumen General</h2>

                <!-- ANI Badge -->
                <div class="bg-gradient-to-r from-${aniColor}-500/10 to-${aniColor}-500/5 rounded-xl p-5 border border-${aniColor}-500/30 mb-4">
                    <div class="flex items-center justify-between">
                        <div>
                            <p class="text-xs text-secondary">Identidad Nucleotidica Promedio (ANI)</p>
                            <p class="text-3xl font-black text-${aniColor}-500">${ani}%</p>
                            <p class="text-sm text-secondary mt-1">${aniData.clasificacion || ''}</p>
                        </div>
                        <div class="text-right text-sm">
                            <p class="text-secondary">Ka/Ks: <span class="font-bold text-primary">${kaKs}</span></p>
                            <p class="text-xs text-${kaKs < 1 ? 'emerald' : kaKs > 1 ? 'red' : 'amber'}-500">${kaKsLabel}</p>
                            <p class="text-secondary mt-1">Ti/Tv: <span class="font-bold text-primary">${mutStats.ratio_ti_tv || 0}</span></p>
                        </div>
                    </div>
                </div>

                <!-- Stats lado a lado -->
                <div class="grid grid-cols-1 md:grid-cols-2 gap-4 mb-4">
                    <div class="bg-card rounded-xl p-5 border border-emerald-500/30">
                        <h3 class="text-lg font-bold text-emerald-500 mb-3">${nombre1}</h3>
                        <div class="space-y-1.5 text-sm">
                            <p><span class="text-secondary">Longitud:</span> <strong class="text-primary">${this.fmt(ec.longitud_genoma_pb)} pb</strong></p>
                            <p><span class="text-secondary">Genes CDS:</span> <strong class="text-primary">${this.fmt(ec.total_genes_cds)}</strong></p>
                            <p><span class="text-secondary">GC:</span> <strong class="text-primary">${ec.contenido_gc_porcentaje}%</strong></p>
                            <p><span class="text-secondary">Densidad:</span> <strong class="text-primary">${ec.densidad_genica_porcentaje}%</strong></p>
                        </div>
                    </div>
                    <div class="bg-card rounded-xl p-5 border border-cyan-500/30">
                        <h3 class="text-lg font-bold text-cyan-500 mb-3">${nombre2}</h3>
                        <div class="space-y-1.5 text-sm">
                            <p><span class="text-secondary">Longitud:</span> <strong class="text-primary">${this.fmt(sal.longitud_genoma_pb)} pb</strong></p>
                            <p><span class="text-secondary">Genes CDS:</span> <strong class="text-primary">${this.fmt(sal.total_genes_cds)}</strong></p>
                            <p><span class="text-secondary">GC:</span> <strong class="text-primary">${sal.contenido_gc_porcentaje}%</strong></p>
                            <p><span class="text-secondary">Densidad:</span> <strong class="text-primary">${sal.densidad_genica_porcentaje}%</strong></p>
                        </div>
                    </div>
                </div>

                <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-4">
                    <div class="bg-card rounded-xl p-5 border border-slate-200">
                        <h3 class="text-sm font-semibold text-primary mb-4">Metricas Generales</h3>
                        <canvas id="chart-compare-metrics" height="250"></canvas>
                    </div>
                    <div class="bg-card rounded-xl p-5 border border-slate-200">
                        <h3 class="text-sm font-semibold text-primary mb-4">Distribucion Tamanos de Genes</h3>
                        <canvas id="chart-compare-sizes" height="250"></canvas>
                    </div>
                </div>
            </div>

            <!-- =============== SECCION: ORTOLOGOS =============== -->
            <div id="sec-ortologos" class="mb-8">
                <h2 class="text-xl font-bold text-primary mb-4 flex items-center gap-2"><span class="w-1.5 h-6 bg-violet-500 rounded-full"></span> Genes Ortologos (Compartidos vs Unicos)</h2>

                <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-4">
                    ${this.statsCard('Compartidos', this.fmt(ortoStats.total_ortologos || 0), 'Genes ortologos', 'emerald')}
                    ${this.statsCard('Unicos ' + n1, this.fmt(ortoStats.unicos_genoma_1 || 0), (ortoStats.porcentaje_compartido_1 || 0) + '% compartido')}
                    ${this.statsCard('Unicos ' + n2, this.fmt(ortoStats.unicos_genoma_2 || 0), (ortoStats.porcentaje_compartido_2 || 0) + '% compartido', 'cyan')}
                    ${this.statsCard('Hipoteticas', (ortoStats.hipoteticos_genoma_1 || 0) + ' / ' + (ortoStats.hipoteticos_genoma_2 || 0), 'G1 / G2', 'amber')}
                </div>

                <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-4">
                    <div class="bg-card rounded-xl p-5 border border-slate-200">
                        <h3 class="text-sm font-semibold text-primary mb-4">Distribucion de Genes</h3>
                        <canvas id="chart-ortologos-dist" height="250"></canvas>
                    </div>
                    <div class="bg-card rounded-xl p-5 border border-slate-200">
                        <h3 class="text-sm font-semibold text-primary mb-4">Conservacion de Hebra en Ortologos</h3>
                        <canvas id="chart-ortologos-hebra" height="250"></canvas>
                    </div>
                </div>

                <!-- Genes unicos G1 -->
                ${(ortologos.unicos_genoma_1 || []).length > 0 ? `
                <div class="bg-card rounded-xl p-5 border border-emerald-500/20 mb-4">
                    <h4 class="text-sm font-semibold text-emerald-500 mb-3">Genes Unicos en ${n1} (top ${ortologos.total_unicos_1 || 0})</h4>
                    <div class="overflow-x-auto max-h-[250px] overflow-y-auto">
                        <table class="w-full text-xs"><thead class="sticky top-0 bg-card"><tr>
                            <th class="text-left px-2 py-1 text-secondary border-b border-slate-200">Locus</th>
                            <th class="text-left px-2 py-1 text-secondary border-b border-slate-200">Gen</th>
                            <th class="text-left px-2 py-1 text-secondary border-b border-slate-200">Producto</th>
                            <th class="text-right px-2 py-1 text-secondary border-b border-slate-200">Longitud</th>
                        </tr></thead><tbody>
                        ${(ortologos.unicos_genoma_1 || []).slice(0, 30).map(g => `<tr class="hover:bg-slate-50 dark:hover:bg-slate-800"><td class="px-2 py-1 font-mono">${g.locus}</td><td class="px-2 py-1 font-medium text-primary">${g.nombre || '-'}</td><td class="px-2 py-1 text-secondary truncate max-w-[200px]">${g.producto}</td><td class="px-2 py-1 text-right">${this.fmt(g.longitud)} pb</td></tr>`).join('')}
                        </tbody></table>
                    </div>
                </div>` : ''}

                <!-- Genes unicos G2 -->
                ${(ortologos.unicos_genoma_2 || []).length > 0 ? `
                <div class="bg-card rounded-xl p-5 border border-cyan-500/20 mb-4">
                    <h4 class="text-sm font-semibold text-cyan-500 mb-3">Genes Unicos en ${n2} (top ${ortologos.total_unicos_2 || 0})</h4>
                    <div class="overflow-x-auto max-h-[250px] overflow-y-auto">
                        <table class="w-full text-xs"><thead class="sticky top-0 bg-card"><tr>
                            <th class="text-left px-2 py-1 text-secondary border-b border-slate-200">Locus</th>
                            <th class="text-left px-2 py-1 text-secondary border-b border-slate-200">Gen</th>
                            <th class="text-left px-2 py-1 text-secondary border-b border-slate-200">Producto</th>
                            <th class="text-right px-2 py-1 text-secondary border-b border-slate-200">Longitud</th>
                        </tr></thead><tbody>
                        ${(ortologos.unicos_genoma_2 || []).slice(0, 30).map(g => `<tr class="hover:bg-slate-50 dark:hover:bg-slate-800"><td class="px-2 py-1 font-mono">${g.locus}</td><td class="px-2 py-1 font-medium text-primary">${g.nombre || '-'}</td><td class="px-2 py-1 text-secondary truncate max-w-[200px]">${g.producto}</td><td class="px-2 py-1 text-right">${this.fmt(g.longitud)} pb</td></tr>`).join('')}
                        </tbody></table>
                    </div>
                </div>` : ''}

                <!-- Ortologos mas divergentes -->
                ${(ortologos.ortologos_divergentes || []).length > 0 ? `
                <div class="bg-card rounded-xl p-5 border border-amber-500/20 mb-4">
                    <h4 class="text-sm font-semibold text-amber-500 mb-3">Ortologos Mas Divergentes (mayor diferencia de tamano)</h4>
                    <div class="overflow-x-auto max-h-[250px] overflow-y-auto">
                        <table class="w-full text-xs"><thead class="sticky top-0 bg-card"><tr>
                            <th class="text-left px-2 py-1 text-secondary border-b border-slate-200">Producto</th>
                            <th class="text-right px-2 py-1 text-secondary border-b border-slate-200">Long. G1</th>
                            <th class="text-right px-2 py-1 text-secondary border-b border-slate-200">Long. G2</th>
                            <th class="text-right px-2 py-1 text-secondary border-b border-slate-200">Diferencia</th>
                        </tr></thead><tbody>
                        ${(ortologos.ortologos_divergentes || []).slice(0, 20).map(o => `<tr class="hover:bg-slate-50 dark:hover:bg-slate-800"><td class="px-2 py-1 text-primary truncate max-w-[250px]">${o.producto}</td><td class="px-2 py-1 text-right">${this.fmt(o.longitud_1)} pb</td><td class="px-2 py-1 text-right">${this.fmt(o.longitud_2)} pb</td><td class="px-2 py-1 text-right font-bold text-amber-500">${this.fmt(o.dif_longitud)} pb</td></tr>`).join('')}
                        </tbody></table>
                    </div>
                </div>` : ''}
            </div>

            <!-- =============== SECCION: MUTACIONES =============== -->
            <div id="sec-mutaciones" class="mb-8">
                <h2 class="text-xl font-bold text-primary mb-4 flex items-center gap-2"><span class="w-1.5 h-6 bg-rose-500 rounded-full"></span> Mutaciones en Ortologos</h2>

                <div class="grid grid-cols-2 md:grid-cols-5 gap-3 mb-4">
                    ${this.statsCard('ANI', mutStats.ani_ortologos + '%', 'Identidad nucleotidica', 'emerald')}
                    ${this.statsCard('SNPs', this.fmt(mutStats.total_snps_codones || 0), 'en codones comparados')}
                    ${this.statsCard('Sinonimas', this.fmt(mutStats.mutaciones_sinonimas || 0), 'misma proteina', 'cyan')}
                    ${this.statsCard('No-sinonimas', this.fmt(mutStats.mutaciones_no_sinonimas || 0), 'cambia proteina', 'rose')}
                    ${this.statsCard('Ka/Ks', kaKs, kaKsLabel, 'amber')}
                </div>

                <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-4">
                    <div class="bg-card rounded-xl p-5 border border-slate-200">
                        <h3 class="text-sm font-semibold text-primary mb-4">Distribucion de Identidad entre Ortologos</h3>
                        <canvas id="chart-mut-identidad" height="250"></canvas>
                    </div>
                    <div class="bg-card rounded-xl p-5 border border-slate-200">
                        <h3 class="text-sm font-semibold text-primary mb-4">Tipo de Mutaciones</h3>
                        <canvas id="chart-mut-tipos" height="250"></canvas>
                    </div>
                </div>

                <!-- Genes mas divergentes -->
                ${(mutaciones.genes_mas_divergentes || []).length > 0 ? `
                <div class="bg-card rounded-xl p-5 border border-rose-500/20 mb-4">
                    <h4 class="text-sm font-semibold text-rose-500 mb-3">Genes con Mayor Divergencia</h4>
                    <div class="overflow-x-auto max-h-[300px] overflow-y-auto">
                        <table class="w-full text-xs"><thead class="sticky top-0 bg-card"><tr>
                            <th class="text-left px-2 py-1 text-secondary border-b border-slate-200">Producto</th>
                            <th class="text-right px-2 py-1 text-secondary border-b border-slate-200">Identidad</th>
                            <th class="text-right px-2 py-1 text-secondary border-b border-slate-200">SNPs</th>
                            <th class="text-right px-2 py-1 text-secondary border-b border-slate-200">Sinonimas</th>
                            <th class="text-right px-2 py-1 text-secondary border-b border-slate-200">No-sinon.</th>
                            <th class="text-right px-2 py-1 text-secondary border-b border-slate-200">Dif. Long.</th>
                        </tr></thead><tbody>
                        ${(mutaciones.genes_mas_divergentes || []).slice(0, 25).map(g => {
                            const idColor = g.identidad_nt >= 95 ? 'emerald' : g.identidad_nt >= 80 ? 'amber' : 'red';
                            return '<tr class="hover:bg-slate-50 dark:hover:bg-slate-800"><td class="px-2 py-1 text-primary truncate max-w-[200px]">' + g.producto + '</td><td class="px-2 py-1 text-right font-bold text-' + idColor + '-500">' + g.identidad_nt + '%</td><td class="px-2 py-1 text-right">' + g.snps_codones + '</td><td class="px-2 py-1 text-right text-emerald-500">' + g.sinonimas + '</td><td class="px-2 py-1 text-right text-rose-500">' + g.no_sinonimas + '</td><td class="px-2 py-1 text-right">' + this.fmt(g.dif_longitud) + ' pb</td></tr>';
                        }).join('')}
                        </tbody></table>
                    </div>
                </div>` : ''}
            </div>

            <!-- =============== SECCION: SINTENIA =============== -->
            <div id="sec-sintenia" class="mb-8">
                <h2 class="text-xl font-bold text-primary mb-4 flex items-center gap-2"><span class="w-1.5 h-6 bg-blue-500 rounded-full"></span> Sintenia (Orden Genico)</h2>

                <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-4">
                    ${this.statsCard('Bloques', sinStats.bloques_sintenicos || 0, 'conservados', 'blue')}
                    ${this.statsCard('Genes en bloques', this.fmt(sinStats.genes_en_bloques || 0), (sinStats.porcentaje_sintenico || 0) + '% sintenico', 'emerald')}
                    ${this.statsCard('Inversiones', sinStats.inversiones_detectadas || 0, 'cambios de orden', 'amber')}
                    ${this.statsCard('Cambio hebra', sinStats.cambios_hebra || 0, 'ortologos en otra hebra', 'rose')}
                </div>

                <div class="bg-card rounded-xl p-5 border border-slate-200 mb-4">
                    <h3 class="text-sm font-semibold text-primary mb-2">Dot Plot de Sintenia</h3>
                    <p class="text-xs text-secondary mb-3">Cada punto = un gen ortologo. Verde = misma hebra, Rojo = diferente hebra. Puntos en diagonal = orden conservado.</p>
                    <canvas id="canvas-dotplot" height="400" style="max-height:500px;"></canvas>
                </div>
            </div>

            <!-- =============== SECCION: ISLAS GENOMICAS =============== -->
            <div id="sec-islas" class="mb-8">
                <h2 class="text-xl font-bold text-primary mb-4 flex items-center gap-2"><span class="w-1.5 h-6 bg-amber-500 rounded-full"></span> Islas Genomicas (Transferencia Horizontal)</h2>

                <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-4">
                    ${this.statsCard('Islas ' + n1, islStats.islas_genoma_1 || 0, this.fmt(islStats.genes_en_islas_1 || 0) + ' genes', 'emerald')}
                    ${this.statsCard('Islas ' + n2, islStats.islas_genoma_2 || 0, this.fmt(islStats.genes_en_islas_2 || 0) + ' genes', 'cyan')}
                    ${this.statsCard('HGT ' + n1, islStats.posibles_hgt_1 || 0, 'GC anormal >3%', 'amber')}
                    ${this.statsCard('HGT ' + n2, islStats.posibles_hgt_2 || 0, 'GC anormal >3%', 'rose')}
                </div>

                <div class="grid grid-cols-1 md:grid-cols-2 gap-4 mb-4">
                    ${(islas.islas_genoma_1 || []).length > 0 ? `
                    <div class="bg-card rounded-xl p-5 border border-emerald-500/20">
                        <h4 class="text-sm font-semibold text-emerald-500 mb-3">Islas en ${n1}</h4>
                        <div class="space-y-2 max-h-[300px] overflow-y-auto">
                        ${(islas.islas_genoma_1 || []).slice(0, 10).map((isla, i) => `
                            <div class="p-3 rounded-lg ${isla.posible_hgt ? 'bg-amber-500/5 border border-amber-500/20' : 'bg-slate-50 dark:bg-slate-800'}">
                                <div class="flex justify-between text-xs mb-1">
                                    <span class="font-bold text-primary">Isla ${i+1}: ${isla.num_genes} genes</span>
                                    <span class="text-secondary">${this.fmt(isla.longitud_pb)} pb</span>
                                </div>
                                <p class="text-xs text-secondary">Pos: ${this.fmt(isla.inicio)}-${this.fmt(isla.fin)} | GC: ${isla.gc_promedio}% (genoma: ${isla.gc_genoma}%, desv: ${isla.desviacion_gc > 0 ? '+' : ''}${isla.desviacion_gc}%)</p>
                                ${isla.posible_hgt ? '<p class="text-xs text-amber-500 font-medium mt-1">Posible transferencia horizontal</p>' : ''}
                                <p class="text-[10px] text-secondary mt-1">${(isla.genes || []).map(g => g.nombre || g.locus).join(', ')}</p>
                            </div>
                        `).join('')}
                        </div>
                    </div>` : '<div></div>'}

                    ${(islas.islas_genoma_2 || []).length > 0 ? `
                    <div class="bg-card rounded-xl p-5 border border-cyan-500/20">
                        <h4 class="text-sm font-semibold text-cyan-500 mb-3">Islas en ${n2}</h4>
                        <div class="space-y-2 max-h-[300px] overflow-y-auto">
                        ${(islas.islas_genoma_2 || []).slice(0, 10).map((isla, i) => `
                            <div class="p-3 rounded-lg ${isla.posible_hgt ? 'bg-amber-500/5 border border-amber-500/20' : 'bg-slate-50 dark:bg-slate-800'}">
                                <div class="flex justify-between text-xs mb-1">
                                    <span class="font-bold text-primary">Isla ${i+1}: ${isla.num_genes} genes</span>
                                    <span class="text-secondary">${this.fmt(isla.longitud_pb)} pb</span>
                                </div>
                                <p class="text-xs text-secondary">Pos: ${this.fmt(isla.inicio)}-${this.fmt(isla.fin)} | GC: ${isla.gc_promedio}% (genoma: ${isla.gc_genoma}%, desv: ${isla.desviacion_gc > 0 ? '+' : ''}${isla.desviacion_gc}%)</p>
                                ${isla.posible_hgt ? '<p class="text-xs text-amber-500 font-medium mt-1">Posible transferencia horizontal</p>' : ''}
                                <p class="text-[10px] text-secondary mt-1">${(isla.genes || []).map(g => g.nombre || g.locus).join(', ')}</p>
                            </div>
                        `).join('')}
                        </div>
                    </div>` : '<div></div>'}
                </div>
            </div>

            <!-- =============== SECCION: PERFIL FUNCIONAL =============== -->
            <div id="sec-funcional" class="mb-8">
                <h2 class="text-xl font-bold text-primary mb-4 flex items-center gap-2"><span class="w-1.5 h-6 bg-cyan-500 rounded-full"></span> Perfil Funcional Comparativo</h2>
                <div class="bg-card rounded-xl p-5 border border-slate-200 mb-4">
                    <canvas id="chart-funcional" height="350"></canvas>
                </div>
            </div>

            <!-- =============== SECCION: RESISTENCIA =============== -->
            <div id="sec-resistencia" class="mb-8">
                <h2 class="text-xl font-bold text-primary mb-4 flex items-center gap-2"><span class="w-1.5 h-6 bg-red-500 rounded-full"></span> Resistencia Antibiotica</h2>

                <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-4">
                    ${this.statsCard('Genes R. ' + n1, resStats.total_genes_resistencia_1 || 0, (resStats.categorias_con_genes_1 || 0) + ' categorias', 'emerald')}
                    ${this.statsCard('Genes R. ' + n2, resStats.total_genes_resistencia_2 || 0, (resStats.categorias_con_genes_2 || 0) + ' categorias', 'cyan')}
                    ${this.statsCard('Virulencia ' + n1, (virulencia.ecoli || {}).total || 0, 'genes', 'amber')}
                    ${this.statsCard('Virulencia ' + n2, (virulencia.salmonella || {}).total || 0, 'genes', 'rose')}
                </div>

                <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-4">
                    <div class="bg-card rounded-xl p-5 border border-slate-200">
                        <h3 class="text-sm font-semibold text-primary mb-4">Resistencia por Categoria</h3>
                        <canvas id="chart-resistencia" height="300"></canvas>
                    </div>
                    <div class="bg-card rounded-xl p-5 border border-slate-200">
                        <h3 class="text-sm font-semibold text-primary mb-4">Virulencia por Categoria</h3>
                        <canvas id="chart-compare-vir-cats" height="300"></canvas>
                    </div>
                </div>
            </div>

            <!-- =============== SECCION: tRNA/rRNA =============== -->
            <div id="sec-arn" class="mb-8">
                <h2 class="text-xl font-bold text-primary mb-4 flex items-center gap-2"><span class="w-1.5 h-6 bg-indigo-500 rounded-full"></span> tRNA y rRNA</h2>

                <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-4">
                    ${this.statsCard('tRNA ' + n1, rnaStats.total_trna_1 || 0, (rnaStats.tipos_trna_1 || 0) + ' tipos', 'emerald')}
                    ${this.statsCard('tRNA ' + n2, rnaStats.total_trna_2 || 0, (rnaStats.tipos_trna_2 || 0) + ' tipos', 'cyan')}
                    ${this.statsCard('rRNA ' + n1, rnaStats.total_rrna_1 || 0, 'operones ribosomales', 'indigo')}
                    ${this.statsCard('rRNA ' + n2, rnaStats.total_rrna_2 || 0, 'operones ribosomales', 'violet')}
                </div>

                <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-4">
                    <div class="bg-card rounded-xl p-5 border border-slate-200">
                        <h3 class="text-sm font-semibold text-primary mb-4">tRNA por Aminoacido</h3>
                        <canvas id="chart-trna" height="300"></canvas>
                    </div>
                    <div class="bg-card rounded-xl p-5 border border-slate-200">
                        <h3 class="text-sm font-semibold text-primary mb-4">rRNA por Subtipo</h3>
                        <canvas id="chart-rrna" height="250"></canvas>
                    </div>
                </div>
            </div>

            <!-- =============== SECCION: GC REGIONAL =============== -->
            <div id="sec-gc" class="mb-8">
                <h2 class="text-xl font-bold text-primary mb-4 flex items-center gap-2"><span class="w-1.5 h-6 bg-teal-500 rounded-full"></span> GC Regional (Ventana Deslizante)</h2>

                <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-4">
                    ${this.statsCard('GC ' + n1, (gcStats.media_gc_1 || 0) + '%', 'std: ' + (gcStats.std_gc_1 || 0) + '%', 'emerald')}
                    ${this.statsCard('GC ' + n2, (gcStats.media_gc_2 || 0) + '%', 'std: ' + (gcStats.std_gc_2 || 0) + '%', 'cyan')}
                    ${this.statsCard('Anomalas ' + n1, gcStats.regiones_anomalas_1 || 0, '>2 std desviacion', 'amber')}
                    ${this.statsCard('Anomalas ' + n2, gcStats.regiones_anomalas_2 || 0, '>2 std desviacion', 'rose')}
                </div>

                <div class="bg-card rounded-xl p-5 border border-slate-200 mb-4">
                    <h3 class="text-sm font-semibold text-primary mb-2">GC% a lo Largo del Genoma</h3>
                    <p class="text-xs text-secondary mb-3">Ventana de ${this.fmt(gcStats.ventana_pb || 50000)} pb. Regiones con GC anormal pueden indicar transferencia horizontal.</p>
                    <canvas id="chart-gc-regional" height="250"></canvas>
                </div>

                <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-4">
                    <div class="bg-card rounded-xl p-5 border border-slate-200">
                        <h3 class="text-sm font-semibold text-primary mb-4">Distribucion GC por Gen</h3>
                        <canvas id="chart-compare-gc" height="250"></canvas>
                    </div>
                    <div class="bg-card rounded-xl p-5 border border-slate-200">
                        <h3 class="text-sm font-semibold text-primary mb-4">Genes de Virulencia</h3>
                        <canvas id="chart-compare-virulence" height="250"></canvas>
                    </div>
                </div>
            </div>

            <!-- =============== SECCION: CODONES =============== -->
            <div id="sec-codones" class="mb-8">
                <h2 class="text-xl font-bold text-primary mb-4 flex items-center gap-2"><span class="w-1.5 h-6 bg-pink-500 rounded-full"></span> Uso de Codones</h2>
                ${usoCodones.ecoli ? `
                <div class="bg-card rounded-xl p-5 border border-slate-200 mb-4">
                    <canvas id="chart-compare-codons" height="300"></canvas>
                </div>` : ''}
            </div>

            <!-- =============== INTERPRETACION IA =============== -->
            ${interpretacionIA ? `
            <div class="bg-gradient-to-br from-violet-500/5 to-fuchsia-500/5 rounded-xl p-6 border border-violet-500/20 mb-6">
                <div class="flex items-center gap-3 mb-4">
                    <div class="w-10 h-10 rounded-xl bg-gradient-to-br from-violet-500 to-fuchsia-500 flex items-center justify-center text-xl shadow-lg">*</div>
                    <div>
                        <h3 class="font-bold text-primary">Interpretacion IA Avanzada</h3>
                        <p class="text-xs text-secondary">Analisis biologico de todos los datos comparativos (Gemini AI)</p>
                    </div>
                </div>
                <div class="text-sm text-primary leading-relaxed space-y-3">
                    ${interpretacionIA.split('\n').filter(p => p.trim()).map(p =>
                        `<p>${p.replace(/\*\*(.*?)\*\*/g, '<strong>$1</strong>').replace(/\*(.*?)\*/g, '<em>$1</em>')}</p>`
                    ).join('')}
                </div>
            </div>` : ''}
        `;

        // ===== CREAR TODOS LOS CHARTS =====
        setTimeout(() => {
            const self = this;
            // Metricas
            this.createChart('chart-compare-metrics', {
                type: 'bar', data: { labels: ['Longitud (Mb)', 'Genes', 'GC %', 'Densidad %'],
                datasets: [
                    { label: n1, data: [(ec.longitud_genoma_pb||0)/1e6, ec.total_genes_cds||0, ec.contenido_gc_porcentaje||0, ec.densidad_genica_porcentaje||0], backgroundColor: '#10b981', borderRadius: 4 },
                    { label: n2, data: [(sal.longitud_genoma_pb||0)/1e6, sal.total_genes_cds||0, sal.contenido_gc_porcentaje||0, sal.densidad_genica_porcentaje||0], backgroundColor: '#06b6d4', borderRadius: 4 }
                ]}, options: { responsive: true, plugins: { legend: { position: 'top' } }, scales: { y: { beginAtZero: true } } }
            });

            // Tamanos
            const r1 = distTamanos.ecoli?.rangos || {}; const r2 = distTamanos.salmonella?.rangos || {};
            if (Object.keys(r1).length > 0) {
                const lb = Object.keys(r1);
                this.createChart('chart-compare-sizes', { type: 'bar', data: { labels: lb.map(l => l + ' pb'),
                    datasets: [{ label: n1, data: lb.map(l => r1[l]||0), backgroundColor: '#10b981', borderRadius: 4 },
                               { label: n2, data: lb.map(l => r2[l]||0), backgroundColor: '#06b6d4', borderRadius: 4 }]},
                    options: { responsive: true, plugins: { legend: { position: 'top' } }, scales: { y: { beginAtZero: true } } }
                });
            }

            // Ortologos distribucion
            this.createChart('chart-ortologos-dist', { type: 'doughnut', data: {
                labels: ['Compartidos', 'Unicos ' + n1, 'Unicos ' + n2, 'Hipoteticas G1', 'Hipoteticas G2'],
                datasets: [{ data: [ortoStats.total_ortologos||0, ortoStats.unicos_genoma_1||0, ortoStats.unicos_genoma_2||0, ortoStats.hipoteticos_genoma_1||0, ortoStats.hipoteticos_genoma_2||0],
                    backgroundColor: ['#10b981', '#f59e0b', '#06b6d4', '#94a3b8', '#64748b'] }]},
                options: { responsive: true, plugins: { legend: { position: 'right' } } }
            });

            // Ortologos hebra
            this.createChart('chart-ortologos-hebra', { type: 'doughnut', data: {
                labels: ['Misma hebra', 'Diferente hebra'],
                datasets: [{ data: [ortoStats.conservacion_hebra||0, ortoStats.cambio_hebra||0],
                    backgroundColor: ['#10b981', '#ef4444'] }]},
                options: { responsive: true, plugins: { legend: { position: 'right' } } }
            });

            // Mutaciones - identidad
            const distId = mutaciones.distribucion_identidad || {};
            if (Object.keys(distId).length > 0) {
                const idLabels = Object.keys(distId);
                this.createChart('chart-mut-identidad', { type: 'bar', data: { labels: idLabels,
                    datasets: [{ label: 'Genes', data: idLabels.map(l => distId[l]||0),
                        backgroundColor: idLabels.map(l => l === '100%' ? '#10b981' : l.includes('99') ? '#34d399' : l.includes('95') ? '#fbbf24' : l.includes('90') ? '#f59e0b' : '#ef4444'),
                        borderRadius: 4 }]},
                    options: { responsive: true, plugins: { legend: { display: false } }, scales: { y: { beginAtZero: true } } }
                });
            }

            // Mutaciones tipos
            this.createChart('chart-mut-tipos', { type: 'doughnut', data: {
                labels: ['Sinonimas', 'No-sinonimas', 'Transiciones', 'Transversiones'],
                datasets: [{ data: [mutStats.mutaciones_sinonimas||0, mutStats.mutaciones_no_sinonimas||0, mutStats.transiciones||0, mutStats.transversiones||0],
                    backgroundColor: ['#10b981', '#ef4444', '#3b82f6', '#f59e0b'] }]},
                options: { responsive: true, plugins: { legend: { position: 'right' } } }
            });

            // Dot plot sintenia (canvas)
            const dotCanvas = document.getElementById('canvas-dotplot');
            if (dotCanvas && (sintenia.puntos_dotplot || []).length > 0) {
                const ctx = dotCanvas.getContext('2d');
                const w = dotCanvas.width = dotCanvas.clientWidth || 600;
                const h = dotCanvas.height = Math.min(dotCanvas.clientHeight || 400, 500);
                const m = 50;
                const len1 = ec.longitud_genoma_pb || 5000000;
                const len2 = sal.longitud_genoma_pb || 5000000;
                ctx.fillStyle = '#0f172a'; ctx.fillRect(0, 0, w, h);
                // Axes
                ctx.strokeStyle = '#475569'; ctx.lineWidth = 1;
                ctx.beginPath(); ctx.moveTo(m, m); ctx.lineTo(m, h - m); ctx.lineTo(w - m, h - m); ctx.stroke();
                // Labels
                ctx.fillStyle = '#94a3b8'; ctx.font = '10px sans-serif'; ctx.textAlign = 'center';
                ctx.fillText(n1, w / 2, h - 5);
                ctx.save(); ctx.translate(12, h / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(n2, 0, 0); ctx.restore();
                // Scale ticks
                for (let i = 0; i <= 4; i++) {
                    const xp = m + (i / 4) * (w - 2 * m);
                    const yp = h - m - (i / 4) * (h - 2 * m);
                    ctx.fillStyle = '#64748b'; ctx.fillText(((i / 4) * len1 / 1e6).toFixed(1) + 'Mb', xp, h - m + 14);
                    ctx.textAlign = 'right'; ctx.fillText(((i / 4) * len2 / 1e6).toFixed(1), m - 5, yp + 3); ctx.textAlign = 'center';
                }
                // Points
                const pts = sintenia.puntos_dotplot;
                for (const p of pts) {
                    const px = m + (p.x / len1) * (w - 2 * m);
                    const py = h - m - (p.y / len2) * (h - 2 * m);
                    ctx.fillStyle = p.misma_hebra ? '#10b981' : '#ef4444';
                    ctx.fillRect(px - 1, py - 1, 2, 2);
                }
            }

            // Funcional
            const funcCats = funcional.categorias || [];
            if (funcCats.length > 0) {
                this.createChart('chart-funcional', { type: 'bar', data: {
                    labels: funcCats.map(c => c.categoria),
                    datasets: [
                        { label: n1, data: funcCats.map(c => c.genoma_1), backgroundColor: '#10b981', borderRadius: 3 },
                        { label: n2, data: funcCats.map(c => c.genoma_2), backgroundColor: '#06b6d4', borderRadius: 3 }
                    ]}, options: { indexAxis: 'y', responsive: true, plugins: { legend: { position: 'top' } }, scales: { x: { beginAtZero: true } } }
                });
            }

            // Resistencia
            const resCmp = (resistencia.comparacion || []).filter(c => c.genoma_1 > 0 || c.genoma_2 > 0);
            if (resCmp.length > 0) {
                this.createChart('chart-resistencia', { type: 'bar', data: {
                    labels: resCmp.map(c => c.categoria),
                    datasets: [
                        { label: n1, data: resCmp.map(c => c.genoma_1), backgroundColor: '#10b981', borderRadius: 3 },
                        { label: n2, data: resCmp.map(c => c.genoma_2), backgroundColor: '#06b6d4', borderRadius: 3 }
                    ]}, options: { indexAxis: 'y', responsive: true, plugins: { legend: { position: 'top' } }, scales: { x: { beginAtZero: true } } }
                });
            }

            // Virulencia categorias
            if (virulencia.ecoli?.categorias && virulencia.salmonella?.categorias) {
                const c1 = virulencia.ecoli.categorias; const c2 = virulencia.salmonella.categorias;
                const allC = [...new Set([...Object.keys(c1), ...Object.keys(c2)])];
                this.createChart('chart-compare-vir-cats', { type: 'bar', data: { labels: allC,
                    datasets: [{ label: n1, data: allC.map(c => c1[c]||0), backgroundColor: '#10b981', borderRadius: 3 },
                               { label: n2, data: allC.map(c => c2[c]||0), backgroundColor: '#06b6d4', borderRadius: 3 }]},
                    options: { indexAxis: 'y', responsive: true, plugins: { legend: { position: 'top' } }, scales: { x: { beginAtZero: true } } }
                });
            }

            // Virulencia totales
            if (virulencia.ecoli && virulencia.salmonella) {
                this.createChart('chart-compare-virulence', { type: 'bar', data: { labels: [n1, n2],
                    datasets: [{ label: 'Genes', data: [virulencia.ecoli.total||0, virulencia.salmonella.total||0], backgroundColor: ['#10b981','#06b6d4'], borderRadius: 6 }]},
                    options: { responsive: true, plugins: { legend: { display: false } }, scales: { y: { beginAtZero: true } } }
                });
            }

            // tRNA
            const trnaData = trnaRrna.trna_por_aminoacido || [];
            if (trnaData.length > 0) {
                this.createChart('chart-trna', { type: 'bar', data: { labels: trnaData.map(t => t.aminoacido),
                    datasets: [{ label: n1, data: trnaData.map(t => t.genoma_1), backgroundColor: '#10b981', borderRadius: 2 },
                               { label: n2, data: trnaData.map(t => t.genoma_2), backgroundColor: '#06b6d4', borderRadius: 2 }]},
                    options: { responsive: true, plugins: { legend: { position: 'top' } }, scales: { y: { beginAtZero: true } } }
                });
            }

            // rRNA
            const rrnaData = trnaRrna.rrna_por_subtipo || [];
            if (rrnaData.length > 0) {
                this.createChart('chart-rrna', { type: 'bar', data: { labels: rrnaData.map(r => r.subtipo),
                    datasets: [{ label: n1, data: rrnaData.map(r => r.genoma_1), backgroundColor: '#6366f1', borderRadius: 4 },
                               { label: n2, data: rrnaData.map(r => r.genoma_2), backgroundColor: '#8b5cf6', borderRadius: 4 }]},
                    options: { responsive: true, plugins: { legend: { position: 'top' } }, scales: { y: { beginAtZero: true } } }
                });
            }

            // GC Regional line chart
            const gc1 = gcRegional.gc_genoma_1 || []; const gc2 = gcRegional.gc_genoma_2 || [];
            if (gc1.length > 0 || gc2.length > 0) {
                this.createChart('chart-gc-regional', { type: 'line', data: {
                    labels: gc1.map(p => (p.posicion / 1e6).toFixed(2) + ' Mb'),
                    datasets: [
                        { label: n1, data: gc1.map(p => p.gc), borderColor: '#10b981', backgroundColor: '#10b98120', fill: true, pointRadius: 0, borderWidth: 1.5, tension: 0.3 },
                        { label: n2, data: gc2.map(p => p.gc), borderColor: '#06b6d4', backgroundColor: '#06b6d420', fill: true, pointRadius: 0, borderWidth: 1.5, tension: 0.3 }
                    ]}, options: { responsive: true, plugins: { legend: { position: 'top' } },
                    scales: { y: { title: { display: true, text: 'GC %' } }, x: { ticks: { maxTicksLimit: 10 } } } }
                });
            }

            // GC por gen
            const gR1 = distGc.ecoli?.rangos || {}; const gR2 = distGc.salmonella?.rangos || {};
            if (Object.keys(gR1).length > 0) {
                const lb = Object.keys(gR1);
                this.createChart('chart-compare-gc', { type: 'bar', data: { labels: lb,
                    datasets: [{ label: n1, data: lb.map(l => gR1[l]||0), backgroundColor: '#10b981', borderRadius: 4 },
                               { label: n2, data: lb.map(l => gR2[l]||0), backgroundColor: '#06b6d4', borderRadius: 4 }]},
                    options: { responsive: true, plugins: { legend: { position: 'top' } }, scales: { y: { beginAtZero: true } } }
                });
            }

            // Codones
            if (usoCodones.ecoli && usoCodones.salmonella) {
                const t1 = usoCodones.ecoli.top_10 || []; const t2 = usoCodones.salmonella.top_10 || [];
                const allC = [...new Set([...t1.map(c=>c.codon),...t2.map(c=>c.codon)])].slice(0,12);
                const f1 = {}; const f2 = {}; t1.forEach(c => f1[c.codon]=c.frecuencia); t2.forEach(c => f2[c.codon]=c.frecuencia);
                this.createChart('chart-compare-codons', { type: 'bar', data: { labels: allC,
                    datasets: [{ label: n1, data: allC.map(c=>f1[c]||0), backgroundColor: '#10b981', borderRadius: 3 },
                               { label: n2, data: allC.map(c=>f2[c]||0), backgroundColor: '#06b6d4', borderRadius: 3 }]},
                    options: { responsive: true, plugins: { legend: { position: 'top' } },
                    scales: { x: { ticks: { font: { family: 'monospace', weight: 'bold' } } }, y: { beginAtZero: true, title: { display: true, text: '%' } } } }
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
                    <button id="btn-3d-spin" onclick="DashboardRenderer._toggle3DSpin()" style="position:absolute; top:10px; right:10px; z-index:10; display:none;" class="px-3 py-1.5 bg-black/60 hover:bg-black/80 text-white text-xs rounded-lg transition backdrop-blur-sm cursor-pointer">\u25B6 Play</button>
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
    _spinning3D: false,

    async _load3DProtein() {
        const select = document.getElementById('protein-3d-select');
        const viewerDiv = document.getElementById('protein-3d-viewer');
        const placeholder = document.getElementById('protein-3d-placeholder');
        const loadingDiv = document.getElementById('protein-3d-loading');
        const infoDiv = document.getElementById('protein-3d-info');
        const infoName = document.getElementById('protein-3d-name');
        const infoDetail = document.getElementById('protein-3d-detail');
        const spinBtn = document.getElementById('btn-3d-spin');

        if (!select || select.value === '') return;

        const idx = parseInt(select.value);
        const proteinData = this._proteinas3DData[idx];
        if (!proteinData) return;

        if (loadingDiv) loadingDiv.classList.remove('hidden');

        try {
            let proxyUrl;
            if (proteinData.source === 'alphafold' && proteinData.url) {
                proxyUrl = `/api/proxy_pdb?pdb_id=${proteinData.pdb_id}&source=alphafold&url=${encodeURIComponent(proteinData.url)}`;
            } else {
                proxyUrl = `/api/proxy_pdb?pdb_id=${proteinData.pdb_id}&source=rcsb`;
            }

            const resp = await fetch(proxyUrl);
            if (!resp.ok) throw new Error('No se pudo descargar la estructura');
            const pdbData = await resp.text();
            if (pdbData.length < 50) throw new Error('Estructura PDB vacia o invalida');
            const format = pdbData.includes('_atom_site') ? 'cif' : 'pdb';

            // Mostrar visor y ocultar placeholder
            if (placeholder) placeholder.style.display = 'none';
            if (viewerDiv) viewerDiv.style.display = 'block';
            if (spinBtn) spinBtn.style.display = 'block';

            // Crear o limpiar visor
            if (this._viewer3D) {
                this._viewer3D.clear();
            } else if (typeof $3Dmol !== 'undefined') {
                viewerDiv.innerHTML = '';
                this._viewer3D = $3Dmol.createViewer(viewerDiv, {
                    backgroundColor: '#1a1a2e',
                    antialias: true
                });

                // Invertir scroll wheel para zoom
                const self = this;
                viewerDiv.addEventListener('wheel', (e) => {
                    e.preventDefault();
                    e.stopPropagation();
                    if (!self._viewer3D) return;
                    const factor = e.deltaY > 0 ? 1.08 : 0.92;
                    self._viewer3D.zoom(factor);
                    self._viewer3D.render();
                }, { passive: false, capture: true });

                // Auto-pausar spin al hacer click
                viewerDiv.addEventListener('mousedown', () => {
                    if (self._spinning3D) {
                        self._spinning3D = false;
                        self._viewer3D.spin(false);
                        const btn = document.getElementById('btn-3d-spin');
                        if (btn) btn.textContent = '\u25B6 Play';
                    }
                });
            } else {
                throw new Error('3Dmol.js no esta cargado. Verifica tu conexion a internet.');
            }

            // Cargar estructura - SIN auto-spin (estatico por defecto)
            this._viewer3D.addModel(pdbData, format);
            this._applyStyle3D(this._current3DStyle);
            this._viewer3D.zoomTo();
            this._viewer3D.render();
            this._spinning3D = false;
            this._viewer3D.spin(false);
            if (spinBtn) spinBtn.textContent = '\u25B6 Play';

            // Mostrar info de la proteina
            const fuente = proteinData.source === 'alphafold' ? 'AlphaFold (prediccion IA)' : 'RCSB PDB (experimental)';
            if (infoDiv) infoDiv.classList.remove('hidden');
            if (infoName) infoName.textContent = (proteinData.nombre_gen || '') + ' - ' + (proteinData.producto || 'Proteina');
            if (infoDetail) infoDetail.textContent = 'Fuente: ' + fuente + ' | ID: ' + (proteinData.pdb_id || 'N/A') + ' | Formato: ' + format.toUpperCase();

        } catch (err) {
            if (placeholder) placeholder.style.display = 'none';
            if (viewerDiv) {
                viewerDiv.style.display = 'flex';
                viewerDiv.innerHTML = '<div class="flex items-center justify-center h-full w-full"><p class="text-red-400 text-sm text-center px-4">' + (err.message || 'Error cargando estructura 3D') + '</p></div>';
            }
            if (spinBtn) spinBtn.style.display = 'none';
        }

        if (loadingDiv) loadingDiv.classList.add('hidden');
    },

    _toggle3DSpin() {
        if (!this._viewer3D) return;
        this._spinning3D = !this._spinning3D;
        if (this._spinning3D) {
            this._viewer3D.spin('y', 0.8);
        } else {
            this._viewer3D.spin(false);
        }
        const btn = document.getElementById('btn-3d-spin');
        if (btn) btn.textContent = this._spinning3D ? '\u23F8 Pause' : '\u25B6 Play';
    },

    _setStyle3D(style) {
        this._current3DStyle = style;
        if (this._viewer3D) {
            this._applyStyle3D(style);
            this._viewer3D.render();
        }
        // Update button active states
        ['cartoon', 'stick', 'sphere', 'line', 'surface'].forEach(s => {
            const btn = document.getElementById(`btn-style-${s}`);
            if (btn) {
                if (s === style) {
                    btn.className = 'px-3 py-1.5 bg-emerald-500 text-white text-xs rounded-lg transition';
                } else {
                    btn.className = 'px-3 py-1.5 bg-slate-200 dark:bg-slate-700 text-xs rounded-lg hover:bg-slate-300 transition';
                }
            }
        });
    },

    _applyStyle3D(style) {
        if (!this._viewer3D) return;
        this._viewer3D.setStyle({}, {});

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
            case 'ss':
                // Color by secondary structure type
                this._viewer3D.setStyle({}, { cartoon: {
                    colorfunc: function(atom) {
                        if (atom.ss === 'h') return '#10b981'; // helix = verde
                        if (atom.ss === 's') return '#3b82f6'; // sheet = azul
                        return '#f59e0b'; // coil = amarillo
                    }
                }});
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
                <p class="text-xs text-secondary mb-3">Vista lineal del cromosoma. Genes en hebra 5'&rarr;3' (+) arriba, hebra 3'&rarr;5' (-) abajo. Pasa el cursor sobre un gen para ver detalles.</p>
                <div class="flex flex-wrap gap-2 mb-3 items-center">
                    <button onclick="DashboardRenderer._zoomGenomeMap(1.5)" class="px-3 py-1 bg-slate-100 dark:bg-slate-800 rounded text-xs hover:bg-slate-200 transition">Zoom +</button>
                    <button onclick="DashboardRenderer._zoomGenomeMap(0.67)" class="px-3 py-1 bg-slate-100 dark:bg-slate-800 rounded text-xs hover:bg-slate-200 transition">Zoom -</button>
                    <span class="text-secondary text-xs">|</span>
                    <button onclick="DashboardRenderer._zoomGenomeMapTo(3)" class="px-2 py-1 bg-slate-100 dark:bg-slate-800 rounded text-xs hover:bg-slate-200 transition">300%</button>
                    <button onclick="DashboardRenderer._zoomGenomeMapTo(10)" class="px-2 py-1 bg-slate-100 dark:bg-slate-800 rounded text-xs hover:bg-slate-200 transition">1000%</button>
                    <button onclick="DashboardRenderer._zoomGenomeMapTo(20)" class="px-2 py-1 bg-slate-100 dark:bg-slate-800 rounded text-xs hover:bg-slate-200 transition">2000%</button>
                    <button onclick="DashboardRenderer._zoomGenomeMapTo(30)" class="px-2 py-1 bg-slate-100 dark:bg-slate-800 rounded text-xs hover:bg-slate-200 transition">3000%</button>
                    <span class="text-secondary text-xs">|</span>
                    <button onclick="DashboardRenderer._zoomGenomeMap(0)" class="px-3 py-1 bg-slate-100 dark:bg-slate-800 rounded text-xs hover:bg-slate-200 transition">Reiniciar</button>
                    <span class="text-xs text-secondary ml-2 self-center" id="genome-map-zoom-label">100%</span>
                </div>
                <div id="genome-map-container" class="overflow-x-auto border border-slate-100 dark:border-slate-700 rounded-lg" style="min-height: 180px; position:relative;">
                    <p class="text-center text-secondary text-sm py-8">Cargando mapa del genoma...</p>
                </div>
                <div id="genome-map-tooltip" style="display:none; position:fixed; z-index:100; pointer-events:none; max-width:320px;" class="bg-slate-900 text-white text-xs rounded-lg px-3 py-2 shadow-xl border border-slate-700"></div>
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
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-center">${op.hebra === '+' ? "5'\u21923'" : "3'\u21925'"}</td>
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
                                        <p class="text-xs text-secondary mt-2">Posicion: ${this.fmt(op.inicio)} - ${this.fmt(op.fin)} | Hebra: ${op.hebra === '+' ? "5'\u21923' (+)" : "3'\u21925' (-)"} | ${op.num_genes} genes en ${this.fmt(op.longitud_total)} pb</p>
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

        // Viewport culling: solo renderizar genes visibles + margen
        const scrollLeft = container.scrollLeft || 0;
        const viewW = container.clientWidth || 800;
        const visibleLeft = scrollLeft - 200;
        const visibleRight = scrollLeft + viewW + 200;
        const usableW = width - 2 * margin;

        let svg = `<svg width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" xmlns="http://www.w3.org/2000/svg" style="min-width: ${width}px">`;
        svg += `<rect width="${width}" height="${height}" fill="transparent"/>`;
        svg += `<line x1="${margin}" y1="${trackCenter}" x2="${width - margin}" y2="${trackCenter}" stroke="#475569" stroke-width="2" stroke-dasharray="4,2"/>`;

        // Marcas de escala
        const scaleStep = Math.pow(10, Math.floor(Math.log10(genomeLength / 10)));
        for (let pos = 0; pos <= genomeLength; pos += scaleStep) {
            const x = margin + (pos / genomeLength) * usableW;
            if (x < visibleLeft || x > visibleRight) continue;
            svg += `<line x1="${x}" y1="${trackCenter - 5}" x2="${x}" y2="${trackCenter + 5}" stroke="#64748b" stroke-width="1"/>`;
            if (pos % (scaleStep * 2) === 0) {
                const label = pos >= 1000000 ? (pos / 1000000).toFixed(1) + ' Mb' : pos >= 1000 ? (pos / 1000).toFixed(0) + ' kb' : pos;
                svg += `<text x="${x}" y="${height - 5}" text-anchor="middle" fill="#94a3b8" font-size="9">${label}</text>`;
            }
        }

        // Genes como rectangulos con data-idx para tooltip
        const fmtFn = this.fmt.bind(this);
        for (let i = 0; i < genes.length; i++) {
            const gene = genes[i];
            const inicio = parseInt(gene.inicio || gene.start || 0);
            const fin = parseInt(gene.fin || gene.end || 0);
            const hebra = gene.hebra || gene.strand || '+';

            const x = margin + (inicio / genomeLength) * usableW;
            const w = Math.max(1, ((fin - inicio) / genomeLength) * usableW);
            // Viewport culling
            if (x + w < visibleLeft || x > visibleRight) continue;

            const y = hebra === '+' || hebra === '1' ? trackFwd : trackRev;
            const color = hebra === '+' || hebra === '1' ? '#10b981' : '#06b6d4';

            svg += `<rect class="gm-gene" data-idx="${i}" x="${x}" y="${y}" width="${w}" height="${geneH}" fill="${color}" opacity="0.6" rx="1" style="cursor:pointer"/>`;
        }

        // Labels de hebras
        svg += `<text x="${margin}" y="25" fill="#10b981" font-size="11" font-weight="bold">Hebra 5'\u21923' (+) - ${genes.filter(g => (g.hebra || g.strand) === '+' || (g.hebra || g.strand) === '1').length} genes</text>`;
        svg += `<text x="${margin}" y="${height - 15}" fill="#06b6d4" font-size="11" font-weight="bold">Hebra 3'\u21925' (-) - ${genes.filter(g => (g.hebra || g.strand) === '-' || (g.hebra || g.strand) === '-1').length} genes</text>`;

        svg += '</svg>';
        container.innerHTML = svg;

        // Custom tooltip via event delegation
        const tooltip = document.getElementById('genome-map-tooltip');
        if (!tooltip) return;
        const self = this;
        container.addEventListener('mousemove', function(e) {
            const target = e.target.closest('.gm-gene');
            if (target && target.dataset.idx != null) {
                const g = self._lastGenomeMapGenes[parseInt(target.dataset.idx)];
                if (!g) return;
                const nombre = g.nombre_gen || g.gene || g.locus_tag || 'Sin nombre';
                const producto = g.producto || 'Sin anotacion';
                const inicio = parseInt(g.inicio || g.start || 0);
                const fin = parseInt(g.fin || g.end || 0);
                const hebra = g.hebra || g.strand || '+';
                const hebraLabel = (hebra === '+' || hebra === '1') ? "5'\u21923' (+)" : "3'\u21925' (-)";
                tooltip.innerHTML = '<strong>' + nombre + '</strong><br>' + producto + '<br><span style="opacity:0.7">Pos: ' + fmtFn(inicio) + ' - ' + fmtFn(fin) + ' | ' + fmtFn(fin - inicio) + ' pb</span><br><span style="opacity:0.7">Hebra: ' + hebraLabel + '</span>';
                tooltip.style.display = 'block';
                tooltip.style.left = (e.clientX + 12) + 'px';
                tooltip.style.top = (e.clientY - 10) + 'px';
            } else {
                tooltip.style.display = 'none';
            }
        });
        container.addEventListener('mouseleave', function() {
            tooltip.style.display = 'none';
        });
    },

    _zoomGenomeMap(factor) {
        const container = document.getElementById('genome-map-container');
        const label = document.getElementById('genome-map-zoom-label');
        if (!container) return;

        if (factor === 0) {
            this._genomeMapZoom = 1;
        } else {
            this._genomeMapZoom = Math.max(0.5, Math.min(30, this._genomeMapZoom * factor));
        }

        if (label) label.textContent = Math.round(this._genomeMapZoom * 100) + '%';

        if (this._lastGenomeMapGenes && this._lastGenomeMapLength) {
            if (this._genomeMapRenderTimer) cancelAnimationFrame(this._genomeMapRenderTimer);
            this._genomeMapRenderTimer = requestAnimationFrame(() => {
                this._renderGenomeMapSVG(this._lastGenomeMapGenes, this._lastGenomeMapLength, container);
            });
        }
    },

    _zoomGenomeMapTo(level) {
        const container = document.getElementById('genome-map-container');
        const label = document.getElementById('genome-map-zoom-label');
        if (!container) return;
        this._genomeMapZoom = Math.max(0.5, Math.min(30, level));
        if (label) label.textContent = Math.round(this._genomeMapZoom * 100) + '%';
        if (this._lastGenomeMapGenes && this._lastGenomeMapLength) {
            this._renderGenomeMapSVG(this._lastGenomeMapGenes, this._lastGenomeMapLength, container);
        }
    },

    // =========================================================================
    // BUSQUEDA DE SECUENCIA
    // =========================================================================

    renderBusquedaSecuencia(genome, container) {
        // Obtener total de genes para el input
        this.fetchResultData(genome, `analisis_genes_${genome}.json`).then(genesData => {
            const totalGenes = genesData?.estadisticas_generales?.total_genes || '?';
            container.innerHTML = `
            <div class="max-w-4xl mx-auto">
                <!-- BUSCAR POR NUMERO DE GEN -->
                <div class="bg-card rounded-xl p-6 border border-slate-200 mb-6">
                    <h3 class="text-lg font-bold text-primary mb-2">Buscar Gen por Numero</h3>
                    <p class="text-sm text-secondary mb-4">Ingresa el numero del gen (1 a ${totalGenes}) para ver su informacion completa: nombre, producto, posicion, hebra y secuencia.</p>

                    <div class="flex gap-3 items-end">
                        <div class="flex-1">
                            <label class="block text-xs font-medium text-secondary mb-1">Numero de Gen</label>
                            <input
                                type="number"
                                id="buscar-gen-numero-input"
                                min="1" max="${totalGenes}" value="1"
                                placeholder="Ej: 5"
                                class="w-full px-4 py-3 bg-slate-50 dark:bg-slate-800 rounded-xl border border-slate-200 dark:border-slate-700 focus:border-emerald-500 focus:ring-2 focus:ring-emerald-500/20 outline-none text-sm text-primary"
                                onkeydown="if(event.key==='Enter') DashboardRenderer._buscarGenPorNumero('${genome}')"
                            >
                        </div>
                        <button
                            onclick="DashboardRenderer._buscarGenPorNumero('${genome}')"
                            class="px-6 py-3 bg-gradient-to-r from-violet-500 to-purple-500 hover:from-violet-600 hover:to-purple-600 text-white rounded-xl font-medium text-sm transition shadow-lg"
                        >
                            Buscar Gen
                        </button>
                    </div>
                    <p class="text-xs text-secondary mt-2">Total de genes en el genoma: <strong class="text-primary">${totalGenes}</strong>. Genes ordenados por posicion genomica.</p>
                </div>

                <!-- Resultado del gen -->
                <div id="buscar-gen-numero-resultado"></div>

                <!-- BUSCAR POR SECUENCIA -->
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
                    <p class="text-xs text-secondary mt-2">Se busca en ambas hebras: 5'&rarr;3' (+) y 3'&rarr;5' (-). Maximo 200 resultados.</p>
                </div>

                <!-- Resultados -->
                <div id="buscar-secuencia-resultados"></div>
            </div>
            `;
        });
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
                        <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-center font-bold ${hebraColor}">${m.hebra === '+' ? "5'\u21923'" : "3'\u21925'"}</td>
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

    async _buscarGenPorNumero(genome) {
        const input = document.getElementById('buscar-gen-numero-input');
        const resultContainer = document.getElementById('buscar-gen-numero-resultado');
        const numGen = parseInt(input?.value);

        if (!numGen || numGen < 1) {
            resultContainer.innerHTML = '<p class="text-amber-500 text-sm text-center py-4">Ingresa un numero de gen valido (>= 1)</p>';
            return;
        }

        resultContainer.innerHTML = `
            <div class="text-center py-6">
                <div class="w-8 h-8 border-2 border-violet-500 border-t-transparent rounded-full animate-spin mx-auto mb-3"></div>
                <p class="text-secondary text-sm">Buscando gen #${numGen}...</p>
            </div>
        `;

        try {
            // Cargar lista de genes del CSV
            const csvFile = 'lista_genes_' + genome + '.csv';
            const resp = await fetch('/api/result_data?genome=' + genome + '&file=' + csvFile);
            const result = await resp.json();

            if (!result.success || !result.data || result.data.length === 0) {
                resultContainer.innerHTML = '<p class="text-red-500 text-sm text-center py-4">No hay datos de genes. Ejecuta el analisis de genes primero.</p>';
                return;
            }

            const genes = result.data;
            if (numGen > genes.length) {
                resultContainer.innerHTML = '<p class="text-amber-500 text-sm text-center py-4">El gen #' + numGen + ' no existe. El genoma tiene ' + genes.length + ' genes (1-' + genes.length + ').</p>';
                return;
            }

            const g = genes[numGen - 1];
            const nombre = g.nombre_gen || g.gene || '-';
            const locus = g.locus_tag || '-';
            const producto = g.producto || 'Sin anotacion';
            const inicio = parseInt(g.inicio || g.start || 0);
            const fin = parseInt(g.fin || g.end || 0);
            const hebra = g.hebra || g.strand || '+';
            const longitud = fin - inicio;
            const gc = g.contenido_gc || g.gc || '-';
            const hebraLabel = (hebra === '+' || hebra === '1') ? "5'\u21923' (+)" : "3'\u21925' (-)";
            const hebraColor = (hebra === '+' || hebra === '1') ? 'emerald' : 'cyan';

            // Calcular posicion relativa en el genoma
            const genomeLen = this._lastGenomeMapLength || 4641652;
            const posRelativa = ((inicio / genomeLen) * 100).toFixed(2);

            // Mini mapa lineal mostrando posicion del gen
            const mapW = 600;
            const mapH = 40;
            const gx = (inicio / genomeLen) * mapW;
            const gw = Math.max(3, (longitud / genomeLen) * mapW);

            let html = `
                <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                    <div class="flex items-center justify-between mb-4">
                        <h4 class="text-sm font-semibold text-primary">Gen #${numGen} de ${genes.length}</h4>
                        <span class="px-3 py-1 bg-${hebraColor}-500/10 text-${hebraColor}-500 rounded-full text-xs font-medium">Hebra ${hebraLabel}</span>
                    </div>

                    <!-- Info cards -->
                    <div class="grid grid-cols-2 md:grid-cols-4 gap-3 mb-4">
                        ${this.statsCard('Gen', nombre, locus)}
                        ${this.statsCard('Longitud', this.fmt(longitud) + ' pb', 'Pos: ' + this.fmt(inicio) + '-' + this.fmt(fin), 'cyan')}
                        ${this.statsCard('GC%', gc !== '-' ? (parseFloat(gc) * 100).toFixed(1) + '%' : gc, 'Contenido GC', 'amber')}
                        ${this.statsCard('Posicion', posRelativa + '%', 'del genoma total', 'violet')}
                    </div>

                    <!-- Producto -->
                    <div class="mb-4 p-3 bg-slate-50 dark:bg-slate-800 rounded-lg">
                        <p class="text-xs font-medium text-secondary mb-1">Producto / Funcion</p>
                        <p class="text-sm text-primary">${producto}</p>
                    </div>

                    <!-- Mini mapa de posicion -->
                    <div class="mb-4">
                        <p class="text-xs font-medium text-secondary mb-2">Ubicacion en el genoma</p>
                        <svg width="100%" height="${mapH}" viewBox="0 0 ${mapW} ${mapH}" style="max-width:${mapW}px">
                            <rect width="${mapW}" height="${mapH}" fill="#1e293b" rx="4"/>
                            <line x1="0" y1="20" x2="${mapW}" y2="20" stroke="#475569" stroke-width="1"/>
                            <rect x="${gx}" y="4" width="${gw}" height="32" fill="${hebraColor === 'emerald' ? '#10b981' : '#06b6d4'}" rx="2" opacity="0.9"/>
                            <text x="${Math.min(gx + gw + 5, mapW - 60)}" y="25" fill="#e2e8f0" font-size="10">${nombre !== '-' ? nombre : locus}</text>
                        </svg>
                    </div>

                    <!-- Acciones -->
                    <div class="flex gap-2">
                        <button onclick="DashboardRenderer._extraerGenIndividual('${genome}', ${numGen})"
                            class="px-4 py-2 bg-gradient-to-r from-emerald-500 to-cyan-500 hover:from-emerald-600 hover:to-cyan-600 text-white text-xs rounded-lg transition shadow-lg">
                            Extraer Secuencia de ADN
                        </button>
                        <button onclick="DashboardRenderer._verGenEnMapa('${genome}', ${inicio}, ${fin})"
                            class="px-4 py-2 bg-violet-500/10 hover:bg-violet-500/20 text-violet-500 text-xs font-medium rounded-lg transition">
                            Ver en Mapa del Genoma
                        </button>
                    </div>
                </div>
            `;

            // Genes vecinos (anterior y siguiente)
            const vecinos = [];
            if (numGen > 1) {
                const prev = genes[numGen - 2];
                vecinos.push({label: 'Gen anterior (#' + (numGen - 1) + ')', gen: prev});
            }
            if (numGen < genes.length) {
                const next = genes[numGen];
                vecinos.push({label: 'Gen siguiente (#' + (numGen + 1) + ')', gen: next});
            }
            if (vecinos.length > 0) {
                html += '<div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">';
                html += '<h4 class="text-sm font-semibold text-primary mb-3">Genes Vecinos</h4>';
                html += '<div class="grid grid-cols-1 md:grid-cols-2 gap-3">';
                for (const v of vecinos) {
                    const vn = v.gen.nombre_gen || v.gen.gene || v.gen.locus_tag || '-';
                    const vp = v.gen.producto || 'Sin anotacion';
                    const vh = v.gen.hebra || v.gen.strand || '+';
                    const vhColor = (vh === '+' || vh === '1') ? 'text-emerald-500' : 'text-cyan-500';
                    html += '<div class="p-3 bg-slate-50 dark:bg-slate-800 rounded-lg">';
                    html += '<p class="text-xs text-secondary">' + v.label + '</p>';
                    html += '<p class="text-sm font-medium text-primary">' + vn + ' <span class="font-bold ' + vhColor + '">' + ((vh === '+' || vh === '1') ? "5'\u21923'" : "3'\u21925'") + '</span></p>';
                    html += '<p class="text-xs text-secondary truncate">' + vp + '</p>';
                    html += '</div>';
                }
                html += '</div></div>';
            }

            resultContainer.innerHTML = html;
        } catch (e) {
            resultContainer.innerHTML = '<p class="text-red-500 text-sm text-center py-4">Error: ' + e.message + '</p>';
        }
    },

    async _extraerGenIndividual(genome, numGen) {
        const resultContainer = document.getElementById('buscar-gen-numero-resultado');
        try {
            const resp = await fetch('/api/extraer_secuencia?genome=' + genome + '&gene_start=' + numGen + '&gene_end=' + numGen + '&mode=cds');
            const data = await resp.json();
            if (!data.success) {
                alert('Error: ' + data.error);
                return;
            }
            const seq = data.secuencia_total || '';
            // Copiar al portapapeles
            const fastaHeader = '>' + genome + '_gen_' + numGen + ' | ' + (data.genes_info?.[0]?.nombre_gen || '') + ' | ' + (data.genes_info?.[0]?.producto || '');
            let fasta = fastaHeader + '\n';
            for (let i = 0; i < seq.length; i += 60) {
                fasta += seq.substring(i, i + 60) + '\n';
            }
            await navigator.clipboard.writeText(fasta);

            // Mostrar secuencia en resultado
            const prev = resultContainer.innerHTML;
            const seqCard = `
                <div class="bg-card rounded-xl p-5 border border-emerald-200 mb-6">
                    <div class="flex items-center justify-between mb-3">
                        <h4 class="text-sm font-semibold text-emerald-500">Secuencia Extraida (copiada al portapapeles)</h4>
                        <span class="px-2 py-1 bg-emerald-500/10 text-emerald-500 text-xs rounded-full">${this.fmt(seq.length)} pb</span>
                    </div>
                    <pre class="text-[10px] font-mono text-secondary bg-slate-50 dark:bg-slate-900 rounded-lg p-3 max-h-[200px] overflow-auto whitespace-pre-wrap break-all">${fasta}</pre>
                </div>
            `;
            resultContainer.innerHTML = prev + seqCard;
        } catch (e) {
            alert('Error extrayendo secuencia: ' + e.message);
        }
    },

    _verGenEnMapa(genome, inicio, fin) {
        // Navegar a la seccion de estructura del gen y hacer scroll al mapa
        const mapContainer = document.getElementById('genome-map-container');
        if (!mapContainer) {
            alert('Primero ve a la seccion "Estructura del Gen" y ejecuta el analisis para ver el mapa del genoma.');
            return;
        }
        // Calcular zoom necesario para ver el gen
        const genomeLen = this._lastGenomeMapLength || 4641652;
        const genLen = fin - inicio;
        const targetZoom = Math.max(1, Math.min(30, (genomeLen / genLen) * 0.05));
        this._genomeMapZoom = targetZoom;
        const label = document.getElementById('genome-map-zoom-label');
        if (label) label.textContent = Math.round(targetZoom * 100) + '%';
        if (this._lastGenomeMapGenes && this._lastGenomeMapLength) {
            this._renderGenomeMapSVG(this._lastGenomeMapGenes, this._lastGenomeMapLength, mapContainer);
        }
        // Scroll to gene position
        const width = Math.max(mapContainer.clientWidth - 20, 800) * targetZoom;
        const scrollTo = (inicio / genomeLen) * width - mapContainer.clientWidth / 2;
        mapContainer.scrollLeft = Math.max(0, scrollTo);
        mapContainer.scrollIntoView({ behavior: 'smooth', block: 'center' });
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
                        <td class="px-2 py-2 text-center font-bold ${g.hebra === '+' ? 'text-emerald-500' : 'text-cyan-500'}">${g.hebra === '+' ? "5'\u21923'" : "3'\u21925'"}</td>
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

    // =========================================================================
    // DASHBOARD DE EVOLUCION - GENES COMPARTIDOS (PANGENOMA)
    // =========================================================================

    renderEvolucionPangenoma(data, container) {
        const pan = data.pangenoma || {};
        const curva = data.curva_pangenoma || [];
        const genomas = data.genomas_analizados || [];
        const genesDetalle = pan.genes_unicos_detalle || {};
        const genesPorGenoma = pan.genes_por_genoma || {};
        const totalConocidos = pan.total_productos_conocidos || pan.total_genes_unicos || 0;
        const corePct = totalConocidos > 0 ? ((pan.core_genome || 0) / totalConocidos * 100).toFixed(1) : '0';
        const accPct = totalConocidos > 0 ? ((pan.accessory_genome || 0) / totalConocidos * 100).toFixed(1) : '0';
        const uniqPct = totalConocidos > 0 ? ((pan.genes_unicos_total || 0) / totalConocidos * 100).toFixed(1) : '0';

        container.innerHTML = `
            <!-- Explicacion general -->
            <div class="bg-emerald-50 dark:bg-emerald-900/15 border border-emerald-200 rounded-xl p-5 mb-6">
                <h3 class="text-lg font-bold text-emerald-700 mb-2">Genes Compartidos entre tus ${data.total_genomas || 0} Bacterias</h3>
                <p class="text-sm text-secondary leading-relaxed">
                    Analizamos <strong>todos los genes</strong> de cada bacteria y los comparamos. Encontramos <strong>${this.fmt(totalConocidos)} tipos de genes diferentes</strong> en total.
                    De esos, <strong class="text-cyan-600">${this.fmt(pan.core_genome || 0)} genes (${corePct}%) son compartidos por TODAS</strong> las bacterias (son esenciales para vivir),
                    <strong class="text-amber-600">${this.fmt(pan.accessory_genome || 0)} (${accPct}%) los tienen algunas pero no todas</strong> (dan ventajas especiales),
                    y <strong class="text-red-600">${this.fmt(pan.genes_unicos_total || 0)} (${uniqPct}%) son unicos</strong> de una sola bacteria (la hacen diferente al resto).
                </p>
            </div>

            <!-- Stat cards interactivos -->
            <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
                <div class="bg-card rounded-xl p-5 border border-slate-200 evo-hover-card" title="Todos los tipos de genes diferentes encontrados entre todas tus bacterias">
                    <p class="text-xs font-medium text-secondary uppercase tracking-wide">Total de Genes</p>
                    <p class="text-2xl font-bold text-emerald-500 mt-1">${this.fmt(totalConocidos)}</p>
                    <p class="text-xs text-secondary mt-1">tipos distintos encontrados</p>
                </div>
                <div class="bg-card rounded-xl p-5 border border-slate-200 evo-hover-card" title="Genes que TODAS las bacterias comparten. Son esenciales para la vida bacteriana">
                    <p class="text-xs font-medium text-secondary uppercase tracking-wide">Genes en Comun</p>
                    <p class="text-2xl font-bold text-cyan-500 mt-1">${this.fmt(pan.core_genome || 0)}</p>
                    <p class="text-xs text-secondary mt-1">compartidos por TODAS (${corePct}%)</p>
                </div>
                <div class="bg-card rounded-xl p-5 border border-slate-200 evo-hover-card" title="Genes que tienen ALGUNAS bacterias pero no todas. Dan habilidades especiales como resistencia">
                    <p class="text-xs font-medium text-secondary uppercase tracking-wide">Genes Parciales</p>
                    <p class="text-2xl font-bold text-amber-500 mt-1">${this.fmt(pan.accessory_genome || 0)}</p>
                    <p class="text-xs text-secondary mt-1">en algunas, no en todas (${accPct}%)</p>
                </div>
                <div class="bg-card rounded-xl p-5 border border-slate-200 evo-hover-card" title="Genes que SOLO tiene una bacteria. La hacen unica respecto a las demas">
                    <p class="text-xs font-medium text-secondary uppercase tracking-wide">Genes Exclusivos</p>
                    <p class="text-2xl font-bold text-red-500 mt-1">${this.fmt(pan.genes_unicos_total || 0)}</p>
                    <p class="text-xs text-secondary mt-1">unicos de 1 sola bacteria (${uniqPct}%)</p>
                </div>
            </div>

            <!-- Charts row -->
            <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h4 class="text-sm font-semibold text-primary mb-1">Como se reparten los genes</h4>
                    <p class="text-xs text-secondary mb-3">La "torta" muestra cuantos genes comparten todas, algunas o solo una bacteria</p>
                    <canvas id="chart-evo-donut" height="280"></canvas>
                </div>
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h4 class="text-sm font-semibold text-primary mb-1">Que pasa al agregar mas bacterias?</h4>
                    <p class="text-xs text-secondary mb-3">Al agregar mas bacterias se descubren mas genes (verde sube), pero los compartidos por TODAS bajan (azul baja)</p>
                    <canvas id="chart-evo-curva" height="280"></canvas>
                </div>
            </div>

            <!-- Tabla detalle por bacteria -->
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h4 class="text-sm font-semibold text-primary mb-1">Detalle por Bacteria</h4>
                <p class="text-xs text-secondary mb-3">Pasa el mouse sobre cada fila para ver mas informacion</p>
                <div class="overflow-x-auto">
                    <table class="w-full text-xs">
                        <thead class="bg-slate-50 dark:bg-slate-800">
                            <tr>
                                <th class="px-3 py-2 text-left text-secondary">Bacteria</th>
                                <th class="px-3 py-2 text-right text-secondary">Total Genes</th>
                                <th class="px-3 py-2 text-right text-secondary" title="Genes que comparte con TODAS las demas">En Comun</th>
                                <th class="px-3 py-2 text-right text-secondary" title="Genes que solo tiene esta bacteria">Exclusivos</th>
                                <th class="px-3 py-2 text-right text-secondary" title="Porcentaje de guanina-citosina en el ADN">GC%</th>
                                <th class="px-3 py-2 text-right text-secondary" title="Tamano del genoma en millones de pares de bases">Tamano</th>
                            </tr>
                        </thead>
                        <tbody>
                            ${genomas.map(g => {
                                const gpg = genesPorGenoma[g.basename] || {};
                                const pctCore = g.total_genes > 0 ? ((gpg.core || 0) / g.total_genes * 100).toFixed(0) : 0;
                                return `
                                    <tr class="border-t border-slate-100 dark:border-slate-800 hover:bg-emerald-50 dark:hover:bg-emerald-900/10 cursor-default"
                                        title="Bacteria: ${(g.nombre || g.basename).substring(0, 80)} | ${gpg.core || 0} genes en comun (${pctCore}% de sus genes) | ${gpg.unicos || 0} exclusivos | ${gpg.hypothetical || 0} sin identificar">
                                        <td class="px-3 py-2.5 font-medium text-primary">${(g.nombre || g.basename).replace(/_/g, ' ').substring(0, 50)}</td>
                                        <td class="px-3 py-2.5 text-right text-primary">${this.fmt(g.total_genes)}</td>
                                        <td class="px-3 py-2.5 text-right"><span class="text-cyan-500 font-bold">${this.fmt(gpg.core || 0)}</span> <span class="text-secondary">(${pctCore}%)</span></td>
                                        <td class="px-3 py-2.5 text-right text-red-500 font-medium">${this.fmt(gpg.unicos || 0)}</td>
                                        <td class="px-3 py-2.5 text-right text-primary">${(g.gc_porcentaje || 0).toFixed(1)}%</td>
                                        <td class="px-3 py-2.5 text-right text-primary">${((g.longitud_pb || 0) / 1e6).toFixed(2)} Mb</td>
                                    </tr>
                                `;
                            }).join('')}
                        </tbody>
                    </table>
                </div>
            </div>

            <!-- Genes exclusivos expandibles -->
            ${Object.keys(genesDetalle).length > 0 ? `
            <div class="bg-card rounded-xl p-5 border border-slate-200">
                <h4 class="text-sm font-semibold text-primary mb-1">Genes Exclusivos de cada Bacteria</h4>
                <p class="text-xs text-secondary mb-3">Haz click en cada bacteria para ver que genes tiene que las demas NO tienen. Estos genes pueden darle habilidades especiales (resistencia a antibioticos, produccion de toxinas, etc.)</p>
                <div class="space-y-2 max-h-[400px] overflow-y-auto">
                    ${Object.entries(genesDetalle).map(([basename, genes]) => `
                        <div class="bg-slate-50 dark:bg-slate-800 rounded-lg p-3 cursor-pointer hover:bg-slate-100 dark:hover:bg-slate-700 transition"
                             onclick="this.nextElementSibling.classList.toggle('hidden')">
                            <div class="flex items-center justify-between">
                                <span class="font-medium text-sm text-primary">${basename.replace(/_/g, ' ').substring(0, 50)}</span>
                                <span class="text-xs px-2 py-0.5 bg-red-500/10 text-red-500 rounded-full">${genes.length} genes exclusivos - click para ver</span>
                            </div>
                        </div>
                        <div class="hidden px-3 pb-2">
                            <div class="bg-white dark:bg-slate-900 rounded-lg p-3 border border-slate-200 text-xs max-h-[200px] overflow-y-auto">
                                ${genes.slice(0, 50).map(g => `
                                    <div class="py-1.5 border-b border-slate-50 dark:border-slate-800 flex justify-between gap-2 hover:bg-emerald-50 dark:hover:bg-emerald-900/10 px-1 rounded"
                                         title="Gen: ${g.nombre_gen || 'sin nombre'} | Codigo: ${g.locus_tag || ''} | Tamano: ${g.longitud_pb || 0} pares de bases | Funcion: ${g.producto || 'desconocida'}">
                                        <span class="text-primary flex-1">${g.producto || 'Funcion desconocida'}</span>
                                        <span class="text-secondary font-mono">${g.locus_tag || ''}</span>
                                    </div>
                                `).join('')}
                                ${genes.length > 50 ? `<p class="text-secondary mt-2 text-center">... y ${genes.length - 50} genes exclusivos mas</p>` : ''}
                            </div>
                        </div>
                    `).join('')}
                </div>
            </div>` : ''}
        `;

        setTimeout(() => {
            this.createChart('chart-evo-donut', {
                type: 'doughnut',
                data: {
                    labels: ['Compartidos por TODAS', 'En algunas (no todas)', 'Exclusivos de una sola'],
                    datasets: [{
                        data: [pan.core_genome || 0, pan.accessory_genome || 0, pan.genes_unicos_total || 0],
                        backgroundColor: ['#06b6d4', '#f59e0b', '#ef4444'],
                        borderWidth: 2,
                        borderColor: document.body.classList.contains('dark') ? '#1a1a1a' : '#ffffff'
                    }]
                },
                options: {
                    responsive: true,
                    plugins: {
                        legend: { position: 'bottom', labels: { padding: 15, font: { size: 11 } } },
                        tooltip: {
                            callbacks: {
                                label: (ctx) => {
                                    const total = totalConocidos || 1;
                                    const pct = ((ctx.raw / total) * 100).toFixed(1);
                                    return `${ctx.label}: ${this.fmt(ctx.raw)} genes (${pct}%)`;
                                }
                            }
                        }
                    }
                }
            });

            if (curva.length > 0) {
                this.createChart('chart-evo-curva', {
                    type: 'line',
                    data: {
                        labels: curva.map(c => c.n_genomas + ' bacterias'),
                        datasets: [
                            {
                                label: 'Total genes descubiertos',
                                data: curva.map(c => c.pan),
                                borderColor: '#10b981',
                                backgroundColor: 'rgba(16,185,129,0.1)',
                                fill: true, tension: 0.3, pointRadius: 5,
                                pointHoverRadius: 8
                            },
                            {
                                label: 'Genes compartidos por TODAS',
                                data: curva.map(c => c.core),
                                borderColor: '#06b6d4',
                                backgroundColor: 'rgba(6,182,212,0.1)',
                                fill: true, tension: 0.3, pointRadius: 5,
                                pointHoverRadius: 8
                            }
                        ]
                    },
                    options: {
                        responsive: true,
                        plugins: {
                            legend: { position: 'bottom' },
                            tooltip: {
                                callbacks: {
                                    label: (ctx) => `${ctx.dataset.label}: ${this.fmt(ctx.raw)} genes`
                                }
                            }
                        },
                        scales: {
                            x: { title: { display: true, text: 'Cantidad de bacterias analizadas' } },
                            y: { title: { display: true, text: 'Cantidad de genes' }, beginAtZero: true }
                        }
                    }
                });
            }
        }, 100);
    },

    // =========================================================================
    // DASHBOARD DE EVOLUCION - ARBOL FAMILIAR
    // =========================================================================

    renderEvolucionArbol(data, container) {
        const arbol = data.arbol_upgma || {};
        const nodos = arbol.nodos || [];
        const genomas = data.genomas_analizados || [];
        const dist = data.distancias_jaccard || {};
        const matriz = dist.matriz || [];
        const raiz = arbol.raiz || 0;

        // Calcular pares mas cercanos y lejanos
        let parCercano = { i: '', j: '', d: 999 };
        let parLejano = { i: '', j: '', d: 0 };
        const dGenomas = dist.genomas || [];
        for (let i = 0; i < dGenomas.length; i++) {
            for (let j = i + 1; j < dGenomas.length; j++) {
                const v = matriz[i]?.[j] || 0;
                if (v < parCercano.d) parCercano = { i: dGenomas[i], j: dGenomas[j], d: v };
                if (v > parLejano.d) parLejano = { i: dGenomas[i], j: dGenomas[j], d: v };
            }
        }
        const cercSimil = ((1 - parCercano.d) * 100).toFixed(0);
        const lejSimil = ((1 - parLejano.d) * 100).toFixed(0);

        container.innerHTML = `
            <!-- Explicacion -->
            <div class="bg-violet-50 dark:bg-violet-900/15 border border-violet-200 rounded-xl p-5 mb-6">
                <h3 class="text-lg font-bold text-violet-700 mb-2">Arbol Familiar de tus Bacterias</h3>
                <p class="text-sm text-secondary leading-relaxed">
                    Este arbol muestra <strong>que tan "parientes" son</strong> tus bacterias entre si, basandose en los genes que comparten.
                    Las bacterias que comparten mas genes estan mas cerca en el arbol (como hermanos).
                    Las que comparten menos genes estan mas lejos (como primos lejanos).
                    <strong>No es un arbol de tiempo real</strong>, sino de similitud genetica: mientras mas parecidos son sus genes, mas "cercanas" estan.
                </p>
            </div>

            <!-- Stats con contexto -->
            <div class="grid grid-cols-1 md:grid-cols-3 gap-4 mb-6">
                <div class="bg-card rounded-xl p-5 border border-slate-200 evo-hover-card" title="Cantidad de bacterias representadas como hojas del arbol">
                    <p class="text-xs font-medium text-secondary uppercase tracking-wide">Bacterias en el Arbol</p>
                    <p class="text-2xl font-bold text-emerald-500 mt-1">${genomas.length}</p>
                    <p class="text-xs text-secondary mt-1">cada una es una "hoja" del arbol</p>
                </div>
                <div class="bg-card rounded-xl p-5 border border-emerald-200 evo-hover-card" title="${parCercano.i.replace(/_/g,' ')} y ${parCercano.j.replace(/_/g,' ')} comparten ${cercSimil}% de sus genes">
                    <p class="text-xs font-medium text-secondary uppercase tracking-wide">Par mas Parecido</p>
                    <p class="text-2xl font-bold text-cyan-500 mt-1">${cercSimil}% similares</p>
                    <p class="text-xs text-secondary mt-1 line-clamp-2">${parCercano.i.replace(/_/g,' ').substring(0,30)} y ${parCercano.j.replace(/_/g,' ').substring(0,30)}</p>
                </div>
                <div class="bg-card rounded-xl p-5 border border-red-200 evo-hover-card" title="${parLejano.i.replace(/_/g,' ')} y ${parLejano.j.replace(/_/g,' ')} solo comparten ${lejSimil}% de sus genes">
                    <p class="text-xs font-medium text-secondary uppercase tracking-wide">Par mas Diferente</p>
                    <p class="text-2xl font-bold text-red-500 mt-1">${lejSimil}% similares</p>
                    <p class="text-xs text-secondary mt-1 line-clamp-2">${parLejano.i.replace(/_/g,' ').substring(0,30)} y ${parLejano.j.replace(/_/g,' ').substring(0,30)}</p>
                </div>
            </div>

            <!-- Estimacion temporal -->
            <div class="bg-amber-50 dark:bg-amber-900/15 border border-amber-200 rounded-xl p-4 mb-6">
                <h4 class="text-sm font-semibold text-amber-700 mb-1">Sobre el tiempo de evolucion</h4>
                <p class="text-xs text-secondary leading-relaxed">
                    Las bacterias evolucionan rapido: una E. coli se reproduce cada ~20 minutos, acumulando mutaciones. Dos cepas de E. coli con ~${cercSimil}% de genes compartidos
                    pudieron divergir hace <strong>miles a decenas de miles de anos</strong>.
                    Si comparas E. coli con otra especie (ej. Salmonella) con ~${lejSimil}% similitud, la divergencia pudo ocurrir hace <strong>cientos de millones de anos</strong>.
                    Estas son estimaciones aproximadas - la velocidad real depende de la presion ambiental y la transferencia horizontal de genes.
                </p>
            </div>

            <!-- Canvas del arbol -->
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6 relative">
                <div class="flex items-center justify-between mb-3">
                    <div>
                        <h4 class="text-sm font-semibold text-primary">Arbol de Parentesco (vertical)</h4>
                        <p class="text-xs text-secondary">El ancestro comun esta arriba. Las bacterias mas parecidas se conectan mas abajo (ramas cortas). Pasa el mouse sobre cada nombre.</p>
                    </div>
                    <div class="flex items-center gap-3 text-[10px] text-secondary">
                        <span class="flex items-center gap-1"><span class="w-8 h-0.5 inline-block" style="background:#10b981;"></span> Ramas</span>
                        <span class="flex items-center gap-1"><span class="w-3 h-3 rounded-full inline-block" style="background:#10b981;"></span> Bacteria</span>
                    </div>
                </div>
                <div class="overflow-x-auto overflow-y-auto" style="position:relative; max-height:700px;">
                    <canvas id="canvas-evo-tree" width="${Math.max(900, genomas.length * 140 + 120)}" height="${Math.max(650, 500)}"></canvas>
                    <div id="evo-tree-tooltip" class="evo-tooltip hidden" style="position:absolute; top:0; left:0;"></div>
                </div>
            </div>

            <!-- Newick -->
            ${arbol.newick ? `
            <div class="bg-card rounded-xl p-5 border border-slate-200">
                <h4 class="text-sm font-semibold text-primary mb-1">Formato Newick (para expertos)</h4>
                <p class="text-xs text-secondary mb-2">Este texto codifica el arbol. Puedes copiarlo y pegarlo en herramientas como iTOL, FigTree o MEGA para verlo en 3D. Haz click para copiar.</p>
                <div class="bg-slate-50 dark:bg-slate-800 rounded-lg p-3 font-mono text-xs break-all max-h-[150px] overflow-y-auto select-all cursor-pointer hover:ring-2 hover:ring-emerald-500/30 transition"
                     onclick="navigator.clipboard.writeText(this.textContent).then(()=>showNotification('Newick copiado al portapapeles','success'))">${arbol.newick}</div>
            </div>` : ''}
        `;

        // Draw tree with hover interactivity
        setTimeout(() => {
            this._drawFamilyTree('canvas-evo-tree', nodos, raiz, genomas, data);
        }, 100);
    },

    _drawFamilyTree(canvasId, nodos, raizId, genomas, data) {
        const canvas = document.getElementById(canvasId);
        if (!canvas || nodos.length === 0) return;
        const ctx = canvas.getContext('2d');
        const dpr = window.devicePixelRatio || 1;

        const hojas = nodos.filter(n => n.hoja);
        const numHojas = hojas.length;

        // Vertical layout: wide for leaves, tall for depth + labels
        const displayW = Math.max(900, numHojas * 140 + 120);
        const displayH = Math.max(650, 500);
        canvas.width = displayW * dpr;
        canvas.height = displayH * dpr;
        canvas.style.width = displayW + 'px';
        canvas.style.height = displayH + 'px';
        ctx.scale(dpr, dpr);

        const isDark = document.body.classList.contains('dark');
        const textColor = isDark ? '#e2e8f0' : '#1e293b';
        const highlightColor = '#06b6d4';
        const branchColors = ['#10b981', '#06b6d4', '#8b5cf6', '#f59e0b', '#ec4899', '#6366f1', '#ef4444', '#14b8a6'];

        const marginLeft = 65;
        const marginRight = 50;
        const marginTop = 55;
        const marginBottom = 130; // space for rotated labels
        const treeWidth = displayW - marginLeft - marginRight;
        const treeHeight = displayH - marginTop - marginBottom;

        let maxDist = 0;
        nodos.forEach(n => { if (n.distancia > maxDist) maxDist = n.distancia; });
        if (maxDist === 0) maxDist = 1;

        const xPos = {};
        const yPos = {};
        let leafIdx = 0;
        const leafPositions = [];

        // Vertical tree: x = leaf spread horizontal, y = distance from top
        // Root (maxDist) at top, leaves (dist~0) at bottom
        function assignPositions(nid) {
            const n = nodos[nid];
            if (!n) return 0;
            if (n.hoja) {
                xPos[nid] = marginLeft + leafIdx * (treeWidth / Math.max(numHojas - 1, 1));
                yPos[nid] = marginTop + treeHeight; // leaves at bottom
                leafIdx++;
                return xPos[nid];
            }
            const childXs = (n.hijos || []).map(cid => assignPositions(cid));
            xPos[nid] = childXs.reduce((a, b) => a + b, 0) / childXs.length;
            yPos[nid] = marginTop + (1 - n.distancia / maxDist) * treeHeight;
            return xPos[nid];
        }
        assignPositions(raizId);

        // Assign colors per subtree from root
        const nodeColors = {};
        function assignColors(nid, colorIdx) {
            nodeColors[nid] = branchColors[colorIdx % branchColors.length];
            const n = nodos[nid];
            if (n && !n.hoja) {
                (n.hijos || []).forEach((cid, i) => {
                    assignColors(cid, nid === raizId ? i : colorIdx);
                });
            }
        }
        assignColors(raizId, 0);

        function drawTree(highlightId) {
            ctx.clearRect(0, 0, displayW, displayH);

            // Title at top center
            ctx.fillStyle = textColor;
            ctx.font = 'bold 13px Space Grotesk, sans-serif';
            ctx.textAlign = 'center';
            ctx.fillText('Ancestro comun (raiz)', displayW / 2, marginTop - 30);
            ctx.textAlign = 'left';

            // Horizontal grid lines with similarity %
            ctx.font = '9px JetBrains Mono, monospace';
            for (let i = 0; i <= 4; i++) {
                const y = marginTop + (treeHeight / 4) * i;
                const dist = maxDist * (1 - i / 4);
                const simPct = ((1 - dist) * 100).toFixed(0);

                ctx.strokeStyle = isDark ? 'rgba(255,255,255,0.06)' : 'rgba(0,0,0,0.06)';
                ctx.lineWidth = 0.5;
                ctx.setLineDash([3, 5]);
                ctx.beginPath();
                ctx.moveTo(marginLeft - 10, y);
                ctx.lineTo(displayW - marginRight + 10, y);
                ctx.stroke();
                ctx.setLineDash([]);

                ctx.fillStyle = isDark ? '#64748b' : '#94a3b8';
                ctx.textAlign = 'right';
                ctx.fillText(`${simPct}%`, marginLeft - 15, y + 3);
                ctx.textAlign = 'left';
            }

            // Vertical axis label
            ctx.save();
            ctx.translate(14, marginTop + treeHeight / 2);
            ctx.rotate(-Math.PI / 2);
            ctx.fillStyle = isDark ? '#64748b' : '#94a3b8';
            ctx.font = '10px Space Grotesk, sans-serif';
            ctx.textAlign = 'center';
            ctx.fillText('Similitud genetica', 0, 0);
            ctx.restore();

            function drawNode(nid) {
                const n = nodos[nid];
                if (!n) return;
                const x = xPos[nid];
                const y = yPos[nid];

                if (n.hoja) {
                    const isHighlight = nid === highlightId;
                    const r = isHighlight ? 9 : 7;
                    const color = isHighlight ? highlightColor : nodeColors[nid] || '#10b981';

                    if (isHighlight) {
                        ctx.shadowColor = highlightColor;
                        ctx.shadowBlur = 16;
                    }
                    ctx.fillStyle = color;
                    ctx.beginPath();
                    ctx.arc(x, y, r, 0, Math.PI * 2);
                    ctx.fill();
                    ctx.shadowBlur = 0;

                    // White border
                    ctx.strokeStyle = isDark ? '#1e293b' : '#ffffff';
                    ctx.lineWidth = 2;
                    ctx.beginPath();
                    ctx.arc(x, y, r, 0, Math.PI * 2);
                    ctx.stroke();

                    // Rotated label below leaf
                    ctx.save();
                    ctx.translate(x, y + r + 8);
                    ctx.rotate(Math.PI / 5);
                    ctx.fillStyle = isHighlight ? highlightColor : textColor;
                    ctx.font = isHighlight ? 'bold 11px Space Grotesk, sans-serif' : '10px Space Grotesk, sans-serif';
                    ctx.textAlign = 'left';
                    const label = (n.nombre || '').replace(/_/g, ' ').substring(0, 32);
                    ctx.fillText(label, 0, 0);
                    ctx.restore();

                    leafPositions.push({ nid, x, y, r: 18, nombre: n.nombre });
                } else {
                    const hijos = n.hijos || [];
                    const childXs = hijos.map(cid => xPos[cid]);
                    const minX = Math.min(...childXs);
                    const maxX = Math.max(...childXs);

                    // Horizontal connector at parent's y level
                    ctx.strokeStyle = nodeColors[nid] || '#10b981';
                    ctx.lineWidth = 3.5;
                    ctx.lineCap = 'round';
                    ctx.beginPath();
                    ctx.moveTo(minX, y);
                    ctx.lineTo(maxX, y);
                    ctx.stroke();

                    // Small dot at junction
                    ctx.fillStyle = nodeColors[nid] || '#10b981';
                    ctx.beginPath();
                    ctx.arc(x, y, 3.5, 0, Math.PI * 2);
                    ctx.fill();

                    hijos.forEach(cid => {
                        const cx = xPos[cid];
                        const cy = yPos[cid];

                        // Vertical branch down
                        ctx.strokeStyle = nodeColors[cid] || '#10b981';
                        ctx.lineWidth = 3;
                        ctx.lineCap = 'round';
                        ctx.beginPath();
                        ctx.moveTo(cx, y);
                        ctx.lineTo(cx, cy);
                        ctx.stroke();

                        drawNode(cid);
                    });
                }
            }

            leafPositions.length = 0;
            drawNode(raizId);

            // Vertical scale bar on right side
            const scaleLen = treeHeight * 0.2;
            ctx.strokeStyle = textColor;
            ctx.lineWidth = 1.5;
            ctx.beginPath();
            ctx.moveTo(displayW - marginRight + 25, marginTop + treeHeight);
            ctx.lineTo(displayW - marginRight + 25, marginTop + treeHeight - scaleLen);
            ctx.stroke();
            // Ticks
            ctx.beginPath();
            ctx.moveTo(displayW - marginRight + 20, marginTop + treeHeight);
            ctx.lineTo(displayW - marginRight + 30, marginTop + treeHeight);
            ctx.stroke();
            ctx.beginPath();
            ctx.moveTo(displayW - marginRight + 20, marginTop + treeHeight - scaleLen);
            ctx.lineTo(displayW - marginRight + 30, marginTop + treeHeight - scaleLen);
            ctx.stroke();
            ctx.fillStyle = textColor;
            ctx.font = '9px JetBrains Mono, monospace';
            ctx.textAlign = 'center';
            const scaleVal = (maxDist * 0.2 * 100).toFixed(0);
            ctx.fillText(`${scaleVal}%`, displayW - marginRight + 25, marginTop + treeHeight + 14);
            ctx.fillText('diferencia', displayW - marginRight + 25, marginTop + treeHeight + 24);
            ctx.textAlign = 'left';
        }

        drawTree(null);

        // Hover interactivity
        const tooltip = document.getElementById('evo-tree-tooltip');
        const genesPorGenoma = data.pangenoma?.genes_por_genoma || {};

        canvas.addEventListener('mousemove', (e) => {
            const rect = canvas.getBoundingClientRect();
            const mx = e.clientX - rect.left;
            const my = e.clientY - rect.top;

            let found = null;
            for (const leaf of leafPositions) {
                const dx = mx - leaf.x;
                const dy = my - leaf.y;
                if (dx * dx + dy * dy < leaf.r * leaf.r * 4) {
                    found = leaf;
                    break;
                }
            }

            if (found) {
                canvas.style.cursor = 'pointer';
                const gInfo = genomas.find(g => g.basename === found.nombre) || {};
                const gpg = genesPorGenoma[found.nombre] || {};
                tooltip.classList.remove('hidden');
                // Position tooltip above the leaf node
                tooltip.style.left = (found.x + 20) + 'px';
                tooltip.style.top = (found.y - 80) + 'px';
                tooltip.innerHTML = `
                    <strong>${(found.nombre || '').replace(/_/g, ' ')}</strong><br>
                    <span style="color:#06b6d4">Genes totales:</span> ${gInfo.total_genes || '?'}<br>
                    <span style="color:#10b981">Genes en comun:</span> ${gpg.core || '?'}<br>
                    <span style="color:#ef4444">Genes exclusivos:</span> ${gpg.unicos || '?'}<br>
                    <span style="color:#94a3b8">GC:</span> ${gInfo.gc_porcentaje || '?'}% | <span style="color:#94a3b8">Tamano:</span> ${((gInfo.longitud_pb || 0) / 1e6).toFixed(2)} Mb
                `;
                drawTree(found.nid);
            } else {
                canvas.style.cursor = 'default';
                tooltip.classList.add('hidden');
                drawTree(null);
            }
        });

        canvas.addEventListener('mouseleave', () => {
            tooltip.classList.add('hidden');
            drawTree(null);
        });
    },

    // =========================================================================
    // DASHBOARD DE EVOLUCION - MAPA DE PARENTESCO (MATRIZ)
    // =========================================================================

    renderEvolucionMatriz(data, container) {
        const dist = data.distancias_jaccard || {};
        const matrizGenomas = dist.genomas || [];
        const matriz = dist.matriz || [];
        const presencia = data.matriz_presencia || {};

        // Store presencia for click-to-expand
        DashboardRenderer._evoPresencia = presencia;

        // Build all pairs sorted, filter out 0% similarity
        const pares = [];
        const paresCero = [];
        for (let i = 0; i < matrizGenomas.length; i++) {
            for (let j = i + 1; j < matrizGenomas.length; j++) {
                const d = matriz[i]?.[j] || 0;
                const sim = ((1 - d) * 100);
                const par = {
                    a: matrizGenomas[i],
                    b: matrizGenomas[j],
                    distancia: d,
                    similitud: sim.toFixed(1)
                };
                if (sim > 0.5) {
                    pares.push(par);
                } else {
                    paresCero.push(par);
                }
            }
        }
        pares.sort((a, b) => a.distancia - b.distancia);

        container.innerHTML = `
            <!-- Explicacion -->
            <div class="bg-cyan-50 dark:bg-cyan-900/15 border border-cyan-200 rounded-xl p-5 mb-6">
                <h3 class="text-lg font-bold text-cyan-700 mb-2">Mapa de Parentesco entre Bacterias</h3>
                <p class="text-sm text-secondary leading-relaxed">
                    Este mapa muestra <strong>que tan parecida es cada bacteria con cada otra</strong>, basandose en cuantos genes comparten.
                    El color <strong class="text-emerald-600">verde</strong> significa "muy parecidas" (comparten muchos genes, como hermanos).
                    El color <strong class="text-red-600">rojo</strong> significa "muy diferentes" (comparten pocos genes, como especies distintas).
                    <strong>Haz click en cualquier par</strong> para ver que genes comparten y en que se diferencian.
                </p>
            </div>

            <!-- Ranking de pares (solo los que comparten algo) -->
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h4 class="text-sm font-semibold text-primary mb-1">Ranking de Similitud</h4>
                <p class="text-xs text-secondary mb-3">Pares que comparten al menos algun gen, ordenados de mas a menos parecidos. <strong>Click para ver genes en detalle.</strong></p>
                <div class="space-y-1 max-h-[450px] overflow-y-auto" id="evo-ranking-list">
                    ${pares.length === 0 ? '<p class="text-sm text-secondary py-4 text-center">No hay pares con genes en comun</p>' : ''}
                    ${pares.map((p, idx) => {
                        const simPct = parseFloat(p.similitud);
                        const barColor = simPct > 80 ? '#10b981' : simPct > 50 ? '#f59e0b' : '#ef4444';
                        const relacion = simPct > 90 ? 'Casi identicos' : simPct > 70 ? 'Muy parecidos' : simPct > 50 ? 'Moderadamente parecidos' : simPct > 30 ? 'Bastante diferentes' : 'Muy diferentes';
                        const emoji = simPct > 90 ? '&#9733;' : simPct > 70 ? '&#9830;' : simPct > 50 ? '&#9679;' : '&#9651;';
                        return `
                            <div class="rounded-lg border border-slate-100 dark:border-slate-800 overflow-hidden">
                                <div class="flex items-center gap-3 p-2.5 hover:bg-slate-50 dark:hover:bg-slate-800 transition cursor-pointer"
                                     onclick="DashboardRenderer._expandParDetail(${idx}, '${p.a}', '${p.b}')">
                                    <span class="text-xs font-mono text-secondary w-6 text-right">${idx + 1}.</span>
                                    <div class="flex-1 min-w-0">
                                        <div class="flex items-center gap-1 text-xs">
                                            <span class="text-primary font-medium truncate">${p.a.replace(/_/g,' ').substring(0,28)}</span>
                                            <span class="text-secondary">vs</span>
                                            <span class="text-primary font-medium truncate">${p.b.replace(/_/g,' ').substring(0,28)}</span>
                                        </div>
                                        <div class="w-full bg-slate-100 dark:bg-slate-800 rounded-full h-2 mt-1.5">
                                            <div class="h-2 rounded-full transition-all" style="width:${simPct}%; background:${barColor};"></div>
                                        </div>
                                    </div>
                                    <div class="text-right shrink-0 flex items-center gap-2">
                                        <div>
                                            <span class="text-sm font-bold" style="color:${barColor}">${p.similitud}%</span>
                                            <p class="text-[10px] text-secondary">${relacion}</p>
                                        </div>
                                        <span class="text-secondary text-xs" id="par-arrow-${idx}">&#9654;</span>
                                    </div>
                                </div>
                                <div id="par-detail-${idx}" class="hidden border-t border-slate-100 dark:border-slate-800 bg-slate-50 dark:bg-slate-800/30 p-3">
                                    <p class="text-xs text-secondary">Calculando genes...</p>
                                </div>
                            </div>
                        `;
                    }).join('')}
                </div>
                ${paresCero.length > 0 ? `
                <div class="mt-3 p-3 bg-slate-50 dark:bg-slate-800/30 rounded-lg border border-slate-200">
                    <p class="text-xs text-secondary"><strong>${paresCero.length} pares</strong> no comparten ningun gen (0% similitud) y no se muestran en el ranking.</p>
                </div>` : ''}
            </div>

            <!-- Heatmap -->
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h4 class="text-sm font-semibold text-primary mb-1">Mapa de Calor (Heatmap)</h4>
                <p class="text-xs text-secondary mb-3">Cada celda muestra el % de genes que comparten dos bacterias. Verde = muchos genes en comun, Rojo = pocos. Pasa el mouse sobre cada celda.</p>
                <div class="overflow-x-auto">
                    <table class="text-xs border-collapse" id="evo-heatmap-table">
                        <thead>
                            <tr>
                                <th class="px-1 py-1 text-secondary text-left sticky left-0 bg-card z-10" style="min-width:30px;"></th>
                                ${matrizGenomas.map(g => `<th class="px-1 py-1 text-secondary text-center" style="writing-mode:vertical-lr; max-height:130px; font-size:10px; cursor:default;" title="${g.replace(/_/g,' ')}">${g.replace(/_/g,' ').substring(0, 22)}</th>`).join('')}
                            </tr>
                        </thead>
                        <tbody>
                            ${matrizGenomas.map((g, i) => `
                                <tr>
                                    <td class="px-2 py-1.5 text-secondary font-medium whitespace-nowrap sticky left-0 bg-card z-10 cursor-default" style="font-size:10px;" title="${g.replace(/_/g,' ')}">${g.replace(/_/g,' ').substring(0, 25)}</td>
                                    ${matrizGenomas.map((g2, j) => {
                                        if (i === j) return `<td class="px-1 py-1.5 text-center" style="background:rgba(16,185,129,0.15); min-width:42px; font-size:9px;">100%</td>`;
                                        const val = (matriz[i] && matriz[i][j] !== undefined) ? matriz[i][j] : 0;
                                        const sim = ((1 - val) * 100).toFixed(0);
                                        const color = this._distToColor(val);
                                        const relacion = sim > 90 ? 'Casi identicos' : sim > 70 ? 'Muy parecidos' : sim > 50 ? 'Moderadamente parecidos' : sim > 30 ? 'Bastante diferentes' : 'Muy diferentes';
                                        return `<td class="px-1 py-1.5 text-center font-mono cursor-default hover:ring-2 hover:ring-emerald-500 hover:z-10 relative transition"
                                                    style="background:${color}; min-width:42px; font-size:9px;"
                                                    title="${g.replace(/_/g,' ')} vs ${g2.replace(/_/g,' ')}: ${sim}% similares. ${relacion}.">${sim}%</td>`;
                                    }).join('')}
                                </tr>
                            `).join('')}
                        </tbody>
                    </table>
                </div>
            </div>

            <!-- Presencia/ausencia -->
            ${presencia.genes && presencia.genes.length > 0 ? `
            <div class="bg-card rounded-xl p-5 border border-slate-200">
                <h4 class="text-sm font-semibold text-primary mb-1">Que genes tiene cada bacteria?</h4>
                <p class="text-xs text-secondary mb-3">Cada fila es un gen y cada columna una bacteria. <span class="inline-block w-3 h-3 rounded" style="background:rgba(16,185,129,0.3)"></span> = lo tiene, vacio = no lo tiene. Mostrando los ${Math.min(100, presencia.genes.length)} genes mas frecuentes.</p>
                <div class="overflow-x-auto max-h-[400px] overflow-y-auto">
                    <table class="text-[9px] border-collapse">
                        <thead class="sticky top-0 z-10">
                            <tr class="bg-card">
                                <th class="px-1 py-1 text-secondary text-left sticky left-0 bg-card z-20">Gen (funcion)</th>
                                ${(presencia.genomas || []).map(g => `<th class="px-1 py-1 text-secondary text-center" style="writing-mode:vertical-lr; max-height:90px;" title="${g.replace(/_/g,' ')}">${g.replace(/_/g,' ').substring(0, 15)}</th>`).join('')}
                            </tr>
                        </thead>
                        <tbody>
                            ${(presencia.genes || []).slice(0, 100).map((gen, gi) => {
                                const count = (presencia.genomas || []).reduce((acc, _, gj) => acc + ((presencia.matriz?.[gi]?.[gj]) || 0), 0);
                                const total = (presencia.genomas || []).length;
                                return `
                                <tr class="border-t border-slate-50 dark:border-slate-900 hover:bg-slate-50 dark:hover:bg-slate-800" title="Gen: ${gen} | Presente en ${count} de ${total} bacterias">
                                    <td class="px-1 py-0.5 text-secondary whitespace-nowrap sticky left-0 bg-card z-10" style="max-width:180px; overflow:hidden; text-overflow:ellipsis;">${gen.substring(0, 30)}</td>
                                    ${(presencia.genomas || []).map((_, gj) => {
                                        const val = presencia.matriz?.[gi]?.[gj] || 0;
                                        return `<td class="px-1 py-0.5 text-center" style="background:${val ? 'rgba(16,185,129,0.25)' : 'transparent'}; min-width:20px;">${val ? '<span style="color:#10b981">&#9679;</span>' : ''}</td>`;
                                    }).join('')}
                                </tr>`;
                            }).join('')}
                        </tbody>
                    </table>
                </div>
                ${(presencia.genes || []).length > 100 ? `<p class="text-xs text-secondary mt-2">Mostrando 100 de ${this.fmt(presencia.genes.length)} genes</p>` : ''}
            </div>` : ''}
        `;
    },

    _expandParDetail(idx, a, b) {
        const el = document.getElementById(`par-detail-${idx}`);
        const arrow = document.getElementById(`par-arrow-${idx}`);
        if (!el) return;

        if (!el.classList.contains('hidden')) {
            el.classList.add('hidden');
            if (arrow) arrow.innerHTML = '&#9654;';
            return;
        }
        el.classList.remove('hidden');
        if (arrow) arrow.innerHTML = '&#9660;';

        // Skip recomputation if already done
        if (el.dataset.computed) return;
        el.dataset.computed = '1';

        const presencia = DashboardRenderer._evoPresencia || {};
        const genes = presencia.genes || [];
        const genomas = presencia.genomas || [];
        const matrizP = presencia.matriz || [];
        const idxA = genomas.indexOf(a);
        const idxB = genomas.indexOf(b);

        if (idxA < 0 || idxB < 0 || genes.length === 0) {
            el.innerHTML = '<p class="text-xs text-secondary p-2">No hay datos de presencia/ausencia disponibles para este par.</p>';
            return;
        }

        const compartidos = [];
        const soloA = [];
        const soloB = [];

        for (let i = 0; i < genes.length; i++) {
            const inA = matrizP[i]?.[idxA] || 0;
            const inB = matrizP[i]?.[idxB] || 0;
            if (inA && inB) compartidos.push(genes[i]);
            else if (inA) soloA.push(genes[i]);
            else if (inB) soloB.push(genes[i]);
        }

        const labelA = a.replace(/_/g, ' ').substring(0, 22);
        const labelB = b.replace(/_/g, ' ').substring(0, 22);
        const maxShow = 40;

        el.innerHTML = `
            <div class="grid grid-cols-1 md:grid-cols-3 gap-3">
                <div class="bg-emerald-50 dark:bg-emerald-900/15 rounded-lg p-3">
                    <p class="text-xs font-bold text-emerald-600 mb-2">Genes compartidos (${this.fmt(compartidos.length)})</p>
                    <p class="text-[10px] text-secondary mb-1">Estos genes estan en ambas bacterias:</p>
                    <div class="max-h-[150px] overflow-y-auto text-[10px] text-secondary leading-relaxed">
                        ${compartidos.slice(0, maxShow).map(g => `<span class="inline-block bg-emerald-100 dark:bg-emerald-900/30 rounded px-1 py-0.5 mr-1 mb-1">${g.substring(0, 35)}</span>`).join('')}
                        ${compartidos.length > maxShow ? `<p class="text-[10px] text-emerald-600 mt-1">...y ${this.fmt(compartidos.length - maxShow)} genes mas</p>` : ''}
                        ${compartidos.length === 0 ? '<p class="text-secondary italic">Ninguno</p>' : ''}
                    </div>
                </div>
                <div class="bg-amber-50 dark:bg-amber-900/15 rounded-lg p-3">
                    <p class="text-xs font-bold text-amber-600 mb-2">Solo en ${labelA} (${this.fmt(soloA.length)})</p>
                    <p class="text-[10px] text-secondary mb-1">Genes que la otra bacteria NO tiene:</p>
                    <div class="max-h-[150px] overflow-y-auto text-[10px] text-secondary leading-relaxed">
                        ${soloA.slice(0, maxShow).map(g => `<span class="inline-block bg-amber-100 dark:bg-amber-900/30 rounded px-1 py-0.5 mr-1 mb-1">${g.substring(0, 35)}</span>`).join('')}
                        ${soloA.length > maxShow ? `<p class="text-[10px] text-amber-600 mt-1">...y ${this.fmt(soloA.length - maxShow)} genes mas</p>` : ''}
                        ${soloA.length === 0 ? '<p class="text-secondary italic">Ninguno</p>' : ''}
                    </div>
                </div>
                <div class="bg-red-50 dark:bg-red-900/15 rounded-lg p-3">
                    <p class="text-xs font-bold text-red-600 mb-2">Solo en ${labelB} (${this.fmt(soloB.length)})</p>
                    <p class="text-[10px] text-secondary mb-1">Genes que la otra bacteria NO tiene:</p>
                    <div class="max-h-[150px] overflow-y-auto text-[10px] text-secondary leading-relaxed">
                        ${soloB.slice(0, maxShow).map(g => `<span class="inline-block bg-red-100 dark:bg-red-900/30 rounded px-1 py-0.5 mr-1 mb-1">${g.substring(0, 35)}</span>`).join('')}
                        ${soloB.length > maxShow ? `<p class="text-[10px] text-red-600 mt-1">...y ${this.fmt(soloB.length - maxShow)} genes mas</p>` : ''}
                        ${soloB.length === 0 ? '<p class="text-secondary italic">Ninguno</p>' : ''}
                    </div>
                </div>
            </div>
        `;
    },

    _matrizMinMax(matriz) {
        let min = Infinity, max = -Infinity;
        for (let i = 0; i < matriz.length; i++) {
            for (let j = i + 1; j < matriz.length; j++) {
                const v = matriz[i][j];
                if (v < min) min = v;
                if (v > max) max = v;
            }
        }
        return { min: min === Infinity ? 0 : min, max: max === -Infinity ? 0 : max };
    },

    _distToColor(val) {
        const clamped = Math.max(0, Math.min(1, val));
        const r = Math.round(239 * clamped + 16 * (1 - clamped));
        const g = Math.round(68 * clamped + 185 * (1 - clamped));
        const b = Math.round(68 * clamped + 129 * (1 - clamped));
        return `rgba(${r},${g},${b},0.3)`;
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
