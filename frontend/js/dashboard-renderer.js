/**
 * GenomeHub - Dashboard Renderer
 *
 * Renderiza dashboards interactivos desde datos JSON usando Chart.js.
 * Cada análisis tiene su propio dashboard con stats cards, gráficos y tablas.
 */

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
                </div>
            </div>

            <!-- Referencias bibliográficas -->
            ${Object.keys(refs).length > 0 ? `
            <div class="bg-card rounded-xl p-5 border border-slate-200 mb-6">
                <h3 class="text-sm font-semibold text-primary mb-3">Referencias Bibliográficas</h3>
                <div class="space-y-2 text-sm text-secondary">
                    ${Object.entries(refs).map(([k, v]) => `<p><span class="font-medium">${k}:</span> ${v}</p>`).join('')}
                </div>
            </div>` : ''}

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
                        labels: ['Forward (+)', 'Reverse (-)'],
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
                                    <td class="px-3 py-2 text-secondary">${info.aminoacido}</td>
                                    <td class="px-3 py-2 text-right text-primary">${this.fmt(info.conteo)}</td>
                                    <td class="px-3 py-2 text-right text-primary">${info.frecuencia_porcentaje}%</td>
                                    <td class="px-3 py-2 text-right text-secondary">${info.densidad_por_kb}</td>
                                </tr>
                            `).join('')}
                        </tbody>
                    </table>
                </div>
            </div>` : ''}

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
                        plugins: { legend: { display: false } },
                        scales: {
                            x: { ticks: { font: { size: 9, family: 'monospace' }, maxRotation: 90 } },
                            y: { beginAtZero: true }
                        }
                    }
                });
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
                <h3 class="text-sm font-semibold text-primary mb-4">Top 20 Regiones Intergénicas Grandes (posibles islas de patogenicidad)</h3>
                <div class="overflow-x-auto max-h-96 overflow-y-auto">
                    <table class="w-full text-sm">
                        <thead class="bg-slate-50 dark:bg-slate-800 sticky top-0">
                            <tr>
                                <th class="px-3 py-2 text-left text-secondary font-medium">#</th>
                                <th class="px-3 py-2 text-left text-secondary font-medium">Gen 1</th>
                                <th class="px-3 py-2 text-left text-secondary font-medium">Gen 2</th>
                                <th class="px-3 py-2 text-right text-secondary font-medium">Distancia (pb)</th>
                                <th class="px-3 py-2 text-center text-secondary font-medium">Hebras</th>
                            </tr>
                        </thead>
                        <tbody>
                            ${top20.map((r, i) => `
                                <tr class="border-t border-slate-100 dark:border-slate-800">
                                    <td class="px-3 py-2 text-secondary">${i + 1}</td>
                                    <td class="px-3 py-2 font-mono text-primary">${r.gen1}</td>
                                    <td class="px-3 py-2 font-mono text-primary">${r.gen2}</td>
                                    <td class="px-3 py-2 text-right font-bold text-emerald-500">${this.fmt(r.distancia_pb)}</td>
                                    <td class="px-3 py-2 text-center text-secondary">${r.hebra_gen1} → ${r.hebra_gen2}</td>
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
                            label: 'Total pares',
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

        const ecoli = metricas.ecoli || {};
        const salmonella = metricas.salmonella || {};

        container.innerHTML = `
            <!-- Stats lado a lado -->
            <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
                <div class="bg-card rounded-xl p-5 border border-emerald-500/30">
                    <h3 class="text-lg font-bold text-emerald-500 mb-3">Genoma 1</h3>
                    <div class="space-y-2 text-sm">
                        <p><span class="text-secondary">Longitud:</span> <span class="font-bold text-primary">${this.fmt(ecoli.longitud_genoma_pb)} pb</span></p>
                        <p><span class="text-secondary">Genes CDS:</span> <span class="font-bold text-primary">${this.fmt(ecoli.total_genes_cds)}</span></p>
                        <p><span class="text-secondary">GC:</span> <span class="font-bold text-primary">${ecoli.contenido_gc_porcentaje}%</span></p>
                        <p><span class="text-secondary">Densidad:</span> <span class="font-bold text-primary">${ecoli.densidad_genica_porcentaje}%</span></p>
                    </div>
                </div>
                <div class="bg-card rounded-xl p-5 border border-cyan-500/30">
                    <h3 class="text-lg font-bold text-cyan-500 mb-3">Genoma 2</h3>
                    <div class="space-y-2 text-sm">
                        <p><span class="text-secondary">Longitud:</span> <span class="font-bold text-primary">${this.fmt(salmonella.longitud_genoma_pb)} pb</span></p>
                        <p><span class="text-secondary">Genes CDS:</span> <span class="font-bold text-primary">${this.fmt(salmonella.total_genes_cds)}</span></p>
                        <p><span class="text-secondary">GC:</span> <span class="font-bold text-primary">${salmonella.contenido_gc_porcentaje}%</span></p>
                        <p><span class="text-secondary">Densidad:</span> <span class="font-bold text-primary">${salmonella.densidad_genica_porcentaje}%</span></p>
                    </div>
                </div>
            </div>

            <!-- Gráficos comparativos -->
            <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Comparación de Métricas</h3>
                    <canvas id="chart-compare-metrics" height="250"></canvas>
                </div>
                <div class="bg-card rounded-xl p-5 border border-slate-200">
                    <h3 class="text-sm font-semibold text-primary mb-4">Genes de Virulencia</h3>
                    <canvas id="chart-compare-virulence" height="250"></canvas>
                </div>
            </div>
        `;

        setTimeout(() => {
            // Métricas comparativas
            this.createChart('chart-compare-metrics', {
                type: 'bar',
                data: {
                    labels: ['Longitud (Mb)', 'Genes CDS', 'GC %', 'Densidad %'],
                    datasets: [
                        {
                            label: 'Genoma 1',
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
                            label: 'Genoma 2',
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

            // Virulencia
            if (virulencia.ecoli && virulencia.salmonella) {
                this.createChart('chart-compare-virulence', {
                    type: 'bar',
                    data: {
                        labels: ['Genoma 1', 'Genoma 2'],
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
        }, 100);
    },

    // =========================================================================
    // VISTA DE ARCHIVOS (lista simple con descargar/eliminar)
    // =========================================================================

    renderArchivosView(results, genome, container) {
        if (!results || results.length === 0) {
            container.innerHTML = `
                <div class="text-center py-12 text-secondary">
                    <div class="text-6xl mb-3">📄</div>
                    <p>No hay archivos de resultados para este genoma</p>
                </div>
            `;
            return;
        }

        container.innerHTML = `
            <div class="space-y-2">
                ${results.map(file => `
                    <div class="flex items-center justify-between bg-card rounded-lg px-5 py-3 border border-slate-200 hover:border-emerald-500/30 transition">
                        <div class="flex items-center gap-3 flex-1 min-w-0">
                            <span class="text-xl">${file.extension === 'json' ? '📋' : file.extension === 'png' ? '📊' : '📄'}</span>
                            <div class="min-w-0">
                                <p class="text-sm font-medium text-primary truncate">${file.filename}</p>
                                <p class="text-xs text-secondary">${file.size_kb} KB</p>
                            </div>
                        </div>
                        <div class="flex items-center gap-2 ml-4">
                            <a href="${file.path}" download class="px-3 py-1.5 bg-emerald-500 hover:bg-emerald-600 text-white text-xs rounded-lg transition">
                                Descargar
                            </a>
                            <button onclick="deleteResult('${file.filename}', 'tablas')" class="px-3 py-1.5 bg-red-500/10 hover:bg-red-500/20 text-red-500 text-xs rounded-lg transition">
                                Eliminar
                            </button>
                        </div>
                    </div>
                `).join('')}
            </div>
        `;
    },

    // =========================================================================
    // UTILIDADES
    // =========================================================================

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
