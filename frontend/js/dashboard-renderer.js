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
                                    <p class="text-center text-slate-400">Cargando preview...</p>
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
                <p class="text-xs text-secondary mb-3">Vista lineal del cromosoma. Genes forward (+) arriba, reverse (-) abajo. Pase el cursor para ver detalles.</p>
                <div class="flex gap-2 mb-3">
                    <button onclick="DashboardRenderer._zoomGenomeMap(1.5)" class="px-3 py-1 bg-slate-100 dark:bg-slate-800 rounded text-xs hover:bg-slate-200 transition">Zoom +</button>
                    <button onclick="DashboardRenderer._zoomGenomeMap(0.67)" class="px-3 py-1 bg-slate-100 dark:bg-slate-800 rounded text-xs hover:bg-slate-200 transition">Zoom -</button>
                    <button onclick="DashboardRenderer._zoomGenomeMap(0)" class="px-3 py-1 bg-slate-100 dark:bg-slate-800 rounded text-xs hover:bg-slate-200 transition">Reset</button>
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
                <p class="text-xs text-secondary mb-3">Grupos de genes adyacentes en la misma hebra con distancia intergenica < 50 pb</p>
                <div class="overflow-x-auto max-h-96 overflow-y-auto">
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
                            ${(operones.operones || []).sort((a, b) => b.num_genes - a.num_genes).slice(0, 20).map(op => `
                                <tr class="hover:bg-slate-50 dark:hover:bg-slate-800">
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700">${op.numero}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 font-mono text-primary">${op.genes.slice(0, 6).join(', ')}${op.genes.length > 6 ? ' ...' : ''}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-center font-bold text-emerald-500">${op.num_genes}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-center">${op.hebra}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-right">${this.fmt(op.inicio)}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-right">${this.fmt(op.fin)}</td>
                                    <td class="px-2 py-2 border-b border-slate-100 dark:border-slate-700 text-right">${this.fmt(op.longitud_total)} pb</td>
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
        svg += `<text x="${margin}" y="25" fill="#10b981" font-size="11" font-weight="bold">Forward (+) - ${genes.filter(g => (g.hebra || g.strand) === '+' || (g.hebra || g.strand) === '1').length} genes</text>`;
        svg += `<text x="${margin}" y="${height - 15}" fill="#06b6d4" font-size="11" font-weight="bold">Reverse (-) - ${genes.filter(g => (g.hebra || g.strand) === '-' || (g.hebra || g.strand) === '-1').length} genes</text>`;

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
                    <p class="text-xs text-secondary mt-2">Se busca en ambas hebras (forward y reverse complement). Maximo 200 resultados.</p>
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
