/**
 * GenomeHub - Analysis Runner + Results Dashboard
 *
 * Ejecuta scripts Python de análisis via servidor.py
 * Resultados como dashboards interactivos con Chart.js
 */

// =============================================================================
// ANÁLISIS RUNNER
// =============================================================================

const AnalysisRunner = {

    async loadGenomeSelector() {
        const select = document.getElementById('analysis-genome');
        if (!select) return;

        const hasServer = await NCBISearch.checkServer();
        if (!hasServer) {
            select.innerHTML = '<option value="">Servidor no disponible</option>';
            return;
        }

        try {
            const resp = await fetch('/api/genomes');
            const data = await resp.json();

            if (data.success && data.genomes && data.genomes.length > 0) {
                select.innerHTML = '<option value="">-- Seleccionar genoma --</option>' +
                    data.genomes.map(g =>
                        `<option value="${g.basename}">${g.organism || g.basename} (${g.accession_id || 'N/A'})</option>`
                    ).join('');

                if (AppState.analysisGenome) {
                    select.value = AppState.analysisGenome;
                    AppState.analysisGenome = null;
                }
            } else {
                select.innerHTML = '<option value="">No hay genomas descargados</option>';
            }
        } catch {
            select.innerHTML = '<option value="">Error al cargar genomas</option>';
        }
    },

    async runSelectedAnalyses() {
        const genomeSelect = document.getElementById('analysis-genome');
        const genome = genomeSelect.value;
        if (!genome) {
            showNotification('Selecciona un genoma', 'warning');
            return;
        }

        const checked = document.querySelectorAll('.analysis-check:checked');
        if (checked.length === 0) {
            showNotification('Selecciona al menos un análisis', 'warning');
            return;
        }

        const scripts = Array.from(checked).map(cb => cb.getAttribute('data-script'));
        const btn = document.getElementById('run-analysis-btn');
        const progress = document.getElementById('analysis-progress');
        const progressText = document.getElementById('analysis-progress-text');
        const progressDetail = document.getElementById('analysis-progress-detail');
        const consolePanel = document.getElementById('analysis-console');
        const consoleOutput = document.getElementById('console-output');

        btn.disabled = true;
        btn.textContent = 'Analizando...';
        progress.classList.remove('hidden');
        consolePanel.classList.remove('hidden');
        consoleOutput.innerHTML = '';

        let completed = 0;
        let errors = 0;

        for (const script of scripts) {
            completed++;
            progressText.textContent = `Ejecutando ${completed}/${scripts.length}...`;
            progressDetail.textContent = script;
            consoleOutput.innerHTML += `\n[${new Date().toLocaleTimeString()}] ${script}...\n`;

            try {
                const payload = { script, genome_basename: genome };

                const response = await fetch('/api/run_analysis', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify(payload)
                });

                const data = await response.json();

                if (!response.ok) {
                    throw new Error(data.error || 'Error en el servidor');
                }

                if (!data.success && data.error) {
                    throw new Error(data.error);
                }

                // Escapar HTML para evitar que <tags> se eliminen
                const safeOutput = (data.output || '').replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
                consoleOutput.innerHTML += safeOutput + '\n';

                if (data.return_code === 0) {
                    consoleOutput.innerHTML += `[OK] Completado en ${data.execution_time}s\n`;
                } else {
                    consoleOutput.innerHTML += `[ERROR] Código: ${data.return_code}\n`;
                    errors++;
                }
            } catch (error) {
                consoleOutput.innerHTML += `[ERROR] ${error.message}\n`;
                errors++;
            }

            consoleOutput.scrollTop = consoleOutput.scrollHeight;
        }

        progress.classList.add('hidden');
        btn.disabled = false;
        btn.textContent = 'Analizar';

        DashboardRenderer.clearCache();
        AppState.resultsGenome = genome;

        if (errors === 0) {
            showNotification('Análisis completado', 'success');
        } else if (errors < scripts.length) {
            showNotification(`${scripts.length - errors} completados, ${errors} con errores`, 'warning');
        } else {
            showNotification(`${errors} análisis con errores`, 'error');
        }

        // Ir a resultados si al menos un analisis fue exitoso
        if (errors < scripts.length) {
            setTimeout(() => showSection('results'), 2000);
        }
    }
};


// =============================================================================
// RESULTADOS - DASHBOARDS
// =============================================================================

let currentResultsGenome = '';
let currentDashboardTab = 'genes';

async function loadResultsGenomeSelector() {
    const select = document.getElementById('results-genome');
    if (!select) return;

    const hasServer = await NCBISearch.checkServer();
    if (!hasServer) {
        select.innerHTML = '<option value="">Servidor no disponible</option>';
        return;
    }

    try {
        const resp = await fetch('/api/results_genomes');
        const data = await resp.json();

        if (data.success && data.genomes && data.genomes.length > 0) {
            select.innerHTML = data.genomes.map(g =>
                `<option value="${g.basename}">${g.label} (${g.tablas} tablas, ${g.figuras} gráficos)</option>`
            ).join('');

            if (AppState.resultsGenome) {
                select.value = AppState.resultsGenome;
                AppState.resultsGenome = null;
            }

            currentResultsGenome = select.value;
            loadDashboard('genes');
        } else {
            select.innerHTML = '<option value="">No hay resultados</option>';
            document.getElementById('dashboard-container').innerHTML = `
                <div class="text-center py-12 text-secondary">
                    <div class="text-6xl mb-3">📊</div>
                    <p>No hay resultados aún</p>
                    <p class="text-sm mt-2">Ejecuta un análisis primero</p>
                </div>
            `;
        }
    } catch {
        select.innerHTML = '<option value="">Error al cargar</option>';
    }
}

async function loadDashboard(analysisType) {
    const container = document.getElementById('dashboard-container');
    const genome = document.getElementById('results-genome')?.value || '';
    currentResultsGenome = genome;
    currentDashboardTab = analysisType;

    if (!genome) {
        container.innerHTML = `
            <div class="text-center py-12 text-secondary">
                <div class="text-6xl mb-3">📊</div>
                <p>Selecciona un genoma para ver sus dashboards</p>
            </div>
        `;
        return;
    }

    // Actualizar tabs activas
    document.querySelectorAll('.dashboard-tab').forEach(t => {
        t.classList.remove('active', 'bg-emerald-500', 'text-white');
        t.classList.add('text-secondary');
    });
    const activeTab = document.querySelector(`.dashboard-tab[data-analysis="${analysisType}"]`);
    if (activeTab) {
        activeTab.classList.add('active', 'bg-emerald-500', 'text-white');
        activeTab.classList.remove('text-secondary');
    }

    DashboardRenderer.destroyCharts();

    // Loading
    container.innerHTML = `
        <div class="flex items-center justify-center py-12">
            <div class="w-8 h-8 border-4 border-emerald-500 border-t-transparent rounded-full animate-spin mr-4"></div>
            <span class="text-secondary">Cargando dashboard...</span>
        </div>
    `;

    if (analysisType === 'archivos') {
        // Vista de archivos - cargar tablas y figuras en paralelo
        try {
            const [respTablas, respFiguras] = await Promise.all([
                fetch(`/api/results?type=tablas&genome=${genome}`),
                fetch(`/api/results?type=figuras&genome=${genome}`)
            ]);
            const dataTablas = await respTablas.json();
            const dataFiguras = await respFiguras.json();
            DashboardRenderer.renderArchivosView(
                dataTablas.results || [],
                dataFiguras.results || [],
                genome,
                container
            );
        } catch {
            container.innerHTML = '<p class="text-red-500 text-center py-8">Error al cargar archivos</p>';
        }
        return;
    }

    // Tab especial: buscar secuencia (no necesita JSON)
    if (analysisType === 'buscar') {
        DashboardRenderer.renderBusquedaSecuencia(genome, container);
        return;
    }

    // Mapear tab a archivo JSON
    const fileMap = {
        genes: `analisis_genes_${genome}.json`,
        codones: `analisis_codones_${genome}.json`,
        distancias: `analisis_distancias_${genome}.json`,
        estructura: `analisis_estructura_${genome}.json`
    };

    const filename = fileMap[analysisType];
    if (!filename) return;

    const data = await DashboardRenderer.fetchResultData(genome, filename);

    if (!data) {
        container.innerHTML = `
            <div class="text-center py-12 text-secondary">
                <div class="text-6xl mb-3">📭</div>
                <p>No hay datos de ${analysisType} para este genoma</p>
                <p class="text-sm mt-2">Ejecuta el análisis correspondiente primero</p>
                <button onclick="explicarErrorAnalisis('${analysisType}', '${genome}')" class="mt-4 px-4 py-2 bg-emerald-500 hover:bg-emerald-600 text-white text-sm rounded-lg transition">
                    Explicar con IA
                </button>
            </div>
        `;
        return;
    }

    // Renderizar dashboard según tipo
    if (analysisType === 'genes') {
        DashboardRenderer.renderGenesDashboard(data, container);
    } else if (analysisType === 'codones') {
        DashboardRenderer.renderCodonesDashboard(data, container);
    } else if (analysisType === 'distancias') {
        DashboardRenderer.renderDistanciasDashboard(data, container);
    } else if (analysisType === 'estructura') {
        DashboardRenderer.renderEstructuraDashboard(data, genome, container);
    }
}

async function deleteResult(filename, type) {
    if (!confirm(`¿Eliminar ${filename}?`)) return;

    try {
        const resp = await fetch('/api/delete_result', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ filename, type, genome: currentResultsGenome })
        });
        const data = await resp.json();

        if (data.success) {
            showNotification('Eliminado', 'success');
            DashboardRenderer.clearCache();
            loadDashboard(currentDashboardTab);
        } else {
            showNotification(data.error || 'Error', 'error');
        }
    } catch {
        showNotification('Error al eliminar', 'error');
    }
}

async function deleteAllResults() {
    const genome = currentResultsGenome;
    if (!genome) {
        showNotification('Selecciona un genoma', 'warning');
        return;
    }

    if (!confirm('¿Eliminar todos los resultados de este genoma?')) return;

    try {
        const resp = await fetch('/api/delete_results_all', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ genome })
        });
        const data = await resp.json();

        if (data.success) {
            showNotification(`${data.deleted} archivos eliminados`, 'success');
            DashboardRenderer.clearCache();
            loadResultsGenomeSelector();
        } else {
            showNotification('Error al eliminar', 'error');
        }
    } catch {
        showNotification('Error al eliminar', 'error');
    }
}


// =============================================================================
// COMPARAR GENOMAS
// =============================================================================

async function loadCompareSelectors() {
    const select1 = document.getElementById('compare-genome-1');
    const select2 = document.getElementById('compare-genome-2');
    if (!select1 || !select2) return;

    const hasServer = await NCBISearch.checkServer();
    if (!hasServer) {
        select1.innerHTML = '<option value="">Servidor no disponible</option>';
        select2.innerHTML = '<option value="">Servidor no disponible</option>';
        return;
    }

    try {
        const resp = await fetch('/api/genomes');
        const data = await resp.json();

        if (data.success && data.genomes && data.genomes.length >= 2) {
            const options = data.genomes.map(g =>
                `<option value="${g.basename}">${g.organism || g.basename}</option>`
            ).join('');

            select1.innerHTML = '<option value="">-- Seleccionar --</option>' + options;
            select2.innerHTML = '<option value="">-- Seleccionar --</option>' + options;

            // Pre-seleccionar los primeros 2
            if (data.genomes.length >= 2) {
                select1.value = data.genomes[0].basename;
                select2.value = data.genomes[1].basename;
            }
        } else {
            select1.innerHTML = '<option value="">Se necesitan al menos 2 genomas</option>';
            select2.innerHTML = '<option value="">Se necesitan al menos 2 genomas</option>';
        }
    } catch {
        select1.innerHTML = '<option value="">Error</option>';
        select2.innerHTML = '<option value="">Error</option>';
    }
}

async function runComparison() {
    const genome1 = document.getElementById('compare-genome-1')?.value;
    const genome2 = document.getElementById('compare-genome-2')?.value;

    if (!genome1 || !genome2) {
        showNotification('Selecciona dos genomas', 'warning');
        return;
    }

    if (genome1 === genome2) {
        showNotification('Selecciona genomas diferentes', 'warning');
        return;
    }

    const btn = document.getElementById('run-compare-btn');
    const progress = document.getElementById('compare-progress');
    const dashboard = document.getElementById('compare-dashboard');

    btn.disabled = true;
    btn.textContent = 'Comparando...';
    progress.classList.remove('hidden');
    dashboard.innerHTML = '';

    try {
        // Ejecutar comparar_genomas con ambos genomas
        const response = await fetch('/api/run_analysis', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                script: 'comparar_genomas',
                genome_basename: genome1,
                genome_basename_2: genome2
            })
        });

        const result = await response.json();

        if (result.return_code === 0) {
            // Cargar datos de comparación con nombre dinámico
            const compFile = `comparacion_${genome1}_vs_${genome2}.json`;
            const data = await DashboardRenderer.fetchResultData(genome1, compFile);

            if (data) {
                DashboardRenderer.destroyCharts();
                DashboardRenderer.renderComparacionDashboard(data, dashboard);
                showNotification('Comparación completada', 'success');
            } else {
                dashboard.innerHTML = `
                    <div class="text-center py-12 text-secondary">
                        <div class="text-6xl mb-3">📭</div>
                        <p>No se pudieron cargar los resultados de comparación</p>
                    </div>
                `;
            }
        } else {
            showNotification('Error en la comparación', 'error');
            const safeCompOutput = (result.output || 'Error desconocido').replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
            dashboard.innerHTML = `
                <div class="bg-card rounded-xl p-6 border border-red-500/30">
                    <h3 class="text-red-500 font-semibold mb-2">Error</h3>
                    <pre class="text-sm text-secondary whitespace-pre-wrap">${safeCompOutput}</pre>
                </div>
            `;
        }
    } catch (error) {
        showNotification('Error al comparar', 'error');
        dashboard.innerHTML = `<p class="text-red-500 text-center py-8">${error.message}</p>`;
    } finally {
        btn.disabled = false;
        btn.textContent = 'Comparar';
        progress.classList.add('hidden');
    }
}


// =============================================================================
// GENERAR INFORME PDF
// =============================================================================

async function generarInforme() {
    const genome = currentResultsGenome;
    if (!genome) {
        showNotification('Selecciona un genoma primero en la seccion Resultados', 'warning');
        return;
    }

    const btn = document.getElementById('generate-report-btn');
    if (!btn) return;

    btn.disabled = true;
    btn.textContent = 'Generando informe...';
    btn.classList.add('opacity-70', 'cursor-wait');

    // Mostrar barra de progreso
    const progressContainer = document.getElementById('report-progress-container');
    const progressBar = document.getElementById('report-progress-bar');
    const progressPercent = document.getElementById('report-progress-percent');
    const progressText = document.getElementById('report-progress-text');
    
    progressContainer.classList.remove('hidden');
    progressBar.style.width = '0%';
    progressPercent.textContent = '0%';
    progressText.textContent = 'Recopilando datos de análisis...';

    // Simular progreso en incrementos
    let progress = 0;
    const progressInterval = setInterval(() => {
        if (progress < 30) {
            progress += Math.random() * 5;
        } else if (progress < 60) {
            progress += Math.random() * 3;
        } else if (progress < 90) {
            progress += Math.random() * 1.5;
        }
        progress = Math.min(progress, 90);
        updateProgressBar(progress);
    }, 500);

    // Actualizar barra de progreso
    function updateProgressBar(percent) {
        percent = Math.min(percent, 100);
        progressBar.style.width = percent + '%';
        progressPercent.textContent = Math.round(percent) + '%';
        
        if (percent < 30) {
            progressText.textContent = 'Recopilando datos de análisis...';
        } else if (percent < 60) {
            progressText.textContent = 'Generando interpretaciones con IA...';
        } else if (percent < 90) {
            progressText.textContent = 'Construyendo documento PDF...';
        } else {
            progressText.textContent = 'Finalizando generación...';
        }
    }

    showNotification('Generando informe PDF con IA... esto puede tardar unos minutos', 'info');

    try {
        const response = await fetch(`/api/generar_informe?genome=${genome}`);

        clearInterval(progressInterval);
        updateProgressBar(95);

        if (response.ok && response.headers.get('Content-Type')?.includes('application/pdf')) {
            // El servidor devolvio el PDF directamente
            const blob = await response.blob();
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = `informe_${genome}.pdf`;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(url);
            
            updateProgressBar(100);
            setTimeout(() => {
                showNotification('Informe PDF descargado exitosamente', 'success');
                progressContainer.classList.add('hidden');
            }, 500);
        } else {
            const data = await response.json();
            const errorMsg = data.error || 'Error desconocido';
            clearInterval(progressInterval);
            progressContainer.classList.add('hidden');
            showNotification(`Error: ${errorMsg}`, 'error');
            if (data.output) {
                console.error('Output del script:', data.output);
            }
        }
    } catch (error) {
        clearInterval(progressInterval);
        progressContainer.classList.add('hidden');
        showNotification(`Error al generar informe: ${error.message}`, 'error');
    } finally {
        btn.disabled = false;
        btn.textContent = 'Generar Informe PDF (IEEE)';
        btn.classList.remove('opacity-70', 'cursor-wait');
    }
}


// =============================================================================
// EXPLICADOR DE ERRORES CON IA
// =============================================================================

async function explicarErrorAnalisis(analysisType, genome) {
    const container = document.getElementById('dashboard-results');
    if (!container) return;

    container.innerHTML = `
        <div class="text-center py-12">
            <p class="text-secondary mb-4">Consultando con IA por qué no hay datos de ${analysisType}...</p>
            <div class="inline-block animate-spin rounded-full h-8 w-8 border-b-2 border-emerald-500"></div>
        </div>
    `;

    try {
        const prompt = `El usuario está intentando ver el análisis de ${analysisType} del genoma "${genome}" en GenomeHub (una aplicación de análisis genómico). 

El análisis no generó datos. Explica de manera clara y educativa:
1. Por qué podría no haber datos para este análisis en este genoma
2. Posibles razones técnicas (errores durante el análisis, formatos de archivo, etc)
3. Qué hacer para resolver el problema
4. Si es normal o esperado en algunos casos

Sé conciso pero informativo. Usa ejemplos relevantes a bioinformática.`;

        const response = await fetch('/api/chat', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ message: prompt, context: 'error_explanation' })
        });

        if (!response.ok) throw new Error('Error en la API');

        const data = await response.json();
        const explanation = data.response || 'No se pudo obtener una explicación en este momento';

        container.innerHTML = `
            <div class="bg-card rounded-xl p-6 border border-amber-200 dark:border-amber-800">
                <div class="flex items-start gap-3 mb-4">
                    <div class="text-3xl">🤖</div>
                    <div>
                        <h3 class="text-lg font-bold text-primary">Explicación de IA</h3>
                        <p class="text-xs text-secondary">Análisis: ${analysisType} | Genoma: ${genome}</p>
                    </div>
                </div>
                <div class="prose prose-sm dark:prose-invert max-w-none text-secondary leading-relaxed">
                    ${explanation.split('\n').map(line => {
                        if (line.trim()) {
                            if (line.match(/^#+\s/)) {
                                // Headers
                                return `<h4 class="font-semibold text-primary mt-3 mb-2">${line.replace(/^#+\s/, '')}</h4>`;
                            } else if (line.match(/^[-•*]/)) {
                                // Bullet points
                                return `<li class="ml-4">${line.replace(/^[-•*]\s/, '')}</li>`;
                            } else {
                                // Paragraphs
                                return `<p class="mb-2">${line}</p>`;
                            }
                        }
                        return '';
                    }).join('')}
                </div>
                <button onclick="loadDashboard('${analysisType}')" class="mt-4 px-4 py-2 bg-emerald-500 hover:bg-emerald-600 text-white text-sm rounded-lg transition">
                    Reintentar análisis
                </button>
            </div>
        `;
    } catch (error) {
        container.innerHTML = `
            <div class="text-center py-12 text-red-500">
                <p class="mb-2">❌ Error al obtener explicación</p>
                <p class="text-sm text-secondary">${error.message}</p>
                <button onclick="loadDashboard('${analysisType}')" class="mt-4 px-4 py-2 bg-slate-500 text-white text-sm rounded-lg">
                    Volver
                </button>
            </div>
        `;
    }
}


// =============================================================================
// EVENT LISTENERS
// =============================================================================

document.addEventListener('DOMContentLoaded', () => {
    // Botón analizar
    const runBtn = document.getElementById('run-analysis-btn');
    if (runBtn) {
        runBtn.addEventListener('click', () => AnalysisRunner.runSelectedAnalyses());
    }

    // Limpiar consola
    const clearConsoleBtn = document.getElementById('clear-console');
    if (clearConsoleBtn) {
        clearConsoleBtn.addEventListener('click', () => {
            document.getElementById('console-output').innerHTML = '';
        });
    }

    // Selector de genoma en resultados
    const resultsGenome = document.getElementById('results-genome');
    if (resultsGenome) {
        resultsGenome.addEventListener('change', () => {
            currentResultsGenome = resultsGenome.value;
            DashboardRenderer.clearCache();
            if (typeof ChatIA !== 'undefined') ChatIA.destroy();
            loadDashboard('genes');
            // Reset tabs
            document.querySelectorAll('.dashboard-tab').forEach(t => {
                t.classList.remove('active', 'bg-emerald-500', 'text-white');
                t.classList.add('text-secondary');
            });
            const tabGenes = document.getElementById('tab-genes');
            if (tabGenes) {
                tabGenes.classList.add('active', 'bg-emerald-500', 'text-white');
                tabGenes.classList.remove('text-secondary');
            }
        });
    }

    // Dashboard tabs
    document.querySelectorAll('.dashboard-tab').forEach(tab => {
        tab.addEventListener('click', () => {
            const analysis = tab.getAttribute('data-analysis');
            loadDashboard(analysis);
        });
    });

    // Eliminar todos los resultados
    const deleteAllBtn = document.getElementById('delete-all-results-btn');
    if (deleteAllBtn) {
        deleteAllBtn.addEventListener('click', deleteAllResults);
    }

    // Botón comparar
    const compareBtn = document.getElementById('run-compare-btn');
    if (compareBtn) {
        compareBtn.addEventListener('click', runComparison);
    }
});
