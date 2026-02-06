/**
 * GenomeHub - Analysis Runner
 *
 * Ejecuta scripts Python de análisis via servidor.py
 * Selector de genoma + checkboxes de análisis + ejecución batch
 */

const AnalysisRunner = {

    /**
     * Cargar selector de genomas desde la biblioteca
     */
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

                // Pre-seleccionar si viene de Biblioteca
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

    /**
     * Ejecutar todos los análisis seleccionados
     */
    async runSelectedAnalyses() {
        const genome = document.getElementById('analysis-genome').value;
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
                const payload = { script };

                // Mapear organismo si el script lo requiere
                if (script === 'analisis_genes') {
                    // Determinar organismo basándose en el basename
                    if (genome.includes('ecoli') || genome.includes('escherichia')) {
                        payload.organism = 'ecoli_k12';
                    } else if (genome.includes('salmonella')) {
                        payload.organism = 'salmonella_lt2';
                    } else {
                        payload.organism = 'ecoli_k12';
                    }
                }

                const response = await fetch('/api/run_analysis', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify(payload)
                });

                const data = await response.json();

                if (!response.ok) {
                    throw new Error(data.error || 'Error en el servidor');
                }

                consoleOutput.innerHTML += data.output + '\n';

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

            // Scroll al final
            consoleOutput.scrollTop = consoleOutput.scrollHeight;
        }

        progress.classList.add('hidden');
        btn.disabled = false;
        btn.textContent = 'Analizar';

        if (errors === 0) {
            showNotification('Análisis completado', 'success');
            setTimeout(() => showSection('results'), 2000);
        } else {
            showNotification(`${errors} análisis con errores`, 'warning');
        }
    }
};

// =============================================================================
// RESULTADOS
// =============================================================================

async function loadResults(type = 'tablas') {
    const tablasContainer = document.getElementById('results-tablas');
    const figurasContainer = document.getElementById('results-figuras');
    const tablasGrid = document.getElementById('tablas-grid');
    const figurasGrid = document.getElementById('figuras-grid');

    const hasServer = await NCBISearch.checkServer();
    if (!hasServer) {
        const grid = type === 'tablas' ? tablasGrid : figurasGrid;
        if (type === 'tablas') {
            tablasContainer.classList.remove('hidden');
            figurasContainer.classList.add('hidden');
        } else {
            tablasContainer.classList.add('hidden');
            figurasContainer.classList.remove('hidden');
        }
        grid.innerHTML = `
            <div class="text-center py-12 text-secondary">
                <div class="text-6xl mb-3">🖥️</div>
                <p class="font-semibold">Servidor no disponible</p>
            </div>
        `;
        return;
    }

    try {
        const response = await fetch(`/api/results?type=${type}`);
        const data = await response.json();

        if (type === 'tablas') {
            tablasContainer.classList.remove('hidden');
            figurasContainer.classList.add('hidden');

            if (data.results.length === 0) {
                tablasGrid.innerHTML = `
                    <div class="text-center py-12 text-secondary">
                        <div class="text-6xl mb-3">📄</div>
                        <p>No hay resultados aún</p>
                        <p class="text-sm mt-2">Ejecuta un análisis primero</p>
                    </div>
                `;
            } else {
                tablasGrid.innerHTML = data.results.map(file => `
                    <div class="flex items-center justify-between bg-card rounded-lg px-5 py-3 border border-slate-200 hover:border-emerald-500/30 transition">
                        <div class="flex items-center gap-3 flex-1 min-w-0">
                            <span class="text-xl">${file.extension === 'json' ? '📋' : '📊'}</span>
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
                                🗑️
                            </button>
                        </div>
                    </div>
                `).join('');
            }
        } else {
            tablasContainer.classList.add('hidden');
            figurasContainer.classList.remove('hidden');

            if (data.results.length === 0) {
                figurasGrid.innerHTML = `
                    <div class="text-center py-12 text-secondary">
                        <div class="text-6xl mb-3">📈</div>
                        <p>No hay gráficos aún</p>
                        <p class="text-sm mt-2">Ejecuta un análisis primero</p>
                    </div>
                `;
            } else {
                figurasGrid.innerHTML = data.results.map(file => `
                    <div class="flex items-center justify-between bg-card rounded-lg px-5 py-3 border border-slate-200 hover:border-emerald-500/30 transition">
                        <div class="flex items-center gap-3 flex-1 min-w-0">
                            <img src="${file.path}" alt="${file.filename}" class="w-12 h-12 object-contain rounded cursor-pointer" onclick="viewImage('${file.path}', '${file.filename}')">
                            <div class="min-w-0">
                                <p class="text-sm font-medium text-primary truncate">${file.filename}</p>
                                <p class="text-xs text-secondary">${file.size_kb} KB</p>
                            </div>
                        </div>
                        <div class="flex items-center gap-2 ml-4">
                            <a href="${file.path}" download class="px-3 py-1.5 bg-emerald-500 hover:bg-emerald-600 text-white text-xs rounded-lg transition">
                                Descargar
                            </a>
                            <button onclick="deleteResult('${file.filename}', 'figuras')" class="px-3 py-1.5 bg-red-500/10 hover:bg-red-500/20 text-red-500 text-xs rounded-lg transition">
                                🗑️
                            </button>
                        </div>
                    </div>
                `).join('');
            }
        }
    } catch (error) {
        showNotification('Error al cargar resultados', 'error');
    }
}

async function deleteResult(filename, type) {
    if (!confirm(`¿Eliminar ${filename}?`)) return;

    try {
        const resp = await fetch('/api/delete_result', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ filename, type })
        });
        const data = await resp.json();

        if (data.success) {
            showNotification('Eliminado', 'success');
            loadResults(type);
        } else {
            showNotification(data.error || 'Error', 'error');
        }
    } catch {
        showNotification('Error al eliminar', 'error');
    }
}

async function deleteAllResults() {
    // Determinar qué tab está activo
    const tablasVisible = !document.getElementById('results-tablas').classList.contains('hidden');
    const type = tablasVisible ? 'tablas' : 'figuras';

    if (!confirm(`¿Eliminar todos los resultados de ${type}?`)) return;

    try {
        const resp = await fetch('/api/delete_results_all', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ type })
        });
        const data = await resp.json();

        if (data.success) {
            showNotification(`${data.deleted} archivos eliminados`, 'success');
            loadResults(type);
        } else {
            showNotification('Error al eliminar', 'error');
        }
    } catch {
        showNotification('Error al eliminar', 'error');
    }
}

function viewImage(path, filename) {
    const modal = document.createElement('div');
    modal.className = 'fixed inset-0 z-50 flex items-center justify-center p-4 bg-black/90 backdrop-blur-sm';
    modal.onclick = (e) => { if (e.target === modal) modal.remove(); };

    modal.innerHTML = `
        <div class="max-w-6xl w-full">
            <div class="flex items-center justify-between mb-4">
                <h3 class="text-white font-medium">${filename}</h3>
                <button onclick="this.closest('.fixed').remove()" class="text-white hover:text-gray-300 text-2xl">✕</button>
            </div>
            <img src="${path}" alt="${filename}" class="w-full h-auto rounded-lg">
        </div>
    `;

    document.body.appendChild(modal);
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

    // Tabs de resultados
    const tabTablas = document.getElementById('tab-tablas');
    const tabFiguras = document.getElementById('tab-figuras');

    if (tabTablas) {
        tabTablas.addEventListener('click', () => {
            document.querySelectorAll('.result-tab').forEach(t => {
                t.classList.remove('active', 'bg-emerald-500', 'text-white');
                t.classList.add('text-secondary');
            });
            tabTablas.classList.add('active', 'bg-emerald-500', 'text-white');
            tabTablas.classList.remove('text-secondary');
            loadResults('tablas');
        });
    }

    if (tabFiguras) {
        tabFiguras.addEventListener('click', () => {
            document.querySelectorAll('.result-tab').forEach(t => {
                t.classList.remove('active', 'bg-emerald-500', 'text-white');
                t.classList.add('text-secondary');
            });
            tabFiguras.classList.add('active', 'bg-emerald-500', 'text-white');
            tabFiguras.classList.remove('text-secondary');
            loadResults('figuras');
        });
    }

    // Eliminar todos los resultados
    const deleteAllBtn = document.getElementById('delete-all-results-btn');
    if (deleteAllBtn) {
        deleteAllBtn.addEventListener('click', deleteAllResults);
    }
});
