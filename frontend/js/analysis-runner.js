/**
 * GenomeHub - Analysis Runner
 *
 * Ejecuta scripts Python de análisis via servidor.py
 * Requiere servidor corriendo en AWS
 */

const AnalysisRunner = {
    async runScript(scriptName, organism = null) {
        const consolePanel = document.getElementById('analysis-console');
        const consoleOutput = document.getElementById('console-output');

        // Verificar servidor
        const hasServer = await NCBISearch.checkServer();
        if (!hasServer) {
            consolePanel.classList.remove('hidden');
            consoleOutput.innerHTML = `[ERROR] Servidor no disponible.\nEjecuta en AWS: python3 servidor.py\n\nLa búsqueda NCBI funciona sin servidor, pero los análisis requieren el servidor corriendo.`;
            showNotification('Se necesita servidor.py para ejecutar análisis', 'warning');
            return;
        }

        consolePanel.classList.remove('hidden');
        consoleOutput.innerHTML = `[${new Date().toLocaleTimeString()}] Ejecutando ${scriptName}...\n`;

        const btn = document.querySelector(`[data-script="${scriptName}"]`);
        if (btn) {
            btn.disabled = true;
            btn.textContent = 'Ejecutando...';
        }

        showNotification(`Ejecutando ${scriptName}...`, 'info');

        try {
            const payload = { script: scriptName };
            if (organism) payload.organism = organism;

            const response = await fetch('/api/run_analysis', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(payload)
            });

            const data = await response.json();

            if (!response.ok) {
                throw new Error(data.error || 'Error en el servidor');
            }

            consoleOutput.innerHTML += `\n${data.output}\n\n`;

            if (data.return_code === 0) {
                consoleOutput.innerHTML += `[OK] Completado en ${data.execution_time}s`;
                showNotification('Análisis completado correctamente', 'success');
                setTimeout(() => showSection('results'), 2000);
            } else {
                consoleOutput.innerHTML += `[ERROR] Finalizó con errores (código: ${data.return_code})`;
                showNotification('El análisis finalizó con errores', 'warning');
            }
        } catch (error) {
            consoleOutput.innerHTML += `\n[ERROR] ${error.message}`;
            showNotification('Error al ejecutar análisis', 'error');
        } finally {
            if (btn) {
                btn.disabled = false;
                btn.textContent = 'Ejecutar';
            }
        }
    }
};

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
            <div class="col-span-full text-center py-12 text-secondary">
                <div class="text-6xl mb-3">🖥️</div>
                <p class="font-semibold">Servidor requerido para ver resultados</p>
                <code class="bg-slate-900 text-emerald-400 px-3 py-1 rounded text-sm font-mono mt-2 inline-block">python3 servidor.py</code>
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
                    <div class="col-span-full text-center py-12 text-secondary">
                        <div class="text-6xl mb-3">📄</div>
                        <p>No hay resultados de tablas aún</p>
                        <p class="text-sm mt-2">Ejecuta un análisis primero</p>
                    </div>
                `;
            } else {
                tablasGrid.innerHTML = data.results.map(file => `
                    <div class="bg-card rounded-lg p-4 border border-slate-200 hover:shadow-md transition">
                        <div class="text-4xl mb-3">${file.extension === 'json' ? '📋' : '📊'}</div>
                        <h4 class="font-medium mb-1 text-sm text-primary line-clamp-2">${file.filename}</h4>
                        <p class="text-xs text-secondary mb-3">${file.size_kb} KB</p>
                        <a href="${file.path}" download class="block text-center px-3 py-2 bg-emerald-500 hover:bg-emerald-600 text-white text-xs rounded-lg transition">
                            Descargar
                        </a>
                    </div>
                `).join('');
            }
        } else {
            tablasContainer.classList.add('hidden');
            figurasContainer.classList.remove('hidden');

            if (data.results.length === 0) {
                figurasGrid.innerHTML = `
                    <div class="col-span-full text-center py-12 text-secondary">
                        <div class="text-6xl mb-3">📈</div>
                        <p>No hay gráficos generados aún</p>
                        <p class="text-sm mt-2">Ejecuta el script de visualizaciones</p>
                    </div>
                `;
            } else {
                figurasGrid.innerHTML = data.results.map(file => `
                    <div class="bg-card rounded-lg p-4 border border-slate-200 hover:shadow-md transition">
                        <img src="${file.path}" alt="${file.filename}" class="w-full h-48 object-contain rounded-lg mb-3 cursor-pointer hover:scale-105 transition" onclick="viewImage('${file.path}', '${file.filename}')">
                        <h4 class="font-medium text-sm mb-1 text-primary line-clamp-1">${file.filename}</h4>
                        <p class="text-xs text-secondary mb-2">${file.size_kb} KB</p>
                        <a href="${file.path}" download class="block text-center px-3 py-2 bg-emerald-500 hover:bg-emerald-600 text-white text-xs rounded-lg transition">
                            Descargar
                        </a>
                    </div>
                `).join('');
            }
        }
    } catch (error) {
        showNotification('Error al cargar resultados', 'error');
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
    document.querySelectorAll('.analysis-btn').forEach(btn => {
        btn.addEventListener('click', (e) => {
            const script = e.currentTarget.getAttribute('data-script');

            if (['analisis_genes'].includes(script)) {
                const organism = prompt('Seleccionar organismo:\n1 = E. coli K-12\n2 = Salmonella LT2');

                if (organism === '1' || organism === '2') {
                    const orgName = organism === '1' ? 'ecoli_k12' : 'salmonella_lt2';
                    AnalysisRunner.runScript(script, orgName);
                } else if (organism !== null) {
                    showNotification('Opción inválida. Usa 1 o 2', 'warning');
                }
            } else {
                AnalysisRunner.runScript(script);
            }
        });
    });

    const clearConsoleBtn = document.getElementById('clear-console');
    if (clearConsoleBtn) {
        clearConsoleBtn.addEventListener('click', () => {
            document.getElementById('console-output').innerHTML = '';
        });
    }

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
});
