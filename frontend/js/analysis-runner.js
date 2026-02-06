/**
 * GenomeHub - Analysis Runner
 *
 * Ejecuta scripts Python de análisis y maneja visualización de resultados
 */

const AnalysisRunner = {
    /**
     * Ejecutar script de análisis
     */
    async runScript(scriptName, organism = null) {
        const consolePanel = document.getElementById('analysis-console');
        const consoleOutput = document.getElementById('console-output');

        consolePanel.classList.remove('hidden');
        consoleOutput.innerHTML = `[${new Date().toLocaleTimeString()}] Ejecutando ${scriptName}...\n`;

        // Deshabilitar botón
        const btn = document.querySelector(`[data-script="${scriptName}"]`);
        if (btn) {
            btn.disabled = true;
            btn.textContent = 'Ejecutando...';
        }

        showNotification(`Ejecutando ${scriptName}...`, 'info');

        try {
            const payload = { script: scriptName };
            if (organism) payload.organism = organism;

            const response = await fetch('/backend/api/run_analysis.php', {
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
                consoleOutput.innerHTML += `[✓ OK] Script completado exitosamente en ${data.execution_time}s`;
                showNotification('Análisis completado correctamente', 'success');

                // Auto-navegar a resultados después de 2s
                setTimeout(() => showSection('results'), 2000);
            } else {
                consoleOutput.innerHTML += `[✕ ERROR] Script finalizó con errores (código: ${data.return_code})`;
                showNotification('El análisis finalizó con errores', 'warning');
            }
        } catch (error) {
            consoleOutput.innerHTML += `\n[✕ ERROR] ${error.message}`;
            showNotification('Error al ejecutar análisis', 'error');
        } finally {
            if (btn) {
                btn.disabled = false;
                btn.textContent = 'Ejecutar';
            }
        }
    }
};

/**
 * Cargar y mostrar resultados
 */
async function loadResults(type = 'tablas') {
    const tablasContainer = document.getElementById('results-tablas');
    const figurasContainer = document.getElementById('results-figuras');
    const tablasGrid = document.getElementById('tablas-grid');
    const figurasGrid = document.getElementById('figuras-grid');

    try {
        const response = await fetch(`/backend/api/get_results.php?type=${type}`);
        const data = await response.json();

        if (!response.ok) {
            throw new Error(data.error || 'Error al cargar resultados');
        }

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
                        <a
                            href="${file.path}"
                            download
                            class="block text-center px-3 py-2 bg-emerald-500 hover:bg-emerald-600 text-white text-xs rounded-lg transition"
                        >
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
                        <img
                            src="${file.path}"
                            alt="${file.filename}"
                            class="w-full h-48 object-contain rounded-lg mb-3 cursor-pointer hover:scale-105 transition"
                            onclick="viewImage('${file.path}', '${file.filename}')"
                        >
                        <h4 class="font-medium text-sm mb-1 text-primary line-clamp-1">${file.filename}</h4>
                        <p class="text-xs text-secondary mb-2">${file.size_kb} KB</p>
                        <a
                            href="${file.path}"
                            download
                            class="block text-center px-3 py-2 bg-emerald-500 hover:bg-emerald-600 text-white text-xs rounded-lg transition"
                        >
                            Descargar
                        </a>
                    </div>
                `).join('');
            }
        }
    } catch (error) {
        showNotification('Error al cargar resultados', 'error');
        console.error(error);
    }
}

/**
 * Ver imagen en tamaño completo
 */
function viewImage(path, filename) {
    const modal = document.createElement('div');
    modal.className = 'fixed inset-0 z-50 flex items-center justify-center p-4 bg-black/90 backdrop-blur-sm';
    modal.onclick = (e) => {
        if (e.target === modal || e.target.tagName === 'BUTTON') modal.remove();
    };

    modal.innerHTML = `
        <div class="max-w-6xl w-full">
            <div class="flex items-center justify-between mb-4">
                <h3 class="text-white font-medium">${filename}</h3>
                <button class="text-white hover:text-gray-300 text-2xl">✕</button>
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
    // Botones de análisis
    document.querySelectorAll('.analysis-btn').forEach(btn => {
        btn.addEventListener('click', (e) => {
            const script = e.currentTarget.getAttribute('data-script');

            // Scripts que requieren seleccionar organismo
            if (['analisis_genes', 'analisis_codones'].includes(script)) {
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

    // Limpiar console
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
});
