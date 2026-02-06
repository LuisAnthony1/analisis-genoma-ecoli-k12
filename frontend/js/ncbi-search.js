/**
 * GenomeHub - NCBI Search Module
 *
 * Maneja búsqueda de genomas en NCBI, selección múltiple y descarga
 */

const NCBISearch = {
    /**
     * Buscar genomas en NCBI
     */
    async search(query, limit = 20) {
        const searchBtn = document.getElementById('search-btn');
        const searchLoading = document.getElementById('search-loading');
        const searchResults = document.getElementById('search-results');

        if (!query || query.trim().length < 3) {
            showNotification('La búsqueda debe tener al menos 3 caracteres', 'warning');
            return;
        }

        // UI: mostrar loading
        searchBtn.disabled = true;
        searchBtn.textContent = 'Buscando...';
        showLoading('search-loading');
        searchResults.innerHTML = '';

        try {
            const response = await fetch(
                `/backend/api/search_ncbi.php?query=${encodeURIComponent(query)}&limit=${limit}`
            );

            const data = await response.json();

            if (data.success && data.results && data.results.length > 0) {
                this.renderResults(data.results);
                showNotification(`${data.count} genomas encontrados`, 'success');
            } else {
                searchResults.innerHTML = `
                    <div class="col-span-full text-center py-12 bg-card rounded-xl border border-slate-200">
                        <div class="text-6xl mb-3">🔍</div>
                        <h3 class="text-xl font-semibold mb-2">No se encontraron resultados</h3>
                        <p class="text-secondary">Intenta con otro término de búsqueda</p>
                        <p class="text-sm text-secondary mt-2">Búsqueda: "${query}"</p>
                    </div>
                `;
            }
        } catch (error) {
            searchResults.innerHTML = `
                <div class="col-span-full bg-red-500/10 border border-red-500/20 rounded-xl p-6">
                    <h3 class="text-red-600 font-semibold mb-2">Error al buscar en NCBI</h3>
                    <p class="text-red-600 text-sm">${error.message}</p>
                </div>
            `;
            showNotification('Error al conectar con el servidor', 'error');
        } finally {
            searchBtn.disabled = false;
            searchBtn.textContent = 'Buscar';
            hideLoading('search-loading');
        }
    },

    /**
     * Renderizar resultados de búsqueda
     */
    renderResults(results) {
        const container = document.getElementById('search-results');

        container.innerHTML = results.map(genome => `
            <div class="bg-card rounded-xl p-5 border border-slate-200 hover:border-emerald-500 transition-all">
                <div class="flex items-start justify-between mb-3">
                    <div class="flex-1 pr-4">
                        <h3 class="font-semibold text-lg mb-1 text-primary">${genome.organism}</h3>
                        <p class="text-sm text-secondary mb-2 line-clamp-2">${genome.title}</p>
                        <div class="flex flex-wrap gap-2 text-xs">
                            <span class="px-2 py-1 bg-emerald-500/10 text-emerald-600 rounded-full">
                                ${genome.accession}
                            </span>
                            <span class="px-2 py-1 bg-cyan-500/10 text-cyan-600 rounded-full">
                                ${genome.length_mb} Mb
                            </span>
                            <span class="px-2 py-1 bg-purple-500/10 text-purple-600 rounded-full">
                                ${formatNumber(genome.length)} pb
                            </span>
                        </div>
                    </div>
                    <label class="flex items-center cursor-pointer">
                        <input
                            type="checkbox"
                            class="genome-select w-5 h-5 text-emerald-500 rounded border-gray-300 focus:ring-emerald-500"
                            data-genome='${JSON.stringify(genome)}'
                            onchange="NCBISearch.toggleSelection(this)"
                        >
                    </label>
                </div>
                <div class="text-xs text-secondary">
                    Actualizado: ${genome.update_date}
                </div>
            </div>
        `).join('');
    },

    /**
     * Manejar selección de genoma
     */
    toggleSelection(checkbox) {
        const genome = JSON.parse(checkbox.getAttribute('data-genome'));

        if (checkbox.checked) {
            AppState.selectedGenomes.add(JSON.stringify(genome));
        } else {
            AppState.selectedGenomes.delete(JSON.stringify(genome));
        }

        this.updateSelectedPanel();
    },

    /**
     * Actualizar panel de seleccionados
     */
    updateSelectedPanel() {
        const panel = document.getElementById('selected-genomes');
        const count = document.getElementById('selected-count');
        const list = document.getElementById('selected-list');

        const selected = Array.from(AppState.selectedGenomes).map(JSON.parse);
        count.textContent = selected.length;

        if (selected.length > 0) {
            panel.classList.remove('hidden');
            list.innerHTML = selected.map((g, idx) => `
                <div class="flex items-center justify-between bg-white dark:bg-slate-800 rounded-lg px-4 py-2 border border-emerald-500/20">
                    <div class="flex-1">
                        <span class="text-sm font-medium">${g.organism}</span>
                        <span class="text-xs text-secondary ml-2">(${g.accession})</span>
                    </div>
                    <button
                        onclick="NCBISearch.removeSelection('${g.accession}')"
                        class="text-red-500 hover:text-red-700 ml-2"
                    >
                        ✕
                    </button>
                </div>
            `).join('');
        } else {
            panel.classList.add('hidden');
        }
    },

    /**
     * Remover selección
     */
    removeSelection(accession) {
        const selected = Array.from(AppState.selectedGenomes).map(JSON.parse);
        const genome = selected.find(g => g.accession === accession);

        if (genome) {
            AppState.selectedGenomes.delete(JSON.stringify(genome));

            // Desmarcar checkbox
            const checkboxes = document.querySelectorAll('.genome-select');
            checkboxes.forEach(cb => {
                const g = JSON.parse(cb.getAttribute('data-genome'));
                if (g.accession === accession) {
                    cb.checked = false;
                }
            });

            this.updateSelectedPanel();
        }
    },

    /**
     * Descargar genomas seleccionados
     */
    async downloadSelected() {
        const selected = Array.from(AppState.selectedGenomes).map(JSON.parse);

        if (selected.length === 0) {
            showNotification('No hay genomas seleccionados', 'warning');
            return;
        }

        const btn = document.getElementById('download-selected-btn');
        btn.disabled = true;
        btn.textContent = 'Descargando...';

        let successCount = 0;
        let errorCount = 0;

        for (const genome of selected) {
            try {
                await this.downloadGenome(genome);
                successCount++;
            } catch (error) {
                errorCount++;
            }
        }

        btn.disabled = false;
        btn.textContent = 'Descargar Seleccionados';

        // Limpiar selección
        AppState.selectedGenomes.clear();
        document.querySelectorAll('.genome-select').forEach(cb => cb.checked = false);
        this.updateSelectedPanel();

        // Mostrar resumen
        if (successCount > 0) {
            showNotification(`${successCount} genomas descargados correctamente`, 'success');
            // Navegar a biblioteca
            setTimeout(() => showSection('library'), 2000);
        }

        if (errorCount > 0) {
            showNotification(`${errorCount} genomas fallaron`, 'error');
        }
    },

    /**
     * Descargar un genoma individual
     */
    async downloadGenome(genome) {
        showNotification(`Descargando ${genome.organism}...`, 'info');

        try {
            const response = await fetch('/backend/api/download_genome.php', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    accession_id: genome.accession,
                    organism_name: genome.organism
                })
            });

            const data = await response.json();

            if (!response.ok || !data.success) {
                throw new Error(data.error || 'Error al descargar');
            }

            showNotification(`✓ ${genome.organism} descargado`, 'success');
            return data;
        } catch (error) {
            showNotification(`✕ Error: ${genome.organism}`, 'error');
            throw error;
        }
    }
};

// =============================================================================
// EVENT LISTENERS
// =============================================================================

document.addEventListener('DOMContentLoaded', () => {
    const searchBtn = document.getElementById('search-btn');
    const searchInput = document.getElementById('search-input');
    const downloadBtn = document.getElementById('download-selected-btn');

    if (searchBtn && searchInput) {
        searchBtn.addEventListener('click', () => {
            const query = searchInput.value.trim();
            NCBISearch.search(query);
        });

        searchInput.addEventListener('keypress', (e) => {
            if (e.key === 'Enter') {
                searchBtn.click();
            }
        });
    }

    if (downloadBtn) {
        downloadBtn.addEventListener('click', () => NCBISearch.downloadSelected());
    }

    // Búsquedas rápidas
    document.querySelectorAll('.quick-search').forEach(btn => {
        btn.addEventListener('click', (e) => {
            const query = e.target.getAttribute('data-query');
            if (searchInput) {
                searchInput.value = query;
                NCBISearch.search(query);
            }
        });
    });
});
