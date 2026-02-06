/**
 * GenomeHub - NCBI Search Module
 *
 * Busca genomas directamente en NCBI E-utilities API desde el navegador.
 * No necesita servidor para buscar. Para descargar a datos/crudo/ usa servidor.py
 */

const NCBI_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';
const NCBI_EMAIL = '194522@unsaac.edu.pe';

const NCBISearch = {
    // Detectar si servidor.py está corriendo
    serverAvailable: null,

    async checkServer() {
        if (this.serverAvailable !== null) return this.serverAvailable;
        try {
            const r = await fetch('/api/genomes', { signal: AbortSignal.timeout(2000) });
            this.serverAvailable = r.ok;
        } catch {
            this.serverAvailable = false;
        }
        return this.serverAvailable;
    },

    /**
     * Buscar genomas en NCBI (directo desde navegador)
     */
    async search(query, limit = 20) {
        const searchBtn = document.getElementById('search-btn');
        const searchResults = document.getElementById('search-results');

        if (!query || query.trim().length < 3) {
            showNotification('La búsqueda debe tener al menos 3 caracteres', 'warning');
            return;
        }

        searchBtn.disabled = true;
        searchBtn.textContent = 'Buscando en NCBI...';
        showLoading('search-loading');
        searchResults.innerHTML = '';

        try {
            // PASO 1: Buscar IDs en NCBI
            // Detectar si es un accession ID (ej: U00096, NC_000913, CP018008) o nombre de organismo
            const isAccession = /^[A-Z]{1,2}_?\d{4,}/i.test(query);
            let searchTerm;
            if (isAccession) {
                // Buscar por accession ID directamente
                searchTerm = `${query}[Accession]`;
            } else {
                // Buscar por nombre de organismo con wildcard
                searchTerm = `${query}*[Organism] AND (complete genome[Title] OR complete sequence[Title])`;
            }
            const searchUrl = `${NCBI_BASE}/esearch.fcgi?db=nucleotide&term=${encodeURIComponent(searchTerm)}&retmax=${limit}&retmode=json&email=${NCBI_EMAIL}&sort=relevance`;

            const searchResp = await fetch(searchUrl);
            const searchData = await searchResp.json();

            const ids = searchData.esearchresult?.idlist || [];

            if (ids.length === 0) {
                searchResults.innerHTML = `
                    <div class="col-span-full text-center py-12 bg-card rounded-xl border border-slate-200">
                        <div class="text-6xl mb-3">🔍</div>
                        <h3 class="text-xl font-semibold mb-2">No se encontraron genomas completos</h3>
                        <p class="text-secondary">Intenta con otro organismo. Ejemplo: "Escherichia coli"</p>
                    </div>
                `;
                hideLoading('search-loading');
                searchBtn.disabled = false;
                searchBtn.textContent = 'Buscar';
                return;
            }

            // PASO 2: Obtener detalles (summary) de cada ID
            const summaryUrl = `${NCBI_BASE}/esummary.fcgi?db=nucleotide&id=${ids.join(',')}&retmode=json&email=${NCBI_EMAIL}`;
            const summaryResp = await fetch(summaryUrl);
            const summaryData = await summaryResp.json();

            const results = [];
            const uids = summaryData.result?.uids || [];

            for (const uid of uids) {
                const s = summaryData.result[uid];
                if (!s) continue;

                const title = s.title || '';
                let organism = title.split(',')[0].trim();
                organism = organism.replace(/ complete genome$/i, '').replace(/ complete sequence$/i, '');

                const length = parseInt(s.slen || s.length || 0);

                results.push({
                    id: uid,
                    accession: s.accessionversion || s.caption || 'N/A',
                    title: title,
                    organism: organism,
                    length: length,
                    length_mb: (length / 1e6).toFixed(2),
                    update_date: s.updatedate || 'N/A',
                    create_date: s.createdate || 'N/A'
                });
            }

            if (results.length > 0) {
                this.renderResults(results);
                showNotification(`${results.length} genomas encontrados`, 'success');
            } else {
                searchResults.innerHTML = `
                    <div class="col-span-full text-center py-12 bg-card rounded-xl border border-slate-200">
                        <div class="text-6xl mb-3">🔍</div>
                        <h3 class="text-xl font-semibold mb-2">No se pudieron obtener detalles</h3>
                        <p class="text-secondary">Intenta de nuevo en unos segundos</p>
                    </div>
                `;
            }

        } catch (error) {
            console.error('Error buscando en NCBI:', error);
            searchResults.innerHTML = `
                <div class="col-span-full bg-red-500/10 border border-red-500/20 rounded-xl p-6">
                    <h3 class="text-red-600 font-semibold mb-2">Error al buscar en NCBI</h3>
                    <p class="text-red-600 text-sm">${error.message}</p>
                    <p class="text-sm text-secondary mt-2">Verifica tu conexión a internet</p>
                </div>
            `;
            showNotification('Error de conexión', 'error');
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

        container.innerHTML = results.map(genome => {
            // Escapar JSON para atributo data
            const genomeJson = JSON.stringify(genome).replace(/'/g, '&#39;');

            return `
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
                            data-genome='${genomeJson}'
                            onchange="NCBISearch.toggleSelection(this)"
                        >
                    </label>
                </div>
                <div class="flex items-center justify-between text-xs text-secondary">
                    <span>Actualizado: ${genome.update_date}</span>
                    <button
                        onclick='NCBISearch.downloadSingle(${genomeJson})'
                        class="px-3 py-1 bg-emerald-500 hover:bg-emerald-600 text-white rounded-lg transition text-xs"
                    >
                        Descargar
                    </button>
                </div>
            </div>
        `}).join('');
    },

    /**
     * Descargar un solo genoma (botón individual)
     */
    async downloadSingle(genome) {
        const hasServer = await this.checkServer();

        if (hasServer) {
            // Servidor disponible: descargar a datos/crudo/
            try {
                await this.downloadGenome(genome);
            } catch (error) {
                showNotification(`Error al descargar: ${error.message}`, 'error');
            }
        } else {
            // Sin servidor: descargar directo desde NCBI al navegador
            this.downloadFromNCBI(genome);
        }
    },

    /**
     * Descargar directo desde NCBI (sin servidor - descarga al navegador)
     */
    downloadFromNCBI(genome) {
        showNotification(`Descargando ${genome.organism}...`, 'info');

        // Abrir GenBank en nueva pestaña
        const gbUrl = `${NCBI_BASE}/efetch.fcgi?db=nucleotide&id=${genome.accession}&rettype=gbwithparts&retmode=text&email=${NCBI_EMAIL}`;
        const fastaUrl = `${NCBI_BASE}/efetch.fcgi?db=nucleotide&id=${genome.accession}&rettype=fasta&retmode=text&email=${NCBI_EMAIL}`;

        // Descargar GenBank
        const a = document.createElement('a');
        a.href = gbUrl;
        a.download = `${genome.accession}.gb`;
        a.target = '_blank';
        document.body.appendChild(a);
        a.click();
        a.remove();

        // Después de 2s, descargar FASTA
        setTimeout(() => {
            const b = document.createElement('a');
            b.href = fastaUrl;
            b.download = `${genome.accession}.fasta`;
            b.target = '_blank';
            document.body.appendChild(b);
            b.click();
            b.remove();
            showNotification(`${genome.organism} descargado`, 'success');
        }, 2000);
    },

    toggleSelection(checkbox) {
        const genome = JSON.parse(checkbox.getAttribute('data-genome'));

        if (checkbox.checked) {
            AppState.selectedGenomes.add(JSON.stringify(genome));
        } else {
            AppState.selectedGenomes.delete(JSON.stringify(genome));
        }

        this.updateSelectedPanel();
    },

    updateSelectedPanel() {
        const panel = document.getElementById('selected-genomes');
        const count = document.getElementById('selected-count');
        const list = document.getElementById('selected-list');

        const selected = Array.from(AppState.selectedGenomes).map(JSON.parse);
        count.textContent = selected.length;

        if (selected.length > 0) {
            panel.classList.remove('hidden');
            list.innerHTML = selected.map(g => `
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

    removeSelection(accession) {
        const selected = Array.from(AppState.selectedGenomes).map(JSON.parse);
        const genome = selected.find(g => g.accession === accession);

        if (genome) {
            AppState.selectedGenomes.delete(JSON.stringify(genome));

            document.querySelectorAll('.genome-select').forEach(cb => {
                const g = JSON.parse(cb.getAttribute('data-genome'));
                if (g.accession === accession) cb.checked = false;
            });

            this.updateSelectedPanel();
        }
    },

    async downloadSelected() {
        const selected = Array.from(AppState.selectedGenomes).map(JSON.parse);

        if (selected.length === 0) {
            showNotification('No hay genomas seleccionados', 'warning');
            return;
        }

        const hasServer = await this.checkServer();
        const btn = document.getElementById('download-selected-btn');
        btn.disabled = true;
        btn.textContent = 'Descargando...';

        if (hasServer) {
            // Con servidor: descargar a datos/crudo/
            let successCount = 0;
            let errorCount = 0;

            for (const genome of selected) {
                try {
                    await this.downloadGenome(genome);
                    successCount++;
                } catch {
                    errorCount++;
                }
            }

            if (successCount > 0) {
                showNotification(`${successCount} genomas descargados`, 'success');
                setTimeout(() => showSection('library'), 2000);
            }
            if (errorCount > 0) {
                showNotification(`${errorCount} descargas fallaron`, 'error');
            }
        } else {
            // Sin servidor: descargar desde NCBI directo
            for (const genome of selected) {
                this.downloadFromNCBI(genome);
                await new Promise(r => setTimeout(r, 3000)); // esperar entre descargas
            }
            showNotification('Descarga en progreso', 'info');
        }

        btn.disabled = false;
        btn.textContent = 'Descargar Seleccionados';

        AppState.selectedGenomes.clear();
        document.querySelectorAll('.genome-select').forEach(cb => cb.checked = false);
        this.updateSelectedPanel();
    },

    /**
     * Descargar genoma via servidor (a datos/crudo/)
     */
    async downloadGenome(genome) {
        showNotification(`Descargando ${genome.organism}...`, 'info');

        const response = await fetch('/api/download', {
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

        showNotification(`${genome.organism} descargado`, 'success');
        return data;
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
            NCBISearch.search(searchInput.value.trim());
        });

        searchInput.addEventListener('keypress', (e) => {
            if (e.key === 'Enter') searchBtn.click();
        });
    }

    if (downloadBtn) {
        downloadBtn.addEventListener('click', () => NCBISearch.downloadSelected());
    }

    document.querySelectorAll('.quick-search').forEach(btn => {
        btn.addEventListener('click', (e) => {
            const query = e.target.getAttribute('data-query');
            if (searchInput) {
                searchInput.value = query;
                NCBISearch.search(query);
            }
        });
    });

    // Detectar servidor al inicio
    NCBISearch.checkServer().then(available => {
        if (available) {
            console.log('[NCBI] Servidor detectado - descargas irán a datos/crudo/');
        } else {
            console.log('[NCBI] Sin servidor - búsqueda directa NCBI, descargas al navegador');
        }
    });
});
