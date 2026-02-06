/**
 * GenomeHub - Genome Library Manager
 *
 * Maneja la biblioteca de genomas descargados.
 * Con servidor: lee datos/crudo/ del servidor
 * Sin servidor: muestra mensaje indicando que se necesita servidor.py
 */

async function loadLibrary() {
    const loading = document.getElementById('library-loading');
    const grid = document.getElementById('library-grid');
    const empty = document.getElementById('library-empty');

    showLoading('library-loading');
    grid.innerHTML = '';
    empty.classList.add('hidden');

    const hasServer = await NCBISearch.checkServer();

    if (!hasServer) {
        hideLoading('library-loading');
        grid.innerHTML = `
            <div class="col-span-full text-center py-12 bg-card rounded-xl border border-slate-200">
                <div class="text-6xl mb-4">🖥️</div>
                <h3 class="text-xl font-semibold mb-2 text-primary">Servidor no disponible</h3>
                <p class="text-secondary mb-4">Se requiere conexión al servidor para ver genomas descargados</p>
            </div>
        `;
        return;
    }

    try {
        const response = await fetch('/api/genomes');
        const data = await response.json();

        if (data.success && data.genomes && data.genomes.length > 0) {
            AppState.libraryGenomes = data.genomes;
            renderLibrary(data.genomes);
        } else {
            empty.classList.remove('hidden');
        }
    } catch (error) {
        grid.innerHTML = `
            <div class="col-span-full bg-red-500/10 border border-red-500/20 rounded-xl p-6">
                <h3 class="text-red-600 font-semibold mb-2">Error al cargar biblioteca</h3>
                <p class="text-red-600 text-sm">${error.message}</p>
            </div>
        `;
    } finally {
        hideLoading('library-loading');
    }
}

function renderLibrary(genomes) {
    const grid = document.getElementById('library-grid');

    grid.innerHTML = genomes.map(genome => `
        <div class="bg-card rounded-xl p-6 border border-slate-200 hover:shadow-lg transition-all">
            <div class="flex items-start justify-between mb-4">
                <div class="w-12 h-12 rounded-xl bg-gradient-to-br from-emerald-500 to-cyan-500 flex items-center justify-center text-2xl">
                    🧬
                </div>
                <div class="flex items-center gap-2">
                    <span class="px-3 py-1 bg-emerald-500/10 text-emerald-600 text-xs font-medium rounded-full">
                        ${genome.size_mb} MB
                    </span>
                    <button
                        onclick="deleteGenome('${genome.basename}')"
                        class="w-8 h-8 flex items-center justify-center rounded-lg bg-red-500/10 hover:bg-red-500/20 text-red-500 transition"
                        title="Eliminar genoma"
                    >
                        🗑️
                    </button>
                </div>
            </div>

            <h3 class="font-semibold text-lg mb-2 text-primary line-clamp-2">
                ${genome.organism || genome.basename}
            </h3>

            <div class="space-y-1 text-sm text-secondary mb-4">
                <p><span class="font-medium">ID:</span> ${genome.accession_id}</p>
                <p><span class="font-medium">Longitud:</span> ${formatNumber(genome.length || 0)} pb</p>
                ${genome.num_features ? `<p><span class="font-medium">Features:</span> ${formatNumber(genome.num_features)}</p>` : ''}
                <p class="text-xs">Descargado: ${formatDate(genome.download_date || genome.modified)}</p>
            </div>

            <div class="flex gap-2">
                <button
                    onclick="analyzeGenome('${genome.basename}')"
                    class="flex-1 px-4 py-2 bg-emerald-500 hover:bg-emerald-600 text-white text-sm rounded-lg transition"
                >
                    Analizar
                </button>
                <button
                    onclick="viewGenomeDetails('${genome.basename}')"
                    class="px-4 py-2 bg-slate-100 dark:bg-slate-800 hover:bg-slate-200 dark:hover:bg-slate-700 text-sm rounded-lg transition"
                >
                    Ver
                </button>
            </div>
        </div>
    `).join('');
}

async function deleteGenome(basename) {
    if (!confirm('¿Eliminar este genoma y todos sus archivos?')) return;

    try {
        const resp = await fetch('/api/delete_genome', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ basename })
        });
        const data = await resp.json();

        if (data.success) {
            showNotification('Genoma eliminado', 'success');
            loadLibrary();
        } else {
            showNotification(data.error || 'Error al eliminar', 'error');
        }
    } catch (error) {
        showNotification('Error al eliminar', 'error');
    }
}

function analyzeGenome(basename) {
    AppState.analysisGenome = basename;
    showSection('analysis');
}

function viewGenomeDetails(basename) {
    const genome = AppState.libraryGenomes.find(g => g.basename === basename);
    if (!genome) {
        showNotification('Genoma no encontrado', 'error');
        return;
    }

    const modal = document.createElement('div');
    modal.className = 'fixed inset-0 z-50 flex items-center justify-center p-4 bg-black/50 backdrop-blur-sm';
    modal.onclick = (e) => { if (e.target === modal) modal.remove(); };

    const metadata = genome.metadata || {};
    const genomaInfo = metadata.genoma || {};

    modal.innerHTML = `
        <div class="bg-card rounded-xl max-w-2xl w-full max-h-[90vh] overflow-y-auto p-6 border border-slate-200">
            <div class="flex items-start justify-between mb-6">
                <h2 class="text-2xl font-bold text-primary">Detalles del Genoma</h2>
                <button onclick="this.closest('.fixed').remove()" class="text-secondary hover:text-primary text-2xl">✕</button>
            </div>

            <div class="space-y-4">
                <div>
                    <h3 class="text-sm font-semibold text-secondary mb-2">Información General</h3>
                    <div class="bg-slate-50 dark:bg-slate-800 rounded-lg p-4 space-y-2 text-sm">
                        <p><span class="font-medium">Organismo:</span> ${genomaInfo.organismo || genome.organism || 'N/A'}</p>
                        <p><span class="font-medium">Accession ID:</span> ${metadata.id_acceso || genome.accession_id}</p>
                        <p><span class="font-medium">Descripción:</span> ${genomaInfo.descripcion || genome.description || 'N/A'}</p>
                        <p><span class="font-medium">Longitud:</span> ${formatNumber(genomaInfo.longitud || genome.length || 0)} pares de bases</p>
                        <p><span class="font-medium">Features:</span> ${formatNumber(genomaInfo.num_features || genome.num_features || 0)} anotaciones</p>
                    </div>
                </div>

                <div>
                    <h3 class="text-sm font-semibold text-secondary mb-2">Archivos</h3>
                    <div class="bg-slate-50 dark:bg-slate-800 rounded-lg p-4 space-y-2 text-sm">
                        <p><span class="font-medium">GenBank:</span> ${genome.filename} (${genome.size_mb} MB)</p>
                        ${genome.has_fasta ? `<p><span class="font-medium">FASTA:</span> ${genome.basename}.fasta</p>` : ''}
                        ${genome.has_metadata ? `<p><span class="font-medium">Metadata:</span> metadata_${genome.basename}.json</p>` : ''}
                    </div>
                </div>

                <div>
                    <h3 class="text-sm font-semibold text-secondary mb-2">Descarga</h3>
                    <div class="bg-slate-50 dark:bg-slate-800 rounded-lg p-4 space-y-2 text-sm">
                        <p><span class="font-medium">Fecha:</span> ${formatDate(metadata.fecha_descarga || genome.download_date)}</p>
                        <p><span class="font-medium">Fuente:</span> ${metadata.fuente || 'NCBI'}</p>
                    </div>
                </div>
            </div>

            <div class="flex gap-2 mt-6">
                <button onclick="this.closest('.fixed').remove(); analyzeGenome('${genome.basename}')" class="flex-1 px-4 py-2 bg-emerald-500 hover:bg-emerald-600 text-white rounded-lg transition">
                    Analizar este Genoma
                </button>
                <button onclick="this.closest('.fixed').remove()" class="px-4 py-2 bg-slate-200 dark:bg-slate-700 hover:bg-slate-300 dark:hover:bg-slate-600 rounded-lg transition">
                    Cerrar
                </button>
            </div>
        </div>
    `;

    document.body.appendChild(modal);
}
